#include "itensor/all.h"
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <ctime>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <iomanip>

#include <vector>
using namespace itensor;



void evolvePlay(int N) {
	clock_t t0, t1;

	// set parameters
	Real energy = NAN, prevEnergy = NAN, tol = 1E-6;
	Args args = Args("Maxm",64,
					 "Minm",32,
					 "Cutoff",1E-15);

	// create sites
	auto sites = SpinHalf(N);

	// initialize Neel state
	auto state = InitState(sites);
	for (int i = 1; i <= N; ++i) {
		if(i%2 == 1)  state.set(i,"Up");
		else          state.set(i,"Dn");
	}
	auto psi = IQMPS(state);

	// create Hamiltonian
	auto ampo = AutoMPO(sites);
	for(int j = 1; j < N; j++) {
		ampo += 0.5,"S+",j,"S-",j+1;
		ampo += 0.5,"S-",j,"S+",j+1;
		ampo +=     "Sz",j,"Sz",j+1;
	}
	auto H = IQMPO(ampo);
	LocalMPO<IQTensor> PH(H);

	energy = overlap(psi,H,psi);
	printfln("energy=%.10f",energy);


	// set psi to be right orthogonolized matrix
	for (int i=N; i>1; i--) {
		auto& A = psi.Aref(i-1);
		auto& B = psi.Aref(i);
		IQTensor D;
		svd(A*B,A,D,B);
		A*=D;
	}

	energy = overlap(psi,H,psi);
	printfln("energy=%.10f",energy);

	auto type = IQGate::tImag;
	auto gates = std::vector<IQGate>();
	Real tstep = 0.1;

	for(int b = 1; b < N; ++b) {
		IQTensor hh = sites.op("Sz",b)*sites.op("Sz",b+1);
		hh += 0.5*sites.op("Sp",b)*sites.op("Sm",b+1);
		hh += 0.5*sites.op("Sm",b)*sites.op("Sp",b+1);
		gates.push_back(IQGate(sites,b,b+1,type,tstep/2.,hh));
	}

	/*printfln("%d\n", gates.size());
	for (int b=0; b < (int)gates.size()/2; b++) {
		Print(gates[b]);
	}

	Print(psi);
	Print(sites);*/
	for (int i=0; i<100; i++) {
		for (int b=1; b < N; b++) {
			auto& A = psi.Aref(b);
			auto& B = psi.Aref(b+1);
			auto AA = A*B * gates[b-1];
			AA.noprime();
			psi.svdBond(b,AA,Fromleft,PH,{"Cutoff",1E-10});
			//IQTensor D;
			//svd(AA,A,D,B,{"Cutoff",1E-10});
			//B*=D;
		}
		for (int b=N-1; b > 0; b--) {
			auto& A = psi.Aref(b);
			auto& B = psi.Aref(b+1);
			auto AA = A*B * gates[b-1];
			AA.noprime();
			psi.svdBond(b,AA,Fromright,PH,{"Cutoff",1E-10});
			//IQTensor D;
			//svd(AA,A,D,B,{"Cutoff",1E-10});
			//A*=D;
		}

		psi.normalize();
		energy = overlap(psi,H,psi);
		printfln("%d energy=%.10f",i,energy);
	}

}

void saveOutput(std::vector< std::vector<Real> > matrix, int n, int m) {
	std::ofstream file;
	file.open("./output.txt");
	for (int i=0; i<n; i++) {
		for (int j=0; j<m; j++)
			file << std::fixed << std::setprecision(15) << matrix[i][j] << "\t";
		file << "\n";
	}
	file.close();
}

/**
 * Create a set of gates for a trotter time step.
 * Supports the following expansions:
 * 	T2 - 2nd order trotter expansion
 * 	T4 - 4th order trotter expansion
 */
std::vector<std::vector<IQGate> > trotterGates(SiteSet sites, double trotterStep, std::string trotterType) {
	auto gates = std::vector<std::vector<IQGate> >();
	std::vector<double> tau;

	if (trotterType == "T2") {
		tau.push_back(trotterStep/2);
		tau.push_back(trotterStep/2);
	}
	else if (trotterType == "T4") {
		double tau1,tau3;
		tau1 = 1/(4-pow(4,1.0/3.0)) * trotterStep;
		tau3 = trotterStep - 4 * tau1;

		for (int t=0; t<10; t++) tau.push_back(tau1/2);
		tau[4] = tau3/2;
		tau[5] = tau3/2;
	}

	for (int t=0; t < (int) tau.size(); t++)
		gates.push_back(std::vector<IQGate>());

	for(int b = 1; b < sites.N(); b++) {
		IQTensor hh = sites.op("Sz",b)*sites.op("Sz",b+1);
		hh += 0.5*sites.op("Sp",b)*sites.op("Sm",b+1);
		hh += 0.5*sites.op("Sm",b)*sites.op("Sp",b+1);
		for (int t=0; t < (int) tau.size(); t++)
			gates[t].push_back(IQGate(sites,b,b+1,IQGate::tReal,tau[t],hh));
	}
	return gates;
}

/**
 * Advance one trotter time step
 */
void trotterAdvance(IQMPS *psi, std::vector<std::vector<IQGate> > gates, LocalMPO<IQTensor> PH, Args args, double *t_gate, double *t_svd, int *maxM) {
	clock_t t0, t1;
	int N = psi->sites().N();

	for (int t=0; t < (int) gates.size()/2; t++) {
		for (int b=1; b < N; b++) {
			auto& A = psi->Aref(b);
			auto& B = psi->Aref(b+1);
			t0 = std::clock();
			auto AA = (A*B * gates[2*t][b-1]).noprime();
			t1 = std::clock();
			*t_gate += double(t1-t0);

			for (int ind=0; ind < AA.inds().r(); ind++)
				if (AA.inds()[ind].m()>*maxM)
					*maxM = AA.inds()[ind].m();

			t0 = std::clock();
			psi->svdBond(b,AA,Fromleft,PH,args);
			t1 = std::clock();
			*t_svd += double(t1-t0);
		}
		for (int b=N-1; b > 0; b--) {
			auto& A = psi->Aref(b);
			auto& B = psi->Aref(b+1);
			t0 = std::clock();
			auto AA = (A*B * gates[2*t+1][b-1]).noprime();
			t1 = std::clock();
			*t_gate += double(t1-t0);

			for (int ind=0; ind < AA.inds().r(); ind++)
				if (AA.inds()[ind].m()>*maxM)
					*maxM = AA.inds()[ind].m();

			t0 = std::clock();
			psi->svdBond(b,AA,Fromright,PH,args);
			t1 = std::clock();
			*t_svd += double(t1-t0);
		}
	}
}

void trotter2order(int N, IQMPS *psi, std::vector<IQGate> gates, LocalMPO<IQTensor> PH, Args args, double *t_gate, double *t_svd, int *maxM) {
	clock_t t0, t1;

	for (int b=1; b < N; b++) {
		auto& A = psi->Aref(b);
		auto& B = psi->Aref(b+1);
		t0 = std::clock();
		auto AA = (A*B * gates[b-1]).noprime();
		t1 = std::clock();
		*t_gate += double(t1-t0);

		for (int ind=0; ind < AA.inds().r(); ind++)
			if (AA.inds()[ind].m()>*maxM)
				*maxM = AA.inds()[ind].m();

		t0 = std::clock();
		psi->svdBond(b,AA,Fromleft,PH,args);
		t1 = std::clock();
		*t_svd += double(t1-t0);
	}
	for (int b=N-1; b > 0; b--) {
		auto& A = psi->Aref(b);
		auto& B = psi->Aref(b+1);
		t0 = std::clock();
		auto AA = (A*B * gates[b-1]).noprime();
		t1 = std::clock();
		*t_gate += double(t1-t0);

		for (int ind=0; ind < AA.inds().r(); ind++)
			if (AA.inds()[ind].m()>*maxM)
				*maxM = AA.inds()[ind].m();

		t0 = std::clock();
		psi->svdBond(b,AA,Fromright,PH,args);
		t1 = std::clock();
		*t_svd += double(t1-t0);
	}
}

void trotterEvolve(IQMPS psi, int totTime, Real timeRes, int stepsPerRes, std::string trotterType) {
	// profiling stuff
	clock_t t0, t1;
	double t_gate = 0, t_svd = 0, t_mes = 0;
	int maxM = 0;

	// setup
	auto args = Args("Cutoff",1e-15	,"Maxm",64, "Minm",32);
	int TSteps = int(totTime / timeRes * stepsPerRes);
	auto sites = psi.sites();
	int N = sites.N();

	// set psi to be right orthogonolized matrix
	for (int i=N; i>1; i--) {
		auto& A = psi.Aref(i-1);
		auto& B = psi.Aref(i);
		IQTensor D;
		svd(A*B,A,D,B);
		A*=D;
	}

	// create Hamiltonian
	auto ampo = AutoMPO(sites);
	auto H = IQMPO(ampo);
	LocalMPO<IQTensor> PH(H);

	// create gates
	std::vector<std::vector<IQGate> > gates = trotterGates(sites, timeRes / stepsPerRes, trotterType);

	// create Sz
	auto SzOp = std::vector<IQMPO>();
	for (int j=1; j<=N; j++) {
		auto ampo = AutoMPO(sites);
		ampo += "Sz",j;
		SzOp.push_back(IQMPO(ampo));
	}

	std::vector< std::vector<Real> > Sz(TSteps / stepsPerRes, std::vector<Real>(N,0));

	for (int i=0; i<TSteps; i++) {
		if (i%(stepsPerRes*10) == 0) {
			printfln("step %d, maxM=%d",i/stepsPerRes,maxM);
			maxM = 0;
		}
		// save
		t0 = std::clock();
		if (i%stepsPerRes == 0)
			for (int j=1; j<=N; j++)
				Sz[i/stepsPerRes][j-1] = overlap(psi,SzOp[j-1],psi);
		t1 = std::clock();
		t_mes += double(t1-t0);

		// advance
		trotterAdvance(&psi, gates, PH, args, &t_gate, &t_svd, &maxM);
		/*trotter2order(N, &psi, gates[0], PH, args, &t_gate, &t_svd, &maxM);
		trotter2order(N, &psi, gates[0], PH, args, &t_gate, &t_svd, &maxM);
		trotter2order(N, &psi, gates[5], PH, args, &t_gate, &t_svd, &maxM);
		trotter2order(N, &psi, gates[0], PH, args, &t_gate, &t_svd, &maxM);
		trotter2order(N, &psi, gates[0], PH, args, &t_gate, &t_svd, &maxM);*/

	}
	saveOutput(Sz,TSteps/stepsPerRes, N);

	printfln("t_gate: %f",t_gate/CLOCKS_PER_SEC);
	printfln("t_svd: %f",t_svd/CLOCKS_PER_SEC);
	printfln("t_mes: %f",t_mes/CLOCKS_PER_SEC);
	printfln("%f\t%f\t%f",t_gate/CLOCKS_PER_SEC,t_svd/CLOCKS_PER_SEC,t_mes/CLOCKS_PER_SEC);
}





void evolveT2(int N, int totTime, Real timeRes, int stepsPerRes) {
	/*auto sites = SpinHalf(N);
	auto state = InitState(sites);
	for (int i = 1; i <= N; ++i) {
		if(i == N/2)  state.set(i,"Dn");
		//if(i%2 == 0)  state.set(i,"Dn"); // Neel State
		else          state.set(i,"Up");
	}
	auto psi = IQMPS(state);

	trotterEvolve(psi, totTime, timeRes, stepsPerRes, "T2");*/

	clock_t t0, t1;
	double t_gate = 0, t_svd = 0, t_mes = 0;
	int maxM = 0;

	auto args = Args("Cutoff",1e-15	,"Maxm",64, "Minm",32);
	int TSteps = int(totTime / timeRes * stepsPerRes);
	Real tstep = timeRes / stepsPerRes;

	// create psi(0)
	auto sites = SpinHalf(N);
	auto state = InitState(sites);
	for (int i = 1; i <= N; ++i) {
		if(i == N/2)  state.set(i,"Dn");
		//if(i%2 == 0)  state.set(i,"Dn"); // Neel State
		else          state.set(i,"Up");
	}
	auto psi = IQMPS(state);
	// set psi to be right orthogonolized matrix
	for (int i=N; i>1; i--) {
		auto& A = psi.Aref(i-1);
		auto& B = psi.Aref(i);
		IQTensor D;
		svd(A*B,A,D,B);
		A*=D;
	}

	// create Hamiltonian
	auto ampo = AutoMPO(sites);
	for(int j = 1; j < N; ++j) {
		ampo += 0.5,"S+",j,"S-",j+1;
		ampo += 0.5,"S-",j,"S+",j+1;
		ampo +=     "Sz",j,"Sz",j+1;
	}
	auto H = IQMPO(ampo);
	LocalMPO<IQTensor> PH(H);

	// create Sz
	auto SzOp = std::vector<IQMPO>();
	for (int j=1; j<=N; j++) {
		auto ampo = AutoMPO(sites);
		ampo += "Sz",j;
		SzOp.push_back(IQMPO(ampo));
	}

	// create gates
	auto gates = std::vector<IQGate>();
	for(int b = 1; b < N; b++) {
		IQTensor hh = sites.op("Sz",b)*sites.op("Sz",b+1);
		hh += 0.5*sites.op("Sp",b)*sites.op("Sm",b+1);
		hh += 0.5*sites.op("Sm",b)*sites.op("Sp",b+1);
		gates.push_back(IQGate(sites,b,b+1,IQGate::tReal,tstep/2.,hh));
	}

	std::vector< std::vector<Real> > Sz(TSteps / stepsPerRes, std::vector<Real>(N,0));

	for (int i=0; i<TSteps; i++) {
		if (i%(stepsPerRes*10) == 0) {
			printfln("step %d, maxM=%d",i/stepsPerRes,maxM);
			maxM = 0;
		}
		// save
		t0 = std::clock();
		if (i%stepsPerRes == 0)
			for (int j=1; j<=N; j++)
				Sz[i/stepsPerRes][j-1] = overlap(psi,SzOp[j-1],psi);
		t1 = std::clock();
		t_mes += double(t1-t0);

		// advance
		trotter2order(N, &psi, gates, PH, args, &t_gate, &t_svd, &maxM);
	}

	saveOutput(Sz,TSteps/stepsPerRes, N);
	printfln("t_gate: %f",t_gate/CLOCKS_PER_SEC);
	printfln("t_svd: %f",t_svd/CLOCKS_PER_SEC);
	printfln("t_mes: %f",t_mes/CLOCKS_PER_SEC);
	printfln("%f\t%f\t%f",t_gate/CLOCKS_PER_SEC,t_svd/CLOCKS_PER_SEC,t_mes/CLOCKS_PER_SEC);

}

void evolveT4(int N, int totTime, Real timeRes, int stepsPerRes) {
	auto sites = SpinHalf(N);
	auto state = InitState(sites);
	for (int i = 1; i <= N; ++i) {
		if(i == N/2)  state.set(i,"Dn");
		//if(i%2 == 0)  state.set(i,"Dn"); // Neel State
		else          state.set(i,"Up");
	}
	auto psi = IQMPS(state);

	trotterEvolve(psi, totTime, timeRes, stepsPerRes, "T4");

	/*
	clock_t t0, t1;
	double t_gate = 0, t_svd = 0, t_mes = 0;
	int maxM = 0;

	auto args = Args("Cutoff",1e-15	,"Maxm",64, "Minm",32);
	int TSteps = int(totTime / timeRes * stepsPerRes);
	Real tstep = timeRes / stepsPerRes;

	// create psi(0)
	auto sites = SpinHalf(N);
	auto state = InitState(sites);
	for (int i = 1; i <= N; ++i) {
		if(i == N/2)  state.set(i,"Dn");
		else          state.set(i,"Up");
	}
	auto psi = IQMPS(state);
	// set psi to be right orthogonolized matrix
	for (int i=N; i>1; i--) {
		auto& A = psi.Aref(i-1);
		auto& B = psi.Aref(i);
		IQTensor D;
		svd(A*B,A,D,B);
		A*=D;
	}

	// create Hamiltonian
	auto ampo = AutoMPO(sites);
	for(int j = 1; j < N; ++j) {
		ampo += 0.5,"S+",j,"S-",j+1;
		ampo += 0.5,"S-",j,"S+",j+1;
		ampo +=     "Sz",j,"Sz",j+1;
	}
	auto H = IQMPO(ampo);
	LocalMPO<IQTensor> PH(H);

	// create Sz
	auto SzOp = std::vector<IQMPO>();
	for (int j=1; j<=N; j++) {
		auto ampo = AutoMPO(sites);
		ampo += "Sz",j;
		SzOp.push_back(IQMPO(ampo));
	}

	// create gates
	auto gates1 = std::vector<IQGate>();
	auto gates3 = std::vector<IQGate>();
	double tau1,tau3;
	tau1 = 1/(4-pow(4,1.0/3.0)) * tstep;
	tau3 = tstep - 4 * tau1;
	for(int b = 1; b < N; b++) {
		IQTensor hh = sites.op("Sz",b)*sites.op("Sz",b+1);
		hh += 0.5*sites.op("Sp",b)*sites.op("Sm",b+1);
		hh += 0.5*sites.op("Sm",b)*sites.op("Sp",b+1);
		gates1.push_back(IQGate(sites,b,b+1,IQGate::tReal,tau1/2,hh));
		gates3.push_back(IQGate(sites,b,b+1,IQGate::tReal,tau3/2,hh));
	}

	std::vector< std::vector<Real> > Sz(TSteps / stepsPerRes, std::vector<Real>(N,0));

	for (int i=0; i<TSteps; i++) {
		if (i%(stepsPerRes*10) == 0) {
			printfln("step %d, maxM=%d",i/stepsPerRes,maxM);
			maxM = 0;
		}
		// save
		t0 = std::clock();
		if (i%stepsPerRes == 0)
			for (int j=1; j<=N; j++)
				Sz[i/stepsPerRes][j-1] = overlap(psi,SzOp[j-1],psi);
		t1 = std::clock();
		t_mes += double(t1-t0);

		// advance
		trotter2order(N, &psi, gates1, PH, args, &t_gate, &t_svd, &maxM);
		trotter2order(N, &psi, gates1, PH, args, &t_gate, &t_svd, &maxM);
		trotter2order(N, &psi, gates3, PH, args, &t_gate, &t_svd, &maxM);
		trotter2order(N, &psi, gates1, PH, args, &t_gate, &t_svd, &maxM);
		trotter2order(N, &psi, gates1, PH, args, &t_gate, &t_svd, &maxM);
	}

	saveOutput(Sz,TSteps/stepsPerRes, N);
	printfln("t_gate: %f",t_gate/CLOCKS_PER_SEC);
	printfln("t_svd: %f",t_svd/CLOCKS_PER_SEC);
	printfln("t_mes: %f",t_mes/CLOCKS_PER_SEC);
	printfln("%f\t%f\t%f",t_gate/CLOCKS_PER_SEC,t_svd/CLOCKS_PER_SEC,t_mes/CLOCKS_PER_SEC);
	*/
}

void evolveT4b(int N, int totTime, Real timeRes, int stepsPerRes) {
	clock_t t0, t1;
	double t_gate = 0, t_svd = 0, t_mes = 0;
	int maxM = 0;

	auto args = Args("Cutoff",1e-15	,"Maxm",64, "Minm",32);
	int TSteps = int(totTime / timeRes * stepsPerRes);
	Real tstep = timeRes / stepsPerRes;

	// create psi(0)
	auto sites = SpinHalf(N);
	auto state = InitState(sites);
	for (int i = 1; i <= N; ++i) {
		if(i == N/2)  state.set(i,"Dn");
		else          state.set(i,"Up");
	}
	auto psi = IQMPS(state);
	// set psi to be right orthogonolized matrix
	for (int i=N; i>1; i--) {
		auto& A = psi.Aref(i-1);
		auto& B = psi.Aref(i);
		IQTensor D;
		svd(A*B,A,D,B);
		A*=D;
	}

	// create Hamiltonian
	auto ampo = AutoMPO(sites);
	for(int j = 1; j < N; ++j) {
		ampo += 0.5,"S+",j,"S-",j+1;
		ampo += 0.5,"S-",j,"S+",j+1;
		ampo +=     "Sz",j,"Sz",j+1;
	}
	auto H = IQMPO(ampo);
	LocalMPO<IQTensor> PH(H);

	// create Sz
	auto SzOp = std::vector<IQMPO>();
	for (int j=1; j<=N; j++) {
		auto ampo = AutoMPO(sites);
		ampo += "Sz",j;
		SzOp.push_back(IQMPO(ampo));
	}

	// create gates
	auto gates = std::vector<std::vector<IQGate> >();
	std::vector<double> tau;
	double tau1,tau3;
	tau1 = 1/(4-pow(4,1.0/3.0)) * tstep;
	tau3 = tstep - 4 * tau1;

	tau.push_back(tau1/2);
	tau.push_back(tau1); //t1+t2
	tau.push_back((tau1+tau3)/2); //t2+t3
	tau.push_back((tau1+tau3)/2); //t3+t2
	tau.push_back(tau1); // t2+t1
	tau.push_back(tau1/2);
	for (int t=0; t < (int) tau.size(); t++)
		gates.push_back(std::vector<IQGate>());

	for(int b = 1; b < N; b++) {
		IQTensor hh = sites.op("Sz",b)*sites.op("Sz",b+1);
		hh += 0.5*sites.op("Sp",b)*sites.op("Sm",b+1);
		hh += 0.5*sites.op("Sm",b)*sites.op("Sp",b+1);
		for (int t=0; t < (int) tau.size(); t++)
			gates[t].push_back(IQGate(sites,b,b+1,IQGate::tReal,tau[t],hh));
	}

	std::vector< std::vector<Real> > Sz(TSteps / stepsPerRes, std::vector<Real>(N,0));

	for (int i=0; i<TSteps; i++) {
		if (i%(stepsPerRes*10) == 0) {
			printfln("step %d, maxM=%d",i/stepsPerRes,maxM);
			maxM = 0;
		}
		// save
		t0 = std::clock();
		if (i%stepsPerRes == 0)
			for (int j=1; j<=N; j++)
				Sz[i/stepsPerRes][j-1] = overlap(psi,SzOp[j-1],psi);
		t1 = std::clock();
		t_mes += double(t1-t0);

		// advance
		for (int t=0; t< (int) tau.size()/2; t++) {
			for (int b=1; b < N; b++) {
				auto& A = psi.Aref(b);
				auto& B = psi.Aref(b+1);
				t0 = std::clock();
				auto AA = (A*B * gates[2*t][b-1]).noprime();
				t1 = std::clock();
				t_gate += double(t1-t0);

				for (int ind=0; ind < AA.inds().r(); ind++)
					if (AA.inds()[ind].m()>maxM)
						maxM = AA.inds()[ind].m();

				t0 = std::clock();
				psi.svdBond(b,AA,Fromleft,PH,args);
				t1 = std::clock();
				t_svd += double(t1-t0);
			}
			for (int b=N-1; b > 0; b--) {
				auto& A = psi.Aref(b);
				auto& B = psi.Aref(b+1);
				t0 = std::clock();
				auto AA = (A*B * gates[2*t+1][b-1]).noprime();
				t1 = std::clock();
				t_gate += double(t1-t0);

				for (int ind=0; ind < AA.inds().r(); ind++)
					if (AA.inds()[ind].m()>maxM)
						maxM = AA.inds()[ind].m();

				t0 = std::clock();
				psi.svdBond(b,AA,Fromright,PH,args);
				t1 = std::clock();
				t_svd += double(t1-t0);
			}
		}
	}

	saveOutput(Sz,TSteps/stepsPerRes, N);
	printfln("t_gate: %f",t_gate/CLOCKS_PER_SEC);
	printfln("t_svd: %f",t_svd/CLOCKS_PER_SEC);
	printfln("t_mes: %f",t_mes/CLOCKS_PER_SEC);
	printfln("%f\t%f\t%f",t_gate/CLOCKS_PER_SEC,t_svd/CLOCKS_PER_SEC,t_mes/CLOCKS_PER_SEC);
}


int main(int argc, char* argv[]) {
	int N = 16;
	int totTime = 10;
	double timeRes = 0.1;
	int stepsPerRes = 1;
	std::string cmnd = "";
	if (argc > 1)
		cmnd = argv[1];
	if (argc > 5) {
		std::sscanf(argv[2],"%d", &N);
		std::sscanf(argv[3],"%d", &totTime);
		std::sscanf(argv[4],"%lf", &timeRes);
		std::sscanf(argv[5],"%d", &stepsPerRes);
	}
	clock_t t0, t1;
	time_t T0, T1;
	t0 = std::clock();
	T0 = std::time(0);
	if (cmnd == "T2")
		evolveT2(N, totTime, timeRes, stepsPerRes);
	else if (cmnd == "T4")
		evolveT4(N, totTime, timeRes, stepsPerRes);
	else if (cmnd == "T4b")
		evolveT4b(N, totTime, timeRes, stepsPerRes);
	t1 = std::clock();
	T1 = std::time(0);
	printfln("TOTAL TIME: %f or %f SECONDS",double(t1-t0)/CLOCKS_PER_SEC, T1-T0);
	printfln("%s\t%d\t%f\t%f\t%f\t%f",cmnd, N, totTime, timeRes, stepsPerRes);
	printfln("%f\t%f",double(t1-t0)/CLOCKS_PER_SEC, T1-T0);
	return 0;
}
