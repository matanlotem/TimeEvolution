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




void saveOutput(std::vector< std::vector<Real> > matrix, char *fname) {
	std::ofstream file;
	file.open(fname);
	for (int i=0; i< (int) matrix.size(); i++) {
		for (int j=0; j< (int) matrix[i].size(); j++)
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
	else if (trotterType == "T4b") {
		double tau1,tau3;
		tau1 = 1/(4-pow(4,1.0/3.0)) * trotterStep;
		tau3 = trotterStep - 4 * tau1;

		tau.push_back(tau1/2);
		tau.push_back(tau1); //t1+t2
		tau.push_back((tau1+tau3)/2); //t2+t3
		tau.push_back((tau1+tau3)/2); //t3+t2
		tau.push_back(tau1); // t2+t1
		tau.push_back(tau1/2);
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

std::vector< std::vector<Real> > trotterEvolve(IQMPS psi, int totTime, Real timeRes, int stepsPerRes, std::string trotterType) {
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
	}

	printfln("t_gate: %f",t_gate/CLOCKS_PER_SEC);
	printfln("t_svd: %f",t_svd/CLOCKS_PER_SEC);
	printfln("t_mes: %f",t_mes/CLOCKS_PER_SEC);
	printfln("%f\t%f\t%f",t_gate/CLOCKS_PER_SEC,t_svd/CLOCKS_PER_SEC,t_mes/CLOCKS_PER_SEC);

	return Sz;
}

void evolve1E(int N, int totTime, Real timeRes, int stepsPerRes, std::string trotterType) {
	auto sites = SpinHalf(N);
	auto state = InitState(sites);
	for (int i = 1; i <= N; ++i) {
		if(i == N/2)  state.set(i,"Dn");
		else          state.set(i,"Up");
	}
	auto psi = IQMPS(state);

	auto Sz = trotterEvolve(psi, totTime, timeRes, stepsPerRes, trotterType);
	saveOutput(Sz,"./output.txt");
}

void evolve2E(int N, int i1, int i2, int totTime, Real timeRes) {
	auto sites = SpinHalf(N);
	auto state = InitState(sites);
	for (int i = 1; i <= N; ++i) {
		if(i == i1 || i == i2)  state.set(i,"Dn");
		else          state.set(i,"Up");
	}
	auto psi = IQMPS(state);

	for(int stepsPerRes = 1; stepsPerRes <=4; stepsPerRes *= 2)
		trotterEvolve(psi, totTime, timeRes, stepsPerRes, "T4");
}

void evolveT2(int N, int totTime, Real timeRes, int stepsPerRes) {
	evolve1E(N, totTime, timeRes, stepsPerRes, "T2");
}

void evolveT4(int N, int totTime, Real timeRes, int stepsPerRes) {
	evolve1E(N, totTime, timeRes, stepsPerRes, "T4");
}

void evolveT4b(int N, int totTime, Real timeRes, int stepsPerRes) {
	evolve1E(N, totTime, timeRes, stepsPerRes, "T4b");
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
