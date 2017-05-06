#include "itensor/all.h"
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <ctime>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <iomanip>

#include <vector>
#include "TnConfig.h"
using namespace itensor;




void saveOutput(Eigen::MatrixXd matrix, const char *fname)  {
	std::ofstream file;
	file.open(fname);
	for (int i=0; i< (int) matrix.rows(); i++) {
		for (int j=0; j< (int) matrix.cols(); j++)
			file << std::fixed << std::setprecision(15) << matrix(i,j) << "\t";
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
	double Jxy = 1, Jz = 0;
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
		IQTensor hh = 0.5*Jxy*sites.op("Sp",b)*sites.op("Sm",b+1)
				    + 0.5*Jxy*sites.op("Sm",b)*sites.op("Sp",b+1)
					+ Jz*sites.op("Sz",b)*sites.op("Sz",b+1);
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

Eigen::MatrixXd trotterEvolve(IQMPS psi, int totTime, Real timeRes, int stepsPerRes, std::string trotterType) {
	// profiling stuff
	clock_t t0, t1;
	double t_gate = 0, t_svd = 0, t_mes = 0;
	int maxM = 0;

	// setup
	auto args = Args("Cutoff",1e-20,"Maxm",64, "Minm",2);
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

	Eigen::MatrixXd Sz(TSteps / stepsPerRes,N);
	for (int i=0; i<TSteps; i++) {
		if (i%(stepsPerRes*10) == 0) {
			printfln("step %d, maxM=%d",i/stepsPerRes,maxM);
			maxM = 0;
		}
		// save
		t0 = std::clock();
		if (i%stepsPerRes == 0)
			for (int j=1; j<=N; j++)
				Sz(i/stepsPerRes,j-1) = overlap(psi,SzOp[j-1],psi);
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

void initFromConfig(TnConfig config) {
	int N = config.getN();

	// create initial state
	auto sites = SpinHalf(N);
	auto state = InitState(sites);
	for (int i = 1; i <= N; ++i) {
		if (config.isUp(i))
			state.set(i,"Dn");
		else
			state.set(i,"Up");
	}
	auto psi = IQMPS(state);
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

	Eigen::Array3d stepsPerRes(1,2,4);
	//
	std::vector<Eigen::MatrixXd> Sz;
	for (int i = 0; i<stepsPerRes.rows(); i++) {
		Sz.push_back(trotterEvolve(psi, totTime, timeRes, int(stepsPerRes(i)) , "T4"));
		Sz[Sz.size()-1] = 0.5 * Eigen::MatrixXd(int(totTime / timeRes),N).setOnes() - Sz[Sz.size()-1];
	}

	// extrapolate
	Eigen::Array3d x = pow(timeRes / stepsPerRes,4);
	std::vector<double> X(5);
	for(int i=0; i<=4; i++)
		X[i] = pow(x,i).sum();
	double den = 2*X[1]*X[2]*X[3] - X[1]*X[1]*X[4] + X[0]*X[2]*X[4] - X[0]*X[3]*X[3] - X[2]*X[2]*X[2];
	double X0Y = (X[2]*X[4] - X[3]*X[3]);
	double X1Y = (X[2]*X[3] - X[1]*X[4]);
	double X2Y = (X[1]*X[3] - X[2]*X[2]);

	Sz.push_back(Eigen::MatrixXd(int(totTime / timeRes),N).setZero());
	for (int i = 0; i<stepsPerRes.rows(); i++)
		Sz[Sz.size()-1] += (X0Y + X1Y*x(i) + X2Y*x(i)*x(i))/den * Sz[i];

	// telemetries
	//abs((Sz[1]-Sz[0]).array())

	// save
	char tag[100];
	sprintf(tag,"Evolve2E_T4_%d_%d_%d__%d_%g",N, i1, i2, totTime, timeRes);
	char fname[1024];
	for (int i = 0; i<stepsPerRes.rows(); i++) {
		 sprintf(fname,"%s_%d.txt",tag,int(stepsPerRes(i)));
		 saveOutput(Sz[i],fname);
	}
	sprintf(fname,"%s_extrap.txt",tag);
	saveOutput(Sz[Sz.size()-1],fname);
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

	TnConfig config("./test.config");
	config.print();

	std::string cmnd = "";
	if (argc > 1)
		cmnd = argv[1];

	clock_t t0, t1;
	time_t T0, T1;
	t0 = std::clock();
	T0 = std::time(0);
	if (cmnd == "T2" || cmnd == "T4" || cmnd == "T4b") {
		int stepsPerRes = 1;
		if (argc > 5) {
			std::sscanf(argv[2],"%d", &N);
			std::sscanf(argv[3],"%d", &totTime);
			std::sscanf(argv[4],"%lf", &timeRes);
			std::sscanf(argv[5],"%d", &stepsPerRes);
		}
		if (cmnd == "T2")
			evolveT2(N, totTime, timeRes, stepsPerRes);
		else if (cmnd == "T4")
			evolveT4(N, totTime, timeRes, stepsPerRes);
		else if (cmnd == "T4b")
			evolveT4b(N, totTime, timeRes, stepsPerRes);
	}
	else if (cmnd == "2E") {
		int i1=0, i2=0;
		if (argc > 6) {
			std::sscanf(argv[2],"%d", &N);
			std::sscanf(argv[3],"%d", &i1);
			std::sscanf(argv[4],"%d", &i2);
			std::sscanf(argv[5],"%d", &totTime);
			std::sscanf(argv[6],"%lf", &timeRes);
		}
		evolve2E(N, i1, i2, totTime, timeRes);
	}

	t1 = std::clock();
	T1 = std::time(0);
	printfln("TOTAL TIME: %f or %f SECONDS",double(t1-t0)/CLOCKS_PER_SEC, T1-T0);
	printfln("%f\t%f",double(t1-t0)/CLOCKS_PER_SEC, T1-T0);
	return 0;
}
