#include "itensor/all.h"
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <ctime>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <sys/stat.h>

#include <vector>
#include "TnConfig.h"
#include "TnMeasurement.h"
using namespace itensor;


class trotterMeasurements {
private:
	int N;
public:
	TnMeasurement<int> maxM = TnMeasurement<int>(MEAS_MAX_VALUE);
	TnMeasurement<Real> maxTruncError = TnMeasurement<double>(MEAS_MAX_VALUE);
	TnMeasurement<Real> totalTruncError = TnMeasurement<double>(MEAS_INCREMENTAL);
	TnMeasurement<Real> maxEntropy = TnMeasurement<double>(MEAS_MAX_VALUE);
	std::vector<TnMeasurement<Real> > Sz;

	trotterMeasurements (int N) : N(N){
		for (int i=0; i<N; i++)
			Sz.push_back(TnMeasurement<double>(MEAS_STATE_LESS));
	}
};

void saveOutput(Eigen::MatrixXd matrix, std::string fname)  {
	std::ofstream file;
	file.open(fname.c_str());
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
std::vector<std::vector<IQGate> > trotterGates(SiteSet sites, double trotterStep, TnConfig config) {

	auto gates = std::vector<std::vector<IQGate> >();
	std::vector<double> tau;

	std::string method = config.getEvMethod();
	if (method == "T2" || method == "ExT2") {
		tau.push_back(trotterStep/2);
		tau.push_back(trotterStep/2);
	}
	else if (method == "T4" || method == "ExT4") {
		double tau1,tau3;
		tau1 = 1/(4-pow(4,1.0/3.0)) * trotterStep;
		tau3 = trotterStep - 4 * tau1;

		for (int t=0; t<10; t++) tau.push_back(tau1/2);
		tau[4] = tau3/2;
		tau[5] = tau3/2;
	}

	for (int t=0; t < (int) tau.size(); t++)
		gates.push_back(std::vector<IQGate>());

	double Jxy = config.getHJxy(), Jz = config.getHJz();//, B = config.getHB();
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
void trotterAdvance(IQMPS *psi, std::vector<std::vector<IQGate> > gates, LocalMPO<IQTensor> PH, Args args, double *t_gate, double *t_svd, trotterMeasurements* measurements) {
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
				measurements->maxM.update(AA.inds()[ind].m());

			t0 = std::clock();
			auto spec = psi->svdBond(b,AA,Fromleft,PH,args);
			t1 = std::clock();
			*t_svd += double(t1-t0);
			measurements->totalTruncError.update(spec.truncerr());
			measurements->maxTruncError.update(spec.truncerr());
			// calculate entropy
            Real S = 0;
            for(auto& p : spec.eigsKept())
                if(p > 1E-13) S -= p*log(p);
			measurements->maxEntropy.update(S);

		}
		for (int b=N-1; b > 0; b--) {
			auto& A = psi->Aref(b);
			auto& B = psi->Aref(b+1);
			t0 = std::clock();
			auto AA = (A*B * gates[2*t+1][b-1]).noprime();
			t1 = std::clock();
			*t_gate += double(t1-t0);

			for (int ind=0; ind < AA.inds().r(); ind++)
				measurements->maxM.update(AA.inds()[ind].m());

			t0 = std::clock();
			auto spec = psi->svdBond(b,AA,Fromright,PH,args);
			t1 = std::clock();
			*t_svd += double(t1-t0);
			measurements->totalTruncError.update(spec.truncerr());
			measurements->maxTruncError.update(spec.truncerr());
			// calculate entropy
			Real S = 0;
			for(auto& p : spec.eigsKept())
				if(p > 1E-13) S -= p*log(p);
			measurements->maxEntropy.update(S);
		}
	}
}

trotterMeasurements trotterEvolve(IQMPS psi, int totTime, Real timeRes, int stepsPerRes, TnConfig config) {
	// profiling stuff
	clock_t t0, t1;
	double t_gate = 0, t_svd = 0, t_mes = 0;

	// setup
	auto args = Args("Cutoff",config.getCutoff(),"Maxm",config.getMaxm(), "Minm",config.getMinm());
	int TSteps = int(totTime / timeRes * stepsPerRes);
	auto sites = psi.sites();
	int N = sites.N();
	double currTime = 0;
	trotterMeasurements measurements(N);

	// set psi to be right orthogonolized matrix
	for (int i=N; i>1; i--) {
		auto& A = psi.Aref(i-1);
		auto& B = psi.Aref(i);
		IQTensor D;
		svd(A*B,A,D,B);
		A*=D;
	}

	// create empty Hamiltonian just for compatibility
	auto ampo = AutoMPO(sites);
	auto H = IQMPO(ampo);
	LocalMPO<IQTensor> PH(H);

	// create gates
	std::vector<std::vector<IQGate> > gates = trotterGates(sites, timeRes / stepsPerRes, config);

	// create Sz
	auto SzOp = std::vector<IQMPO>();
	for (int j=1; j<=N; j++) {
		auto ampo = AutoMPO(sites);
		ampo += "Sz",j;
		SzOp.push_back(IQMPO(ampo));
	}

	for (int i=0; i<TSteps; i++) {
		if (i%(stepsPerRes*10) == 0) {
			std::cout << "step " << i/stepsPerRes << ", maxM=" << measurements.maxM.getLastValue()
					  << ", truncErr=" << measurements.maxTruncError.getLastValue()
					  << ", totalErr=" << measurements.totalTruncError.getLastValue() << std::endl;
		}

		// save
		t0 = std::clock();
		if (i%stepsPerRes == 0)
			for (int j=1; j<=N; j++)
				measurements.Sz[j-1].updateAndMeasure(0.5 - overlap(psi,SzOp[j-1],psi), currTime);
		t1 = std::clock();
		t_mes += double(t1-t0);
		measurements.maxM.measure(currTime);
		measurements.maxTruncError.measure(currTime);
		measurements.totalTruncError.measure(currTime);
		measurements.maxEntropy.measure(currTime);

		// advance
		trotterAdvance(&psi, gates, PH, args, &t_gate, &t_svd, &measurements);
		currTime += (timeRes / stepsPerRes);
	}

	printfln("t_gate: %f",t_gate/CLOCKS_PER_SEC);
	printfln("t_svd: %f",t_svd/CLOCKS_PER_SEC);
	printfln("t_mes: %f",t_mes/CLOCKS_PER_SEC);
	printfln("%f\t%f\t%f",t_gate/CLOCKS_PER_SEC,t_svd/CLOCKS_PER_SEC,t_mes/CLOCKS_PER_SEC);

	return measurements;
}


IQMPS stateFromConfig(TnConfig config) {
	auto sites = SpinHalf(config.getN());
	auto state = InitState(sites);
	for (int i = 1; i <= config.getN(); ++i) {
		if (config.isUp(i))
			state.set(i,"Dn");
		else
			state.set(i,"Up");
	}
	return IQMPS(state);
}

trotterMeasurements evolveAndSave(TnConfig config, std::string runName, int stepsPerRes) {
	// evolve
	auto psi = stateFromConfig(config);
	int totTime = config.getEvTime();
	double timeRes = config.getEvRes();
	trotterMeasurements measurements =  trotterEvolve(psi, totTime, timeRes, stepsPerRes , config);

	// save
	measurements.maxM.saveToFile(config.getOutputPath() + "/" + runName + "_maxM.txt");
	measurements.maxTruncError.saveToFile(config.getOutputPath() + "/" + runName + "_maxTruncError.txt");
	measurements.totalTruncError.saveToFile(config.getOutputPath() + "/" + runName + "_totalTruncError.txt");
	measurements.maxEntropy.saveToFile(config.getOutputPath() + "/" + runName + "_maxEntropy.txt");
	saveToFile(config.getOutputPath() + "/" + runName + "_Sz.txt", measurements.Sz);

	return measurements;
}


void evolveExtrap(TnConfig config, const char* runName) {
	// copy configuration file
	mkdir(config.getOutputPath().c_str(),0777);
	std::ifstream f1 (config.getConfigPath().c_str(), std::fstream::binary);
	std::ofstream f2 ((config.getOutputPath() + "/" +  runName + ".config").c_str(), std::fstream::trunc|std::fstream::binary);
	f2 << f1.rdbuf();

	// evolve
	Eigen::Array3d stepsPerRes(1,2,4);
	std::vector<trotterMeasurements> measurements;
	for (int i = 0; i<stepsPerRes.rows(); i++)
		measurements.push_back(evolveAndSave(config, std::string(runName) + "_Ex" + std::to_string(i), (int) stepsPerRes(i)));

	// extrapolate - bad code - converts to Eigen:VectorXd and then back to std::vector because I am lazy
	int N = config.getN();
	double timeRes = config.getEvRes();
	// prepare extrapolation variables
	Eigen::Array3d x = pow(timeRes / stepsPerRes,4);
	std::vector<double> X(5);
	for(int i=0; i<=4; i++)
		X[i] = pow(x,i).sum();
	double den = 2*X[1]*X[2]*X[3] - X[1]*X[1]*X[4] + X[0]*X[2]*X[4] - X[0]*X[3]*X[3] - X[2]*X[2]*X[2];
	double X0Y = (X[2]*X[4] - X[3]*X[3]);
	double X1Y = (X[2]*X[3] - X[1]*X[4]);
	double X2Y = (X[1]*X[3] - X[2]*X[2]);
	// calculate extrapolated vector
	std::vector<TnMeasurement<double> > SzExtrap;
	for (int j=0; j<N; j++) {
		Eigen::VectorXd tmpExtrapValue = Eigen::VectorXd::Zero(measurements[0].Sz[j].size());
		for (int i = 0; i<stepsPerRes.rows(); i++)
			tmpExtrapValue += (X0Y + X1Y*x(i) + X2Y*x(i)*x(i))/den * Eigen::VectorXd::Map(measurements[i].Sz[j].getValues().data(),tmpExtrapValue.size());

		std::vector<double> tmpExtrapVec(tmpExtrapValue.data(), tmpExtrapValue.data() + tmpExtrapValue.size());
		SzExtrap.push_back(TnMeasurement<double>(tmpExtrapVec, measurements[0].Sz[j].getTimes()));
	}

	// save
	saveToFile(config.getOutputPath() + "/" + runName + "_extrap_Sz.txt", SzExtrap);

}


trotterMeasurements evolveOnce(TnConfig config, const char* runName) {
	// copy configuration file
	mkdir(config.getOutputPath().c_str(),0777);
	std::ifstream f1 (config.getConfigPath().c_str(), std::fstream::binary);
	std::ofstream f2 ((config.getOutputPath() + "/" +  runName + ".config").c_str(), std::fstream::trunc|std::fstream::binary);
	f2 << f1.rdbuf();

	return evolveAndSave(config, runName, 1);
}

int main(int argc, char* argv[]) {
	int N = 16;
	int totTime = 10;
	double timeRes = 0.1;

	clock_t t0, t1;
	time_t T0, T1;
	t0 = std::clock();
	T0 = std::time(0);

	if (argc > 2) {
		TnConfig config(argv[1]);
		printf("configuration\n=============\n");
		config.print();
		printf("==================================\n\n");
		std::string method = config.getEvMethod();
		if (method == "T2" || method == "T4")
			evolveOnce(config, argv[2]);

		else if (method == "ExT2" || method == "ExT4")
			evolveExtrap(config, argv[2]);
	}

	t1 = std::clock();
	T1 = std::time(0);
	printfln("TOTAL TIME: %f or %f SECONDS",double(t1-t0)/CLOCKS_PER_SEC, T1-T0);
	printfln("%f\t%f",double(t1-t0)/CLOCKS_PER_SEC, T1-T0);
	return 0;
}
