/*
 * TnConfig.h
 *
 *  Created on: Apr 11, 2017
 *      Author: Matan
 */

#ifndef TNCONFIG_H_
#define TNCONFIG_H_

#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cstring>

#define STR_LEN 1025

class TnConfig {
private:
	// state
	int N = 0;
	std::vector<bool> state;

	// Hamiltonian
	double HJxy = 1;
	double HJz = 1;
	double HB = 0;

	// Time Evolution
	int EvTime = 0;
	double EvRes = 0;
	std::string EvMethod;

	// DMRG
	int Maxm = 64;
	int Minm = 2;
	double Cutoff = 1e-20;

	// Paths
	std::string outputPath;
	std::string outputName;
	std::string configPath;


public:
	TnConfig(const char* path);
	void print();

	// State
	int getN();
	bool isUp(int i);

	// Hamiltonian
	double getHJxy();
	double getHJz();
	double getHB();

	// Time evolution
	int getEvTime();
	double getEvRes();
	std::string getEvMethod();

	// DMRG
	int getMaxm();
	int getMinm();
	double getCutoff();

	// Paths
	std::string getOutputPath();
	std::string getOutputName();
	std::string getConfigPath();

};



TnConfig::TnConfig(const char* path) : configPath(path) {
	char line[STR_LEN];
	std::string substr,param, valstr;
	std::ifstream file;
	file.open(configPath.c_str());

	while (file.getline(line,STR_LEN)) {
		// comment or empty line
		if (strlen(line) == 0 || line[0] == '#')
			continue;

		std::stringstream ss(line);
		getline(ss, param, ' ');
		getline(ss, valstr, ' ');

		if (param == "N") {
			N = atoi(valstr.c_str());
			state.resize(N,false);
		}
		else if (param == "Jxy")	HJxy = atof(valstr.c_str());
		else if (param == "Jz")		HJz = atof(valstr.c_str());
		else if (param == "B")		HB = atof(valstr.c_str());
		else if (param == "EvTime")	EvTime = atoi(valstr.c_str());
		else if (param == "EvRes")	EvRes = atof(valstr.c_str());
		else if (param == "EvMethod")	EvMethod = valstr;
		else if (param == "Output")	outputPath = valstr;
		else if (param == "Name")	outputName= valstr;
		else if (param == "Maxm")	Maxm = atoi(valstr.c_str());
		else if (param == "Minm")	Minm = atoi(valstr.c_str());
		else if (param == "Cutoff")	Cutoff = atof(valstr.c_str());
		// get state
		else if (param == "State") {
			// Neel state
			if (valstr == "neel") {
				for (int i=0; i<N; i++)
					state[i] = (i%2 == 0);
			}
			// only up spins specified
			else if (valstr ==  "up") {
				while (getline(ss, substr, ' ')) {
					state[atoi(substr.c_str())-1] = true;
				}
			}
			// all spins specified
			else if (valstr == "full") {
				getline(ss, substr);
				for (int i=0; i < (int) substr.length(); i++)
					state[i] = (substr[i] == '1');
			}
		}
	}
	file.close();
}

void TnConfig::print() {
	printf ("N = %d\n", N);
	printf ("State: ");
	for (int i=0; i<N; i++)
		printf("%d", (int) state[i]);
	printf("\n");
	printf ("H: Jxy=%.1f Jz=%.1f B=%.1f\n", HJxy, HJz, HB);
	printf ("Evolve: Method=%s Time=%d Res=%.2f\n", EvMethod.c_str(), EvTime, EvRes);
	printf ("DMRG: Maxm=%d Minm=%d Cutoff=%f\n", Maxm, Minm, Cutoff);
}

int TnConfig::getN() {
	return N;
}

bool TnConfig::isUp(int i) {
	return (i>0 && i<=N && state[i-1]);
}

// Hamiltonian
double TnConfig::getHJxy() {
	return HJxy;
}

double TnConfig::getHJz() {
	return HJz;
}

double TnConfig::getHB() {
	return HB;
}

// Time evolution
int TnConfig::getEvTime() {
	return EvTime;
}

double TnConfig::getEvRes() {
	return EvRes;
}

std::string TnConfig::getEvMethod() {
	return EvMethod;
}

// DMRG
int TnConfig::getMaxm() {
	return Maxm;
}

int TnConfig::getMinm() {
	return Minm;
}

double TnConfig::getCutoff() {
	return Cutoff;
}

// Paths
std::string TnConfig::getOutputPath() {
	return outputPath;
}

std::string TnConfig::getOutputName() {
	return outputName;
}

std::string TnConfig::getConfigPath() {
	return configPath;
}

#endif /* TNCONFIG_H_ */
