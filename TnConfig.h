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
	int N;
	std::vector<bool> state;

	// Hamiltonian
	double HJxy;
	double HJz;
	double HB;

	// Time Evolution
	int EvTime;
	double EvRes;

	// Paths
	std::string outputPath;


public:
	TnConfig(const char* configPath);
	void print();

	int getN();
	bool isUp(int i);

	// Hamiltonian
	double getHJxy();
	double getHJz();
	double getHB();

	int getEvTime();
	double getEvRes();
};



TnConfig::TnConfig(const char* configPath) {
	char line[STR_LEN];
	std::string substr,param, valstr;
	std::ifstream file;
	file.open(configPath);

	// set default values
	N = 0;
	HJxy = 0;
	HJz = 0;
	HB = 0;
	EvTime = 0;
	EvRes = 0;

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
		else if (param == "Jxy")
			HJxy = atof(valstr.c_str());
		else if (param == "Jz")
			HJz = atof(valstr.c_str());
		else if (param == "B")
			HB = atof(valstr.c_str());
		else if (param == "EvTime")
			EvTime = atoi(valstr.c_str());
		else if (param == "EvRes")
			EvRes = atof(valstr.c_str());
		else if (param == "Output")
			outputPath = valstr;
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
	printf ("Evolve: Time=%d Res=%.1f\n", EvTime, EvRes);
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

int TnConfig::getEvTime() {
	return EvTime;
}

double TnConfig::getEvRes() {
	return EvRes;
}



#endif /* TNCONFIG_H_ */
