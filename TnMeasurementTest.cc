#include <stdio.h>
#include <stdlib.h>

#include "TnMeasurement.h"

class trotterMeasurements {
private:
	int N;
public:
	TnMeasurement<int> maxM = TnMeasurement<int>(MEAS_MAX_VALUE);
	TnMeasurement<double> maxTruncError = TnMeasurement<double>(MEAS_MAX_VALUE);
	TnMeasurement<double> totalTruncError = TnMeasurement<double>(MEAS_INCREMENTAL);
	TnMeasurement<double> maxEntropy = TnMeasurement<double>(MEAS_MAX_VALUE);
	std::vector<TnMeasurement<double> > Sz;

	trotterMeasurements (int N) : N(N){
		for (int i=0; i<N; i++)
			Sz.push_back(TnMeasurement<double>(MEAS_STATE_LESS));
	}
};


int main() {
	TnMeasurement<double> measD(MEAS_INCREMENTAL);
	measD.update(5.5);
	measD.update(5.324324325);
	measD.measure(0.1);
	measD.update(0.000001111);
	measD.measure(0.2);
	measD.saveToFile("./Tests/test1.txt");

	TnMeasurement<int> measI(MEAS_MAX_VALUE);
	measI.update(5);
	measI.update(2);
	measI.measure(0.1);
	measI.update(3);
	measI.measure(0.2);
	measI.saveToFile("./Tests/test2.txt");

	std::vector<TnMeasurement<double> > measVD;
	measVD.push_back(TnMeasurement<double>(MEAS_STATE_LESS));
	measVD.push_back(TnMeasurement<double>(MEAS_STATE_LESS));
	measVD[0].updateAndMeasure(0.3, 0.1);
	measVD[1].updateAndMeasure(0.315, 0.1);

	measVD[0].updateAndMeasure(0.113, 0.2);
	measVD[1].updateAndMeasure(0.352, 0.2);
	saveToFile("./Tests/test3.txt", measVD);

	int N = 10;
	trotterMeasurements measurements(N);
	for (double t=0; t<1; t+= 0.1)
		for (int i=0; i<N; i++)
			measurements.Sz[i].updateAndMeasure(i*t,t);
	saveToFile("./Tests/test4.txt", measurements.Sz);
	measurements.Sz[0].saveToFile("./Tests/test5.txt");
}
