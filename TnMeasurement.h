/*
 * TnMeasurement.h
 *
 *  Created on: May 17, 2017
 *      Author: Matan
 */

#ifndef TNMEASUREMENT_H_
#define TNMEASUREMENT_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

typedef enum measure_mode {
	MEAS_MAX_VALUE,
	MEAS_INCREMENTAL,
	MEAS_STATE_LESS,
	MEAS_ALL_DONE
} MeasureMode;

template <class T>
class TnMeasurement {
private:
	MeasureMode mode;
	T defaultValue;
	T currentValue;
	std::vector<T> values;
	std::vector<double> times;
public:
	TnMeasurement(MeasureMode mode);
	TnMeasurement(MeasureMode mode, T defaultValue);
	TnMeasurement(std::vector<T> values, std::vector<double> times);

	void update(T value);
	void measure(double time);
	void updateAndMeasure(T value, double time);
	int size();

	T getCurrentValue();
	T getLastValue();
	double getLastTime();
	T getValue(int i);
	double getTime(int i);
	std::vector<T> getValues();
	std::vector<double> getTimes();

	void saveToFile(std::string filename);
};

template <class T> void saveToFile(std::string filename, std::vector<TnMeasurement<T>> vec);



template <class T> TnMeasurement<T>::TnMeasurement(MeasureMode mode) : mode(mode) {
	defaultValue = 0;
	currentValue = defaultValue;
}

template <class T>  TnMeasurement<T>::TnMeasurement(MeasureMode mode, T defaultValue) :
		mode(mode), defaultValue(defaultValue) {
	currentValue = defaultValue;
}

template <class T> TnMeasurement<T>::TnMeasurement(std::vector<T> values, std::vector<double> times) :
	values(values), times(times) {
	mode = MEAS_ALL_DONE;
	defaultValue = 0;
	currentValue = 0;
}

template <class T> void TnMeasurement<T>::update(T value) {
	switch (mode) {
	case MEAS_MAX_VALUE:
		if (value > currentValue)
			currentValue = value;
		break;
	case MEAS_INCREMENTAL:
		currentValue += value;
		break;
	case MEAS_STATE_LESS:
		currentValue = value;
	default: break;
	}
}

template <class T> void TnMeasurement<T>::measure(double time) {
	times.push_back(time);
	values.push_back(currentValue);

	// reset currentValue
	switch (mode) {
	case MEAS_MAX_VALUE:
		currentValue = defaultValue;
		break;
	default: break;
	}
}

template <class T> void TnMeasurement<T>::updateAndMeasure(T value, double time) {
	update(value);
	measure(time);
}

template <class T> int TnMeasurement<T>::size() {
	return values.size();
}

template <class T> T TnMeasurement<T>::getCurrentValue() {
	return currentValue;
}

template <class T> T TnMeasurement<T>::getLastValue() {
	if (values.size() > 0)
		return values[values.size()-1];
	return defaultValue;
}

template <class T> double TnMeasurement<T>::getLastTime() {
	if (times.size() > 0)
		return times[times.size()-1];
	return 0;
}

template <class T> T TnMeasurement<T>::getValue(int i) {
	return values[i];
}

template <class T> double TnMeasurement<T>::getTime(int i) {
	return times[i];
}

template <class T> std::vector<T> TnMeasurement<T>::getValues() {
	return values;
}

template <class T> std::vector<double> TnMeasurement<T>::getTimes() {
	return times;
}

template <class T> void TnMeasurement<T>::saveToFile(std::string filename) {
	std::ofstream file;
    std::streamsize ss = std::cout.precision();
	file.open(filename.c_str());
	for (int i=0; i< (int) times.size(); i++) {
		file << std::defaultfloat << std::setprecision(ss) << times[i] << "\t";
		file << std::fixed << std::setprecision(15) << values[i] << std::endl;
	}
	file.close();

}

template <class T> void saveToFile(std::string filename, std::vector<TnMeasurement<T>> vec) {
	// validate dimensions
	if (vec.size() == 0) return;
	std::vector<double> times = vec[0].getTimes();
	for (int j=0; j< (int) vec.size(); j++) {
		if (vec[j].getValues().size() != times.size()) {
			std::cout << "Not all measurements of same dimension" << std::endl;
			return;
		}
	}

	// write to file
	std::ofstream file;
	std::streamsize ss = std::cout.precision();
	file.open(filename.c_str());
	for (int i=0; i < (int) times.size(); i++) {
		file << std::defaultfloat << std::setprecision(ss) << times[i] << "\t";
		for (int j=0; j < (int) vec.size(); j++)
			file << std::fixed << std::setprecision(15) << vec[j].getValue(i) << "\t";
		file << std::endl;
	}
	file.close();
}

#endif /* TNMEASUREMENT_H_ */
