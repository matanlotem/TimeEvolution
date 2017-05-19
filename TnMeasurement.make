CPP = g++
CPP_COMP_FLAG = -std=c++11

EXEC = TnMeasurementTest.exe
OBJS = TnMeasurementTest.o

$(EXEC): $(OBJS)
	$(CPP) $(OBJS) -o $@

TnMeasurementTest.o: TnMeasurementTest.cc TnMeasurement.h
	$(CPP) $(CPP_COMP_FLAG) -c $*.cc