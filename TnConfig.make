CPP = g++
CPP_COMP_FLAG = -std=c++11

EXEC = TnConfigTest.exe
OBJS = TnConfigTest.o

$(EXEC): $(OBJS)
	$(CPP) $(OBJS) -o $@

TnConfigTest.o: TnConfigTest.cc TnConfig.h
	$(CPP) $(CPP_COMP_FLAG) -c $*.cc