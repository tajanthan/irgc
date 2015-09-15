# Makefile

CPP = g++
CPPFLAGS= -DIL_STD -msse4 -O3 -s -ffast-math -funroll-loops -ansi -Wall -DTRWS_DEBUG_OUTPUT -DNDEBUG

all: irgc

irgc: *.cpp *.h
	$(CPP) $(CPPFLAGS) *.cpp -o irgc

clean:
	rm -f irgc
