# Makefile

CPP = g++
CPPFLAGS= -DIL_STD -msse4 -O3 -s -ffast-math -funroll-loops -ansi -Wall -DTRWS_DEBUG_OUTPUT -DNDEBUG

all: irgc

irgc: maxflow/*.cpp maxflow/*.h alpha_exapnsion_west/*.cpp alpha_exapnsion_west/*.h *.cpp *.h
	$(CPP) $(CPPFLAGS)  maxflow/*.cpp alpha_exapnsion_west/*.cpp *.cpp -o irgc

clean:
	rm -f irgc
