.PHONY: all clean O2 O3

all: O2 O3

O2: *.cpp
	g++ -O2 -o trjreadO2 -std=c++11 *.cpp

O3: *.cpp
	g++ -O3 -o trjreadO3 -std=c++11 *.cpp

clean:
	rm trjreadO2 trjreadO3
