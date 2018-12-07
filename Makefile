metacomm_evo : metacomm_evo.cpp include/classes.h include/procedures.h
	g++ -Wall -g -o metacomm_evo metacomm_evo.cpp -lgsl -lgslcblas -lm
	
all:	metacomm_evo

clean:
	rm -f metacomm_evo
