//============================================================================
// Name        : Metacommunity evolution
// Author      : Emanuel A. Fronhofer & Peter Kamal
// Version     : v2
// Date	       : November 2024
//============================================================================

/*
	Copyright (C) 2024  Emanuel A. Fronhofer & Peter Kamal

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
	
*/

//============================================================================

#include <iostream>
#include <cstdlib>								//standard C library
#include <ctime>								//access system time library
#include <fstream>								//file streaming library
#include <string>								//string library included
#include <sstream>								//string streaming for reading numbers from

#include <vector>
#include <cmath>								//standard math library

#include <gsl/gsl_rng.h>						//gsl random number generator
#include <gsl/gsl_randist.h>					//gsl random distributions
#include <gsl/gsl_statistics.h>					//gsl some statistical methods
#include <gsl/gsl_statistics_double.h> 			//gsl some double statistical methods
#include <gsl/gsl_sort_double.h> 				//gsl sort double arrays

#include <algorithm>

using namespace std;

#include "include/procedures.h"					//procedure simplifications
#include "include/classes.h"					//class definitions

//_____________________________________________________________________________
//------------------------------------------------------------ global variables
unsigned int sim_time;															// actual time in simulation
int max_runs;																	// no. of repeats
float mut_sd_disp;																// variance used for mutations
float mut_rate_disp;															// mutation rate
float mut_sd_abiot;																// variance used for mutations
float mut_rate_abiot;															// mutation rate
float lambda_null;																// fertility
float mu0;																		// dispersal mortality
float sigma;																	// environmental stochasticity
float epsilon;																	// random patch extinctions
bool nnd;																		// nearest neighbour dispersal (yes=true, no=false)

TPatch world[WORLDDIM_X][WORLDDIM_Y];											// simulated world

unsigned int metacommunity_size;												// relative metapopulation size
float occupancy;																// metapopulation occupancy
float rel_emigrants;															// relative number of emigrants

float sigma_env;																// with of local adaptation curve
float start_disp;																// dispersal value at simultion start

float mean_local_adaptation;													// globally tracking level of adaptation (from 1 best to 0 worst)

float community_matrix[NO_SPECIES][NO_SPECIES];									// community matrix (holds all pairwise competition coefficients)

unsigned int env_change_start;													// start time of environmental change
unsigned int env_change_end;													// end time of environmental change
float env_change_delta;															// magnitude of environmental change

float gamma_diversity;															// gamma diversity (total no. of species in metacommunity)
float alpha_diversity;															// alpha diversity (mean local species number)
float gamma_diversity_analog;													// gamma diversity in analog patches
float alpha_diversity_analog;													// alpha diversity in analog patches
float gamma_diversity_nonanalog;												// gamma diversity in non-analog patches
float alpha_diversity_nonanalog;												// alpha diversity in non-analog patches

//_____________________________________________________________________________
//------------------------------------------------------------------ procedures

//------------------------------------------------------------- read parameters
void readParameters(){
	ifstream parinfile("input/parameters.in");							//parameter input file
	string buffer;
	istringstream is;

	getline(parinfile,buffer); getline(parinfile,buffer);
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> sim_time;																		//simulation time
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> max_runs;																		//no. of repeats
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> mut_sd_disp;																	//variance used for mutations dispersal tait
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> mut_rate_disp;																//mutation rate dispersal trait
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> start_disp;																	//dispersal rate at simulatipon start
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> mut_sd_abiot;																	//variance used for mutations abiotic trait
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> mut_rate_abiot;																//mutation rate abiotic trait
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> lambda_null;																	//fertility
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> sigma;																		//uncorrelated evironmental stochasticity
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> epsilon;																		//random patch extinction probability
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> mu0;																			//dispersal mortality
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> sigma_env;																	//with of local adaptation curve
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> env_change_start;																//start of environmental change
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> env_change_delta;																//magnitude of environmental change
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> env_change_end;																//end of environmental change
	parinfile.close();
}

//------------------------------------------------------------- read community matrix from file
void readCommunityMatrix(unsigned int actrun){
	stringstream community_matrix_path_stream;
	community_matrix_path_stream << "input/community_matrices/community_matrix_" << actrun + 1 << "_competition_symmetric.txt";
	string community_matrix_path = community_matrix_path_stream.str();
	ifstream matinfile(community_matrix_path.c_str());

	for (int i = 0; i < NO_SPECIES; ++i) {
		for (int j = 0; j < NO_SPECIES; ++j) {
			matinfile >> community_matrix[i][j]; // this is [row][column] or [to][from]
		}
	}

	matinfile.close();
}

//------------------------------------------------------- initialize simulation
void Initialize(){
	
	float init_abiot_trait[NO_SPECIES];
	
	// get initial values for all species
	for (unsigned int s = 0; s < NO_SPECIES; ++s){
		init_abiot_trait[s] = ran()*(float(WORLDDIM_X)/2-1);
	}
	
	// initialize patches and individuals in patches
	for (unsigned int x = 0; x < WORLDDIM_X; ++x) {
			for (unsigned int y = 0; y < WORLDDIM_Y; ++y) {
			
			// initialize abiotic conditions
			if (float(x) < float(WORLDDIM_X)/2) {
				world[x][y].abiotCond = float(x);
			} else {
				world[x][y].abiotCond = float(WORLDDIM_X) - float(x) - 1;
			}
			
			if((world[x][y].abiotCond + (env_change_end - env_change_start)*env_change_delta) < float(WORLDDIM_X)/2){
				world[x][y].analog = 1;
			}else{
				world[x][y].analog = 0;
			}
			
			//cout << x << "    " << world[x].abiotCond << "    " << "    " << world[x].analog <<endl;
				
			for (unsigned int s = 0; s < NO_SPECIES; ++s){

				// clear the world
				world[x][y].species[s].females.clear();
				world[x][y].species[s].males.clear();
				world[x][y].species[s].newFemales.clear();
				world[x][y].species[s].newMales.clear();

				// initialize individuals in this patch
				// females
				for (int f = 0; f < round((lambda_null-1)/community_matrix[s][s]/2); ++f) {
					TIndiv newfemale;
					for (int a = 0; a < N_ALLELES; ++a) {
						for(int act_loc = 0; act_loc < N_LOCI; ++act_loc) {

							// initialize dispersal randomly, i.e. with max standing genetic variation if evolution is modelled
							if (mut_sd_disp > 0 && mut_rate_disp > 0) {
								newfemale.dispRate[a][act_loc] = ran()/(float(N_LOCI)*float(N_ALLELES));
							} else {
								// otherwise set dispersal rate to fixed value
								newfemale.dispRate[a][act_loc] = start_disp/(float(N_LOCI)*float(N_ALLELES));
							}

							//newfemale.abiotTrait[a] = ran()*(float(WORLDDIM)/2-1);
							newfemale.abiotTrait[a][act_loc] = init_abiot_trait[s]/(float(N_LOCI)*float(N_ALLELES));
						}
					}
					
					// init pre change location with 0
					newfemale.preChangeLocation = 0;
					
					world[x][y].species[s].females.push_back(newfemale);
				}

				// males
				for (int m = 0; m < round((lambda_null-1)/community_matrix[s][s]/2); ++m) {
					TIndiv newmale;
					for (int a = 0; a < N_ALLELES; ++a) {
						for(int act_loc = 0; act_loc < N_LOCI; ++act_loc) {
							// initialize dispersal randomly, i.e. with max standing genetic variation if evolution is modelled
							if (mut_sd_disp > 0 && mut_rate_disp > 0) {
								newmale.dispRate[a][act_loc] = ran()/(float(N_LOCI)*float(N_ALLELES));
							} else {
								// otherwise set dispersal rate to fixed value
								newmale.dispRate[a][act_loc] = start_disp/(float(N_LOCI)*float(N_ALLELES));
							}

							//newmale.abiotTrait[a] = ran()*float(WORLDDIM)/2;
							newmale.abiotTrait[a][act_loc] = init_abiot_trait[s]/(float(N_LOCI)*float(N_ALLELES));
						}
					}
					
					// init pre change location with 0
					newmale.preChangeLocation = 0;
					
					world[x][y].species[s].males.push_back(newmale);
				}
			}
		}
	}
}
// ------------------------------------------------ analyze population dynamics
void Analyze(unsigned int acttime, int actrun){
	//reset metacommunity size and occupancy
	metacommunity_size = 0;
	occupancy = 0;

	unsigned int numberoccupied = 0;

	// help array to calculate gamma diversity (total number of species in metacommunity)
	float help_array_species[NO_SPECIES];
	float help_array_species_analog[NO_SPECIES];
	float help_array_species_nonanalog[NO_SPECIES];
	for (int s = 0; s < NO_SPECIES; ++s) {
		help_array_species[s] = 0;
		help_array_species_analog[s] = 0;
		help_array_species_nonanalog[s] = 0;
	}

	//help array to calculate mean alpha diversity (mean number of species per patch in the metacommunity)
	float help_array_alphadiv[WORLDDIM_X*WORLDDIM_Y];
	int alpha_cnt = 0;
	for (int x = 0; x < WORLDDIM_X; ++x) {
		for (int y = 0; y < WORLDDIM_Y; ++y) {
			help_array_alphadiv[alpha_cnt] = 0;
			alpha_cnt = alpha_cnt + 1;
		}
	}
	
	alpha_cnt = 0;
	for (int x = 0; x < WORLDDIM_X; ++x) {
		for (int y = 0; y < WORLDDIM_Y; ++y) {
			int help_act_occ = 0;
			for (int s = 0; s < NO_SPECIES; ++s) {
				unsigned int localpopsize = world[x][y].species[s].females.size() + world[x][y].species[s].males.size();
				metacommunity_size += localpopsize;
				if (localpopsize > 0) {
					help_act_occ = 1;
					// save that there is this species
					help_array_species[s] = 1;
					
					// get analog and non-analog counts of species
					if(world[x][y].analog == 1){
						help_array_species_analog[s] = 1;
					}else{
						help_array_species_nonanalog[s] = 1;
					}
					
					// save that this species is there for alpha diversity
					help_array_alphadiv[alpha_cnt] = help_array_alphadiv[alpha_cnt] + 1;
				}
			}

			if (help_act_occ == 1) {
				++numberoccupied;
			}
			alpha_cnt = alpha_cnt + 1;
		}
	}
	// calculate occupancy
	occupancy = float(numberoccupied) / (float(WORLDDIM_X)*float(WORLDDIM_Y));
	// calculate gamma diversity
	gamma_diversity = 0;
	gamma_diversity_analog = 0;
	gamma_diversity_nonanalog = 0;
	for (int s = 0; s < NO_SPECIES; ++s) {
		gamma_diversity = gamma_diversity + help_array_species[s];
		gamma_diversity_analog = gamma_diversity_analog + help_array_species_analog[s];
		gamma_diversity_nonanalog = gamma_diversity_nonanalog + help_array_species_nonanalog[s];
	}
	
	// total alpha diversity
	alpha_diversity = mean(help_array_alphadiv, WORLDDIM_X*WORLDDIM_Y);
	
	// get analog and non-analog alpha diversity
	float cnt_analog = 0;
	float cnt_nonanalog = 0;
	float help_alpha_analog = 0;
	float help_alpha_nonanalog = 0;
	alpha_cnt = 0;
	for (int x = 0; x < WORLDDIM_X; ++x) {
		for (int y = 0; y < WORLDDIM_Y; ++y) {
			if(world[x][y].analog == 1){
				cnt_analog = cnt_analog + 1;
				help_alpha_analog = help_alpha_analog + help_array_alphadiv[alpha_cnt];
			}else{
				cnt_nonanalog = cnt_nonanalog + 1;
				help_alpha_nonanalog = help_alpha_nonanalog + help_array_alphadiv[alpha_cnt];
			}
			alpha_cnt = alpha_cnt + 1;
		}
	}
	alpha_diversity_analog = help_alpha_analog / cnt_analog;
	alpha_diversity_nonanalog = help_alpha_nonanalog / cnt_nonanalog;
	
}

// ---------------------------------------------------- save individual results
void saveResults(int actrun, int acttime){
	// output file: individuals
	stringstream outputindiv_path_stream;
	outputindiv_path_stream << "output/output_individuals_run" << actrun << "_t" << acttime << ".out";
	string outputindiv_path = outputindiv_path_stream.str();
	ofstream outputindiv(outputindiv_path.c_str());

	// headers
	outputindiv << "x" << "    " << "y" << "    " << "species" << "    " << "sex"; 
	
	for (int a = 0; a < N_ALLELES; ++a) {
		for(int act_loc = 0; act_loc < N_LOCI; ++act_loc) {
			outputindiv << "    " << "dispRate_a" << a <<"l"<<act_loc << "    "  << "abiotTrait_a" << a <<"l"<<act_loc;
		}
	}
	
	outputindiv << "    " << "preChangePos" << endl;

	for (int x = 0; x < WORLDDIM_X; ++x) {
		for (int y = 0; y < WORLDDIM_Y; ++y) {
			for (int s = 0; s < NO_SPECIES; ++s){
				for (unsigned int f = 0; f < world[x][y].species[s].females.size(); ++f) {
					// write metapop results to file
					outputindiv << x << "    " << y << "    " << s << "    " << "f" << "    ";
					
					for (int a = 0; a < N_ALLELES; ++a) {
						for(int act_loc = 0; act_loc < N_LOCI; ++act_loc) {
							outputindiv << world[x][y].species[s].females.at(f).dispRate[a][act_loc] << "    " << world[x][y].species[s].females.at(f).abiotTrait[a][act_loc] << "    " ;
						}
					}
					
					 outputindiv << world[x][y].species[s].females.at(f).preChangeLocation << endl;
				}
				for (unsigned int m = 0; m < world[x][y].species[s].males.size(); ++m) {
					// write metapop results to file
					outputindiv << x << "    " << y << "    " << s << "    " << "m" << "    ";
					for (int a = 0; a < N_ALLELES; ++a) {
						for(int act_loc = 0; act_loc < N_LOCI; ++act_loc) {
							outputindiv << world[x][y].species[s].males.at(m).dispRate[a][act_loc] << "    " << world[x][y].species[s].males.at(m).abiotTrait[a][act_loc] << "    " ;
						}
					}
					 outputindiv << world[x][y].species[s].males.at(m).preChangeLocation << endl;
				}
			}
		}
	}
	// close indiv output file
	outputindiv.close();
}

// --------------------------------------------- find new patch during dispersal
vector<int> findNewPatch(int x, int y){

	vector<int> res;
	int dir;

	// nearest neighbour dispersal (nnd8)
	dir = floor(ran()*8);

	switch (dir) {
	case 0:
		res.push_back(x+1);
		res.push_back(y+1);
		break;
	case 1:
		res.push_back(x+1);
		res.push_back(y);
		break;
	case 2:
		res.push_back(x+1);
		res.push_back(y-1);
		break;
	case 3:
		res.push_back(x);
		res.push_back(y-1);
		break;
	case 4:
		res.push_back(x-1);
		res.push_back(y-1);
		break;
	case 5:
		res.push_back(x-1);
		res.push_back(y);
		break;
	case 6:
		res.push_back(x-1);
		res.push_back(y+1);
		break;
	case 7:
		res.push_back(x);
		res.push_back(y+1);
	break;
	default:
		cout << "Error in NND" << endl;
		break;
	}

	if (res.at(0) == WORLDDIM_X){
		res.at(0) = 0;
	}
	if (res.at(0) == -1){
		res.at(0) = WORLDDIM_X-1;
	}
	if (res.at(1) == WORLDDIM_Y){
		res.at(1) = 0;
	}
	if (res.at(1) == -1){
		res.at(1) = WORLDDIM_Y-1;
	}
	return(res);

}

// -------------------------------------------------------- phenotype calculation
// calculate phenotypes from 2D array
float calc_phenotype(float data[N_ALLELES][N_LOCI]){
	float allele_sum = 0;
	
	for (int a = 0; a < N_ALLELES; ++a) {
		for(int act_loc = 0; act_loc < N_LOCI; ++act_loc) {
			allele_sum = allele_sum + data[a][act_loc];
			}
		}
	return(allele_sum);
}

// -------------------------------------------------------- dispersal procedure
void Dispersal(){

	unsigned int no_emigrants = 0;
	unsigned int metapopsize = 0;
	rel_emigrants = 0;

	for (int x = 0; x < WORLDDIM_X; ++x) {
		for (int y = 0; y < WORLDDIM_Y; ++y) {
			for (int s = 0; s < NO_SPECIES; ++s) {
				// counter for metapopsize
				metapopsize += world[x][y].species[s].females.size() + world[x][y].species[s].males.size();

				// start with females
				for (unsigned int f = 0; f < world[x][y].species[s].females.size(); ++f) {
					// should the individual disperse?
					if (ran() < calc_phenotype(world[x][y].species[s].females.at(f).dispRate)){
						// this individual will disperse
						// increase counter
						++no_emigrants;
						// check whether this emigrant survives the dispersal process
						if (ran() > mu0){
							// find new patch (global dispersal)
							vector<int> newPatch = findNewPatch(x, y);
							// copy disperser into new patch
							TIndiv Disperser = world[x][y].species[s].females.at(f);
							world[newPatch.at(0)][newPatch.at(1)].species[s].newFemales.push_back(Disperser);
						}
						// delete emigrant from natal patch
						world[x][y].species[s].females.at(f) = world[x][y].species[s].females.back();
						world[x][y].species[s].females.pop_back();
						// decrease female loop counter
						--f;
					}
				}
				// continue with males
				for (unsigned int m = 0; m < world[x][y].species[s].males.size(); ++m) {
					// should the individual disperse?
					if (ran() < calc_phenotype(world[x][y].species[s].males.at(m).dispRate)){
						// this individual will disperse
						// increase counter
						++no_emigrants;
						// check whether this emigrant survives the dispersal process
						if (ran() > mu0){
							// find new patch (global dispersal)
							vector<int> newPatch = findNewPatch(x, y);
							// copy disperser into new patch
							TIndiv Disperser = world[x][y].species[s].males.at(m);
							world[newPatch.at(0)][newPatch.at(1)].species[s].newMales.push_back(Disperser);
						}
						// delete emigrant from natal patch
						world[x][y].species[s].males.at(m) = world[x][y].species[s].males.back();
						world[x][y].species[s].males.pop_back();
						// decrease female loop counter
						--m;
					}
				}
			}
		}
	}

	// now that dispersal is over, merge philopatrics and residents
	for (int x = 0; x < WORLDDIM_X; ++x) {
		for (int y = 0; y < WORLDDIM_Y; ++y) {
			for (int s = 0; s < NO_SPECIES; ++s) {
				// first copy the females
				for (unsigned int f = 0; f < world[x][y].species[s].newFemales.size(); ++f) {
					world[x][y].species[s].females.push_back(world[x][y].species[s].newFemales.at(f));
				}
				// erase the "old" immigrants from newFemales
				world[x][y].species[s].newFemales.clear();
				// then copy the males
				for (unsigned int m = 0; m < world[x][y].species[s].newMales.size(); ++m) {
					world[x][y].species[s].males.push_back(world[x][y].species[s].newMales.at(m));
				}
				// erase the "old" immigrants from newFemales
				world[x][y].species[s].newMales.clear();
			}		
		}
	}
	
	rel_emigrants = float(no_emigrants) / float(metapopsize);
}

// ------------------------------------------------------------------ mutations for dispersal rate
float mutate_disp(float allele){
	if(ran()< mut_rate_disp){
		float newallele = allele + gauss(mut_sd_disp/(float(N_LOCI)*float(N_ALLELES)));
		return(newallele);
	} else {
		return(allele);
	}
}

// ------------------------------------------------------------------ mutations for abiotic trait
float mutate_abiot(float allele){
	if(ran()< mut_rate_abiot){
		float newallele = allele + gauss(mut_sd_abiot/(float(N_LOCI)*float(N_ALLELES)));
		return(newallele);
	} else {
		return(allele);
	}
}


// ------------------------------------------------------------ larval survival
float larvalSurvival(unsigned int x, unsigned int y, unsigned int s){
	// following the Beverton-Holt model

	// total competition for species s can be calculated once per species here
	// reset total competition
	float total_competition = 0;

	// loop over all species to get total competition
	for (int s1 = 0; s1 < NO_SPECIES; ++s1) {
		total_competition = total_competition + ((world[x][y].species[s1].males.size() + world[x][y].species[s1].females.size()) * community_matrix[s][s1]);
	}

	// this version assumes a neutral community matrix
	float survival = 1/(1+total_competition);				//survival-probability of newborns based on theory of logistic growth

	return(survival);
}

float localAdaptation (unsigned int x, unsigned int y, unsigned int s, unsigned int f){
	// local adaptation: reduction of lambda_null due to mismatch between mother trait mean and environment
	// values range between 0 and 1
	// sigma is the width of the distribution

	float trait_mean = calc_phenotype(world[x][y].species[s].females.at(f).abiotTrait);
	float local_env = world[x][y].abiotCond;

	float adapt_fact = exp(-pow((local_env - trait_mean) / (2* sigma_env),double(2)));

	return(adapt_fact);

}

// --------------------------------------------------------------- reproduction
void Reproduction(){

	float sum_local_adaptation = 0;
	int sum_females = 0;

	for (int x = 0; x < WORLDDIM_X; ++x) {
		for (int y = 0; y < WORLDDIM_Y; ++y) {
			for (int s = 0; s < NO_SPECIES; ++s) {
				// just to be sure: resize new females and males vectors
				world[x][y].species[s].newFemales.clear();
				world[x][y].species[s].newMales.clear();

				// counter for all females
				sum_females = sum_females + world[x][y].species[s].females.size();


				// for each patch check whether there are females and males
				if (world[x][y].species[s].females.size() > 0 && world[x][y].species[s].males.size() > 0) {

					// calculate larval survival for all individuals of species s in patch x
					float survival = larvalSurvival(x, y, s);

					// local lambda
					float lambda_local = lognorm(lambda_null, sigma);

					// females choose their mates
					for (unsigned int f = 0; f < world[x][y].species[s].females.size(); ++f) {

						// randomly choose male
						unsigned int m = floor(ran()*world[x][y].species[s].males.size());
						// calculate local adaptation for current female
						float loc_adapt = localAdaptation(x,y,s,f);
						sum_local_adaptation = sum_local_adaptation + loc_adapt;
						// calculate number of offspring from this mating
						int no_offspring = poisson(2*lambda_local*loc_adapt*survival);
						// loop over offspring
						for (int o = 0; o < no_offspring; ++o) {
							// sex ratio is 0.5
							if (ran() < 0.5) {
								// females
								// initialize new individual
								TIndiv newOffspring;
								
								int allele_dm = floor(ran()*N_ALLELES);
								int allele_af = floor(ran()*N_ALLELES);
								int allele_am = floor(ran()*N_ALLELES);
								int allele_df = floor(ran()*N_ALLELES);
								for(int act_loc = 0; act_loc < N_LOCI; ++act_loc) {
									// inherit emigration rate alleles from mother
									newOffspring.dispRate[0][act_loc] = mutate_disp(world[x][y].species[s].females.at(f).dispRate[allele_dm][act_loc]);
									// inherit emigration rate allele from father
									newOffspring.dispRate[1][act_loc] = mutate_disp(world[x][y].species[s].males.at(m).dispRate[allele_df][act_loc]);
									// inherit abiotic trait allele from mother
									newOffspring.abiotTrait[0][act_loc] = mutate_abiot(world[x][y].species[s].females.at(f).abiotTrait[allele_am][act_loc]);
									// inherit abiotic trait allele from father
									newOffspring.abiotTrait[1][act_loc] = mutate_abiot(world[x][y].species[s].males.at(m).abiotTrait[allele_af][act_loc]);
								}
							
								// inherit pre change location marker maternally
								newOffspring.preChangeLocation = world[x][y].species[s].females.at(f).preChangeLocation;

								// add new individual to new females vector
								world[x][y].species[s].newFemales.push_back(newOffspring);

							} else {
								//males
								// initialize new individual
								TIndiv newOffspring;
								int allele_dm = floor(ran()*N_ALLELES);
								int allele_af = floor(ran()*N_ALLELES);
								int allele_am = floor(ran()*N_ALLELES);
								int allele_df = floor(ran()*N_ALLELES);
								for(int act_loc = 0; act_loc < N_LOCI; ++act_loc) {
									// inherit emigration rate allele from mother
									newOffspring.dispRate[0][act_loc] = mutate_disp(world[x][y].species[s].females.at(f).dispRate[allele_dm][act_loc]);
									// inherit dispersal rate allele from father
									newOffspring.dispRate[1][act_loc] = mutate_disp(world[x][y].species[s].males.at(m).dispRate[allele_df][act_loc]);
									// inherit abiotic trait allele from mother
									newOffspring.abiotTrait[0][act_loc] = mutate_abiot(world[x][y].species[s].females.at(f).abiotTrait[allele_am][act_loc]);
									// inherit abiotic trait allele from father
									newOffspring.abiotTrait[1][act_loc] = mutate_abiot(world[x][y].species[s].males.at(m).abiotTrait[allele_af][act_loc]);
								}

								// inherit pre change location marker maternally
								newOffspring.preChangeLocation = world[x][y].species[s].females.at(f).preChangeLocation;

								// add new individual to new males vector
								world[x][y].species[s].newMales.push_back(newOffspring);
							}
						}
					}
				}
			}
		}
	}
	// calculate mean local adaptation
	mean_local_adaptation = sum_local_adaptation / float(sum_females);

}



// -------------------------------------------------- death of annual organisms
void Death(){
	for (int x = 0; x < WORLDDIM_X; ++x) {
		for (int y = 0; y < WORLDDIM_Y; ++y) {
			// include local patch extinctions
			if (ran() > epsilon){
				// no local patch extinction
				for (int s = 0; s < NO_SPECIES; ++s) {
					int local_offspring_no = world[x][y].species[s].newFemales.size()+world[x][y].species[s].newMales.size();
					//cout << local_offspring_no << endl;
					if (local_offspring_no > 0) {
						// now clear adult vectors
						world[x][y].species[s].males.clear();
						world[x][y].species[s].females.clear();
						// now copy new females into females
						for (unsigned int nf = 0; nf < world[x][y].species[s].newFemales.size(); ++nf) {
							world[x][y].species[s].females.push_back(world[x][y].species[s].newFemales.at(nf));
						}
						// now copy surviving new males into males
						for (unsigned int nm = 0; nm < world[x][y].species[s].newMales.size(); ++nm) {
							world[x][y].species[s].males.push_back(world[x][y].species[s].newMales.at(nm));
					}
					// clear new females vector
					world[x][y].species[s].newFemales.clear();
					// clear new males vector
					world[x][y].species[s].newMales.clear();
				} else {
					world[x][y].species[s].males.clear();
					world[x][y].species[s].females.clear();
					}
				}
			}else{
				// local patch extinction: empty all vectors
				for (int s = 0; s < NO_SPECIES; ++s) {
					world[x][y].species[s].newFemales.clear();
					world[x][y].species[s].newMales.clear();
					world[x][y].species[s].males.clear();
					world[x][y].species[s].females.clear();
				}
			}
		}
	}
}

// -------------------------------------------------- implement environmental change
void EnvironmentalChange(){
	for (int x = 0; x < WORLDDIM_X; ++x) {
		for (int y = 0; y < WORLDDIM_Y; ++y) {
			world[x][y].abiotCond = world[x][y].abiotCond + env_change_delta;
		}
	}
}


// -------------------------------------------------- calculate mean trait of one species in one patch at one time
float meantrait(int s, int x, int y){
	
	float abiot_tait_sum = 0;
	for (unsigned int m = 0; m < world[x][y].species[s].males.size(); ++m){
		abiot_tait_sum = abiot_tait_sum + calc_phenotype(world[x][y].species[s].males.at(m).abiotTrait);
	}
	for (unsigned int f = 0; f < world[x][y].species[s].females.size(); ++f){
		abiot_tait_sum = abiot_tait_sum + calc_phenotype(world[x][y].species[s].females.at(f).abiotTrait);
	}

	float meantrait = abiot_tait_sum / float((world[x][y].species[s].females.size() + world[x][y].species[s].males.size()));

	return meantrait;
}


//_____________________________________________________________________________
//------------------------------------------------------------------------ main

int main() {
	// random number generator
	//specify_rng(time(NULL));
	specify_rng(RS);

	//read parameters for all simulation runs
	readParameters();

	// output file: sim
	stringstream outputsim_path_stream;
	outputsim_path_stream << "output/output_sim.out";
	string outputsim_path = outputsim_path_stream.str();
	ofstream outputsim(outputsim_path.c_str());

	// outputfile header
	outputsim << "patch_x"  << "    " << "patch_y"  << "    " << "species" << "    " << "popsize"<< "    " << "meantrait"<< "    " << "time"<<   "    " << "replicate" << "    " << "emirate" << "    " << "mutrate" <<  endl;

	// repeat loop
	for (int actrun = 0; actrun < max_runs; ++actrun) {

		// output file: metapop
		stringstream outputmetapop_path_stream;
		outputmetapop_path_stream << "output/output_metapop_run" << actrun << ".out";
		string outputmetapop_path = outputmetapop_path_stream.str();
		ofstream outputmetapop(outputmetapop_path.c_str());

		// outputfile header
		outputmetapop << "time"  << "    " << "occupancy" << "    " << "metacommunity_size"<< "    " << "alpha_diversity"<< "    " << "alpha_diversity_analog" << "    " << "alpha_diversity_nonanalog"<< "    " << "gamma_diversity"<<  "    " << "gamma_diversity_analog" << "    " << "gamma_diversity_nonanalog"<< "    " << "emirate" << "    " << "mean_local_adaptation" <<  endl;

		// read the community matrix
		readCommunityMatrix(actrun);

		// initialize
		Initialize();

		// time loop
		for (unsigned int acttime = 0; acttime < sim_time; ++acttime) {

			// environmental change happens after a given time
			if (acttime > env_change_start && acttime < env_change_end) {
				EnvironmentalChange();
			}

			// natal dispersal
			Dispersal();
			
			// reproduction
			Reproduction();

			// density regulation and death of adults
			Death();

			// analyze metapopulation
			Analyze(acttime, actrun);

			// write metapop results to file
			outputmetapop << acttime  << "    " << occupancy << "    " << metacommunity_size <<  "    " << alpha_diversity <<"    " << alpha_diversity_analog <<"    " << alpha_diversity_nonanalog <<"    " << gamma_diversity <<"    " << gamma_diversity_analog <<"    " << gamma_diversity_nonanalog << "    " << rel_emigrants << "    " << mean_local_adaptation << endl;
			// write run and time to console
			//cout << actrun << "   " << acttime  << "    " << occupancy << "    " << metacommunity_size << "    " << alpha_diversity <<"    " << alpha_diversity_analog <<"    " << alpha_diversity_nonanalog <<"    " << gamma_diversity <<"    " << gamma_diversity_analog <<"    " << gamma_diversity_nonanalog <<  "   " << rel_emigrants << "    " << mean_local_adaptation << endl;

			// save individual results before environmental change starts
			if (acttime == env_change_start-1){
				
				// generate simulation level output
				for (int x = 0; x < WORLDDIM_X; ++x) {
					for (int y = 0; y < WORLDDIM_Y; ++y) {
						for (int s = 0; s < NO_SPECIES; ++s) {						
							// output to file
							outputsim << x  << "    " << y  << "    " << s << "    " <<  world[x][y].species[s].females.size() + world[x][y].species[s].males.size() << "    " <<  meantrait(s, x, y)<< "    " << acttime<<   "    " << actrun << "    " << start_disp << "    " << mut_sd_abiot <<  endl;

							// assign marker for pre change position in landscape
							for (unsigned int f = 0; f < world[x][y].species[s].females.size(); ++f) {
								world[x][y].species[s].females.at(f).preChangeLocation = x;
							}
							for (unsigned int m = 0; m < world[x][y].species[s].males.size(); ++m) {
								world[x][y].species[s].males.at(m).preChangeLocation = x;
							}		
						}
					}
				}
				
				saveResults(actrun, acttime);
				
			}			
			// save individual results after environmental change ends
			if (acttime == env_change_end-1){
				saveResults(actrun, acttime);
				
				// gerate simulation level output
				for (int x = 0; x < WORLDDIM_X; ++x) {
					for (int y = 0; y < WORLDDIM_Y; ++y) {
						for (int s = 0; s < NO_SPECIES; ++s) {
							// output to file
							outputsim << x  << "    " << y  << "    " << s << "    " <<  world[x][y].species[s].females.size() + world[x][y].species[s].males.size() << "    " <<  meantrait(s, x, y)<< "    " << acttime<<   "    " << actrun << "    " << start_disp << "    " << mut_sd_abiot <<  endl;					
						}
					}
				}
			}
			// save individual results at simulation end
			if (acttime == sim_time-1){
				saveResults(actrun, acttime);
				
				// gerate simulation level output
				for (int x = 0; x < WORLDDIM_X; ++x) {
					for (int y = 0; y < WORLDDIM_Y; ++y) {
						for (int s = 0; s < NO_SPECIES; ++s) {
							// output to file
							outputsim << x  << "    " << y  << "    " << s << "    " <<  world[x][y].species[s].females.size() + world[x][y].species[s].males.size() << "    " <<  meantrait(s, x, y)<< "    " << acttime<<   "    " << actrun << "    " << start_disp << "    " << mut_sd_abiot <<  endl;					
						}
					}
				}
			}
		}
		// close metapop output file
		outputmetapop.close();
	}
	
	// close sim output file
	outputsim.close();
		
	return 0;
}
