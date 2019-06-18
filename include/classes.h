//============================================================================
// Name        : Metacommunity evolution
// Author      : Emanuel A. Fronhofer
// Version     : v1
// Date	       : June 2019
//============================================================================

/*
	Copyright (C) 2019  Emanuel A. Fronhofer

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

//____________________________________________________________________________
//----------------------------------------------------------- define constants

const int WORLDDIM_X = 30;														// world dimensions, always quadratic
const int WORLDDIM_Y = 1;														// world dimensions, always quadratic
const int RS = 2;															// random seed
const int N_ALLELES = 2;
const int N_LOCI = 1;
const int NO_SPECIES = 80;

//____________________________________________________________________________
//------------------------------------------------------------- define classes

// one individual ------------------------------------------------------------
class TIndiv {
public:
	TIndiv();
	float dispRate[N_ALLELES][N_LOCI];
	float abiotTrait[N_ALLELES][N_LOCI];
	int preChangeLocation;
};

TIndiv::TIndiv() { //constructor for TIndiv
	for (int i = 0; i < N_ALLELES; ++i) {
		for (int act_loc = 0; act_loc < N_LOCI; ++act_loc) {
			dispRate[i][act_loc] = 0;
			abiotTrait[i][act_loc] = 0;
		}
	}
	preChangeLocation = 0;
}

// one population -----------------------------------------------------------------
class TPop {
public:
	TPop();
	vector<TIndiv> males;
	vector<TIndiv> newMales;
	vector<TIndiv> females;
	vector<TIndiv> newFemales;
};

TPop::TPop() {
	males.clear();
	newMales.clear();
	females.clear();
	newFemales.clear();
}

// one patch with a community-------------------------------------------------
class TPatch {
public:
	TPatch();
	TPop species[NO_SPECIES];
	float abiotCond;
	bool analog;
};

TPatch::TPatch() {
	for (int s = 0; s < NO_SPECIES; ++s) {
		species[s].males.clear();
		species[s].females.clear();
		species[s].newMales.clear();
		species[s].newFemales.clear();
	}
	abiotCond = 0;	
	analog = 0;
}
