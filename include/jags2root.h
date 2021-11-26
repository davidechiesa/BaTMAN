/*===============================================================================
 * Name of the project: BaTMAN 
 *     (Bayesian-unfolding Toolkit for Multi-foil Activation with Neutrons)
 *
 * Copyright (C) 2021 
 *
 * Author: Davide Chiesa    (University and INFN of Milano - Bicocca)
 *
 * Address: Piazza della Scienza 3, Edificio U2, Milano, Italy 20126
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ================================================================================*/
#ifndef jags2root_h
#define jags2root_h

#include <string>
#include <vector>
#include <fstream>
#include <string>
#include "TNtuple.h"
using namespace std;

double fitfunz(double *x, double *p);
double fitfunzExpo(double *x, double *p);

void SetCanvasOptions ();

void Jags2Root (string OutputRootFile, vector <string> &VarList, vector <string> &NtupleList, TNtuple* NtuSum);
void CalcNtupleStat (TNtuple *ntu, vector<string> &VarList, vector<double> &mean, vector<double> &RMS);
void PrintOutputSimple (vector<string> &VarList, vector<double> &mean, vector<double> &RMS, ofstream &out);

void TracePlots (string OutputRootFile, vector <string> &VarList, string NtuName);

#endif
