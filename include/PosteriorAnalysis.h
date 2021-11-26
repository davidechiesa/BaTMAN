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

#ifndef PosteriorAnalysis_h
#define PosteriorAnalysis_h

#include <string>
#include <vector>
#include <fstream>
#include <string>
#include "TNtuple.h"
#include "TH1F.h"

using namespace std;

void ReadInput (string filename, vector<string> &ReactionName, vector<double> &obsR, vector<double> &sigmaR, vector<double> &AtomMass, vector<double> &IsotAb);
vector<double> FillEgroups (string filename);

double fitfunz(double *x, double *p);
double fitfunzExpo(double *x, double *p);

void SetCanvasOptions ();

void MarginalPDF (string OutputRootFile, vector <string> &VarList, TNtuple* ntu, ofstream &out);
void DrawMarginalPDF(TFile* &file, vector <TH1F*> MargPDF, ofstream &out);
TH1F* HighlightedRegion (TH1F* Original, double prob, int& binlow, int& binup);
void CorrelationPlots (string OutputRootFile, vector <string> &VarList, TNtuple* ntu, double** CorrCoeffMatrix, ofstream &out);

void CorrelRatePHI(string OutputRootFile, vector <string> &PHIList, vector <string> &RateList,  TNtuple* ntu, double** CorrCoeffMatrix, ofstream &out);
void CorrelRates(string OutputRootFile, vector <string> &RateList,  TNtuple* ntu, double** CorrCoeffMatrix, ofstream &out);

#endif
