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

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip> 
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "TROOT.h"
#include "TStyle.h"
#include "TGaxis.h"

#include "TCanvas.h"
#include "TLegend.h"

#include "TTree.h"
#include "TNtuple.h"

#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"

#include "TF1.h"
#include "TGraph.h"

#include "TFile.h"
#include "TDirectory.h"

#include "TApplication.h"

#include "jags2root.h"

using namespace std;


int main(int argc, char** argv)
{	
	cout << "\nThis program is part of the BaTMAN software package\n"
	     << "(Bayesian-unfolding Toolkit for Multi-foil Activation with Neutrons)\n"
	     << "Copyright (C) 2021 Davide Chiesa (University and INFN of Milano - Bicocca)\n"
	     << "This program comes with ABSOLUTELY NO WARRANTY\n"<<endl;
	     
	cout << "\tPress any key to START\n";
	cin.get();
	
	string OutputTextFile = "JagsOutput.txt";
	string OutputRootFile = "JagsOutput.root";
	
	ofstream out(OutputTextFile.c_str());	
	
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

	
	vector <string> VarList;
	vector <string> NtupleList;
	TNtuple* NtuSum;
	Jags2Root (OutputRootFile, VarList, NtupleList, NtuSum);	
	int nv=VarList.size();
	
	// TRACE PLOTS	
	SetCanvasOptions();	
	TracePlots (OutputRootFile, VarList, "NtuSum");
		
	// MEAN AND STANDARD DEVIATION CALCULATION
	TFile *file = new TFile(OutputRootFile.c_str(),"update");
	for (int i=0; i<NtupleList.size(); i++)
	{
		TNtuple *ntu = (TNtuple*) file->Get(NtupleList[i].c_str());		
		cout<<"\n"<< NtupleList[i] <<endl;
		cout<<"Calculating MEAN and STANDARD DEVIATION "<<endl;
		vector<double> mean (nv,0.);
		vector<double> RMS (nv,0.);
		CalcNtupleStat (ntu, VarList, mean, RMS);
	}
	
	cout<<"\nNtuSum"<< endl;	
	cout<<"Calculating MEAN and STANDARD DEVIATION "<<endl;
	vector<double> mean (nv,0.);
	vector<double> RMS (nv,0.);

	NtuSum = (TNtuple*) file->Get("NtuSum"); 
	CalcNtupleStat (NtuSum, VarList, mean, RMS);
	
	// Create a TTree
	TTree *t = new TTree("Statistics","Tree with vectors");
	t->Branch("VarList",&VarList);
	t->Branch("mean",&mean);
	t->Branch("RMS",&RMS);
	t->Fill();
	t->Write();
	file->Close();

	PrintOutputSimple (VarList, mean, RMS, out);	

	cout << "\nOutput files are: "<< OutputTextFile << " and " << OutputRootFile <<endl;
	cout << "Bye!\n" <<endl;
	
	return 0; 
}

