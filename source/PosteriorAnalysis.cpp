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

#include "TGaxis.h"

#include "TApplication.h"

#include "PosteriorAnalysis.h"
#include "MyLib.h"

using namespace std;


int main(int argc, char** argv)
{
    cout << "\nThis program is part of the BaTMAN software package\n"
         << "(Bayesian-unfolding Toolkit for Multi-foil Activation with Neutrons)\n"
         << "Copyright (C) 2021 Davide Chiesa (University and INFN of Milano - Bicocca)\n"
         << "This program comes with ABSOLUTELY NO WARRANTY\n"<<endl;
    cout << "\tPress any key to START\n";
    cin.get();

	if (argc != 3)	{
		cout<<"\nInstructions to launch PosteriorAnalysis:"<<endl<<endl;
		cout<< "PosteriorAnalysis <Input_data.txt> <Energy_Groups.txt>" << endl << endl;
		return 1;
	}
	
	for (int i=1; i<argc; i++)	if (CheckFileExists(argv[i]) == 0) return 1;
			
	string InputData=argv[1];
	string En_group=argv[2];
	
	string OutputTextFile = "PosteriorAnalysis.txt";
	string OutputRootFile = "PosteriorAnalysis.root";
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	SetCanvasOptions();
	ofstream out(OutputTextFile.c_str());	
	TFile *file = new TFile(OutputRootFile.c_str(),"RECREATE");
	file->Close();
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	cout<<"Loading the data saved in JagsOutput.root"<<endl;
	
	TFile *f = TFile::Open("JagsOutput.root","READ");
	if (!f) { 
		cout << "Error! The file JagsOutput.root is not in this folder. \n"
		     << "       Run jags2root to produce it.\n";
		return 1;
	}
	
	TTree *t; 
	f->GetObject("Statistics",t);
	
	vector<double> *mean_prt = 0;
	vector<double> *RMS_prt = 0;
	vector<string> *VarList_prt = 0;
		
	t->SetBranchAddress("mean",&mean_prt);
	t->SetBranchAddress("RMS",&RMS_prt);
	t->SetBranchAddress("VarList",&VarList_prt);
	t->GetEntry(0);
	t->ResetBranchAddresses();
	
	vector<double> mean;
	vector<double> RMS;
	vector<string> VarList;
	
	for (int i=0; i < VarList_prt->size(); i++)
	{
		VarList.push_back (VarList_prt->at(i));
		mean.push_back (mean_prt->at(i));
		RMS.push_back (RMS_prt->at(i));		
		//cout<< VarList[i] << "\t" << mean[i] << "\t" << RMS[i] << endl;
	}
	
	TNtuple* NtuSum = (TNtuple*) f->Get("NtuSum");
	cout<<"OK!"<<endl;	
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	cout<<"Subdividing the monitored variables"<<endl;
	
	string mask = "phi";	
	string varNameCopy;
	
	vector <string> ListPHI; 
	vector<string> ListGroups; 
	vector<double> meanPHI;
	vector<double> RMSPHI;
	
	for (int i=0; i < VarList.size(); i++)	{
		varNameCopy = VarList[i];
		varNameCopy.resize(mask.length());		
		if (varNameCopy.compare(mask)==0) {
			ListPHI.push_back(VarList[i]);
			meanPHI.push_back(mean[i]);
			RMSPHI.push_back(RMS[i]);
		}
	}
	
	mask = "R";
	vector <string> ListRate; 
	vector<double> meanRate;
	vector<double> RMSRate;
	
	for (int i=0; i < VarList.size(); i++)	{
		varNameCopy = VarList[i];
		varNameCopy.resize(mask.length());		
		if (varNameCopy.compare(mask)==0) {
			ListRate.push_back(VarList[i]);
			meanRate.push_back(mean[i]);
			RMSRate.push_back(RMS[i]);
		}
	}
	
	mask = "random";
	vector <string> ListRandom; 
	vector<double> meanRandom;
	vector<double> RMSRandom;
	
	for (int i=0; i < VarList.size(); i++)	{
		varNameCopy = VarList[i];
		varNameCopy.resize(mask.length());		
		if (varNameCopy.compare(mask)==0) {
			ListRandom.push_back(VarList[i]);
			meanRandom.push_back(mean[i]);
			RMSRandom.push_back(RMS[i]);
		}
	}	
		
	cout<<"OK!"<<endl;	
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	vector<double> EGroups = FillEgroups (En_group);
	int Ngroups=ListPHI.size()-1;
	if (Ngroups != EGroups.size()-1) {
		cout << "Error: Ngroups != EGroups.size()-1" <<endl;
		return 1;
	}
	
	double FluxUnit;
	cout<<"Set order of magnitude of neutron flux used for unfolding (e.g. 1e12): ";
	cin>>FluxUnit;
	
	out << "Integral neutron flux per group (n/cm2/s)"<<endl;	
	out << left << setw(15) <<"E_low (MeV)"
		    << setw(15) <<"E_up (MeV)" 
		    << setw(15) << "Mean" 
		    << setw(15) << "StdDev" 
		    << setw(15) <<right<< "RelErr(\%)" << endl;

	for (int i=0; i<Ngroups; i++) {
		ListGroups.push_back(ListPHI[i]);
		out<< scientific <<left << setw(15) <<setprecision(4)<< EGroups[i] 
				<< setw(15) << EGroups[i+1] 
				<< setw(15) << meanPHI[i]*FluxUnit
				<< setw(15) << RMSPHI[i]*FluxUnit 
				<<setw(15) << fixed<<right<<setprecision(2) << 100*RMSPHI[i]/meanPHI[i]<<"\%" << endl;
	}
	out<< scientific <<left << setw(30) <<setprecision(3)<< "Phi_Tot" 
				<< setw(15) << meanPHI[Ngroups]*FluxUnit
				<< setw(15) << RMSPHI[Ngroups]*FluxUnit 
				<< setw(15) << fixed<<right<<setprecision(2) << 100*RMSPHI[Ngroups]/meanPHI[Ngroups]<<"\%" << endl;
	out<<endl;
	
	cout<<"\nANALYSIS OF ACTIVATION RATES\n";
	cout<<"Reading the file: "<<InputData<<endl;
	vector<string> ReactionName;
	vector<double> obsR;
	vector<double> sigmaR;
	vector<double> AM;
	vector<double> IA;
	ReadInput (InputData, ReactionName, obsR, sigmaR, AM, IA);
	
	if (ListRate.size() != ReactionName.size()) {
		cout << "Error: ListRate != ReactionName" <<endl;
		return 1;
	}
	
	out<<"\nANALYSIS OF ACTIVATION RATES [Bq/g]\n";
	out << left << setw(34) <<"Reaction"
	    << setw(15) << "ActRate" 
	    << setw(15) << "Err" 
	    << setw(15) << "ExpRate"
	    << setw(15) << "ExpErr"
	    << setw(10) << right << "Nsigma"
	    << setw(10) << "\% diff"<< endl;
	double N_avog=6.0221413e23;
	double nSigma, PercDiff;
	for (int i=0; i<ReactionName.size(); i++) {
		meanRate[i] /= (1e24*AM[i]/(FluxUnit*N_avog*IA[i]));
		RMSRate[i] /= (1e24*AM[i]/(FluxUnit*N_avog*IA[i]));
		nSigma = (meanRate[i]-obsR[i])/sqrt(pow(RMSRate[i],2)+pow(sigmaR[i],2));
		PercDiff = 100.*(meanRate[i]-obsR[i])/obsR[i];
		
		out 	<< setw(2) <<i+1<<setw(2) <<") "
			<< left<<setw(30) <<ReactionName[i]
			<< setw(15) <<scientific << setprecision(3) << meanRate[i]
			<< setw(15) << RMSRate[i]
			<< setw(15) << obsR[i]
			<< setw(15) << sigmaR[i]
			<< setw(10) << right << fixed << setprecision(2) << nSigma
			<< setw(10) << PercDiff<<"\%"<< endl;		
	}
	
	cout<<"OK!"<<endl;
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	TGaxis::SetMaxDigits(2);
	// MARGINAL DISTRIBUTIONS
	cout<< "\nDrawing Marginal PDF " << endl;
	MarginalPDF (OutputRootFile, VarList, NtuSum, out);
	
	TGaxis::SetMaxDigits();
	// CORRELATION PLOTS
	double **CorrCoeffMatrix;
	int nv=Ngroups+1;
	cout<< "\nDrawing correlation scatter plots" << endl;
	CorrCoeffMatrix = new double*[nv];
	for (int j=0; j<nv; j++){
		CorrCoeffMatrix[j] = new double [nv];
	}
	CorrelationPlots (OutputRootFile, ListGroups, NtuSum, CorrCoeffMatrix, out);
	
	
	double **CorrMatrix2;
	int Nreactions = ListRate.size();
	cout<< "\nAnalyzing correlations between PHI and ActRATE" << endl;
	CorrMatrix2 = new double*[Nreactions];
	for (int k=0; k<nv; k++){
		CorrMatrix2[k] = new double [Nreactions];
	}	
	
	CorrelRatePHI(OutputRootFile, ListGroups, ListRate,  NtuSum, CorrMatrix2, out);
	
	
	/*double **CorrMatrix3;
	cout<< "\nAnalyzing correlations between ActRATEs" << endl;
	CorrMatrix3 = new double*[Nreactions];
	for (int k=0; k<Nreactions; k++){
		CorrMatrix3[k] = new double [Nreactions];
	}	
	
	CorrelRates(OutputRootFile, ListRate,  NtuSum, CorrMatrix3, out);*/

	
	cout << "\nOutput files are: "<< OutputTextFile << " and " << OutputRootFile <<endl;
	
	cout << "Bye!\n" <<endl;

	return 0; 
}

