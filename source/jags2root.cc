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


#ifndef jags2root_cc
#define jags2root_cc

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

#include "TApplication.h"
#include "TGaxis.h"
#include "TColor.h"
#include "TCut.h"

#include "jags2root.h"
#include "MyLib.h"

using namespace std;

double fitfunz(double *x, double *p)  {
	double PI=4*atan(1);
	double arg=(x[0]-p[0])/p[1];
	return (exp(-0.5*arg*arg)/sqrt(2*PI*p[1]*p[1]));
}

double fitfunzExpo(double *x, double *p)  { 
	 return p[0]*exp(-p[0]*x[0]);
}


// CANVAS OPTIONS
int n_lines;
int n_columns;
int n_pad;	
int CanvasPx;
int CanvasPy;
int CanvasNx;
int CanvasNy;	
int SleepSeconds;

const int Ncolors = 22;
//Int_t MyPalette[Ncolors];
Color_t ColorList[] = {
	kRed,
	kGreen,
	kBlue,
	kMagenta,
	kOrange,
	kCyan,
	kRed+2,
	kGreen+3,
	kBlue-7,
	kMagenta+3,
	kOrange-7,
	kCyan+3,
	kYellow+2,
	kGray,
	kRed-9,
	kGreen-10,
	kAzure+7,
	kMagenta-8,
	kOrange+2,
	kCyan-8,
	kYellow-7,
	kGray+2
};
	
void SetCanvasOptions ()
{
	gROOT->SetBatch();
	gROOT->SetStyle("Plain");
	gStyle->SetOptFit(0);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetOptStat(kFALSE); 
	//TGaxis::SetMaxDigits(2);
	
	n_lines=2;
	n_columns=3;
	n_pad=n_lines*n_columns;
	
	CanvasPx = 1200 - (n_columns*250);
	CanvasPy = 0;
	CanvasNx = n_columns*250;
	CanvasNy = n_lines*350;
	
	SleepSeconds = 0;
	
	return;	
}


void Jags2Root (string OutputRootFile, vector <string> &VarList, vector <string> &NtupleList, TNtuple* NtuSum)
{
	TFile *file = new TFile(OutputRootFile.c_str(),"RECREATE");
	
	cout<<"Reading index file: CODAindex.txt..."<<endl;
	int nv=0;
	vector <int> varIniz;
	vector <int> varFine;
	int varIniz1, varFine1;
	string varName1;  

	ifstream in("CODAindex.txt"); 
	if (in.good()==0)
	{
		cout<<"ERROR: The JAGS output file CODAindex.txt is not here. "<<endl;
		cin.get();
	}
	while(1)
	{
		in >> varName1 >> varIniz1>> varFine1;
		if (in.eof()==1) break;
		
		// erasing square brackets from variable names
		for (int i=0; i<varName1.size();i++)
		{
			if (varName1[i]=='[') varName1[i]='_';
			if (varName1[i]==']') varName1[i]='_';
		}
		cout << varName1 <<endl;
				
		VarList.push_back(varName1);
		if (varName1.c_str()!="")
		{
			varIniz.push_back(varIniz1);
			varFine.push_back(varFine1);
			nv++;
		}
		else break;
	}
	in.close();

	const int N_DATA=varFine[0]-varIniz[0]+1;
	cout << "Number of variables = "<<nv<<endl;
	cout << "Number of entries = "<<N_DATA<<endl<<endl;
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	cout<<"Defining the list of variables for Ntuples:"<<endl;
	string ListVarNtu;
	for (int i=0; i<nv;i++)
	{
		if (i<(nv-1)) ListVarNtu+=VarList[i]+":";
		else ListVarNtu+=VarList[i];
	}
	cout<<"ListVarNtu = "<<ListVarNtu<<endl<<endl;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	NtuSum = new TNtuple("NtuSum","NtuSum",ListVarNtu.c_str());

	int filecount=1;
	string FileName;
	while(1)
	{
		stringstream ss;
		ss << "CODAchain" << filecount << ".txt";
		FileName=ss.str();
		ifstream in (FileName.c_str(), ios::in);
		if (in.good())
		{
			cout<<"Reading file: "<<ss.str()<<endl;

			string NtuName = "NtuChain"+Int2String(filecount);
			NtupleList.push_back(NtuName);

			TNtuple *ntu = new TNtuple (NtuName.c_str(),NtuName.c_str(),ListVarNtu.c_str());
			float *array = new float[nv];
			cout<<"Filling Ntuple: "<<NtuName<<endl<<endl;

			vector <double> var;
			int indx1;
			double var1;
			while (1) {
				in >> indx1 >>  var1;
				if (in.eof()==1) break;
				var.push_back(var1);
			}

			for (int i=0; i<N_DATA; i++) {
				
				for (int j=0; j<nv; j++) {
					  array[j]=var[i+varIniz[j]-1];
				}				
				ntu->Fill(array);
				NtuSum->Fill(array);
				double perc = 100*i/N_DATA;
				cout << "\r" << perc << "%\t\t";
			}
			cout<<endl;
			    
			filecount++;
			in.close();
			delete [] array;
		}
		else break;
	}

	cout<<"Writing ntuples in file: "<<file->GetName()<<endl<<endl;
	file->Write();
	file->Close();
	return;
}

void TracePlots (string OutputRootFile, vector <string> &VarList, string NtuName)
{
    TFile *file = new TFile(OutputRootFile.c_str(),"UPDATE");
	int nv = VarList.size();
		
	TDirectory *TracePlots = file->mkdir("Trace Plots");
	TracePlots->cd(); 

	TNtuple *ntu = (TNtuple*) file->Get(NtuName.c_str());
	float *array = new float[nv];

	vector <TGraph*> TracePlot;        
	for (int k=0; k<nv; k++)
	{
		ntu->SetBranchAddress(VarList[k].c_str(),&array[k]);
		TGraph* gr1 = new TGraph();
		TracePlot.push_back(gr1);
	}

	for (int j=0; j<ntu->GetEntries(); j++)
	{
		ntu -> GetEntry(j);
		for (int k=0; k<nv; k++)
		{
			TracePlot[k]->SetPoint(j,j,array[k]);
		}     
	}
	
	TCanvas *cT;
	int PadNum=0;
	for (int k=0; k<TracePlot.size(); k++)
	{
		PadNum++;			
		if ((k%n_pad)==0)
		{
			string CanvasName = "TracePlot";
			cT= new TCanvas(CanvasName.c_str(),CanvasName.c_str(),CanvasPx,CanvasPy,CanvasNx,CanvasNy);
			cT->Divide(n_columns,n_lines);
		}
		cT->cd(PadNum);
		
		TGraph *gr1 = TracePlot[k];
		string title="Trace Plot "+VarList[k];
		gr1->SetTitle(title.c_str());
		gr1->GetXaxis()->SetTitle("Extraction Sequence");
		gr1->GetYaxis()->SetTitle(VarList[k].c_str());
		gr1->GetXaxis()->CenterTitle();
		gr1->GetYaxis()->CenterTitle();
		gr1->Draw("AL");
		gPad->Update();
		
		if (((k+1)%n_pad)==0 || (k+1)==TracePlot.size())
		{
			sleep(SleepSeconds);
			cT->Write();
			delete cT;
			PadNum=0;
		}			
	}

	for (int k=0; k<TracePlot.size(); k++) delete TracePlot[k]; 
	file->Close();
	return;
}


void CalcNtupleStat (TNtuple *ntu, vector<string> &VarList, vector<double> &mean, vector<double> &RMS)
{
	const int nv = VarList.size();
	float *array = new float[nv];

	for (int k=0; k<nv; k++)
	{
		mean[k]=0;
		RMS[k]=0;
		ntu->SetBranchAddress(VarList[k].c_str(),&array[k]);
	}

	for (int j=0; j<ntu->GetEntries(); j++)
	{
		ntu -> GetEntry(j);
		for (int k=0; k<nv; k++)
		{
			mean[k] += array[k];
			RMS[k] += array[k]*array[k];
		}     
	}

	cout<<"VAR\tMEAN\t\tSTD.DEV\t\tErr(%)\n";
	cout.precision(2);

	for (int k=0; k<nv; k++)
	{
		mean[k] /= ntu->GetEntries();
		RMS[k] /= ntu->GetEntries();
		RMS[k] -= mean[k]*mean[k];
		RMS[k] = sqrt (RMS[k]);

		cout.precision(3);
		cout <<scientific<<VarList[k]<<"\t"<<mean[k]<<"\t"<<RMS[k] <<"\t";
		cout.precision(2);
		cout<<fixed<< 100*(RMS[k]/mean[k])<<"%"<< endl;      
	}  
	delete [] array;
}

void PrintOutputSimple (vector<string> &VarList, vector<double> &mean, vector<double> &RMS, ofstream &out)
{
	const int nv = VarList.size();

	out<<"VAR\tMEAN\t\tSTD.DEV\t\tErr(%)\n";
	out.precision(2);
	
	for (int k=0; k<nv; k++)
	{
		out.precision(3);
		out <<scientific<<VarList[k]<<"\t"<<mean[k]<<"\t"<<RMS[k] <<"\t";
		out.precision(2);
		out<<fixed<< 100*(RMS[k]/mean[k])<<"%"<< endl;        
	}	
	return;
}

#endif

