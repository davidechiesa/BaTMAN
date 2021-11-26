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
#include "TGraphErrors.h"

#include "TFile.h"
#include "TDirectory.h"

#include "TApplication.h"
#include "TGaxis.h"
#include "TColor.h"
#include "TCut.h"

#include "PosteriorAnalysis.h"
#include "MyLib.h"

using namespace std;

void ReadInput (string filename, vector<string> &ReactionName, vector<double> &obsR, vector<double> &sigmaR, vector<double> &AtomMass, vector<double> &IsotAb)  {
	CheckFileExists(filename);
	double AM,IA,ASS,err;
	double R,errR; // Activation rate divided by number of target isotopes (R/N)
	double N_avog=6.0221413e23;

	string reaction, SFfile;
	string line;
	int Nreactions=0;

	ifstream in (filename.c_str()); 

	cout<<"\nReading input file: "<<filename<<endl;
  
	cout.width(15);
	cout<<"Reaction";
	cout.width(15);
	cout<<"AtomicMass";
	cout.width(15);
	cout<<"Isot.Abund.";
	cout.width(15);
	cout<<"S.S.A.[Bq/g]";
	cout.width(15);
	cout<<"AbsErr";
	cout.width(15);
	cout<<"obsR";
	cout.width(15);
	cout<<"sigmaR";
	cout.width(20);
	cout<<"SelfShield"<<endl;
    
	while(getline(in,line))
	{
		if (line[0]=='#'||line[0]=='\0') continue;
		stringstream str;
		str<<line;
		str >> reaction >> ASS >> err >> AM >> IA >> SFfile;
		
		ReactionName.push_back(RemoveExtension(reaction));
		obsR.push_back(ASS);
		sigmaR.push_back(err);
		AtomMass.push_back(AM);
		IsotAb.push_back(IA);

		cout.width(15);
		cout<<reaction;
		cout.width(15);
		cout<<fixed<<AM;
		cout.width(15);
		cout<<IA;
		cout.width(15);
		cout<<scientific<<ASS;
		cout.width(15);
		cout<<err<<endl;

		Nreactions++;
	}
	in.close();
	cout<<"Number of activation reactions: "<<Nreactions<<endl;

	return;
}

vector<double> FillEgroups (string filename) {
	CheckFileExists(filename);
	int M=0;
	double x;	
	string line;
	vector<double> EGroups;	
	ifstream in(filename.c_str()); 
  
	cout<<"\nReading file of energy groups: "<<filename<<endl;
 
	while(getline(in,line))	{
		if (line[0]=='#'||line[0]=='\0') continue;
		stringstream(line) >> x;
		EGroups.push_back(x);
		M++;
  	}
  	in.close();
  	cout<<"\tN_groups: "<<M-1<<endl;
  	return EGroups;
}

double fitfunz(double *x, double *p)
{
	double PI=4*atan(1);
	double arg=(x[0]-p[0])/p[1];
	return (exp(-0.5*arg*arg)/sqrt(2*PI*p[1]*p[1]));
}


double fitfunzExpo(double *x, double *p)
{ 
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
	gStyle->SetOptFit(0);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetOptStat(kFALSE); 
	//TGaxis::SetMaxDigits(3);
	
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

void MarginalPDF (string OutputRootFile, vector <string> &VarList, TNtuple* ntu, ofstream &out)
{
    TFile *file = new TFile(OutputRootFile.c_str(),"UPDATE");
	int nv = VarList.size();
	TDirectory *MargPDFdir = file->mkdir("Marginal PDF"); 
	MargPDFdir->cd();
	
	int nBins=50;
	double xMin, xMax;
	bool PosDef=true;
	
	vector <TH1F*> MargPDF;
	TH1F *Marg_Dstr;
	TCanvas *c1 = new TCanvas();
       
	for (int k=0; k<nv; k++)
	{
		xMin= ntu->GetMinimum(VarList[k].c_str());
		xMax= ntu->GetMaximum(VarList[k].c_str());

		if (xMin < 0) PosDef=false;
		
		double range = xMax - xMin;
		xMin -= 0.2*range;
		xMax += 0.2*range;
		if (xMin < 0 && PosDef==true) xMin=0;
		
		string TH1F_Name = VarList[k];
		Marg_Dstr = new TH1F(TH1F_Name.c_str(), TH1F_Name.c_str(), nBins,xMin,xMax);

		string DrawOpt = VarList[k] +">>" + TH1F_Name;
		ntu->Draw(DrawOpt.c_str());
		
		MargPDF.push_back(Marg_Dstr);
	}
	
	DrawMarginalPDF(file, MargPDF, out);
	delete c1;	
	file->Close();	
	return;
}


void DrawMarginalPDF(TFile* &file, vector <TH1F*> MargPDF, ofstream &out)
{
	
	TCanvas *cT;
	int PadNum=0;	
	
	double probability = 0.90;
	out.precision(0);
	out << "\n\n";

	out<<left<<setw(20)<<" "<<setw(64)<<"GAUS FIT RESULTS"<<fixed<<probability*100<<"% CONFIDENCE RANGE\n";
	out.width(24);
	out <<left<<setw(20)<< "Variable Name"<<setw(8)<<"Func"<<setw(16)<<"p0"<<setw(16)<<"p1"<<setw(8)<<"Chi2-r"<<setw(16)<<"p-value"<<setw(16)<<"Minimum"<<"Maximum"<<endl;
	
	TF1 *f1=new TF1("gaus_norm",fitfunz,-100,100,2);
	f1->SetLineColor(kBlue);
	f1->SetLineWidth(2);
	TF1 *f2=new TF1("expo_norm",fitfunzExpo,0,100,1);
	f2->SetLineColor(kRed);
	f2->SetLineWidth(2);
	int GausFlag, ExpFlag;
	
	for (int k=0; k<MargPDF.size(); k++)
	{
		GausFlag=0;
		ExpFlag=0;
		PadNum++;			
		if ((k%n_pad)==0)
		{
			string CanvasName = "Marginal-PDF";
			cT= new TCanvas(CanvasName.c_str(),CanvasName.c_str(),CanvasPx,CanvasPy,CanvasNx,CanvasNy);
			cT->Divide(n_columns,n_lines);
		}
		cT->cd(PadNum);
		TH1F *Marg_Dstr = MargPDF[k];
		double norm=Marg_Dstr->Integral("width");
		double bincontent, binerror;
		for (int i=0; i<Marg_Dstr->GetNbinsX(); i++)
		{
			bincontent = Marg_Dstr -> GetBinContent(i+1);
			binerror = Marg_Dstr -> GetBinError(i+1);
			Marg_Dstr -> SetBinContent(i+1, bincontent/norm);
			Marg_Dstr -> SetBinError(i+1, binerror/norm);
		}
		int binlow, binup;
		TH1F *yellow = HighlightedRegion (Marg_Dstr, probability, binlow, binup);
		double LowEdge = yellow->GetBinLowEdge(binlow+1);
		double UpEdge = yellow->GetBinLowEdge(binup+2);
				
		f1->SetParameter(0, Marg_Dstr->GetMean());
		f1->SetParameter(1, Marg_Dstr->GetRMS());
	    
		Marg_Dstr->Fit("gaus_norm","0");
		double prob = f1 -> GetProb ();
		double redchisq = f1 -> GetChisquare() / f1->GetNDF();
		double mu = f1->GetParameter (0);
		double sigma = f1->GetParameter (1);
		if (prob>=0.01) GausFlag=1;			
		
		yellow -> Draw("H");
		yellow -> SetFillColor(kYellow);
		yellow->GetXaxis()->SetTitle(Marg_Dstr->GetName());
		yellow->GetYaxis()->SetTitle("Marginal p.d.f.");
		yellow->GetXaxis()->CenterTitle();
		yellow->GetYaxis()->CenterTitle();
		yellow->GetXaxis()->SetTitleSize(0.05);
		yellow->GetYaxis()->SetTitleSize(0.05);
		yellow->GetXaxis()->SetLabelSize(0.05);
		yellow->GetYaxis()->SetLabelSize(0.05);
		yellow->GetXaxis()->SetNdivisions(506);		
					
		Marg_Dstr->SetFillStyle(0);
		Marg_Dstr->Draw("sameHE");
		
		out <<setw(20)<<left<<Marg_Dstr->GetTitle();
		
		if (prob >= 0.01 && GausFlag==1) {
			Marg_Dstr->Fit("gaus_norm");
			out <<left<<setw(8)<< "gaus_n";
			out.precision(3);
			out <<left<<scientific<<setw(16)<< mu << setw(16) << sigma;

		}
		else  {
			Marg_Dstr->Fit("gaus");
			TF1 *gausfunc = Marg_Dstr->GetFunction("gaus");
			out <<left<<setw(8)<< "*Gaus";
			out.precision(3);
			mu=gausfunc->GetParameter(1);
			sigma=gausfunc->GetParameter(2);
			out <<left<<scientific<<setw(16)<< mu << setw(16) << sigma;
			prob = gausfunc -> GetProb ();
			redchisq = gausfunc -> GetChisquare() / gausfunc->GetNDF();
		}
		
		out.precision(2);
		out <<fixed<<setw(8)<< redchisq <<setw(16)<< prob;
		out <<scientific<<setw(16)<< LowEdge<<setw(16)<< UpEdge <<endl;   			

		gPad->Update();
		
		if (((k+1)%n_pad)==0 || (k+1)==MargPDF.size())
		{
			sleep(SleepSeconds);
			cT->Write();
			string PrintName;
			if (((k+1)%n_pad)==0) PrintName = "MargPDF"+Int2String((k+1)/n_pad)+".png";
			else PrintName = "MargPDF"+Int2String(1+(k+1)/n_pad)+".png";
			cT->Print(PrintName.c_str(),"png");
			delete cT;
			PadNum=0;
		}			
	}
	out << "\n\n";
	for (int k=0; k<MargPDF.size(); k++) delete MargPDF[k]; 
	
	return;
}

TH1F* HighlightedRegion (TH1F* Original, double prob, int& binlow, int& binup)
{
	TH1F *yellow = (TH1F*)Original->Clone("yellow");
	vector<double> BinArea;
	for (int i=0; i<Original->GetNbinsX(); i++)
		BinArea.push_back(Original -> GetBinContent(i+1) * Original -> GetBinWidth(i+1));
	int nbins=2;
	double sum=0;
	double MaxSum=0;
	while (1)
	{
		MaxSum=0;
		for (int j=0; j<=(Original->GetNbinsX() - nbins); j++)
		{
			sum=0;
			for (int i=j; i<(j+nbins); i++)
				sum +=BinArea[i];
			if ((sum>=prob) && (sum>MaxSum))
			{
				binlow=j;
				binup=j+nbins-1;
				MaxSum=sum; 
				
			}  
		}
		if (MaxSum > prob) break;
		nbins++;
	}
	//cout << nbins << "\t" <<binlow << "\t"<<binup << "\t" << MaxSum << endl; 
	
	for (int i=0; i<yellow->GetNbinsX(); i++)
	{
		if ((i<binlow) || (i>binup)) 
		{
			yellow ->SetBinContent(i+1, 0.);
			yellow ->SetBinError(i+1, 0.);
		}
	}	
	//cin.get();
	return yellow;
}

void CorrelationPlots (string OutputRootFile, vector <string> &VarList, TNtuple* ntu, double** CorrCoeffMatrix, ofstream &out)
{
    TFile *file = new TFile(OutputRootFile.c_str(),"UPDATE");
    //TNtuple* ntu = (TNtuple*) file->Get("NtuSum"); 
	int nv = VarList.size();	
		
	TDirectory *ScatPlots = file->mkdir("Correlation Plots");
	ScatPlots->cd();
	
	double xMin, yMin, xMax, yMax;
	int NBin_SC=100;
	string ScatPL;
	double comb=0.5*nv*(nv-1);  
	gStyle->SetOptStat(kFALSE);
	//TGaxis::SetMaxDigits(3);
 
	//TNtuple *ntu = (TNtuple*) file->Get(NtuName.c_str());
	vector <TH2F*> ScatPDF;
	
	//double CorrCoeffMatrix [nv][nv];
	for (int i=0; i<nv; i++) CorrCoeffMatrix[i][i]=1.;
	
	TCanvas *c1;
	for (int k=0; k<nv; k++)
	{
		for (int j=0; j<nv; j++)
		{
			if (j>k){
				c1 = new TCanvas();
				xMin=ntu->GetMinimum(VarList[k].c_str());
				xMax=ntu->GetMaximum(VarList[k].c_str());
				yMin=ntu->GetMinimum(VarList[j].c_str());
				yMax=ntu->GetMaximum(VarList[j].c_str());
				string TH2F_Name = "Scat_" + VarList[k] + "vs" + VarList[j];
				TH2F *Scat_plot=new TH2F(TH2F_Name.c_str(), TH2F_Name.c_str(), 
							NBin_SC, xMin, xMax, NBin_SC, yMin, yMax);
										
				string DrawOpt = VarList[j]+":"+VarList[k]+">>" + TH2F_Name;
				ntu->Draw(DrawOpt.c_str());
				ScatPDF.push_back(Scat_plot);
				CorrCoeffMatrix[j][k] = Scat_plot->GetCorrelationFactor();
				CorrCoeffMatrix[k][j] = Scat_plot->GetCorrelationFactor();
				delete c1;
			}
		}				
	}

	TCanvas *cT;
	int PadNum=0;
	for (int k=0; k<ScatPDF.size(); k++)
	{
		PadNum++;			
		if ((k%n_pad)==0)
		{
			string CanvasName = "ScatterPlot";
			cT= new TCanvas(CanvasName.c_str(),CanvasName.c_str(),CanvasPx,CanvasPy,CanvasNx,CanvasNy);
			cT->Divide(n_columns,n_lines);
		}
		cT->cd(PadNum);
		TH2F *Scat_plot = ScatPDF[k];
		Scat_plot->Draw("COLZ");
		gPad->Update();
		
		if (((k+1)%n_pad)==0 || (k+1)==ScatPDF.size())
		{
			sleep(SleepSeconds);
			cT->Write();
			delete cT;
			PadNum=0;
		}			
	}

	for (int k=0; k<ScatPDF.size(); k++) delete ScatPDF[k]; 
	
	// draw correlation matrix 
		
	file->cd();
	string CanvasName = "Correlation-Matrix";
	cT= new TCanvas(CanvasName.c_str(),CanvasName.c_str(),0,0,1200,1200);
	cT->SetRightMargin(0.12); 
	//cT->SetLeftMargin(0.1); 
	cT->SetTopMargin(0.12); 
	//cT->SetBottomMargin(0.1); 
	
	TH2F *CorrMatrix = new TH2F("CorrMatrix", "Correlation Factors", nv, 0.5, nv+0.5, nv, 0.5, nv+0.5);
	out<<"\n////////////////////////////////////////////////////\n";
	out<<"/////////     Correlation Coefficients     /////////\n";
	out<<"////////////////////////////////////////////////////\n\n";
	out.width(10);
	out<<" ";
	for (int i=0; i<nv; i++)
	{
		out.width(10);
		out<<right<<VarList[i];
	}
	out<<"\n";

	out.precision(3);
	for (int i=0; i<nv; i++)
	{
		out.width(10);
		out << VarList[i];
		for (int j=0; j<nv; j++)
		{
			out.width(10);
			if (i==j) 
			{
				CorrMatrix->Fill(i+0.5, j+0.5, CorrCoeffMatrix[i][j]);
				CorrCoeffMatrix[i][j]=1.;			
			}
			else CorrMatrix->Fill(i+0.5, j+0.5, CorrCoeffMatrix[i][j]);
			out<<fixed<<right<<CorrCoeffMatrix[i][j];			
		}
		out<<endl;
	}
	out<<endl;
	
	Int_t nb=100;
	Int_t Palette[nb];
	const Int_t Number = 3;
	Double_t Red[Number]   = {0., 1.0, 1.0};
	Double_t Green[Number] = {0., 1.0, 0.};
	Double_t Blue[Number]  = {1.0, 1.0, 0.};
	Double_t Stops[Number] = {0., .50, 1.0};
	Int_t FI = TColor::CreateGradientColorTable(Number,Stops,Red,Green,Blue,nb);
	for (int i=0;i<nb;i++) Palette[i] = FI+i;
	gStyle->SetPalette(nb, Palette);

	CorrMatrix->Draw("colz");
	CorrMatrix->GetXaxis()->SetNdivisions(115);
	CorrMatrix->GetYaxis()->SetNdivisions(115);
	CorrMatrix->GetZaxis()->SetRangeUser(-1.,1.1);
	CorrMatrix->GetXaxis()->SetTitle("Flux Group");
	CorrMatrix->GetYaxis()->SetTitle("Flux Group");
	cT->SetGridx();
	cT->SetGridy();
	
	cT->Write();
	cT->Print("CorrelationMatrix.pdf", "pdf");
	delete cT;	
	file->Close();
	return;
}



void CorrelRatePHI(string OutputRootFile, vector <string> &PHIList, vector <string> &RateList,  TNtuple* ntu, double** CorrCoeffMatrix, ofstream &out)
{
    TFile *file = new TFile(OutputRootFile.c_str(),"UPDATE");

	int nRows = PHIList.size();	
	int nColumns = RateList.size();	
	
	double xMin, yMin, xMax, yMax;
	int NBin_SC=100;
	string ScatPL;
	
	gStyle->SetOptStat(kFALSE);
	//vector <TH2F*> ScatPDF;
	TCanvas *c1;
	for (int k=0; k<nRows; k++)
	{
		for (int j=0; j<nColumns; j++)
		{
			c1 = new TCanvas();
			xMin=ntu->GetMinimum(PHIList[k].c_str());
			xMax=ntu->GetMaximum(PHIList[k].c_str());
			yMin=ntu->GetMinimum(RateList[j].c_str());
			yMax=ntu->GetMaximum(RateList[j].c_str());
			string TH2F_Name = "Scat_" + PHIList[k] + "vs" + RateList[j];
			TH2F *Scat_plot=new TH2F(TH2F_Name.c_str(), TH2F_Name.c_str(), 
						NBin_SC, xMin, xMax, NBin_SC, yMin, yMax);
									
			string DrawOpt = RateList[j]+":"+PHIList[k]+">>" + TH2F_Name;
			ntu->Draw(DrawOpt.c_str());
			//ScatPDF.push_back(Scat_plot);
			CorrCoeffMatrix[k][j] = Scat_plot->GetCorrelationFactor();
			delete c1;
		}				
	}
		
	file->cd();
	TCanvas *cT;
	string CanvasName = "Correlation-ActRate_vs_PHI";
	cT= new TCanvas(CanvasName.c_str(),CanvasName.c_str(),0,0,100*nColumns,100*nRows);
	cT->SetRightMargin(0.12); 
	//cT->SetLeftMargin(0.1); 
	cT->SetTopMargin(0.12); 
	//cT->SetBottomMargin(0.1); 
	
	TH2F *CorrMatrix = new TH2F("CorrMatrix", "Correlation Factors", nColumns, 0.5, nColumns+0.5, nRows, 0.5, nRows+0.5);
	out<<"\n//////////////////////////////////////////////////\n";
	out<<"/////////     Correlation ActRate_vs_PHI    /////////\n";
	out<<"////////////////////////////////////////////////////\n\n";
	out.width(10);
	out<<" ";
	for (int j=0; j<nColumns; j++)
	{
		out.width(10);
		out<<right<<RateList[j];
	}
	out<<"\n";

	out.precision(2);
	for (int k=0; k<nRows; k++)
	{
		out.width(10);
		out << PHIList[k];
		for (int j=0; j<nColumns; j++)
		{
			out.width(10);
			CorrMatrix->Fill(j+0.5, k+0.5, CorrCoeffMatrix[k][j]);
			out<<fixed<<right<<CorrCoeffMatrix[k][j];			
		}
		out<<endl;
	}
	out<<endl;
	
	Int_t nb=100;
	Int_t Palette[nb];
	const Int_t Number = 3;
	Double_t Red[Number]   = {0., 1.0, 1.0};
	Double_t Green[Number] = {0., 1.0, 0.};
	Double_t Blue[Number]  = {1.0, 1.0, 0.};
	Double_t Stops[Number] = {0., .50, 1.0};
	Int_t FI = TColor::CreateGradientColorTable(Number,Stops,Red,Green,Blue,nb);
	for (int i=0;i<nb;i++) Palette[i] = FI+i;
	gStyle->SetPalette(nb, Palette);

	CorrMatrix->Draw("colz");
	CorrMatrix->GetXaxis()->SetNdivisions(511);
	CorrMatrix->GetYaxis()->SetNdivisions(511);
	CorrMatrix->GetZaxis()->SetRangeUser(-1.,1.1);
	CorrMatrix->GetXaxis()->SetTitle("Reaction Number");
	CorrMatrix->GetYaxis()->SetTitle("Flux Group");
	cT->SetGridx();
	cT->SetGridy();
	
	cT->Write();
	cT->Print("CorrMatrix-R_vs_PHI.pdf", "pdf");
	delete cT;	
	file->Close();
	return;
}

void CorrelRates(string OutputRootFile, vector <string> &RateList,  TNtuple* ntu, double** CorrCoeffMatrix, ofstream &out)
{
        TFile *file = new TFile(OutputRootFile.c_str(),"UPDATE");

	int nRows = RateList.size();	
	int nColumns = RateList.size();	
	
	double xMin, yMin, xMax, yMax;
	int NBin_SC=100;
	string ScatPL;
	
	gStyle->SetOptStat(kFALSE);
	TCanvas *c1;
	for (int k=0; k<nRows; k++)
	{
		for (int j=0; j<nColumns; j++)
		{
			c1 = new TCanvas();
			xMin=ntu->GetMinimum(RateList[k].c_str());
			xMax=ntu->GetMaximum(RateList[k].c_str());
			yMin=ntu->GetMinimum(RateList[j].c_str());
			yMax=ntu->GetMaximum(RateList[j].c_str());
			string TH2F_Name = "Scat_" + RateList[k] + "vs" + RateList[j];
			TH2F *Scat_plot=new TH2F(TH2F_Name.c_str(), TH2F_Name.c_str(), 
						NBin_SC, xMin, xMax, NBin_SC, yMin, yMax);
									
			string DrawOpt = RateList[j]+":"+RateList[k]+">>" + TH2F_Name;
			ntu->Draw(DrawOpt.c_str());
			CorrCoeffMatrix[k][j] = Scat_plot->GetCorrelationFactor();
			delete c1;
		}				
	}
		
	file->cd();
	TCanvas *cT;
	string CanvasName = "Correlation-ActRates";
	cT= new TCanvas(CanvasName.c_str(),CanvasName.c_str(),0,0,100*nColumns,100*nRows);
	cT->SetRightMargin(0.12); 
	//cT->SetLeftMargin(0.1); 
	cT->SetTopMargin(0.12); 
	//cT->SetBottomMargin(0.1); 
	
	TH2F *CorrMatrix = new TH2F("CorrMatrix", "Correlation Factors", nColumns, 0.5, nColumns+0.5, nRows, 0.5, nRows+0.5);
	out<<"\n//////////////////////////////////////////////////\n";
	out<<"/////////       Correlation ActRates       /////////\n";
	out<<"////////////////////////////////////////////////////\n\n";
	out.width(10);
	out<<" ";
	for (int j=0; j<nColumns; j++)
	{
		out.width(10);
		out<<right<<RateList[j];
	}
	out<<"\n";

	out.precision(2);
	for (int k=0; k<nRows; k++)
	{
		out.width(10);
		out << RateList[k];
		for (int j=0; j<nColumns; j++)
		{
			out.width(10);
			CorrMatrix->Fill(j+0.5, k+0.5, CorrCoeffMatrix[k][j]);
			out<<fixed<<right<<CorrCoeffMatrix[k][j];			
		}
		out<<endl;
	}
	out<<endl;
	
	Int_t nb=100;
	Int_t Palette[nb];
	const Int_t Number = 3;
	Double_t Red[Number]   = {0., 1.0, 1.0};
	Double_t Green[Number] = {0., 1.0, 0.};
	Double_t Blue[Number]  = {1.0, 1.0, 0.};
	Double_t Stops[Number] = {0., .50, 1.0};
	Int_t FI = TColor::CreateGradientColorTable(Number,Stops,Red,Green,Blue,nb);
	for (int i=0;i<nb;i++) Palette[i] = FI+i;
	gStyle->SetPalette(nb, Palette);

	CorrMatrix->Draw("colz");
	CorrMatrix->GetXaxis()->SetNdivisions(511);
	CorrMatrix->GetYaxis()->SetNdivisions(511);
	CorrMatrix->GetZaxis()->SetRangeUser(-1.,1.1);
	CorrMatrix->GetXaxis()->SetTitle("Reaction Number");
	CorrMatrix->GetYaxis()->SetTitle("Reaction Number");
	cT->SetGridx();
	cT->SetGridy();
	
	cT->Write();
	cT->Print("CorrMatrix-ActRates.pdf", "pdf");
	delete cT;	
	file->Close();
	return;
}

