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

#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#include "TGraphErrors.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"

#include "NucReaction.h"
#include "MyLib.h"


using namespace std;

//ctor
NucReaction::NucReaction (string name) {
	name_p=name;
	grXS_p=new TGraphErrors();
	Flux_p=new TH1F();
	SelfShieding_p=new TH1F();

	NormPhi_p=0.;
	FluxTOT_p=0.;
	Ngroups_p=0.;
	
	XSeff_p=0.;
	ErrCorr_p=0.;
	ErrUncorr_p=0.;
	PHI_p=0.;
	ErrPHI_p=0.;
}

//dtor
NucReaction::~NucReaction ()  {
	//cout<<"Dtor!"<<endl;
	delete grXS_p;
	delete Flux_p;	
	delete SelfShieding_p;
	
	delete[] PHIgr_p;
	delete[] errPHIgr_p;
	delete[] XSgr_p;	
	delete[] errXSgr_p;
	delete[] errCorrXSgr_p;	
	delete[] errUncorrXSgr_p;
	delete[] GroupFrac_p;	
}

void NucReaction::FillXS (string filename)  {
	
	CheckFileExists(filename);
	
	int NdataXS=0;
	double x,y, xPre;
	string line;
	ifstream in(filename.c_str()); 

	cout<<"\nReading XS file: "<<filename<<endl;
	cout<<"\tConversion eV -> MeV"<<endl;

	while(getline(in,line))	{
		if (line[0]=='#'||line[0]=='\0') continue;
		stringstream(line) >> x >> y;
		//cout<<x<<"\t"<<y<<"\n";
		if ((NdataXS!=0) && (x==xPre)) continue;

		grXS_p -> SetPoint (NdataXS, x*1.0e-6 , y);
		grXS_p -> SetPointError (NdataXS, 0. , 0.);

		NdataXS++;
		xPre=x;
	}
	in.close();
	cout<<"\tN XS data: "<<NdataXS<<endl;
	return;	
}


void NucReaction::FillXSwithErr (string filename)  {
	
	CheckFileExists(filename);
	
	int NdataXS=0;
	double x,y,ey, xPre;
	string line;
	ifstream in(filename.c_str()); 

	cout<<"\nReading XS file: "<<filename<<endl;
	cout<<"\tConversion eV -> MeV"<<endl;

	while(getline(in,line))	{
		if (line[0]=='#'||line[0]=='\0') continue;
		stringstream(line) >> x >> y >> ey;
		//cout<<x<<"\t"<<y<<"\n";
		if ((NdataXS!=0) && (x==xPre)) continue;

		grXS_p -> SetPoint (NdataXS, x*1.0e-6 , y);
		grXS_p -> SetPointError (NdataXS, 0. , ey);

		NdataXS++;
		xPre=x;
	}
	in.close();
	cout<<"\tN XS data: "<<NdataXS<<endl;
	return;	
}

void NucReaction::SetXSErrors (double errRelTHERM, double errRelRES)  {
	double CdCutoff = 0.5e-6; // 0.5 eV
	double E, XS;
	for (int i=0; i < grXS_p->GetN(); i++) {
		grXS_p -> GetPoint (i, E, XS);
		if (E <= CdCutoff) grXS_p -> SetPointError (i, 0., errRelTHERM*XS);
		else 	           grXS_p -> SetPointError (i, 0., errRelRES*XS);	
	}
	return;
}

void NucReaction::FillFlux (string filename, double NormPhi)  {
	NormPhi_p=NormPhi;
	FluxTOT_p=0.;
	CheckFileExists(filename);

	int NdataFlux=0.;
	double x,y,z;
	string line;
	ifstream in (filename.c_str()); 
	
	vector <double> Ebin, phi, err;
	cout<<"\nReading flux file: "<<filename<<endl;
	cout<<"\tFlux norm factor: "<<NormPhi<<endl;
	while(getline(in,line))  {
		if (line[0]=='#'||line[0]=='\0') continue;
			stringstream(line) >> x >> y >> z;
			Ebin.push_back(x);
			phi.push_back(y*NormPhi_p);
			err.push_back(z);
			FluxTOT_p += (y*NormPhi_p);
			NdataFlux++;
	}
	in.close();
	
	double *E = new double [NdataFlux];
	for (int i=0; i<Ebin.size(); i++) E[i]=Ebin[i];
	
	int Nbins = NdataFlux-1 ;
	string hname = "h"+name_p;
	string hnameSF = "SF"+name_p;

	TH1F* h1 = new TH1F(hname.c_str(), hname.c_str(), Nbins ,E); 
	TH1F* h2 = new TH1F(hnameSF.c_str(), hnameSF.c_str(), Nbins ,E); 
	for (int i=1; i<=Nbins; i++) {
		h1->SetBinContent(i, phi[i]);
		h1->SetBinError(i, err[i]*phi[i]);
		h2->SetBinContent(i, 1.);
		h2->SetBinError(i, 0.);
	}
	
	*Flux_p = *h1;
	*SelfShieding_p = *h2;
		
	cout<<"\tN bins flux: "<< Nbins <<"\tFlux total: "<<FluxTOT_p<<endl;
	delete [] E;
	delete h1;
	delete h2;
	return;
}

void NucReaction::SetSelfShielding (string filename)  {
	CheckFileExists(filename);
	int Nbins = Flux_p->GetNbinsX();
	string line;
	if (Nbins<2) {
		cout << "Error: FillFlux method must be applied before SetSelfShielding."<< endl;
		cout << "Nbins flux = "<<Nbins<< endl;
		getline (cin, line);
		getline (cin, line);
		return;
	}
	
	int Ndata=0.;
	double x,y,z;
	
	ifstream in (filename.c_str()); 
	
	vector <double> Ebin, SF, err;
	cout<<"\nReading self-shielding file: "<<filename<<endl;
	
	while(getline(in,line))  {
		if (line[0]=='#'||line[0]=='\0') continue;
			stringstream(line) >> x >> y; // >> z;
			if (x>=1) break; // self-shielding correction not applied for E > 1 MeV 
			Ebin.push_back(x);
			SF.push_back(y);
			//err.push_back(z);
			Ndata++;
	}
	in.close();
	
	double LowEgde;
	double epsilon=1e-4;
	for (int i=0; i<Ndata; i++) {
		LowEgde=Flux_p->GetBinLowEdge(i+1);
		if ( fabs(LowEgde-Ebin[i]) > epsilon*LowEgde) {
			cout << "Error: Flux binning is different with respect to Self-shielding binning:"<< endl;
			cout << "E_Flux["<<i<<"] = "<<LowEgde<< endl;
			cout << "Self-sh["<<i<<"] = "<<Ebin[i]<< endl;
			getline (cin, line);
			getline (cin, line);
			return;
		}
	}
	
	for (int i=1; i<=Nbins; i++)  {
		if (i<SF.size()) SelfShieding_p->SetBinContent(i, SF[i]);
		else SelfShieding_p->SetBinContent(i, 1.);	
	}
	return;	 
}

double NucReaction::CalcXSeff(double min, double max, double& ErrCorr, double& ErrUncorr, double& Flux, double& ErrFlux) {
	cout.precision(2);
	cout<<scientific<<"\nCalculating effective cross section in the range: [" << min << ", " << max <<"] MeV"<< endl;
	
	///////////////////////////////////////////////////////////////////////////////
	// Compatibility check of energy ranges of XSeff, FLUX and XSDATA
	///////////////////////////////////////////////////////////////////////////////
	
	int Nbins = Flux_p->GetNbinsX();
	//cout<< "Nbins = " <<Nbins <<endl;
	double FluxMin= Flux_p->GetBinLowEdge(1);
	double FluxMax= Flux_p->GetBinLowEdge(Nbins+1);
	//cout<< "Flux spectrum range [" << Flux_p->GetBinLowEdge(1) << ", " << Flux_p->GetBinLowEdge(Nbins+1) <<"] MeV"<< endl;
	double eps = 1e-13;
	if (FluxMin>(min+eps) || FluxMax<(max-eps)) { 
		cout << "Error! Effective cross section range exceeding flux range!"<<endl;
		cout << "Press ANY key to continue, assuming flux = 0 in the absence of data ";
		cin.ignore(256,'\n');
		cin.get();
	}
	int binMin = Flux_p->Fill(min, 0.);
	int binMax = Flux_p->Fill(max, 0.);	
	if (FluxMax<=(max+eps) ) binMax=Nbins; // FluxMax<=max
	if (FluxMin>=(min-eps) ) binMin=1; // FluxMin>=min
	cout<< "\tRange of flux bins: min = " <<binMin << "\tmax = "<<binMax<<endl;
	
		
	int XSpoints = grXS_p->GetN();
	double EXSmin, EXSmax, XS;
	grXS_p->GetPoint(0, EXSmin, XS);
	grXS_p->GetPoint(XSpoints-1, EXSmax, XS);
	//cout<< "XS data range [" << EXSmin << ", " << EXSmax <<"] MeV"<< endl;	
	
	
	if (EXSmin>(min+eps) || EXSmax<(max-eps)) { // maggiore e minore "stretti"	
		cout << "Warning! Effective cross section range exceeding XS data range" << endl;
		cout << "\t => assuming XS = 0 in the absence of XS data" << endl;
		//cout << "Press ANY key to continue, assuming XS = 0 in the absence of XS data ";
		//cin.ignore(256,'\n');
		//cin.get();
		if (EXSmax<(max-eps)) { // adding two points equal to zero at the end of the cross section
			grXS_p->SetPoint(XSpoints, EXSmax + 1e-9*(max-EXSmax), 0.);
			grXS_p->SetPoint(XSpoints+1, max, 0.);			
			grXS_p->SetPointError(XSpoints, 0., 0.);
			grXS_p->SetPointError(XSpoints+1, 0., 0.);
		}
		if (EXSmin>(min+eps)) { // adding two points equal to zero at the beginning of the cross section
			double *energia =new double [XSpoints];
			double *sezUrto =new double [XSpoints];
			double *errSezUrto =new double [XSpoints];
			for (int i=0; i>XSpoints; i++) {
				grXS_p->GetPoint(i, energia[i], sezUrto[i]);
				errSezUrto[i] = grXS_p->GetErrorY(i);
			}			
			grXS_p->SetPoint (0, min, 0.);
			grXS_p->SetPoint (1, EXSmin - 1e-9*(EXSmin-min), 0.);			
			grXS_p->SetPointError (0, 0., 0.);
			grXS_p->SetPointError (1, 0., 0.);
			for (int i=0; i>XSpoints; i++)  {
				grXS_p->SetPoint (i+2, energia[i], sezUrto[i]);
				grXS_p->SetPointError (i+2, 0., errSezUrto[i]);
			}			
			delete [] energia;
			delete [] sezUrto;
			delete [] errSezUrto;
		}
	}
	
	
	///////////////////////////////////////////////////////////////////////////////
	// Effective cross section calculation between MIN and MAX
	///////////////////////////////////////////////////////////////////////////////
	
	int point;
	double x, y, xPre, yPre, errY, errYpre=0.;
	
	double phi, errPhi, Elow, Eup, xsLow, xsUp, integral, errIntXS, selfSh;
	double XSeff=0; 
	ErrCorr=0.;
	Flux=0;
	ErrFlux=0.;
	
	vector<double> XSmean;	// average XS value in each bin of the flux spectrum
	vector<double> DeltaPhi;  // flux error in each bin of the spectrum
	
	for (int i=binMin; i<=binMax; i++) {  // calculating the subtended area to the cross section graph in each bin of the flux spectrum
		phi = Flux_p->GetBinContent(i);
		errPhi = Flux_p->GetBinError(i);
		
		selfSh = SelfShieding_p->GetBinContent(i);
		
		Elow = Flux_p->GetBinLowEdge(i);
		Eup = Flux_p->GetBinLowEdge(i+1);
		phi /= (Eup-Elow); // flux per unit energy
		errPhi /= (Eup-Elow);
		
		if (fabs(Elow-max) < eps) continue; // Elow==max -> stop integrating because I'm at the bin with LowEdge = Max
		// In the first and last bin of the range [min, max], Elow and Eup must coincide with min and max
		if (Elow<=(min+eps)) Elow=min; // Elow<=min
		if (Eup>=(max-eps)) Eup=max; // Eup>=max
		// Interpolating the cross section at min and max 
		xsLow = grXS_p->Eval(Elow);
		xsUp = grXS_p->Eval(Eup);
		
		//cout << scientific << Elow <<"\t" << Eup << endl;
		
		if (i==binMin) point=0;		
		integral=0.;  // subtended area to the cross section graph
		errIntXS=0.; // uncertainty associated to integral
		xPre=Elow;
		yPre=xsLow;

		while (point < grXS_p->GetN()) {
			grXS_p->GetPoint(point, x, y);
			errY = grXS_p->GetErrorY(point);
			
			//cycling over the cross section points until arriving at the first bin low edge
			if ( x<=(Elow+eps) ) {  // x<=Elow  
				if (y!=0) errYpre = xsLow*errY/y;
				point++;
				continue;
			}
			
			// break the cycle when passing the up edge
			// --> the "point" value will be correctly initialized for the next bin)
			else if ( x>=(Eup-eps) ) break; // x>=Eup 

			// integral calculated with the trapezoid method
			else if ( (x>Elow+eps) && (x<Eup-eps) ) { // Elow < x < Eup
				if ((xPre/Elow-1.)<eps && y!=0 ) {  // xPre==Elow
					errYpre=xsLow*errY/y;
					//cout << "      a " << 100*errYpre/xsLow << endl;
				}
				integral += (0.5*(yPre+y)*(x-xPre)); // 
				errIntXS += (0.5*(errYpre+errY)*(x-xPre)); // correlation Y errors = 1
				xPre=x;
				yPre=y;
				errYpre = errY;
				point++;
			}
			else {
				cout <<"Error in CalcXSeff!" << endl;
				return -1;
			}			

		} 
		// adding the subtended area in the energy range between the last XS point ans the bin UpEdge 
		if (y!=0) {
			errY = xsUp*errY/y;
			//cout << "      b " << 100*errY/xsUp << endl;
		}
		integral += (0.5*(yPre+xsUp)*(Eup-xPre));
		errIntXS += (0.5*(errYpre+errY)*(Eup-xPre)); 
		errYpre=errY;
		
		//cout.precision(4);
		//cout <<fixed<< Eup <<"\t"<< 100*errIntXS/integral<<"\t" <<selfSh << endl;
		
		XSmean.push_back(integral/((Eup-Elow)));
		DeltaPhi.push_back( errPhi*(Eup-Elow) );
		
		//    For each bin, integral is filled with XS(E)*dE
		//    ---> multiply by Phi(E) (i.e. flux per unit energy) to get xsEff numerator: Phi(E)*XS(E)*dE
		//    ---> multiply by selfSh to introduce SelfShielding correction
		integral *= phi*selfSh;
		errIntXS *= phi*selfSh;
				
		XSeff += integral;		
		ErrCorr += errIntXS;  // XS uncertainty calculated assuming correlation = 1 among bins (rigidly moved up and down)
		
		// calculating xsEff denominator: Phi(E)*dE
		Flux += phi*(Eup-Elow);
		ErrFlux += pow(errPhi*(Eup-Elow),2);		
	}
	ErrFlux = sqrt(ErrFlux);
	
	//calculating XSeff
	XSeff /= Flux;
	ErrCorr /= Flux;
	
	// calculating the uncertainty associated to XSeff due to the guess flux (e.g evaluated throug a MC)
	ErrUncorr=0.; 
	double derivata;
	//cout << "\tXSmean.size() = " << XSmean.size() << endl;
	for (int i=0; i<XSmean.size(); i++) {
		derivata = (XSmean[i]-XSeff)/Flux;
		ErrUncorr += pow( (derivata*DeltaPhi[i]) , 2);
	}
	
	ErrUncorr = sqrt(ErrUncorr);
	
	// total XSeff uncertainty
	double ErrXSeff = sqrt(pow(ErrCorr,2)+pow(ErrUncorr,2)); 
	
	//cout << "Energy range: [" << min << ", " << max <<"]"<< endl;
	cout.unsetf(ios_base::scientific);
	cout.precision(2);
	cout << scientific << "Flux Intensity (a.u.): "<<Flux<<" +- "<<ErrFlux<< "\t("<< fixed <<100*ErrFlux/Flux<<"%)"<< endl;
	cout << scientific << "Effective XS (barn):   "<<XSeff<<" +- "<<ErrXSeff<< "\t("<<fixed <<100*ErrXSeff/XSeff<<"%)"<< endl;	
	cout << scientific <<"\t\tErrCorr   = "<<ErrCorr<< "\t("<<fixed <<100*ErrCorr/XSeff<<"%)"<< endl;
	cout << scientific <<"\t\tErrUncorr = "<<ErrUncorr<< "\t("<<fixed <<100*ErrUncorr/XSeff<<"%)"<< endl;
	
	return XSeff;
}

void NucReaction::CalcXSgroups (vector<double> &EGroup){	
	Ngroups_p = EGroup.size()-1;
	
	PHIgr_p=new double[Ngroups_p];
	errPHIgr_p=new double[Ngroups_p];
	XSgr_p=new double[Ngroups_p];
	errXSgr_p=new double[Ngroups_p];
	errCorrXSgr_p=new double[Ngroups_p];
	errUncorrXSgr_p=new double[Ngroups_p];
	GroupFrac_p=new double[Ngroups_p];
	
	double integral=0.;
	
	for (int i=0; i<Ngroups_p; i++)  {
		XSgr_p[i]=CalcXSeff(EGroup[i], EGroup[i+1], errCorrXSgr_p[i], errUncorrXSgr_p[i], PHIgr_p[i], errPHIgr_p[i]);
		errXSgr_p[i] = sqrt(pow(errCorrXSgr_p[i],2)+pow(errUncorrXSgr_p[i],2));
		
		GroupFrac_p[i] = XSgr_p[i]*PHIgr_p[i];
		integral += GroupFrac_p[i];
	}
	
	for (int i=0; i<Ngroups_p; i++)  GroupFrac_p[i] /= integral;
	
	XSeff_p = CalcXSeff(EGroup[0], EGroup[Ngroups_p], ErrCorr_p, ErrUncorr_p, PHI_p, ErrPHI_p);
	
	return;
}




