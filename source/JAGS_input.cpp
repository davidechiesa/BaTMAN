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
#include <iomanip>
#include <math.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include "TGraph.h"
#include "TAxis.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"

#include "NucReaction.h"
#include "MyLib.h"

#define FluxNorm 1

using namespace std;

void ReadInputData (string filename, vector<string> &ReactionName, vector<double> &obsR, vector<double> &sigmaR, vector<string> &SelfShield);
vector<double> FillEgroups (string filename);

int main (int numArg, char *listArg[])
{
    cout << "\nThis program is part of the BaTMAN software package\n"
         << "(Bayesian-unfolding Toolkit for Multi-foil Activation with Neutrons)\n"
         << "Copyright (C) 2021 Davide Chiesa (University and INFN of Milano - Bicocca)\n"
         << "This program comes with ABSOLUTELY NO WARRANTY\n"<<endl;
    cout << "\tPress any key to START\n";
    cin.get();

	if (numArg!=6)	{
		cout<<"\nTo launch program: JAGS_input InputData  FluxSpectrum  XSpath  EnergyGroups SelfShieldingPATH"<<endl;
    		return 1;
 	}
	string InputData=listArg[1];
	string Guess_Flux=listArg[2];
	string XSpath=listArg[3];
	string En_group=listArg[4];
	string pathSF=listArg[5];
	
	string OutputName=RemovePath(RemoveExtension(InputData));
	
	ofstream out((OutputName+"_XSgr.txt").c_str());	
	for (int i=0; i<numArg; i++) out << listArg[i]<<" "; 
	out<<endl<<endl;
  		
	double fatt_norm=FluxNorm; 
	
  	vector<string> ReactionName;
	vector<double> obsR; 
	vector<double> sigmaR;	
	vector<string> SelfShield;
	
	ReadInputData (InputData, ReactionName, obsR, sigmaR, SelfShield);
	int Nreactions=ReactionName.size();
	
	
	vector<double> EGroups = FillEgroups (En_group);
	int Ngr = EGroups.size()-1;
	
	
	// Reading cross section data and calculate effective cross sections
	string XSfile;

  	double* PHIgr;
  	double* ErrPHIgr;
  	double FluxTot=0., ErrFluxTot;    	
  	
	vector <double*> XSgr;
	vector <double*> ErrCorr;
	vector <double*> GroupFrac;
	
	vector<double> XSeff;
	
	NucReaction* a;
		
	for (int j=0; j<Nreactions; j++)  {
  		XSfile=XSpath+ReactionName[j];
  		  		
  		a = new NucReaction (ReactionName[j]);
  		
  		a->FillXSwithErr(XSfile.c_str());
  		a->FillFlux(Guess_Flux.c_str(),fatt_norm);  		
  		
  		if (SelfShield[j].compare("NoSelfShield")!=0) a->SetSelfShielding (pathSF+SelfShield[j]);		
  		
  		a->CalcXSgroups (EGroups);	
  		 
  		XSgr.push_back(a->GetXSgr());
  		ErrCorr.push_back(a->GetErrCorrXSgr());
  		GroupFrac.push_back(a->GetGroupFraction());
  		XSeff.push_back(a->GetXSeff());
  		
  		if (j==0)  {
  			PHIgr=a->GetFluxGr();
  			ErrPHIgr=a->GetErrFluxGr();
  			FluxTot=a->GetPHI();  	
  			ErrFluxTot=a->GetErrPHI(); 		
  		}
 	}
  	
  	cout<<"\nSize vector XSgr: "<<XSgr.size()<<endl;
  	cout.precision(2);
  	out.precision(2);
  	
  	cout<<"\nGuess Flux"<<endl<<"E_low[MeV]\tE_up[MeV]\tPHI[n/cm2/s]\tErr[n/cm2/s]\terr%\tGroup(%)\n";
  	out<<"\nGuess Flux"<<endl<<"E_low[MeV]\tE_up[MeV]\tPHI[n/cm2/s]\tErr[n/cm2/s]\terr%\tGroup(%)\n";
	for (int i=0; i<Ngr; i++)  {
		cout<<scientific<<EGroups[i]<<"\t"<<EGroups[i+1]<<"\t"<<PHIgr[i]<<"\t"<< ErrPHIgr[i] <<"\t";
		cout<<fixed<< 100*ErrPHIgr[i]/PHIgr[i] <<"%\t"<<100*PHIgr[i]/FluxTot<<"%"<<endl;
		out<<scientific<<EGroups[i]<<"\t"<<EGroups[i+1]<<"\t"<<PHIgr[i]<<"\t"<< ErrPHIgr[i] <<"\t";
		out<<fixed<< 100*ErrPHIgr[i]/PHIgr[i] <<"%\t"<<100*PHIgr[i]/FluxTot<<"%"<<endl;
	}
	cout << "Total"<<endl;
	cout<<scientific<<EGroups[0]<<"\t"<<EGroups[Ngr]<<"\t"<<FluxTot<<"\t"<< ErrFluxTot <<"\t";
	cout<<fixed<< 100*ErrFluxTot/FluxTot <<"%\n";
	out << "Total"<<endl;
	out<<scientific<<EGroups[0]<<"\t"<<EGroups[Ngr]<<"\t"<<FluxTot<<"\t"<< ErrFluxTot <<"\t";
	out<<fixed<< 100*ErrFluxTot/FluxTot <<"%\n";
	
	cout<<"\nEffective cross section (global)"<<endl;
	out<<"\nEffective cross section (global)"<<endl;
	for ( int i=0; i<Nreactions; i++ )    	{	
	    	  cout<<left<<setw(20)<<ReactionName[i]<<"  "<<scientific<<XSeff[i]<<"\n"; //<<err_xs_eff_tot[i]<<endl;	
	    	  out<<left<<setw(20)<<ReactionName[i]<<"  "<<scientific<<XSeff[i]<<"\n"; //<<err_xs_eff_tot[i]<<endl;	
	}
		
	cout<<"\nContribution percentage to total activation by neutron flux groups"<<endl;
	out<<"\nContribution percentage to total activation by neutron flux groups"<<endl;
	for ( int j=0; j<Nreactions; j++ )  {
		cout <<left<<setw(20)<<ReactionName[j]<< "  ";
		for (int i=0; i<Ngr; i++) cout<<fixed<<"\t"<<100*GroupFrac[j][i]<<"%";
		cout<<endl;	
		out <<left<<setw(20)<<ReactionName[j]<< "  ";
		for (int i=0; i<Ngr; i++) out<<fixed<<"\t"<<100*GroupFrac[j][i]<<"%";
		out<<endl;
	}
	
	cout<<"\nGroup effective cross sections"<<endl;
	out<<"\nGroup effective cross sections"<<endl;
	for ( int j=0; j<Nreactions; j++ ) {
		cout <<left<<setw(20)<<ReactionName[j]<< "  ";
		for (int i=0; i<Ngr; i++) cout<<scientific<< XSgr[j][i] <<"  ";	
		cout<<endl;
		out <<left<<setw(20)<<ReactionName[j]<< "  ";
		for (int i=0; i<Ngr; i++) out<<scientific<< XSgr[j][i] <<"  ";	
		out<<endl;	
	}
	
	cout<<"\nGroup effective cross section uncertainties"<<endl;
	out<<"\nGroup effective cross section uncertainties"<<endl;
	for ( int j=0; j<Nreactions; j++ ) {
		cout <<left<<setw(20)<<ReactionName[j]<< "  ";
		for (int i=0; i<Ngr; i++) cout<<scientific<< ErrCorr[j][i] <<"  ";	
		cout<<endl;
		out <<left<<setw(20)<<ReactionName[j]<< "  ";
		for (int i=0; i<Ngr; i++) out<<scientific<< ErrCorr[j][i] <<"  ";	
		out<<endl;	
	}
	out.close();

	
	///////////////////////////////////////////
	// JAGS input file creation
	///////////////////////////////////////////
	
	string Question = "Do you want cross section uncertainties to be included in the unfolding? (y/N) ";
	bool answer = YesNo(Question);
	
	out.open("Data4Jags.dat");
	out.setf(ios_base::scientific);
	out.precision(4);
  
	out<<"Ngr <- "<<Ngr<<endl;
	out << endl;
	for (int i=0; i<Nreactions; i++) {
		out<<"obsR"<<i+1<<" <- "<<obsR[i]<<endl;
		out<<"sigmaR"<<i+1<<" <- "<<sigmaR[i]<<endl;
		if (answer==true) out<<"XSmu"<<i+1<<" <-c( ";
		else out<<"XS"<<i+1<<" <-c( ";
		for (int j=0; j<Ngr; j++)  {
			if (j<(Ngr-1)) out<< XSgr[i][j] <<", ";
			else out<< XSgr[i][j] <<")\n";
		}
		if (answer==true) {
		    out<<"ErrCorr"<<i+1<<" <-c( ";
		    for (int j=0; j<Ngr; j++)  {
			    if (j<(Ngr-1)) out<< ErrCorr[i][j] <<", ";
			    else out<< ErrCorr[i][j] <<")\n";
		    }
		}
		out << endl;
	}
	out.close(); 
	
	out.open("Model4Jags.bug");
	out.setf(ios_base::scientific);
	out.precision(4);
  	out<<"model {"<<endl;
	for (int i=0; i<Nreactions; i++) {
	    if (answer==true) {
		    out<<"\trandom"<<i+1<<"~dnorm(0,1)"<<endl;
		    out<<"\tXS"<<i+1<<"<- XSmu"<<i+1<<" + random"<<i+1<<"*ErrCorr"<<i+1<<endl; 
		}
		out<<"\tR"<<i+1<<"<- inprod(XS"<<i+1<<",phi)"<<endl; 
		out<<"\tobsR"<<i+1<<" ~ dnorm ( R"<<i+1<<", pow(sigmaR"<<i+1<<",-2))"<<endl;
		out << endl;
	}
	out<<"\tphiTOT<-sum(phi)"<<endl;
	out<<"\tfor (i in 1:Ngr) {"<<endl;
	out<<"\t\tphi[i]~dunif(0,10)"<<endl;
	out<<"\t}"<<endl;
	out<<"}"<<endl;
	out.close();
	
	out.open("LaunchJags.cmd");
	out<<"data in Data4Jags.dat"<<endl;
	out<<"model in Model4Jags.bug"<<endl;
	out<<"compile, nchains(4)"<<endl;
	out<<"initialize"<<endl;
	out<<"update 25000"<<endl;
	out<<"monitor set phi, thin(10)"<<endl;
	out<<"monitor set phiTOT, thin(10)"<<endl;
	for (int i=0; i<Nreactions; i++) out<<"monitor set R"<<i+1<<", thin(10)"<<endl;
	if (answer==true) 
	    for (int i=0; i<Nreactions; i++) 
	        out<<"monitor set random"<<i+1<<", thin(10)"<<endl;
	out<<"update 100000, by(1000)"<<endl;
	out<<"coda *"<<endl;
	out<<"data to useddata"<<endl;
	out.close();
	
	
	cout << "\nOutput written to files:"<< endl;
	cout << "--> "<<OutputName+"_XSgr.txt"<< endl;
	cout << "--> Data4Jags.dat"<< endl;
	cout << "--> Model4Jags.bug"<< endl;
	cout << "--> LaunchJags.cmd"<< endl;
	 
	return 0;
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

void ReadInputData (string filename, vector<string> &ReactionName, vector<double> &obsR, vector<double> &sigmaR, vector<string> &SelfShield)  {
	CheckFileExists(filename);
	double AM,IA,ASS,err;
	double R,errR; // variables for the activation rates divided by the nubler of target isotopes (R/N)
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
		ReactionName.push_back(reaction);

		R=(ASS*AM)/(IA*N_avog);
		errR=R*err/ASS;

		obsR.push_back(R);
		sigmaR.push_back(errR);
		
		SelfShield.push_back(SFfile);

		cout.width(15);
		cout<<reaction;
		cout.width(15);
		cout<<fixed<<AM;
		cout.width(15);
		cout<<IA;
		cout.width(15);
		cout<<scientific<<ASS;
		cout.width(15);
		cout<<err;
		cout.width(15);
		cout<<scientific<<R;
		cout.width(15);
		cout<<errR;
		cout.width(20);
		cout<<SFfile<<endl;

		Nreactions++;
	}
	in.close();
	cout<<"Number of activation reactions: "<<Nreactions<<endl;

	double FluxUnit;
	cout<<"Set order of magnitude of total neutron flux (e.g. 1e12): ";
	cin>>FluxUnit;
	double factor = 1e24/FluxUnit;
	for (int i=0; i<obsR.size(); i++) {
		obsR[i] *= factor;
		sigmaR[i] *= factor;
	}	
	
	return;
}

