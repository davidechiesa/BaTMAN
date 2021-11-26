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


#ifndef NucReaction_h
#define NucReaction_h

#include<TH1F.h>
#include<TGraphErrors.h>

using namespace std;

class NucReaction
{
	public:

	NucReaction (string name);
	~NucReaction ();	
	
	void FillXS (string filename);	
	void FillXSwithErr (string filename);	
	void FillFlux (string filename, double NormPhi);
	void SetXSErrors (double errRelTHERM, double errRelRES);
	void SetSelfShielding (string filename);
	
	double CalcXSeff(double min, double max, double& ErrCorr, double& ErrUncorr, double& Flux, double& ErrFlux);
	void CalcXSgroups (vector<double> &EGroups);
	
	string		GetName() {return name_p;}
	TH1F* 		GetFlux() {return Flux_p;}
	TGraphErrors* 	GetXS() {return grXS_p;}
	int		GetNgroups() {return Ngroups_p;}
	double 		GetFluxTOT() {return FluxTOT_p;}
	double 		GetNormPhi() {return NormPhi_p;}
	
	double 		GetXSeff() {return XSeff_p;}
	double 		GetErrCorr() {return ErrCorr_p;}
	double 		GetErrUncorr() {return ErrUncorr_p;}
	double 		GetPHI() {return PHI_p;}
	double 		GetErrPHI() {return ErrPHI_p;}
	
	double*		GetFluxGr() {return PHIgr_p;}
	double*		GetErrFluxGr() {return errPHIgr_p;}
	double*		GetXSgr() {return XSgr_p;}
	double*		GetErrXSgr() {return errXSgr_p;}
	double*		GetErrCorrXSgr() {return errCorrXSgr_p;}
	double*		GetErrUncorrXSgr() {return errUncorrXSgr_p;}
	double*		GetGroupFraction() {return GroupFrac_p;}
	
	private:
	
	string name_p;
	
	TGraphErrors* grXS_p;
	
	TH1F* Flux_p;
	TH1F* SelfShieding_p;

	double NormPhi_p;
	double FluxTOT_p;
	
	double XSeff_p;
	double ErrCorr_p;
	double ErrUncorr_p;
	double PHI_p;
	double ErrPHI_p;	
	
	int Ngroups_p;
	double* PHIgr_p;
	double* errPHIgr_p;
	double* XSgr_p;	
	double* errXSgr_p;
	double* errCorrXSgr_p;	
	double* errUncorrXSgr_p;
	double* GroupFrac_p;	
	
};

#endif
