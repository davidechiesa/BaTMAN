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
#include <fstream>
#include <ctime>
#include <sstream>
#include "date.h"

using namespace std;

void PrintTM (time_t rawtime) {
	tm *date = localtime(&rawtime);

	cout << "\tdate.tm_sec = "<<date->tm_sec<<"\n";
	cout << "\tdate.tm_min = "<<date->tm_min<<"\n";
	cout << "\tdate.tm_hour = "<<date->tm_hour<<"\n";
	cout << "\tdate.tm_mday = "<<date->tm_mday<<"\n";
	cout << "\tdate.tm_mon = "<<date->tm_mon<<"\n";
	cout << "\tdate.tm_year = "<<date->tm_year<<"\n";
	cout << "\tdate.tm_wday = "<<date->tm_wday<<"\n";
	cout << "\tdate.tm_yday = "<<date->tm_yday<<"\n";
	cout << "\tdate.tm_isdst = "<<date->tm_isdst<<"\n";
	cout << "\trawtime = "<<rawtime<<"\n";
	return;
}

void PrintDATE (time_t rawtime) {

	tm *date = localtime(&rawtime);

	//cout << "Date: ";
	cout <<setfill('0')<<setw(2)<< date->tm_mday<<"/";
	cout <<setfill('0')<<setw(2)<< date->tm_mon+1 << "/" << date->tm_year +1900 << " ";
	if (date->tm_isdst == 1) {
		cout <<setfill('0')<<setw(2) << date->tm_hour-1 <<":";
		cout <<setfill('0')<<setw(2) << date->tm_min << ":";
		cout <<setfill('0')<<setw(2) << date->tm_sec << "\n";
	}
	else {
		cout <<setfill('0')<<setw(2) << date->tm_hour <<":";
		cout <<setfill('0')<<setw(2) << date->tm_min << ":";
		cout <<setfill('0')<<setw(2) << date->tm_sec << "\n";
	}

	return;
}


time_t CreateDATE( string date, string time )
{
	int* d=ReadDate(date, "/");
	int* t=ReadDate(time, ":");
		
	int Day = d[0];
	int Month = d[1];
	int Year = d[2];
	int Hour = t[0];
	int Minutes = t[1];
	int Seconds = t[2];	

	if (Month>12) {
		cout << "Error: Month = "<< Month <<"\tPress ANY key to continue"<<endl;
		cin.get();
		return 0;
	}
	
	bool gg30=0;
	bool gg31=0;
	bool feb=0;
	if (Month==11 || Month==4 || Month==6 || Month==9) gg30=true;
	else if (Month==2) feb=true;
	else gg31=true;
	
	if (gg30==true && Day>30) {
		cout << "Error: Month = "<< Month <<"\tDay = "<< Day <<"\tPress ANY key to continue"<<endl;
		cin.get();
		return 0;
	}
	
	if (gg31==true && Day>31) {
		cout << "Error: Month = "<< Month <<"\tDay = "<< Day <<"\tPress ANY key to continue"<<endl;
		cin.get();
		return 0;
	}
	
	if (feb==true && Day>29) {
		cout << "Error: Month = "<< Month <<"\tDay = "<< Day <<"\tPress ANY key to continue"<<endl;
		cin.get();
		return 0;
	}
	
	if (Hour>23) {
		cout << "Error: Hour = "<< Hour <<"\tPress ANY key to continue"<<endl;
		cin.get();
		return 0;
	}
	
	if (Minutes>59) {
		cout << "Error: Minutes = "<< Minutes <<"\tPress ANY key to continue"<<endl;
		cin.get();
		return 0;
	}
	
	if (Seconds>59) {
		cout << "Error: Seconds = "<< Seconds <<"\tPress ANY key to continue"<<endl;
		cin.get();
		return 0;
	}

	tm data; 
	data.tm_sec = Seconds;
	data.tm_min = Minutes;
	data.tm_hour = Hour;
	data.tm_mday = Day;
	data.tm_mon = Month - 1;
	data.tm_year = Year - 1900;
	data.tm_isdst = 0;

	time_t rawtime =  mktime ( &data );
	/*cout << rawtime << endl;
	PrintDATE (rawtime );
	PrintTM (rawtime); 
	rawtime =  mktime ( &data );
	cout << rawtime << endl;
	PrintDATE (rawtime );
	PrintTM (rawtime); */

	return rawtime;
}

int* ReadDate (string date, string separator){
	for (int i=0; i<2; i++) {
		size_t f = date.find(separator);
		if (f==date.npos) cout <<"Error! No matches found for separator '"<<separator<<"' in string '"<<date<<"'"<< endl;
		date.replace(f, std::string(separator).length(), " ");
	}
	//cout << date << endl;
	stringstream ss;
	ss<<date;
	
	int *d = new int (3);
	ss>>d[0]>>d[1]>>d[2];
	return d;
}

