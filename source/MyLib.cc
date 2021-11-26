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
#include <algorithm> 

using namespace std;

bool CheckFileExists (string filename)  {
	ifstream in(filename.c_str());
	if (in.good()==0) 
	{
		cout << "Error! The file "<<filename<<" is not in this folder!\n";
		cin.get();
		return 0;
	}
	else return 1;
}

string RemoveExtension (string filename)  {
	int pos = filename.find_last_of('.');
	if (pos>0) filename.erase (filename.begin()+pos, filename.end()); 
	return filename;
}

string GetExtension (string filename)  {
	int pos = filename.find_last_of('.');
	if (pos>0) filename.erase (filename.begin(), filename.begin()+pos); 
	return filename;
}

string RemovePath (string filename)  {
	int pos = filename.find_last_of('/');
	if (pos>0) filename.erase (filename.begin(), filename.begin()+pos+1); 
	return filename;
}

string EraseSubstr (string original, string substr)  {
	int length = substr.length();
	int pos = original.find(substr);
	if (pos>=0) 
	{
		//cout<<"found at pos "<<pos<<" with legth "<<length<<endl;
		original.erase (original.begin()+pos, original.begin()+pos+length); 
	}
	return original;
}

string Int2String (int number)  {
	ostringstream convert;
	convert<<number;
	string word=convert.str();
	return word;
}

string Double2String (double number)
{
    ostringstream convert;
    convert<<number;
    string word=convert.str();
    return word;
}

int String2Int (string word)
{
    return atoi (word.c_str());
}

float String2Float (string word)
{
    return atof (word.c_str());
}

double String2Double (string word)
{
    char *ptr;
    return strtod(word.c_str(), &ptr);
}

string DelBlanks(string input)
{
  int i=0;
  while (i<input.size())
    {
      if (input[0]==' ') input.erase(input.begin()+i);
      if (input[i]==' ' && input[i-1]==' ') input.erase(input.begin()+i);
      else i++;
    }
  return input;
}


vector<string> Line2VecString (string &line)
{
      vector<string> VecString;
      istringstream iss(line);
      string word;
      while (iss >> word)
      {
        VecString.push_back(word);
      }
      return VecString;
}

double WeightedAverage (vector<double> &x, vector<double> &err, double &errWghtAvg)  {
	double mean=0, w, sumw=0;
	if(x.size()!=err.size()) {
		cout<<"Error: x.size()!= err.size()"<< endl;
		cin.get();
		return 1.;
	}
	
	for (int i=0; i<x.size(); i++) {
		if (err[i]<=0) {
			cout<<"Error: err["<<i<<"] = "<<err[i]<< endl;
			cin.get();
			return 1.;
		}
		
		w = pow(err[i], -2);
		mean += w*x[i];
		sumw += w;
	}
	mean = mean/sumw;
	errWghtAvg = sqrt(1./sumw);
	return mean;
}

string Num2String (double number)
{
	ostringstream convert;
	convert<<scientific<<setprecision(2)<<number;
	string word=convert.str();
	return word;
}

bool YesNo(string Question)
{
	char answer;
	string line;
	cout << Question;
	cin.ignore(1,'\n');
	getline (cin, line);
	if (line.length()==0) answer='N';
	else 	answer = line[0];
	if (answer == 'y') return true;
	return false;
}


