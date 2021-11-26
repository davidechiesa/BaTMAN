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

#ifndef mylib_h
#define mylib_h

#include <string>
#include <vector>
using namespace std;


string RemoveExtension (string filename);
string GetExtension (string filename);
bool CheckFileExists (string filename);
string RemovePath (string filename);
string EraseSubstr (string original, string substr);

string Num2String (double number);
string Int2String (int number);
int String2Int (string word);
float String2Float (string word);
double String2Double (string word);
string Double2String (double number);
string DelBlanks(string input);
vector<string> Line2VecString (string &line);

double WeightedAverage (vector<double> &x, vector<double> &err, double &errWghtAvg);

bool YesNo(string Question);

#endif
