#pragma once
//區隔不同實驗的儲存file
#include<iostream>
#include<vector>
#include<math.h>
#include<fstream>

using namespace std;

// type1 : 區隔simulation and analytical, type2 : 
void Exp1_Record(vector<vector<double>> res, int type1, int type2) {// 模擬或數學/第幾個值/第幾個實驗
	fstream file;
	string s;
	string e = ".csv";
	string value1 = to_string(type2);
	if (type1 == 1) {
		s = "C:\\Users\\User\\Documents\\Github  liberary\\Matlab\\Result data\\exp1\\SimRes_exp1_";
	}
	else if (type1 == 2) {
		s = "C:\\Users\\User\\Documents\\Github  liberary\\Matlab\\Result data\\exp1\\AnaRes_exp1_";
	}
	s = s + value1 + e;
	file.open(s, ios::out);
	for (int i = 0; i < res.size(); i++) {
		for (int j = 0; j < res[0].size(); j++) {
			file << res[i][j] << ",";
			//cout << res[i][j] << " ";
		}
		file << endl;
		//cout << endl;
	}
}

void Exp2_Record(vector<vector<double>> res, int type1, int type2) {// 模擬或數學/第幾個值/第幾個實驗
	fstream file;
	string s;
	string e = ".csv";
	string value1 = to_string(type2);
	if (type1 == 1) {
		s = "C:\\Users\\User\\Documents\\Github  liberary\\Matlab\\Result data\\exp2\\SimRes_exp2_";
	}
	else if (type1 == 2) {
		s = "C:\\Users\\User\\Documents\\Github  liberary\\Matlab\\Result data\\exp2\\AnaRes_exp2_";
	}
	s = s + value1 + e;
	file.open(s, ios::out);
	for (int i = 0; i < res.size(); i++) {
		for (int j = 0; j < res[0].size(); j++) {
			file << res[i][j] << ",";
			//cout << res[i][j] << " ";
		}
		file << endl;
		//cout << endl;
	}
}