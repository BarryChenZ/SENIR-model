// Final_paper_simulation.cpp : 此檔案包含 'main' 函式。程式會於該處開始執行及結束執行。
// Main parｔcall header
#include "pch.h"
#include <iostream>
#include <random>
#include <cmath>
#include <time.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <vector>"
#include <omp.h>
#include "ProcessModel.h"
#include "OtherE_model.h"
#include "Distribution.h"
#include "Mobility.h"

using namespace std;
/*
double infected_rate = 0.32;
double received_rate = 0.228;
double iNsidious_rate = 0.8;
double recovered_rate = 0.083; 
double infected_sim_array[4] = { 0.2, 0.3, 0.4, 0.5 };
double bornRate = 5;
double DeathRate = 1;
*/

int n = 4039;
int iter_times = 100;
int time_stamp = 12;
int interval = 1; // can run 3 times a day or one time a day 

vector<double> Montecarlo();
void record_result_to_file(vector<vector<double>>& result, int type);

int main()
{	
	srand((unsigned)time(NULL));
	/*
	int** adj_matrix = new int*[n];
	for (int i = 0; i < n; i++) {
		adj_matrix[i] = new int[n];
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			adj_matrix[i][j] = 0;
		}
	}
	vector<int> id;
	int tmp;
	ifstream input("C:\\Users\\User\\Desktop\\Dataset\\facebook_combined.txt", ios::in);
	while (input >> tmp) {
		id.push_back(tmp);
	}
	for (int i = 0; i < id.size(); i = i + 2) {
		int a = id[i];
		int b = id[i + 1];
		adj_matrix[a][b] = 1;
		adj_matrix[b][a] = 1;
	}
	input.close();
	*/

	//MonteCarlo Simulation
	//call node and network(S,P) -> distribution node on the network -> doing simulation -> record result
	vector<Node> NODES;
	for (int i = 0; i < n; i++) NODES.push_back(Node(i));
	Social_network S(n);//BA model
	Physical_network P(n, NODES);//only uniform distribution, 詳細的走法在mobility.h
	
	Levy_Walk LEVY(n);

	vector<double> record = Montecarlo(NODES, P, S);//Running simulation
	
	//vector<double> Result;

	record_result_to_file(record, 2);
	
	return 0;
	
}

//SENIRD for single infected_rate
vector<double> Montecarlo(vector<Node>& NODES, Physical_network& P, Social_network& S) {
	vector<double> ans;
	double t = 0;
	//Setting initial infected node

	while (t < iter_times) {
		vector<int> temp(n, 0);
		//依據前一step的state來執行
		exposedProcess(NODES, temp, S.adj_matrix);
		iNsidiousProcess(NODES, temp, P);
		infectedProcess_1(NODES, temp);
		infectedProcess_2(NODES, temp);
		recoveredProcess(NODES, temp);

		wakeupProcess(NODES, temp);
		deathProcess(NODES, temp);
		
		update_state(NODES, temp);
		record(NODES, ans);
		t += interval;
	}
	return ans;
};

void record_result_to_file(vector<double>& result, int type) {
	fstream file;
	if (type == 1)file.open("C:\\Users\\User\\Desktop\\Matlab Code\\MonteCarlo_Simulation\\sim_result.txt", ios::out);
	else         file.open("C:\\Users\\User\\Desktop\\Matlab Code\\MonteCarlo_Simulation\\sim_result_2.txt", ios::out);
	for (int i = 0; i < result.size(); i++) {
		file << result[i] << endl;
	}
};
