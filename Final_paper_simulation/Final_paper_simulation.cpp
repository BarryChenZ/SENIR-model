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
#include <vector>
#include <omp.h>
#include "ProcessModel.h"
#include "OtherE_model.h"
#include "Distribution.h"
#include "Mobility.h"
#include "FinalModel.h"

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

int num = 1000;//4039;
int iter_times = 10;
int time_stamp = 12;
int interval = 1; // can run 3 times a day or one time a day 

vector<double> Montecarlo();
void record_result_to_file(vector<vector<double>>& result, int type);
void test();

int main()
{	
	//srand((unsigned)time(NULL));
	test();
	return 0;
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
	*/
}

//SENIRD for single infected_rate
vector<double> Montecarlo(vector<Node>& NODES, Physical_network& P, Social_network& S) {
	vector<double> ans;
	double t = 0;
	//Setting initial infected node
	while (t < iter_times) {
		vector<int> temp(num, 0);
		//依據前一step的state來執行
		exposedProcess(NODES, temp, S.adj_matrix, num);
		iNsidiousProcess(NODES, temp, P, num);
		infectedProcess_1(NODES, temp, num);
		infectedProcess_2(NODES, temp, num);
		recoveredProcess(NODES, temp, num);

		wakeupProcess(NODES, temp, num);
		deathProcess(NODES, temp, num);

		update_state(NODES, temp, num);
		record(NODES, ans, num);
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

//test network
void test() {
	srand((unsigned)time(NULL));
	vector<Node> NODES;
	for (int i = 0; i < num; i++) NODES.push_back(Node(i)); //ID: 0 to n-1
	Social_network S(num);
	Physical_network P(num, NODES);
	/*
	for (int i = 0; i < num; i++) {
		cout << i << ":" << NODES[i].get_position_x() << " " << NODES[i].get_position_y() << endl;
	}
	//position test ok
	*/


	//for (int i = 0; i < num; i++) {
		//cout << i << ":" << S.get_degree(i) << " " <<  P.get_degree(i, NODES) << endl;
	//}
	
	
	//Mobility test
	//Orient ok
	/*
	vector<Identity> Pos;
	for (int i = 0; i < num; i++) Pos.push_back(Identity(i,30)); // 30個位置
	int t = 0;
	while (t < 10) {
		for (int i = 0; i < num; i++) {
			Pos[i].Process();
		}
		cout << Pos[3].At_position << endl;
		t++;
	}
	*/
	//Gauss-Markov moblity model ok
	/*
	vector<Gauss_Markov> GM_model;
	for (int i = 0; i < num; i++) GM_model.push_back(Gauss_Markov(i));
	cout << GM_model[5].get_velocity() << endl;
	for (int i = 0; i < 40; i++) {
		GM_model[5].Set_v();
		cout << GM_model[5].get_velocity() << endl;
	}
	
	int t = 0;
	while (t < 10) {
		for (int i = 0; i < num; i++) {
			GM_model[i].Set_v(), GM_model[i].Set_d();
			GM_model[i].Process(NODES);
		}
		cout << NODES[3].get_position_x() << " " << NODES[3].get_position_y() << endl;
		t++;
	}
	
	*/

	/*
	//ProcessModel.h test for no mobility node
	vector<double> ans_I; // record one state
	vector<vector< double>> ans_all(iter_times, vector<double>(0,0)); // 6 state
	double t = 0;
	//Setting initial infected node
	initial_infected(NODES, num);
	while (t < iter_times) {
		vector<int> temp(num, -1);//-1 
		//依據前一step的state來執行
		
		exposedProcess(NODES, temp, S.adj_matrix, num);
		iNsidiousProcess(NODES, temp, P, num);
		infectedProcess_1(NODES, temp, num);
		infectedProcess_2(NODES, temp, num);
		recoveredProcess(NODES, temp, num);
		loseImProcess(NODES, temp, num);

		wakeupProcess(NODES, temp, num);
		deathProcess(NODES, temp, num);
		unchanged(NODES, temp, num)
		
		Judgement(temp, NODES, S.adj_matrix, P, num);
		update_state(NODES, temp, num);
		record_all_state(NODES, ans_all, num, t);
		//record(NODES, ans_I, num);
		t += interval;
		
	}
	*/

	// Analytical model test
	vector<Node_a> NODES_A;
	process_a(NODES_A, P);
	return;
}
