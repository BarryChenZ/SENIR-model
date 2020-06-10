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
#include "ExperimentRecord.h"

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
int iter_times = 50;
int time_stamp = 10;
int Montecarlo_times = 1;
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

//SEIRD for single infected_rate 
vector<vector<double>> Montecarlo(vector<Node> NODES, vector<Gauss_Markov> GM, Physical_network P, Social_network S, vector<vector<double>> record) {
	vector<vector<double>> res;
	res.resize(iter_times, vector<double>(6, 0.0));
	int t = 0;
	//Setting initial infected node
	cout << 0 << endl;
	initial_infected(NODES, NODES.size());
	record_all_state(NODES, res, NODES.size(), t);
	t = 1;
	while (t < iter_times) {
		// Movement
		for (int j = 0; j < NODES.size(); j++) {
			GM[j].Set_d(), GM[j].Set_d();
			GM[j].Process(NODES);
		}
		//Calculating neighbor of P and S(no change)
		for (int j = 0; j < NODES.size(); j++) {
			P.get_degree(j, NODES);//short range
		}
		//Change probability 有點刻意
		for (int j = 0; j < NODES.size(); j++) {
			//4 value set by average
			NODES[j].exposed_rate = record[t][0];
			NODES[j].iNsidious_rate = record[t][1];
			NODES[j].infected_rate = record[t][2];
			NODES[j].infected_rate_2 = record[t][3];
			
		}
		cout << record[t][0] << " " << record[t][1] << endl;
		//Calculating state change
		vector<int> temp(NODES.size(), -1);
		exposedProcess(NODES, temp, S.adj_matrix, num);
		//iNsidiousProcess(NODES, temp, P, num);
		infectedProcess_1(NODES, temp, num);
		//infectedProcess_2(NODES, temp, num);
		recoveredProcess(NODES, temp, num);
		//wakeupProcess(NODES, temp, num);
		deathProcess(NODES, temp, num);
		loseImProcess(NODES, temp, num);
		unchanged(NODES, temp, num);
		//Node leaving and new node coming
		Leaving_and_Coming(temp, num);
		//Change
		update_state(NODES, temp, NODES.size());
		//Store
		cout << t << endl;
		record_all_state(NODES, res, NODES.size(), t);
		t += interval;
	}
	return res;
};
//Node, Mobility model, networks
vector<vector<double>> Montecarlo_modified(vector<Node> NODES, vector<Gauss_Markov> GM, Physical_network P, Social_network S) {//return the fraction of state along with time, run 100 times
	vector<vector<double>> res;
	res.resize(iter_times, vector<double>(6, 0.0));
	for (int i = 0; i < Montecarlo_times; i++) {
		int t = 0;
		//initial infection
		cout << 0 << endl;
		initial_infected(NODES, NODES.size());
		record_all_state(NODES, res, NODES.size(), t);
		t = 1;
		while (t < iter_times) {
			//Movement
			for (int j = 0; j < NODES.size(); j++) {
				GM[j].Set_d(), GM[j].Set_d();
				GM[j].Process(NODES);
			}
			//Calculating neighbor of P and S(no change)
			for (int j = 0; j < NODES.size(); j++) {
				P.get_degree(j, NODES);//short range
			}
			//Calculating state change
			vector<int> temp(NODES.size(), 0);
			Judgement(temp,NODES,S.adj_matrix, P, NODES.size());
			//Node leaving and new node coming
			Leaving_and_Coming(temp, num);
			//Change
			update_state(NODES, temp, NODES.size());
			//Store
			cout << t << endl;
			record_all_state(NODES, res, NODES.size(), t);
			t += interval;
		}

	}
	return res;
}
void record_result_to_file(vector<double>& result, int type) {
	fstream file;
	if (type == 1)file.open("C:\\Users\\User\\Desktop\\Matlab Code\\MonteCarlo_Simulation\\sim_result.txt", ios::out);
	else         file.open("C:\\Users\\User\\Desktop\\Matlab Code\\MonteCarlo_Simulation\\sim_result_2.txt", ios::out);
	for (int i = 0; i < result.size(); i++) {
		file << result[i] << endl;
	}
};
void record_result_and_print(vector<vector<double>> res, int type) {
	fstream file;
	//file.open("C:\\Users\\User\\Desktop\\Code\\Matlab Code\\SENIRD result\\SimRes.txt",ios::out);
	if (type == 1) {
		file.open("C:\\Users\\User\\Documents\\Github  liberary\\Matlab\\Result data\\SimRes_SEIRD_1_exp1.txt", ios::out);
	}//Modified
	else if (type == 2) {
		file.open("C:\\Users\\User\\Documents\\Github  liberary\\Matlab\\Result data\\SimRes_SEIRD_2_exp1.txt", ios::out);
	}//origin
	else if (type == 3) {
		file.open("C:\\Users\\User\\Documents\\Github  liberary\\Matlab\\Result data\\AnaRes_SEIRD_exp1.txt", ios::out);
	
	}//Anal
	for (int i = 0; i < res.size(); i++) {
		for (int j = 0; j < res[0].size(); j++) {
			file << res[i][j] << " ";
			//cout << res[i][j] << " ";
		}
		file << endl;
		//cout << endl;
	}
}


//Let the parameters of simulation and analytical model be the same.
void setting(Node i) {
	omega = i.wake_up_rate;
	//received_rate = i.received_rate;
	lambda = i.lose_immunity_rate;
	delta = i.death_rate;
	ex_delta = i.ex_death_rate;
}

//test network
void test() {
	srand((unsigned)time(NULL));
	vector<Node> NODES;
	vector<Gauss_Markov> GM;
	
	for (int i = 0; i < num; i++) NODES.push_back(Node(i)); //ID: 0 to n-1
	for (int i = 0; i < num; i++) GM.push_back(Gauss_Markov(i)); // assign initial speed
	Social_network S(num);
	Physical_network P(num, NODES);
	//vector<vector<double>> res1;
	//cout << "First:" << endl;
	//res1 = Montecarlo_modified(NODES, GM, P, S);
	//record_result_and_print(res1, 1);
	//vector<vector<double>> res2;
	//cout << "Second:" << endl;
	//res2 = Montecarlo(NODES, GM, P, S);
	//record_result_and_print(res2, 2);
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

	for (int k = 2; k <= 3; k++) {
		cout << "Third: Analytical model" << endl;
		vector<vector<double>> record;
		vector<vector<double>> res3;
		vector<Node_a> NODES_A;
		res3 = process_a(NODES_A, P, NODES, k, record);
		Exp1_Record(res3, 2, k);

		//第一種模擬方式
		/*
		vector<vector<double>> res2(iter_times, vector<double>(6, 0.0));
		for (int i = 0; i < time_stamp; i++) {
			//cout << "i: " << i << endl;
			vector<vector<double>> tmp;
			tmp = Montecarlo_modified(NODES, GM, P, S);
			//store
			for (int j = 0; j < tmp.size(); j++) {
				for (int q = 0; q < tmp[j].size(); q++) {
					res2[j][q] += tmp[j][q] / (double)time_stamp;
				}
			}
		}
		Exp1_Record(res2, 1, k);
		*/
		//第二種
		
		vector<vector<double>> res2(iter_times, vector<double>(6, 0.0));
		for (int i = 0; i < time_stamp; i++) {
			//cout << "i: " << i << endl;
			vector<vector<double>> tmp;
			tmp = Montecarlo(NODES, GM, P, S, record);
			//store
			for (int j = 0; j < tmp.size(); j++) {
				for (int q = 0; q < tmp[j].size(); q++) {
					res2[j][q] += tmp[j][q] / (double)time_stamp;
				}
			}
		}
		Exp1_Record(res2, 1, k);
		
	}
	

	//第二種模擬方式
	//vector<vector<double>> res2;
	//cout << "Second:" << endl;
	//res2 = Montecarlo(NODES, GM, P, S);
	//record_result_and_print(res2, 2);

	//Threshold
	//bool res = Calculating_Threshold(NODES_A[0]);
	//cout << res << endl;
	return;
}

