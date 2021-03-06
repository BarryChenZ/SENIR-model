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
#include <unordered_map>
#include <omp.h>
#include "ProcessModel.h"
#include "OtherE_model.h"
#include "Distribution.h"
#include "Mobility.h"
#include "FinalModel.h"
#include "ExperimentRecord.h"
//2020 08 10 資料處理
#include "DataProcess.h" 
//平行化
#include <ppl.h>
#include <algorithm>

using namespace std;
using namespace concurrency;

/*
double infected_rate = 0.32;
double received_rate = 0.228;
double iNsidious_rate = 0.8;
double recovered_rate = 0.083; 
double infected_sim_array[4] = { 0.2, 0.3, 0.4, 0.5 };
double bornRate = 5;
double DeathRate = 1;
*/

int num = 1000;//1000,4039;
int num_array[3] = { 3000, 5000, 10000 };
int iter_times = 50;//50
//一種要跑1000次 100*1000次
int time_stamp = 100;
int ScatterTime = 1;
int Montecarlo_times = 1;
int interval = 1; // can run 3 times a day or one time a day 

//vector<double> Montecarlo();
void record_result_to_file(vector<vector<double>>& result, int type);
void test();

int main()
{	
	//srand((unsigned)time(NULL));
	test();
	return 0;
}

//SEIRD for single infected_rate 
vector<vector<double>> Montecarlo(vector<Node> NODES, vector<Gauss_Markov> GM, Physical_network P, Social_network S, unordered_map<int, vector<vector<double>>> record) {
	vector<vector<double>> res;
	res.resize(iter_times, vector<double>(6, 0.0));
	int t = 0;
	//Setting initial infected node
	//cout << 0 << endl;
	initial_infected(NODES, NODES.size());
	record_all_state(NODES, res, NODES.size(), t);
	t = 1;

	//vector<int> order;
	//for (int j = 0; j < 5; j++) order.push_back(j);
	//random_shuffle(order.begin(), order.end());
	while (t < iter_times) {
		// Movement
		for (int j = 0; j < NODES.size(); j++) {
			GM[j].Set_v(), GM[j].Set_d();
			GM[j].Process(NODES);
		}
		//Calculating neighbor of P and S(no change)
		for (int j = 0; j < NODES.size(); j++) {
			P.get_degree(j, NODES);//short range
		}
		//Change probability 有點刻意
		for (int j = 0; j < NODES.size(); j++) {
			//4 value set by average
			NODES[j].exposed_rate = record[NODES[j].degree_S][t][0];
			NODES[j].iNsidious_rate = record[NODES[j].degree_S][t][1];
			NODES[j].infected_rate = record[NODES[j].degree_S][t][2];
			NODES[j].infected_rate_2 = record[NODES[j].degree_S][t][3];
		}
		//Calculating state change 0E 1I 2R 3D 4L
		vector<int> temp(NODES.size(), -1);
		//打亂順序
		/*
		for (int j = 0; j < order.size(); j++) {
			int numOrder = order[j];
			cout << numOrder << endl;
			switch (numOrder) {
				case 0:
					exposedProcess(NODES, temp, S.adj_matrix, num);
					break;
				case 1:
					infectedProcess_1(NODES, temp, num);
					break;
				case 2:
					recoveredProcess(NODES, temp, num);
					break;
				case 3:
					deathProcess(NODES, temp, num);
					break;
				case 4:
					loseImProcess(NODES, temp, num);
					break;
			}
		*/

		recoveredProcess(NODES, temp, NODES.size());
		deathProcess(NODES, temp, NODES.size());
		loseImProcess(NODES, temp, NODES.size());
		exposedProcess(NODES, temp, S.adj_matrix, NODES.size());
		infectedProcess_1(NODES, temp, NODES.size());
		//recoveredProcess(NODES, temp, NODES.size());
		//deathProcess(NODES, temp, NODES.size());
		//loseImProcess(NODES, temp, NODES.size());
		
		
		unchanged(NODES, temp, NODES.size());
		//Node leaving and new node coming
		Leaving_and_Coming(NODES, temp, NODES.size());
		//Change
		update_state(NODES, temp, NODES.size());
		//Store
		//cout << t << endl;
		record_all_state(NODES, res, NODES.size(), t);
		t += interval;
	}
	return res;
};
//Node, Mobility model, networks
vector<vector<double>> Montecarlo_modified(vector<Node> NODES, vector<Gauss_Markov> GM, Physical_network P, Social_network S, unordered_map<int, vector<vector<double>>> record) {//return the fraction of state along with time, run 100 times
	vector<vector<double>> res;
	res.resize(iter_times, vector<double>(6, 0.0));
	for (int i = 0; i < Montecarlo_times; i++) {
		int t = 0;
		//initial infection
		//cout << 0 << endl;
		initial_infected(NODES, NODES.size());
		record_all_state(NODES, res, NODES.size(), t);
		t = 1;
		while (t < iter_times) {
			//Movement
			for (int j = 0; j < NODES.size(); j++) {
				GM[j].Set_v(), GM[j].Set_d();
				GM[j].Process(NODES);
			}
			//Calculating neighbor of P and S(no change)
			for (int j = 0; j < NODES.size(); j++) {
				P.get_degree(j, NODES);//short range
			}
			//Change probability 有點刻意
			for (int j = 0; j < NODES.size(); j++) {
				//4 value set by average
				NODES[j].exposed_rate = record[NODES[j].degree_S][t][0];
				NODES[j].iNsidious_rate = record[NODES[j].degree_S][t][1];
				NODES[j].infected_rate = record[NODES[j].degree_S][t][2];
				NODES[j].infected_rate_2 = record[NODES[j].degree_S][t][3];
			}
			//Calculating state change
			vector<int> temp(NODES.size(), 0);
			Judgement(temp,NODES,S.adj_matrix, P, NODES.size());

			unchanged(NODES, temp, num);
			//Node leaving and new node coming
			Leaving_and_Coming(NODES, temp, NODES.size());
			//Change
			update_state(NODES, temp, NODES.size());
			//Store
			//cout << t << endl;
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
/*
void setting(Node i) {
	omega = i.wake_up_rate;
	//received_rate = i.received_rate;
	lambda = i.lose_immunity_rate;
	delta = i.death_rate;
	ex_delta = i.ex_death_rate;
}
*/
//test network
void test() {
	
	//num = num_array[2];
	//cout << num << endl;
	
	srand((unsigned)time(NULL));
	vector<Node> NODES;
	vector<Gauss_Markov> GM;
	
	for (int i = 0; i < num; i++) NODES.push_back(Node(i)); //ID: 0 to n-1
	for (int i = 0; i < num; i++) GM.push_back(Gauss_Markov(i)); // assign initial speed
	Social_network S(num, NODES);
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
		//cout << i << ":" << NODES[i].degree_S << " " << S.get_degree(i) << endl;
	//}
	//BA model
	/*
	char buffer[80];
	fstream file;
	file.open("C:\\Users\\User\\Documents\\Github  liberary\\Matlab\\Result data\\degree.txt", ios::in);
	do {
		file.getline(buffer,sizeof(buffer));
		S.degree_individual.push_back(atoi(buffer));
		S.degree_count[atoi(buffer)]++;

	} while (!file.eof());
	map<int, int>::iterator it;
	for (it = S.degree_count.begin(); it != S.degree_count.end(); it++) {
		cout << it->first << " " << it->second << endl;
	}
	*/
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
	// Analytical model test

	//covid19 dataset
	DataStore(P, S);
	

	return;
	//FB dataset
	vector<Node> NODES_FB;
	for (int i = 0; i < 4039; i++) NODES_FB.push_back(Node(i)); //ID: 0 to n-1
	Social_network S_fb(4039, NODES_FB, 1);
	vector<Gauss_Markov> GM_FB;
	for (int i = 0; i < 4039; i++) GM_FB.push_back(Gauss_Markov(i)); // assign initial speed
	/*
	map<int, int>::iterator it;
	for (it = S_fb.degree_count.begin(); it != S_fb.degree_count.end(); it++) {
		cout << it->first << " " << it->second << endl;
	}
	*/

	for (int k = 1; k <= 1; k++) {
		//P.setting_range(k);
		cout << "Third: Analytical model" << endl;
		unordered_map<int, vector<vector<double>>> record;
		vector<vector<double>> res3;
		vector<Node_a> NODES_A;
		res3 = process_a(NODES_A, S_fb, P, NODES_FB, k, record);
		//res3 = process_a(NODES_A, S, P, NODES, k, record);
		Exp_Record(res3, 2, k, 18); //要修改

		//第一種模擬方式
		vector<vector<double>> res2(iter_times, vector<double>(6, 0.0));
		
		/*
		vector<vector<double>> tmp;
		for (int i = 0; i < time_stamp; i++) {
			cout << i << endl;
			//tmp = Montecarlo(NODES, GM, P, S, record);
			tmp = Montecarlo_modified(NODES, GM, P, S, record);
			//store
			for (int j = 0; j < tmp.size(); j++) {
				for (int q = 0; q < tmp[j].size(); q++) {
					res2[j][q] += tmp[j][q] / (double)time_stamp;
				}
			}
		}
		*/
		
		for (int z = 0; z < ScatterTime; z++) {
			//Physical_network P1(num, NODES);
			Physical_network P1(4039, NODES_FB);
			//每次都重新撒點
			cout << z << endl;
			parallel_for(int(0), time_stamp, [&](int i) {
				cout << i << endl;
				vector<vector<double>> tmp;
				tmp = Montecarlo(NODES_FB, GM_FB, P1, S_fb, record);
				//tmp = Montecarlo(NODES, GM, P1, S, record);
				//tmp = Montecarlo_modified(NODES, GM, P1, S, record);
				//tmp = Montecarlo_modified(NODES_FB, GM_FB, P1, S_fb, record);
				//store
				for (int j = 0; j < tmp.size(); j++) {
					for (int q = 0; q < tmp[j].size(); q++) {
						res2[j][q] += tmp[j][q] / (double)(time_stamp*ScatterTime);
					}
				}
			});
		}
		
		Exp_Record(res2, 1, k, 17);
		
		
		//Threshold
		//bool res = Calculating_Threshold(NODES_A[0], num);
		//cout << res << endl;
		
		//test other model
		//cout << "liu2017web: " << endl;
		//liu2017web A(num, 50);
		//vector<vector<double>> res_liu2017web = A.Process(record, NODES_A[0]);
		//Exp_Record(res_liu2017web, 2, 3, 15);
		/*
		zhao2019virus B(num, 50);
		vector<vector<double>> res_zhao2019virus = B.Process(S);
		Exp_Record(res_zhao2019virus, 2, 4, 15);
		shen2019snird C(num, 50);
		vector<vector<double>> res_shen2019snird = C.Process();
		Exp_Record(res_shen2019snird, 2, 5, 15);
		shen2020epidemiology D(num, 50);
		vector<vector<double>> res_shen2020epidemiology = D.Process();//未完成
		Exp_Record(res_shen2020epidemiology, 2, 6, 15);
		*/
		//cout << "acarali2019modelling: " << endl;
		//acarali2019modelling E(num, 50);
		//vector<vector<double>> res_acarali2019modelling = E.Process(record, NODES_A[0]);
		//Exp_Record(res_acarali2019modelling, 2, 4, 15);
	}
	
	
	
	//第二種模擬方式
	//vector<vector<double>> res2;
	//cout << "Second:" << endl;
	//res2 = Montecarlo(NODES, GM, P, S);
	//record_result_and_print(res2, 2);

	return;
}

