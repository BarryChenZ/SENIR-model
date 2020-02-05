// Final_paper_simulation.cpp : 此檔案包含 'main' 函式。程式會於該處開始執行及結束執行。
//
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

using namespace std;

int n = 4039;
double infected_rate = 0.32;
double received_rate = 0.228;
double iNsidious_rate = 0.8;
double recovered_rate = 0.083; 
double infected_sim_array[4] = { 0.2, 0.3, 0.4, 0.5 };

int iter_times = 100;
int time_stamp = 12;

vector<vector<double>> Montecarlo_2(int** adj_matrix);
void record_result_to_file(vector<vector<double>>& result, int type);
int main()
{
	srand((unsigned)time(NULL));
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
	vector<vector<double>> rescor2 = Montecarlo_2(adj_matrix);
	vector<vector<double>> result2;
	for (int i = 0; i < 4; i++) {
		vector<double> store2;
		cout << infected_sim_array[i] << " : " << endl;
		for (int k = 0; k < time_stamp + 1; k++) {
			double x = 0.0;
			for (int j = i * iter_times; j < iter_times * (i + 1); j++) {
				x += rescor2[j][k];
			}
			store2.push_back(x / iter_times);
			cout << x / iter_times << endl;
		}
		result2.push_back(store2);
		cout << endl;
	}
	record_result_to_file(result2, 2);

	return 0;
}

vector<vector<double>> Montecarlo_2(int** adj_matrix) {
	vector<vector<double>> ans;
	for (int x = 0; x < 4; x++) {
		int t = 0;
		while (t < iter_times) {
			cout << infected_sim_array[x] << " " << t << " : ";
			vector<int> state(n, 0);
			vector<double> tmp;
			//起始感染點
			int a = (rand() % n);
			state[a] = 2;

			//SEIS model
			tmp.push_back(1 / n); // 初始化只有一個點 + 跑12次
			for (int i = 0; i < time_stamp; i++) {
				vector<int> temp(n, 0);

				//依據上一次的step做
				exposedProcess(state, temp, adj_matrix);
				iNsidiousProcess(state, temp);
				infectedProcess_2(state, temp, x);
				recoveredProcess(state, temp);

				update_state(state, temp);
				record(state, tmp);
			}
			cout << endl;
			t++;
			ans.push_back(tmp);
		}
	}
	return ans;
};

void record_result_to_file(vector<vector<double>>& result, int type) {
	fstream file;
	if (type == 1)file.open("C:\\Users\\User\\Desktop\\Matlab Code\\MonteCarlo_Simulation\\sim_result.txt", ios::out);
	else         file.open("C:\\Users\\User\\Desktop\\Matlab Code\\MonteCarlo_Simulation\\sim_result_2.txt", ios::out);
	for (int i = 0; i < result.size(); i++) {
		for (int j = 0; j < result[0].size(); j++) {
			file << result[i][j] << " ";
		}
		file << endl;
	}
};
