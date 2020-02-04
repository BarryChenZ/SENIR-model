#pragma once

#include <iostream>
#include <random>
#include <cmath>
#include <time.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <omp.h>
using namespace std;

int n = 4039;
double infected_rate = 0.32;
double received_rate = 0.228;
double recovered_rate = 0.083; 
double infected_sim_array[4] = { 0.2, 0.3, 0.4, 0.5 };

bool guessTrue(double probability) {
	static default_random_engine e;
	static uniform_real_distribution<double> u(0, 1);
	double guess = u(e);
	if (guess <= probability) {
		return true;
	}
	else {
		return false;
	}
};

void exposedProcess(vector<int>& state, vector<int>& temp, int** adj_matrix) {
#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		int node_state = state[i];
		if (node_state == 0) {
			int infected_neighbor_num = 0;
#pragma omp parallel for
			for (int j = 0; j < n; j++) {
				if (adj_matrix[j][i] == 1 && state[j] == 2)
					infected_neighbor_num++;
			}
			double infected_probaility = 1 - pow(1 - received_rate, infected_neighbor_num);
			if (guessTrue(infected_probaility)) node_state = 1;
			temp[i] = node_state;
		}
	}
};

void infectedProcess_2(vector<int>& state, vector<int>& temp, int x) {
	for (int i = 0; i < n; i++) {
		int node_state = state[i];
		if (node_state == 1) {
			if (guessTrue(infected_sim_array[x])) node_state = 2;
			temp[i] = node_state;
		}
	}
};

void recoveredProcess(vector<int>& state, vector<int>& temp) {
	for (int i = 0; i < n; i++) {
		int node_state = state[i];
		if (node_state == 2) {
			if (guessTrue(recovered_rate)) node_state = 0;
			temp[i] = node_state;
		}
	}
};

void update_state(vector<int>& state, vector<int>& temp) {
#pragma omp parallel for
	for (int i = 0; i < n; i++) state[i] = temp[i];
};

void record(vector<int>& state, vector<double>& tmp) {
	//¥u°O¿ýI state
	double num = 0;
#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		if (state[i] == 2) num++;
	}
	tmp.push_back(double(num / n));
	cout << double(num / n) << " ";
};