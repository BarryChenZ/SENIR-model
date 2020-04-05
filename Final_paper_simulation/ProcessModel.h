#pragma once
//Simulation and analytical model part
#include <iostream>
#include <random>
#include <cmath>
#include <time.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <vector>
#include "Mobility.h"
#include "Distribution.h"
using namespace std;

int n = 4039;
double infected_rate = 0.32;
double received_rate = 0.228;
double iNsidious_rate = 0.8;
double recovered_rate = 0.083; 
double infected_sim_array[4] = { 0.2, 0.3, 0.4, 0.5 };
double death_rate = 0.2;

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


//0 : sus 
//1 : expose 
//2 : iNsidious 
//3 : infected 
//4 : recovered
//5:  death

//S -> E for long range
void exposedProcess(vector<Node>& NODES, vector<int>& temp, vector<vector<int>> adj_matrix) { // S->E
#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		int node_state = NODES[i].get_state();
		if (node_state == 0) {
			int infected_neighbor_num = 0;
#pragma omp parallel for
			for (int j = 0; j < n; j++) {
				if (adj_matrix[j][i] == 1 && NODES[j].get_state() == 3)
					infected_neighbor_num++;
			}
			double infected_probaility = 1 - pow(1 - received_rate, infected_neighbor_num);
			if (guessTrue(infected_probaility)) node_state = 1;
			temp[i] = node_state;
		}
	}
};

//S->N for short range
void iNsidiousProcess(vector<Node>& NODES, vector<int>& temp) { // E->N
	for (int i = 0; i < n; i++) {
		int node_state = NODES[i].get_state();
		if (node_state == 0) {
			if (guessTrue(iNsidious_rate)) node_state = 2;
			temp[i] = node_state;
		}
	}
}

//E->I and N->I
void infectedProcess_1(vector<Node>& NODES, vector<int>& temp) { // N->I
	for (int i = 0; i < n; i++) {
		int node_state = NODES[i].get_state();
		if (node_state == 1) {
			if (guessTrue(infected_rate)) node_state = 3;
			temp[i] = node_state;
		}
	}
};
void infectedProcess_2(vector<Node>& NODES, vector<int>& temp) {
	for (int i = 0; i < n; i++) {
		int node_state = NODES[i].get_state();
		if (node_state == 2) {
			if (guessTrue(infected_rate)) node_state = 3;
			temp[i] = node_state;
		}
	}
}
//I->R
void recoveredProcess(vector<Node>& NODES, vector<int>& temp) { // I->S
	for (int i = 0; i < n; i++) {
		int node_state = NODES[i].get_state();
		if (node_state == 3) {
			if (guessTrue(recovered_rate)) node_state = 4;
			temp[i] = node_state;
		}
	}
};
//any state -> D
void deathProcess(vector<Node>& NODES, vector<int>& temp) {
	for (int i = 0; i < n; i++) {
		if(guessTrue(death_rate)) temp[i] = 5;
	}
}
//D -> S wake up and lose malware
void wakeupProcess(vector<Node>& NODES, vector<int>& temp) {
	for (int i = 0; i < n; i++) {
		int node_state = NODES[i].get_state();
		if (node_state == 5) {
			if (guessTrue(NODES[i].wake_up_rate)) node_state = 0;
			temp[i] = node_state;
		}
	}
}
void loseImProcess(vector<Node>& NODES, vector<int>& temp) {
	for (int i = 0; i < n; i++) {
		int node_state = NODES[i].get_state();
		if (node_state == 4) {
			if (guessTrue(NODES[i].lose_immunity_rate)) node_state = 0;
			temp[i] = node_state;
		}
	}
}

void update_state(vector<Node>& NODES, vector<int>& temp) {
#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		NODES[i].change_state(temp[i]);
	}
	return;
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