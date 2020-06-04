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

//int n = 4039;
double infected_rate = 0.32;
double received_rate = 0.228;
double iNsidious_rate = 0.8;
double recovered_rate = 0.083; 
double infected_sim_array[4] = { 0.2, 0.3, 0.4, 0.5 };
double death_rate = 0.2;
double new_rate = 0.01;// new rate = leave rate
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

double infected_fraction = 0.05; //5%
void initial_infected(vector<Node>& NODES, int n) {
	double tmp = n * infected_fraction;
	int i = 0;
	while (i < tmp) {
		int x = (int)(rand() % n);
		if (NODES[x].get_state() == 0) {
			NODES[x].change_state(3);
			i++;
		}
	}
	return;
}

//0 : sus 
//1 : expose 
//2 : iNsidious 
//3 : infected 
//4 : recovered
//5:  death
//6:  recovered to death

//S -> E for long range SEIRD
void exposedProcess(vector<Node>& NODES, vector<int>& temp, vector<vector<int>> adj_matrix, int n) { // S->E
	for (int i = 0; i < n; i++) {
		int node_state = NODES[i].get_state();
		if (node_state == 0 && temp[i] == -1) {
			int infected_neighbor_num_e = 0;
			for (int j = 0; j < n; j++) {
				//cout << adj_matrix[i][j];
				if (adj_matrix[j][i] == 1 && NODES[j].get_state() == 3 && j != i)
					infected_neighbor_num_e++;
			}
			double infected_probability_e = 1 - pow(1 - NODES[i].exposed_rate, infected_neighbor_num_e);

			int infected_neighbor_num_n = 0;
			for (int j = 0; j < NODES[i].neighbor_set.size(); j++) {
				if (NODES[NODES[i].neighbor_set[j]].get_state() == 3)
					infected_neighbor_num_n++;
			}
			double infected_probability_n = 1 - pow(1 - NODES[i].iNsidious_rate, infected_neighbor_num_n);
			//cout << infected_probability << endl;
			if (guessTrue(infected_probability_e + infected_probability_n)) {
				temp[i] = 1;
			}
		}
	}
};

//S->N for short range
void iNsidiousProcess(vector<Node>& NODES, vector<int>& temp, Physical_network& P, int n) { // E->N
	for (int i = 0; i < n; i++) {
		int node_state = NODES[i].get_state();
		if (node_state == 0 && temp[i] == -1) {
			int infected_neighbor_num = 0;
			for (int j = 0; j < NODES[i].neighbor_set.size(); j++) {
				if (NODES[NODES[i].neighbor_set[j]].get_state() == 3)
					infected_neighbor_num++;
			}
			double infected_probability = 1 - pow(1 - NODES[i].iNsidious_rate, infected_neighbor_num);
			//cout << infected_neighbor_num << endl;
			if (guessTrue(infected_probability)) temp[i] = 2;
		}
	}
}

//E->I and N->I
void infectedProcess_1(vector<Node>& NODES, vector<int>& temp, int n) { // N->I
	for (int i = 0; i < n; i++) {
		int node_state = NODES[i].get_state();
		if (node_state == 1 && temp[i] == -1) {
			if (guessTrue(NODES[i].infected_rate + NODES[i].infected_rate_2)) temp[i] = 3;
		}
	}
};
void infectedProcess_2(vector<Node>& NODES, vector<int>& temp, int n) {
	for (int i = 0; i < n; i++) {
		int node_state = NODES[i].get_state();
		if (node_state == 2 && temp[i] == -1) {
			if (guessTrue(infected_rate)) temp[i] = 3;
		}
	}
}
//I->R || D->R
void recoveredProcess(vector<Node>& NODES, vector<int>& temp, int n) { // I->S
	for (int i = 0; i < n; i++) {
		int node_state = NODES[i].get_state();
		if (node_state == 3 || node_state == 0 || node_state == 1 && temp[i] == -1) {
			if (guessTrue(NODES[i].recovered_rate)) temp[i] = 4;
		}
	}
};
//any state -> D
void deathProcess(vector<Node>& NODES, vector<int>& temp, int n) {
	for (int i = 0; i < n; i++) {
		int node_state = NODES[i].get_state();
		if (node_state == 0 || node_state == 3 || node_state == 4 && temp[i] == -1) {
			if (guessTrue(NODES[i].death_rate)) {
				temp[i] = 5;
			}
		}
		else if(node_state == 1 && node_state == 3 && temp[i] == -1){
			if (guessTrue(NODES[i].death_rate+ NODES[i].ex_death_rate)) temp[i] = 5;
		}
	}
}
//D -> S wake up and lose malware
void wakeupProcess(vector<Node>& NODES, vector<int>& temp, int n) {
	for (int i = 0; i < n; i++) {
		int node_state = NODES[i].get_state();
		if (node_state == 5) {
			if (guessTrue(NODES[i].wake_up_rate)) temp[i] = 0;
		}
	}
}
void loseImProcess(vector<Node>& NODES, vector<int>& temp, int n) {
	for (int i = 0; i < n; i++) {
		int node_state = NODES[i].get_state();
		if (node_state == 4 && temp[i] == -1) {
			if (guessTrue(NODES[i].lose_immunity_rate)) temp[i] = 0;
		}
	}
}

void unchanged(vector<Node>& NODES, vector<int>& temp, int n) {
	for (int i = 0; i < n; i++) {
		int node_state = NODES[i].get_state();
		if (temp[i] == -1) temp[i] = node_state;
	}
}

void update_state(vector<Node>& NODES, vector<int>& temp, int n) {
#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		NODES[i].change_state(temp[i]);
	}
	return;
};

void record(vector<Node>& NODES, vector<double>& tmp, int n) {
	//只記錄I state, tmp 紀錄各時間的I比率
	double num = 0.0;
#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		if (NODES[i].get_state() == 3 || NODES[i].get_state() == 2) num++;
	}
	tmp.push_back(double(num / (double)n));
	cout << double(num / (double)n) << " ";
};

void record_all_state(vector<Node>& NODES, vector<vector<double>>& tmp, int n, int time) {
	vector<double> num_state(6, 0.0);
	for (auto i : NODES) {
		num_state[(i.get_state() == 6)? 5 : i.get_state()]++;
	}
	for (int i = 0; i < num_state.size(); i++) {
		tmp[time][i] = num_state[i] / double(n);
		cout << num_state[i] / double(n) << " ";
	}
	cout << endl;
}

//Modified simulaiton idea.
//focus on node current state. Each node in the network will decided that which state to go. 
//Case 1 : SENIRD, Case 2 : SEIRD
int modified_sus_process(vector<Node>& NODES, vector<int>& temp, vector<vector<int>> adj_matrix, Physical_network& P, int n, int i, double guess) {
	//to E N D S
	double eff_expose;
	double eff_iN;
	int infected_neighbor_num_e = 0;
	int infected_neighbor_num_n = 0;
	for (int j = 0; j < n; j++) {
		//cout << adj_matrix[i][j];
		if (adj_matrix[j][i] == 1 && NODES[j].get_state() == 3 && j != i)
			infected_neighbor_num_e++;
	}
	for (int j = 0; j < NODES[i].neighbor_set.size(); j++) {
		//cout << NODES[i].neighbor_set[j] << endl;
		if (NODES[NODES[i].neighbor_set[j]].get_state() == 3 && j != i) {
			infected_neighbor_num_n++;
			//cout << 1;
		}
	}
	/*
	for (int j = 0; j < infected_neighbor_num_e; j++) {
		if (guessTrue(NODES[i].exposed_rate)) return 1;
	}
	for (int j = 0; j < infected_neighbor_num_n; j++) {
		if (guessTrue(NODES[i].iNsidious_rate)) return 2;
	}
	if (guess < NODES[i].death_rate) return 5;
	return 0;
	*/
	
	eff_expose = 1 - pow(1 - NODES[i].exposed_rate, infected_neighbor_num_e);
	eff_iN = 1 - pow(1 - NODES[i].iNsidious_rate, infected_neighbor_num_n);
	
	//cout << eff_expose << " " << eff_iN << endl;
	/*
	if (guess <= eff_expose) return 1;
	else if (guess <= eff_expose + eff_iN) return 2;
	else if (guess <= eff_expose + eff_iN + NODES[i].death_rate) return 5;
	else return 0;
	*/
	
	if (guess <= eff_expose + eff_iN) return 1;//case 2
	else if (guess <= eff_expose + eff_iN + NODES[i].recovered_rate) return 4;
	else if (guess <= eff_expose + eff_iN + NODES[i].recovered_rate + NODES[i].death_rate) return 5;
	else return 0;
	
}
// Particular node i
int modified_ex_process(Node i, double guess) {

	if (guess <= i.infected_rate + i.infected_rate_2) return 3; // case 2
	else if (guess <= i.infected_rate + i.infected_rate_2 + i.recovered_rate) return 4;
	else if (guess <= i.infected_rate + i.infected_rate_2 + i.recovered_rate + i.death_rate) return 5;
	else return 1;
	/*
	if (guess <= i.infected_rate) return 3;
	else if (guess <= i.infected_rate + i.death_rate) return 5;
	else return 1;
	*/
}

int modified_iN_process(Node i, double guess) {//case 2 不需要
	
	if (guess <= i.infected_rate_2) return 3;
	else if (guess <= i.infected_rate_2 + i.death_rate + i.ex_death_rate) return 5;
	else return 2;
}

int modified_inf_process(Node i, double guess) {
	
	if (guess <= i.recovered_rate) return 4;
	else if (guess <= i.recovered_rate + i.death_rate + i.ex_death_rate) return 5;
	else return 3;
}

int modified_dea_process(Node i, double guess) {//case 2 不需要
	
	if (guess <= i.wake_up_rate) {
		if (i.get_state() == 6) return 4;
		return 0;
	}
	else {
		return 5;
	}
	
}

int modified_rec_process(Node i, double guess) {
	
	if (guess <= i.death_rate) return 5;//5 without omega_r, 6 with omega_r
	else if (guess <= i.death_rate + i.lose_immunity_rate) return 0;
	else return 4;
}
void Leaving_and_Coming(vector<int>& temp, int num) {
	int number = new_rate * num;
	vector<int>  record_D;
	for (int i = 0; i < temp.size(); i++) {
		if (temp[i] == 5) record_D.push_back(i);
	}
	random_shuffle(record_D.begin(), record_D.end());
	if (number > record_D.size()) for (int i = 0; i < record_D.size(); i++) temp[record_D[i]] = 0;
	else for (int i = 0; i < number; i++) temp[record_D[i]] = 0;

	return;
}
//0 : sus 
//1 : expose 
//2 : iNsidious 
//3 : infected 
//4 : recovered
//5:  death
//6:  recover to death

void Judgement(vector<int>& temp, vector<Node>& NODES, vector<vector<int>> adj_matrix, Physical_network& P, int n) {
	static default_random_engine e;
	static uniform_real_distribution<double> u(0, 1);
	for (int i = 0; i < n; i++) {
		//cout << i << endl;
		double guess = u(e);
		int state = NODES[i].get_state();
		switch (state) {
		case 0:
			temp[i] = modified_sus_process(NODES, temp, adj_matrix, P, n, i, guess);
			break;
		case 1:
			temp[i] = modified_ex_process(NODES[i], guess);
			break;
		case 2:
			temp[i] = modified_iN_process(NODES[i], guess);
			break;
		case 3:
			temp[i] = modified_inf_process(NODES[i], guess);
			break;
		case 4:
			temp[i] = modified_dea_process(NODES[i], guess);
			break;
		case 5:
			temp[i] = modified_rec_process(NODES[i], guess);
			break;
		case 6:
			temp[i] = modified_rec_process(NODES[i], guess);
			break;
		}
	}
}
