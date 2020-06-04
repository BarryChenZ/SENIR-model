#pragma once
//Model function I proposed
//Directly computing
#include<iostream>
#include<vector>
#include<math.h>
#include<vector>
#include "Distribution.h"
#include "Mobility.h"

#define PI acos(-1)
using namespace std;

struct Node_a {//for analytical
	int num, degree;
	double Area_i, Op_k;
	vector<double> state;
	vector<double> v;//velocity for node
};
//parameters
int number = 1000; // 4
int total_time = 50;
double max_x = 1000, max_y = 1000;

double contact_rate = 0.5;//Rate_c
double success_prob = 0.5;//P_sucess
double open_rate = 0.8;//Rate_o
double scan_rate = 0.05;//
double collision_cost = 0.2;
double prob_NI = 0.8;
//wake_up_rate, Loss_imu_rate, Recover_rate
//death_rare, extra_death_rate
double omega = 0.3;
double omega_r = 0.2;
double gamma = 0.1;
double lambda = 0.1;
double delta = 0.05;
double ex_delta = 0.01;

double new_node = 0.01;

void computeArea_OpK(vector<Node_a>& NODES_A, Physical_network P, Node_a& i) {
	double mean_v = 0.0;
	for (int j = 0; j < i.v.size(); j++) {
		mean_v += i.v[j];
	}
	mean_v = mean_v / (double)i.num;
	//cout << mean_v << " " << P.get_range() << endl;
	i.Area_i = PI * pow(P.get_range(), 2) + 2 * P.get_range()*mean_v * 1;//t = 1;
	double total_link = 0.0;
	for (int j = 0; j < NODES_A.size(); j++){
		total_link += NODES_A[j].num * NODES_A[j].degree;
	}
	i.Op_k = (double)(i.degree * i.num) / total_link;
	return;
}
/*
double prob_SE(Node_a i, Social_network& S, vector<Node_a> NODES_A) {//long
	double E_I = 0;
	double degree = 0;
	for (int j = 0; j < number; j++) {
		if (S.adj_matrix[i.ID][j] == 1 && i.ID != j) {
			degree++;
			if (NODES_A[j].state == 1 || NODES_A[j].state == 3) E_I++;
		}
	}
	return contact_rate * success_prob*E_I / degree;
}
double prob_EI(Node_a i, Social_network& S, vector<Node_a> NODES_A) {
	double I = 0;
	double degree = 0;
	for (int j = 0; j < number; j++) {
		if (S.adj_matrix[i.ID][j] == 1 && i.ID != j) {
			degree++;
			if (NODES_A[j].state == 3) I++;
		}
	}
	return open_rate*I / degree;
}
double prob_NI(Node_a i, vector<Node_a> NODES_A, vector<Node> NODES) {//short
	double tmp = 0.0;
	for (int j = 0; j < NODES[i.ID].neighbor_set.size(); j++) {
		if (NODES_A[NODES[i.ID].neighbor_set[j]].state == 3) tmp++;
	}
	return (1 - collision_cost * tmp)*prob_NIC;
}
double prob_SN(Node_a i, Physical_network P, vector<Node_a> NODES_A) {
	double tmp = 0.0;
	for (int j = 0; j < NODES_A.size(); j++) {
		if(NODES_A[j].state == 3) tmp++;
	}
	return i.Area_i *number*scan_rate / P.area();
}
*/
/*
int Compute_S(Node_a i, double guess) {//Desicision a node for one time
	if (guess <= i.prob_E) return 1;
	else if (guess <= i.prob_E + i.prob_N) return 2;
	else if (guess <= i.prob_E + i.prob_N + i.prob_D) return 5;
	else return 0;
}

int Compute_E(Node_a i, double guess) {
	if (guess <= i.prob_I) return 3;
	else if (guess <= i.prob_I + i.prob_D) return 5;
	else return 1;
}

int Compute_N(Node_a i, double guess) {
	if (guess <= i.prob_I) return 3;
	else if (guess <= i.prob_I + i.prob_D) return 5;
	else return 2;
}

int Compute_I(Node_a i, double guess) {
	if (guess <= i.prob_R) return 4;
	else if (guess <= i.prob_R + i.prob_D) return 5;
	else return 3;
}

int Compute_R(Node_a i, double guess) {
	if (guess <= i.prob_D) return 5;
	else if (guess <= i.prob_D + i.prob_S) return 0;
	else return 4;
}

int Compute_D(Node_a i, double guess) {
	if (i.state == 5) {
		if (guess <= i.prob_S) return 0;
		return 5;
	}
	else {//6
		if (guess <= i.prob_R) return 4;
		return 6;
	}
}
*/
void setting(Node i) {
	omega = i.wake_up_rate;
	//received_rate = i.received_rate;
	lambda = i.lose_immunity_rate;
	delta = i.death_rate;
	ex_delta = i.ex_death_rate;
}
/*
void Compute_prob_up(Node_a& i, Social_network& S, Physical_network P, vector<Node_a>& NODES_A, vector<Node> NODES, vector<int>& temp) {//total
	for (int j = 0; j < number; j++) {
		Node_a tmp;
		tmp.prob_S = i.prob_D * wake_up_rate + i.prob_R * recovered_rate - prob_SE(i,S,NODES_A) * i.prob_S - prob_SN(i, P, NODES_A) * i.prob_S - i.prob_S * death_rate;
		tmp.prob_E = prob_SE(i, S, NODES_A) *i.prob_S - prob_EI(i, S, NODES_A) * i.prob_E - death_rate * i.prob_E;
		tmp.prob_N = prob_SN(i, P, NODES_A) *i.prob_S - prob_NI(i,  NODES_A, NODES)*i.prob_N - (death_rate + ex_death_rate) * i.prob_N;
		tmp.prob_I = prob_EI(i, S, NODES_A)*i.prob_E + prob_NI(i, NODES_A, NODES)*i.prob_N - (death_rate + ex_death_rate) * i.prob_I - recovered_rate * i.prob_I;
		tmp.prob_R = recovered_rate * i.prob_I - death_rate * i.prob_R - lose_immunity_rate * i.prob_R + wake_up_rate*i.prob_D;
		tmp.prob_D = (i.prob_S+ i.prob_E+ i.prob_R)*death_rate + (i.prob_N+i.prob_I)*(death_rate + ex_death_rate) - wake_up_rate * i.prob_D;
		//store
		i.prob_S = tmp.prob_S;
		i.prob_E = tmp.prob_E;
		i.prob_N = tmp.prob_N;
		i.prob_I = tmp.prob_I;
		i.prob_R = tmp.prob_R;
		i.prob_D = tmp.prob_D;
		//updating state
		i.state = temp[i.ID];
	}
}
*/
/*
void Judgement(vector<int>& temp, vector<Node_a>& NODES_A) {
	static default_random_engine e;
	static uniform_real_distribution<double> u(0, 1);
	for (int i = 0; i < number; i++) {
		//cout << i << endl;
		double guess = u(e);
		int state = NODES_A[i].state;
		switch (state) {
		case 0:
			temp[i] = Compute_S(NODES_A[i], guess);
			break;
		case 1:
			temp[i] = Compute_E(NODES_A[i], guess);
			break;
		case 2:
			temp[i] = Compute_N(NODES_A[i], guess);
			break;
		case 3:
			temp[i] = Compute_I(NODES_A[i], guess);
			break;
		case 4:
			temp[i] = Compute_R(NODES_A[i], guess);
			break;
		case 5:
			temp[i] = Compute_D(NODES_A[i], guess);
			break;
		case 6:
			temp[i] = Compute_D(NODES_A[i], guess);
			break;
		}
	}
	return;
}
*/
//SENIRD or SEIRD model
void Compute_prob_frac(double E, double N, double I, double Op_k, double area, double area_i, vector<double>& tmp, vector<double> state) {
	//case 1 delta_ex != 0 without omega_r
	//case 2 with omega_r
	//case 3 delta_ex = 0
	//case 4 SENIRD -> SEIRD model and D不會回溯 + alpha(有新點)
	//cout << Op_k << " " << area_i << endl;
	//S
	tmp[0] = -1 * contact_rate*success_prob*Op_k*(E + I)*state[0] - (scan_rate*(area_i / area * ((E + I)*number)))*state[0] - delta * state[0] + lambda * state[4] - gamma * state[0] + new_node;//case 4
	//tmp[0] = -1 * contact_rate*success_prob*Op_k*(E + I)*state[0] - (scan_rate*(area_i/area *((N+I)*number)))*state[0] - delta * state[0] + omega * state[5] + lambda * state[4];//area_i * density *scan rate
	//E
	tmp[1] = contact_rate * success_prob*Op_k*(E + I)*state[0]+ (scan_rate*(area_i / area * ((E + I)*number)))*state[0] - Op_k * I * open_rate*state[1] - (1 - collision_cost * (I))*prob_NI*state[1] - delta * state[1] - gamma * state[1];//case 4
	//tmp[1] = contact_rate*success_prob*Op_k*(E + I)*state[0] - Op_k*I*open_rate*state[1] - delta * state[1];
	//N
	//tmp[2] = (scan_rate*(area_i/area *((N+I)*number)))*state[0] - (1-collision_cost*(I))*prob_NI*state[2] - (delta+ex_delta)*state[2];
	// tmp[2] = (scan_rate*(area_i*I) / area)*state[0] - (1 - collision_cost * (I / area))*prob_NI*state[2] - (delta)*state[2];//case 3
	//I
	tmp[3] = Op_k * I*open_rate*state[1] + (1 - collision_cost * (I))*prob_NI*state[1] - (delta + ex_delta + gamma)*state[3];//case 4
	//tmp[3] = Op_k*I*open_rate*state[1] + (1-collision_cost*(I))*prob_NI*state[2] - (delta + ex_delta + gamma)*state[3];
	// tmp[3] = Op_k * I*open_rate*state[1] + (1 - collision_cost * (I / area))*prob_NI*state[2] - (delta + gamma)*state[3];//case 3
	//R
	tmp[4] = gamma * state[3] +gamma * state[0] + gamma * state[1] - (delta + lambda)*state[4]; // case 4
	//tmp[4] = gamma*state[3]-(delta+lambda)*state[4];
	//tmp[4] = gamma * state[3] - (delta + lambda)*state[4] + omega_r * state[5]; //case 2
	//D
	tmp[5] = delta * (state[0] + state[1] + state[4]) + (delta + ex_delta)*state[3] - new_node;//case 4
	//tmp[5] = delta*(state[0]+state[1]+state[4])+(delta+ex_delta)*(state[2]+state[3])-omega*state[5];
	//tmp[5] = delta * (state[0] + state[1] + state[4]) + (delta + ex_delta)*(state[2] + state[3])-omega*state[5] -omega_r * state[5];//case 2
	//tmp[5] = delta * (state[0] + state[1] + state[4]) + (delta)*(state[2] + state[3])-omega*state[5] * state[5];//case 3

}	
void update(vector<double>& tmp, vector<double>& state) {
	for (int i = 0; i < state.size(); i++) {
		state[i] = (state[i] + tmp[i] <= 0.0) ? 0.0 : state[i] + tmp[i];
		//cout << tmp[i] << " ";
	}
	//cout << endl;
	return;
}
void Printing(vector<Node_a> NODES_A, vector<double>& res1) {
	
	vector<double> res(6, 0.0);//SENIRD
	for (int i = 0; i < NODES_A.size(); i++) {
		for (int j = 0; j < NODES_A[i].state.size(); j++) {//fraction
			res[j] += ((NODES_A[i].state[j])* (double)NODES_A[i].num)/(double)number;
		}
	}
	for (int i = 0; i < res.size(); i++) {
		res1[i] = res[i];
		cout << res[i] << " ";
	}
	cout << endl;
}

vector<vector<double>> process_a(vector<Node_a> NODES_A, Physical_network P) {
	vector<vector<double>> res;
	res.resize(total_time, vector<double>(6, 0.0));
	double E = 0.0, N = 0.0, I = 0.0, TE = 0.0, TN = 0.0, TI = 0.0;//frac and total
	Gauss_Markov GM = Gauss_Markov();
	//initial test only 3種degree 50 20 10/10% 30% 60%
	NODES_A.resize(3);
	vector<vector<double>> tmp(NODES_A.size(), vector<double>(6, 0.0));//SENIRD:012345
	NODES_A[0].v.resize(number * 0.1), NODES_A[0].num = number * 0.1, NODES_A[0].degree = 50;
	NODES_A[1].v.resize(number * 0.3), NODES_A[1].num = number * 0.3, NODES_A[1].degree = 20;
	NODES_A[2].v.resize(number * 0.6), NODES_A[2].num = number * 0.6, NODES_A[2].degree = 10;
	for(int i = 0; i < NODES_A.size(); i++) GM.initial_v(NODES_A[i].v);
	//state initialize 10% infection at t = 0
	for (int i = 0; i < NODES_A.size(); i++) {
		NODES_A[i].state.resize(6, 0.0);
		NODES_A[i].state[0] = 1.0 - 0.1, NODES_A[i].state[3] = 0.1;
	}
	cout << 0 << endl;
	Printing(NODES_A, res[0]);
	int t = 1;
	double area = max_x * max_y;
	
	// start
	while (t < total_time) {
		E = 0.0, N = 0.0, I = 0.0;
		TE = 0.0, TN = 0.0, TI = 0.0;
		//Mobility model
		for (int i = 0; i < NODES_A.size(); i++) {
			for (int j = 0; j < NODES_A[i].v.size(); j++) {
				NODES_A[i].v[j] = GM.change_velocity(NODES_A[i].v[j]);
			}
		}

		//first calculating total fraction of ENI
		for (int i = 0; i < NODES_A.size(); i++) {
			//cout << NODES_A[i].state[1] * (double)NODES_A[i].num /number << endl;
			E += (double)(NODES_A[i].state[1] * (double)NODES_A[i].num) / (double)number;
			N += (double)(NODES_A[i].state[2] * (double)NODES_A[i].num) / (double)number;
			I += (double)(NODES_A[i].state[3] * (double)NODES_A[i].num) / (double)number;
			TE += NODES_A[i].state[1] * (double)NODES_A[i].num;
			TN += NODES_A[i].state[2] * (double)NODES_A[i].num;
			TI += NODES_A[i].state[3] * (double)NODES_A[i].num;
		}
		//cout << "TE: " << TE << " TN:" << TN << " TI:" << TI << endl;
		//cout << "E:" << E << " N:" << N << " I:" << I << endl;
		//Calculating the fraction change of degree group
		for (int i = 0; i < NODES_A.size(); i++) {
			computeArea_OpK(NODES_A, P, NODES_A[i]);//Area會變 但是Op_k不變
			Compute_prob_frac(E, N, I, NODES_A[i].Op_k, area, NODES_A[i].Area_i, tmp[i], NODES_A[i].state);
		}
		//updating
		for (int i = 0; i < NODES_A.size(); i++) {
			update(tmp[i], NODES_A[i].state);
		}
		//print fraction of state
		cout << t << endl;
		Printing(NODES_A, res[t]);
		t++;
		cout << contact_rate * success_prob*NODES_A[0].Op_k*(E + I) << " ";
		cout << (scan_rate*(NODES_A[0].Area_i / area * ((E + I)*number))) << " ";
		cout << NODES_A[0].Op_k*I*open_rate << " ";
		cout << (1 - collision_cost * (I))*prob_NI << endl;
	}
	cout << NODES_A[0].Area_i << " " << NODES_A[0].Op_k << endl;
	return res;
}