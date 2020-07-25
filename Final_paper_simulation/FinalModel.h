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
	//experiment 1 
	double contact_rate_array[3] = {4, 2, 6};
	double contact_rate = 20;//Rate_c 
	double success_prob = 0.8;//P_sucess
	double open_rate = 90;//Rate_o prob
	double open_rate_array[3] = {0.1, 0.3, 0.8};
	double scan_rate = 0.0;// 0.005;//1-100不會爆掉
	double collision_cost = 0.2;
	double collision_cost_array[3] = {2000, 5000, 10000};//Decay
	double prob_NI = 0.0;// 0.6;
	double prob_NI_array[3] = {0.2, 0.5, 0.9};
	//wake_up_rate, Loss_imu_rate, Recover_rate
	//death_rare, extra_death_rate
	double omega = 0.3;
	double omega_r = 0.2;
	double gamma = 0.02;//0.01;
	double lambda = 0.1;//0.1;
	double delta = 0.05;
	double delta_array[3] = { 0.05, 0.04, 0.4 };
	double ex_delta = 0.01;
	double ex_delta_array[3] = { 0.01, 0.08, 0.1 };
	double effect_delta;
	double new_node = 0.04;
	double new_node_array[3] = { 0.01, 0.03, 0.06 };
};
//cout << NODES_A[0].contact_rate * NODES_A[0].success_prob*NODES_A[0].Op_k*(E + I) << " ";
//cout << (NODES_A[0].scan_rate*(NODES_A[0].Area_i * ((E + I) / area))) << " ";
//cout << NODES_A[0].Op_k*I*NODES_A[0].open_rate << " ";
//cout << (1 - NODES_A[0].collision_cost * (I))*NODES_A[0].prob_NI << endl;
//parameters
int number = 1000; // 4
int FB_number = 4039;
int number_array[3] = { 3000,5000,10000 };
int total_time = 50;
double max_x = 1000, max_y = 1000;

double contact_rate = 5;//Rate_c
double success_prob = 0.0005;//P_sucess
double open_rate = 0.08;//Rate_o
double scan_rate = 0.005;//
double collision_cost = 0.02;
double prob_NI = 0.2;
//wake_up_rate, Loss_imu_rate, Recover_rate
//death_rare, extra_death_rate
double omega = 0.3;
double omega_r = 0.2;
double gamma = 0.1;
double lambda = 0.1;
double delta = 0.05;//0.05
double ex_delta = 0.01;
double ex_delta_array[3] = { 0.01, 0.05, 0.09};
double new_node = 0.04;

//contact_rate* success_prob* NODES_A[0].Op_k* (E + I) = exposed rate
//(scan_rate * (NODES_A[0].Area_i / area * ((E + I) * number))) = iNsidious rate
//NODES_A[0].Op_k * I * open_rate = infected_rate
//(1 - collision_cost * (I)) * prob_NI = infected_rate_2
// gamma = recovered_rate
// lambda = lose_immunity_rate
// delta = death_rate
// ex_delta = ex_death_rate
// new_node = new rate

void setting(vector<Node>& NODES, Node_a j) {
	for (int i = 0; i < NODES.size(); i++) {
		//4 value set by average
		NODES[i].recovered_rate = j.gamma;
		NODES[i].death_rate = j.delta;
		NODES[i].ex_death_rate = j.ex_delta;
		NODES[i].lose_immunity_rate = j.lambda;
		NODES[i].new_rate = j.new_node;
	}
	return;
}
void computeArea_OpK(vector<Node_a>& NODES_A, Physical_network P, Node_a& i) {
	double mean_v = 0.0;
	for (int j = 0; j < i.v.size(); j++) {
		mean_v += i.v[j];
	}
	mean_v = mean_v / (double)i.num;
	i.Area_i = PI * pow(P.get_range(), 2) + 2 * P.get_range()*mean_v * 1;//t = 1;
	//cout << i.Area_i << endl;
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
void Compute_prob_frac(double S, double E, double N, double I, double TE, double TN, double TI, double area, vector<double>& tmp, Node_a i) {
	//case 1 delta_ex != 0 without omega_r
	//case 2 with omega_r
	//case 3 delta_ex = 0
	//case 4 SENIRD -> SEIRD model and D不會回溯 + alpha(有新點)
	//cout << Op_k << " " << area_i << endl;
	//S
	tmp[0] = -1 * i.contact_rate*i.success_prob*i.Op_k*(E + I)*i.state[0] - (i.scan_rate*(i.Area_i*((TE + TI))/area))*i.state[0] - i.delta * i.state[0] + i.lambda * i.state[4] - i.gamma * i.state[0] + i.new_node;//case 4
	//tmp[0] = -1 * contact_rate*success_prob*Op_k*(E + I)*state[0] - (scan_rate*(area_i/area *((N+I)*number)))*state[0] - delta * state[0] + omega * state[5] + lambda * state[4];//area_i * density *scan rate
	//E
	tmp[1] = i.contact_rate * i.success_prob*i.Op_k*(E + I)*i.state[0]+ (i.scan_rate*(i.Area_i*((TE + TI)) / area))*i.state[0] - i.Op_k * I*i.open_rate*i.state[1] - (1 - i.collision_cost * (I))*i.prob_NI*i.state[1] - i.delta * i.state[1] - i.gamma * i.state[1];//case 4
	//tmp[1] = contact_rate*success_prob*Op_k*(E + I)*state[0] - Op_k*I*open_rate*state[1] - delta * state[1];
	//N
	//tmp[2] = (scan_rate*(area_i/area *((N+I)*number)))*state[0] - (1-collision_cost*(I))*prob_NI*state[2] - (delta+ex_delta)*state[2];
	//tmp[2] = (scan_rate*(area_i*I) / area)*state[0] - (1 - collision_cost * (I / area))*prob_NI*state[2] - (delta)*state[2];//case 3
	//I
	tmp[3] = i.Op_k*I*i.open_rate*i.state[1] + (1 - i.collision_cost * (I))*i.prob_NI*i.state[1] - (i.delta + i.ex_delta + i.gamma)*i.state[3];//case 4
	//tmp[3] = Op_k*I*open_rate*state[1] + (1-collision_cost*(I))*prob_NI*state[2] - (delta + ex_delta + gamma)*state[3];
	// tmp[3] = Op_k * I*open_rate*state[1] + (1 - collision_cost * (I / area))*prob_NI*state[2] - (delta + gamma)*state[3];//case 3
	//R
	//tmp[4] = i.gamma * i.state[3] + i.gamma * i.state[0] + i.gamma * i.state[1] - (i.delta + i.lambda)*i.state[4]; // case 4
	tmp[4] = i.gamma * i.state[3] +i.gamma * i.state[0] + i.gamma * i.state[1] - (i.effect_delta + i.lambda)*i.state[4]; // case 4
	//tmp[4] = gamma*state[3]-(delta+lambda)*state[4];
	//tmp[4] = gamma * state[3] - (delta + lambda)*state[4] + omega_r * state[5]; //case 2
	//D
	tmp[5] = i.delta * (i.state[0] + i.state[1]) + (i.delta + i.ex_delta)*i.state[3] + i.effect_delta*i.state[4] - i.new_node;//case 4
	//tmp[5] = i.delta * (i.state[0] + i.state[1]) + (i.delta + i.ex_delta)*i.state[3] + i.delta*i.state[4] - i.new_node;//case 4
	//tmp[5] = delta*(state[0]+state[1]+state[4])+(delta+ex_delta)*(state[2]+state[3])-omega*state[5];
	//tmp[5] = delta * (state[0] + state[1] + state[4]) + (delta + ex_delta)*(state[2] + state[3])-omega*state[5] -omega_r * state[5];//case 2
	//tmp[5] = delta * (state[0] + state[1] + state[4]) + (delta)*(state[2] + state[3])-omega*state[5] * state[5];//case 3

}	
void update(vector<double>& tmp, vector<double>& state) {
	for (int i = 0; i < state.size(); i++) {
		state[i] = (state[i] + tmp[i] < 0) ? (double)0.0 : state[i] + tmp[i];
		//state[i] = state[i] + tmp[i];
		//cout << tmp[i] << " ";
	}
	//cout << endl;
	return;
}
void Printing(vector<Node_a> NODES_A, vector<double>& res1) {
	
	vector<double> res(6, 0.0);//SENIRD
	for (int i = 0; i < NODES_A.size(); i++) {
		for (int j = 0; j < NODES_A[i].state.size(); j++) {//fraction
			//res[j] += ((double)(NODES_A[i].state[j])* (double)NODES_A[i].num)/(double)number;
			res[j] += ((double)(NODES_A[i].state[j])* (double)NODES_A[i].num) / (double)FB_number;
		}
	}
	for (int i = 0; i < res.size(); i++) {
		res1[i] = res[i];
		cout << res[i] << " ";
	}
	cout << endl;
}

void storeParameters(vector<Node_a> NODES_A, unordered_map<int, vector<vector<double>>>& record, double E, double N, double I,double TE, double TN, double TI, int t) {
	for (int j = 0; j < NODES_A.size(); j++) {
		//record[NODES_A[j].degree][t][0] = (NODES_A[j].contact_rate * NODES_A[j].success_prob*NODES_A[j].Op_k*(TE + TI));
		record[NODES_A[j].degree][t][0] = (NODES_A[j].contact_rate * NODES_A[j].success_prob*NODES_A[j].Op_k*(E + I));
		//record[NODES_A[j].degree][t][1] = NODES_A[j].scan_rate*(NODES_A[j].Area_i * ((E + I) / (max_x * max_y)));
		record[NODES_A[j].degree][t][1] = NODES_A[j].scan_rate*(NODES_A[j].Area_i * ((TE + TI) / (max_x * max_y)));
		record[NODES_A[j].degree][t][2] = NODES_A[j].Op_k * I * NODES_A[j].open_rate;
		//record[NODES_A[j].degree][t][2] = NODES_A[j].Op_k * TI * NODES_A[j].open_rate;
		record[NODES_A[j].degree][t][3] = (1 - NODES_A[j].collision_cost * (I))*NODES_A[j].prob_NI;
		//record[NODES_A[j].degree][t][3] = (1 - NODES_A[j].collision_cost * (TI) / (max_x * max_y))*NODES_A[j].prob_NI;

		//如果公式有改 這邊也要改
		cout << NODES_A[j].contact_rate * NODES_A[j].success_prob*NODES_A[j].Op_k*(E + I) << " ";//This
		//cout << NODES_A[j].contact_rate * NODES_A[j].success_prob*NODES_A[j].Op_k*(TE + TI) << " ";
		//cout << (NODES_A[j].scan_rate*(NODES_A[j].Area_i * ((E + I) / (max_x * max_y)))) << " ";
		//cout << (NODES_A[j].scan_rate*(NODES_A[j].Area_i * ((TE + TI) / (max_x * max_y)))) << " ";//This
		//cout << NODES_A[j].scan_rate << " " << (NODES_A[j].Area_i/ (max_x * max_y) )* (TE + TI) << endl;
		cout << NODES_A[j].Op_k*I*NODES_A[j].open_rate << " " << endl;//This
		//cout << NODES_A[j].Op_k*TI*NODES_A[j].open_rate << " ";
		//cout << (1 - NODES_A[j].collision_cost * (I))*NODES_A[j].prob_NI << endl;//This
		//cout << (1 - NODES_A[j].collision_cost * (TI) / (max_x * max_y))*NODES_A[j].prob_NI << endl;
		//cout << NODES_A[j].effect_delta << endl;//This
	}
}

//Model threshold and stationary point
bool Calculating_Threshold(Node_a i, int n) {
	double stationary_S, stationary_E, stationary_I, stationary_R, stationary_D, threshold;
	//case 1: Non-constant value for E->I // test ok
	//stationary_S = (i.new_node * (i.delta + i.lambda)) / (i.delta*(i.lambda + i.gamma + i.delta));
	//double C1 = i.contact_rate*i.success_prob*i.Op_k + (i.scan_rate*(i.Area_i / (max_x*max_y) * number / (max_x*max_y)));
	//double C3 = i.Op_k * i.open_rate + i.prob_NI;
	//threshold = (C1 * stationary_S * (i.delta + i.ex_delta + i.gamma)) / (i.delta*i.ex_delta+2*i.delta*i.gamma+pow(i.delta,2)+pow(i.gamma,2));
	//case 2: Constant value for E->I //test ok
	stationary_S = (i.new_node * (i.delta + i.lambda))/(i.delta*(i.lambda+i.gamma+i.delta));//算malware free的
	stationary_R = (i.new_node * i.gamma)/(i.delta*(i.lambda + i.gamma + i.delta));
	stationary_D = 1 - stationary_S - stationary_R;
	stationary_E = 0, stationary_I = 0;
	double C1 = i.contact_rate*i.success_prob*i.Op_k + (i.scan_rate*n*(i.Area_i / (max_x*max_y)));
	double C3 = i.Op_k * i.open_rate + i.prob_NI;
	//threshold = (C1 * stationary_S * (C3 + i.delta + i.ex_delta + i.gamma))/((C3 +i.delta + i.gamma)*(i.delta + i.ex_delta + i.gamma));
	threshold = (C1 * stationary_S * (i.prob_NI + i.delta + i.ex_delta + i.gamma)) / ((i.prob_NI + i.delta + i.gamma)*(i.delta + i.ex_delta + i.gamma));
	cout << "Threshold: " << threshold << endl; // if > 1 I E large , else if i < 1 I E small.
	if (threshold > 1) return true;
	else {
		cout << "stationary_S: " << stationary_S << endl;
		cout << "stationary_R: " << stationary_R << endl;
		cout << "stationary_D: " << stationary_D << endl;
		return false;
	}
}

vector<vector<double>> process_a(vector<Node_a>& NODES_A, Social_network S, Physical_network P, vector<Node>& NODES, int k, unordered_map<int, vector<vector<double>>>& record) {
	//number = number_array[k - 1];
	vector<vector<double>> res;
	res.resize(total_time, vector<double>(6, 0.0));
	double E = 0.0, N = 0.0, I = 0.0, TE = 0.0, TN = 0.0, TI = 0.0;//frac and total
	double S_frac = 0.0;
	//record.resize(total_time, vector<double>(4,0.0)); // record 4 value, prapare to transmit to the simulation parameters.
	Gauss_Markov GM = Gauss_Markov();
	//initial test only 3種degree 50 20 10/10% 30% 60% 舊的
	/*
	NODES_A.resize(3);
	//cout << NODES_A[0].contact_rate_array[k - 1] << endl;//test 要修改1
	NODES_A[0].v.resize(number * 0.1), NODES_A[0].num = number * 0.1, NODES_A[0].degree = 50;
	NODES_A[1].v.resize(number * 0.3), NODES_A[1].num = number * 0.3, NODES_A[1].degree = 20;
	NODES_A[2].v.resize(number * 0.6), NODES_A[2].num = number * 0.6, NODES_A[2].degree = 10;
	*/
	//FB dataset
	NODES_A.resize(S.degree_count.size());
	map<int, int>::iterator it = S.degree_count.begin();
	for (int i = 0; i < NODES_A.size(); i++, it++) {
		NODES_A[i].v.resize(it->second,0);
		NODES_A[i].num = it->second;
		NODES_A[i].degree = it->first;
	}
	//新的BA model模式 會有部分Degree只有一個點
	/*
	NODES_A.resize(S.degree_count.size());
	cout << NODES_A[0].contact_rate_array[k - 1] << endl;//test 要修改1
	vector<vector<double>> tmp(NODES_A.size(), vector<double>(6, 0.0));//SENIRD:012345
	map<int, int>::iterator it = S.degree_count.begin();
	for (int i = 0; i < S.degree_count.size(); i++, it++) {
		NODES_A[i].v.resize(it->second), NODES_A[i].num = it->second, NODES_A[i].degree = it->first;
		cout << NODES_A[i].degree << " " << NODES_A[i].num << endl;
	}
	*/

	//GM.initial_v(NODES_A[i].v);
	//GM.initial_v(NODES_A[i].v, k);
	for (int i = 0; i < NODES_A.size(); i++) GM.initial_v(NODES_A[i].v);//改速度的實驗 Exp7
	//state initialize 10% infection at t = 0
	for (int i = 0; i < NODES_A.size(); i++) {
		NODES_A[i].state.resize(6, 0.0);
		//NODES_A[i].state[0] = 1.0 - 0.1, NODES_A[i].state[3] = 0.1;//fb 不適用
		NODES_A[i].state[0] = 1.0;
		NODES_A[i].effect_delta = NODES_A[i].delta;
		if (k != 0) {//修改參數
			//NODES_A[i].contact_rate = NODES_A[i].contact_rate_array[k - 1]; //Exp1
			//NODES_A[i].open_rate = NODES_A[i].open_rate_array[k - 1]; // Exp2
			//NODES_A[i].prob_NI = NODES_A[i].prob_NI_array[k - 1]; // Exp 3
			//NODES_A[i].collision_cost = NODES_A[i].collision_cost_array[k - 1]; // Exp 4
			//NODES_A[i].new_node = NODES_A[i].new_node_array[k - 1]; // Exp 5
			//NODES_A[i].delta = NODES_A[i].delta_array[k - 1];// Exp 9
			//NODES_A[i].ex_delta = NODES_A[i].ex_delta_array[k - 1];// Exp 9
		}

		if (record.find(NODES_A[i].degree) == record.end()) {
			record[NODES_A[i].degree] = {vector<vector<double>>(total_time, vector<double>(4, 0.0))};
		}
	}
	//分開做 random找
	vector<int> random_infection_initial;
	for (int i = 0; i < 4039; i++) random_infection_initial.push_back(i);
	random_shuffle(random_infection_initial.begin(), random_infection_initial.end());
	for (int i = 0; i < 404; i++) {
		
		for (int j = 0; j < NODES_A.size(); j++) {
			if (NODES_A[j].degree == NODES[random_infection_initial[i]].degree_S) {
				NODES_A[j].state[0] = (NODES_A[j].state[0] * NODES_A[j].num - 1) / (double)NODES_A[j].num;
				NODES_A[j].state[3] = (NODES_A[j].state[3] * NODES_A[j].num + 1) / (double)NODES_A[j].num;
				//cout << NODES_A[j].state[0] << " " << NODES_A[j].state[3] << endl;
				break;
			}
		}
	}
	cout << 0 << endl;
	//for(int i = 0; i < NODES_A.size(); i++) computeArea_OpK(NODES_A, P, NODES_A[i]);
	Printing(NODES_A, res[0]);
	int t = 1;
	double area = max_x * max_y;
	
	// start
	while (t < total_time) {
		E = 0.0, N = 0.0, I = 0.0, S_frac = 0.0;
		TE = 0.0, TN = 0.0, TI = 0.0;
		//Mobility model
		for (int i = 0; i < NODES_A.size(); i++) {
			for (int j = 0; j < NODES_A[i].v.size(); j++) {
				NODES_A[i].v[j] = GM.change_velocity(NODES_A[i].v[j]);
			}
		}
		
		//first calculating total fraction of ENI
		for (int i = 0; i < NODES_A.size(); i++) {
			S_frac += (double)(NODES_A[i].state[0] * (double)NODES_A[i].num) / (double)number;
			E += (double)(NODES_A[i].state[1] * (double)NODES_A[i].num) / (double)number;
			N += (double)(NODES_A[i].state[2] * (double)NODES_A[i].num) / (double)number;
			I += (double)(NODES_A[i].state[3] * (double)NODES_A[i].num) / (double)number;
			TE += (double)(NODES_A[i].state[1] * (double)NODES_A[i].num);
			TN += (double)(NODES_A[i].state[2] * (double)NODES_A[i].num);
			TI += (double)(NODES_A[i].state[3] * (double)NODES_A[i].num); 
		}
		//Calculating the fraction change of degree group
		vector<vector<double>> tmp(NODES_A.size(), vector<double>(6, 0.0));//SENIRD:012345
		for (int i = 0; i < NODES_A.size(); i++) {
			computeArea_OpK(NODES_A, P, NODES_A[i]);//Area會變 但是Op_k不變
			Compute_prob_frac(S_frac, E, N, I, TE, TN, TI, area, tmp[i], NODES_A[i]);
			NODES_A[i].effect_delta = NODES_A[i].delta + (NODES_A[i].ex_delta)*(I)*NODES_A[i].gamma;
			//由前一次的I來決定 所以後算，雖然有算 但是機率變化非常小
		}
		//store rate parameters
		cout << t << endl;
		storeParameters(NODES_A, record, E, N, I, TE, TN, TI, t);
		//cout << NODES_A[0].Area_i << " ";
		//cout << NODES_A[0].contact_rate * NODES_A[0].success_prob*NODES_A[0].Op_k*(E + I) << " ";
		//cout << (NODES_A[0].scan_rate*(NODES_A[0].Area_i * ((E + I) / area))) << " ";
		//cout << NODES_A[0].Op_k*I*NODES_A[0].open_rate << " ";
		//cout << (1 - NODES_A[0].collision_cost * (I))*NODES_A[0].prob_NI << endl;

		//updating
		for (int i = 0; i < NODES_A.size(); i++) {
			update(tmp[i], NODES_A[i].state);
		}
		//print fraction of state and store
		Printing(NODES_A, res[t]);
		t++;
	}
	//cout << NODES_A[0].Area_i << " " << NODES_A[0].Op_k << endl;
	//cout << record[0] << " " << record[1] << " " << record[2] << " " << record[3] << endl;
	setting(NODES, NODES_A[0]);
	return res;
}