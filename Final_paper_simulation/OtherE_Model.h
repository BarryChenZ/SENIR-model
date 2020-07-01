#pragma once
//3 types epidemic model from other papers

//OLD
//1. A Two layer Model of Malware Propagation in a Search Engine Context
//2. Modelling the Spread of Botnet Malware in IoT-Based Wireless Sensor NetworksModelling the Spre
//3. Optimal Impulse Control of Bi-Virus SIR Epidemics with Application to Heterogeneous Internet of Things
//4. SNIRD  Disclosing Rules of Malware Spread in Heter wireless Sensor networks
//5. A cloud-assisted malware detection and suppression framework for wireless multimedia system in IoT based on dynamic differential game

//NEW

#include <iostream>
#include <vector>
#include <math.h>
#include <time.h>
#include <fstream>
#include "Distribution.h"

using namespace std;
//1 old
/*
class Community {
private:
	double alpha = 0.1; //search engine prob.
	double beta = 0.2;//0.15 best fit
	double initial_infected = 1;//average
	double infected = 0;
	double suscepted;  
	int number_of_nodes;
public:
	int time = 1;
	Community(int size) {
		number_of_nodes = size;
		suscepted = number_of_nodes - initial_infected;
		infected = initial_infected;
	}
	void infected_process() {
		// Equation (9),(17)
		double a = (2 * initial_infected + alpha * number_of_nodes*(time - 1))*number_of_nodes*exp(beta * time);
		double b = 2 * number_of_nodes + (2 * initial_infected + alpha * number_of_nodes * (time - 1)) * (exp(beta * time) - 1);
		time++;
		if (a / b > number_of_nodes) infected = number_of_nodes;
		else infected = a / b;
	}
	double get_infected() {
		return infected;
	}
};
bool probability(double prob) {
	srand((unsigned)time(NULL));
	double x = (double)rand() / (RAND_MAX + 1.0);
	if (x > prob) return true;
	return false;
}
void infected_community_process(vector<pair<int, Community>>& total_community) {
	// Equation (17)
	int C = 0.9;
	int delta = 2;
	for (int i = 0; i < total_community.size(); i++) {
		double prob = C / (pow(i + 1, delta));
		if (probability(prob) && total_community[i].first == 0) {
			total_community[i].second.time = 1;
			total_community[i].first = 1;
		}
	}
}
//1
void SI_model_with_search_engine() {
	//first, computing the number of infected host in a community
	//Second. computing the number of infected community in the networks                                             
	fstream file_1;
	file_1.open("C:\\Users\\User\\Desktop\\Code\\Matlab Code\\Final_paper_simulation\\analytical_result_1.txt", ios::out);
	int array[10] = {1034, 786, 747, 534, 333, 224, 168, 150, 61 ,52};
	int number_of_community = 10;
	vector<pair<int, Community>> total_community;
	for (int i = 0; i < number_of_community; i++) {
		total_community.push_back(make_pair(0, Community(array[i])));
	}
	file_1 << 1 << endl;
	int t = 1;
	while (t < 35) {
		cout << t << endl;
		infected_community_process(total_community);
		for (int i = 0; i < total_community.size(); i++) {
			if (total_community[i].first == 1) total_community[i].second.infected_process();
			cout<< i << ":" << (int)total_community[i].second.get_infected() << " ";
		}
		cout << endl;
		int total_infected_node = 0;
		for (int i = 0; i < total_community.size(); i++) {
			total_infected_node += (int)total_community[i].second.get_infected();
		}
		
		cout << total_infected_node << endl;
		file_1 << total_infected_node << endl;
		++t;
	}
	return;
}
*/
//2
double S_node_2 = 0.998, I_node_2 = 0.002, node_2 = 500; //S =  random scanning full network
double alpha_2 = 0.25;
double I_r_2 = 0.002, I_loc_2 = 0, I_nhb_2 = 0;
double beta_r = 0.3, beta_loc = 0.4, beta_nhb = 0.6;
double dth_r = 0.01, dth = 0.01, dth_l = 0.018, dth_p = 0.01;
double number_of_node_2 = 500;
double sum_I() {
	return I_r_2 + I_loc_2 + I_nhb_2;
}
double sus_process_2() {
	float tmp = -beta_r * S_node_2*I_node_2 - beta_loc * S_node_2*I_node_2 - beta_nhb * S_node_2*I_node_2 - dth * S_node_2 + alpha_2 * I_node_2;
	return tmp + S_node_2;
}
double infected_r_process_2() {
	float tmp = beta_r * S_node_2*I_node_2 - alpha_2 * I_r_2 - dth * I_r_2 - dth_r * I_r_2;
	return tmp + I_r_2;
}
double infected_l_process_2() {
	float tmp = beta_loc * S_node_2*I_node_2 - alpha_2 * I_loc_2 - dth * I_loc_2 - dth_l * I_loc_2;
	return tmp + I_loc_2;
}
double infected_p_process_2() {
	float tmp = beta_nhb * S_node_2*I_node_2 - alpha_2 * I_nhb_2 - dth * I_nhb_2 - dth_p * I_nhb_2;
	return tmp + I_nhb_2;
}
void BOT_SIS_model() {
	fstream file_2;
	file_2.open("C:\\Users\\User\\Desktop\\Code\\Matlab Code\\Final_paper_simulation\\analytical_result_2.txt", ios::out);
	//Deployment area = 50->100m^2, transmission range = 10->100m (幾乎是全部)
	int WSN_size = 5;//1,5,10
	int t = 1; 
	file_2 << S_node_2 << " " << I_r_2 << " " << I_loc_2 << " " << I_nhb_2 << endl;
	cout << S_node_2 << " " << I_r_2 << " " << I_loc_2 << " " << I_nhb_2 << endl;
	while (t < 20) {
		vector<float> tmp(4, 0);
		I_node_2 = sum_I();
		cout << I_node_2 << endl;
		tmp[0] = sus_process_2();
		tmp[1] = infected_r_process_2();
		tmp[2] = infected_l_process_2();
		tmp[3] = infected_p_process_2();
		//store
		S_node_2 = tmp[0];
		I_r_2 = tmp[1];
		I_loc_2 = tmp[2];
		I_nhb_2 = tmp[3];
		cout << S_node_2 << " " << I_r_2 << " " << I_loc_2 << " " << I_nhb_2 << endl;
		file_2 << S_node_2 << " " << I_r_2 << " " << I_loc_2 << " " << I_nhb_2 << endl;
		++t;
	}

}
//3 SIIR
//fraction
float S_node = 0.4, I1_node = 0.2, I2_node = 0.3, I12_node = 0.1, R_node = 0.0;
//論文I1,I2比例寫反
float beta_1 = 0.3, beta_2 = 0.4, epsilon = 0.5, v1 = 0.36 , v2 = 0.43, v3 = 0.25;
//v1 = v2 = v3 = 0.0, v1 = 0.36, v2 = 0.43, v3 = 0.25
float sus_process() {
	float tmp = -beta_1 * S_node*(I1_node + I12_node) - beta_2 * S_node*(I2_node + I12_node);
	float res = S_node + tmp;
	return res;
}
float infect_process_1() {
	float tmp = beta_1 * S_node*(I1_node + I12_node) - epsilon * beta_2 * I1_node*(I2_node + I12_node) - v1 * I1_node;
	float res = I1_node + tmp;
	return res;
}
float infect_process_2() {
	float tmp = beta_2 * S_node*(I2_node + I12_node) - epsilon * beta_1 * I2_node*(I1_node + I12_node) - v2 * I2_node;
	float res = I2_node + tmp;
	return res;
}
float infect_process_12() {
	float tmp = epsilon * beta_2 * I1_node*(I2_node + I12_node) + epsilon * beta_1 * I2_node*(I1_node + I12_node) - v3 * I12_node;
	float res = I12_node + tmp;
	return res;
}
float recover_process() {
	float tmp = v1 * I1_node + v2 * I2_node + v3 * I12_node;
	float res = R_node + tmp;
	return res;
}
void Bi_Virus_model() {
	fstream file_3;
	//file_3.open("C:\\Users\\User\\Desktop\\Code\\Matlab Code\\Final_paper_simulation\\analytical_result_3_without_impulse.txt", ios::out);
	file_3.open("C:\\Users\\User\\Desktop\\Code\\Matlab Code\\Final_paper_simulation\\analytical_result_3_with_impulse.txt", ios::out);
	int t = 1;
	file_3 << S_node << " " << I1_node << " " << I2_node << " " << I12_node << " " << R_node << endl;
	while (t < 40) {
		vector<float> tmp(5, 0);
		tmp[0] = sus_process();
		tmp[1] = infect_process_1();
		tmp[2] = infect_process_2();
		tmp[3] = infect_process_12();
		tmp[4] = recover_process();
		//store
		S_node = tmp[0], I1_node = tmp[1], I2_node = tmp[2], I12_node = tmp[3], R_node = tmp[4];
		cout << t << ": " << S_node << " " << I1_node << " " << I2_node << " " << I12_node << " " << R_node << endl;
		file_3 << S_node << " " << I1_node << " " << I2_node << " " << I12_node << " " << R_node << endl;
		t++;
	}
}

//4 SNIRD model, simple assunming fully connected
float S_node_4 = 0.8, N_node_4 = 0.0, I_node_4 = 0.2, R_node_4 = 0.0, D_node_4 = 0.0;
float q_RS = 0.25, q_RD = 0.0125;
float q_SN = 0.4, q_SR = 0.15, q_SD = 0.0125;
float q_NR = 0.1, q_NI = 0.2, q_ND = 0.0125;
float q_IR = 0.1, q_ID = 0.05;
float omega_4, mu = 0.01;
int k = 20; float alpha_k = 1, delta_k;//delta_k = k, C, function, k = 20 結果比較準,如果k都相同
void omega_function() {
	omega_4 = (1 / (k*alpha_k)) * alpha_k * delta_k * I_node_4;
	return;
}
void delta_k_function() {
	delta_k = (5 * (float)pow(k, 0.5)) / (1 + 1 * (float)pow(k, 0.5));
	return;
}
float S_process_4() {
	float tmp = S_node_4 + q_RS * R_node_4 - q_SN * omega_4*S_node_4 - q_SR * S_node_4 - q_SD * S_node_4;
	return tmp;
}
float N_process_4() {
	float tmp = N_node_4 + q_SN * omega_4*S_node_4 - q_NI * N_node_4 - q_ND * N_node_4;
	return tmp;
}
float I_process_4() {
	float tmp = I_node_4 + q_NI * N_node_4 - q_IR * I_node_4 - q_ID * I_node_4;
	return tmp;
}
float R_process_4() {
	float tmp = R_node_4 + mu + q_SR * S_node_4 + q_NR * N_node_4 + q_IR * I_node_4 - q_RS * R_node_4 - q_RD * R_node_4;
	return tmp;
}
float D_process_4() {
	float tmp = D_node_4 + q_SD * S_node_4 + q_ND * N_node_4 + q_ID * I_node_4 + q_RD * R_node_4 - mu;
	return tmp;
}
void SNIRD_model() {
	fstream file_4;
	file_4.open("C:\\Users\\User\\Desktop\\Code\\Matlab Code\\Final_paper_simulation\\analytical_result_4.txt", ios::out);
	file_4  << S_node_4 << " " << N_node_4 << " " << I_node_4 << " " << R_node_4 << " " << D_node_4 << endl;
	int t = 1;
	delta_k_function();
	while (t < 100) {
		vector<float> temp(5, 0);
		omega_function();
		temp[0] = S_process_4(), temp[1] = N_process_4(), temp[2] = I_process_4();
		temp[3] = R_process_4(), temp[4] = D_process_4();
		//updating
		S_node_4 = temp[0], N_node_4 = temp[1], I_node_4 = temp[2];
		R_node_4 = temp[3], D_node_4 = temp[4];
		cout << t << ": "  << S_node_4 << " " << N_node_4 << " " << I_node_4 << " ";
		cout << R_node_4 << " " << D_node_4 << endl;
		file_4 << S_node_4 << " " << N_node_4 << " " << I_node_4 << " " << R_node_4 << " " << D_node_4 << endl;
		t++;
	}
}

//5 SEIQRD model
float S_node_5 = 1-0.05, I_node_5 = 0.05, R_node_5 = 0, S1_node_5 = 0, I1_node_5 = 0, R1_node_5 = 0;
float E_node_5 = 0, D_node_5 = 0 , Q_node_5 = 0;
float p_Nw = 0.1, p_Ns = 0.02, p_Nr = 0.02, p_Nq = 0.2, p_Mt = 0.6, p_Md = 0.04;//pNq = 0.9 pNs = 0.4
float a, b, alpha = 0.3, beta = 0.3, gamma_5 = 0.6, xi = 0.3, p = 0.3;//3個參數自行假設
void a_function() {
	float omega_1 = (1 - I_node_5) / (3.14 * 10000);
	a = 2 * sqrt(omega_1*3.14) * 50;
	return;
}
void b_function() {
	float omega_1 = (1 - I_node_5) / (3.14 * 10000);
	b = -omega_1*3.14*2500;
	return;
}
float S_process_5() {
	float tmp = p_Nw * S1_node_5 - p_Ns * S_node_5 - beta * p_Nr*S_node_5 - alpha * p_Mt*(a*sqrt(I_node_5) + b)*S_node_5*xi - p * S_node_5;
	float res = S_node_5 + tmp;
	return res;
}
float S1_process_5() {
	float tmp = p_Ns * S_node_5 - p_Nw * S1_node_5;
	float res = S1_node_5 + tmp;
	return res;
}
float I_process_5() {
	float tmp = gamma_5 * E_node_5 - (p_Md+p_Ns)*I_node_5 + p_Nw * I1_node_5 - p_Nq * I_node_5 - beta * p_Nr*I_node_5;
	float res = I_node_5 + tmp;
	return res;
}
float I1_process_5() {
	float tmp = p_Ns * I_node_5 - p_Nw * I1_node_5;
	float res = I1_node_5 + tmp;
	return res;
}
float R_process_5() {
	float tmp = p_Nw * R1_node_5 - p_Ns * I1_node_5 + beta * p_Nr*(S_node_5 + I_node_5 + Q_node_5);
	float res = R_node_5 + tmp;
	return res;
}
float R1_process_5() {
	float tmp = p_Ns * R_node_5 - p_Nw * R1_node_5;
	float res = R1_node_5 + tmp;
	return res;
}
float D_process_5() {
	float tmp = p * S_node_5 + p * R_node_5 + p * E_node_5 + p_Md * I_node_5;
	float res = D_node_5 + tmp;
	return res;
}
float Q_process_5() {
	float tmp = p_Nq * I_node_5 - beta * p_Nr*Q_node_5;
	float res = Q_node_5 + tmp;
	return res;
}
float E_process_5() {
	float tmp = alpha * p_Mt*(a*sqrt(I_node_5) + b)*S_node_5 - p * E_node_5 - gamma_5 * E_node_5;
	float res = E_node_5 + tmp;
	return res;
} 
void SEIQRD_model() {
	fstream file_5;
	file_5.open("C:\\Users\\User\\Desktop\\Code\\Matlab Code\\Final_paper_simulation\\analytical_result_5.txt", ios::out);
	file_5 << S_node_5 << " " << I_node_5 << " " << R_node_5 << " " << D_node_5 << endl;
	int t = 1;
	while (t < 20) {
		a_function(), b_function();
		//cout << a << " " << b << endl;
		vector<float> temp(9, 0.0);
		temp[0] = S_process_5();
		temp[1] = S1_process_5();
		temp[2] = I_process_5();
		temp[3] = I1_process_5();
		temp[4] = R_process_5();
		temp[5] = R1_process_5();
		temp[6] = D_process_5();
		temp[7] = Q_process_5();
		temp[8] = E_process_5();
		//store
		S_node_5 = temp[0];
		S1_node_5 = temp[1];
		I_node_5 = temp[2];
		I1_node_5 = temp[3];
		R_node_5 = temp[4];
		R1_node_5 = temp[5];
		D_node_5 = temp[6];
		Q_node_5 = temp[7];
		E_node_5 = temp[8];
		cout << S_node_5 << " " << I_node_5 << " " << R_node_5 << " " << D_node_5 << endl;
		file_5 << S_node_5 << " " << I_node_5 << " " << R_node_5 << " " << D_node_5  << endl;
		++t;
	}
}


//NEW comparing papers
//Idea: Comparing only long, only short, hybrid
/*
Direct:
1.Web malware spread modelling and optimal control strategies %liu2017web test ok
2.Virus Propagation and Patch Distribution in Multiplex Networks: Modeling, Analysis, and Optimal Allocation
*/
/*
Short:
3.SNIRD Disclosing Rules of Malware Spread in Heterogeneous Wireless Sensor Networks(without mobility)
4.An Epidemiology-Based Model for Disclosing Dynamics of Malware Propagation in Heterogeneous and Mobile WSNs(with mobility)
*/
/*
Hybrid:
5.Modelling the Spread of Botnet Malware in IoT-Based Wireless Sensor Networks(without mobility)
*/

class liu2017web{//SDIR (susceptible,delitescent,infected,recovered)
private:
	int n, RunT;
	double infected_start = 0.1;
	//parameters
	double lambda = 0.00005, mu = 0.5, epsilon = 0.2, gamma = 0.2, d = 0.01, b = 100, zeta = 0.1;//論文數據
public:
	liu2017web(int num, int time) {
		n = num, RunT = time;
	}
	vector<vector<double>> Process() {
		vector<vector<double>> res(RunT, vector<double>(4, 0.0));
		res.resize(RunT, vector<double>(4, 0.0));
		vector<double> state(4, 0.0);
		//initial
		state[0] = n - infected_start * n, state[1] = 0, state[2] = infected_start * n, state[3] = 0;
		cout << 0 << ": ";
		Printing(state, res[0]);
		int t = 1;
		while (t < RunT) {
			vector<double> tmp(4,0.0);
			Computing(tmp, state);
			Update(tmp, state);
			cout << t << ": ";
			Printing(state, res[t]);
			++t;
		}
		return res;
	}
	//0S 1D 2I 3R
	void Computing(vector<double>& tmp, vector<double> state) {
		tmp[0] = b - lambda*state[0]*state[2]+zeta*state[3]-d*state[0];
		tmp[1] = lambda * state[0] * state[2]-mu*state[1]-epsilon*state[1]-d*state[1];
		tmp[2] = epsilon * state[1] -gamma*state[2]-d*state[2];
		tmp[3] = mu * state[1]+ gamma * state[2]- zeta * state[3]-d*state[3];
	}
	void Update(vector<double> tmp, vector<double>& state) {
		for (int i = 0; i < state.size(); i++) {
			state[i] = (state[i] + tmp[i] < 0) ? (double)0.0 : state[i] + tmp[i];
		}
	}
	void Printing(vector<double> state, vector<double>& res) {
		for (int i = 0; i < state.size(); i++) {
			res[i] = state[i]/(double)n;
			cout << res[i] << " ";
		}
		cout << endl;
	}
};

class zhao2019virus {//Markov chain 建好尚未測試
private:
	int n, RunT;
	double infected_start = 0.1;
	vector<vector<double>> prob_matrix;//store probility
	double m = 0.03, beta_v = 0.3, beta_p = 0.3 , delta_v = 0.2, delta_p = 0.4;
public:
	zhao2019virus(int num, int time) {
		n = num, RunT = time;
	}
	vector<vector<double>> Process(Social_network S) {//need adj network
		vector<vector<double>> res(RunT, vector<double>(4, 0.0));
		prob_matrix.resize(n, vector<double>(4, 0.0));
		vector<double> state(n, 0.0);
		//initial 初始化感染點
		for (int i = 0; i < prob_matrix.size(); i++) {
			prob_matrix[i][0] = 1.0;
		}
		for (int i = 0; i < n*infected_start; i++) {
			prob_matrix[i][0] = 0.0;
			prob_matrix[i][1] = 1.0;
			state[i] = 1.0;
		}

		cout << 0 << ": ";
		Printing(state, res[0]);
		int t = 1;
		while (t < RunT) {
			vector<double> tmp(n, 0.0);
			//因為是markov chain 每個點要各自判斷
			ComputingProb(prob_matrix, S);//compute t + 1, 0 -> 1, 1 -> 2
			Determine_state(state, prob_matrix);
			cout << t << ": ";
			Printing(state, res[t]);
			++t;
		}
		return res;
	}
	void ComputingProb(vector<vector<double>>& prob_matrix, Social_network S) {
		vector<vector<double>> tmp(n, vector<double>(4, 0.0));
		for (int i = 0; i < prob_matrix.size(); i++) {
			double q_i = 1.0, r_i = 1.0;
			for (int j = 0; j < n; j++) {
				if (S.adj_matrix[i][j] == 1 && i != j) {
					q_i *= 1.0 - prob_matrix[j][1] * beta_v;
					r_i *= 1.0 - prob_matrix[j][2] * beta_p;
				}
			}
			tmp[i][0] = prob_matrix[i][0] * (1 - m)*r_i*q_i + prob_matrix[i][1] * (1 - m)*r_i*delta_v;
			tmp[i][1] = prob_matrix[i][0] * (1 - m)*r_i*(1 - q_i) + prob_matrix[i][1] * (1 - m)*r_i*(1 - delta_v);
			tmp[i][2] = prob_matrix[i][0]*(1-((1 - m)*r_i))+ prob_matrix[i][1]*(1-((1 - m)*r_i))+ prob_matrix[i][2]*(1-delta_v);
			tmp[i][3] = prob_matrix[i][3] + prob_matrix[i][2] * delta_p;
		}
		//Updating
		for (int i = 0; i < prob_matrix.size(); i++) {
			prob_matrix[i][0] = tmp[i][0];
			prob_matrix[i][1] = tmp[i][1];
			prob_matrix[i][2] = tmp[i][2];
			prob_matrix[i][3] = tmp[i][3];
		}
	}
	// 0 US 1 UI 2 PV 3 DV
	void Determine_state(vector<double> & state, vector<vector<double>> prob_matrix) {
		for (int i = 0; i < n; i++) {
			double r = (double)rand() / (RAND_MAX);
			if (r <= prob_matrix[i][0]) state[i] = 0;
			r -= prob_matrix[i][0];
			if (r <= prob_matrix[i][1]) state[i] = 1;
			r -= prob_matrix[i][0];
			if (r <= prob_matrix[i][2]) state[i] = 2;
			r -= prob_matrix[i][0];
			if (r <= prob_matrix[i][3]) state[i] = 3;
		}
	}
	void Update(vector<double> tmp, vector<double>& state) {//Markov 
		for (int i = 0; i < state.size(); i++) {
			state[i] = tmp[i];
		}
	}
	void Printing(vector<double> state, vector<double>& res) {
		vector<double> temp(4,0.0);
		for (int i = 0; i < state.size(); i++) {
			temp[state[i]] += 1 / (double)n;
		}
		for (int i = 0; i < temp.size(); i++) {
			cout << temp[i] << " ";
		}
		cout << endl;
		res = temp;
	}
};

class shen2019snird {// 修改上面的 建好
private:
	int n, RunT;
	double q_RS = 0.25, q_RD = 0.0125;
	double q_SN = 0.4, q_SR = 0.15, q_SD = 0.0125;
	double q_NR = 0.1, q_NI = 0.2, q_ND = 0.0125;
	double q_IR = 0.1, q_ID = 0.05;
	double omega_k, mu = 0.01;
	double alpha_k = 0.333;//delta_k = k, C, function, k = 20 結果比較準,如果k都相同
	//degree 設50 20 10
	double infected_start = 0.1;
public:
	shen2019snird(int num, int time) {
		n = num, RunT = time;
	}
	vector<vector<double>> Process() {
		vector<vector<double>> res(RunT, vector<double>(5, 0.0));
		vector<double> state(5, 0.0);
		//initial
		state[0] = 1 - infected_start , state[1] = 0, state[2] = infected_start, state[3] = 0, state[4] = 0;
		cout << 0 << ": ";
		Printing(state, res[0]);
		int t = 1;
		while (t < RunT) {
			vector<double> tmp(5, 0.0);
			omega_function(state);

			Computing(tmp, state);
			Update(tmp, state);
			cout << t << ": ";
			Printing(state, res[t]);
			++t;
		}
		return res;
	}
	//0S 1N 2I 3R 4D
	void Computing(vector<double>& tmp, vector<double> state) {
		tmp[0] = q_RS * state[3] - q_SN * omega_4*state[0] - q_SR * state[0] - q_SD * state[0];
		tmp[1] = q_SN * omega_4*state[0] - q_NI * state[1] - q_ND * state[1];
		tmp[2] = q_NI * state[1] - q_IR * state[2] - q_ID * state[2];
		tmp[3] = mu + q_SR * state[0] + q_NR * state[1] + q_IR * state[2] - q_RS * state[3] - q_RD * state[3];
		tmp[4] = q_SD * state[0] + q_ND * state[1] + q_ID * state[2] + q_RD * state[3] - mu;
	}
	void Update(vector<double> tmp, vector<double>& state) {
		for (int i = 0; i < state.size(); i++) {
			state[i] = (state[i] + tmp[i] < 0) ? (double)0.0 : state[i] + tmp[i];
		}
	}
	void Printing(vector<double> state, vector<double>& res) {
		for (int i = 0; i < state.size(); i++) {
			res[i] = state[i];
			cout << res[i] << " ";
		}
		cout << endl;
	}
	void omega_function(vector<double> state) {
		double m = 50 * alpha_k + 20 * alpha_k + 10 * alpha_k;
		double tmp = 0.0;
		omega_k = (1 / m)*(alpha_k*delta_k(50)*state[2] + alpha_k * delta_k(20)*state[2] + alpha_k * delta_k(10)*state[2]);
		return;
	}
	double delta_k(int k) {
		double delta_k = (5 * (double)pow(k, 0.5)) / (1 + 1 * (double)pow(k, 0.5));
		return delta_k;
	}
};

class shen2020epidemiology {//VCQPS
private:
	int n, RunT;
	double infected_start = 0.1;

	double alpha;
	double phi_QV, phi_QP, phi_QS;
	double phi_PV, phi_PS;
	double phi_VC, phi_VP, phi_VS;
	double phi_CQ, phi_CP, phi_CS;
	double varphi, varsigma, R, v_mean, chi_n;
public:
	shen2020epidemiology(int num, int time) {
		n = num, RunT = time;
	}
	//0V 1C 2Q 3P 4S
	vector<vector<double>> Process() {
		vector<vector<double>> res(RunT, vector<double>(5, 0.0));
		vector<double> state(5, 0.0);
		//initial
		state[0] = 1 - infected_start, state[1] = infected_start;
		cout << 0 << ": ";
		Printing(state, res[0]);
		int t = 1;
		while (t < RunT) {
			vector<double> tmp(5, 0.0);
			//還要計算

			Computing(tmp, state);
			Update(tmp, state);
			cout << t << ": ";
			Printing(state, res[t]);
			++t;
		}
		return res;
	}
	//0V 1C 2Q 3P 4S
	void Computing(vector<double>& tmp, vector<double> state) {
		tmp[0] = phi_QV*state[2]+phi_PV*state[3]-phi_VC*varsigma*varphi*state[0]*chi_n-phi_VP*state[0]-phi_VS*state[0];
		tmp[1] = phi_VC * varsigma*varphi*state[0] - phi_CQ * state[1] - phi_CP * state[1] - phi_CS * state[1];
		tmp[2] = phi_CQ * state[1] - phi_QV * state[2] - phi_QP * state[2] - phi_QS * state[2];
		tmp[3] = alpha + phi_VP * state[0] + phi_CP * state[1] + phi_QP * state[2] - phi_PS * state[3];
		tmp[4] = phi_VS * state[0] + phi_CS * state[1] + phi_QS * state[2] + phi_PS * state[3] - alpha;
	}
	void Update(vector<double> tmp, vector<double>& state) {
		for (int i = 0; i < state.size(); i++) {
			state[i] = (state[i] + tmp[i] < 0) ? (double)0.0 : state[i] + tmp[i];
		}
	}
	void Printing(vector<double> state, vector<double>& res) {
		for (int i = 0; i < state.size(); i++) {
			res[i] = state[i];
			cout << res[i] << " ";
		}
		cout << endl;
	}
};
class acarali2019modelling {//SIIIS 問題1:S感覺都一樣 建好
private:
	int n, RunT, WSN = 5;
	double alpha = 0.25;
	double I_r_2 = 0.002, I_loc_2 = 0, I_nhb_2 = 0;
	double beta_random = 0.3 , beta_local = 0.4 , beta_p2p = 0.6;
	double dth_random  = 0.01, dth_local = 0.018, dth_p2p = 0.01, dth = 0.01;
	double infected_start = 0.1;
public:
	acarali2019modelling(int num, int time) {
		n = num, RunT = time;
	}
	vector<vector<double>> Process() {
		vector<vector<double>> res(RunT, vector<double>(2, 0.0));//0 S 1 I
		vector<double> state(4, 0.0);//比例
		//0S 1I_random 2I_local 3I_p2p
		//一開始透過random擴散
		state[0] = 1 - infected_start, state[1] = infected_start;
		
		cout << 0 << ": ";
		Printing(state, res[0]);
		int t = 1;
		while (t < RunT) {
			double I = 0.0;
			I = state[1] + state[2] + state[3];
			vector<double> tmp(4, 0.0);
			Computing(tmp, state, I);
			Update(tmp, state);
			cout << t << ": ";
			Printing(state, res[t]);
			++t;
		}
		return res;
	}
	void Computing(vector<double>& tmp, vector<double> state,double I) {// S = S_loc = S_nhb
		tmp[0] = -beta_random * state[0] * I - beta_local * state[0] * I - beta_p2p * state[0] * I - dth * state[0] + alpha * I;
		tmp[1] = beta_random * state[0] * I - alpha * state[1] - dth * state[1] - dth_random * state[1];
		tmp[2] = beta_local * state[0] * I - alpha * state[2] - dth * state[2] - dth_local * state[2];
		tmp[3] = beta_p2p * state[0] * I - alpha * state[3] - dth * state[3] - dth_p2p * state[3];
	}
	void Update(vector<double> tmp, vector<double>& state) {
		for (int i = 0; i < state.size(); i++) {
			state[i] = (state[i] + tmp[i] < 0) ? (double)0.0 : state[i] + tmp[i];
		}
	}
	void Printing(vector<double> state, vector<double>& res) {
		res[0] = state[0];
		res[1] = (state[1] + state[2] + state[3]);
		
		for(int i = 0; i < res.size(); i++) cout << res[i] << " ";//S I
		cout << endl;
	}

};