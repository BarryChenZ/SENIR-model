#pragma once
//3 types epidemic model from other papers

//1. A Two layer Model of Malware Propagation in a Search Engine Context
//2. Modelling the Spread of Botnet Malware in IoT-Based Wireless Sensor NetworksModelling the Spre
//3. Optimal Impulse Control of Bi-Virus SIR Epidemics with Application to Heterogeneous Internet of Things
//4. SNIRD  Disclosing Rules of Malware Spread in Heter wireless Sensor networks
//5. A cloud-assisted malware detection and suppression framework for wireless multimedia system in IoT based on dynamic differential game

#include <iostream>
#include <vector>
#include <math.h>
#include <time.h>
#include <fstream>

using namespace std;
//1
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
	//Deployment area = 50->100m^2, transmission range = 10->100m (碭琌场)
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
//阶ゅI1,I2ゑㄒ糶は
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
float omega, mu = 0.01;
int k = 20; float alpha_k = 1, delta_k;//delta_k = k, C, function, k = 20 挡狦ゑ耕非,狦k常
void omega_function() {
	omega = (1 / (k*alpha_k)) * alpha_k * delta_k * I_node_4;
	return;
}
void delta_k_function() {
	delta_k = (5 * (float)pow(k, 0.5)) / (1 + 1 * (float)pow(k, 0.5));
	return;
}
float S_process_4() {
	float tmp = S_node_4 + q_RS * R_node_4 - q_SN * omega*S_node_4 - q_SR * S_node_4 - q_SD * S_node_4;
	return tmp;
}
float N_process_4() {
	float tmp = N_node_4 + q_SN * omega*S_node_4 - q_NI * N_node_4 - q_ND * N_node_4;
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
float a, b, alpha = 0.3, beta = 0.3, gamma = 0.6, xi = 0.3, p = 0.3;//3把计︽安砞
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
	float tmp = gamma * E_node_5 - (p_Md+p_Ns)*I_node_5 + p_Nw * I1_node_5 - p_Nq * I_node_5 - beta * p_Nr*I_node_5;
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
	float tmp = alpha * p_Mt*(a*sqrt(I_node_5) + b)*S_node_5 - p * E_node_5 - gamma * E_node_5;
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
