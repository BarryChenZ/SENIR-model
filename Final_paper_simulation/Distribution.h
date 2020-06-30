#pragma once
//Initial node position distribution
// Class node
// Node includes:
// 1. Position (initial, mobility)
// 2. Hetergeneous infected, recovered rate etc.
// 3. State
// 4. Useful information, example degree (Reducing complexity), velocity, direction etc.

// Full network
// Social network (degree)
// 
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <map>
#include <math.h>
//#include "Mobility.h"
using namespace std;

double Uniform_distri(double min, double max);

class Node {
private:
	int number;//絪腹
	int status = 0;
	//transmission
	vector<double> rate;//把计琌rateだ秨安砞//既﹚
	double degree_P;
	//mobility
	double x_pos, y_pos;
public:
	double degree_S;
	vector<int> neighbor_set;
	double infected_rate = 0.04;
	double infected_rate_2 = 0.72;
	double received_rate = 0.02;
	double exposed_rate = 0.015;
	double iNsidious_rate = 0.015;
	double recovered_rate = 0.1;
	double death_rate = 0.05;
	double ex_death_rate = 0.01;
	double wake_up_rate = 0.3;
	double wake_up_rate_r = 0.2;
	double lose_immunity_rate = 0.1;
	double new_rate = 0.01;// new rate = leave rate

	Node(int ID) {
		number = ID;
	};
	double countNeigh(vector<vector<int>>& adj) {//on the social network
		double tmp = 0;
		for (int i = 0; i < adj.size(); i++) {
			if (adj[number][i] == 1 && i != number) tmp++;
		}
		return tmp;
	}
	void set_position(double x, double y) {
		x_pos = x, y_pos = y;
		return;
	}
	double get_position_x() {
		return x_pos;
	}
	double get_position_y() {
		return y_pos;
	}
	void set_degreeP(int degree) {
		degree_P = degree;
	}
	int get_state() {
		return status;
	}
	int get_number() {
		return number;
	}
	void change_state(int state) {
		status = state;
		return;
	}
};

class Physical_network {
	//can be one or >1 group
private:
	int nums; //Τ场だ┪翴
	double range = 100.0;//10 average degree:1, 30 average degree: 5, 50 average degree: 8
	double range_array[3] = {10, 100, 200};
	double max_x = 1000, max_y = 1000;
public:
	//type 1 : based on range
	//type 2 : based on group
	Physical_network(int NUM, vector<Node>& NODES) {
		nums = NUM;
		for (int i = 0; i < nums; i++) {
			double tmp1  = Uniform_distri(0.0, max_x);
			double tmp2  = Uniform_distri(0.0, max_y);
			
			NODES[i].set_position(tmp1, tmp2);
		}
		for (int i = 0; i < nums; i++) {
			NODES[i].set_degreeP(get_degree(i,NODES));
		}
	}
	double get_degree(int ID, vector<Node>& NODES) {
		double res = 0;
		NODES[ID].neighbor_set.clear();
		for (int i = 0; i < nums; i++) {
			double tmp = pow(NODES[i].get_position_x() - NODES[ID].get_position_x(), 2) + pow(NODES[i].get_position_y() - NODES[ID].get_position_y(), 2);//distance
			if (sqrt(tmp) <= range && i != ID) {
				NODES[ID].neighbor_set.push_back(i);
				res++;
			}
		}
		return res;
	}
	double get_range() {//for mathmatics model
		return range;
	}
	void setting_range(int k) {
		range = range_array[k - 1];
	}
	double area() {
		return max_x * max_y;
	}
};

class Social_network {
private:
	int nodes;
	double average_degree;
public:
	vector<vector<int>> adj_matrix;
	map<int, int> degree_count;//参璸ぃdegree翴计
	vector<int> degree_individual;
	Social_network(int num, vector<Node>& NODES){
		nodes = num;
		//initial();

		//fixed degree
		
		vector<int> degree(num, 0);
		for (int i = 0; i < nodes*0.1; i++) degree[i] = 50;
		for (int i = nodes * 0.1; i < nodes*0.4; i++) degree[i] = 20;
		for (int i = nodes * 0.4; i < nodes; i++) degree[i] = 10;
		random_shuffle(degree.begin(), degree.end());
		for (int i = 0; i < nodes; i++) {
			NODES[i].degree_S = degree[i];
		}
		//initial_fixedDegree(degree);
		
	}
	void initial() {
		//srand((unsigned)time(NULL));
		//average_degree = AVG_D;
		adj_matrix.resize(nodes, vector<int>(nodes));
		vector<double> degree(nodes, 0); 
		double total_link = 0.0;//total > start*start
		//承full connected small network
		int start = 10;
		for (int i = 0; i < start; i++) {
			for (int j = 0; j < start; j++) {
				adj_matrix[i][j] = 1;
			}
			degree[i] = start-1;//ぃ
			total_link += start-1;
		}
		cout << total_link << endl;
		
		//だ皌link by BA model(Scare-free)
		start++;
		while (start < nodes) {
			double tmp = 0.0;
			for (int i = 0; i < start-1; i++) {
				//cout << degree[i] << endl;
				double prob = degree[i]/total_link;
				//cout << prob << endl;;
				if (isValid(prob)) {
					adj_matrix[start][i] = 1, adj_matrix[i][start] = 1;
					degree[i]++;
					tmp++;
				}
			}
			degree[start] = tmp;
			total_link += degree[start] * 2;
			start++;
		}
		cout << total_link << endl;

		for (int i = 0; i < nodes; i++) {
			int temp = get_degree(i);
			degree_count[temp]++;
		}
		return;
	}
	void initial_fixedDegree(vector<int> degree) {
		adj_matrix.resize(nodes, vector<int>(nodes));
		for (int i = 0; i < nodes; i++) {
			for (int j = 0; j < nodes; j++) {
				if (degree[i] > 0 && degree[j] > 0 && i != j) {
					degree[i]--, degree[j]--;
					adj_matrix[i][j] = 1, adj_matrix[j][i] = 1;
				}
			}
		}
	}

	bool isValid(double prob) {
		double x = (double)rand() / (RAND_MAX + 1.0);
		if (x < prob) return true;
		return false;
	}
	int get_degree(int ID) {
		int res = 0;
		for (int i = 0; i < nodes; i++) {
			if (adj_matrix[ID][i] == 1 && ID != i) res++;
		}
		return res;
	}
};