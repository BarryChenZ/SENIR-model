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
#include "Mobility.h"
#include <vector>
#include <math.h>
using namespace std;

class Node {
private:
	int number;//編號
	int status = 0;
	//transmission
	vector<double> rate;//參數可以是rate也可以分開假設//暫定
	double degree_P;
	double degree_S;
	//mobility
	double x_pos, y_pos;
public:
	vector<int> neighbor_set;
	double infected_rate = 0.32;
	double received_rate = 0.228;
	double exposed_rate;
	double iNsidious_rate = 0.8;
	double recovered_rate = 0.083;
	double death_rate = 0.2;
	double wake_up_rate;
	double lose_immunity_rate;

	Node(int ID) {
		number = ID;
	};
	double countNeigh(vector<vector<int>> adj) {//on the social network
		double tmp = 0;
		for (int i = 0; i < adj.size; i++) {
			if (adj[number][i] == 1 && i != number) tmp++;
		}
		return tmp;
	}
	void set_position(double x, double y) {
		x_pos = x, y_pos = y;
	}
	double get_position_x() {
		return x_pos;
	}
	double get_position_y() {
		return y_pos;
	}
	int get_state() {
		return status;
	}
	void change_state(int state) {
		status = state;
		return;
	}
};

class Physical_network {
	//can be one or >1 group
private:
	int nums; //可以有部分或者很多點
	double range;
	double max_x = 1000, max_y = 1000;
public:
	//type 1 : based on range
	//type 2 : based on group
	Physical_network(int NUM, vector<Node> NODES) {
		nums = NUM;
		for (int i = 0; i < nums; i++) {
			double tmp1 = Uniform_distri(0.0, max_x);
			double tmp2 = Uniform_distri(0.0, max_y);
			NODES[i].set_position(tmp1, tmp2);
		}
	}
	double get_degree_s(int ID, vector<Node> NODES) {
		double res = 0;
		for (int i = 0; i < nums; i++) {
			NODES[i].neighbor_set.clear();
			double tmp = pow(NODES[i].get_position_x() - NODES[ID].get_position_x, 2) + pow(NODES[i].get_position_y() - NODES[ID].get_position_y, 2);//distance
			if (sqrt(tmp) <= range && i != ID) {
				NODES[i].neighbor_set.push_back(i);
				res++;
			}
		}
		return res;
	}
};

class Social_network {
private:
	int nodes;
	double average_degree;
public:
	vector<vector<int>> adj_matrix;

	Social_network(int num){
		nodes = num;
		initial();
	}
	void initial() {
		srand(time(NULL));
		//average_degree = AVG_D;
		adj_matrix.resize(nodes, vector<int>(nodes));
		vector<double> degree(nodes, 0); 
		double total_link = 0;//total > start*start
		//先創一個full connected small network
		int start = 10;
		for (int i = 0; i < start; i++) {
			for (int j = 0; j < start; j++) {
				adj_matrix[i][j] = 1;
			}
			degree[i] = start;
			total_link += start;
		}
		total_link -= start * start;
		//分配link by BA model(Scare-free)
		while (start <= nodes) {
			double tmp = 0;
			for (int i = 0; i < start; i++) {
				double prob = degree[i]/total_link;
				if (isValid(prob)) {
					adj_matrix[start][i] = 1, adj_matrix[i][start] = 1;
					degree[i]++, degree[start]++;
					tmp++;
				}
				total_link += tmp * 2;
			}
			start++;
		}
		return;
	}
	bool isValid(double prob) {
		double x = (double)rand() / (RAND_MAX + 1.0);
		if (x > prob) return true;
		return false;
	}
	double get_degree_s(int ID) {
		double res = 0;
		for (int i = 0; i < nodes; i++) {
			if (adj_matrix[ID][i] == 1) res++;
		}
		return res;
	}
};