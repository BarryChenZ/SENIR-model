#pragma once
// 1. Implement human mobility for discrete time
// 2. Caculating effective transmission node based on the range

#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <random>
#include "Distribution.h"

using namespace std;

#define PI acos(-1)

double Uniform_distri(double min, double max) {
	double x = (double)rand() / (RAND_MAX + 1.0);
	double res = x * (max - min) + min;
	return res;
}
double PowerLaw_distri(double ALPHA, double min, double max) {
	double x = (double)rand() / (RAND_MAX + 1.0);
	double res = pow(((pow(max, (ALPHA + 1.0)) - 1.0) * x + 1.0), 1.0 / (ALPHA + 1.0));
	return res;
}

class Levy_Walk {
private:
	//only exist in this class
	int num_of_nodes;
	int max_x = 100, max_y = 100;
	int ndim = 2; int dimension[2] = { 100,100 };
	double FL_EXP = 2.6, FL_MAX = 50.0, WT_EXP = -1.8, WT_MAX = 100.0;
public:
	Levy_Walk(int number) { num_of_nodes = number; }
	double FL_distri() {
		return PowerLaw_distri(FL_EXP, 1.0, FL_MAX);
	}
	double Velocity_distri(double value) {
		return sqrt(value) / 10;
	}
	void check(int border_policy, vector<double>& x_pos, vector<double>& y_pos, vector<double>& movement_x, vector<double>& movement_y) {
		if (border_policy) reflect(x_pos, y_pos, movement_x, movement_y);
		else               wrap(x_pos, y_pos, movement_x, movement_y);
	}

	void reflect(vector<double>& x_pos, vector<double>& y_pos, vector<double>& movement_x, vector<double>& movement_y) {
		for (int i = 0; i < num_of_nodes; i++) {
			if (x_pos[i] < 0) {
				x_pos[i] = -x_pos[i];
				movement_x[i] = -movement_x[i];
			}
			else if (x_pos[i] > max_x) {
				x_pos[i] = 2 * max_x - x_pos[i];
				movement_x[i] = -movement_x[i];
			}
			if (y_pos[i] < 0) {
				y_pos[i] = -y_pos[i];
				movement_y[i] = -movement_y[i];
			}
			else if (y_pos[i] > max_y) {
				y_pos[i] = 2 * max_y - y_pos[i];
				movement_y[i] = -movement_y[i];
			}
		}
	}

	void wrap(vector<double>& x_pos, vector<double>& y_pos, vector<double>& movement_x, vector<double>& movement_y) {
		for (int i = 0; i < num_of_nodes; i++) {
			if (x_pos[i] < 0) x_pos[i] += max_x;
			else if (x_pos[i] > max_x) x_pos[i] -= max_x;
			if (y_pos[i] < 0) y_pos[i] += max_y;
			else if (y_pos[i] > max_y) y_pos[i] -= max_y;
		}
	}

	void StochasticWalk(int policy, vector<double>& x_pos, vector<double>& y_pos) {
		//assign flight lengths, velocities
		vector<double> fl(num_of_nodes, 0);
		vector<double> velocity(num_of_nodes, 0);
		for (int i = 0; i < num_of_nodes; i++) {
			//x_pos[i] = Uniform_distri(0.0, dimension[0]); // (0. max_x)
			//y_pos[i] = Uniform_distri(0.0, dimension[1]); // (0, max_y)
			fl[i] = FL_distri();
			velocity[i] = Velocity_distri(fl[i]);
		}
		//assign nodes' movements
		vector<double> direction_x(num_of_nodes, 0);
		vector<double> direction_y(num_of_nodes, 0);
		vector<double> movement_x(num_of_nodes, 0);
		vector<double> movement_y(num_of_nodes, 0);
		for (int i = 0; i < num_of_nodes; i++) {
			direction_x[i] = Uniform_distri(0.0, 1.0) - 0.5;
			direction_y[i] = Uniform_distri(0.0, 1.0) - 0.5;
			//calculating norm
			double norm = sqrt(pow(direction_x[i], 2) + pow(direction_y[i], 2));
			movement_x[i] = direction_x[i] / norm * velocity[i];
			movement_y[i] = direction_y[i] / norm * velocity[i];
		}
		//Starting with no waiting time
		vector<double> waiting(num_of_nodes, 0);
		int t = 0;
		while (t < 100) { //true
			vector<double> arrived;
			for (int i = 0; i < num_of_nodes; i++) {
				x_pos[i] += movement_x[i];
				y_pos[i] += movement_y[i];
				fl[i] -= velocity[i];
				if (velocity[i] > 0 && fl[i] <= 0.0) arrived.push_back(i);
			}
			//step back (超過飛行距離)
			if (arrived.size() > 0) {
				for (int i = 0; i < arrived.size(); i++) {
					double diff = fl[arrived[i]] / velocity[i];
					x_pos[arrived[i]] += diff * movement_x[arrived[i]];
					y_pos[arrived[i]] += diff * movement_y[arrived[i]];
				}
			}
			//apply border
			check(policy, x_pos, y_pos, movement_x, movement_y);
			//updating the moving nodes,if nodes complete fl then apply new fl, v, direc, movement
			if (arrived.size() > 0) {
				for (int i = 0; i < arrived.size(); i++) {
					fl[arrived[i]] = FL_distri();
					velocity[arrived[i]] = Velocity_distri(fl[arrived[i]]);
					double v = velocity[arrived[i]];
					direction_x[arrived[i]] = Uniform_distri(0.0, 1.0) - 0.5;
					direction_y[arrived[i]] = Uniform_distri(0.0, 1.0) - 0.5;
					//calculating norm
					double norm = sqrt(pow(direction_x[arrived[i]], 2) + pow(direction_y[arrived[i]], 2));
					movement_x[arrived[i]] = direction_x[arrived[i]] / norm * velocity[arrived[i]];
					movement_y[arrived[i]] = direction_y[arrived[i]] / norm * velocity[arrived[i]];
				}
			}
			//output the position(pos_x,pos_y) at t;
			cout << "Time: " << t << endl;
			for (int i = 0; i < num_of_nodes; i++) {
				cout << i << " : " << x_pos[i] << " " << y_pos[i] << endl;
			}
			t++;
		}
	}
};

class Gauss_Markov {
private:
	int number;
	double velocity;
	double direction;
	double chi_v = 0.85, mu_v = 10.0;//chi = alpha,mu = mean
	double chi_d = 0.85, mu_d = (1/3)* PI;//degree 60
	double chi_s = 1;
	double mu = 0.0, sigma = 5.0; // sigma越大浮動越大
	int max_x = 1000, max_y = 1000;


public:
	Gauss_Markov(int ID) {
		number = ID;
		velocity = mu_v;
	}
	Gauss_Markov() {}

	double change_velocity(double v) {
		std::random_device rd;
		std::default_random_engine generator(rd());
		std::normal_distribution<double> distribution(mu, sigma);
		//t-1 to t
		double ra = distribution(generator);
		v = chi_v * v + (1 - chi_v)*mu_v + chi_s * (ra*sqrt(1 - pow(chi_v, 2)));
		return v;
	}
	void initial_v(vector<double> v) {
		for (int i = 0; i < v.size(); i++) {
			v[i] = mu_v;
		}
		return;
	}
	void Set_v() {
		std::random_device rd;
		std::default_random_engine generator(rd());
		std::normal_distribution<double> distribution(mu, sigma);
		//t-1 to t
		double ra = distribution(generator);
		//cout << ra << endl;
		velocity = chi_v * velocity + (1 - chi_v)*mu_v + chi_s * (ra*sqrt(1 - pow(chi_v, 2)));
	}
	void Set_d() {
		std::random_device rd;
		std::default_random_engine generator(rd());
		std::normal_distribution<double> distribution(mu, sigma);
		direction = chi_d * direction + (1 - chi_d)*mu_v + chi_s * (distribution(generator)*sqrt(1 - pow(chi_d, 2)));
	}
	void Process(vector<Node>& NODES) {//Moving one time
		double x = NODES[number].get_position_x() + velocity * cos(direction);
		double y = NODES[number].get_position_y() + direction * sin(direction);
		if (x > max_x || x < 0) x = NODES[number].get_position_x() - velocity * cos(direction);
		if (y > max_y || y < 0) y = NODES[number].get_position_y() - velocity * cos(direction);
		NODES[number].set_position(x, y);
	}
	double get_velocity() {
		return velocity;
	}
};

//Identiyty mobility model
class Identity {
private:
	int number; //for ID
	int positions;// number of positions
	double rho = 0.7, gamma = 0.7;//Test value
	vector<int> times;
	double total_time = 1;
	double visited = 1;
public:
	int At_position;
	Identity(int ID, int n) {//initial
		number = ID;
		positions = n;
		times.resize(positions, 0);
		At_position = rand() % (positions - 0) + 0;
		times[At_position]++;
	}
	void Process() {
		double exploraiton = rho * pow(visited, -1 * gamma);
		//double Pre_return = 1 - exploraiton;
		double tmp = (double)rand() / (RAND_MAX + 1.0);
		if (isValid(exploraiton)) {
			int unvisited = positions - visited;
			double temp1 = 0.0;
			for (int i = 0; i < positions; i++) {
				if (times[i] == 0) {
					temp1 += pow(unvisited, -1);// 1/unvisited 永遠=0
					if (temp1 > tmp) {
						At_position = i;
						times[i]++;
						break;
					}
				}
			}
			visited++;
		}
		else {
			double temp = 0.0;
			for (int i = 0; i < positions; i++) {
				temp += times[i] / total_time;
				//cout << temp << endl;
				if (times[i] > 0 && temp > tmp) {
					At_position = i;
					times[i]++;
					break;
				}
			}
		}
		total_time++;
		return;
	}
	bool isValid(double prob) {
		double x = (double)rand() / (RAND_MAX + 1.0);
		if (x < prob) return true;
		return false;
	}
};