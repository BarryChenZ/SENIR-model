#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include <cctype>

#include "Distribution.h"
#include "ProcessModel.h"
#include "ExperimentRecord.h"
#include "FinalModel.h"

using namespace std;


int times = 1;

struct transaction_ID {
	int ID;
	string tx_muuid;
	
	int data;
	//int unix_timestamp;
	int upload_time_hour;
	int upload_time_minute;
	int upload_time_second;
	double distan; //方便收集限制範圍內的資料

	string rx_muuid_short;
	string tx_muuid_short;
};

string mmuid_to_short(string s) {
	string res;
	bool capital = false;
	for (int i = 0; i < s.size(); i++) {
		if (s[i] >= 65 && s[i] <= 90) {
			capital = true;
			break;
		}
	}
	stringstream ss;
	if (capital != true) {
		//AND
		ss << 'A' << s[0] << s[14] << s[s.size() - 1];
	}
	else {
		//IOS
		ss << 'B' << s[0] << s[14] << s[s.size() - 1];
	}
	res = ss.str();
	return res;
}

bool infectedProcess(double probability) { // N->I
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

string to_lower(string s) {
	for (int i = 0; i < s.size(); i++) {
		s[i] = tolower(s[i]);
	}
	return s;
}

void iNsidiousProcess(vector<Node>& NODES, vector<int>& temp, unordered_map<string, string> mmuid, vector<transaction_ID> at_time) { // E->N
	for (int i = 0; i < at_time.size(); i++) {
		//如何找到對應?
		int pos_tx = -1, pos_rx = -1;
		//cout << at_time[i].tx_muuid_short << " " << at_time[i].rx_muuid_short << endl;
		if (at_time[i].tx_muuid_short.size() <= 4) {
			pos_tx = distance(mmuid.begin(), mmuid.find(at_time[i].tx_muuid_short));
			pos_rx = distance(mmuid.begin(), mmuid.find(at_time[i].rx_muuid_short));

		}
		else {//ios 不會有mmuid IOS->tx and or IOS ->rx
			string up_s = at_time[i].rx_muuid_short;
			string lo_s = to_lower(up_s);
			for (auto its = mmuid.begin(); its != mmuid.end(); its++) {
				if (its->second == at_time[i].tx_muuid_short) {
					pos_tx = distance(mmuid.begin(), its);
				}
				if (its->second == up_s || its->second == lo_s) {
					pos_rx = distance(mmuid.begin(), its);
				}
			}

			//rx要轉大小寫去找
		}
		//cout << pos_tx << " " << pos_rx << endl;

		//例外處理
		if (pos_tx == -1 || pos_rx == -1) continue;

		
		if (NODES[pos_tx].get_state() == 3 && NODES[pos_rx].get_state() == 0) {
			if (infectedProcess(NODES[pos_rx].iNsidious_rate) && temp[pos_rx] == -1) {
				temp[pos_rx] = 1;
			}
		}
		
	}
}

void DataStore(Physical_network P, Social_network S) {
	//vector 存, total 58 devices
	ifstream  file_COVID19("C:\\Users\\User\\Desktop\\Dataset\\users_tcn_rx.csv", ios::in);

	int num = 0;
	vector<vector<string>> record;
	string line;
	while (getline(file_COVID19, line)){

		//cout << line << endl;

		istringstream delim(line);
		string token;
		vector<string> tmp;
		while (getline(delim, token, ',')) {
			//cout << token << endl;
			tmp.push_back(token);
		}
		record.push_back(tmp);

		//num++;
	} 		
	cout << record.size() << endl;
	//過濾無用的data 與 分割
	vector<transaction_ID> trans_record(record.size());
	for (int i = 1; i < record.size(); i++) {
		//transaction_ID temp;
		trans_record[i].tx_muuid = record[i][1];
		trans_record[i].tx_muuid_short = record[i][10];
		trans_record[i].rx_muuid_short = record[i][9];
		trans_record[i].distan = stod(record[i][5]);
		//轉型 time 處理
		int current = 0; //最初由 0 的位置開始找
		int next;
		//放日期
		//cout << record[i][7] << endl;
		next = record[i][7].find_first_of(" ", current);
		string data = record[i][7].substr(next - current - 2, next - current);
		trans_record[i].data = stoi(data);
		//放時間
		current = next + 1;
		next = record[i][7].find_first_of(":", current);
		data = record[i][7].substr(current, next - current);
		trans_record[i].upload_time_hour = stoi(data);

		current = next + 1;
		next = record[i][7].find_first_of(":", current);
		data = record[i][7].substr(current, next - current);
		trans_record[i].upload_time_minute = stoi(data);

		current = next + 1;
		data = record[i][7].substr(current, record[i][7].size() - current);
		trans_record[i].upload_time_second = round(stod(data));

		//cout << trans_record[i - 1].tx_muuid << " " << trans_record[i - 1].tx_muuid_short << " " << trans_record[i - 1].rx_muuid_short << endl;
		//cout << trans_record[i - 1].data << " " << trans_record[i - 1].upload_time_hour << " " << trans_record[i - 1].upload_time_minute << " " << trans_record[i].upload_time_second << endl;
	}
	cout << trans_record.size() << endl;
	//技術mmuid對應map
	unordered_map<string, string> mmuid;
	for (int i = 1; i < record.size(); i++) {
		string s_1 = record[i][1];
		string s_2 = mmuid_to_short(s_1);
		if (mmuid.find(s_2) == mmuid.end()) {
			mmuid[s_2] = s_1;
		}
	}

	record.clear();

	//cout << mmuid.size() << endl;
	//直接從使用者資料建對應表 補足原本資料的不足
	ifstream  file_user("C:\\Users\\User\\Desktop\\Dataset\\users_user.csv", ios::in);
	vector<vector<string>> record_user;
	while (getline(file_user, line)) {

		//cout << line << endl;

		istringstream delim(line);
		string token;
		vector<string> tmp;
		while (getline(delim, token, ',')) {
			//cout << token << endl;
			tmp.push_back(token);
		}
		record_user.push_back(tmp);
	}
	cout << record_user.size() << endl;
	for (int i = 1; i < record_user.size(); i++) {
		string s_1 = record_user[i][16];
		//cout << s_1 << endl;
		if (s_1.size() <= 30) {
			continue;
		}
		string s_2 = mmuid_to_short(s_1);
		if (mmuid.find(s_2) == mmuid.end()) {
			mmuid[s_2] = s_1;
		}
	}

	record_user.clear();

	cout << mmuid.size() << endl;
	unordered_map<string, string>::iterator it;
	//for (it = mmuid.begin(); it != mmuid.end(); it++) {
		//cout << it->first << " " << it->second << endl;
	//}
	
	//start for simulation
	//不用跑social network也不用跑mobility model，但可以設定range，點數只有51


	vector<vector<double>> res;//沒有確定時間不能寫死
	
	//先蒐集要用的資料
	//27日
	vector < transaction_ID> record_27;
	//28日
	vector < transaction_ID> record_28;
	//29日
	vector < transaction_ID> record_29;
	//466611
	//cout << trans_record[466611].data;
	//cout << trans_record[466611].data;
	//631971
	//for (int i = 0; i < 631971; i++) {
		//cout << i << " " << trans_record[i].data << endl;
		//if (trans_record[i].data == 27) record_27.push_back(trans_record[i]);
		//if (trans_record[i].data == 28) record_28.push_back(trans_record[i]);
		//if (trans_record[i].data == 29) record_29.push_back(trans_record[i]);
	//}

	//cout << record_27.size() << endl;

	srand((unsigned)time(NULL));
	for (int z = 0; z < times; z++) {
		vector<Node> NODES;
		for (int i = 0; i < mmuid.size(); i++) NODES.push_back(Node(i));

		vector<vector<double>> time_record;

		int t_hour, t_min, t_sec;
		int t = 0;
		//initial infection, only one node infection
		int k = 0;
		while (k < 5) {
			int x = (int)(rand() % NODES.size());
			if (NODES[x].get_state() == 0) {
				NODES[x].change_state(3);
				k++;
			}
		}
		//infection 
		int j = 1;


		int count_28 = 0;
		while (j <= trans_record.size() - 1) {
			if (trans_record[j].data != 29) {
				//27 
				//28 1000000
				//29 530000
				//cout << j << endl;
				//break;
				j++;
				continue;
			}
			vector<double> time_temp;
			//t_hour = record_27[j].upload_time_hour, t_min = record_27[j].upload_time_minute, t_sec = record_27[j].upload_time_second;
			t_hour = trans_record[j].upload_time_hour, t_min = trans_record[j].upload_time_minute, t_sec = trans_record[j].upload_time_second;
			time_temp.push_back(t_hour), time_temp.push_back(t_min), time_temp.push_back(t_sec);
			//cout << record_27[j].upload_time_hour << " " << record_27[j].upload_time_minute << " " << record_27[j].upload_time_second << endl;
			time_record.push_back(time_temp);

			vector<transaction_ID> at_time;
			while (j <= trans_record.size() - 1) {
				//cout << record_27[j].upload_time_hour << " " << record_27[j].upload_time_minute << " " << record_27[j].upload_time_second << endl;
				//if (record_27[j].upload_time_hour == t_hour && record_27[j].upload_time_minute == t_min && record_27[j].upload_time_second == t_sec) {
				if (trans_record[j].upload_time_hour == t_hour && trans_record[j].upload_time_minute == t_min && trans_record[j].upload_time_second == t_sec) {
					at_time.push_back(trans_record[j]);
				}
				else {
					//cout << record_27[j].upload_time_hour << " " << record_27[j].upload_time_minute << " " << record_27[j].upload_time_second << endl;
					break;
				}
				j++;
				count_28++;
			}
			//cout << at_time.size() << endl;

			vector<int> temp(NODES.size(), -1);

			iNsidiousProcess(NODES, temp, mmuid, at_time);
			infectedProcess_1(NODES, temp, NODES.size());
			recoveredProcess(NODES, temp, NODES.size());
			deathProcess(NODES, temp, NODES.size());
			loseImProcess(NODES, temp, NODES.size());
			unchanged(NODES, temp, NODES.size());
			//Leaving_and_Coming(NODES, temp, NODES.size());
			update_state(NODES, temp, NODES.size());
			//record_all_state(NODES, res, NODES.size(), t);

			/*
			for (int i = 0; i < at_time.size(); i++) {
			//如何找到對應?
				int pos_tx = -1, pos_rx = -1;
				//cout << at_time[i].tx_muuid_short << " " << at_time[i].rx_muuid_short << endl;
				if (at_time[i].tx_muuid_short.size() <= 4) {
					pos_tx = distance(mmuid.begin(), mmuid.find(at_time[i].tx_muuid_short));
					pos_rx = distance(mmuid.begin(), mmuid.find(at_time[i].rx_muuid_short));

				}
				else {//ios 不會有mmuid IOS->tx and or IOS ->rx
					string up_s = at_time[i].rx_muuid_short;
					string lo_s = to_lower(up_s);
					for (auto its = mmuid.begin(); its != mmuid.end(); its++) {
						if (its->second == at_time[i].tx_muuid_short) {
							pos_tx = distance(mmuid.begin(), its);
						}
						if (its->second == up_s || its->second == lo_s) {
							pos_rx = distance(mmuid.begin(), its);
						}
					}

					//rx要轉大小寫去找
				}
				//cout << pos_tx << " " << pos_rx << endl;

				//例外處理
				if (pos_tx == -1 || pos_rx == -1) continue;


				if (NODES[pos_tx].get_state() == 1) {
					if (infectedProcess(NODES[pos_tx].infected_rate)) {
						NODES[pos_rx].change_state(1);
					}
				}

			}
			*/

			if(t >= res.size()) res.push_back(vector<double>(6, 0.0));
			for (int i = 0; i < NODES.size(); i++) {
				int state_value = NODES[i].get_state();
				res[t][state_value] += 1 / ((double)NODES.size()*times);
			}
			t++;
		}

		//cout << res.size() << endl;

		//for (int i = 0; i < res.size(); i++) {
			//cout << time_record[i][0] << ":" << time_record[i][1] << ":" << time_record[i][2] << endl;
			//for (int g = 0; g < res[i].size(); g++) {
				//cout << res[i][g] << " ";
			//}
			//cout << endl;
		//}

		//cout << z << endl;
		cout << count_28 << endl;
	}
	//Exp_Record(res, 1, 1, 29);//simulation

	unordered_map<int, vector<vector<double>>> record_ana;
	vector<vector<double>> res3;
	vector<Node_a> NODES_A;
	vector<Node> NODES_no_use;
	for (int i = 0; i < mmuid.size(); i++) NODES_no_use.push_back(Node(i));

	res3 = process_a(NODES_A, S, P, NODES_no_use, 1, record_ana);
	//res3 = process_a(NODES_A, S, P, NODES, k, record);
	Exp_Record(res3, 2, 1, 29); //要修改

}



