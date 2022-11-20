#include <bits/stdc++.h>
#include <iostream>
#include <stack>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include "Functions\exact.h"
#include "Functions\metrics.h"

using namespace std;

int main(){
    // double match_search = match_searching();  
    // double** trace = trace_back();  
    // double exact_repr = exact_repair();  
    string version = "-test";
    string datasets = "energy";
    string methods = "exact"; //aproximate
    bool data_characteristic = false;
    vector<vector<string>> content;
    vector<string> row;
    string line, word;
    map<string, double> result_dfs = {};

    result_dfs["rmse"] = {0};
    result_dfs["time"] = {0};

    int file_counts = 5;
    int truth_col = 1;
    string truth_dir = "../data/energy";
    int original_col = 1;
    string original_dir = "../data/dirty_energy";
    int start_point_granularity = 1;
    int interval_granularity = 1;
    int lmd_a = 100;
    int lmd_d = 100;
    double eps_t_e, s_0_e, m_e;
    double eps_t_a, s_0_a, m_a;
    int bias_d = 1;
    int bias_s = 3;
    // map<string, double> result_map = {};
    // result_map["exact-rmse"] = 0;
    // result_map["exact-time"] = 0;
    double result_rmse = 0;
    string dataset_path = "../result/energy";
    for(int ts = 0; ts<file_counts; ts++){
        string file_name = "../data/dirty_energy/series_0.csv";
        string data_truth = "../data/energy/series_0.csv";
        double* original_seq;
        double* file_rows_2;
        double* ground_truth_seq;
        double* data_rows_2;
        string metric = "cost";
        double** data_rows;
        int cnt = 0;

        fstream file(file_name, ios::in);
        fstream data(data_truth, ios::in);

        if (file.is_open())
        {
            while (getline(file, line))
            {
                row.clear();

                stringstream str(line);

                while (getline(str, word, ','))
                    row.push_back(word);
                    cout << row[0] << endl;
                    cout << row[1] << endl;
                    cout << row[2] << endl;
                    original_seq[cnt] = stol(row[1]);
                    file_rows_2[cnt] = stol(row[2]);
                    cnt += 1;
                content.push_back(row);
            }
        }
        else {
            cout << "Could not open the file\n";
        }
        cnt = 0;
        if (data.is_open())
        {
            while (getline(data, line))
            {
                row.clear();

                stringstream str(line);

                while (getline(str, word, ','))
                    row.push_back(word);
                    cout << row[0] << endl;
                    cout << row[1] << endl;
                    cout << row[2] << endl;
                    ground_truth_seq[cnt] = stol(row[1]);
                    file_rows_2[cnt] = stol(row[2]);
                    cnt += 1;
                content.push_back(row);
            }
        }
        else {
            cout << "Could not open the file\n";
        }

        double source_values = 0;
        double time_scale;
        if(time_scale){
            original_seq = time2ts(original_seq, time_scale);
            ground_truth_seq = time2ts(ground_truth_seq, time_scale);
        }

        // if (data_characteristic){
        //     eps_t_e, s_0_e, m_e = exact_repair_v(original_seq, source_values, lmd_a, lmd_d, interval_granularity, start_point_granularity);
        // } else {
        eps_t_e, s_0_e, m_e = exact_repair(original_seq, lmd_a, lmd_d, interval_granularity, start_point_granularity, bias_d, bias_s);
        // }
        
        // eps_t_a, s_0_a, m_a = median_approximation_all(original_seq, lmd_a, lmd_d, interval_granularity);

        double* exact_res = equal_series_generate(eps_t_e, s_0_e, m_e);
        // double* appro_res = equal_series_generate(eps_t_a, s_0_a, m_a);

        result_rmse = metric_res(exact_res, ground_truth_seq, original_seq, metric);
        // result_map["exact-time"] = exact_time;
        // result_map["approximate-rmse"] = metric_res(appro_res, truth, original, metric);
        // result_map["approximate-time"] = appro_time;
    }
    cout << result_rmse <<endl;      
    // cout << result_time <<endl;   
    return 0;
}