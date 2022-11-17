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
#include <exact.h>
#include <metrics.h>

long exact_repair();
long match_searching();
long** trace_back();
double cal_cost();
double calDTW();
double calAccuracy();
double cal_rmse();

using namespace std;

// long* time2ts(long* seq, long time_scale){
//     long* ts_list;
//     for(int i = 0; i <= sizeof(seq); i++){
//         long timeArray = datetime.strptime(seq[i], "%Y-%m-%d %H:%M:%S.%f");
//         long timeStamp = float(timeArray.timestamp()) * time_scale;
//         ts_list[i] = timeStamp;
//     }
//     return ts_list;
// }

double* equal_series_generate(long eps_t, long s_0, long m){
    double* ret;
    for(int i = 0; i <= m; i++){
        ret[i] = s_0 + i*eps_t;
    }
    return ret;
}

double metric_res(double* repair, double* truth, double* fault, string metric_name="cost"){
    double calCost = cal_cost();  
    double cal_DTW = calDTW();  
    double cal_Accuracy = calAccuracy();  
    double calRmse = cal_rmse();
    if(metric_name == "cost"){
        long lmd_a = 5 * (truth[1] - truth[0]);
        long lmd_d = 5 * (truth[1] - truth[0]);
        return calCost(truth, repair, lmd_a, lmd_d);
    } else if(metric_name == "dtw"){
        return cal_DTW(truth, repair);
    } else if(metric_name == "accuracy"){
        return cal_Accuracy(truth, fault, repair);
    } else{
        return calRmse(truth, repair);
    }
}

void main(){
    long match_search = match_searching();  
    long** trace = trace_back();  
    long exact_repr = exact_repair();  
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
    long eps_t_e, s_0_e, m_e;
    long eps_t_a, s_0_a, m_a;
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
        long** data_rows;
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

        long source_values = 0;
        // if "time_scale" in param:
        //     original = time2ts(original_seq, param["time_scale"])
        //     truth = time2ts(ground_truth_seq, param["time_scale"])

        // if (data_characteristic){
        //     eps_t_e, s_0_e, m_e = exact_repair_v(original_seq, source_values, lmd_a, lmd_d, interval_granularity, start_point_granularity);
        // } else {
        eps_t_e, s_0_e, m_e = exact_repair(original_seq, lmd_a, lmd_d, interval_granularity, start_point_granularity);
        // }
        
        // eps_t_a, s_0_a, m_a = median_approximation_all(original_seq, lmd_a, lmd_d, interval_granularity);

        double* exact_res;
        exact_res = equal_series_generate(eps_t_e, s_0_e, m_e);
        // double* appro_res = equal_series_generate(eps_t_a, s_0_a, m_a);

        result_rmse = metric_res(exact_res, ground_truth_seq, original_seq, metric);
        // result_map["approximate-rmse"] = metric_res(appro_res, truth, original, metric);
        // result_map["approximate-time"] = appro_time;

        cout << result_rmse <<endl;

        // for metric in (metrics + ["time"]):
        //     result_dfs["rmse"].at[dataset, "exact"] = np.mean(result_map[f"exact-{metric}"])
        //     np.savetxt(os.path.join(dataset_path, f"exact-{metric}{version}.txt"), result_map[f"exact-{metric}"])
        //     result_dfs[metric].at[dataset, "approximate"] = np.mean(result_map[f"approximate-{metric}"])
        //     np.savetxt(os.path.join(dataset_path, f"approximate-{metric}{version}.txt"), result_map[f"approximate-{metric}"])

    // for metric in (metrics + ["time"]):
    //     result_dfs[metric].to_csv(os.path.join("result", f"exp1-{metric}{version}.csv"))
        
    return;
}