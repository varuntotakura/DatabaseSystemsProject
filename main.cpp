#include <bits/stdc++.h>
#include <iostream>
#include <stack>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>

using namespace std;

void main(){
    string version = "-test";
    string datasets = "energy";
    string methods = "exact"; //aproximate
    bool data_characteristic = false;
    map<string, double> result_dfs = {};

    result_dfs["rmse"] = pd.DataFrame(0, columns=methods, index=datasets);
    result_dfs["time"] = pd.DataFrame(0, columns=methods, index=datasets);

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
    map<string, double> result_map = {};
    result_map["exact-rmse"] = 0;
    result_map["exact-time"] = 0;
    string dataset_path = "../result/energy";
    for(int ts = 0; ts<file_counts; ts++){
        string file_name = "../data/dirty_energy/series_"+ts+".csv";
        string data_truth = "../data/energy/series_"+ts+".csv";

        data = pd.read_csv(file_name);
        original_seq = data.iloc[:, param["original_col"]]
        source_values = None
        # source_values = data.iloc[:, param["original_col"] + 1] # values
        ground_truth_seq = data_truth.iloc[:, param["truth_col"]]

        double* original = original_seq;
        double* truth = ground_truth_seq;

        start = time.time()
        
        if (data_characteristic){
            eps_t_e, s_0_e, m_e = exact_repair_v(original, source_values, lmd_a, lmd_d, interval_granularity, start_point_granularity);
        } else {
            eps_t_e, s_0_e, m_e = exact_repair(original, lmd_a, lmd_d, interval_granularity, start_point_granularity);
        }
        
        end = time.time()
        double exact_time = end - start;

        start = time.time()
        
        eps_t_a, s_0_a, m_a = median_approximation_all(original, lmd_a, lmd_d, interval_granularity);
        
        end = time.time()

        double appro_time = end - start;

        double* exact_res = equal_series_generate(eps_t_e, s_0_e, m_e);
        double* appro_res = equal_series_generate(eps_t_a, s_0_a, m_a);

        result_map["exact-rmse"] = metric_res(exact_res, truth, original, metric);
        // result_map["approximate-rmse"] = metric_res(appro_res, truth, original, metric);
        result_map["exact-time"] = exact_time;
        // result_map["approximate-time"] = appro_time;

        for metric in (metrics + ["time"]):
            result_dfs["rmse"].at[dataset, "exact"] = np.mean(result_map[f"exact-{metric}"])
            np.savetxt(os.path.join(dataset_path, f"exact-{metric}{version}.txt"), result_map[f"exact-{metric}"])
            result_dfs[metric].at[dataset, "approximate"] = np.mean(result_map[f"approximate-{metric}"])
            np.savetxt(os.path.join(dataset_path, f"approximate-{metric}{version}.txt"), result_map[f"approximate-{metric}"])

    for metric in (metrics + ["time"]):
        result_dfs[metric].to_csv(os.path.join("result", f"exp1-{metric}{version}.csv"))
        
    return;
}


long* time2ts(long* seq, long time_scale){
    long* ts_list;
    for(int i = 0; i <= sizeof(seq); i++){
        long timeArray = datetime.strptime(seq[i], "%Y-%m-%d %H:%M:%S.%f");
        long timeStamp = float(timeArray.timestamp()) * time_scale;
        ts_list[i] = timeStamp;
    }
    return ts_list;
}


double* equal_series_generate(long eps_t, long s_0, long m){
    double* ret;
    for(int i = 0; i <= m; i++){
        ret[i] = s_0 + i*eps_t;
    }
    return ret;
}

double metric_res(long* repair, long* truth, fault, string metric_name="cost"){
    if(metric_name == "cost"){
        long lmd_a = 5 * (truth[1] - truth[0]);
        long lmd_d = 5 * (truth[1] - truth[0]);
        return cal_cost(truth, repair, lmd_a, lmd_d);
    } else if(metric_name == "dtw"){
        return calDTW(truth, repair);
    } else if(metric_name == "accuracy"){
        return calAccuracy(truth, fault, repair);
    } else{
        return cal_rmse(truth, repair);
    }
}