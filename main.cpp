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
#include <ctime>
// #include "Functions\exact.h"
// #include "Functions\metrics.h"

using namespace std;

double median(double* arr, int size){
   sort(arr, arr+size);
   if (size % 2 != 0) {
      return (double)arr[size/2];
   }
   return (double)(arr[(size-1)/2] + arr[size/2])/2.0;
}

double determine_interval(double* t){
    double eps[100];
    int i;
    for(i = 1; i<sizeof(t); i++){
        eps[i-1] = t[i] - t[i-1];
    }
    return median(eps, sizeof(eps));
}

// bool find(double* arr, double elem){
bool find(stack<double> arr, double elem){
    int n = sizeof(arr)/sizeof(arr);
    int i = n;
    bool temp = false;
    while (i >= 0)
    {
        if (arr.top() == elem) {
            arr.pop();
            temp = true;
            break;
        }
        i-=1;
    }
    return temp;
}

bool check_interval_lb(double interval, double min_cost, double* eps_list){
    double c = 0;
    for(int i = 0; i<=sizeof(eps_list); i++){
        c += abs(interval - eps_list[i]);
    }
    return (c <= min_cost);
}

bool check_st_lb(int d, double* eps_list, double min_cost, double lmd_d, double eps_t){
    double c = d * lmd_d;
    for(int i=d; i<=sizeof(eps_list); i++){
        c += abs(eps_t - eps_list[i]);
    }
    return (c < min_cost);
}

double** trace_back(double** op, double* t, double s_0, double eps_t, double m_best){
    int n = sizeof(t);
    // vector<int> M;
    double** M;
    int i = n;
    int j = m_best;
    int count = 0;
    while(i > 0 && j > 0){
        if(op[i][j] == 0){
            M[count][0] = i-1;
            M[count][1] = j-1;
            i = i - 1;
            j = j - 1;
        } else if(op[i][j] == 1){
            M[count][0] = -1;
            M[count][1] = j-1;
            j = j - 1;
        } else{
            M[count][0] = i-1;
            M[count][1] = -1;
            i = i - 1;
        }  
        count += 1;
    }
    // M = reverse(M.begin(), M.end());
    return M;
}

double match_searching(double* t, double eps_t, double s_0, double lmd_a, double lmd_d){
    // vector<double> dp; // double dp[99]; // stack<double> dp;
    // vector<double> op; // double op[99]; // stack<double> op;
    // for(int i=1; i<(n+1); i++){
    //     dp.push_back({});
    //     op.push_back({});
    // }
    int n = sizeof(t);
    double** dp;
    double** op;
    dp[0][0] = 0;
    for(int i=1; i<n+1; i++){
        dp[i][0] = i*lmd_a;
        op[i][0] = 1;
    }
    for(int i=1; i<n+1; i++){
        dp[0][i] = i*lmd_d;
        op[0][i] = 2;
    }
    float m_best = 10e8;
    float m_ub = 10e8;
    float min_cost = 10e8;
    int m = 1;
    while(m <= m_ub){
        for(int i = 1; i <= n; i++){
            double s_m = s_0 + (m-1) * eps_t;
            double move_res = dp[i-1][m-1] + abs(t[i-1]-s_m);
            double add_res = dp[i][m - 1] + lmd_a;
            double del_res = dp[i-1][m]+lmd_d;
            if(move_res <= add_res && move_res <= del_res){
                dp[i][m] = move_res;
                op[i][m] = 0;
            } else if(add_res <= move_res && add_res <= del_res){
                dp[i][m] = add_res;
                op[i][m] = 1;
            } else {
                dp[i][m] = del_res;
                op[i][m] = 2;
            }
        }
        if(dp[n][m] < min_cost){
            double min_cost = dp[n][m];
            int m_best = m;
            double m_ub = floor(min_cost/lmd_a)+n;
        }
        m += 1;
    };
    double** M = trace_back(op, t, s_0, eps_t, m_best);
    return M, min_cost, m_best;
}

double round_to_granularity(double value, double granularity){
    return round(value / granularity) * granularity;
}

double exact_repair(double* t, double lmd_a, double lmd_d, int interval_granularity, int start_point_granularity, int bias_d, int bias_s){
    double* eps_list;
    int n = sizeof(t);
    for(int i = 0; i <= n; i++){
        eps_list[i] = t[i]-t[i-1];
    }
    double eps_md = median(eps_list, n);
    double eps_t = round_to_granularity(eps_md, interval_granularity);
    stack<double> eps_t_traverse_range;// double* eps_t_traverse_range;
    stack<double> eps_t_traversed;// double* eps_t_traversed;
    float min_cost = 10e8;
    double** M;
    float cost;
    float m;
    double min_eps_t;
    double min_s_0;
    bool flag_increase = false;
    bool flag_decrease = false;
    bool min_cost_change = false;
    float m_best;
    while(true){
        int d = 0;
        while(((d == 0) || check_st_lb(d, eps_list, min_cost, lmd_d, eps_t)) && d < n && d < bias_d){
            double s_0 = t[d];
            while(s_0 <= t[d] + bias_s){
                M, cost, m = match_searching(t, eps_t, s_0, lmd_a, lmd_d);
                if(cost < min_cost){
                    min_cost, m_best, min_eps_t, min_s_0 = cost, m, eps_t, s_0;
                    min_cost_change = true;
                    s_0 += start_point_granularity;
                } else {
                    flag_increase = true;
                    break;
                }
            }
            s_0 = t[d]-1;
            while(s_0 >= t[d] - bias_s){
                s_0 -= start_point_granularity;
                M, cost, m = match_searching(t, eps_t, s_0, lmd_a, lmd_d);
                if(cost < min_cost){
                    min_cost, m_best, min_eps_t, min_s_0 = cost, m, eps_t, s_0;
                    min_cost_change = true;
                } else {
                    flag_decrease = true;
                    break;
                }
            }
            if(flag_increase && flag_decrease){
                break;
            }
            d += 1;
        }
        if(!min_cost_change || (!check_interval_lb(eps_t, min_cost, eps_list))){
            break;
        }
        if(find(eps_t_traversed, (eps_t + interval_granularity)) && (eps_t + interval_granularity) <= round_to_granularity(eps_md, interval_granularity) + interval_granularity){
            double temp = eps_t + interval_granularity;
            eps_t_traverse_range.push(temp);
        }
        if(find(eps_t_traversed, (eps_t - interval_granularity)) && (eps_t - interval_granularity) >= round_to_granularity(eps_md, interval_granularity) + interval_granularity){
            double temp = eps_t - interval_granularity;
            eps_t_traverse_range.push(temp);
        }
        eps_t_traversed.push(eps_t);
        if(sizeof(eps_t_traverse_range) == 0){
            break;
        }
        eps_t = eps_t_traverse_range.top();
        eps_t_traverse_range.pop();
    }
    return min_eps_t, min_s_0, m_best;
}

double determine_interval(double* t, double interval_granularity){
    double* eps_list;
    double eps_md;
    for(int i=0; i<=sizeof(t); i++){
        eps_list[i] = t[i] - t[i-1];
    }
    if (sizeof(eps_list) % 2 == 0)
        eps_md = (eps_list[sizeof(eps_list)/2 - 1] + eps_list[sizeof(eps_list)/2]) / 2;
    else
        eps_md = eps_list[sizeof(eps_list)/2];
    return round(eps_md / interval_granularity) * interval_granularity;
}

double start_point_approximation(double* t, double lmd_a=5, double lmd_d=5, double interval_granularity=1000){
    int n = sizeof(t);
    double eps_t = determine_interval(t, interval_granularity);
    double s_0 = t[0];
    double** dp;
    double** op;
    dp[0][0] = 0;
    op[0][0] = 0;
    int count = 0;
    for(int i=0; i<=n; i++){
        dp[i][count] = i * lmd_d;
        op[i][count] = 2;
        count+=1;
    } 
    float m_best = 10e8;
    float m_ub = 10e8;
    float min_cost = 10e8;
    int m = 1;
    double s_m, move_res, add_res, del_res;
    while(m <= m_ub){
        dp[0][count] = m * lmd_a;
        op[0][count] = 1;
        for(int i=0; i<=n; i++){
            s_m = s_0 + (m - 1) * eps_t;
            move_res = dp[i - 1][m - 1] + abs(t[i - 1] - s_m);
            add_res = dp[i][m - 1] + lmd_a;
            del_res = dp[i - 1][m] + lmd_d;
            int c = count;
            if(move_res <= add_res and move_res <= del_res){
                dp[i][c] = move_res;
                op[i][c] = 0;
            }
            else if(add_res <= move_res and add_res <= del_res){
                dp[i][c] = add_res;
                op[i][c] = 1;
            }
            else{
                dp[i][c] = del_res;
                op[i][c] = 2;
            }
            c+=1;
        }
        if(dp[n][m] < min_cost){
            min_cost = dp[n][m];
            m_best = m;
            m_ub = floor(min_cost / lmd_a) + n;
        }
        m += 1;
        count += 1;
    }
    return min_cost, eps_t, s_0, m_best;
}

double median_approximation(double* t, double lmd_a=5, double lmd_d=5, double interval_granularity=1){
    int n = sizeof(t);
    double eps_t = determine_interval(t, interval_granularity);
    int n_md = floor(n/2);
    double s_md ;
    if (sizeof(t) % 2 == 0)
        s_md = (t[sizeof(t)/2 - 1] + t[sizeof(t)/2]) / 2;
    else
        s_md = t[sizeof(t)/2];
    double** dp_l;
    double** dp_r;
    double** op_l;
    double** op_r;
    int count = 0;
    dp_l[0][0] = 0;
    op_l[0][0] = 0;
    dp_r[0][0] = 0;
    op_r[0][0] = 0;
    for(int i=0; i<=n_md; i++){
        dp_l[i][count] = i* lmd_d;
        op_l[i][count] = 2;
        dp_r[i][count] = i * lmd_d;
        op_r[i][count] = 2;
        count += 1;
    }
    float m_best = 10e8;
    float m_ub = 10e8;
    float min_cost = 10e8;
    int m = 1;
    double s_m_l, s_m_r, t_i_l, t_i_r;
    while(m <= m_ub){
        dp_l[0][count] = m*lmd_a;
        op_l[0][count] = 1;
        dp_r[0][count] = m * lmd_a;
        op_r[0][count] = 1;
        for(int i=0; i<=n_md; i++){
            if(n % 2 == 1){
                s_m_l = s_md - m * eps_t;
                s_m_r = s_md + m * eps_t;
                t_i_l = t[(int)((n-1)/2)-i];
                t_i_r = t[(int)((n+1)/2)+(i-1)];
            } else {
                s_m_l = s_md - (m - 0.5)*eps_t;
                s_m_r = s_md + (m - 0.5)*eps_t;
                t_i_l = t[(int)(n / 2)-i];
                t_i_r = t[(int)(n / 2)+i-1];
            }
            double move_res_l = dp_l[i - 1][m - 1] + abs(t_i_l - s_m_l);
            double move_res_r = dp_r[i - 1][m - 1] + abs(t_i_r - s_m_r);
            double add_res_l = dp_l[i][m - 1] + lmd_a;
            double add_res_r = dp_r[i][m - 1] + lmd_a;
            double del_res_l = dp_l[i - 1][m] + lmd_d;
            double del_res_r = dp_r[i - 1][m] + lmd_d;
            double min_res_l = min(move_res_l, add_res_l, del_res_l);
            if(move_res_l == min_res_l){
                dp_l[i][count] = move_res_l;
                op_l[i][count] = 0;
            }
            else if(add_res_l == min_res_l){
                dp_l[i][count] = add_res_l;
                op_l[i][count] = 1;
            } else {
                dp_l[i][count] = del_res_l;
                op_l[i][count] = 2;
            }
            double min_res_r = min(move_res_r, add_res_r, del_res_r);
            if(move_res_r == min_res_r){
                dp_r[i][count] = move_res_r;
                op_r[i][count] = 0;
            }
            else if(add_res_r == min_res_r){
                dp_r[i][count] = add_res_r;
                op_r[i][count] = 1;
            }
            else {
                dp_r[i][count] = del_res_r;
                op_r[i][count] = 2;
            }
        }
        if(dp_r[n_md][m] + dp_l[n_md][m] < min_cost){
            min_cost = dp_r[n_md][m] + dp_l[n_md][m];
            int m_best = m;
            if(n % 2 == 1)
                m_ub = ((int)(floor(min_cost/lmd_a + n)) - 1) / 2;
            else
                m_ub = ((int)(floor(min_cost / lmd_a + n))) / 2;
        }
        m += 1;
    }
    double s_0;
    if(n % 2 == 1){
        s_0 = s_md - m_best * eps_t;
        m = m_best * 2 + 1;
    }
    else{
        s_0 = s_md - (m_best - 0.5) * eps_t;
        m = m_best * 2;
    }
    return min_cost, eps_t, s_0, m;
}

double** trace_back(double** op, double* t, double s_0, double eps_t, double m_best){
    int n = sizeof(t);
    // vector<int> M;
    double** M;
    int i = n - 1;
    int j = m_best - 1;
    int count = 0;
    while(i >= 0 && j >= 0){
        if(op[i][j] == 0){
            M[count][0] = i;
            M[count][1] = j;
            i = i - 1;
            j = j - 1;
        } else if(op[i][j] == 1){
            M[count][0] = -1;
            M[count][1] = j;
            j = j - 1;
        } else{
            M[count][0] = i;
            M[count][1] = -1;
            i = i - 1;
        }  
        count += 1;
    }
    // M = reverse(M.begin(), M.end());
    return M;
}

double median_approximation_all(double* t, double  lmd_a=5, double lmd_d=5, double interval_granularity=1){
    double median_min_cost, median_eps_t, median_s_0, median_m = median_approximation(t, lmd_a, lmd_d, interval_granularity);
    double sp_min_cost, sp_eps_t, sp_median_s_0, sp_m = start_point_approximation(t, lmd_a, lmd_d, interval_granularity);
    if(median_min_cost <= sp_min_cost)
        return median_eps_t, median_s_0, median_m;
    else
        return sp_eps_t, sp_median_s_0, sp_m;
}

double* time2ts(double* seq, double time_scale){
    double* ts_list;
    for(int i = 0; i <= sizeof(seq); i++){
        time_t timeArray = datetime.strptime(seq[i], "%Y-%m-%d %H:%M:%S.%f");
        time_t timeStamp = float(timeArray.timestamp()) * time_scale;
        ts_list[i] = timeStamp;
    }
    return ts_list;
}

double* equal_series_generate(double eps_t, double s_0, double m) {
    double* ret;
    for(int i = 0; i <= m; i++){
        ret[i] = s_0 + i*eps_t;
    }
    return ret;
}

double cal_cost(double* truth, double* repair, double lmd_a, double lmd_d) {
    lmd_a = 5;
    lmd_d = 5;
    double* s1 = repair;
    double* s2 = truth;
    int n = sizeof(s1);
    int m = sizeof(s2);
    double** dp;
    lmd_a = lmd_a * (truth[1] - truth[0]);
    lmd_d = lmd_d * (truth[1] - truth[0]);
    for(int i=1; i<=n+1; i++){
        dp[i][0] = i*lmd_d;
    }
    for(int j=1; j<=m+1; j++){
        dp[0][0] = j*lmd_a;
        for(int i=1; i<=n+1; i++){
            double s_m = s2[j-1];
            double move_res = dp[i-1][j-1] + abs(s1[i-1]-s_m);
            double add_res = dp[i][j-1] + lmd_a;
            double del_res = dp[i-1][j] + lmd_d;
            if(move_res <= add_res and move_res <= del_res) {
                dp[i][0] = move_res;
            }
            else if(add_res <= move_res and add_res <= del_res) {
                dp[i][0] = add_res;
            }
            else {
                dp[i][0] = del_res;
            }
        }
    }
    double res = dp[n][m];
    return res;
}

double cal_rmse(double* tuth, double* repir) {
    int min_len = min(sizeof(tuth), sizeof(repir));
    double* truth;
    double* repair;
    double* diff;
    double res;
    double sum = 0;
    for(int i=0; i<=min_len; i++){
        truth[i] = tuth[i];
        repair[i] = repir[i];
    }
    for(int i=0; i<=min_len; i++){
        diff[i] = pow(abs(truth[i] - repair[i]), 2);
        sum += diff[i];
    }
    res = sqrt(sum / sizeof(diff));
    return res;
}

double calAccuracy(double* tuth, double* falt, double* repir) {
    int min_len = min(sizeof(tuth), sizeof(falt), sizeof(repir));
    double* truth;
    double* repair;
    double* fault;
    double* diff;
    double error = 0;
    double cost = 0;
    double inject;
    for(int i=0; i<=min_len; i++){
        truth[i] = tuth[i];
        repair[i] = repir[i];
        fault[i] = falt[i];
    }
    for(int i=0; i<=min_len; i++){
        diff[i] = abs(truth[i] - repair[i]);
        error += diff[i];
    }
    for(int i=0; i<=min_len; i++){
        diff[i] = abs(fault[i] - repair[i]);
        cost += diff[i];
    }
    for(int i=0; i<=min_len; i++){
        diff[i] = abs(truth[i] - fault[i]);
        inject += diff[i];
    }
    if(error == 0){
        return 1;
    }
    double ret = (1 - (error / (cost + inject)));
    return ret;
}

double calDTW(double* s1, double* s2){
    int m = sizeof(s1);
    int n = sizeof(s2);
    double** dp;
    for(int i=0; i<=m; i++){
        dp[i][0] = abs(s1[i] - s2[0]);
    }
    for(int j=0; j<=n; j++){
        dp[0][j] = abs(s1[0] - s2[j]);
    }
    for(int j=1; j<=n; j++){
        for(int i=1; i<=m; i++){
            dp[i][j] = min(dp[i - 1][j - 1], dp[i - 1][j], dp[i][j - 1]) + abs(s1[i] - s2[j]);
        }
    }
    double ret = dp[-1][-1];
    return ret;
}

double metric_res(double* repair, double* truth, double* fault, string metric_name){
    string metric_name="cost";
    if(metric_name == "cost"){
        double lmd_a = 5 * (truth[1] - truth[0]);
        double lmd_d = 5 * (truth[1] - truth[0]);
        return cal_cost(truth, repair, lmd_a, lmd_d);
    } else if(metric_name == "dtw") {
        return calDTW(truth, repair);
    } else if(metric_name == "accuracy") {
        return calAccuracy(truth, fault, repair);
    } else {
        return cal_rmse(truth, repair);
    }
}

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

        if (data_characteristic){
            eps_t_e, s_0_e, m_e = exact_repair_v(original_seq, source_values, lmd_a, lmd_d, interval_granularity, start_point_granularity);
        } else {
        eps_t_e, s_0_e, m_e = exact_repair(original_seq, lmd_a, lmd_d, interval_granularity, start_point_granularity, bias_d, bias_s);
        }
        
        eps_t_a, s_0_a, m_a = median_approximation_all(original_seq, lmd_a, lmd_d, interval_granularity);

        double* exact_res = equal_series_generate(eps_t_e, s_0_e, m_e);
        double* appro_res = equal_series_generate(eps_t_a, s_0_a, m_a);

        result_rmse = metric_res(exact_res, ground_truth_seq, original_seq, metric);
        // result_map["exact-time"] = exact_time;
        // result_map["approximate-rmse"] = metric_res(appro_res, truth, original, metric);
        // result_map["approximate-time"] = appro_time;
    }
    cout << result_rmse <<endl;      
    // cout << result_time <<endl;   
    return 0;
}