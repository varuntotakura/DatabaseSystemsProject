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
#include <time.h>
#include <array> 
// #include "Functions\exact.h"
// #include "Functions\metrics.h"

using namespace std;

long median(vector<long> arr, long size){
    // cout << "median" <<endl;
    if (size % 2 != 0) {
        return (long)arr[size/2];
    }
    return (long)(arr[(size-1)/2] + arr[size/2])/2.0;
}

long determine_interval(vector<long> t){
    // cout << "determine_interval" <<endl;
    vector<long> eps;
    int i;
    for(i = 1; i<t.size(); i++){
        eps.push_back((long)t[i] - t[i-1]);
    }
    return median(eps, eps.size());
}

// bool find(long* arr, long elem){
bool find(vector<long> arr, long elem){
    // cout << "find" <<endl;
    int n = arr.size();
    int i = n;
    bool temp = false;
    while (i > 0)
    {
        if (arr.back() == elem) {
            arr.pop_back();
            temp = true;
            break;
        }
        i-=1;
    }
    return temp;
}

bool check_interval_lb(long interval, long min_cost, vector<long> eps_list){
    // cout << "check_interval_lb" <<endl;
    long c = 0;
    for(int i = 0; i<=eps_list.size(); i++){
        c += abs(interval - eps_list[i]);
    }
    return (c <= min_cost);
}

bool check_st_lb(int d, vector<long> eps_list, long min_cost, long lmd_d, long eps_t){
    // cout << "check_st_lb" <<endl;
    long c = d * lmd_d;
    for(int i=d; i<=eps_list.size(); i++){
        c += abs(eps_t - eps_list[i]);
    }
    return (c < min_cost);
}

vector<vector<long>> trace_back(vector<vector<long>> op, vector<long> t, long s_0, long eps_t, long m_best){
    // cout << "trace_back" <<endl; 
    int n = t.size();
    // vector<int> M;
    vector<vector<long>> M;
    long i = n;
    long j = m_best;
    int count = 0;
    while(i > 0 && j > 0){
        if(op[i][j] == 0){
            M.push_back({i-1, j-1});
            i = i - 1;
            j = j - 1;
        } else if(op[i][j] == 1){
            M.push_back({-1, j-1});
            j = j - 1;
        } else{
            M.push_back({i-1, -1});
            i = i - 1;
        }  
        count += 1;
    }
    // M = reverse(M.begin(), M.end());
    return M;
}

long match_searching(vector<long> t, long eps_t, long s_0, long lmd_a, long lmd_d){
    // cout << "match_searching" <<endl; 
    // vector<long> dp; // long dp[99]; // stack<long> dp;
    // vector<long> op; // long op[99]; // stack<long> op;
    // for(int i=1; i<(n+1); i++){
    //     dp.push_back({});
    //     op.push_back({});
    // }
    int n = t.size();
    vector<vector<long>> dp;
    vector<vector<long>> op;
    for(int i=0; i<=n; i++){
        dp.push_back(vector<long>());
        op.push_back(vector<long>());
    }
    dp[0].push_back(0);
    op[0].push_back(0);
    for(int i=1; i<n+1; i++){
        dp[i].push_back(i*lmd_a);
        op[i].push_back(1);
    }
    for(int i=1; i<n+1; i++){
        dp[0].push_back(i*lmd_d);
        op[0].push_back(2);
    }
    float m_best = 10e8;
    float m_ub = 10e8;
    float min_cost = 10e8;
    int m = 1;
    while(m <= m_ub){
        for(int i = 1; i <= n; i++){
            long s_m = s_0 + (m-1) * eps_t;
            long move_res = dp[i-1][m-1] + abs(t[i-1]-s_m);
            long add_res = dp[i][m - 1] + lmd_a;
            long del_res = dp[i-1][m]+lmd_d;
            if(move_res <= add_res && move_res <= del_res){
                dp[i].push_back(move_res);
                op[i].push_back(0);
            } else if(add_res <= move_res && add_res <= del_res){
                dp[i].push_back(add_res);
                op[i].push_back(1);
            } else {
                dp[i].push_back(del_res);
                op[i].push_back(2);
            }
        }
        if(dp[n][m] < min_cost){
            min_cost = dp[n][m];
            m_best = m;
            m_ub = floor(min_cost/lmd_a)+n;
        }
        m += 1;
    };
    vector<vector<long>> M = trace_back(op, t, s_0, eps_t, m_best);
    return M, min_cost, m_best;
}

long round_to_granularity(long value, long granularity){
    return round(value / granularity) * granularity;
}

long exact_repair(vector<long> t, long lmd_a, long lmd_d, int interval_granularity, int start_point_granularity, int bias_d, int bias_s){
    // cout << "exact_repair" <<endl; 
    vector<long> eps_list;
    int n = t.size();
    for(int i = 0; i <= n; i++){
        eps_list.push_back((long)t[i]-t[i-1]);
    }
    long eps_md = median(eps_list, n);
    long eps_t = round_to_granularity(eps_md, interval_granularity);
    vector<long> eps_t_traverse_range;// long* eps_t_traverse_range;
    vector<long> eps_t_traversed;// long* eps_t_traversed;
    float min_cost = 10e8;
    vector<vector<long>> M;
    float cost;
    float m;
    long min_eps_t;
    long min_s_0;
    bool flag_increase = false;
    bool flag_decrease = false;
    bool min_cost_change = false;
    float m_best;
    int d = 0;
    while(true){
        while(((d == 0) || check_st_lb(d, eps_list, min_cost, lmd_d, eps_t)) && d < n && d < bias_d){
            long s_0 = t[d];
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
            long temp = eps_t + interval_granularity;
            eps_t_traverse_range.push_back(temp);
        }
        if(find(eps_t_traversed, (eps_t - interval_granularity)) && (eps_t - interval_granularity) >= round_to_granularity(eps_md, interval_granularity) + interval_granularity){
            long temp = eps_t - interval_granularity;
            eps_t_traverse_range.push_back(temp);
        }
        eps_t_traversed.push_back(eps_t);
        if(eps_t_traverse_range.size() == 0){
            break;
        }
        eps_t = eps_t_traverse_range.back();
        eps_t_traverse_range.pop_back();
    }
    return min_eps_t, min_s_0, m_best;
}

long determine_interval(vector<long> t, long interval_granularity){
    // cout << "detrimine_interval" <<endl; 
    vector<long> eps_list;
    long eps_md;
    for(int i=0; i<=t.size(); i++){
        eps_list.push_back((long)t[i] - t[i-1]);
    }
    int len = eps_list.size();//*(&eps_list + 1) - eps_list;
    if (len % 2 == 0)
        eps_md = (eps_list[len/2 - 1] + eps_list[len/2]) / 2;
    else
        eps_md = eps_list[len/2];
    return round(eps_md / interval_granularity) * interval_granularity;
}

long start_point_approximation(vector<long> t, long lmd_a=5, long lmd_d=5, long interval_granularity=1000){
    // cout << "start_point_approximation" <<endl; 
    int n = t.size();
    long eps_t = determine_interval(t, interval_granularity);
    long s_0 = t[0];
    vector<vector<long>> dp;
    vector<vector<long>> op;
    for(int i=0; i<=n; i++){
        dp.push_back(vector<long>());
        op.push_back(vector<long>());
    }
    dp[0].push_back(0);
    op[0].push_back(0);
    for(int i=0; i<=n; i++){
        dp[i].push_back(i * lmd_d);
        op[i].push_back(2);
    } 
    float m_best = 10e8;
    float m_ub = 10e8;
    float min_cost = 10e8;
    int m = 1;
    long s_m, move_res, add_res, del_res;
    while(m <= m_ub){
        dp[0].push_back(m * lmd_a);
        op[0].push_back(1);
        for(int i=1; i<n; i++){
            s_m = s_0 + (m - 1) * eps_t;
            move_res = dp[i - 1][m - 1] + abs(t[i - 1] - s_m);
            add_res = dp[i][m - 1] + lmd_a;
            del_res = dp[i - 1][m] + lmd_d;
            if(move_res <= add_res and move_res <= del_res){
                dp[i].push_back(move_res);
                op[i].push_back(0);
            }
            else if(add_res <= move_res and add_res <= del_res){
                dp[i].push_back(add_res);
                op[i].push_back(1);
            }
            else{
                dp[i].push_back(del_res);
                op[i].push_back(2);
            }
        }
        if(dp[n][m] < min_cost){
            min_cost = dp[n][m];
            m_best = m;
            m_ub = floor(min_cost / lmd_a) + n;
        }
        m += 1;
    }
    return min_cost, eps_t, s_0, m_best;
}

long median_approximation(vector<long> t, long lmd_a=5, long lmd_d=5, long interval_granularity=1){
    // cout << "median approximation" <<endl; 
    int n = t.size();
    long eps_t = determine_interval(t, interval_granularity);
    int n_md = floor(n/2);
    long s_md ;
    if (n % 2 == 0)
        s_md = (t[n/2 - 1] + t[n/2]) / 2;
    else
        s_md = t[n/2];
    vector<vector<long>> dp_l;
    vector<vector<long>> dp_r;
    vector<vector<long>> op_l;
    vector<vector<long>> op_r;
    for(int i=0; i<=n; i++){
        dp_l.push_back(vector<long>());
        op_l.push_back(vector<long>());
        dp_r.push_back(vector<long>());
        op_r.push_back(vector<long>());
    }
    dp_l[0].push_back(0);
    op_l[0].push_back(0);
    dp_r[0].push_back(0);
    op_r[0].push_back(0);
    for(int i=0; i<=n_md; i++){
        dp_l[i].push_back(i * lmd_d);
        op_l[i].push_back(2);
        dp_r[i].push_back(i * lmd_d);
        op_r[i].push_back(2);
    }
    float m_best = 10e8;
    float m_ub = 10e8;
    float min_cost = 10e8;
    int m = 1;
    long s_m_l, s_m_r, t_i_l, t_i_r;
    while(m <= m_ub){
        dp_l[0].push_back(m*lmd_a);
        op_l[0].push_back(1);
        dp_r[0].push_back(m * lmd_a);
        op_r[0].push_back(1);
        for(int i=1; i<n_md; i++){
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
            long move_res_l = dp_l[i-1][m-1] + abs(t_i_l - s_m_l);
            long move_res_r = dp_r[i-1][m-1] + abs(t_i_r - s_m_r);
            long add_res_l = dp_l[i][m - 1] + lmd_a;
            long add_res_r = dp_r[i][m - 1] + lmd_a;
            long del_res_l = dp_l[i-1][m] + lmd_d;
            long del_res_r = dp_r[i-1][m] + lmd_d;
            long min_res_l = min({move_res_l, add_res_l, del_res_l});
            // long temp = min(move_res_l, add_res_l);
            // long min_res_l = min(temp, del_res_l);
            if(move_res_l == min_res_l){
                dp_l[i].push_back(move_res_l);
                op_l[i].push_back(0);
            }
            else if(add_res_l == min_res_l){
                dp_l[i].push_back(add_res_l);
                op_l[i].push_back(1);
            } else {
                dp_l[i].push_back(del_res_l);
                op_l[i].push_back(2);
            }
            long min_res_r = min({move_res_r, add_res_r, del_res_r});
            if(move_res_r == min_res_r){
                dp_r[i].push_back(move_res_r);
                op_r[i].push_back(0);
            }
            else if(add_res_r == min_res_r){
                dp_r[i].push_back(add_res_r);
                op_r[i].push_back(1);
            }
            else {
                dp_r[i].push_back(del_res_r);
                op_r[i].push_back(2);
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
    long s_0;
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

long median_approximation_all(vector<long> t, long  lmd_a=5, long lmd_d=5, long interval_granularity=1){
    // cout << "median_approximation_all" <<endl; 
    long median_min_cost, median_eps_t, median_s_0, median_m = median_approximation(t, lmd_a, lmd_d, interval_granularity);
    long sp_min_cost, sp_eps_t, sp_median_s_0, sp_m = start_point_approximation(t, lmd_a, lmd_d, interval_granularity);
    if(median_min_cost <= sp_min_cost)
        return median_eps_t, median_s_0, median_m;
    else
        return sp_eps_t, sp_median_s_0, sp_m;
}

vector<long> equal_series_generate(long eps_t, long s_0, long m=100) {
    // cout << "equal_series_generate" <<endl; 
    vector<long> ret;
    for(int i = 0; i < m; i++){
        ret.push_back((long)(s_0 + i*eps_t));
    }
    return ret;
}

long cal_cost(vector<long> truth, vector<long> repair, long lmd_a, long lmd_d) {
    // cout << "cal_cost" <<endl; 
    lmd_a = 5;
    lmd_d = 5;
    vector<long> s1;
    for(int i=0; i<=repair.size(); i++){
        s1.push_back(repair[i]);
    }
    vector<long> s2;
    for(int i=0; i<=truth.size(); i++){
        s2.push_back(truth[i]);
    }
    int n = s1.size();
    int m = s2.size();
    vector<vector<long>> dp;
    for(int i=0; i<=n+m; i++){
        dp.push_back(vector<long>());
    }
    lmd_a = lmd_a * (truth[1] - truth[0]);
    lmd_d = lmd_d * (truth[1] - truth[0]);
    for(int i=1; i<=n+1; i++){
        dp[i].push_back(i*lmd_d);
    }
    for(int j=1; j<=m+1; j++){
        dp[0].push_back(j*lmd_a);
        for(int i=1; i<=n+1; i++){
            long s_m = s2[j-1];
            long move_res = dp[i-1][j-1] + abs(s1[i-1]-s_m);
            long add_res = dp[i][j-1] + lmd_a;
            long del_res = dp[i-1][j] + lmd_d;
            if(move_res <= add_res and move_res <= del_res) {
                dp[i].push_back(move_res);
            }
            else if(add_res <= move_res and add_res <= del_res) {
                dp[i].push_back(add_res);
            }
            else {
                dp[i].push_back(del_res);
            }
        }
    }
    long res = dp[n][m];
    return res;
}

long cal_rmse(vector<long> truth, vector<long> repair) {
    // cout << "cal_rmse" <<endl; 
    int min_len = min(truth.size(), repair.size());
    vector<long> diff;
    long res;
    long sum = 0;
    for(int i=0; i<=min_len; i++){
        diff.push_back((long)pow(abs(truth[i] - repair[i]), 2));
        sum += diff[i];
    }
    res = sqrt(sqrt(sum / diff.size()));
    return res;
}

long calAccuracy(vector<long> truth, vector<long> fault, vector<long> repair) {
    // cout << "calAccuracy" <<endl; 
    int min_len;
    min_len = min({truth.size(), fault.size(), repair.size()});
    vector<long> diff;
    long error = 0;
    long cost = 0;
    long inject = 0;
    for(int i=0; i<=min_len; i++){
        diff.push_back((long)abs(truth[i] - repair[i]));
        error += diff[i];
    }
    for(int i=0; i<=min_len; i++){
        diff.push_back((long)abs(fault[i] - repair[i]));
        cost += diff[i];
    }
    for(int i=0; i<=min_len; i++){
        diff.push_back((long)abs(truth[i] - fault[i]));
        inject += diff[i];
    }
    if(error == 0){
        return 1;
    }
    long ret = (1 - (error / (cost + inject)))*100;
    return ret;
}


long metric_res(vector<long> repair, vector<long> truth, vector<long> fault, string metric_name){
    // cout << "metric_res" <<endl;
    repair = fault;
    if(metric_name == "cost"){
        long lmd_a = 5 * (truth[1] - truth[0]);
        long lmd_d = 5 * (truth[1] - truth[0]);
        return cal_cost(truth, repair, lmd_a, lmd_d);
    } else if(metric_name == "accuracy") {
        return calAccuracy(truth, fault, repair);
    } else if(metric_name == "rmse") {
        return cal_rmse(truth, repair);
        // return cal_rmse(truth, fault);
    }
    return 0;
}

int main(){ 
    string version = "-test";
    string datasets = "energy";
    string methods = "exact"; //aproximate
    vector<vector<string>> content;
    vector<string> row;
    string line, word;
    map<string, long> result_dfs = {};

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
    int bias_d = 1;
    int bias_s = 3;
    long result_exact_rmse = 0;
    long result_approx_rmse = 0;
    for(int ts = 0; ts<file_counts; ts++){
        string file_name = "D:\\FSU\\Semester 1\\Database Systems\\Project\\DatabaseSystemsProject\\data\\dirty_energy_test\\series_0.csv";
        string data_truth = "D:\\FSU\\Semester 1\\Database Systems\\Project\\DatabaseSystemsProject\\data\\energy_test\\series_0.csv";
        vector<long> original_seq;
        vector<long> ground_truth_seq;
        vector<vector<long>> data_rows;

        int cnt = 0;
        long* original_sq;
        ifstream inputFile1;
        inputFile1.open(file_name);
        getline(inputFile1, line);
        line = "";
        while(getline(inputFile1, line)){
            string first_column;
            string second_column;
            stringstream inputString(line);
            getline(inputString, first_column, ',');
            getline(inputString, second_column);
            original_seq.push_back(stol(second_column));
            cnt += 1;
            line = "";
        }
        cnt = 0;
        long* truth_sq;
        ifstream inputFile2;
        inputFile2.open(data_truth);
        getline(inputFile2, line);
        line = "";
        while(getline(inputFile2, line)){
            string first_column;
            string second_column;
            string thrid_column;
            stringstream inputString(line);
            getline(inputString, first_column, ',');
            getline(inputString, second_column, ',');
            getline(inputString, thrid_column);
            ground_truth_seq.push_back(stol(second_column));
            cnt += 1;
            line = "";
        }

        // cout << original_seq[23];

        long source_values = 0;
        long time_scale;
        
        eps_t_e, s_0_e, m_e = exact_repair(original_seq, lmd_a, lmd_d, interval_granularity, start_point_granularity, bias_d, bias_s);
        eps_t_a, s_0_a, m_a = median_approximation_all(original_seq, lmd_a, lmd_d, interval_granularity);

        vector<long> exact_res = equal_series_generate(eps_t_e, s_0_e, m_e);
        // vector<long> appro_res = equal_series_generate(eps_t_a, s_0_a, m_a);

        string metric = "rmse"; //rmse cost accuracy
        result_exact_rmse = metric_res(exact_res, ground_truth_seq, original_seq, metric);
        // result_approx_rmse = metric_res(appro_res, ground_truth_seq, original_seq, metric);
    }
    cout << result_exact_rmse <<endl;      
    // cout << result_approx_rmse <<endl;   
    return 0;
}