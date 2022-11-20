// #include "exact.h"
#include <bits/stdc++.h>
#include <iostream>
#include <stack>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <vector>

using namespace std;

double median(double* arr, int size){
   sort(arr, arr+size);
   if (size % 2 != 0)
      return (double)arr[size/2];
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
        i--;
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


