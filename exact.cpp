#include <bits/stdc++.h>
#include <iostream>
#include <stack>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <vector>

using namespace std;

double median(double arr[], int size){
   sort(arr, arr+size);
   if (size % 2 != 0)
      return (double)arr[size/2];
   return (double)(arr[(size-1)/2] + arr[size/2])/2.0;
}

double determine_interval(double t[]){
    double eps[100];
    int i;
    for(i = 1; i<sizeof(t); i++){
        eps[i-1] = t[i] - t[i-1];
    }
    return median(eps, sizeof(eps));
}

double match_searching(double t[], double eps_t, double s_0, double lmd_a, double lmd_d){
    int n = sizeof(t);
    double dp[99]; // stack<double> dp;
    double op[99]; // stack<double> op;
    double temp[99];
    for(int i=1; i<(n+1); i++){
        dp[i].push(temp);
        op[i].push(temp);
    }
    dp[0].push(0);
    op[0].push(0);
    for(int i=1; i<n+1; i++){
        dp[i].push(i*lmd_d);
        op[i].push(2);
    }
    double m_best = 10e8;
    double m_ub = 10e8;
    double min_cost = 10e8;
    int m = 1;
    while(m <= m_ub){
        dp[0].push(m*lmd_a);
        op[0].push(1);
        for(int i = 1; i<n+1; i++){
            double s_m = s_0 + (m-1) * eps_t;
            double move_res = dp[i-1][m-1]+abs(t[i-1]-s_m);
            double add_res = dp[i][m - 1] + lmd_a;
            double del_res = dp[i-1][m]+lmd_d;
            if(move_res <= add_res && move_res <= del_res){
                dp[i].push(move_res);
                op[i].push(0);
            } else if(add_res <= move_res && add_res <= del_res){
                dp[i].push(add_res);
                op[i].push(1);
            } else {
                dp[i].push(del_res);
                op[i].push(2);
            }
        }
        if(dp[n][m] < min_cost){
            double min_cost = dp[n][m];
            double m_best = m;
            double m_ub = floor(min_cost/lmd_a)+n;
        }

        m += 1;
    }
    double M = trace_back(op, t, (s_0, eps_t, m_best));
    return M, min_cost, m_best;
}

double trace_back(double op, int t, double s){
    double s_0, eps_t, m_best;
    m_best = s;
    int n = sizeof(t);
    vector<int> M;
    int i = n;
    int j = m_best;
    while(i > 0 && j > 0){
        if(op[i][j] == 0){
            M.push_back([i-1, j-1]);
            i = i - 1;
            j = j - 1;
        } else if(op[i][j] == 1){
            M.push_back([-1, j-1]);
            j = j - 1;
        } else{
            M.push_back([i-1, -1]);
            i = i - 1;
        }  
    }
    M = reverse(M.begin(), M.end());
    double ret[];
    ret = M;
    return ret;
}


double round_to_granularity(double value, double granularity){
    return round(value / granularity) * granularity;
}

double exact_repair(int t, long lmd_a=10, long lmd_d=10, int interval_granularity=1, int start_point_granularity=1, int bias_d=1, int bias_s=3){
    double eps_list[];
    for(int i=0; i<sizeof(t)-1; i++){
        eps_list.push(t[i]-t[i-1]);
    }
    double eps_md = median(eps_list);
    double eps_t = round_to_granularity(eps_md, interval_granularity);
    double eps_t_traverse_range[];
    double eps_t_traversed[];
    float min_cost = 10e8;
    while(1){
        int d = 0;
        while(d == 0 || check_st_lb(d, eps_list, min_cost, lmd_d, eps_t)) && d < len(eps_list) && d < bias_d){
            double s_0 = t[d];
            bool flag_increase = false;
            bool flag_decrease = false;
            bool min_cost_change = false;
            while(s_0 <= t[d] + bias_s){
                double M, cost, m = match_searching(t, eps_t, s_0, lmd_a, lmd_d);
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
        if(not min_cost_change || (not check_interval_lb(eps_t, min_cost, eps_list))){
            break;
        }
        bool a = find(begin(eps_t_traversed), end(eps_t_traversed), (eps_t + interval_granularity)) != end(eps_t_traversed);
        if(a && (eps_t + interval_granularity) <= round_to_granularity(eps_md, interval_granularity) + interval_granularity){
            eps_t_traverse_range.push([eps_t + interval_granularity]);
        }
        bool b = find(begin(eps_t_traversed), end(eps_t_traversed), (eps_t - interval_granularity)) != end(eps_t_traversed);
        if(b && (eps_t - interval_granularity) >= round_to_granularity(eps_md, interval_granularity) + interval_granularity){
            eps_t_traverse_range.push([eps_t - interval_granularity]);
        }
        eps_t_traversed.push(eps_t);
        if(len(eps_t_traverse_range) == 0){
            break;
        }
        eps_t = eps_t_traverse_range.pop();
    }
    return min_eps_t, min_s_0, m_best;
}

bool check_interval_lb(long interval, long min_cost, double eps_list[]){
    double c = 0;
    for(int eps = 0; i<=sizeof(eps_list); i++){
        c += abs(interval - eps_list[eps]);
    }
    return (c <= min_cost);
}

bool check_st_lb(int d, double eps_list[], long min_cost, long lmd_d, double eps_t){
    double c = d * lmd_d;
    for(int i=d; i<=sizeof(eps_list); i++){
        c += abs(eps_t - eps_list[eps]);
    }
    return (c < min_cost);
}












