#include <bits/stdc++.h>
#include <iostream>
#include <stack>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <vector>

using namespace std;

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