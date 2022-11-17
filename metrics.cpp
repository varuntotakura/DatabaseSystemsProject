#include "metrics.h"
#include <bits/stdc++.h>
#include <iostream>
#include <stack>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <vector>
#include <cmath>
#include <cstdlib>

using namespace std;

double cal_rmse(double* tuth, double* repir)  {
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

double calDTW(double* s1, double* s2) {
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
    double ret = dp[-1][-1];
    return ret;
}

double cal_cost(double* truth, double* repair, long lmd_a=5, long lmd_d=5) {
    double* s1 = repair;
    double* s2 = truth;
    int n = sizeof(s1);
    int m = sizeof(s2);
    double** dp;
    long lmd_a = lmd_a * (truth[1] - truth[0]);
    long lmd_d = lmd_d * (truth[1] - truth[0]);
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
