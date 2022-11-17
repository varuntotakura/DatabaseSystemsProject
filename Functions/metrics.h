// #include "metrics.h"
#include <bits/stdc++.h>
#include <iostream>
#include <stack>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <vector>
#include <cstdlib>

using namespace std;

double* time2ts(double*, long);

double* equal_series_generate(long, long, long);

double cal_cost(double*, double*, long, long);

double cal_rmse(double*, double*);

double calAccuracy(double*, double*, double*);

double calDTW(double*, double*);

double metric_res(double*, double*, double*, string);