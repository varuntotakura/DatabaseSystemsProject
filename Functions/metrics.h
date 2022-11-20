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

double* time2ts(double*, double);

double* equal_series_generate(double, double, double);

double cal_cost(double*, double*, double, double);

double cal_rmse(double*, double*);

double calAccuracy(double*, double*, double*);

double calDTW(double*, double*);

double metric_res(double*, double*, double*, string);