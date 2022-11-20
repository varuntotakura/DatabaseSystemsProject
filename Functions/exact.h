// #ifndef EXACT_H    // To make sure you don't declare the function more than once by including the header multiple times.
// #define EXACT_H

// #include "exact.h"
#include <bits/stdc++.h>
#include <iostream>
#include <stack>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <vector>

using namespace std;

double median(double*, int);

double determine_interval(double*);

// bool find(double*, double);
bool find(stack<double>, double);

bool check_interval_lb(double, double, double*);

bool check_st_lb(int, double*, double, double, double);

double** trace_back(double**, double*, double, double, double);

double match_searching(double*, double, double, double, double);

double round_to_granularity(double, double);

double exact_repair(double*, double, double, int, int, int, int);