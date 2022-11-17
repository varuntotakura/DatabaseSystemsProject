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

long median(long*, int);

long determine_interval(long*);

// bool find(long*, long);
bool find(stack<long>, long);

bool check_interval_lb(long, long, long*);

bool check_st_lb(int, long*, long, long, long);

long** trace_back(double**, double*, long, long, long);

long match_searching(double*, long, long, long, long);

long round_to_granularity(long, long);

long exact_repair(double*, long, long, int, int, int, int);