#include <iostream>

using namespace std;

void main() {
    char parameters = {
        "energy":{
            "file_counts": 5,
            "truth_col": 1,
            "truth_dir": "../data/energy",
            "original_col": 1,
            "original_dir": "../data/dirty_energy",
            "start_point_granularity": 1,
            "interval_granularity": 1,
            "lmd_a": 100,
            "lmd_d": 100,
        },
    };
    char version = "-test";
    char datasets = ["energy"];
    char metrics = ["rmse"];
    char methods = ["exact", "approximate"];
    bool data_characteristic = false; 

    char result_dfs = {};
    return;
}