/*
 * main_lkf_final.cpp - Simple wrapper for LKF final assembly
 * Team ChaiGPT - Milestone 4
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>

using namespace std;

const int TOTAL_STATES = 276;
const int TOTAL_MEAS = 69;

extern "C" {
    void run_lkf_asm(double* meas, int num_frames, double* F, double* Q, 
                     double* H, double* R, double* all_states);
    void buildF(double* F, double dt);
    void buildQ(double* Q, double dt, double sigma);
    void buildH(double* H);
}

int load_measurements(const string& filename, vector<double>& meas, int& num_frames) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "ERROR: Cannot open " << filename << endl;
        return -1;
    }
    
    string line;
    getline(file, line); // Skip header
    
    vector<vector<double>> data;
    while (getline(file, line)) {
        if (line.empty()) continue;
        stringstream ss(line);
        vector<double> row;
        string value;
        while (getline(ss, value, ',')) {
            row.push_back(stod(value));
        }
        if (row.size() == TOTAL_MEAS) {
            data.push_back(row);
        }
    }
    file.close();
    
    num_frames = data.size();
    meas.resize(num_frames * TOTAL_MEAS);
    for (int k = 0; k < num_frames; k++) {
        for (int i = 0; i < TOTAL_MEAS; i++) {
            meas[k * TOTAL_MEAS + i] = data[k][i];
        }
    }
    return 0;
}

int save_states(const string& filename, const vector<double>& states, int num_frames) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "ERROR: Cannot write " << filename << endl;
        return -1;
    }
    
    file << fixed << setprecision(8);
    for (int k = 0; k < num_frames; k++) {
        for (int i = 0; i < TOTAL_STATES; i++) {
            file << states[k * TOTAL_STATES + i];
            if (i < TOTAL_STATES - 1) file << ",";
        }
        file << "\n";
    }
    file.close();
    return 0;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <input.csv> <output.csv>" << endl;
        return 1;
    }
    
    cout << "LINEAR KALMAN FILTER - Milestone 4" << endl;
    
    // Load measurements
    vector<double> measurements;
    int num_frames;
    if (load_measurements(argv[1], measurements, num_frames) != 0) {
        return 1;
    }
    cout << "Loaded " << num_frames << " frames" << endl;
    
    // Build system matrices
    double dt = 0.01;
    double sigma_jerk = 0.1;
    double sigma_meas = 0.05;
    
    double* F = new double[TOTAL_STATES * TOTAL_STATES];
    double* Q = new double[TOTAL_STATES * TOTAL_STATES];
    double* H = new double[TOTAL_MEAS * TOTAL_STATES];
    double* R = new double[TOTAL_MEAS * TOTAL_MEAS];
    
    buildF(F, dt);
    buildQ(Q, dt, sigma_jerk);
    buildH(H);
    
    for (int i = 0; i < TOTAL_MEAS * TOTAL_MEAS; i++) R[i] = 0.0;
    double r_diag = sigma_meas * sigma_meas;
    for (int i = 0; i < TOTAL_MEAS; i++) {
        R[i * TOTAL_MEAS + i] = r_diag;
    }
    
    cout << "Running LKF..." << endl;
    
    // Run filter
    vector<double> all_states(num_frames * TOTAL_STATES);
    run_lkf_asm(measurements.data(), num_frames, F, Q, H, R, all_states.data());
    
    // Save results
    if (save_states(argv[2], all_states, num_frames) != 0) {
        return 1;
    }
    
    cout << "Saved " << num_frames << " frames to " << argv[2] << endl;
    cout << "SUCCESS" << endl;
    
    delete[] F;
    delete[] Q;
    delete[] H;
    delete[] R;
    
    return 0;
}