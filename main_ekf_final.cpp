/*
 * main_ekf_final.cpp - Simple wrapper for EKF final assembly
 * Team ChaiGPT - Milestone 3
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>

using namespace std;

const int NUM_JOINTS = 23;
const int TOTAL_STATES = 276;
const int TOTAL_MEAS = 69;

extern "C" {
    void run_ekf_asm(double* meas_cart, double* meas_sph, int num_frames,
                     double* F, double* Q, double* R, double* all_states);
    void buildF(double* F, double dt);
    void buildQ(double* Q, double dt, double sigma);
    double arctan2_manual(double y, double x);
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

void cartesian_to_spherical(vector<double>& meas_sph, const vector<double>& meas_cart, int num_frames) {
    meas_sph.resize(num_frames * TOTAL_MEAS);
    
    for (int k = 0; k < num_frames; k++) {
        for (int j = 0; j < NUM_JOINTS; j++) {
            int idx = k * TOTAL_MEAS + j * 3;
            double px = meas_cart[idx + 0];
            double py = meas_cart[idx + 1];
            double pz = meas_cart[idx + 2];
            
            double r = sqrt(px*px + py*py + pz*pz);
            double theta = arctan2_manual(py, px);
            double rho_xy = sqrt(px*px + py*py);
            double phi = arctan2_manual(pz, rho_xy);
            
            meas_sph[idx + 0] = r;
            meas_sph[idx + 1] = theta;
            meas_sph[idx + 2] = phi;
        }
    }
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
    
    cout << "EXTENDED KALMAN FILTER - Milestone 3" << endl;
    
    // Load Cartesian measurements
    vector<double> meas_cart;
    int num_frames;
    if (load_measurements(argv[1], meas_cart, num_frames) != 0) {
        return 1;
    }
    cout << "Loaded " << num_frames << " frames" << endl;
    
    // Convert to spherical
    cout << "Converting to spherical coordinates..." << endl;
    vector<double> meas_sph;
    cartesian_to_spherical(meas_sph, meas_cart, num_frames);
    
    // Build system matrices
    double dt = 0.01;
    double sigma_jerk = 0.1;
    double sigma_r = 0.05;
    double sigma_ang = 0.005;
    
    double* F = new double[TOTAL_STATES * TOTAL_STATES];
    double* Q = new double[TOTAL_STATES * TOTAL_STATES];
    double* R = new double[TOTAL_MEAS * TOTAL_MEAS];
    
    buildF(F, dt);
    buildQ(Q, dt, sigma_jerk);
    
    for (int i = 0; i < TOTAL_MEAS * TOTAL_MEAS; i++) R[i] = 0.0;
    for (int j = 0; j < NUM_JOINTS; j++) {
        R[(j*3 + 0) * TOTAL_MEAS + (j*3 + 0)] = sigma_r * sigma_r;
        R[(j*3 + 1) * TOTAL_MEAS + (j*3 + 1)] = sigma_ang * sigma_ang;
        R[(j*3 + 2) * TOTAL_MEAS + (j*3 + 2)] = sigma_ang * sigma_ang;
    }
    
    cout << "Running EKF..." << endl;
    
    // Run filter
    vector<double> all_states(num_frames * TOTAL_STATES);
    run_ekf_asm(meas_cart.data(), meas_sph.data(), num_frames, F, Q, R, all_states.data());
    
    // Save results
    if (save_states(argv[2], all_states, num_frames) != 0) {
        return 1;
    }
    
    cout << "Saved " << num_frames << " frames to " << argv[2] << endl;
    cout << "SUCCESS" << endl;
    
    delete[] F;
    delete[] Q;
    delete[] R;
    
    return 0;
}
