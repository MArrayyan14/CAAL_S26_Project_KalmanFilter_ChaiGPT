/*
 * Linear Kalman Filter - Milestone 2
 * Team ChaiGPT: Laiba (31736), Ameer (30535), Usman (29126), Arrayyan (30557)
 * 
 * UPDATED: Handles CSV with header row
 * - Input: 3041 rows (1 header + 3040 data)
 * - Output: 3040 rows (no header, data only)
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <Eigen/Dense>
#include <chrono>

using namespace Eigen;
using namespace std;

const int NUM_JOINTS = 23;
const int STATES_PER_JOINT = 12;
const int TOTAL_STATES = 276;
const int TOTAL_MEAS = 69;

// ============================================================================
// CSV LOADING WITH HEADER SKIP
// ============================================================================

vector<VectorXd> loadMeasurements(const string& filename, int& num_frames) {
    vector<VectorXd> measurements;
    ifstream file(filename);
    
    if (!file.is_open()) {
        cerr << "ERROR: Cannot open " << filename << endl;
        exit(1);
    }
    
    string line;
    int line_number = 0;
    
    // Skip header line
    if (getline(file, line)) {
        line_number++;
        cout << "  Skipped header row" << endl;
    }
    
    // Read data rows
    while (getline(file, line)) {
        line_number++;
        if (line.empty()) continue;
        
        stringstream ss(line);
        VectorXd meas(TOTAL_MEAS);
        string value;
        int idx = 0;
        
        while (getline(ss, value, ',') && idx < TOTAL_MEAS) {
            meas(idx++) = stod(value);
        }
        
        if (idx == TOTAL_MEAS) {
            measurements.push_back(meas);
        }
    }
    
    file.close();
    num_frames = measurements.size();
    return measurements;
}

// ============================================================================
// MATRICES
// ============================================================================

MatrixXd buildF(double dt) {
    MatrixXd F = MatrixXd::Zero(TOTAL_STATES, TOTAL_STATES);
    double dt2 = dt * dt;
    double dt3 = dt2 * dt;
    
    Matrix4d F_axis;
    F_axis << 1.0,  dt,  dt2/2.0,  dt3/6.0,
              0.0, 1.0,      dt,   dt2/2.0,
              0.0, 0.0,     1.0,        dt,
              0.0, 0.0,     0.0,       1.0;
    
    for (int joint = 0; joint < NUM_JOINTS; ++joint) {
        int base = joint * STATES_PER_JOINT;
        for (int axis = 0; axis < 3; ++axis) {
            int start = base + axis * 4;
            F.block<4, 4>(start, start) = F_axis;
        }
    }
    return F;
}

MatrixXd buildQ(double dt, double sigma) {
    MatrixXd Q = MatrixXd::Zero(TOTAL_STATES, TOTAL_STATES);
    
    double dt2 = dt * dt;
    double dt3 = dt2 * dt;
    double dt4 = dt3 * dt;
    double dt5 = dt4 * dt;
    double dt6 = dt5 * dt;
    double sigma2 = sigma * sigma;
    
    // DISCRETE-TIME Q matrix (corrected from Milestone 1 feedback)
    Matrix4d Q_axis;
    Q_axis << dt6/36.0,  dt5/12.0,  dt4/6.0,  dt3/6.0,
              dt5/12.0,   dt4/4.0,  dt3/2.0,  dt2/2.0,
               dt4/6.0,   dt3/2.0,      dt2,       dt,
               dt3/6.0,   dt2/2.0,       dt,      1.0;
    Q_axis *= sigma2;
    
    for (int joint = 0; joint < NUM_JOINTS; ++joint) {
        int base = joint * STATES_PER_JOINT;
        for (int axis = 0; axis < 3; ++axis) {
            int start = base + axis * 4;
            Q.block<4, 4>(start, start) = Q_axis;
        }
    }
    return Q;
}

MatrixXd buildH() {
    MatrixXd H = MatrixXd::Zero(TOTAL_MEAS, TOTAL_STATES);
    
    for (int joint = 0; joint < NUM_JOINTS; ++joint) {
        int meas_base = joint * 3;
        int state_base = joint * STATES_PER_JOINT;
        
        H(meas_base + 0, state_base + 0) = 1.0;  // px
        H(meas_base + 1, state_base + 4) = 1.0;  // py
        H(meas_base + 2, state_base + 8) = 1.0;  // pz
    }
    return H;
}

// ============================================================================
// KALMAN FILTER
// ============================================================================

void runLKF(const vector<VectorXd>& measurements,
            const MatrixXd& F, const MatrixXd& Q,
            const MatrixXd& H, const MatrixXd& R,
            vector<VectorXd>& states) {
    
    int num_frames = measurements.size();
    states.reserve(num_frames);
    
    // Initialize
    VectorXd x = VectorXd::Zero(TOTAL_STATES);
    VectorXd z0 = measurements[0];
    for (int joint = 0; joint < NUM_JOINTS; ++joint) {
        x(joint * STATES_PER_JOINT + 0) = z0(joint * 3 + 0);
        x(joint * STATES_PER_JOINT + 4) = z0(joint * 3 + 1);
        x(joint * STATES_PER_JOINT + 8) = z0(joint * 3 + 2);
    }
    
    MatrixXd P = MatrixXd::Identity(TOTAL_STATES, TOTAL_STATES) * 1.0;
    MatrixXd Ht = H.transpose();
    
    for (int k = 0; k < num_frames; ++k) {
        if ((k + 1) % 500 == 0 || k == num_frames - 1) {
            cout << "\r  Frame " << (k + 1) << " / " << num_frames << flush;
        }
        
        // Predict
        x = F * x;
        P = F * P * F.transpose() + Q;
        
        // Update
        VectorXd z = measurements[k];
        VectorXd y = z - H * x;
        MatrixXd PHt = P * Ht;
        MatrixXd S = H * PHt + R;
        
        LLT<MatrixXd> llt(S);
        MatrixXd K = llt.solve(PHt.transpose()).transpose();
        
        x = x + K * y;
        
        MatrixXd I_KH = MatrixXd::Identity(TOTAL_STATES, TOTAL_STATES) - K * H;
        P = I_KH * P * I_KH.transpose() + K * R * K.transpose();
        
        states.push_back(x);
    }
    cout << endl;
}

// ============================================================================
// SAVE (NO HEADER)
// ============================================================================

void saveStates(const string& filename, const vector<VectorXd>& states) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "ERROR: Cannot write to " << filename << endl;
        exit(1);
    }
    
    file << fixed;
    file.precision(8);
    
    for (const auto& state : states) {
        for (int i = 0; i < state.size(); ++i) {
            file << state(i);
            if (i < state.size() - 1) file << ",";
        }
        file << "\n";
    }
    file.close();
}

// ============================================================================
// MAIN
// ============================================================================

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <input_csv> <output_csv>" << endl;
        return 1;
    }
    
    cout << "\n========================================" << endl;
    cout << "  LINEAR KALMAN FILTER - Milestone 2  " << endl;
    cout << "========================================\n" << endl;
    
    // Load (skips header)
    int num_frames = 0;
    vector<VectorXd> measurements = loadMeasurements(argv[1], num_frames);
    cout << "  Loaded " << num_frames << " data frames" << endl;
    
    if (num_frames != 3040) {
        cerr << "  ⚠ WARNING: Expected 3040 frames, got " << num_frames << endl;
    }
    
    // Parameters
    double dt = 0.01;
    double process_noise = 0.1;
    double meas_noise = 0.05;
    
    // Build matrices
    MatrixXd F = buildF(dt);
    MatrixXd Q = buildQ(dt, process_noise);
    MatrixXd H = buildH();
    MatrixXd R = MatrixXd::Identity(TOTAL_MEAS, TOTAL_MEAS) * (meas_noise * meas_noise);
    
    // Run filter
    cout << "\nRunning filter..." << endl;
    vector<VectorXd> states;
    
    auto start = chrono::high_resolution_clock::now();
    runLKF(measurements, F, Q, H, R, states);
    auto end = chrono::high_resolution_clock::now();
    
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "  Time: " << duration.count() << " ms" << endl;
    
    // Save (no header)
    saveStates(argv[2], states);
    cout << "  Saved " << states.size() << " rows to " << argv[2] << endl;
    cout << "\n=== SUCCESS ===\n" << endl;
    
    return 0;
}
