/*
 * Extended Kalman Filter - Milestone 2
 * Team ChaiGPT: Laiba (31736), Ameer (30535), Usman (29126), Arrayyan (30557)
 
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <chrono>

using namespace Eigen;
using namespace std;

const int    NUM_JOINTS       = 23;
const int    STATES_PER_JOINT = 12;
const int    TOTAL_STATES     = 276;
const int    TOTAL_MEAS       = 69;
const double PI = 3.14159265358979323846;

// MANUAL ARCTAN2  
double arctan2_manual(double y, double x) {
    if (x == 0.0 && y == 0.0) return 0.0;

    double abs_y = (y < 0) ? -y : y;
    double abs_x = (x < 0) ? -x : x;
    double angle;

    // Fast, highly accurate Minimax polynomial for atan(t) on [0, 1]
    auto atan_approx = [](double t) {
        double t2 = t * t;
        return t * (0.999866 - t2 * (0.3302995 - t2 * (0.180141 - t2 * (0.085133 - t2 * 0.0208351))));
    };

    if (abs_x >= abs_y) {
        double t = abs_y / abs_x;
        angle = atan_approx(t);
        if (x < 0) angle = PI - angle;
    } else {
        double t = abs_x / abs_y;
        angle = PI/2.0 - atan_approx(t);
        if (x < 0) angle = PI - angle;
    }
    
    return (y < 0) ? -angle : angle;
}

// CSV LOADING
vector<VectorXd> loadMeasurements(const string& fn, int& N) {
    vector<VectorXd> data;
    ifstream f(fn);
    if (!f.is_open()) { cerr << "ERROR: Cannot open " << fn << "\n"; exit(1); }
    string line;
    if (getline(f, line)) cout << "  Skipped header\n";
    while (getline(f, line)) {
        if (line.empty()) continue;
        stringstream ss(line);
        VectorXd m(TOTAL_MEAS); string v; int i = 0;
        while (getline(ss, v, ',') && i < TOTAL_MEAS) m(i++) = stod(v);
        if (i == TOTAL_MEAS) data.push_back(m);
    }
    N = (int)data.size();
    return data;
}

// CARTESIAN → SPHERICAL  (elevation phi = arctan2(pz, rho_xy))
VectorXd cartToSph(const VectorXd& c) {
    VectorXd s(TOTAL_MEAS);
    for (int j = 0; j < NUM_JOINTS; ++j) {
        int i = j*3;
        double px=c(i), py=c(i+1), pz=c(i+2);
        s(i)   = sqrt(px*px + py*py + pz*pz);
        s(i+1) = arctan2_manual(py, px);
        s(i+2) = arctan2_manual(pz, sqrt(px*px + py*py)); 
    }
    return s;
}

// SYSTEM MATRICES
MatrixXd buildF(double dt) {
    MatrixXd F = MatrixXd::Zero(TOTAL_STATES, TOTAL_STATES);
    double dt2=dt*dt, dt3=dt2*dt;
    Matrix4d Fa;
    Fa << 1,dt,dt2/2,dt3/6, 0,1,dt,dt2/2, 0,0,1,dt, 0,0,0,1;
    for (int j=0; j<NUM_JOINTS; ++j)
        for (int a=0; a<3; ++a)
            F.block<4,4>(j*12+a*4, j*12+a*4) = Fa;
    return F;
}

MatrixXd buildQ(double dt, double sig) {
    MatrixXd Q = MatrixXd::Zero(TOTAL_STATES, TOTAL_STATES);
    double dt2=dt*dt,dt3=dt2*dt,dt4=dt3*dt,dt5=dt4*dt,dt6=dt5*dt,s2=sig*sig;
    Matrix4d Qa;
    Qa << dt6/36,dt5/12,dt4/6,dt3/6,
          dt5/12,dt4/4, dt3/2,dt2/2,
          dt4/6, dt3/2, dt2,  dt,
          dt3/6, dt2/2, dt,   1.0;
    Qa *= s2;
    for (int j=0; j<NUM_JOINTS; ++j)
        for (int a=0; a<3; ++a)
            Q.block<4,4>(j*12+a*4, j*12+a*4) = Qa;
    return Q;
}

// NONLINEAR MEASUREMENT FUNCTION h(x)
VectorXd h_func(const VectorXd& x) {
    VectorXd z(TOTAL_MEAS);
    for (int j=0; j<NUM_JOINTS; ++j) {
        int sb=j*12, mb=j*3;
        double px=x(sb), py=x(sb+4), pz=x(sb+8);
        z(mb)   = sqrt(px*px + py*py + pz*pz);
        z(mb+1) = arctan2_manual(py, px);
        z(mb+2) = arctan2_manual(pz, sqrt(px*px+py*py)); 
    }
    return z;
}

// JACOBIAN ∂h/∂x
MatrixXd computeJacobian(const VectorXd& x) {
    MatrixXd H = MatrixXd::Zero(TOTAL_MEAS, TOTAL_STATES);

    // eps_rxy = 0.01 m instead of 1e-8.
    // theta entries = ±1/rho_xy^2.
    //   eps=1e-8  → max entry = 1e16  (causes enormous Kalman gain → spikes)
    //   eps=0.01  → max entry = 1e4   (large but bounded)
    const double eps_r   = 1e-3;
    const double eps_rxy = 1e-2;

    for (int j=0; j<NUM_JOINTS; ++j) {
        int sb=j*12, mb=j*3;
        double px=x(sb), py=x(sb+4), pz=x(sb+8);
        double rho_xy = sqrt(px*px + py*py);
        double r      = sqrt(rho_xy*rho_xy + pz*pz);
        if (r      < eps_r  ) r      = eps_r;
        if (rho_xy < eps_rxy) rho_xy = eps_rxy;
        double r2   = r*r;
        double rxy2 = rho_xy*rho_xy;  

        // ∂r/∂[px,py,pz]
        H(mb,   sb  ) = px/r;
        H(mb,   sb+4) = py/r;
        H(mb,   sb+8) = pz/r;

        // ∂θ/∂[px,py,pz]
        H(mb+1, sb  ) = -py/rxy2;   
        H(mb+1, sb+4) =  px/rxy2;
        H(mb+1, sb+8) =  0.0;

        // ∂φ/∂[px,py,pz]   for φ = arctan2(pz, rho_xy)
        H(mb+2, sb  ) = -(px*pz)/(r2*rho_xy);   // NEGATIVE
        H(mb+2, sb+4) = -(py*pz)/(r2*rho_xy);   // NEGATIVE
        H(mb+2, sb+8) =   rho_xy/r2;             // POSITIVE
    }
    return H;
}

// EKF
void runEKF(const vector<VectorXd>& meas_cart,
            const vector<VectorXd>& meas_sph,
            const MatrixXd& F, const MatrixXd& Q, const MatrixXd& R,
            vector<VectorXd>& states) {

    int N = (int)meas_cart.size();
    states.reserve(N);

    // Initialise from first Cartesian measurement
    VectorXd x = VectorXd::Zero(TOTAL_STATES);
    for (int j=0; j<NUM_JOINTS; ++j) {
        x(j*12  ) = meas_cart[0](j*3  );
        x(j*12+4) = meas_cart[0](j*3+1);
        x(j*12+8) = meas_cart[0](j*3+2);
    }
    MatrixXd P = MatrixXd::Identity(TOTAL_STATES, TOTAL_STATES);

    const double skip_threshold = 1e6; // Only catch absolute divergence

    int skipped = 0;

    for (int k=0; k<N; ++k) {
        if ((k+1)%500==0 || k==N-1)
            cout << "\r  Frame " << (k+1) << "/" << N << flush;

        x = F * x;
        P = F * P * F.transpose() + Q;

   
        VectorXd y = meas_sph[k] - h_func(x);

        // Wrap angles to [-pi, pi]
        for (int j=0; j<NUM_JOINTS; ++j) {
            int mb = j*3;
            while (y(mb+1) >  PI) y(mb+1) -= 2*PI;
            while (y(mb+1) < -PI) y(mb+1) += 2*PI;
            while (y(mb+2) >  PI) y(mb+2) -= 2*PI;
            while (y(mb+2) < -PI) y(mb+2) += 2*PI;
        }

        
        double innov_norm = y.norm();
        if (innov_norm > skip_threshold) {
            ++skipped;
            states.push_back(x);   // output the predict-only state
            continue;
        }

   
        MatrixXd Hk  = computeJacobian(x);
        MatrixXd PHt = P * Hk.transpose();
        MatrixXd S   = Hk * PHt + R;

        LDLT<MatrixXd> ldlt(S);

        if (ldlt.info() != Success) {
            ++skipped;
            states.push_back(x);
            continue;
        }

        // K = (S^{-1} * H * P)^T  solved without explicit inversion
        MatrixXd K = ldlt.solve(PHt.transpose()).transpose();  // 276×69
        x = x + K * y;

        // Joseph form: guarantees P stays positive semi-definite
        MatrixXd IKH = MatrixXd::Identity(TOTAL_STATES, TOTAL_STATES) - K*Hk;
        P = IKH * P * IKH.transpose() + K * R * K.transpose();

        states.push_back(x);
    }

    cout << "\n";
    if (skipped > 0)
        cout << "  Skipped " << skipped << " frames ("
             << (100.0*skipped/N) << "%) due to large innovations or decomp failure\n";
}

// SAVE
void saveStates(const string& fn, const vector<VectorXd>& states) {
    ofstream f(fn);
    if (!f.is_open()) { cerr << "ERROR: Cannot write " << fn << "\n"; exit(1); }
    f << fixed; f.precision(8);
    for (const auto& s : states) {
        for (int i=0; i<s.size(); ++i) { f<<s(i); if(i<s.size()-1) f<<","; }
        f << "\n";
    }
}

// MAIN
int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <input.csv> <output.csv>\n";
        return 1;
    }

    cout<< " EXTENDED KALMAN FILTER - Milestone 2\n";
       

    int N = 0;
    auto meas_cart = loadMeasurements(argv[1], N);
    cout << "  Loaded " << N << " frames\n";
    if (N != 3040) cerr << "  WARNING: Expected 3040, got " << N << "\n";

    cout << "Converting to spherical...\n";
    vector<VectorXd> meas_sph;
    meas_sph.reserve(N);
    for (auto& m : meas_cart) meas_sph.push_back(cartToSph(m));

    double dt  = 0.01;   // 100 Hz
    double sig = 0.1;    // process noise sigma_jerk

    // block-diagonal R — r in metres, angles in radians
    // sigma_r   = 0.05 m   (50 mm range noise)
    // sigma_ang = 0.005 rad (≈ 5 mm at r=1 m, i.e. 0.29 degrees)
    double sigma_r   = 0.05;
    double sigma_ang = 0.005;
    MatrixXd R = MatrixXd::Zero(TOTAL_MEAS, TOTAL_MEAS);
    for (int j=0; j<NUM_JOINTS; ++j) {
        R(j*3,   j*3  ) = sigma_r   * sigma_r;
        R(j*3+1, j*3+1) = sigma_ang * sigma_ang;
        R(j*3+2, j*3+2) = sigma_ang * sigma_ang;
    }

    MatrixXd F = buildF(dt);
    MatrixXd Q = buildQ(dt, sig);

    cout << "Running EKF...\n";
    vector<VectorXd> states;
    auto t0 = chrono::high_resolution_clock::now();
    runEKF(meas_cart, meas_sph, F, Q, R, states);
    auto t1 = chrono::high_resolution_clock::now();
    cout << "  Time: "
         << chrono::duration_cast<chrono::milliseconds>(t1-t0).count()
         << " ms\n";

    saveStates(argv[2], states);
    cout << "  Saved " << states.size() << " rows to " << argv[2] << "\n"
         << "SUCCESS\n";
    return 0;
}
