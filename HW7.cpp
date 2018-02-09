//Preston Peck
//CS 365
//November 20, 2017
//HW7

#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
using namespace std;

int binomial_simple(double S, double K, double r, double q, double sigma, double T, double t0, bool call, bool American, int n, double & V, bool print);
void printTree(double** stock_nodes, double** option_nodes, int n);

int main() {
    double S, K, r, q, v, T, t0, fair = 0.0;
    int n = 0;

    cout << "7.10 Tests" << endl;
    S = 100;
    cout << "S: $" << S << endl;
    K = 100;
    cout << "K: $" << K << endl;
    r = 0.1;
    cout << "r: " << r << endl;
    q = 0.0;
    cout << "q: " << q << endl;
    v = 0.5;
    cout << "v: " << v << endl;
    T = 0.3;
    cout << "T: " << T << endl;
    t0 = 0.0;
    cout << "t0: " << t0 << endl;
    n = 3;
    cout << "n: " << n << endl << endl;

    binomial_simple(S, K, r, q, v, T, t0, true, false, n, fair, true);
    binomial_simple(S, K, r, q, v, T, t0, false, false, n, fair, true);
    binomial_simple(S, K, r, q, v, T, t0, true, true, n, fair, true);
    binomial_simple(S, K, r, q, v, T, t0, false, true, n, fair, true);

    cout << endl << "7.11 New calculations" << endl;
    r = 0.1;
    cout << "r: " << r << endl;
    q = 0.1;
    cout << "q: " << q << endl;
    T = 1;
    cout << "T: " << T << endl;
    n = 100;
    cout << "n: " << n << endl;

    binomial_simple(S, K, r, q, v, T, t0, true, false, n, fair, false);
    binomial_simple(S, K, r, q, v, T, t0, false, false, n, fair, false);
    binomial_simple(S, K, r, q, v, T, t0, true, true, n, fair, false);
    binomial_simple(S, K, r, q, v, T, t0, false, true, n, fair, false);

    cout << endl << "T (> t0): ";
    cin >> T;

    while (S <= 0) {
        cout << "T > t0: ";
        cin >> T;
    }

    cout << "S (> 0): ";
    cin >> S;

    while (S <= 0) {
        cout << "S > 0: ";
        cin >> S;
    }

    K = S;
    cout << "K: " << K << endl;

    cout << "r (> 0): ";
    cin >> r;

    while (r <= 0) {
        cout << "r > 0: ";
        cin >> r;
    }

    q = r;
    cout << "q: " << q << endl << endl;

    binomial_simple(S, K, r, q, v, T, t0, true, false, n, fair, false);
    binomial_simple(S, K, r, q, v, T, t0, false, false, n, fair, false);
    binomial_simple(S, K, r, q, v, T, t0, true, true, n, fair, false);
    binomial_simple(S, K, r, q, v, T, t0, false, true, n, fair, false);
}

//7.1 Function signature
int binomial_simple(double S, double K, double r, double q, double sigma, double T, double t0, bool call, bool American, int n, double & V, bool print) {
    //7.2 Validation tests
    if (n < 1 || S <= 0 || T <= t0 || sigma <= 0.0) {
        return 1;
    }

    //7.3 Parameters
    double dt = (T - t0) / double(n);
    double df = exp(-r * dt);
    double growth = exp((r - q) * dt);

    double u = exp(sigma * sqrt(dt));
    double d = 1.0 / u;
    double p_prob = (growth - d) / (u - d);
    double q_prob = 1.0 - p_prob;

    if (p_prob < 0.0 || p_prob > 1.0) {
        return 1;
    }

    //7.4 Allocate memory/ set up arrays
    //array
    //allocate memory
    double** stock_nodes = new double* [n + 1];
    double** option_nodes = new double* [n + 1];
    double* S_tmp = NULL;
    double* V_tmp = NULL;

    for (int i = 0; i <= n; ++i) {
        stock_nodes[i] = new double[n + 1];
        option_nodes[i] = new double[n + 1];

        S_tmp = stock_nodes[i];
        V_tmp = option_nodes[i];

        for (int j = 0; j <= n; ++j) {
            S_tmp[j] = 0;
            V_tmp[j] = 0;
        }
    }

    //7.5 Set up stock prices in nodes
    S_tmp = stock_nodes[0];
    S_tmp[0] = S;

    for (int i = 1; i <= n; ++i) {
        double* prev = stock_nodes[i - 1];
        S_tmp = stock_nodes[i];
        S_tmp[0] = prev[0] * d;

        for (int j = 1; j <= n; ++j) {
            S_tmp[j] = S_tmp[j - 1] * u * u;
        }
    }
    
    //7.6 Terminal payoff
    int i = n;
    S_tmp = stock_nodes[i];
    V_tmp = option_nodes[i];

    for (int j = 0; j <= n; ++j) {
        double intrinsic = 0;
        
        if (call) {
            if (S_tmp[j] > K) {
                intrinsic = S_tmp[j] - K;
            }
        }
        
        else {
            if (S_tmp[j] < K) {
                intrinsic = K - S_tmp[j];
            }
        }

        V_tmp[j] = intrinsic;
    }
    
    //7.7 Main valuation loop
    for (int i = n - 1; i >= 0; --i) {
        S_tmp = stock_nodes[i];
        V_tmp = option_nodes[i];
        double* V_next = option_nodes[i + 1];

        for (int j = 0; j <= i; ++j) {
            V_tmp[j] = df * (p_prob * V_next[j+1] + q_prob * V_next[j]);
            
            // early exercise test
            if (American) {
                if (call) {
                    V_tmp[j] = fmax(V_tmp[j], fmax(S_tmp[j] - K, 0));
                }

                else {
                    V_tmp[j] = fmax(V_tmp[j], fmax(K - S_tmp[j], 0));
                }
            }
        } 
    }

    if (print) {
        printTree(stock_nodes, option_nodes, n);
    }

    //7.8 option fair value
    i = 0;
    V_tmp = option_nodes[i];
    V = V_tmp[0];

    if (American == false && call == true) {
        cout << "European CALL: $";
    }

    else if (American == false && call == false) {
        cout << "European PUT: $"; 
    }

    else if (American == true && call == true) {
        cout << "American CALL: $";
    }

    else {
        cout << "American PUT: $";
    }

    cout << V << endl;

    if (print) {
        cout << endl;
    }
        
    //7.9 Memory deallocation
    for (int i = 0; i <= n; ++i) {
        delete[] stock_nodes[i];
        delete[] option_nodes[i];
    }

    delete[] stock_nodes;
    delete[] option_nodes;
    return 0;
}

void printTree(double** stock_nodes, double** option_nodes, int n) {
    cout << "[stock]" << endl;
    cout << "(option)" << endl << endl;

    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= i; j++) {
            cout << "[" << stock_nodes[i][j] << "]" << '\t';

            ostringstream ss;
            ss << stock_nodes[i][j];
            string st = ss.str();

            if (st.length() <= 3) {
                cout << '\t';
            }
        }
        cout << endl;

        for (int k = 0; k <= i; k++) {
            cout << "(" << option_nodes[i][k] << ")" << '\t';

            ostringstream ss;
            ss << option_nodes[i][k];
            string st = ss.str();

            if (st.length() <= 5) {
                cout << '\t';
            }
        }

        if (i != n) {
            cout << endl;
            cout << '\t' << "|-------.";

            for (int l = 0; l < i; l++) {
                cout << '\t' << "|-------.";
            }
        }

        cout << endl;
    }
}
