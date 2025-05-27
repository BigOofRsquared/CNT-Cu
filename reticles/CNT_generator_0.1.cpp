#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <armadillo>
#include <string>
#include <iomanip>
#include <numeric>

// C++17 needed

using namespace std;
using namespace arma;

int main(int argc, char** argv) {

    if (argc != 6) {
        cout << "Usage: " << argv[0] << " <n> <m> <repetitions> <output file name> <a>" << endl;
        return -1;
    }

    string filename{argv[4]};
    int n{atoi(argv[1])};
    int m{atoi(argv[2])};
    int repetitions{atoi(argv[3])};
    double a{atof(argv[5])};

    // actual nanotube 3d points
    vector<vec> points(0);

    // Bravais and basis vectors
    vec a1 = {1., 0.};
    vec a2 = {0.5, sqrt(3.) / 2.};
    vec b = (a1 + a2) * 2. / 3.;

    // axial indexes
    int nc{-(n + 2*m) / gcd(n, m) * repetitions};
    int mc{(m + 2*n) / gcd(n, m) * repetitions};

    // tricky condition for (n, n) nanotubes
    if (n == m) {
        nc/=3;
        mc/=3;
    }

    vec A = n * a1 + m * a2;
    double A2{arma::dot(A, A)};
    double sqrtA2{sqrt(A2)};

    cout << "width = " << sqrtA2 / 2. / M_PI << endl;

    vec C = nc * a1 + mc * a2;
    double C2{arma::dot(C, C)};
    double sqrtC2{sqrt(C2)};

    cout << "lenght = " << sqrtC2 << endl;
    
    vec candidateOnGraphene(2); // point on the graphene sheet, normalized
    vec pointOnTube(3); // the selected point, wrapped in the tube

    // grid extremal points
    int minN {min(min(0, n), min(nc, n + nc))},
        maxN {max(max(0, n), max(nc, n + nc))},
        minM {min(min(0, m), min(mc, m + mc))},
        maxM {max(max(0, m), max(mc, m + mc))};

    // points counter
    int counter{};

    // loop on possible m and n
    for (int n_i{minN}; n_i < maxN; n_i++) {
        for (int m_i{minM}; m_i < maxM; m_i++) {

            candidateOnGraphene = n_i * a1 + m_i * a2;
            double relativeProjectionA{arma::dot(candidateOnGraphene, A) / A2};
            double relativeProjectionC{arma::dot(candidateOnGraphene, C) / C2};

            // checking the presence
            if (relativeProjectionA >= -0.00001 && relativeProjectionA < 0.99999
                && relativeProjectionC >= -0.00001 && relativeProjectionC < 0.99999) 
            {

                relativeProjectionC += 0.0001;
                
                // pure Bravais point
                pointOnTube(0) = cos(relativeProjectionA * 2. * M_PI) * sqrtA2 / 2. / M_PI * a;
                pointOnTube(1) = sin(relativeProjectionA * 2. * M_PI) * sqrtA2 / 2. / M_PI * a;
                pointOnTube(2) = (relativeProjectionC - floor(relativeProjectionC) - 0.5) * sqrtC2 * a;
                points.push_back(pointOnTube);
                counter++;

                // base shifted point
                
                candidateOnGraphene = n_i * a1 + m_i * a2 + b;
                relativeProjectionA = arma::dot(candidateOnGraphene, A) / A2;
                relativeProjectionC = arma::dot(candidateOnGraphene, C) / C2 + 0.00001;
                

                pointOnTube(0) = cos(relativeProjectionA * 2. * M_PI) * sqrtA2 / 2. / M_PI * a;
                pointOnTube(1) = sin(relativeProjectionA * 2. * M_PI) * sqrtA2 / 2. / M_PI * a;
                pointOnTube(2) = (relativeProjectionC - floor(relativeProjectionC) - 0.5) * sqrtC2 * a;
                points.push_back(pointOnTube);
                counter++;
            }
        }
    }


    // feedback
    cout << "Total points generated: " << counter << endl;

    // printing on file
    ofstream fout(filename, ios::out);
    fout << counter << endl << "#nanotube (" << n << "," << m << "), pbc friendly" << endl;
    fout << scientific;
    fout << setprecision(10);
    for (auto v : points) {
        fout << "C" << setw(20) << v(0) << setw(20) << v(1) << setw(20) << v(2) << endl;
    }

    fout.close();

    return 0;
}