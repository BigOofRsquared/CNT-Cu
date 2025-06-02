#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <armadillo>
#include <string>
#include <iomanip>
#include <numeric>



//   "\n\nCNT+METAL RETICLE GENERATOR - Riccardo Piazza\n\n"
//   "This code generates carbon nanotubes .xyz PBC-friendly ions positioning files.\n"
//   "It is possible to generate every carbon nanotube and include extra atoms"
//   "inside or outside the surface.\n\n\n\n"

//   "Usage:" << argv[0] << "\n"
//   "[file name] [n chiral index] [m chiral index] [number of cell repetition]"
//   "[Bravais lattice vector lenght] [x-y simulation cell width] [metal position]\n\n"
//   "and for every metal atom position specify\n\n"
//    "[Xx element name] [distance from the graphene sheet (negative is inside, positive is outside)] [chiral index n] [chiral index m]\n\n"
//    Metal (n, m) indexes can also be non-integer (they are general points on the sheet).

// C++17 needed

using namespace std;
using namespace arma;

int main(int argc, char** argv) {

    vector<vec> metalIndexes;
    int nMetals{argc - 8};
    vec distanceFromSurface;
    vector<string> names;

    if (argc < 7 || nMetals % 4 != 0 || (nMetals) / 4 != atoi(argv[7])) {
       cout << "+     \n\n      CNT+METAL RETICLE GENERATOR - Riccardo Piazza\n\n"
            << "+     This code generates carbon nanotubes .xyz PBC-friendly ions positioning files.\n"
            << "+     It is possible to generate every carbon nanotube and include extra atoms"
            << " inside or outside the surface.\n\n\n\n"
            << "+     Usage:" << argv[0] << "\n"
            << "+     [file name] [n chiral index] [m chiral index] [number of cell repetition]\n"
            << "+     [Bravais lattice vector lenght] [x-y simulation cell width] [number of metal ions]\n\n"
            << "+     and for every metal atom position specify\n\n"
            << "+     [Xx element name] [distance from the graphene sheet (negative is inside, positive is outside)]\n"
            << "+     [chiral index n] [chiral index m]\n\n"
            << "+     Metal (n, m) indexes can also be non-integer (they are general points on the sheet).\n\n\n" << endl;        
        return -1;
    } else if ((nMetals) / 4 == atoi(argv[7]) && (nMetals) % 4 == 0){
        nMetals /= 4;
        distanceFromSurface.set_size(nMetals);

        // metal site name, distance and indexes on graphene
        for (uint i{}; i < nMetals; i++) {
            names.push_back(argv[8 + 2*i]);
            distanceFromSurface(i) = atof(argv[8 + 2*i + 1]);

            vec index(2);
            index(0) = atof(argv[8 + 2*i + 2]);
            index(1) = atof(argv[8 + 2*i + 3]);
            metalIndexes.push_back(index);
        }
    }

    // file
    string filename{argv[1]};

    // tube
    int n{atoi(argv[2])};
    int m{atoi(argv[3])};
    int repetitions{atoi(argv[4])};
    double a{atof(argv[5])};

    // pbc cell width
    double width{atof(argv[6])};    

    // actual nanotube 3d points
    vector<vec> points(0);
    // metal sites 3d points
    vector<vec> extraSites(0);

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

    cout << "Tube width = " << sqrtA2 / 2. / M_PI * a<< endl;

    vec C = nc * a1 + mc * a2;
    double C2{arma::dot(C, C)};
    double sqrtC2{sqrt(C2)};

    cout << "Tube lenght = " << sqrtC2 * a << endl;
    
    vec candidateOnGraphene(2); // point on the graphene sheet, normalized
    vec pointOnTube(3); // the selected, temporary, point wrapped in the tube

    // grid extremal points
    int minN {min(min(0, n), min(nc, n + nc))},
        maxN {max(max(0, n), max(nc, n + nc))},
        minM {min(min(0, m), min(mc, m + mc))},
        maxM {max(max(0, m), max(mc, m + mc))};

    // points counter
    int counter{};
    double relativeProjectionA{}, relativeProjectionC{};

    // loop on possible m and n
    for (int n_i{minN}; n_i < maxN; n_i++) {
        for (int m_i{minM}; m_i < maxM; m_i++) {

            candidateOnGraphene = n_i * a1 + m_i * a2;
            relativeProjectionA = arma::dot(candidateOnGraphene, A) / A2;
            relativeProjectionC = arma::dot(candidateOnGraphene, C) / C2;

            // checking the presence
            if (relativeProjectionA >= -0.00001 && relativeProjectionA < 0.99999
                && relativeProjectionC >= -0.00001 && relativeProjectionC < 0.99999) 
            {
                relativeProjectionC += 0.00001;
                
                // pure Bravais point
                pointOnTube(0) = (cos(relativeProjectionA * 2. * M_PI)) * sqrtA2 / 2. / M_PI * a + width/2.;
                pointOnTube(1) = (sin(relativeProjectionA * 2. * M_PI)) * sqrtA2 / 2. / M_PI * a + width/2.;
                pointOnTube(2) = (relativeProjectionC - floor(relativeProjectionC)) * sqrtC2 * a;
                points.push_back(pointOnTube);
                counter++;

                // base shifted point                
                candidateOnGraphene = n_i * a1 + m_i * a2 + b;
                relativeProjectionA = arma::dot(candidateOnGraphene, A) / A2;
                relativeProjectionC = arma::dot(candidateOnGraphene, C) / C2 + 0.00001;

                pointOnTube(0) = cos(relativeProjectionA * 2. * M_PI) * sqrtA2 / 2. / M_PI * a + width/2.;
                pointOnTube(1) = sin(relativeProjectionA * 2. * M_PI) * sqrtA2 / 2. / M_PI * a + width/2.;
                pointOnTube(2) = (relativeProjectionC - floor(relativeProjectionC)) * sqrtC2 * a;
                points.push_back(pointOnTube);
                counter++;
            }
        }
    }

    // adding metal atoms
    for (uint i{}; i < nMetals; i++) {

        candidateOnGraphene = metalIndexes[i](0) * a1 + metalIndexes[i](1) * a2;
        relativeProjectionA = arma::dot(candidateOnGraphene, A) / A2;
        relativeProjectionC = arma::dot(candidateOnGraphene, C) / C2;
        relativeProjectionC += 0.00001;

        pointOnTube(0) = (cos(relativeProjectionA * 2. * M_PI)) * ( sqrtA2 * a / 2. / M_PI + distanceFromSurface(i) ) + width/2.;
        pointOnTube(1) = (sin(relativeProjectionA * 2. * M_PI)) * ( sqrtA2 * a / 2. / M_PI + distanceFromSurface(i) ) + width/2.;
        pointOnTube(2) = (relativeProjectionC - floor(relativeProjectionC)) * sqrtC2 * a;
        extraSites.push_back(pointOnTube);
    }


    // feedback
    cout << "Total points generated: " << counter << " + " << nMetals << endl;

    // printing on file
    ofstream fout(filename, ios::out);
    fout << counter + nMetals << endl << "#nanotube (" << n << "," << m << ")  with copper atoms, pbc friendly" << endl;
    fout << scientific;
    fout << setprecision(10);

    // nanotube
    for (auto v : points) {
        fout << "C" << setw(20) << v(0) << setw(20) << v(1) << setw(20) << v(2) << endl;
    }

    // atoms
    for (uint i{}; i < nMetals; i++) {
        fout << names[i] << setw(20) << extraSites[i](0)
             << setw(20) << extraSites[i](1)
             << setw(20) << extraSites[i](2) << endl;
    }

    fout.close();

    return 0;
}