#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream>
#include "aoa.hpp"

// How to use AOA-NM code:
// First, define the problem, and create a function for the objective. Make sure to penalize out of bounds solutions.
// Second, initialize a Solver object with template arguments for the population size and solution dimension.
// - The arguments for the solver will be in order:
//   - the objective function
//   - alpha: controlling parameter
//   - mu: controlling parameter
//   - min MOP
//   - max MOP
//   - AOA iterations
//   - Lower bounds of variables
//   - Upper bounds of variables
//   - alpha (for NM)
//   - gamma (for NM)
//   - rho (for NM)
//   - sigma (for NM)
//   - NM iterations
// Finally, call the solve() method on the Solver object to get the AOA-NM solution

// a, b, c, e, f, l, delta
const int solution_dimension = 7;
const int a_idx = 0;
const int b_idx = 1;
const int c_idx = 2;
const int e_idx = 3;
const int f_idx = 4;
const int l_idx = 5;
const int delta_idx = 6;

// Constants set by paper
const double Y_MIN = 50;
const double Y_G = 150;
const double Y_MAX = 100;
const int Z_MAX = 100;
const double P = 100;

// Bounds set by paper
const std::array<double, solution_dimension> LB = {10, 10, 100, 0, 10, 100, 1.0};
const std::array<double, solution_dimension> UB = {150, 150, 200, 50, 150, 300, 3.14};

// F_k(x,z) from the paper
double fk(const std::array<double, solution_dimension> &x, double z) {
    double a = x[a_idx];
    double b = x[b_idx];
    double c = x[c_idx];
    double e = x[e_idx];
    double l = x[l_idx];
    double g = sqrt((l+z)*(l+z) + e*e);
    double phi = atan(e/(l-z));
    double alpha = acos((a*a + g*g - b*b)/(2*a*g)) - phi;
    double beta = acos((b*b + g*g - a*a)/(2*b*g)) - phi;
    return (P*b*sin(alpha + beta))/(2*c*cos(alpha));
}

// y(x,z) from the paper
double y(const std::array<double, solution_dimension> &x, double z) {
    double a = x[a_idx];
    double b = x[b_idx];
    double c = x[c_idx];
    double e = x[e_idx];
    double f = x[f_idx];
    double l = x[l_idx];
    double delta = x[delta_idx];
    double g = sqrt((l+z)*(l+z) + e*e);
    double phi = atan(e/(l-z));
    double beta = acos((b*b + g*g - a*a)/(2*b*g)) - phi;
    return 2*(e + f + c*sin(beta + delta));
}

// The fitness function which weighs between two objectives
double fitness(const std::array<double, solution_dimension> &x, double w1, double w2) {
    double a = x[a_idx];
    double b = x[b_idx];
    double c = x[c_idx];
    double e = x[e_idx];
    double f = x[f_idx];
    double l = x[l_idx];
    double delta = x[delta_idx];

    // Apply fuzzy rule
    // if ((a < 4*b) && (c < a + b)) {
        // f = 2*e + 10;
    // }
    // if ((a < 4*b) && (c > a + b)) {
        // f = e + 50;
    // }

    std::array<double, solution_dimension> new_x = {a, b, c, e, f, l, delta};

    // Penalize out of bounds
    if (
        a < LB[a_idx] || a > UB[a_idx] ||
        b < LB[b_idx] || b > UB[b_idx] ||
        c < LB[c_idx] || c > UB[c_idx] ||
        e < LB[e_idx] || e > UB[e_idx] ||
        f < LB[f_idx] || f > UB[f_idx] ||
        l < LB[l_idx] || l > UB[l_idx] ||
        delta < LB[delta_idx] || delta > UB[delta_idx]
    ) {
        return -1.0e10;
    }

    // Penalize breaking constraints
    double y_x_zmax = y(new_x, Z_MAX);
    double y_x_0 = y(new_x, 0);
    if (
        (Y_MIN - y_x_zmax) < 0 ||
        (y_x_zmax) < 0 ||
        (y_x_0 - Y_MAX) < 0 ||
        (Y_G - y_x_0) < 0 ||
        ((a + b)*(a + b) - l*l - e*e) < 0 ||
        ((l - Z_MAX)*(l - Z_MAX) + (a - e)*(a - e) - b*b) < 0 ||
        (l - Z_MAX) < 0
    ) {
        return -1.0e10;
    }

    // Calculate f1, f2 and then return f from the weighted sum
    double min_fk = INFINITY;
    double max_fk = -INFINITY;
    for (int z = 0; z <= Z_MAX; z++) {
        double fk_xz = fk(new_x, z);
        if (fk_xz < min_fk) {
            min_fk = fk_xz;
        }
        if (fk_xz > max_fk) {
            max_fk = fk_xz;
        }
    }
    double f1 = fabs(max_fk - min_fk);
    double f2 = P/min_fk;
    // We negate it since our solver is a maximizer and the fitness is a min
    return -(w1*f1/1000 + w2*f2/10);
}

int main() {
    const int population_size = 1000;
    assert(population_size > solution_dimension); // We need at least n + 1 for the population size

    std::ofstream file1("gripper_sol.csv");

    file1 << "w1,w2,fitness,a,b,c,e,f,l,delta\n";

    // Testing on robot gripper problem and printing to csv
    double w1_s[] = {1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0};
    for (double w1 : w1_s) {
        double w2 = 1.0 - w1;
        
        auto fit = [&](const std::array<double, solution_dimension> &x) {return fitness(x, w1, w2);};
        Solver<population_size, solution_dimension> solver(
            fit, // objective function to maximize
            10.0,    // alpha
            0.5,     // mu
            0.2,     // min MOP
            0.8,     // max MOP
            1000,    // max iterations
            LB,      // lower bounds
            UB,      // upper bounds
            1.0,     // NM alpha
            2.0,     // NM gamma
            0.5,     // NM rho
            0.5,     // NM sigma
            1000      // NM max iterations
        );

        // Get best solution and add it to csv
        auto solution = solver.solve();

        std::cout << "Finished w1 = " << w1 << '\n';
        file1 << w1 << ',' << w2 << ',' << fit(solution) << ',' << solution[0] << ',' << solution[1] << ',' << solution[2] << ',' << solution[3] << ',' << solution[4] << ',' << solution[5] << ',' << solution[6] << '\n';
    }

    std::ofstream file2("beale_sol.csv");
    std::ofstream file3("beale_eval.csv");

    // Testing on Beale function and printing results
    auto fit = [](const std::array<double, 2> &v) {
        double x = v[0];
        double y = v[1];
        double a = 1.5 - x + x*y;
        double b = 2.25 - x + x*y*y;
        double c = 2.625 - x + x*y*y*y;
        return -(a*a + b*b + c*c); // Negative since its a min
    };
    Solver<1000, 2> solver(
        fit,     // objective function to maximize
        10.0,    // alpha
        0.5,     // mu
        0.2,     // min MOP
        0.8,     // max MOP
        10,     // max iterations
        {-2, 2}, // lower bounds
        {-2, 2}, // upper bounds
        1.0,     // NM alpha
        2.0,     // NM gamma
        0.5,     // NM rho
        0.5,     // NM sigma
        1000     // NM max iterations
    );

    auto solution = solver.solve();
    file2 << "fitness,x,y\n";
    file2 << fit(solution) << ',' << solution[0] << ',' << solution[1] << '\n';

    file3 << "best_fitness,n_fitness_evals\n";
    for (auto pair : solver.getFitnessEvaluation()) {
        file3 << pair.first << ',' << pair.second << '\n';
    }

    return 0;
}
