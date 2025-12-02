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
// Optionally, call the getFitnessEvaluation() method on the Solver object to get
// the fitness of best solution vs. number of fitness calls evaluation metric.

int main() {
    std::ofstream file1("beale_sol.csv");
    std::ofstream file2("beale_eval.csv");

    // Testing on Beale function and printing results
    auto fit = [](const std::array<double, 2> &v) {
        double x = v[0];
        double y = v[1];
        double a = 1.5 - x + x*y;
        double b = 2.25 - x + x*y*y;
        double c = 2.625 - x + x*y*y*y;
        return -(a*a + b*b + c*c); // Negative since its a min
    };
    Solver<1000, 2> solver( // 1000 population size, 2 solution dimension
        fit,     // objective function to maximize
        10.0,    // alpha
        0.5,     // mu
        0.2,     // min MOP
        0.8,     // max MOP
        10,      // max iterations
        {-2, 2}, // lower bounds
        {-2, 2}, // upper bounds
        1.0,     // NM alpha
        2.0,     // NM gamma
        0.5,     // NM rho
        0.5,     // NM sigma
        1000     // NM max iterations
    );

    auto solution = solver.solve();
    file1 << "fitness,x,y\n";
    file1 << fit(solution) << ',' << solution[0] << ',' << solution[1] << '\n';

    file2 << "best_fitness,n_fitness_evals\n";
    for (auto pair : solver.getFitnessEvaluation()) {
        file2 << pair.first << ',' << pair.second << '\n';
    }

    return 0;
}
