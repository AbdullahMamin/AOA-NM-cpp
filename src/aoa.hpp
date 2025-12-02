#ifndef AOA_HPP
#define AOA_HPP

#include <array>
#include <functional>
#include <random>

// Some small value to avoid divisions by 0
#define EPSILON (1.0e-3) 

// AOA-NM solver
template <int N, int n>
class Solver {
public:
    Solver(
        std::function<double(const std::array<double, n>&)> objective_function,
        double alpha,
        double mu,
        double mop_min,
        double mop_max,
        int max_iterations,
        const std::array<double, n> &lower_bounds,
        const std::array<double, n> &upper_bounds,
        double reflection_coeff,
        double expansion_coeff,
        double contraction_coeff,
        double shrinkage_coeff,
        int nm_max_iterations
    );
    std::array<double, n> solve();
    std::vector<std::pair<double, int>> getFitnessEvaluation();

private:
    // Objective function to maximize
    std::function<double(const std::array<double, n>&)> m_objective_function;

    // Parameters to AOA
    double m_alpha;
    double m_mu;
    double m_mop_min;
    double m_mop_max;
    int m_max_iterations;
    std::array<double, n> m_lower_bounds;
    std::array<double, n> m_upper_bounds;

    // Paramters to NM
    double m_reflection_coeff;
    double m_expansion_coeff;
    double m_contraction_coeff;
    double m_shrinkage_coeff;
    int m_nm_max_iterations;

    // Candidate solution list (array of solution-fitness pairs)
    std::array<std::pair<std::array<double, n>, double>, N> m_X;

    // Evaluation metrics:
    // - anytime objective_function is evaluated, it is tallied.
    // - best fitness per AOA iteration is kept track of alongside
    //   how many objective evaluations it took to get there.
    std::vector<std::pair<double, int>> m_best_fitness_vs_fitness_evaluations;
    int m_n_fitness_evaluations;
};

template <int N, int n>
Solver<N, n>::Solver(
    std::function<double(const std::array<double, n>&)> objective_function,
    double alpha,
    double mu,
    double mop_min,
    double mop_max,
    int max_iterations,
    const std::array<double, n> &lower_bounds,
    const std::array<double, n> &upper_bounds,
    double reflection_coeff,
    double expansion_coeff,
    double contraction_coeff,
    double shrinkage_coeff,
    int nm_max_iterations
) {
    m_objective_function = objective_function;
    m_alpha = alpha;
    m_mu = mu;
    m_mop_min = mop_min;
    m_mop_max = mop_max;
    m_max_iterations = max_iterations;
    m_lower_bounds = lower_bounds;
    m_upper_bounds = upper_bounds;
    m_reflection_coeff = reflection_coeff;
    m_expansion_coeff = expansion_coeff;
    m_contraction_coeff = contraction_coeff;
    m_shrinkage_coeff = shrinkage_coeff;
    m_nm_max_iterations = nm_max_iterations;
    m_n_fitness_evaluations = 0;
}

template <int N, int n>
std::array<double, n> Solver<N, n>::solve() {
    // Initialize RNG
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis_0_1(0.0, 1.0);

    // Initialize random population
    for (int j = 0; j < n; j++) {
        std::uniform_real_distribution<> dis_l_u(m_lower_bounds[j], m_upper_bounds[j]);
        for (int i = 0; i < N; i++) {
            m_X[i].first[j] = dis_l_u(gen);
        }
    }

    // Calculate fitness of the population
    for (int i = 0; i < N; i++) {
        m_X[i].second = m_objective_function(m_X[i].first);
        m_n_fitness_evaluations++;
    }

    // AOA-NM loop
    auto aoa_best = m_X[0];
    double inv_alpha = 1.0/m_alpha;
    for (int c_iter = 1; c_iter <= m_max_iterations; c_iter++) {
        // Find best solution
        for (int i = 0; i < N; i++) {
            if (m_X[i].second > aoa_best.second) {
                aoa_best = m_X[i];
            }
        }
        m_best_fitness_vs_fitness_evaluations.push_back(std::make_pair(aoa_best.second, m_n_fitness_evaluations));

        // Calculate MOA and MOP
        double MOA = m_mop_min + c_iter*(m_mop_max - m_mop_min)/m_max_iterations;
        double MOP = 1.0 - pow(c_iter, inv_alpha)/pow(m_max_iterations, inv_alpha);

        // AOA exploration/exploitation
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < n; j++) {
                double r1 = dis_0_1(gen);
                double r2 = dis_0_1(gen);
                double r3 = dis_0_1(gen);
                double lb = m_lower_bounds[j];
                double ub = m_upper_bounds[j];
                if (r1 > MOA) {
                    // Exploration
                    if (r2 > 0.5) {
                        // Division
                        m_X[i].first[j] = aoa_best.first[j]/(MOP + EPSILON)*((ub - lb)*m_mu + lb);
                    } else {
                        // Multiplication
                        m_X[i].first[j] = aoa_best.first[j]*MOP*((ub - lb)*m_mu + lb);
                    }
                } else {
                    // Exploitation
                    if (r3 > 0.5) {
                        // Subtraction
                        m_X[i].first[j] = aoa_best.first[j] - MOP*((ub - lb)*m_mu + lb);
                    } else {
                        // Addition
                        m_X[i].first[j] = aoa_best.first[j] + MOP*((ub - lb)*m_mu + lb);
                    }
                }
            }

            // Recalculte the fitness
            m_X[i].second = m_objective_function(m_X[i].first);
            m_n_fitness_evaluations++;
        }

        // The NM part of AOA-NM

        // First initialize the simplex from some of the candidate solutions
        std::array<std::pair<std::array<double, n>, double>, n + 1> simplex;
        for (int i = 0; i < n + 1; i++) {
            simplex[i] = m_X[i];
        }

        // NM iterations
        for (int k = 1; k <= m_nm_max_iterations; k++) {
            // Step 1: sort (get worst, best and second worst)
            int worst_idx = 0;
            auto worst = simplex[0];
            auto second_worst = simplex[0];
            int best_idx = 0;
            auto best = simplex[0];
            for (int i = 1; i < n + 1; i++) {
                auto current = simplex[i];
                if (current.second > best.second) {
                    best_idx = i;
                    best = current;
                }
                if (current.second < worst.second) {
                    second_worst = worst;
                    worst_idx = i;
                    worst = current;
                }
            }

            // Step 2: calculate centroid of all except worst
            std::pair<std::array<double, n>, double> centroid{{0.0}, 0.0};
            for (int i = 0; i < n + 1; i++) {
                if (i == worst_idx) {
                    continue;
                }
                for (int j = 0; j < n; j++) {
                    centroid.first[j] += simplex[i].first[j];
                }
            }
            for (int j = 0; j < n; j++) {
                centroid.first[j] /= (n + 1);
            }
            centroid.second = m_objective_function(centroid.first);
            m_n_fitness_evaluations++;

            // Step 3: reflection
            auto reflection = centroid;
            for (int j = 0; j < n; j++) {
                reflection.first[j] += m_reflection_coeff*(centroid.first[j] - worst.first[j]);
            }
            reflection.second = m_objective_function(reflection.first);
            m_n_fitness_evaluations++;

            // If reflection is better than second worst, but not better than the best, replace the worst with the reflection and go to step 1
            if (reflection.second > second_worst.second && reflection.second <= best.second) {
                simplex[worst_idx] = reflection;

                continue;
            }

            // Step 4: expansion
            // If the reflection is the best, we try to expand it and see if it is better
            if (reflection.second > best.second) {
                auto expansion = centroid;
                for (int j = 0; j < n; j++) {
                    expansion.first[j] += m_expansion_coeff*(reflection.first[j] - centroid.first[j]);
                }
                expansion.second = m_objective_function(expansion.first);
                m_n_fitness_evaluations++;

                // If it is indeed better, we replace the worst with it and go to step 1
                if (expansion.second > reflection.second) {
                    simplex[worst_idx] = expansion;

                    continue;
                }
            }

            // Step 5: contraction
            if (reflection.second > worst.second) {
                // Contracted point is on the outside
                auto contracted = centroid;
                for (int j = 0; j < n; j++) {
                    contracted.first[j] += m_contraction_coeff*(reflection.first[j] - centroid.first[j]);
                }
                contracted.second = m_objective_function(contracted.first);
                m_n_fitness_evaluations++;

                // If it is better than the reflected, we replace the worst with it and go to step 1
                if (contracted.second > reflection.second) {
                    simplex[worst_idx] = contracted;

                    continue;
                }
            } else if (reflection.second <= worst.second) {
                // Contracted point is on the inside
                auto contracted = centroid;
                for (int j = 0; j < n; j++) {
                    contracted.first[j] += m_contraction_coeff*(worst.first[j] - centroid.first[j]);
                }
                contracted.second = m_objective_function(contracted.first);
                m_n_fitness_evaluations++;

                // If it is better than the worst, we replace the worst with it and go to step 1
                if (contracted.second > worst.second) {
                    simplex[worst_idx] = contracted;

                    continue;
                }
            }

            // Step 6: shrinkage
            // Shrink all points towards the best, except for the best itself and go to step 1
            for (int i = 0; i < n + 1; i++) {
                if (i == best_idx) {
                    continue;
                }

                for (int j = 0; j < n; j++) {
                    simplex[i].first[j] = best.first[j] + m_shrinkage_coeff*(simplex[i].first[j] - best.first[j]);
                }
                simplex[i].second = m_objective_function(simplex[i].first);
                m_n_fitness_evaluations++;
            }
        }

        // If the best simplex value is better than the worst AOA value, replace the AOA one with the simplex one
        auto best_simplex = simplex[0];
        for (int i = 1; i < n + 1; i++) {
            if (simplex[i].second > best_simplex.second) {
                best_simplex = simplex[i];
            }
        }

        int worst_aoa_idx = 0;
        auto worst_aoa = m_X[0];
        for (int i = 1; i < N; i++) {
            if (m_X[i].second < worst_aoa.second) {
                worst_aoa_idx = i;
                worst_aoa = m_X[i];
            }
        }

        if (best_simplex.second > worst_aoa.second) {
            m_X[worst_aoa_idx] = best_simplex;
        }
    }

    // Return the best solution
    return aoa_best.first;
}


template <int N, int n>
std::vector<std::pair<double, int>> Solver<N, n>::getFitnessEvaluation() {
    return m_best_fitness_vs_fitness_evaluations;
}

#endif // AOA_HPP
