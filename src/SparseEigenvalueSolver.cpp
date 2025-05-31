#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/ArpackSupport>
#include <vector>
#include <complex>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <iostream>
#include <cmath>

struct EigenResults {
    Eigen::VectorXd eigenvalues;
    Eigen::MatrixXd eigenvectors;
    std::vector<double> residual_norms;
    int iterations;
    bool converged;
    double computation_time;
};

class SparseEigenvalueSolver {
private:
    Eigen::SparseMatrix<double> system_matrix;
    Eigen::SparseMatrix<double> mass_matrix;
    int krylov_dimension;
    double tolerance;
    int max_iterations;
    double spectral_shift;
    bool use_preconditioning;
    
    // Preconditioner for improved convergence
    Eigen::IncompleteILUT<double> preconditioner;
    
public:
    SparseEigenvalueSolver(int krylov_dim = 150, double tol = 1e-12, int max_iter = 1000) 
        : krylov_dimension(krylov_dim), tolerance(tol), max_iterations(max_iter), 
          spectral_shift(0.0), use_preconditioning(true) {}
    
    void setMatrices(const Eigen::SparseMatrix<double>& A, const Eigen::SparseMatrix<double>& M) {
        system_matrix = A;
        mass_matrix = M;
        
        if (use_preconditioning) {
            setupPreconditioner();
        }
    }
    
    void setSpectralShift(double sigma) {
        spectral_shift = sigma;
    }
    
    EigenResults solveGeneralizedEigenvalue(int num_eigenvalues) {
        auto start_time = std::chrono::high_resolution_clock::now();
        
        EigenResults results;
        
        if (spectral_shift != 0.0) {
            results = solveWithSpectralTransformation(num_eigenvalues);
        } else {
            results = solveDirect(num_eigenvalues);
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        results.computation_time = std::chrono::duration<double>(end_time - start_time).count();
        
        return results;
    }
    
    EigenResults solveWithSpectralTransformation(int num_eigenvalues) {
        // Shift-invert transformation: solve (A - σM)^{-1}M v = θ v
        // where θ = 1/(λ - σ) and λ are the original eigenvalues
        
        Eigen::SparseMatrix<double> shifted_matrix = system_matrix - spectral_shift * mass_matrix;
        
        // Set up iterative solver for shift-invert
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteILUT<double>> solver;
        solver.setTolerance(tolerance * 0.1);
        solver.setMaxIterations(max_iterations / 10);
        
        if (use_preconditioning) {
            solver.setPreconditioner(preconditioner);
        }
        
        solver.compute(shifted_matrix);
        
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Failed to factorize shifted matrix");
        }
        
        // Arnoldi iteration for the transformed problem
        return arnoldiIteration(solver, num_eigenvalues);
    }
    
    EigenResults solveDirect(int num_eigenvalues) {
        // Direct eigenvalue solution using ARPACK
        Eigen::ArpackGeneralizedSelfAdjointEigenSolver<Eigen::SparseMatrix<double>> solver;
        
        solver.compute(system_matrix, mass_matrix, num_eigenvalues, "SM"); // Smallest magnitude
        
        EigenResults results;
        results.eigenvalues = solver.eigenvalues();
        results.eigenvectors = solver.eigenvectors();
        results.converged = (solver.info() == Eigen::Success);
        results.iterations = solver.getNbrIterations();
        
        // Compute residuals
        computeResiduals(results);
        
        return results;
    }
    
    EigenResults arnoldiIteration(Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteILUT<double>>& solver, 
                                 int num_eigenvalues) {
        
        int n = system_matrix.rows();
        int k = std::min(num_eigenvalues, krylov_dimension);
        
        // Initialize random starting vector
        Eigen::VectorXd v0 = Eigen::VectorXd::Random(n);
        v0.normalize();
        
        // Krylov subspace construction
        Eigen::MatrixXd V = Eigen::MatrixXd::Zero(n, k + 1);
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(k + 1, k);
        
        V.col(0) = v0;
        
        for (int j = 0; j < k; ++j) {
            // Apply operator: w = (A - σM)^{-1} M v_j
            Eigen::VectorXd Mv_j = mass_matrix * V.col(j);
            Eigen::VectorXd w = solver.solve(Mv_j);
            
            if (solver.info() != Eigen::Success) {
                throw std::runtime_error("Linear solver failed in Arnoldi iteration");
            }
            
            // Modified Gram-Schmidt orthogonalization
            for (int i = 0; i <= j; ++i) {
                H(i, j) = V.col(i).dot(w);
                w -= H(i, j) * V.col(i);
            }
            
            H(j + 1, j) = w.norm();
            
            if (H(j + 1, j) < tolerance) {
                // Lucky breakdown - exact invariant subspace found
                k = j + 1;
                break;
            }
            
            if (j < k - 1) {
                V.col(j + 1) = w / H(j + 1, j);
            }
        }
        
        // Solve reduced eigenvalue problem
        Eigen::MatrixXd H_reduced = H.topLeftCorner(k, k);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> small_solver(H_reduced);
        
        EigenResults results;
        results.converged = (small_solver.info() == Eigen::Success);
        
        if (results.converged) {
            // Transform back to original eigenvalues
            Eigen::VectorXd theta = small_solver.eigenvalues();
            results.eigenvalues.resize(k);
            
            for (int i = 0; i < k; ++i) {
                results.eigenvalues(i) = spectral_shift + 1.0 / theta(i);
            }
            
            // Compute eigenvectors
            Eigen::MatrixXd Y = small_solver.eigenvectors();
            results.eigenvectors = V.leftCols(k) * Y;
            
            // Sort eigenvalues and eigenvectors
            sortEigenResults(results);
            
            // Keep only requested number
            if (num_eigenvalues < k) {
                results.eigenvalues = results.eigenvalues.head(num_eigenvalues);
                results.eigenvectors = results.eigenvectors.leftCols(num_eigenvalues);
            }
            
            computeResiduals(results);
        }
        
        results.iterations = k;
        return results;
    }
    
    void computeResiduals(EigenResults& results) {
        int num_eigs = results.eigenvalues.size();
        results.residual_norms.resize(num_eigs);
        
        for (int i = 0; i < num_eigs; ++i) {
            Eigen::VectorXd v = results.eigenvectors.col(i);
            Eigen::VectorXd residual = system_matrix * v - results.eigenvalues(i) * mass_matrix * v;
            results.residual_norms[i] = residual.norm() / v.norm();
        }
    }
    
    void sortEigenResults(EigenResults& results) {
        int num_eigs = results.eigenvalues.size();
        std::vector<int> indices(num_eigs);
        std::iota(indices.begin(), indices.end(), 0);
        
        // Sort by eigenvalue magnitude
        std::sort(indices.begin(), indices.end(), 
                 [&](int i, int j) { return results.eigenvalues(i) < results.eigenvalues(j); });
        
        // Reorder eigenvalues and eigenvectors
        Eigen::VectorXd sorted_eigenvalues(num_eigs);
        Eigen::MatrixXd sorted_eigenvectors(results.eigenvectors.rows(), num_eigs);
        
        for (int i = 0; i < num_eigs; ++i) {
            sorted_eigenvalues(i) = results.eigenvalues(indices[i]);
            sorted_eigenvectors.col(i) = results.eigenvectors.col(indices[i]);
        }
        
        results.eigenvalues = sorted_eigenvalues;
        results.eigenvectors = sorted_eigenvectors;
    }
    
    bool verifyOrthogonality(const EigenResults& results) {
        int num_eigs = results.eigenvalues.size();
        
        for (int i = 0; i < num_eigs; ++i) {
            for (int j = i + 1; j < num_eigs; ++j) {
                Eigen::VectorXd vi = results.eigenvectors.col(i);
                Eigen::VectorXd vj = results.eigenvectors.col(j);
                
                double orthogonality_error = std::abs(vi.transpose() * mass_matrix * vj);
                
                if (orthogonality_error > tolerance * 100) {
                    return false;
                }
            }
        }
        return true;
    }
    
    double estimateSpectralShift() {
        // Estimate shift parameter for better convergence
        // Use Gershgorin circle theorem for rough estimate
        
        Eigen::VectorXd diagonal = system_matrix.diagonal();
        Eigen::VectorXd row_sums(system_matrix.rows());
        
        for (int i = 0; i < system_matrix.outerSize(); ++i) {
            double sum = 0.0;
            for (Eigen::SparseMatrix<double>::InnerIterator it(system_matrix, i); it; ++it) {
                if (it.row() != it.col()) {
                    sum += std::abs(it.value());
                }
            }
            row_sums(i) = sum;
        }
        
        double min_eigenvalue_estimate = (diagonal.array() - row_sums.array()).minCoeff();
        
        return min_eigenvalue_estimate - 1.0; // Shift slightly below estimated minimum
    }
    
    void setupPreconditioner() {
        preconditioner.setDroptol(0.01);
        preconditioner.setFillfactor(10);
        preconditioner.compute(system_matrix);
        
        if (preconditioner.info() != Eigen::Success) {
            use_preconditioning = false;
            std::cerr << "Warning: Preconditioner setup failed, proceeding without preconditioning" << std::endl;
        }
    }
    
    void setTolerances(double solver_tol, int max_iter) {
        tolerance = solver_tol;
        max_iterations = max_iter;
    }
    
    void enablePreconditioning(bool enable) {
        use_preconditioning = enable;
    }
};

// Fallback solvers for robustness
class LOBPCGSolver {
private:
    double tolerance;
    int max_iterations;
    int block_size;
    
public:
    LOBPCGSolver(double tol = 1e-12, int max_iter = 1000, int bs = 3) 
        : tolerance(tol), max_iterations(max_iter), block_size(bs) {}
    
    EigenResults solve(const Eigen::SparseMatrix<double>& A, 
                      const Eigen::SparseMatrix<double>& M, 
                      int num_eigenvalues) {
        
        int n = A.rows();
        int block_size_actual = std::min(block_size, num_eigenvalues);
        
        // Initialize random block of vectors
        Eigen::MatrixXd X = Eigen::MatrixXd::Random(n, block_size_actual);
        
        // Mass-orthogonalize initial block
        massOrthogonalize(X, M);
        
        Eigen::MatrixXd P = Eigen::MatrixXd::Zero(n, block_size_actual); // Search directions
        Eigen::MatrixXd R = Eigen::MatrixXd::Zero(n, block_size_actual); // Residuals
        
        EigenResults results;
        results.eigenvalues.resize(num_eigenvalues);
        results.eigenvectors.resize(n, num_eigenvalues);
        
        for (int iter = 0; iter < max_iterations; ++iter) {
            // Compute residuals: R = AX - MX*Λ
            Eigen::MatrixXd AX = A * X;
            Eigen::MatrixXd MX = M * X;
            Eigen::MatrixXd XtMX = X.transpose() * MX;
            Eigen::MatrixXd XtAX = X.transpose() * AX;
            
            // Solve small eigenvalue problem
            Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> small_solver;
            small_solver.compute(XtAX, XtMX);
            
            Eigen::VectorXd lambda = small_solver.eigenvalues();
            Eigen::MatrixXd eigvecs = small_solver.eigenvectors();
            
            // Update X
            X = X * eigvecs;
            MX = MX * eigvecs;
            AX = AX * eigvecs;
            
            // Compute residuals
            for (int i = 0; i < block_size_actual; ++i) {
                R.col(i) = AX.col(i) - lambda(i) * MX.col(i);
            }
            
            // Check convergence
            double max_residual = 0.0;
            for (int i = 0; i < block_size_actual; ++i) {
                double res_norm = R.col(i).norm() / X.col(i).norm();
                max_residual = std::max(max_residual, res_norm);
            }
            
            if (max_residual < tolerance) {
                // Converged
                results.converged = true;
                results.iterations = iter;
                
                for (int i = 0; i < std::min(num_eigenvalues, block_size_actual); ++i) {
                    results.eigenvalues(i) = lambda(i);
                    results.eigenvectors.col(i) = X.col(i);
                }
                
                break;
            }
            
            // Precondition residuals (simplified - could use actual preconditioner)
            Eigen::MatrixXd W = R; // In practice, apply preconditioner: W = K^{-1} R
            
            // Orthogonalize W against X and P
            massOrthogonalizeAgainst(W, X, M);
            if (iter > 0) {
                massOrthogonalizeAgainst(W, P, M);
            }
            
            // Update search directions
            if (iter > 0) {
                // Rayleigh-Ritz on [X, W, P]
                int total_size = 3 * block_size_actual;
                Eigen::MatrixXd combined(n, total_size);
                combined.leftCols(block_size_actual) = X;
                combined.middleCols(block_size_actual, block_size_actual) = W;
                combined.rightCols(block_size_actual) = P;
                
                // Mass-orthogonalize combined block
                massOrthogonalize(combined, M);
                
                // Rayleigh-Ritz
                Eigen::MatrixXd A_combined = A * combined;
                Eigen::MatrixXd M_combined = M * combined;
                Eigen::MatrixXd H = combined.transpose() * A_combined;
                Eigen::MatrixXd S = combined.transpose() * M_combined;
                
                Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> rr_solver;
                rr_solver.compute(H, S);
                
                // Extract best approximations
                Eigen::MatrixXd ritz_vectors = rr_solver.eigenvectors().leftCols(block_size_actual);
                X = combined * ritz_vectors;
                
                // Update P for next iteration
                P = combined * rr_solver.eigenvectors().middleCols(block_size_actual, block_size_actual);
            } else {
                P = W;
            }
        }
        
        if (!results.converged) {
            results.iterations = max_iterations;
        }
        
        return results;
    }
    
private:
    void massOrthogonalize(Eigen::MatrixXd& X, const Eigen::SparseMatrix<double>& M) {
        int k = X.cols();
        
        for (int i = 0; i < k; ++i) {
            // Orthogonalize against previous vectors
            for (int j = 0; j < i; ++j) {
                double proj = X.col(j).transpose() * M * X.col(i);
                X.col(i) -= proj * X.col(j);
            }
            
            // Normalize
            double norm = std::sqrt(X.col(i).transpose() * M * X.col(i));
            X.col(i) /= norm;
        }
    }
    
    void massOrthogonalizeAgainst(Eigen::MatrixXd& W, const Eigen::MatrixXd& X, 
                                 const Eigen::SparseMatrix<double>& M) {
        int k_w = W.cols();
        int k_x = X.cols();
        
        for (int i = 0; i < k_w; ++i) {
            for (int j = 0; j < k_x; ++j) {
                double proj = X.col(j).transpose() * M * W.col(i);
                W.col(i) -= proj * X.col(j);
            }
        }
    }
};

class PhysicalQuantityExtractor {
private:
    double hbar_coherence;
    double coherence_light_speed;
    
public:
    PhysicalQuantityExtractor(double hbar_phi = 1.0e-34, double c_phi = 3.0e8) 
        : hbar_coherence(hbar_phi), coherence_light_speed(c_phi) {}
    
    std::vector<double> computeMassSpectrum(const EigenResults& eigenvalue_results) {
        std::vector<double> masses;
        
        for (int i = 0; i < eigenvalue_results.eigenvalues.size(); ++i) {
            double lambda_k = eigenvalue_results.eigenvalues(i);
            double mass_k = hbar_coherence * lambda_k / (2.0 * M_PI);
            masses.push_back(mass_k);
        }
        
        return masses;
    }
    
    std::vector<double> computeRecurrenceTimes(const EigenResults& eigenvalue_results) {
        std::vector<double> recurrence_times;
        
        for (int i = 0; i < eigenvalue_results.eigenvalues.size(); ++i) {
            double lambda_k = eigenvalue_results.eigenvalues(i);
            double tau_rec = 2.0 * M_PI / lambda_k;
            recurrence_times.push_back(tau_rec);
        }
        
        return recurrence_times;
    }
    
    double estimateBECCoherenceForce(const CoherenceField& field, 
                                   const CoherenceField::Point& location,
                                   double atom_mass = 1.4e-25) { // Rb-87 mass in kg
        
        Eigen::Vector3d coherence_gradient = field.gradientCoherenceDensity(location);
        double force_magnitude = (hbar_coherence / atom_mass) * coherence_gradient.norm();
        
        return force_magnitude;
    }
    
    double computeGravitationalWaveModification(double lambda_c, double kappa_c) {
        return std::sqrt(lambda_c / kappa_c);
    }
    
    int computeCMBCoherencePeakPosition(double lambda_c, double kappa_c, double redshift = 1100.0) {
        double base_peak = 200.0 * std::sqrt(lambda_c / kappa_c);
        return static_cast<int>(base_peak * std::sqrt(1.0 + redshift));
    }
    
    void generateExperimentalPredictions(const EigenResults& results,
                                       const CoherenceField& field,
                                       double lambda_c = 1e-4,
                                       double kappa_c = 1e-6) {
        
        auto masses = computeMassSpectrum(results);
        auto recurrence_times = computeRecurrenceTimes(results);
        
        std::cout << "=== LCFT Experimental Predictions ===" << std::endl;
        std::cout << "Fundamental mass scale: " << masses[0] << " kg" << std::endl;
        std::cout << "Fundamental recurrence time: " << recurrence_times[0] << " s" << std::endl;
        
        // BEC predictions
        CoherenceField::Point bec_location(0.0, 0.0, 0.0);
        double bec_force = estimateBECCoherenceForce(field, bec_location);
        std::cout << "Expected BEC coherence force: " << bec_force << " N" << std::endl;
        
        // Gravitational wave predictions
        double gw_modification = computeGravitationalWaveModification(lambda_c, kappa_c);
        std::cout << "Gravitational wave polarization ratio h_phi/h_+: " << gw_modification << std::endl;
        
        // CMB predictions
        int cmb_peak = computeCMBCoherencePeakPosition(lambda_c, kappa_c);
        std::cout << "Expected CMB coherence peak at l = " << cmb_peak << std::endl;
        
        std::cout << "===============================" << std::endl;
    }
};

