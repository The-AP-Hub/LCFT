#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <memory>

// Main LCFT computational framework
class LCFTFramework {
private:
    std::unique_ptr<CoherenceField> coherence_field;
    std::unique_ptr<AdaptiveMesh> mesh;
    std::unique_ptr<SystemMatrices> matrices;
    std::unique_ptr<SparseEigenvalueSolver> primary_solver;
    std::unique_ptr<LOBPCGSolver> backup_solver;
    std::unique_ptr<PhysicalQuantityExtractor> physics_extractor;
    
    // Physical parameters
    double hbar_phi;
    double lambda_c;
    double kappa_c;
    double xi_c;
    
    // Computational parameters
    double mesh_refinement_tolerance;
    int num_eigenvalues_requested;
    double eigenvalue_tolerance;
    
    // Results storage
    EigenResults current_results;
    std::vector<double> mass_spectrum;
    std::vector<double> recurrence_times;
    bool computation_complete;

public:
    LCFTFramework() : computation_complete(false) {
        // Initialize physical parameters with theoretical estimates
        hbar_phi = 1.054571817e-34; // Coherence action quantum
        lambda_c = 1e-4;  // Coherence coupling constant
        kappa_c = 1e-6;   // Correlation strength
        xi_c = 1e-16;     // Coherence correlation length
        
        // Initialize computational parameters
        mesh_refinement_tolerance = 1e-8;
        num_eigenvalues_requested = 50;
        eigenvalue_tolerance = 1e-12;
        
        // Initialize components
        coherence_field = std::make_unique<CoherenceField>(1e-12, 0.5);
        physics_extractor = std::make_unique<PhysicalQuantityExtractor>(hbar_phi, 3e8);
        primary_solver = std::make_unique<SparseEigenvalueSolver>(150, eigenvalue_tolerance, 1000);
        backup_solver = std::make_unique<LOBPCGSolver>(eigenvalue_tolerance, 1000, 5);
    }
    
    void setPhysicalParameters(double hbar, double lambda, double kappa, double correlation_length) {
        hbar_phi = hbar;
        lambda_c = lambda;
        kappa_c = kappa;
        xi_c = correlation_length;
        
        physics_extractor = std::make_unique<PhysicalQuantityExtractor>(hbar_phi, 3e8);
    }
    
    void setComputationalParameters(double mesh_tol, int num_eigs, double eig_tol) {
        mesh_refinement_tolerance = mesh_tol;
        num_eigenvalues_requested = num_eigs;
        eigenvalue_tolerance = eig_tol;
        
        primary_solver->setTolerances(eigenvalue_tolerance, 1000);
    }
    
    bool initializeCoherenceField(const std::vector<CoherenceField::Point>& domain_bounds,
                                 std::function<double(const CoherenceField::Point&)> initial_condition) {
        
        std::cout << "Initializing coherence field..." << std::endl;
        
        try {
            // Create initial mesh
            mesh = std::make_unique<AdaptiveMesh>(coherence_field.get(), 0.1, 1.0);
            mesh->generateMesh(domain_bounds);
            
            // Initialize coherence field on mesh
            coherence_field->initializeField(mesh->getNodes(), initial_condition);
            
            // Refine mesh near coherence boundaries
            std::cout << "Refining mesh near coherence boundaries..." << std::endl;
            mesh->refineNearBoundaries();
            
            std::cout << "Coherence field initialized with " << mesh->getNodes().size() 
                      << " nodes and " << mesh->getElements().size() << " elements." << std::endl;
            
            return true;
            
        } catch (const std::exception& e) {
            std::cerr << "Error initializing coherence field: " << e.what() << std::endl;
            return false;
        }
    }
    
    bool assembleSystemMatrices() {
        std::cout << "Assembling system matrices..." << std::endl;
        
        try {
            size_t num_nodes = mesh->getNodes().size();
            matrices = std::make_unique<SystemMatrices>(num_nodes);
            
            // Assemble stiffness matrix
            std::cout << "  Assembling stiffness matrix..." << std::endl;
            matrices->assembleStiffnessMatrix(*coherence_field, *mesh);
            
            // Assemble mass matrix
            std::cout << "  Assembling mass matrix..." << std::endl;
            matrices->assembleMassMatrix(*coherence_field, *mesh);
            
            // Apply boundary conditions
            std::cout << "  Applying boundary conditions..." << std::endl;
            matrices->applyBoundaryConditions(*coherence_field, *mesh);
            
            // Verify matrix properties
            if (matrices->verifyMatrixProperties()) {
                std::cout << "Matrix assembly successful - symmetry and positive definiteness verified." << std::endl;
                return true;
            } else {
                std::cerr << "Warning: Matrix verification failed." << std::endl;
                return false;
            }
            
        } catch (const std::exception& e) {
            std::cerr << "Error assembling matrices: " << e.what() << std::endl;
            return false;
        }
    }
    
    bool solveEigenvalueProblem() {
        std::cout << "Solving eigenvalue problem..." << std::endl;
        std::cout << "Requesting " << num_eigenvalues_requested << " eigenvalues with tolerance " 
                  << eigenvalue_tolerance << std::endl;
        
        try {
            // Set up matrices for solver
            auto& stiffness = matrices->getStiffnessMatrix();
            auto& mass = matrices->getMassMatrix();
            auto& boundary = matrices->getBoundaryMatrix();
            
            // Combined system matrix
            Eigen::SparseMatrix<double> system_matrix = stiffness + boundary;
            
            primary_solver->setMatrices(system_matrix, mass);
            
            // Estimate optimal spectral shift
            double shift = primary_solver->estimateSpectralShift();
            primary_solver->setSpectralShift(shift);
            std::cout << "Using spectral shift: " << shift << std::endl;
            
            // Attempt primary solution
            current_results = primary_solver->solveGeneralizedEigenvalue(num_eigenvalues_requested);
            
            if (!current_results.converged) {
                std::cout << "Primary solver failed, attempting backup LOBPCG solver..." << std::endl;
                current_results = backup_solver->solve(system_matrix, mass, num_eigenvalues_requested);
            }
            
            if (current_results.converged) {
                std::cout << "Eigenvalue problem solved successfully!" << std::endl;
                std::cout << "Iterations: " << current_results.iterations << std::endl;
                std::cout << "Computation time: " << current_results.computation_time << " seconds" << std::endl;
                
                // Verify solution quality
                bool orthogonal = primary_solver->verifyOrthogonality(current_results);
                std::cout << "Eigenvector orthogonality check: " << (orthogonal ? "PASSED" : "FAILED") << std::endl;
                
                // Display residual norms
                std::cout << "Residual norms for first 10 eigenvalues:" << std::endl;
                for (int i = 0; i < std::min(10, static_cast<int>(current_results.residual_norms.size())); ++i) {
                    std::cout << "  \u03bb_" << i << ": " << std::scientific << current_results.residual_norms[i] << std::endl;
                }
                
                return true;
            } else {
                std::cerr << "Both primary and backup solvers failed to converge." << std::endl;
                return false;
            }
            
        } catch (const std::exception& e) {
            std::cerr << "Error solving eigenvalue problem: " << e.what() << std::endl;
            return false;
        }
    }
    
    void extractPhysicalQuantities() {
        if (!current_results.converged) {
            std::cerr << "Cannot extract physical quantities - eigenvalue problem not solved." << std::endl;
            return;
        }
        
        std::cout << "\nExtracting physical quantities..." << std::endl;
        
        // Compute mass spectrum
        mass_spectrum = physics_extractor->computeMassSpectrum(current_results);
        
        // Compute recurrence times
        recurrence_times = physics_extractor->computeRecurrenceTimes(current_results);
        
        // Generate experimental predictions
        physics_extractor->generateExperimentalPredictions(current_results, *coherence_field, lambda_c, kappa_c);
        
        computation_complete = true;
    }
    
    void outputResults(const std::string& filename_prefix) {
        if (!computation_complete) {
            std::cerr << "Cannot output results - computation not complete." << std::endl;
            return;
        }
        
        // Output eigenvalues
        std::ofstream eigen_file(filename_prefix + "_eigenvalues.dat");
        eigen_file << std::scientific << std::setprecision(16);
        eigen_file << "# Eigenvalue index, Eigenvalue, Residual norm" << std::endl;
        
        for (int i = 0; i < current_results.eigenvalues.size(); ++i) {
            eigen_file << i << "\t" << current_results.eigenvalues(i) << "\t" 
                      << current_results.residual_norms[i] << std::endl;
        }
        eigen_file.close();
        
        // Output mass spectrum
        std::ofstream mass_file(filename_prefix + "_mass_spectrum.dat");
        mass_file << std::scientific << std::setprecision(16);
        mass_file << "# Mass index, Mass (kg), Recurrence time (s)" << std::endl;
        
        for (size_t i = 0; i < mass_spectrum.size(); ++i) {
            mass_file << i << "\t" << mass_spectrum[i] << "\t" << recurrence_times[i] << std::endl;
        }
        mass_file.close();
        
        // Output eigenvectors (first few)
        std::ofstream eigvec_file(filename_prefix + "_eigenvectors.dat");
        eigvec_file << std::scientific << std::setprecision(16);
        eigvec_file << "# Node coordinates and first 5 eigenvector components" << std::endl;
        
        auto& nodes = mesh->getNodes();
        int num_vecs_output = std::min(5, static_cast<int>(current_results.eigenvectors.cols()));
        
        for (size_t i = 0; i < nodes.size(); ++i) {
            eigvec_file << nodes[i].x << "\t" << nodes[i].y << "\t" << nodes[i].z;
            for (int j = 0; j < num_vecs_output; ++j) {
                eigvec_file << "\t" << current_results.eigenvectors(i, j);
            }
            eigvec_file << std::endl;
        }
        eigvec_file.close();
        
        // Output computational summary
        std::ofstream summary_file(filename_prefix + "_summary.txt");
        summary_file << "LCFT Computational Results Summary" << std::endl;
        summary_file << "=================================" << std::endl;
        summary_file << "Physical Parameters:" << std::endl;
        summary_file << "  hbar_phi = " << hbar_phi << std::endl;
        summary_file << "  lambda_c = " << lambda_c << std::endl;
        summary_file << "  kappa_c = " << kappa_c << std::endl;
        summary_file << "  xi_c = " << xi_c << std::endl;
        summary_file << std::endl;
        
        summary_file << "Computational Parameters:" << std::endl;
        summary_file << "  Mesh nodes: " << mesh->getNodes().size() << std::endl;
        summary_file << "  Mesh elements: " << mesh->getElements().size() << std::endl;
        summary_file << "  Eigenvalues computed: " << current_results.eigenvalues.size() << std::endl;
        summary_file << "  Solver iterations: " << current_results.iterations << std::endl;
        summary_file << "  Computation time: " << current_results.computation_time << " seconds" << std::endl;
        summary_file << std::endl;
        
        summary_file << "Key Results:" << std::endl;
        summary_file << "  Fundamental eigenvalue: " << current_results.eigenvalues(0) << std::endl;
        summary_file << "  Fundamental mass: " << mass_spectrum[0] << " kg" << std::endl;
        summary_file << "  Fundamental recurrence time: " << recurrence_times[0] << " s" << std::endl;
        
        double gw_ratio = physics_extractor->computeGravitationalWaveModification(lambda_c, kappa_c);
        summary_file << "  GW polarization ratio: " << gw_ratio << std::endl;
        
        int cmb_peak = physics_extractor->computeCMBCoherencePeakPosition(lambda_c, kappa_c);
        summary_file << "  CMB coherence peak: l = " << cmb_peak << std::endl;
        
        summary_file.close();
        
        std::cout << "Results written to files with prefix: " << filename_prefix << std::endl;
    }
    
    void runConvergenceStudy() {
        std::cout << "\nRunning mesh convergence study..." << std::endl;
        
        std::vector<double> mesh_sizes = {0.2, 0.1, 0.05, 0.025};
        std::vector<double> fundamental_eigenvalues;
        
        for (double h : mesh_sizes) {
            std::cout << "Testing mesh size h = " << h << std::endl;
            
            // Reinitialize with new mesh size
            mesh = std::make_unique<AdaptiveMesh>(coherence_field.get(), h, 1.0);
            
            // Simple domain for convergence test
            std::vector<CoherenceField::Point> domain = {{-1, -1, -1}, {1, 1, 1}};
            
            // Gaussian coherence field for testing
            auto gaussian_field = [](const CoherenceField::Point& p) {
                double r2 = p.x*p.x + p.y*p.y + p.z*p.z;
                return std::exp(-r2);
            };
            
            if (initializeCoherenceField(domain, gaussian_field) &&
                assembleSystemMatrices() &&
                solveEigenvalueProblem()) {
                
                fundamental_eigenvalues.push_back(current_results.eigenvalues(0));
                std::cout << "  Fundamental eigenvalue: " << current_results.eigenvalues(0) << std::endl;
            }
        }
        
        // Analyze convergence
        if (fundamental_eigenvalues.size() >= 3) {
            std::cout << "\nConvergence analysis:" << std::endl;
            for (size_t i = 1; i < fundamental_eigenvalues.size(); ++i) {
                double error = std::abs(fundamental_eigenvalues[i] - fundamental_eigenvalues.back());
                double ratio = error / (mesh_sizes[i] * mesh_sizes[i]); // Expecting h^2 convergence
                std::cout << "  h = " << mesh_sizes[i] << ", error = " << error 
                          << ", error/h^2 = " << ratio << std::endl;
            }
        }
    }
    
    // Main execution function
    bool executeFullComputation(const std::string& output_prefix = "lcft_results") {
        auto total_start = std::chrono::high_resolution_clock::now();
        
        std::cout << "===== LCFT Framework Execution =====" << std::endl;
        std::cout << "Coherence action quantum: " << hbar_phi << " J\u22c5s" << std::endl;
        std::cout << "Coherence coupling: " << lambda_c << std::endl;
        std::cout << "Correlation strength: " << kappa_c << std::endl;
        std::cout << "Correlation length: " << xi_c << " m" << std::endl;
        std::cout << "=====================================" << std::endl;
        
        // Define computational domain
        std::vector<CoherenceField::Point> domain_bounds = {{-2, -2, -2}, {2, 2, 2}};
        
        // Define coherence field - example with localized coherence structure
        auto coherence_profile = [this](const CoherenceField::Point& p) {
            double r = std::sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
            double coherence_scale = this->xi_c * 1e15; // Scale up for computational convenience
            return 0.8 * std::exp(-r*r / (coherence_scale*coherence_scale)) + 0.2;
        };
        
        // Execute computation pipeline
        bool success = true;
        success &= initializeCoherenceField(domain_bounds, coherence_profile);
        success &= assembleSystemMatrices();
        success &= solveEigenvalueProblem();
        
        if (success) {
            extractPhysicalQuantities();
            outputResults(output_prefix);
            
            auto total_end = std::chrono::high_resolution_clock::now();
            double total_time = std::chrono::duration<double>(total_end - total_start).count();
            
            std::cout << "\n===== Computation Complete =====" << std::endl;
            std::cout << "Total execution time: " << total_time << " seconds" << std::endl;
            std::cout << "=================================" << std::endl;
            
            return true;
        } else {
            std::cerr << "Computation failed at some stage." << std::endl;
            return false;
        }
    }
};

// Example usage and testing
int main() {
    try {
        LCFTFramework framework;
        
        // Set physical parameters based on theoretical estimates
        framework.setPhysicalParameters(
            1.054571817e-34,  // hbar_phi
            1e-4,             // lambda_c
            1e-6,             // kappa_c
            1e-16             // xi_c
        );
        
        // Set computational parameters for high precision
        framework.setComputationalParameters(
            1e-8,   // mesh refinement tolerance
            50,     // number of eigenvalues
            1e-12   // eigenvalue solver tolerance
        );
        
        // Execute full computation
        bool success = framework.executeFullComputation("lcft_baseline");
        
        if (success) {
            std::cout << "\nPerforming convergence study..." << std::endl;
            framework.runConvergenceStudy();
        }
        
        return success ? 0 : 1;
        
    } catch (const std::exception& e) {
        std::cerr << "Fatal error: " << e.what() << std::endl;
        return 1;
    }
}

