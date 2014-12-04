#ifndef DPJ_GLMNET_H
#define DPJ_GLMNET_H
#include <vector>

class GLMnet
{
  public:
    GLMnet();
    void spelnet(
                 std::vector<double>& y,
                 std::vector<double>& w,
                 const std::vector<double>& xx,
                 const std::vector<int>& xi,
                 const std::vector<int>& xp
                );

    void flags(const std::vector<int>& flgs);
    const std::vector<int>& flags() const;
    const std::vector<double>& coeffs() const;
    const std::vector<int>& coeff_counts() const;
    const std::vector<int>& coeff_ptrs() const;
    const std::vector<double>& intercepts() const;
    const std::vector<double>& r_squareds() const;
    const std::vector<double>& lambdas_used() const;

    bool standardize;
    double lambda_min_ratio;
    double alpha;
    double converge_thresh;

    int max_final_features;
    int max_path_features;
    int num_lambdas;
    int cov_updating;
    int max_iterations;
    int num_passes;    
    int error_flag;    
    int num_fits;

    std::vector<double> penalties;
    std::vector<double> lambdas;

  private:

    std::vector<int> m_flags, m_coeff_ptrs, m_coeff_counts;
    std::vector<double> m_r_squareds, m_lambdas_used, m_intercepts, m_coeffs;

    void set_defaults(int rows, int cols);
    void process_outputs();
};

extern "C" {
void spelnet_(
            const int& cov_updating,
            const double& alpha,
            const int& num_observations,
            const int& num_predictors,

            const double* vals,         // sparse matrix
            const int* outer_pointers,  // sparse matrix
            const int* inner_indices,   // sparse matrix

            double* observations,              // vector, overwritten
            double* weights,                   // vector, overwritten
            const int* flags,                  // vector
            const double* predictor_penalties, // vector
            const int& max_final_features,
            const int& max_path_features,
            const int& num_lambdas,
            const double& lambda_min_ratio,
            const double* user_defined_lambdas, // vector
            const double& convergence_threshold,
            const int& standardize_flag, // standardize observations
            const int& maxit,           // maximum iterations
            // Outputs
            int& num_fits,              // actual number of lambdas used 
            double* intercepts,         // vector of intercepts 
            double* coeffs,             // compressed matrix values
            int* coeffs_ptrs,           // compressed matrix ptrs
            int* coeff_counts,          // number non-zero per sol
            double* r_squareds,         // r squared statistic vector
            double* lambdas_used,       // actual lambda values vector
            int& num_passes,            // total iterations
            int& error_flag             // errors
           );
}
            

#endif // DPJ_GLMNET_H