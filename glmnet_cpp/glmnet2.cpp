#include <iostream>
#include <algorithm>
#include <numeric>

#include "dpj_utils.hpp"
#include "glmnet.h"

GLMnet::GLMnet() 
    :       standardize(true),
            lambda_min_ratio(0),
            alpha(1.0),
            converge_thresh(0.0000001),
            max_final_features(0),
            max_path_features(0),
            num_lambdas(100),
            cov_updating(2),
            max_iterations(100000) { }

const std::vector<int>& GLMnet::flags() const { return m_flags; }
const std::vector<double>& GLMnet::coeffs() const { return m_coeffs; }
const std::vector<int>& GLMnet::coeff_counts() const { return m_coeff_counts; }
const std::vector<int>& GLMnet::coeff_ptrs() const { return m_coeff_ptrs; }
const std::vector<double>& GLMnet::intercepts() const { return m_intercepts; }
const std::vector<double>& GLMnet::r_squareds() const { return m_r_squareds; }
const std::vector<double>& GLMnet::lambdas_used() const {return m_lambdas_used;}

// Turns an indicator vector into a compressed vector for the fortran
// routine.
void GLMnet::flags(const std::vector<int>& flags)
{
    std::vector<int> tmp;
    for (unsigned i = 0; i < flags.size(); ++i)
        if (flags[i] > 0) tmp.push_back(i + 1);
    m_flags.resize(tmp.size() + 1);
    std::copy(tmp.begin(), tmp.end(), m_flags.data() + 1);
    m_flags.front() = tmp.size();
    std::cerr << tmp.size() << " flagged probes excluded.\n";
}

// To be called directly before fortran routines.
void GLMnet::set_defaults(int rows, int cols)
{
    // Set these if unset.
    if (m_flags.size() == 0) m_flags.push_back(0);
    if (penalties.size() == 0) penalties = std::vector<double>(cols, 1);

    // TODO: make the following configurable

    if (max_final_features == 0) max_final_features = cols + 1;
    max_path_features = std::min(max_final_features * 2, cols);
    lambda_min_ratio = rows < cols ? 0.01 : 0.0001;    
    m_intercepts.resize(num_lambdas);
    m_coeffs.resize(max_path_features * num_lambdas);
    m_coeff_ptrs.resize(max_path_features);
    m_coeff_counts.resize(num_lambdas);
    m_r_squareds.resize(num_lambdas);
    m_lambdas_used.resize(num_lambdas);
}

void GLMnet::process_outputs()
{
    // We resize everything to num_lambdas_used rather than num_lambdas.
    m_coeffs.resize(num_fits * max_path_features);
    m_coeff_counts.resize(num_fits);
    m_intercepts.resize(num_fits);
    m_r_squareds.resize(num_fits);
    m_lambdas_used.resize(num_fits);

    // hack: fix.lam from R package "lam[1]=exp(2*llam[2]-llam[3])"
    m_lambdas_used[0] = ::exp(2*::log(m_lambdas_used[1]) -
        ::log(m_lambdas_used[2]));
}

// TODO
// Am disregarding the warnings on Eigen::SparseMatrix API stability here. And
// then some.
void GLMnet::spelnet(
                     std::vector<double>& y,
                     std::vector<double>& w,
                     const std::vector<double>& xx,
                     const std::vector<int>& xi,
                     const std::vector<int>& xp
                    )
{
    set_defaults(y.size(), xp.size() - 1);
    std::cerr << "calling spelnet...\n";
    spelnet_(
             cov_updating,
             alpha,
             (int)y.size(),
             (int)(xp.size() - 1),
             // Sparse matrix representation. newGLMnet_f.f90 source says row
             // major but all other evidence points to col major. I use col
             // major here. This must be correct or there would be 'major'
             // errors, groan.

             xx.data(),                 // values
             xp.data(),                 // pointers to blocks of col indices
             xi.data(),                 // col indices

             // 
             y.data(),
             w.data(),
             m_flags.data(),
             penalties.data(),
             max_final_features,
             max_path_features,
             num_lambdas,
             lambda_min_ratio,
             lambdas.data(),
             converge_thresh,
             (int)standardize,
             max_iterations,
             // outputs
             num_fits,
             m_intercepts.data(),
             m_coeffs.data(),
             m_coeff_ptrs.data(),
             m_coeff_counts.data(), // The number of non-zeros for each lambda.
             m_r_squareds.data(),
             m_lambdas_used.data(),
             num_passes,
             error_flag
            );
    fprintf(stderr, "spelnet done with error code: %d\n", error_flag);
    process_outputs();    
}