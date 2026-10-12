// Copyright @ Chun Shen 2020

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include "MomentumSamplerBase.h"

MomentumSamplerBase::MomentumSamplerBase() : tol_(1e-8), eps_(1e-16),
    CDFc_m_tilde_(std::numeric_limits<double>::quiet_NaN()),
    CDFc_mu_tilde_(std::numeric_limits<double>::quiet_NaN()),
    CDFc_m_minus_mu_tilde_(std::numeric_limits<double>::quiet_NaN()) {}


void MomentumSamplerBase::update_cache(const double m_minus_mu_tilde) {
    m_minus_mu_tilde_ = m_minus_mu_tilde;
    CDF_0_mtilde_ = CDF_0(m_minus_mu_tilde);
    CDF_1_mtilde_ = CDF_1(m_minus_mu_tilde);
    CDF_2_mtilde_ = CDF_2(m_minus_mu_tilde);
}


void MomentumSamplerBase::update_combined_CDF(
        const double m_tilde, const double mu_tilde, const double m_term) {
// inverse_CDF bisects the combined CDF, which it evaluated from the three
// tables at every step. On an isothermal surface (m, mu, T) stay the same
// for a species, so it is tabulated once and kept while m_tilde, mu_tilde
// and the cached CDF_x(m_minus_mu_tilde) stay exactly the same: the same
// expression gives the same doubles.
    if (m_tilde == CDFc_m_tilde_ && mu_tilde == CDFc_mu_tilde_
        && m_minus_mu_tilde_ == CDFc_m_minus_mu_tilde_) {
        return;
    }
    CDFc_.resize(tb_len_);
    for (int i = 0; i < tb_len_; i++) {
        CDFc_[i] = (CDF_2_[i] + 2.*mu_tilde*CDF_1_[i]
                    + (mu_tilde*mu_tilde - m_tilde*m_tilde/2.)*CDF_0_[i]
                    - m_term);
    }
    CDFc_m_tilde_ = m_tilde;
    CDFc_mu_tilde_ = mu_tilde;
    CDFc_m_minus_mu_tilde_ = m_minus_mu_tilde_;
}


double MomentumSamplerBase::Sample_a_momentum(double m, double T, double mu) {
    T = std::max(eps_, T);
    const double m_tilde = m/T;
    const double mu_tilde = mu/T;
    if (std::abs((m_tilde - mu_tilde) - m_minus_mu_tilde_) > tol_)
        update_cache(m_tilde - mu_tilde);
    double m_term = (CDF_2_mtilde_ + 2.*mu_tilde*CDF_1_mtilde_
                     + (mu_tilde*mu_tilde - m_tilde*m_tilde/2.)*CDF_0_mtilde_);

    update_combined_CDF(m_tilde, mu_tilde, m_term);
    const int idx_max = tb_len_ - 1;
    double CDF_max = CDFc_[idx_max];

    const int idx_min = static_cast<int>(
            (m_tilde - mu_tilde - Etilde_[0])/(Etilde_[1] - Etilde_[0]));
    if (idx_min < 0 || idx_min >= idx_max) {
        std::cout << "[Error: MomentumSampler] out of range, "
                  << "Etilde_min = " << Etilde_[0]
                  << ", Etilde_max = " << Etilde_[idx_max]
                  << ", m = " << m << " GeV, T = " << T << " GeV, mu = "
                  << mu << " GeV. mtilde - mutilde = " << m_tilde - mu_tilde
                  << std::endl;
        exit(1);
    }

    double p_sample = -1.;
    double E_sample = -1.;
    double accept_ratio = 1.;
    do {
        double r = ran_gen_ptr_->rand_uniform()*CDF_max;
        double Etilde_sample = inverse_CDF(r, idx_min, CDF_max,
                                           m_tilde, mu_tilde);
        // convert back to particle momentum
        E_sample = T*Etilde_sample + mu;
        p_sample = sqrt(E_sample*E_sample - m*m);
        accept_ratio = (p_sample/E_sample)/(1. - m*m/(2.*E_sample*E_sample));
    } while (ran_gen_ptr_->rand_uniform() > accept_ratio);
    return(p_sample);
}


double MomentumSamplerBase::inverse_CDF(
                        const double r, const int idx_0,
                        const double CDF_max,
                        const double m_tilde, const double mu_tilde) const {
// bisection of the combined CDF (update_combined_CDF) between idx_0 and the
// end of the table
    int idx_min = idx_0;
    int idx_max = tb_len_ - 1;
    double r_min = CDFc_[idx_min];
    double r_max = CDF_max;
    while ((idx_max - idx_min) > 1) {
        const int idx_mid = static_cast<int>((idx_max + idx_min)/2 + 0.1);
        const double r_mid = CDFc_[idx_mid];
        // the same bisection, as selects: the branch is unpredictable
        const bool go_left = (r < r_mid);
        idx_max = go_left ? idx_mid : idx_max;
        r_max   = go_left ? r_mid : r_max;
        idx_min = go_left ? idx_min : idx_mid;
        r_min   = go_left ? r_min : r_mid;
    }
    double E0 = Etilde_[idx_min];
    if (E0 < m_tilde - mu_tilde) {
        E0 = m_tilde - mu_tilde;
        r_min = 0.;
    }
    double res = (E0 + (Etilde_[idx_max] - E0)
                       /std::max(eps_, (r_max - r_min))*(r - r_min));
    return(res);
}
