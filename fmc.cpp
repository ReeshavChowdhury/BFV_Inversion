// fmc.cpp -- Frobenius Monomial Circuit, depth L_F + s.
//
//   m^{-1} = T_0^{p-2} · ∏_{k=1}^{d-1} T_k^{p-1}
//
// T_k^{p-1} = κ_{p^k}(m^{p-1}). Both powers come from one square chain
// when p-1 = (2^k+1) 2^v. The d factors are then multiplied in a tree
// of depth s.
#include "bfv_inversion.h"
#include <iostream>
#include <iomanip>
#include <thread>
#include <vector>
#include <chrono>
#include <exception>
#include <mutex>

namespace bfv_inv {

seal::Ciphertext
fmc_invert(const seal::Ciphertext&  ct_m,
           const seal::SEALContext& context,
           const seal::RelinKeys&   rk,
           const seal::GaloisKeys&  gk,
           TimingBreakdown*         stats)
{
    using Clock = std::chrono::high_resolution_clock;
    auto ms = [](auto t0, auto t1) {
        return std::chrono::duration<double, std::milli>(t1 - t0).count();
    };

    if (verbose()) {
        std::cout << "=== FMC (target depth " << fmc_depth()
                  << ", 2-exp + Frob, d=" << D << ") ===\n";
    }

    auto t_all = Clock::now();

    seal::Ciphertext E0, A;
    std::exception_ptr eptr;
    std::mutex eptr_mu;

    auto t1s = Clock::now();
    bool fused = false;
    try {
        fused = fused_p_exps(ct_m, context, rk, E0, A);
    } catch (...) {
        eptr = std::current_exception();
    }
    if (!fused && !eptr) {
        std::thread th_e0([&]() {
            try {
                E0 = binary_exp(ct_m, P - 2, context, rk);
            } catch (...) {
                std::lock_guard<std::mutex> g(eptr_mu);
                if (!eptr) eptr = std::current_exception();
            }
        });
        std::thread th_a([&]() {
            try {
                A = binary_exp(ct_m, P - 1, context, rk);
            } catch (...) {
                std::lock_guard<std::mutex> g(eptr_mu);
                if (!eptr) eptr = std::current_exception();
            }
        });
        th_a.join();
        th_e0.join();
    }
    if (eptr) std::rethrow_exception(eptr);

    // κ_{p^k}(m^{p-1}) for k = 1..d-1. Depth 0.
    auto t0 = Clock::now();
    std::vector<seal::Ciphertext> frob_tail;
    std::exception_ptr frob_err;
    std::thread th_frob([&]() {
        try {
            if (D > 1) {
                auto F = frobenius_images(A, context, gk);
                frob_tail.reserve(static_cast<size_t>(D - 1));
                for (int k = 1; k < D; ++k)
                    frob_tail.push_back(std::move(F[k]));
            }
        } catch (...) {
            frob_err = std::current_exception();
        }
    });
    th_frob.join();
    if (eptr) std::rethrow_exception(eptr);
    if (frob_err) std::rethrow_exception(frob_err);
    auto t1e = Clock::now();
    if (verbose()) {
        std::cout << "  Phase 1 (m^{p-2} || m^{p-1}, Frob overlapped): "
                  << std::fixed << std::setprecision(0) << ms(t1s, t1e)
                  << " ms\n";
    }

    std::vector<seal::Ciphertext> fac;
    fac.reserve(static_cast<size_t>(D));
    fac.push_back(std::move(E0));
    for (auto& c : frob_tail)
        fac.push_back(std::move(c));
    auto t1 = t1e;
    if (verbose()) {
        std::cout << "  Frob of m^{p-1}: "
                  << std::fixed << std::setprecision(0) << ms(t0, t1)
                  << " ms\n";
    }

    auto t2s = Clock::now();
    seal::Ciphertext result = product_tree_parallel(fac, context, rk);
    auto t2e = Clock::now();
    if (verbose()) {
        std::cout << "  Phase 2 (product tree): "
                  << std::fixed << std::setprecision(0) << ms(t2s, t2e)
                  << " ms\n";
    }

    auto t_end = Clock::now();
    if (stats) {
        stats->frobenius_ms = ms(t0, t1);
        stats->phase1_ms    = ms(t1s, t1e);
        stats->phase2_ms    = ms(t2s, t2e);
        stats->total_ms     = ms(t_all, t_end);
        // Critical path: L_F levels for the exponent, then s tree levels.
        stats->mul_ops      = LF + S;
    }
    return result;
}

} // namespace bfv_inv
