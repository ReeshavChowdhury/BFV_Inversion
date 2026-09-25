// norm_tower.cpp -- Galois Norm Tower, depth L_F+s+1 = 22.
// Correction product tree runs concurrently with the scalar inversion.
#include "bfv_inversion.h"
#include <iostream>
#include <thread>
#include <chrono>
#include <exception>
#include <mutex>

namespace bfv_inv {

seal::Ciphertext
norm_tower_invert(const seal::Ciphertext&  ct_m,
                  const seal::SEALContext& context,
                  const seal::RelinKeys&   rk,
                  const seal::GaloisKeys&  gk,
                  TimingBreakdown*         stats)
{
    using Clock = std::chrono::high_resolution_clock;
    auto ms = [](auto t0, auto t1) {
        return std::chrono::duration<double, std::milli>(t1 - t0).count();
    };

    if (verbose())
        std::cout << "=== Norm Tower (target depth " << nt_depth() << ") ===\n";

    auto t_all = Clock::now();

    auto t0 = Clock::now();
    auto T = frobenius_images(ct_m, context, gk);
    auto t1 = Clock::now();
    if (verbose()) {
        std::cout << "  T_0..T_" << D - 1 << " generated (depth 0).\n";
    }

    auto t_norm_s = Clock::now();
    // Correction uses T_1..T_{d-1}; the full norm tree consumes a
    // moved copy of T so we do not keep two full d-vectors live.
    std::vector<seal::Ciphertext> corr(T.begin() + 1, T.end());
    seal::Ciphertext nm = product_tree_parallel(T, context, rk);
    if (verbose())
        std::cout << "  Nm computed (depth s=" << S << ").\n";

    seal::Ciphertext nm_inv, correction;
    std::exception_ptr eptr;
    std::mutex eptr_mu;
    std::thread th_inv([&]() {
        try {
            nm_inv = binary_exp(nm, P - 2, context, rk);
        } catch (...) {
            std::lock_guard<std::mutex> g(eptr_mu);
            if (!eptr) eptr = std::current_exception();
        }
    });
    std::thread th_corr([&]() {
        try {
            correction = product_tree_parallel(corr, context, rk);
        } catch (...) {
            std::lock_guard<std::mutex> g(eptr_mu);
            if (!eptr) eptr = std::current_exception();
        }
    });
    th_inv.join();
    th_corr.join();
    if (eptr) std::rethrow_exception(eptr);
    auto t_par_e = Clock::now();
    if (verbose()) {
        std::cout << "  Nm^{-1} + correction (overlapped): depth s+L_F="
                  << S + LF << ".\n";
    }

    auto t_fin_s = Clock::now();
    seal::Evaluator ev(context);
    // Align + multiply + drop, same helper path as the product tree:
    // reuse a one-pair tree so modulus switching stays consistent.
    std::vector<seal::Ciphertext> last{std::move(correction), std::move(nm_inv)};
    seal::Ciphertext result = product_tree(last, ev, rk);
    auto t_fin_e = Clock::now();
    if (verbose())
        std::cout << "  Final multiply (depth=" << nt_depth() << ").\n";

    auto t_end = Clock::now();
    if (stats) {
        stats->frobenius_ms = ms(t0, t1);
        stats->phase1_ms    = ms(t_norm_s, t_par_e); // tree + invert
        stats->phase2_ms    = ms(t_fin_s, t_fin_e);
        stats->total_ms     = ms(t_all, t_end);
        // Critical path: s (norm tree) + L_F (exponent) + 1 final multiply.
        stats->mul_ops = S + LF + 1;
    }
    return result;
}

} // namespace bfv_inv
