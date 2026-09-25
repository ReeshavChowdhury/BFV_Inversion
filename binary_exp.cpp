// binary_exp.cpp -- addition chains for x^{p-1} and x^{p-2} at depth L_F.
// Odd parts 1, 3, 5, 7, 9 use a fixed lockstep. Any other odd part uses
// a streaming product of the required 2-powers.
#include "bfv_inversion.h"
#include <stdexcept>
#include <thread>
#include <vector>
#include <exception>
#include <mutex>
#include <algorithm>
#include <cstddef>

namespace {

seal::MemoryPoolHandle thread_pool()
{
    thread_local seal::MemoryPoolHandle pool = seal::MemoryPoolHandle::New();
    return pool;
}

} // namespace

namespace bfv_inv {

static void drop_level(seal::Ciphertext& a, seal::Evaluator& ev)
{
    if (a.coeff_modulus_size() > 1)
        ev.mod_switch_to_next_inplace(a, thread_pool());
}

static void align_levels(seal::Ciphertext& a, seal::Ciphertext& b,
                         seal::Evaluator& ev)
{
    if (a.parms_id() == b.parms_id())
        return;
    auto switch_to = [&](seal::Ciphertext& src,
                         const seal::parms_id_type& target) -> bool {
        try {
            ev.mod_switch_to_inplace(src, target);
            return src.parms_id() == target;
        } catch (const std::exception&) {
            return false;
        }
    };
    if (a.coeff_modulus_size() > b.coeff_modulus_size()) {
        if (!switch_to(a, b.parms_id())) switch_to(b, a.parms_id());
    } else if (b.coeff_modulus_size() > a.coeff_modulus_size()) {
        if (!switch_to(b, a.parms_id())) switch_to(a, b.parms_id());
    } else {
        if (!switch_to(a, b.parms_id())) switch_to(b, a.parms_id());
    }
    if (a.parms_id() != b.parms_id())
        throw std::invalid_argument("align_levels failed");
}

static seal::Ciphertext
mul_relin(const seal::Ciphertext& a_in, const seal::Ciphertext& b_in,
          seal::Evaluator& ev, const seal::RelinKeys& rk, bool drop = true)
{
    const seal::Ciphertext *ap = &a_in, *bp = &b_in;
    seal::Ciphertext a_tmp, b_tmp;
    if (a_in.parms_id() != b_in.parms_id()) {
        a_tmp = a_in; b_tmp = b_in;
        align_levels(a_tmp, b_tmp, ev);
        ap = &a_tmp; bp = &b_tmp;
    }
    seal::Ciphertext out;
    auto pool = thread_pool();
    ev.multiply(*ap, *bp, out, pool);
    ev.relinearize_inplace(out, rk, pool);
    if (drop) drop_level(out, ev);
    return out;
}

static void square_relin_inplace(seal::Ciphertext& a, seal::Evaluator& ev,
                                 const seal::RelinKeys& rk, bool drop = true)
{
    auto pool = thread_pool();
    ev.square_inplace(a, pool);
    ev.relinearize_inplace(a, rk, pool);
    if (drop) drop_level(a, ev);
}

static seal::Ciphertext
exp_generic(const seal::Ciphertext& ct, uint64_t exp,
            seal::Evaluator& ev, const seal::RelinKeys& rk)
{
    if (exp == 0) throw std::invalid_argument("exponent must be >= 1");
    if (exp == 1) return ct;
    seal::Ciphertext result = ct;
    int msb = 63;
    while (msb > 0 && !((exp >> msb) & 1)) --msb;
    for (int bit = msb - 1; bit >= 0; --bit) {
        square_relin_inplace(result, ev, rk);
        if ((exp >> bit) & 1)
            result = mul_relin(result, ct, ev, rk);
    }
    return result;
}

// A = x^{2^i} and B = x^{2^i-1} advance together. The square of A does not
// depend on the multiply that updates B, so those two ciphertext
// multiplications run at the same time. Each still drops one level, so the
// SEAL depth of the chain is unchanged.
static seal::Ciphertext
lockstep_pminus2(const seal::Ciphertext& ct, int nA, int nB, bool keep_mid,
                 const seal::SEALContext& context, const seal::RelinKeys& rk)
{
    seal::Ciphertext A = ct, B = ct, C;
    bool have_C = false;
    for (int i = 1; i <= nA; ++i) {
        seal::Ciphertext A2;
        std::exception_ptr square_err;
        std::thread th([&]() {
            try {
                seal::Evaluator ev(context);
                A2 = A;
                square_relin_inplace(A2, ev, rk, true);
            } catch (...) {
                square_err = std::current_exception();
            }
        });
        std::exception_ptr mul_err;
        if (i >= 2 && i <= nB) {
            try {
                seal::Evaluator ev(context);
                B = mul_relin(A, B, ev, rk, true);
            } catch (...) {
                mul_err = std::current_exception();
            }
        }
        th.join();
        if (square_err) std::rethrow_exception(square_err);
        if (mul_err) std::rethrow_exception(mul_err);
        if (keep_mid && i == nB + 1) { C = A2; have_C = true; }
        A = std::move(A2);
    }
    if (nA == 0) return B;
    if (have_C) {
        // odd=7: A=x^{2^{v+2}} (deepest), C=x^{2^{v+1}}, B=x^{2^v-1}.
        // Multiply the two shallow factors first so the last mul is
        // only +1 on A and lands at L_F, not L_F+1.
        seal::Evaluator ev(context);
        seal::Ciphertext bc = mul_relin(B, C, ev, rk, true);
        return mul_relin(A, bc, ev, rk, true);
    }
    seal::Evaluator ev(context);
    return mul_relin(std::move(A), std::move(B), ev, rk, true);
}

static seal::Ciphertext
exp_odd_small(const seal::Ciphertext& ct, int odd,
              seal::Evaluator& ev, const seal::RelinKeys& rk)
{
    if (odd == 1) return ct;
    if (odd == 3) {
        seal::Ciphertext x2 = ct;
        square_relin_inplace(x2, ev, rk);
        return mul_relin(x2, ct, ev, rk);
    }
    if (odd == 5) {
        seal::Ciphertext x2 = ct;
        square_relin_inplace(x2, ev, rk);
        seal::Ciphertext x4 = x2;
        square_relin_inplace(x4, ev, rk);
        return mul_relin(x4, ct, ev, rk);
    }
    if (odd == 7) {
        seal::Ciphertext x2 = ct;
        square_relin_inplace(x2, ev, rk);
        seal::Ciphertext x4 = x2;
        square_relin_inplace(x4, ev, rk);
        seal::Ciphertext x3 = mul_relin(x2, ct, ev, rk);
        return mul_relin(x4, x3, ev, rk);
    }
    if (odd == 9) {
        seal::Ciphertext x2 = ct;
        square_relin_inplace(x2, ev, rk);
        seal::Ciphertext x4 = x2;
        square_relin_inplace(x4, ev, rk);
        seal::Ciphertext x8 = x4;
        square_relin_inplace(x8, ev, rk);
        return mul_relin(x8, ct, ev, rk);
    }
    return seal::Ciphertext{};
}

static seal::Ciphertext
huffman_product(std::vector<seal::Ciphertext> fac,
                seal::Evaluator& ev, const seal::RelinKeys& rk)
{
    if (fac.empty()) throw std::invalid_argument("empty product");
    while (fac.size() > 1) {
        size_t i = 0, j = 1;
        if (fac[j].coeff_modulus_size() > fac[i].coeff_modulus_size())
            std::swap(i, j);
        for (size_t k = 2; k < fac.size(); ++k) {
            auto s = fac[k].coeff_modulus_size();
            if (s > fac[i].coeff_modulus_size()) { j = i; i = k; }
            else if (s > fac[j].coeff_modulus_size()) j = k;
        }
        auto prod = mul_relin(fac[i], fac[j], ev, rk);
        if (i < j) std::swap(i, j);
        fac.erase(fac.begin() + static_cast<std::ptrdiff_t>(i));
        fac.erase(fac.begin() + static_cast<std::ptrdiff_t>(j));
        fac.push_back(std::move(prod));
    }
    return std::move(fac[0]);
}

// Streaming ladder: one running square, keep only the 2-powers we need.
static seal::Ciphertext
exp_pminus2_stream(const seal::Ciphertext& ct,
                   seal::Evaluator& ev, const seal::RelinKeys& rk)
{
    int maxk = 0;
    for (int t = ODD; t > 1; t >>= 1) ++maxk;
    const int maxe = V2 + maxk;
    seal::Ciphertext A = ct, B = ct;
    std::vector<seal::Ciphertext> fac;
    for (int i = 1; i <= maxe; ++i) {
        seal::Ciphertext A2 = A;
        square_relin_inplace(A2, ev, rk, true);
        if (i >= 2 && i <= V2)
            B = mul_relin(A, B, ev, rk, true);
        if (i == V2)
            fac.push_back(B);
        const int k = i - V2;
        if (k >= 1 && (ODD & (1 << k)))
            fac.push_back(A2);
        A = std::move(A2);
    }
    if (fac.empty()) return B;
    return huffman_product(std::move(fac), ev, rk);
}

static seal::Ciphertext
exp_odd_stream(const seal::Ciphertext& ct, int odd,
               seal::Evaluator& ev, const seal::RelinKeys& rk)
{
    if (odd == 1) return ct;
    int maxe = 0;
    for (int t = odd; t > 1; t >>= 1) ++maxe;
    seal::Ciphertext A = ct;
    std::vector<seal::Ciphertext> fac;
    if (odd & 1) fac.push_back(ct);
    for (int i = 1; i <= maxe; ++i) {
        square_relin_inplace(A, ev, rk, true);
        if (odd & (1 << i))
            fac.push_back(A);
    }
    return huffman_product(std::move(fac), ev, rk);
}

static seal::Ciphertext
exp_pminus1(const seal::Ciphertext& ct, seal::Evaluator& ev,
            const seal::RelinKeys& rk)
{
    seal::Ciphertext y;
    if (ODD == 1 || ODD == 3 || ODD == 5 || ODD == 7 || ODD == 9)
        y = exp_odd_small(ct, ODD, ev, rk);
    else
        y = exp_odd_stream(ct, ODD, ev, rk);
    for (int i = 0; i < V2; ++i)
        square_relin_inplace(y, ev, rk);
    return y;
}

static seal::Ciphertext
exp_pminus2(const seal::Ciphertext& ct, seal::Evaluator& ev,
            const seal::SEALContext& context, const seal::RelinKeys& rk)
{
    // Lean lockstep (2–3 live ciphertexts) for the common primes.
    if (ODD == 1) {
        seal::Ciphertext A = ct, B = ct;
        for (int i = 1; i <= V2; ++i) {
            seal::Ciphertext A2 = A;
            square_relin_inplace(A2, ev, rk, true);
            if (i >= 2) B = mul_relin(A, B, ev, rk, true);
            A = std::move(A2);
        }
        return B;
    }
    if (ODD == 3)
        return lockstep_pminus2(ct, V2 + 1, V2, false, context, rk);
    if (ODD == 5)
        return lockstep_pminus2(ct, V2 + 2, V2, false, context, rk);
    if (ODD == 7)
        return lockstep_pminus2(ct, V2 + 2, V2, true, context, rk);
    if (ODD == 9)
        return lockstep_pminus2(ct, V2 + 3, V2, false, context, rk);
    return exp_pminus2_stream(ct, ev, rk);
}

bool fused_p_exps(const seal::Ciphertext& ct,
                  const seal::SEALContext& context,
                  const seal::RelinKeys& rk,
                  seal::Ciphertext& pm2,
                  seal::Ciphertext& pm1)
{
    // p - 1 = odd * 2^{V2}. When odd = 2^k + 1,
    //   p - 1 = 2^{V2+k} + 2^{V2}
    //   p - 2 = 2^{V2+k} + 2^{V2} - 1
    // so both powers are products of values the p-2 lockstep already has:
    //   A_final = x^{2^{V2+k}},  A_mid = x^{2^{V2}},  B_mid = x^{2^{V2}-1}.
    if (ODD < 3) return false;
    const int even = ODD - 1;
    if (even & (even - 1)) return false; // not a power of two
    int k = 0;
    for (int t = even; t > 1; t >>= 1) ++k;

    const int nA = V2 + k;
    const int nB = V2;
    seal::Ciphertext A = ct, B = ct, A_mid;
    for (int i = 1; i <= nA; ++i) {
        seal::Ciphertext A2;
        std::exception_ptr square_err;
        std::thread th([&]() {
            try {
                seal::Evaluator ev(context);
                A2 = A;
                square_relin_inplace(A2, ev, rk, true);
            } catch (...) {
                square_err = std::current_exception();
            }
        });
        std::exception_ptr mul_err;
        if (i >= 2 && i <= nB) {
            try {
                seal::Evaluator ev(context);
                B = mul_relin(A, B, ev, rk, true);
            } catch (...) {
                mul_err = std::current_exception();
            }
        }
        th.join();
        if (square_err) std::rethrow_exception(square_err);
        if (mul_err) std::rethrow_exception(mul_err);
        A = std::move(A2);
        if (i == nB) A_mid = A;
    }

    std::exception_ptr e0, e1;
    std::thread th_pm2([&]() {
        try {
            seal::Evaluator ev(context);
            pm2 = mul_relin(A, B, ev, rk, true);
        } catch (...) { e0 = std::current_exception(); }
    });
    std::thread th_pm1([&]() {
        try {
            seal::Evaluator ev(context);
            pm1 = mul_relin(A, A_mid, ev, rk, true);
        } catch (...) { e1 = std::current_exception(); }
    });
    th_pm2.join();
    th_pm1.join();
    if (e0) std::rethrow_exception(e0);
    if (e1) std::rethrow_exception(e1);
    return true;
}

seal::Ciphertext
binary_exp(const seal::Ciphertext& ct, uint64_t exp,
           const seal::SEALContext& context, const seal::RelinKeys& rk)
{
    seal::Evaluator ev(context);
    if (exp == P - 1) return exp_pminus1(ct, ev, rk);
    if (exp == P - 2) return exp_pminus2(ct, ev, context, rk);
    return exp_generic(ct, exp, ev, rk);
}

seal::Ciphertext
product_tree(std::vector<seal::Ciphertext>& v,
             seal::Evaluator& ev, const seal::RelinKeys& rk)
{
    if (v.empty()) throw std::invalid_argument("empty input");
    while (v.size() > 1) {
        std::vector<seal::Ciphertext> next;
        next.reserve((v.size() + 1) / 2);
        for (size_t i = 0; i + 1 < v.size(); i += 2)
            next.push_back(mul_relin(std::move(v[i]), std::move(v[i + 1]), ev, rk));
        if (v.size() % 2 == 1) next.push_back(std::move(v.back()));
        v = std::move(next);
    }
    return std::move(v[0]);
}

seal::Ciphertext
product_tree_parallel(std::vector<seal::Ciphertext>& v,
                      const seal::SEALContext& context,
                      const seal::RelinKeys& rk)
{
    if (v.empty()) throw std::invalid_argument("empty input");
    while (v.size() > 1) {
        const size_t npair = v.size() / 2;
        const bool oddp = (v.size() % 2) == 1;
        std::vector<seal::Ciphertext> next(npair + (oddp ? 1 : 0));
        std::exception_ptr eptr;
        std::mutex eptr_mu;
        // While a ciphertext still has many RNS limbs, run two pairs at a time.
        // After the exponentiation the ciphertexts are short and every pair
        // in the level runs together.
        const size_t width = (v[0].coeff_modulus_size() >= 16) ? 2 : npair;
        if (oddp) next.back() = std::move(v.back());
        for (size_t base = 0; base < npair; base += width) {
            const size_t end = std::min(npair, base + width);
            std::vector<std::thread> ts;
            ts.reserve(end - base);
            for (size_t i = base; i < end; ++i) {
                ts.emplace_back([&, i]() {
                    try {
                        seal::Evaluator ev(context);
                        next[i] = mul_relin(v[2 * i], v[2 * i + 1], ev, rk);
                    } catch (...) {
                        std::lock_guard<std::mutex> g(eptr_mu);
                        if (!eptr) eptr = std::current_exception();
                    }
                });
            }
            for (auto& t : ts) t.join();
            if (eptr) std::rethrow_exception(eptr);
        }
        v = std::move(next);
    }
    return std::move(v[0]);
}

} // namespace bfv_inv
