// galois_ops.cpp -- Frobenius images via a single κ_p key (depth 0).
#include "bfv_inversion.h"
#include <iomanip>
#include <stdexcept>
#include <iostream>

namespace bfv_inv {

uint64_t frobenius_galois_elt(int k)
{
    const uint64_t mod = 2 * static_cast<uint64_t>(N);
    uint64_t val = 1;
    for (int i = 0; i < k; ++i)
        val = (val * P) % mod;
    if ((val & 1) == 0)
        throw std::logic_error("Frobenius galois element is even");
    return val;
}

std::vector<seal::Ciphertext>
frobenius_images(const seal::Ciphertext& ct_m,
                 const seal::SEALContext& context,
                 const seal::GaloisKeys& galois_keys)
{
    std::vector<seal::Ciphertext> T(D);
    T[0] = ct_m;
    seal::Evaluator ev(context);
    const uint32_t frobenius = static_cast<uint32_t>(frobenius_galois_elt(1));
    for (int k = 1; k < D; ++k)
        ev.apply_galois(T[k - 1], frobenius, galois_keys, T[k]);
    return T;
}

int current_depth(const seal::Ciphertext& ct, const seal::SEALContext& ctx)
{
    auto first = ctx.first_context_data()->chain_index();
    auto curr  = ctx.get_context_data(ct.parms_id())->chain_index();
    return static_cast<int>(first - curr);
}

void print_depth(const std::string& label,
                 const seal::Ciphertext& ct,
                 const seal::SEALContext& ctx)
{
    int depth = current_depth(ct, ctx);
    std::cout << "[depth " << std::setw(3) << depth << "]  " << label << "\n";
}

} // namespace bfv_inv
