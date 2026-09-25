// main.cpp -- BFV slot inversion demo, SEAL 4.1
#include "bfv_inversion.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <sstream>
#include <cctype>

using namespace bfv_inv;
using namespace seal;

static uint64_t mod_inv(uint64_t m, uint64_t p)
{
    uint64_t r = 1, b = m % p, e = p - 2;
    for (; e; e >>= 1) { if (e & 1) r = r * b % p; b = b * b % p; }
    return r;
}

static Plaintext encode_const(uint64_t m)
{
    std::ostringstream oss;
    oss << std::hex << m;
    return Plaintext(oss.str());
}

static uint64_t decode_const(const Plaintext& pt)
{
    return pt.is_zero() ? 0 : pt[0];
}

// Returns 0 if FMC, NT, and m*inv==1 all succeed.
static int run_one(uint64_t requested_val)
{
    const uint64_t test_val = std::min(requested_val, P - 2);
    if (test_val == 0) {
        std::cerr << "value must be in [1," << P - 1 << "]\n";
        return 1;
    }
    const uint64_t expected = mod_inv(test_val, P);
    std::cout << "m         = " << test_val << "\n"
              << "m^{-1}    = " << expected << "\n"
              << "m*m^{-1}  = " << (test_val * expected % P) << "\n\n";

    std::cout << "Setting up BFV context (p=" << P << ", N=" << N << ")...\n";
    SEALContext ctx  = make_context();
    std::cout << "Generating keys...\n";
    Keys        keys = generate_keys(ctx);
    Evaluator   ev(ctx);
    Decryptor   dec(ctx, keys.sec);
    Encryptor   enc(ctx, keys.pub);
    std::cout << "Ready.\n\n";

    Plaintext  pt_m = encode_const(test_val);
    Ciphertext ct_m;
    enc.encrypt(pt_m, ct_m);

    std::cout << "Running FMC (depth " << fmc_depth() << ")...\n";
    TimingBreakdown fmc_tb;
    Ciphertext ct_inv;
    try {
        ct_inv = fmc_invert(ct_m, ctx, keys.relin, keys.galois, &fmc_tb);
    } catch (const std::exception& e) {
        std::cerr << "FMC failed: " << e.what() << "\n";
        return 1;
    }

    Plaintext pt_inv;
    dec.decrypt(ct_inv, pt_inv);
    uint64_t fmc_result = decode_const(pt_inv);
    bool ok = (fmc_result == expected);
    std::cout << "FMC result = " << fmc_result << (ok ? "  CORRECT\n" : "  WRONG\n");
    print_depth("FMC", ct_inv, ctx);
    std::cout << "FMC wall " << std::fixed << std::setprecision(0)
              << fmc_tb.total_ms << " ms\n";

    uint64_t prod_val = 0;
    try {
        Ciphertext ct_m_at_inv = ct_m;
        ev.mod_switch_to_inplace(ct_m_at_inv, ct_inv.parms_id());
        Ciphertext ct_prod;
        ev.multiply(ct_m_at_inv, ct_inv, ct_prod);
        ev.relinearize_inplace(ct_prod, keys.relin);
        Plaintext pt_prod;
        dec.decrypt(ct_prod, pt_prod);
        prod_val = decode_const(pt_prod);
        std::cout << "m * FMC(m) mod p = " << prod_val
                  << (prod_val == 1 ? "  (=1 OK)\n\n" : "  WRONG\n\n");
    } catch (const std::exception& e) {
        std::cerr << "check multiply failed: " << e.what() << "\n\n";
        prod_val = 0;
    }

    std::cout << "Running Norm Tower (depth " << nt_depth() << ")...\n";
    TimingBreakdown nt_tb;
    Ciphertext ct_nt;
    try {
        ct_nt = norm_tower_invert(ct_m, ctx, keys.relin, keys.galois, &nt_tb);
    } catch (const std::exception& e) {
        std::cerr << "NT failed: " << e.what() << "\n";
        return 1;
    }

    Plaintext pt_nt;
    dec.decrypt(ct_nt, pt_nt);
    uint64_t nt_result = decode_const(pt_nt);
    bool ok_nt = (nt_result == expected);
    std::cout << "NT  result = " << nt_result << (ok_nt ? "  CORRECT\n" : "  WRONG\n");
    print_depth("NT ", ct_nt, ctx);
    std::cout << "NT  wall " << std::fixed << std::setprecision(0)
              << nt_tb.total_ms << " ms\n";

    std::cout << "\nFMC depth=" << fmc_depth()
              << "   NT depth=" << nt_depth() << "\n";
    return (ok && ok_nt && prod_val == 1) ? 0 : 1;
}

int main(int argc, char** argv)
{
    std::cout << std::unitbuf;
    std::cerr << std::unitbuf;
    std::string preset = "main";
    uint64_t test_val = 12345;
    bool run_all = false;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-h" || a == "--help" || a == "--list") {
            list_presets();
            return 0;
        }
        if (a == "--all") { run_all = true; continue; }
        if (a == "--preset" && i + 1 < argc) { preset = argv[++i]; continue; }
        if (!a.empty() && std::isalpha(static_cast<unsigned char>(a[0]))) {
            preset = a;
            continue;
        }
        if (!a.empty() && a[0] == '-') {
            std::cerr << "unknown option: " << a << "\n";
            list_presets();
            return 2;
        }
        try {
            test_val = std::stoull(a);
        } catch (const std::exception&) {
            std::cerr << "not a number or preset: " << a << "\n";
            list_presets();
            return 2;
        }
    }

    if (run_all) {
        const auto& tags = all_128bit_tags();
        int fails = 0;
        std::cout << "Running " << tags.size()
                  << " 128-bit sets, smallest N first\n\n";
        std::cout << std::left << std::setw(36) << "tag"
                  << std::right << std::setw(8) << "FMC"
                  << std::setw(8) << "NT"
                  << std::setw(8) << "ok?" << "\n";
        for (size_t i = 0; i < tags.size(); ++i) {
            apply_preset(tags[i]);
            set_verbose(false);
            std::cout << "\n======== [" << (i + 1) << "/" << tags.size()
                      << "] " << tags[i] << " ========\n";
            int rc = 1;
            try {
                rc = run_one(test_val);
            } catch (const std::exception& e) {
                std::cerr << "set failed: " << e.what() << "\n";
            }
            if (rc) ++fails;
            std::cout << std::left << std::setw(36) << tags[i]
                      << std::right << (rc ? "  FAIL\n" : "  ok\n");
        }
        std::cout << "\n" << (tags.size() - fails) << "/" << tags.size()
                  << " sets passed.  " << fails << " failed.\n";
        return fails ? 1 : 0;
    }

    apply_preset(preset);
    return run_one(test_val);
}
