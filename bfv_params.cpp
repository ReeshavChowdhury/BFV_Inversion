// bfv_params.cpp -- presets for several ring dimensions.
// Security figures: HomEncStd-2021 classical 128-bit cap, and the
// Albrecht rough usvp / ADPS16 number (same model as LWE.estimate.rough).
#include "bfv_inversion.h"
#include <stdexcept>
#include <cstring>
#include <iostream>
#include <atomic>
#include <cmath>
#include <algorithm>
#include <cstdio>

namespace bfv_inv {

uint64_t P  = 40961;
size_t   N  = 32768;
int      D  = 8;
int      LF = 16;
int      S  = 3;
int      V2 = 13;
int      ODD = 5;
std::vector<int> PRIME_BITS(22, 32);

static std::string g_preset = "n15s";
static std::atomic<bool> g_verbose{true};

void set_verbose(bool v) { g_verbose.store(v, std::memory_order_relaxed); }
bool verbose()           { return g_verbose.load(std::memory_order_relaxed); }

const char* preset_name() { return g_preset.c_str(); }

int logq_bits()
{
    int s = 0;
    for (int b : PRIME_BITS) s += b;
    return s;
}

int homenc128_cap()
{
    switch (N) {
    case 8192:   return 218;
    case 16384:  return 438;
    case 32768:  return 881;
    case 65536:  return 1761;
    case 131072: return 3523;
    default:     return 0;
    }
}

// Cached Albrecht-rough usvp bits for the built-in chains (estimate_security.py).
double albrecht_usvp_bits()
{
    const int lq = logq_bits();
    if (N == 32768 && lq == 880)  return 92.9;
    if (N == 32768 && lq == 704)  return 128.2;
    if (N == 32768 && lq == 1900) return 24.5;
    if (N == 65536 && lq == 920)  return 186.0;
    if (N == 65536 && lq == 1150) return 168.0;
    if (N == 65536 && lq == 1200) return 160.9;
    if (N == 131072 && lq == 1200) return 405.3;
    if (N == 131072 && lq == 1300) return 365.6;
    if (N == 131072 && lq == 1750) return 247.3;
    return -1.0;
}

static std::vector<int> repeats(int n, int bits)
{
    return std::vector<int>(n, bits);
}

const std::vector<std::string>& all_128bit_tags()
{
    // Smallest N first. N<=2^17 only.
    static const std::vector<std::string> tags = {
        "N32768-p12289-d16-21x33",
        "N32768-p40961-d8-22x32",
        "N32768-p114689-d4-22x32",
        "N32768-p163841-d2-22x32",
        "N65536-p40961-d16-23x60",
        "N65536-p65537-d2-20x60",
        "N65536-p114689-d8-23x60",
        "N65536-p147457-d8-24x58",
        "N65536-p163841-d4-23x60",
        "N65536-p188417-d16-25x56",
        "N65536-p557057-d4-25x56",
        "N65536-p1376257-d2-25x56",
        "N65536-p1769473-d2-25x56",
        "N131072-p65537-d4-21x60",
        "N131072-p114689-d16-24x60",
        "N131072-p147457-d16-25x60",
        "N131072-p163841-d8-24x60",
        "N131072-p557057-d8-26x60",
        "N131072-p1179649-d2-25x60",
        "N131072-p1376257-d4-26x60",
        "N131072-p1769473-d4-26x60",
    };
    return tags;
}

void list_presets()
{
    std::cout
        << "Usage:\n"
        << "  ./bfv_inversion                 paper parameters, test value 12345\n"
        << "  ./bfv_inversion 12345           same parameters, chosen test value\n"
        << "  ./bfv_inversion 12345 small     short test (smaller ring)\n"
        << "  ./bfv_bench 3                   time both circuits, 3 runs\n"
        << "  ./bfv_bench 3 small             time the short test, 3 runs\n"
        << "  ./bfv_inversion --list          this message\n\n"
        << "Paper parameters (the default):\n"
        << "  prime 163841, polynomial degree 131072, 16384 slots\n"
        << "  24 primes of 50 bits in the ciphertext modulus\n\n"
        << "Short test:\n"
        << "  prime 40961, polynomial degree 32768, 4096 slots\n"
        << "  22 primes of 32 bits in the ciphertext modulus\n";
}

bool apply_tag(const std::string& tag)
{
    unsigned long long n = 0, p = 0, d = 0, np = 0, pb = 0;
    if (std::sscanf(tag.c_str(), "N%llu-p%llu-d%llu-%llux%llu",
                    &n, &p, &d, &np, &pb) != 5)
        return false;
    if (n == 0 || p < 3 || d < 1 || np < 2 || pb < 20 || pb > 60)
        return false;
    g_preset = tag;
    P = static_cast<uint64_t>(p);
    N = static_cast<size_t>(n);
    D = static_cast<int>(d);
    uint64_t pm1 = P - 1;
    V2 = 0;
    while ((pm1 & 1) == 0) { pm1 >>= 1; ++V2; }
    ODD = static_cast<int>(pm1);
    LF = static_cast<int>(std::ceil(std::log2(static_cast<double>(P - 1))));
    S = 0;
    for (int t = 1; t < D; t <<= 1) ++S;
    PRIME_BITS.assign(static_cast<size_t>(np), static_cast<int>(pb));
    return true;
}

void apply_preset(const std::string& name)
{
    std::string id = name;
    std::transform(id.begin(), id.end(), id.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    auto set = [&](const char* tag, uint64_t p, size_t n, int d,
                   std::vector<int> bits) {
        g_preset = tag;
        P = p; N = n; D = d;
        uint64_t pm1 = p - 1;
        V2 = 0;
        while ((pm1 & 1) == 0) { pm1 >>= 1; ++V2; }
        ODD = static_cast<int>(pm1);
        LF = static_cast<int>(std::ceil(std::log2(static_cast<double>(p - 1))));
        S  = 0;
        for (int t = 1; t < d; t <<= 1) ++S;
        PRIME_BITS = std::move(bits);
    };

    if (apply_tag(name))
        return;

    if (id == "small" || id == "short" || id == "quick"
        || id == "n15s" || id == "n15-128" || id == "fast") {
        // log2 Q = 704. HomEncStd classical cap at N=2^15 is 881.
        set("small", 40961, 32768, 8, repeats(22, 32));
    } else if (id.empty() || id == "main" || id == "default" || id == "paper"
               || id == "n17s" || id == "paper128") {
        // 23 data levels: depth 22 for the Norm Tower plus one check
        // multiply, and one special key-switch prime. log2 Q = 1200.
        set("main", 163841, 131072, 8, repeats(24, 50));
    } else if (id == "n15") {
        // Wider 40-bit primes, HomEncStd 128, Albrecht usvp ≈ 93.
        set("n15", 40961, 32768, 8, repeats(22, 40));
    } else if (id == "n16") {
        // First power-of-two N where Albrecht usvp is >= 128 with enough levels.
        set("n16", 114689, 65536, 8, repeats(23, 40));
    } else if (id == "n16w") {
        set("n16w", 114689, 65536, 8, repeats(23, 50));
    } else if (id == "n17") {
        set("n17", 163841, 131072, 8, repeats(26, 50));
    } else {
        list_presets();
        throw std::invalid_argument("unknown parameter set: " + name);
    }
}

seal::SEALContext make_context()
{
    seal::EncryptionParameters parms(seal::scheme_type::bfv);
    parms.set_poly_modulus_degree(N);
    parms.set_plain_modulus(P);
    parms.set_coeff_modulus(seal::CoeffModulus::Create(N, PRIME_BITS));

    // SEAL's built-in table only knows N <= 2^15.
    const auto sec = (N <= 32768 && logq_bits() <= homenc128_cap())
                         ? seal::sec_level_type::tc128
                         : seal::sec_level_type::none;

    seal::SEALContext ctx(parms, /*expand_mod_chain=*/true, sec);

    const char* err = ctx.parameter_error_message();
    if (err && std::strlen(err) > 0 && std::strcmp(err, "valid") != 0)
        throw std::runtime_error(std::string("SEAL error: ") + err);

    if (verbose()) {
        auto key_primes  = ctx.key_context_data()->parms().coeff_modulus().size();
        auto data_primes = ctx.first_context_data()->parms().coeff_modulus().size();
        std::cout << "Plaintext prime p = " << P << "\n"
                  << "Polynomial degree N = " << N << "\n"
                  << "Slots in one ciphertext = " << (N / static_cast<size_t>(D)) << "\n"
                  << key_primes << " primes in the ciphertext modulus, "
                  << data_primes << " available for multiplications, log2 Q = "
                  << logq_bits() << "\n"
                  << "Frobenius Monomial Circuit depth = " << fmc_depth() << "\n"
                  << "Norm Tower depth = " << nt_depth() << "\n"
                  << "Rough security estimate ≈ ";
        double b = albrecht_usvp_bits();
        if (b > 0) std::cout << b << " bits\n";
        else       std::cout << "(see estimate_security.py)\n";
    }
    return ctx;
}

Keys generate_keys(const seal::SEALContext& ctx)
{
    Keys k;
    seal::KeyGenerator kg(ctx);
    k.sec = kg.secret_key();
    kg.create_public_key(k.pub);
    kg.create_relin_keys(k.relin);

    // One generator key (Frobenius κ_p) is enough: T_{k} = κ_p(T_{k-1}).
    // Smaller key material, same depth-0 Galois cost.
    std::vector<uint32_t> elts = {
        static_cast<uint32_t>(frobenius_galois_elt(1))
    };
    kg.create_galois_keys(elts, k.galois);
    return k;
}

} // namespace bfv_inv
