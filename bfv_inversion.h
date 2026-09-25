#pragma once
// bfv_inversion.h -- SEAL 4.1 BFV slot inversion (parameterised)
#include <seal/seal.h>
#include <cstdint>
#include <vector>
#include <string>

namespace bfv_inv {

// Active parameter set. Call apply_preset() before make_context().
extern uint64_t P;
extern size_t   N;
extern int      D;
extern int      LF;
extern int      S;
extern int      V2;     // v2(p-1); p-1 = ODD · 2^V2
extern int      ODD;
extern std::vector<int> PRIME_BITS;

inline int fmc_depth() { return LF + S; }
inline int nt_depth()  { return LF + S + 1; }

void apply_preset(const std::string& name);
bool apply_tag(const std::string& tag);  // N{N}-p{p}-d{d}-{nprimes}x{pbits}
void list_presets();
// All inversion-feasible 128-bit tags, smallest N first (N <= 2^17).
const std::vector<std::string>& all_128bit_tags();
const char* preset_name();
int  logq_bits();          // sum of prime bit-lengths
int  homenc128_cap();      // HomEncStd 128-bit max log2 Q for this N
double albrecht_usvp_bits(); // cached rough usvp figure for this chain

struct TimingBreakdown {
    double frobenius_ms = 0;
    double phase1_ms    = 0;
    double phase2_ms    = 0;
    double total_ms     = 0;
    int    mul_ops      = 0;
};

void set_verbose(bool v);
bool verbose();

seal::SEALContext make_context();

struct Keys {
    seal::PublicKey  pub;
    seal::SecretKey  sec;
    seal::RelinKeys  relin;
    seal::GaloisKeys galois;
};
Keys generate_keys(const seal::SEALContext& ctx);

uint64_t frobenius_galois_elt(int k);

std::vector<seal::Ciphertext>
frobenius_images(const seal::Ciphertext&   ct_m,
                 const seal::SEALContext&  context,
                 const seal::GaloisKeys&   galois_keys);

seal::Ciphertext
binary_exp(const seal::Ciphertext& ct,
           uint64_t                exp,
           const seal::SEALContext& context,
           const seal::RelinKeys&  relin_keys);

// One square-chain for both m^{p-2} and m^{p-1}, when p-1 = 2^a + 2^b.
// Returns false for odd parts that are not one more than a power of two
// (the caller falls back to two separate exponentiations).
bool fused_p_exps(const seal::Ciphertext& ct,
                  const seal::SEALContext& context,
                  const seal::RelinKeys& relin_keys,
                  seal::Ciphertext& pm2,
                  seal::Ciphertext& pm1);

seal::Ciphertext
product_tree(std::vector<seal::Ciphertext>& v,
             seal::Evaluator&               evaluator,
             const seal::RelinKeys&         relin_keys);

seal::Ciphertext
product_tree_parallel(std::vector<seal::Ciphertext>& v,
                      const seal::SEALContext&       context,
                      const seal::RelinKeys&         relin_keys);

seal::Ciphertext
fmc_invert(const seal::Ciphertext&  ct_m,
           const seal::SEALContext& context,
           const seal::RelinKeys&   relin_keys,
           const seal::GaloisKeys&  galois_keys,
           TimingBreakdown*         stats = nullptr);

seal::Ciphertext
norm_tower_invert(const seal::Ciphertext&  ct_m,
                  const seal::SEALContext& context,
                  const seal::RelinKeys&   relin_keys,
                  const seal::GaloisKeys&  galois_keys,
                  TimingBreakdown*         stats = nullptr);

int  current_depth(const seal::Ciphertext& ct, const seal::SEALContext& ctx);
void print_depth(const std::string& label, const seal::Ciphertext& ct,
                 const seal::SEALContext& ctx);

} // namespace bfv_inv
