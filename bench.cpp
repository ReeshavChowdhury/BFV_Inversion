// bench.cpp -- benchmark FMC (parallel) vs Norm Tower
// Prints the breakdown used by Table tab:timing in the experimental section.
#include "bfv_inversion.h"
#include <chrono>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <sstream>
#include <thread>
#include <cmath>
#include <cctype>
#include <string>

using namespace bfv_inv;
using namespace seal;

static Plaintext encode_const(uint64_t m)
{
    std::ostringstream o; o << std::hex << m;
    return Plaintext(o.str());
}

int main(int argc, char** argv)
{
    std::cout << std::unitbuf;
    std::cerr << std::unitbuf;
    int iters = 3;
    std::string preset = "main";
    bool run_all = false;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-h" || a == "--help" || a == "--list") { list_presets(); return 0; }
        if (a == "--all") { run_all = true; continue; }
        if (a == "--preset" && i + 1 < argc) { preset = argv[++i]; continue; }
        if (!a.empty() && std::isalpha(static_cast<unsigned char>(a[0]))) {
            preset = a; continue;
        }
        iters = std::atoi(a.c_str());
    }

    if (run_all) {
        const auto& tags = all_128bit_tags();
        std::cout << "bfv_bench --all  (" << tags.size()
                  << " sets, 1 run each, smallest N first)\n\n"
                  << std::left << std::setw(36) << "tag"
                  << std::right << std::setw(10) << "FMC_s"
                  << std::setw(10) << "NT_s" << "\n";
        set_verbose(false);
        for (size_t i = 0; i < tags.size(); ++i) {
            apply_preset(tags[i]);
            std::cout << "[" << (i + 1) << "/" << tags.size() << "] "
                      << tags[i] << std::flush;
            try {
                SEALContext ctx = make_context();
                Keys keys = generate_keys(ctx);
                Encryptor enc(ctx, keys.pub);
                uint64_t m = std::min<uint64_t>(12345, P - 2);
                Ciphertext ct; enc.encrypt(encode_const(m), ct);
                TimingBreakdown fmc, nt;
                Ciphertext c = ct;
                fmc_invert(c, ctx, keys.relin, keys.galois, &fmc);
                c = ct;
                norm_tower_invert(c, ctx, keys.relin, keys.galois, &nt);
                std::cout << "  FMC " << std::fixed << std::setprecision(1)
                          << fmc.total_ms / 1000.0 << "s  NT "
                          << nt.total_ms / 1000.0 << "s\n";
            } catch (const std::exception& e) {
                std::cout << "  FAIL " << e.what() << "\n";
            }
        }
        return 0;
    }
    apply_preset(preset);

    std::cout << "Hardware threads available: "
              << std::thread::hardware_concurrency() << "\n"
              << "Preset " << preset_name()
              << "  p=" << P << " N=" << N << " d=" << D
              << " L_F=" << LF << " s=" << S
              << " slots=" << (N / D)
              << " log2Q=" << logq_bits() << "\n\n";

    SEALContext ctx  = make_context();
    std::cout << "Generating keys...\n";
    Keys        keys = generate_keys(ctx);
    Encryptor   enc(ctx, keys.pub);
    std::cout << "Keys ready.\n\n";

    Plaintext  pt = encode_const(12345);
    Ciphertext ct_m; enc.encrypt(pt, ct_m);

    std::cout << "Running " << iters << " iterations each.\n\n";

    struct Row {
        std::vector<double> total, ph1, ph2, fro;
        int ops = 0;
    };

    auto run_fmc = [&]() {
        Row r;
        for (int i = 0; i < iters; ++i) {
            TimingBreakdown tb;
            Ciphertext c = ct_m;
            fmc_invert(c, ctx, keys.relin, keys.galois, &tb);
            r.total.push_back(tb.total_ms);
            r.ph1.push_back(tb.phase1_ms);
            r.ph2.push_back(tb.phase2_ms);
            r.fro.push_back(tb.frobenius_ms);
            r.ops = tb.mul_ops;
            std::cout << "  FMC run " << i + 1 << ": "
                      << std::fixed << std::setprecision(0) << tb.total_ms
                      << " ms  (Ph1 " << tb.phase1_ms
                      << "  Ph2 " << tb.phase2_ms << ")\n";
        }
        return r;
    };

    auto run_nt = [&]() {
        Row r;
        for (int i = 0; i < iters; ++i) {
            TimingBreakdown tb;
            Ciphertext c = ct_m;
            norm_tower_invert(c, ctx, keys.relin, keys.galois, &tb);
            r.total.push_back(tb.total_ms);
            r.ph1.push_back(tb.phase1_ms);
            r.ph2.push_back(tb.phase2_ms);
            r.fro.push_back(tb.frobenius_ms);
            r.ops = tb.mul_ops;
            std::cout << "  NT  run " << i + 1 << ": "
                      << std::fixed << std::setprecision(0) << tb.total_ms
                      << " ms\n";
        }
        return r;
    };

    // Quiet the per-phase banners; bench prints its own lines.
    set_verbose(false);

    Row fmc = run_fmc();
    std::cout << "\n";
    Row nt  = run_nt();

    auto med = [](std::vector<double> v) {
        std::sort(v.begin(), v.end());
        return v[v.size() / 2];
    };
    auto at = [](const std::vector<double>& v, int i) {
        return i < (int)v.size() ? v[i] : 0.0;
    };

    const double slots = static_cast<double>(N / D);
    auto thru = [&](double ms) {
        return slots / (ms / 1000.0);
    };

    auto print_row = [&](const char* name, const std::vector<double>& v,
                         int ops) {
        double m = med(v);
        std::cout << std::left << std::setw(32) << name
                  << std::right << std::fixed << std::setprecision(1)
                  << std::setw(10) << at(v, 0) / 1000.0
                  << std::setw(10) << at(v, 1) / 1000.0
                  << std::setw(10) << at(v, 2) / 1000.0
                  << std::setw(10) << m / 1000.0
                  << std::setw(8)  << ops
                  << std::setw(14) << std::setprecision(0) << thru(m)
                  << " slots/s\n";
    };

    std::cout << "\n--- Table tab:timing (seconds) ---\n"
              << std::left << std::setw(32) << "Operation"
              << std::right
              << std::setw(10) << "R1"
              << std::setw(10) << "R2"
              << std::setw(10) << "R3"
              << std::setw(10) << "Med.(s)"
              << std::setw(8)  << "Ops"
              << std::setw(14) << "Throughput\n";

    print_row("FMC Ph.1: exponentiation", fmc.ph1, fmc.ops - S);
    print_row("FMC Ph.2: product tree (s=3)",  fmc.ph2, S);
    print_row("FMC total",                     fmc.total, fmc.ops);
    print_row("NT total",                      nt.total,  nt.ops);

    std::cout << "\nFMC (depth " << fmc_depth() << ") median: "
              << std::fixed << std::setprecision(1) << med(fmc.total) / 1000.0 << " s\n"
              << "NT  (depth " << nt_depth()  << ") median: "
              << med(nt.total) / 1000.0 << " s\n"
              << "Ratio FMC/NT: " << std::setprecision(3)
              << med(fmc.total) / med(nt.total) << "x\n";
    return 0;
}
