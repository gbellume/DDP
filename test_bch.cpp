/**
 * test_bch.cpp – Stress test for the BCH/CRC codec (supports large words)
 *
 * Compile:
 *   g++ -std=c++17 -O2 -DBCH_NO_MAIN -o test_bch test_bch.cpp gf2_poly.cpp
 * Run:
 *   ./test_bch
 * Plot:
 *   python plot_EDAC_C.py
 */

#include <cstdint>
#include <cstdio>
#include <map>
#include <random>
#include <vector>

// ── Forward declarations (defined in gf2_poly.cpp) ───────────────────────────

struct BigPoly;   // full definition in gf2_poly.cpp – needed for declarations

// Re-declare BigPoly inline so this TU can use it without a shared header.
struct BigPoly {
    std::vector<uint64_t> w;
    BigPoly() = default;
    explicit BigPoly(uint64_t v) { if (v) w.push_back(v); }
    bool empty() const { return w.empty(); }
    int bit_length() const {
        for (int i = (int)w.size()-1; i >= 0; --i)
            if (w[i]) return i*64 + (64 - __builtin_clzll(w[i]));
        return 0;
    }
    void trim() { while (!w.empty() && w.back() == 0) w.pop_back(); }
    bool operator==(const BigPoly& o) const { return w == o.w; }
    bool operator!=(const BigPoly& o) const { return w != o.w; }
};

// Functions implemented in gf2_poly.cpp
int              generate_gf_tables(int m, int prim_poly = -1);
uint64_t         get_g_polynomial(int m, int t = 2);
BigPoly          random_bigpoly(int bits, std::mt19937_64& rng);
BigPoly          CRC_encoding(const BigPoly& word);
bool             CRC_check(const BigPoly& received);
BigPoly          encoding(const BigPoly& word, uint64_t g);
std::vector<int> calculate_syndromes(const BigPoly& msg, int t);
std::vector<int> bch_berlekamp_massey(const std::vector<int>& syndromes, int t);
std::vector<int> bch_chien_search(const std::vector<int>& lambda_poly,
                                  int message_length, int m);
void             flip_bit(BigPoly& p, int pos);

// ─────────────────────────────────────────────────────────────────────────────

int main() {
    constexpr int DATA_BYTES = 590;    // <── cambia qui per testare altre dimensioni
    constexpr int m          = 13;
    constexpr int t          = 2;
    constexpr int MAX_FLIPS  = 4;
    constexpr int N          = 10000;   // parole per scenario

    generate_gf_tables(m);
    const uint64_t g = get_g_polynomial(m, t);

    std::mt19937_64 rng(std::random_device{}());

    struct Stats { int total = 0, corrected = 0, failed = 0;
                      long long effective_sum = 0; };
    std::vector<Stats> stats(MAX_FLIPS + 1);

    std::printf("BCH stress test — %d bytes/parola, %d parole/scenario\n\n",
                DATA_BYTES, N);

    for (int flip_count = 0; flip_count <= MAX_FLIPS; ++flip_count) {
        Stats& s = stats[flip_count];

        for (int trial = 0; trial < N; ++trial) {
            // ── Genera parola non nulla ──────────────────────────────────────
            BigPoly data;
            do { data = random_bigpoly(DATA_BYTES * 8, rng); }
            while (data.empty());

            // ── Codifica ─────────────────────────────────────────────────────
            BigPoly word_crc = CRC_encoding(data);
            BigPoly codeword = encoding(word_crc, g);
            int     cw_len   = codeword.bit_length();

            // ── Corrompi: flip_count volte (duplicati ammessi, come Python) ────
            BigPoly corrupted = codeword;
            int effective_flips = 0;
            {
                std::uniform_int_distribution<int> bit_dist(0, cw_len - 1);
                std::map<int,int> flip_count_map;
                for (int k = 0; k < flip_count; ++k)
                    flip_count_map[bit_dist(rng)]++;
                for (auto& [pos, cnt] : flip_count_map) {
                    if (cnt % 2 == 1) { flip_bit(corrupted, pos); ++effective_flips; }
                }
            }

            // ── Decodifica ───────────────────────────────────────────────────
            auto syndromes = calculate_syndromes(corrupted, t);
            auto lambda    = bch_berlekamp_massey(syndromes, t);
            auto err_pos   = bch_chien_search(lambda, cw_len, m);

            BigPoly corrected = corrupted;
            for (int p : err_pos) flip_bit(corrected, p);

            ++s.total;
            s.effective_sum += effective_flips;
            if (corrected == codeword) ++s.corrected;
            else                       ++s.failed;
        }
    }

    // ── Tabella risultati ────────────────────────────────────────────────────
    std::printf("%-6s  %8s  %10s  %8s  %8s  %s\n",
                "Flips", "Total", "Corrected", "Failed", "Avg.eff", "Result");
    std::printf("%-6s  %8s  %10s  %8s  %8s  %s\n",
                "-----", "-----", "---------", "------", "-------", "------");
    for (int f = 0; f <= MAX_FLIPS; ++f) {
        const Stats& s = stats[f];
        double avg_eff = s.total > 0 ? (double)s.effective_sum / s.total : 0.0;
        const char* res = (f <= t)
            ? (s.failed == 0 ? "OK" : "UNEXPECTED FAILURE")
            : (s.corrected == 0 ? "OK (uncorrectable)" : "accidental (eff.flips < t)");
        std::printf("%-6d  %8d  %10d  %8d  %8.3f  %s\n",
                    f, s.total, s.corrected, s.failed, avg_eff, res);
    }

    // ── CSV per plot.py ──────────────────────────────────────────────────────
    FILE* csv = std::fopen("bch_results.csv", "w");
    if (csv) {
        std::fprintf(csv, "flips,success_rate\n");
        for (int f = 0; f <= MAX_FLIPS; ++f) {
            const Stats& s = stats[f];
            double rate = s.total > 0
                          ? static_cast<double>(s.corrected) / s.total
                          : 0.0;
            std::fprintf(csv, "%d,%.6f\n", f, rate);
        }
        std::fclose(csv);
        std::printf("\nCSV salvato: bch_results.csv\n");
    }

    return 0;
}