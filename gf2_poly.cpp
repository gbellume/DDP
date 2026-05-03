/**
 * gf2_poly.cpp – BCH/CRC codec over GF(2^m)
 *
 * Supports arbitrarily large data words via BigPoly (vector<uint64_t>).
 * GF-level operations (gf_mult, gf_inv, BM, Chien) use plain int/uint64_t
 * since field elements are always < 2^m.
 *
 * Optimisations:
 *  1. Hot GF functions force-inlined.
 *  2. GF tables in a single contiguous heap allocation.
 *  3. poly_div_gf2: bit_length maintained incrementally.
 *  4. calculate_syndromes_fast: iterates only set bits (O(popcount)).
 *  5. bch_berlekamp_massey: O(1) shift via integer offset.
 */

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <map>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// =============================================================================
// BigPoly – polynomial over GF(2) of arbitrary degree
//
// Stored as a little-endian vector of uint64_t words:
//   w[0] holds bits  0..63
//   w[1] holds bits 64..127
//   ...
// The representation is always trimmed (no trailing zero words).
// =============================================================================

struct BigPoly {
    std::vector<uint64_t> w;   // little-endian limbs

    BigPoly() = default;
    explicit BigPoly(uint64_t v) { if (v) w.push_back(v); }

    bool empty() const { return w.empty(); }

    int bit_length() const {
        for (int i = static_cast<int>(w.size()) - 1; i >= 0; --i)
            if (w[i]) return i * 64 + (64 - __builtin_clzll(w[i]));
        return 0;
    }

    int bit(int pos) const {
        int wi = pos / 64, bi = pos % 64;
        if (wi >= static_cast<int>(w.size())) return 0;
        return static_cast<int>((w[wi] >> bi) & 1);
    }

    void trim() { while (!w.empty() && w.back() == 0) w.pop_back(); }

    BigPoly& operator^=(const BigPoly& o) {
        if (o.w.size() > w.size()) w.resize(o.w.size(), 0);
        for (int i = 0; i < static_cast<int>(o.w.size()); ++i) w[i] ^= o.w[i];
        trim();
        return *this;
    }
    BigPoly operator^(const BigPoly& o) const { BigPoly r = *this; r ^= o; return r; }

    // XOR in (val << shift_bits) without building a full BigPoly temporary.
    // val is a small uint64_t (e.g. the divisor in poly division).
    void xor_shifted_u64(uint64_t val, int shift) {
        int wi = shift / 64, bi = shift % 64;
        int needed = wi + (bi ? 2 : 1);
        if (static_cast<int>(w.size()) < needed) w.resize(needed, 0);
        w[wi] ^= val << bi;
        if (bi) w[wi + 1] ^= val >> (64 - bi);
        trim();
    }

    BigPoly operator<<(int n) const {
        if (w.empty() || n == 0) { BigPoly r = *this; return r; }
        int wi = n / 64, bi = n % 64;
        BigPoly res;
        res.w.resize(w.size() + wi + 1, 0);
        for (int i = 0; i < static_cast<int>(w.size()); ++i) {
            res.w[i + wi] |= w[i] << bi;
            if (bi) res.w[i + wi + 1] |= w[i] >> (64 - bi);
        }
        res.trim();
        return res;
    }

    bool operator==(const BigPoly& o) const { return w == o.w; }
    bool operator!=(const BigPoly& o) const { return w != o.w; }
};

/// Generate a random BigPoly with exactly `bits` random bits.
BigPoly random_bigpoly(int bits, std::mt19937_64& rng) {
    BigPoly r;
    if (bits <= 0) return r;
    int nwords = (bits + 63) / 64;
    r.w.resize(nwords, 0);
    for (auto& x : r.w) x = rng();
    int rem = bits % 64;
    if (rem) r.w.back() &= (uint64_t(1) << rem) - 1;
    r.trim();
    return r;
}

// =============================================================================
// GF table storage
// =============================================================================

static int* exp_table = nullptr;
static int* log_table = nullptr;
static int  gf_n      = 0;

static std::unique_ptr<int[]> gf_storage;

// =============================================================================
// Force-inlined GF helpers
// =============================================================================

[[gnu::always_inline]] static inline int gf_mult(int a, int b) {
    if (__builtin_expect(a == 0 || b == 0, 0)) return 0;
    return exp_table[(log_table[a] + log_table[b]) % gf_n];
}

[[gnu::always_inline]] static inline int gf_inv(int a) {
    if (a == 0) throw std::domain_error("0 has no inverse in GF");
    return exp_table[(gf_n - log_table[a]) % gf_n];
}

// =============================================================================
// GF(2) polynomial arithmetic  (small polys stay uint64_t)
// =============================================================================

/// Multiply two small GF(2) polynomials (used only for building g).
uint64_t poly_mult_gf2(uint64_t p1, uint64_t p2) {
    uint64_t result = 0;
    while (p2) {
        if (p2 & 1) result ^= p1;
        p1 <<= 1; p2 >>= 1;
    }
    return result;
}

/// Divide BigPoly dividend by a small uint64_t divisor; returns remainder.
/// The divisor (g polynomial) has degree ≤ 26, always fits in uint64_t.
BigPoly poly_div_gf2(BigPoly dividend, uint64_t divisor) {
    int div_len = 64 - __builtin_clzll(divisor);   // bit_length of divisor
    int dd_len  = dividend.bit_length();
    while (dd_len >= div_len) {
        dividend.xor_shifted_u64(divisor, dd_len - div_len);
        dd_len = dividend.bit_length();
    }
    return dividend;
}

/// Convenience overload: small dividend, small divisor (used for CRC on
/// uint64_t values and for the CRC key itself).
uint64_t poly_div_gf2_u64(uint64_t dividend, uint64_t divisor) {
    int div_len = 64 - __builtin_clzll(divisor);
    int dd_len  = dividend ? 64 - __builtin_clzll(dividend) : 0;
    while (dd_len >= div_len) {
        dividend ^= divisor << (dd_len - div_len);
        dd_len = dividend ? 64 - __builtin_clzll(dividend) : 0;
    }
    return dividend;
}

// =============================================================================
// GF(2^m) table generation
// =============================================================================

int generate_gf_tables(int m, int prim_poly = -1) {
    static const std::map<int,int> standard_polys = {
        {3, 0b1011}, {4, 0x13}, {5, 0b100101},
        {8, 0x11D }, {13, 0x201B}, {16, 0x1100B},
    };
    if (prim_poly == -1) {
        auto it = standard_polys.find(m);
        if (it == standard_polys.end())
            throw std::invalid_argument("No primitive polynomial for m=" + std::to_string(m));
        prim_poly = it->second;
    }
    int field_size = 1 << m;
    gf_n           = field_size - 1;
    gf_storage     = std::make_unique<int[]>(field_size * 2);
    exp_table      = gf_storage.get();
    log_table      = gf_storage.get() + field_size;
    int x = 1;
    for (int i = 0; i < gf_n; ++i) {
        exp_table[i] = x; log_table[x] = i;
        x <<= 1;
        if (x & field_size) x ^= prim_poly;
    }
    exp_table[gf_n] = 1;
    return field_size;
}

// =============================================================================
// BCH generator polynomial
// =============================================================================

uint64_t get_g_polynomial(int m, int t = 2) {
    static const std::map<int, std::array<uint64_t,2>> min_polys = {
        {3,  {0b1011,   0b1101  }},
        {4,  {0b10011,  0b11111 }},
        {5,  {0b100101, 0b111101}},
        {13, {0x201B,   0x26B1  }},
    };
    auto it = min_polys.find(m);
    if (it == min_polys.end())
        throw std::invalid_argument("No data for m=" + std::to_string(m));
    if (t == 1) return it->second[0];
    if (t == 2) return poly_mult_gf2(it->second[0], it->second[1]);
    throw std::invalid_argument("Only t=1 or t=2 supported.");
}

// =============================================================================
// Systematic BCH encoding
// =============================================================================

BigPoly encoding(const BigPoly& word, uint64_t g) {
    int     parity_bits = 64 - __builtin_clzll(g) - 1;   // degree of g
    BigPoly word_s      = word << parity_bits;
    BigPoly remainder   = poly_div_gf2(word_s, g);
    return word_s ^ remainder;
}

// =============================================================================
// Syndrome calculation
// =============================================================================

/// Horner evaluation of BigPoly at alpha^root_idx.
int calculate_horner(const BigPoly& msg, int root_idx) {
    int root = exp_table[root_idx];
    int s    = 0;
    int len  = msg.bit_length();
    for (int bp = len - 1; bp >= 0; --bp)
        s = gf_mult(s, root) ^ msg.bit(bp);
    return s;
}

std::vector<int> calculate_syndromes(const BigPoly& msg, int t) {
    int s1 = calculate_horner(msg, 1);
    int s2 = gf_mult(s1, s1);
    if (t == 1) return {s1, s2};
    int s3 = calculate_horner(msg, 3);
    int s4 = gf_mult(s2, s2);
    return {s1, s2, s3, s4};
}

/// Fast syndrome: iterates only over set bits of each limb using ctzll.
std::vector<int> calculate_syndromes_fast(const BigPoly& msg, int t = 2, int m = 13) {
    const int n = (1 << m) - 1;
    int s1 = 0, s3 = 0;
    for (int wi = 0; wi < static_cast<int>(msg.w.size()); ++wi) {
        for (uint64_t tmp = msg.w[wi]; tmp; tmp &= tmp - 1) {
            int bit_idx = wi * 64 + __builtin_ctzll(tmp);
            s1 ^= exp_table[bit_idx % n];
            if (t == 2)
                s3 ^= exp_table[static_cast<int>((static_cast<long long>(bit_idx) * 3) % n)];
        }
    }
    int s2 = s1 ? exp_table[(log_table[s1] * 2) % n] : 0;
    if (t == 1) return {s1, s2};
    int s4 = s2 ? exp_table[(log_table[s2] * 2) % n] : 0;
    return {s1, s2, s3, s4};
}

// =============================================================================
// Berlekamp-Massey  (operates on GF elements – no BigPoly needed)
// =============================================================================

std::vector<int> bch_berlekamp_massey(const std::vector<int>& syndromes, int t = 2) {
    std::vector<int> lambda_poly = {1};
    std::vector<int> b_poly      = {1};
    int b_offset = 0, L = 0;

    for (int k = 1; k <= 2 * t; ++k) {
        int d = syndromes[k - 1];
        for (int i = 1; i <= L; ++i)
            if (i < static_cast<int>(lambda_poly.size()))
                d ^= gf_mult(lambda_poly[i], syndromes[k - 1 - i]);
        ++b_offset;
        if (d != 0) {
            std::vector<int> T = lambda_poly;
            int needed = static_cast<int>(b_poly.size()) + b_offset;
            if (needed > static_cast<int>(lambda_poly.size()))
                lambda_poly.resize(needed, 0);
            for (int i = 0; i < static_cast<int>(b_poly.size()); ++i)
                lambda_poly[i + b_offset] ^= gf_mult(d, b_poly[i]);
            if (2 * L < k) {
                L = k - L;
                int d_inv = gf_inv(d);
                b_poly.resize(T.size());
                for (std::size_t i = 0; i < T.size(); ++i)
                    b_poly[i] = gf_mult(T[i], d_inv);
                b_offset = 0;
            }
        }
    }
    lambda_poly.resize(t + 1, 0);
    return lambda_poly;
}

// =============================================================================
// Chien Search  (operates on GF elements – no BigPoly needed)
// =============================================================================

std::vector<int> bch_chien_search(const std::vector<int>& lambda_poly,
                                  int message_length, int m = 13) {
    std::vector<int> out;
    const int n     = (1 << m) - 1;
    int alpha_inv_1 = exp_table[n - 1];
    int alpha_inv_2 = exp_table[n - 2];
    int term1 = lambda_poly[1];
    int term2 = lambda_poly[2];
    for (int i = 0; i < message_length; ++i) {
        if ((1 ^ term1 ^ term2) == 0) out.push_back(i);
        term1 = gf_mult(term1, alpha_inv_1);
        term2 = gf_mult(term2, alpha_inv_2);
    }
    return out;
}

std::vector<int> bch_chien_search_fast(const std::vector<int>& lambda_poly,
                                       int message_length, int m = 13) {
    const int n = (1 << m) - 1;
    int L1 = lambda_poly[1], L2 = lambda_poly[2];
    if (L1 == 0 && L2 == 0) return {};
    if (L2 == 0) {
        int pos = log_table[L1];
        return (pos < message_length) ? std::vector<int>{pos} : std::vector<int>{};
    }
    std::vector<int> out;
    int log_t1 = log_table[L1];
    int log_t2 = log_table[L2];
    for (int i = 0; i < message_length; ++i) {
        if ((1 ^ exp_table[log_t1] ^ exp_table[log_t2]) == 0)
            out.push_back(i);
        log_t1 = (log_t1 - 1 + n) % n;
        log_t2 = (log_t2 - 2 + n) % n;
    }
    return out;
}

// =============================================================================
// CRC-32  (key = 0xBADB17CC, 32 parity bits appended)
// =============================================================================

BigPoly CRC_encoding(const BigPoly& word) {
    constexpr uint64_t key = 0xBADB17CC;
    BigPoly word_s = word << 32;
    return word_s ^ poly_div_gf2(word_s, key);
}

bool CRC_check(const BigPoly& received) {
    constexpr uint64_t key = 0xBADB17CC;
    return poly_div_gf2(received, key).empty();
}

// =============================================================================
// Helpers for BigPoly bit-flip (used in test)
// =============================================================================

void flip_bit(BigPoly& p, int pos) {
    int wi = pos / 64, bi = pos % 64;
    if (wi >= static_cast<int>(p.w.size())) p.w.resize(wi + 1, 0);
    p.w[wi] ^= uint64_t(1) << bi;
    p.trim();
}

// =============================================================================
// Entry point  (single-word smoke test)
// =============================================================================

#ifndef BCH_NO_MAIN
int main() {
    constexpr int DATA_BYTES = 590;
    constexpr int BIT_FLIPS  = 2;
    constexpr int m          = 13;

    generate_gf_tables(m);
    const uint64_t g = get_g_polynomial(m);

    std::mt19937_64 rng(std::random_device{}());

    BigPoly original_data = random_bigpoly(DATA_BYTES * 8, rng);
    BigPoly word_crc      = CRC_encoding(original_data);
    BigPoly codeword      = encoding(word_crc, g);

    BigPoly corrupted = codeword;
    int cw_len = codeword.bit_length();
    std::uniform_int_distribution<int> dist(0, cw_len - 1);
    for (int i = 0; i < BIT_FLIPS; ++i) flip_bit(corrupted, dist(rng));

    auto syndromes  = calculate_syndromes(corrupted, 2);
    auto lambda     = bch_berlekamp_massey(syndromes, 2);
    auto error_pos  = bch_chien_search(lambda, cw_len, m);

    BigPoly corrected = corrupted;
    for (int p : error_pos) flip_bit(corrected, p);

    std::printf("%s\n", corrected == codeword ? "CORRECTED" : "UNCORRECTABLE");
    std::printf("Data size:      %d bytes (%d bits)\n", DATA_BYTES, DATA_BYTES * 8);
    std::printf("Codeword size:  %d bits\n", cw_len);
    std::printf("Error positions: [");
    for (std::size_t i = 0; i < error_pos.size(); ++i)
        std::printf("%s%d", i ? ", " : "", error_pos[i]);
    std::printf("]\n");

    return 0;
}
#endif // BCH_NO_MAIN