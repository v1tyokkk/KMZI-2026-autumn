#include "bigint.hpp"

#include <cstring>

namespace {

int hex_digit(char c) {
    if (c >= '0' && c <= '9') {
        return c - '0';
    }
    if (c >= 'a' && c <= 'f') {
        return c - 'a' + 10;
    }
    if (c >= 'A' && c <= 'F') {
        return c - 'A' + 10;
    }
    return -1;
}

char hex_char(int v) { return v < 10 ? static_cast<char>('0' + v) : static_cast<char>('a' + v - 10); }

}

namespace u256 {

void zero(U256& v) { std::memset(v.w, 0, sizeof(v.w)); }

void set_word(U256& v, uint64_t x) {
    zero(v);
    v.w[0] = x;
}

bool from_hex(U256& out, const char* hex) {
    if (hex == nullptr) {
        return false;
    }
    if (hex[0] == '0' && (hex[1] == 'x' || hex[1] == 'X')) {
        hex += 2;
    }
    const std::size_t len = std::strlen(hex);
    if (len == 0 || len > 64) {
        return false;
    }
    U256 r;
    zero(r);
    for (std::size_t i = 0; i < len; ++i) {
        const int d = hex_digit(hex[i]);
        if (d < 0) {
            return false;
        }
        uint64_t carry = static_cast<uint64_t>(d);
        for (int j = 0; j < 4; ++j) {
            const unsigned __int128 t = (static_cast<unsigned __int128>(r.w[j]) << 4) | carry;
            r.w[j] = static_cast<uint64_t>(t);
            carry = static_cast<uint64_t>(t >> 64);
        }
        if (carry != 0) {
            return false;
        }
    }
    out = r;
    return true;
}

std::string to_hex(const U256& v) {
    std::string out(64, '0');
    for (int limb = 3; limb >= 0; --limb) {
        uint64_t x = v.w[limb];
        for (int i = 15; i >= 0; --i) {
            out[static_cast<std::size_t>((3 - limb) * 16 + i)] = hex_char(static_cast<int>(x & 0xF));
            x >>= 4;
        }
    }
    return out;
}

bool from_bytes(U256& out, const uint8_t* be, int len) {
    if (be == nullptr || len < 0 || len > 32) {
        return false;
    }
    zero(out);
    for (int i = 0; i < len; ++i) {
        for (int j = 3; j > 0; --j) {
            out.w[j] = (out.w[j] << 8) | (out.w[j - 1] >> 56);
        }
        out.w[0] = (out.w[0] << 8) | be[i];
    }
    return true;
}

void to_bytes(const U256& v, uint8_t* be) {
    for (int i = 0; i < 32; ++i) {
        const int limb = 3 - i / 8;
        const int shift = 8 * (7 - i % 8);
        be[i] = static_cast<uint8_t>(v.w[limb] >> shift);
    }
}

bool is_zero(const U256& v) { return (v.w[0] | v.w[1] | v.w[2] | v.w[3]) == 0; }

int compare(const U256& a, const U256& b) {
    for (int i = 3; i >= 0; --i) {
        if (a.w[i] != b.w[i]) {
            return a.w[i] < b.w[i] ? -1 : 1;
        }
    }
    return 0;
}

bool bit(const U256& v, int i) { return ((v.w[i / 64] >> (i % 64)) & 1u) != 0; }

U256 select(const U256& a, const U256& b, uint64_t bit_value) {
    const uint64_t mask = 0 - (bit_value & 1u);
    U256 out;
    for (int i = 0; i < 4; ++i) {
        out.w[i] = (a.w[i] & ~mask) | (b.w[i] & mask);
    }
    return out;
}

}

Mod::Mod(const char* modulus_hex) : n0_(0) {
    if (!u256::from_hex(m_, modulus_hex)) {
        u256::zero(m_);
    }
    uint64_t inv = 1;
    for (int i = 0; i < 6; ++i) {
        inv *= 2 - m_.w[0] * inv;
    }
    n0_ = 0 - inv;

    // R mod m и R^2 mod m, где R = 2^256.
    u256::set_word(one_, 1);
    for (int i = 0; i < 256; ++i) {
        one_ = add_raw(one_, one_);
    }
    r2_ = one_;
    for (int i = 0; i < 256; ++i) {
        r2_ = add_raw(r2_, r2_);
    }
}

U256 Mod::zero() const {
    U256 z;
    u256::zero(z);
    return z;
}

U256 Mod::from_small(uint64_t x) const {
    U256 v;
    u256::set_word(v, x);
    return to_mont(v);
}

U256 Mod::conditional_sub(const U256& a, uint64_t high) const {
    U256 d;
    uint64_t borrow = 0;
    for (int i = 0; i < 4; ++i) {
        const unsigned __int128 s =
            static_cast<unsigned __int128>(a.w[i]) - m_.w[i] - borrow;
        d.w[i] = static_cast<uint64_t>(s);
        borrow = static_cast<uint64_t>(s >> 64) & 1u;
    }
    const uint64_t ge = ((high != 0) || (borrow == 0)) ? 1u : 0u;
    return u256::select(a, d, ge);
}

U256 Mod::add_raw(const U256& a, const U256& b) const {
    uint64_t s[4];
    uint64_t carry = 0;
    for (int i = 0; i < 4; ++i) {
        const unsigned __int128 t = static_cast<unsigned __int128>(a.w[i]) + b.w[i] + carry;
        s[i] = static_cast<uint64_t>(t);
        carry = static_cast<uint64_t>(t >> 64);
    }
    U256 r{{s[0], s[1], s[2], s[3]}};
    return conditional_sub(r, carry);
}

U256 Mod::add(const U256& a, const U256& b) const { return add_raw(a, b); }

U256 Mod::sub(const U256& a, const U256& b) const {
    uint64_t d[4];
    uint64_t borrow = 0;
    for (int i = 0; i < 4; ++i) {
        const unsigned __int128 s =
            static_cast<unsigned __int128>(a.w[i]) - b.w[i] - borrow;
        d[i] = static_cast<uint64_t>(s);
        borrow = static_cast<uint64_t>(s >> 64) & 1u;
    }
    const uint64_t mask = 0 - borrow;
    uint64_t carry = 0;
    U256 r;
    for (int i = 0; i < 4; ++i) {
        const unsigned __int128 s =
            static_cast<unsigned __int128>(d[i]) + (m_.w[i] & mask) + carry;
        r.w[i] = static_cast<uint64_t>(s);
        carry = static_cast<uint64_t>(s >> 64);
    }
    return r;
}

// t = a * b * R^{-1} mod m, R = 2^256.
void Mod::mont_mul(U256& out, const U256& a, const U256& b) const {
    uint64_t t[9] = {};

    for (int i = 0; i < 4; ++i) {
        uint64_t carry = 0;
        for (int j = 0; j < 4; ++j) {
            const unsigned __int128 uv =
                static_cast<unsigned __int128>(a.w[i]) * b.w[j] + t[i + j] + carry;
            t[i + j] = static_cast<uint64_t>(uv);
            carry = static_cast<uint64_t>(uv >> 64);
        }
        for (int k = i + 4; k < 9; ++k) {
            const unsigned __int128 uv = static_cast<unsigned __int128>(t[k]) + carry;
            t[k] = static_cast<uint64_t>(uv);
            carry = static_cast<uint64_t>(uv >> 64);
        }
    }

    for (int i = 0; i < 4; ++i) {
        const uint64_t u = t[i] * n0_;
        uint64_t carry = 0;
        for (int j = 0; j < 4; ++j) {
            const unsigned __int128 uv =
                static_cast<unsigned __int128>(u) * m_.w[j] + t[i + j] + carry;
            t[i + j] = static_cast<uint64_t>(uv);
            carry = static_cast<uint64_t>(uv >> 64);
        }
        for (int k = i + 4; k < 9; ++k) {
            const unsigned __int128 uv = static_cast<unsigned __int128>(t[k]) + carry;
            t[k] = static_cast<uint64_t>(uv);
            carry = static_cast<uint64_t>(uv >> 64);
        }
    }

    U256 r{{t[4], t[5], t[6], t[7]}};
    out = conditional_sub(r, t[8]);
}

U256 Mod::mul(const U256& a, const U256& b) const {
    U256 r;
    mont_mul(r, a, b);
    return r;
}

U256 Mod::to_mont(const U256& a) const {
    U256 r;
    mont_mul(r, a, r2_);
    return r;
}

U256 Mod::from_mont(const U256& a) const {
    U256 one;
    u256::set_word(one, 1);
    U256 r;
    mont_mul(r, a, one);
    return r;
}

U256 Mod::inv(const U256& a) const {
    U256 exp = m_;
    exp.w[0] -= 2;
    U256 result = one_;
    const U256 base = a;
    for (int i = 255; i >= 0; --i) {
        result = mul(result, result);
        const U256 prod = mul(result, base);
        result = u256::select(result, prod, u256::bit(exp, i) ? 1u : 0u);
    }
    return result;
}

U256 Mod::reduce(const U256& a) const {
    U256 r = a;
    for (int i = 0; i < 2; ++i) {
        r = conditional_sub(r, 0);
    }
    return r;
}

bool Mod::equal(const U256& a, const U256& b) const { return u256::compare(a, b) == 0; }

bool Mod::is_zero(const U256& a) const { return u256::is_zero(a); }
