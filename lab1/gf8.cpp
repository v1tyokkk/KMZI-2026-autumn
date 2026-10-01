#include "gf8.hpp"

Gf8::Gf8(uint16_t poly) : poly_(poly), p_byte_(static_cast<uint8_t>(poly & 0xFF)) {}

uint8_t Gf8::mul4(uint8_t a, uint8_t b) {
    uint8_t r = 0;
    r ^= static_cast<uint8_t>(a) & static_cast<uint8_t>(-static_cast<int8_t>(b & 1));
    r ^= static_cast<uint8_t>(a << 1) &
         static_cast<uint8_t>(-static_cast<int8_t>((b >> 1) & 1));
    r ^= static_cast<uint8_t>(a << 2) &
         static_cast<uint8_t>(-static_cast<int8_t>((b >> 2) & 1));
    r ^= static_cast<uint8_t>(a << 3) &
         static_cast<uint8_t>(-static_cast<int8_t>((b >> 3) & 1));
    return r;
}

uint8_t Gf8::reduce(uint16_t value) const {
    for (int i = 14; i >= 8; --i) {
        const uint16_t bit = static_cast<uint16_t>((value >> i) & 1u);
        const uint16_t mask = static_cast<uint16_t>(-static_cast<int16_t>(bit));
        value ^= static_cast<uint16_t>(poly_ << (i - 8)) & mask;
    }
    return static_cast<uint8_t>(value);
}

uint8_t Gf8::multiply(uint8_t a, uint8_t b) const {
    const uint8_t a0 = static_cast<uint8_t>(a & 0x0F);
    const uint8_t a1 = static_cast<uint8_t>(a >> 4);
    const uint8_t b0 = static_cast<uint8_t>(b & 0x0F);
    const uint8_t b1 = static_cast<uint8_t>(b >> 4);

    const uint8_t p0 = mul4(a0, b0);
    const uint8_t p2 = mul4(a1, b1);
    const uint8_t p1 = mul4(static_cast<uint8_t>(a0 ^ a1),
                            static_cast<uint8_t>(b0 ^ b1));

    const uint16_t unreduced = static_cast<uint16_t>(
        (static_cast<uint16_t>(p2) << 8) ^
        (static_cast<uint16_t>(p0 ^ p1 ^ p2) << 4) ^
        p0);
    return reduce(unreduced);
}

uint8_t Gf8::xtime(uint8_t a) const {
    const uint8_t hi = static_cast<uint8_t>((a >> 7) & 1);
    const uint8_t mask = static_cast<uint8_t>(-static_cast<int8_t>(hi));
    return static_cast<uint8_t>((a << 1) ^ (p_byte_ & mask));
}

uint8_t Gf8::inverse(uint8_t a) const {
    // a^254 = a^128 * a^64 * a^32 * a^16 * a^8 * a^4 * a^2; 0^n = 0.
    uint8_t sq = a;
    uint8_t acc = 1;
    // sq проходит a^2, a^4, ..., a^128; перемножаем все, кроме a^1.
    sq = multiply(sq, sq);  // a^2
    acc = multiply(acc, sq);
    sq = multiply(sq, sq);  // a^4
    acc = multiply(acc, sq);
    sq = multiply(sq, sq);  // a^8
    acc = multiply(acc, sq);
    sq = multiply(sq, sq);  // a^16
    acc = multiply(acc, sq);
    sq = multiply(sq, sq);  // a^32
    acc = multiply(acc, sq);
    sq = multiply(sq, sq);  // a^64
    acc = multiply(acc, sq);
    sq = multiply(sq, sq);  // a^128
    acc = multiply(acc, sq);
    return acc;
}

uint8_t Gf8::affine_aes(uint8_t x) {
    const auto rotl = [](uint8_t v, int n) -> uint8_t {
        return static_cast<uint8_t>((v << n) | (v >> (8 - n)));
    };
    return static_cast<uint8_t>(x ^ rotl(x, 1) ^ rotl(x, 2) ^ rotl(x, 3) ^ rotl(x, 4) ^ 0x63);
}

void Gf8::generate_sbox(std::array<uint8_t, 256>& sbox) const {
    for (int i = 0; i < 256; ++i) {
        sbox[static_cast<std::size_t>(i)] = affine_aes(inverse(static_cast<uint8_t>(i)));
    }
}

void Gf8::generate_inv_sbox(const std::array<uint8_t, 256>& sbox,
                            std::array<uint8_t, 256>& inv_sbox) const {
    for (int i = 0; i < 256; ++i) {
        inv_sbox[sbox[static_cast<std::size_t>(i)]] = static_cast<uint8_t>(i);
    }
}
