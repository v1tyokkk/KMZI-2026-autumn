#pragma once

#include <array>
#include <cstdint>

class Gf8 {
public:
    static constexpr uint16_t kVariantPoly = 0x14D;  // x^8 + x^6 + x^3 + x^2 + 1

    explicit Gf8(uint16_t poly = kVariantPoly);

    uint16_t polynomial() const { return poly_; }
    uint8_t p_byte() const { return p_byte_; }

    static constexpr uint8_t add(uint8_t a, uint8_t b) { return static_cast<uint8_t>(a ^ b); }

    uint8_t multiply(uint8_t a, uint8_t b) const;

    uint8_t inverse(uint8_t a) const;
    uint8_t xtime(uint8_t a) const;

    void generate_sbox(std::array<uint8_t, 256>& sbox) const;
    void generate_inv_sbox(const std::array<uint8_t, 256>& sbox,
                           std::array<uint8_t, 256>& inv_sbox) const;

    static uint8_t affine_aes(uint8_t x);

private:
    uint16_t poly_;
    uint8_t p_byte_;

    static uint8_t mul4(uint8_t a, uint8_t b);
    uint8_t reduce(uint16_t value) const;
};
