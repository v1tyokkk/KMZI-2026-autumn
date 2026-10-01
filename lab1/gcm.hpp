#pragma once

#include "aes128.hpp"

#include <array>
#include <cstdint>
#include <cstddef>

// AES-GCM (AEAD): CTR + GHASH в GF(2^128) по модулю
// f(x) = x^128 + x^7 + x^2 + x + 1.
class AesGcm {
public:
    static constexpr int kBlockSize = 16;
    static constexpr int kTagSize = 16;

    AesGcm(const std::array<uint8_t, Aes128::kKeySize>& key, const Gf8& gf);

    void encrypt(const uint8_t* iv, std::size_t iv_len,
                 const uint8_t* aad, std::size_t aad_len,
                 const uint8_t* pt, std::size_t pt_len,
                 uint8_t* ct,
                 std::array<uint8_t, kTagSize>& tag) const;

    bool decrypt(const uint8_t* iv, std::size_t iv_len,
                 const uint8_t* aad, std::size_t aad_len,
                 const uint8_t* ct, std::size_t ct_len,
                 const uint8_t tag[kTagSize],
                 uint8_t* pt) const;

    static bool ct_equal(const uint8_t* a, const uint8_t* b, std::size_t n);

private:
    Aes128 aes_;
    std::array<uint8_t, kBlockSize> h_{};

    void compute_j0(const uint8_t* iv, std::size_t iv_len, uint8_t j0[kBlockSize]) const;
    void ctr_crypt(const uint8_t j0[kBlockSize],
                   const uint8_t* in, std::size_t len, uint8_t* out) const;
    void ghash(const uint8_t* aad, std::size_t aad_len,
               const uint8_t* ct, std::size_t ct_len,
               uint8_t out[kBlockSize]) const;

    static void gf128_mul(uint8_t z[kBlockSize],
                          const uint8_t x[kBlockSize],
                          const uint8_t y[kBlockSize]);
    static void inc32(uint8_t block[kBlockSize]);
    static void xor_block(uint8_t dst[kBlockSize], const uint8_t src[kBlockSize]);
};
