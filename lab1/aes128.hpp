#pragma once

#include "gf8.hpp"

#include <array>
#include <cstdint>
class Aes128 {
public:
    static constexpr int kBlockSize = 16;
    static constexpr int kKeySize = 16;
    static constexpr int kRounds = 10;

    Aes128(const std::array<uint8_t, kKeySize>& key, const Gf8& gf);

    const Gf8& gf() const { return gf_; }

    void encrypt_block(const uint8_t in[kBlockSize], uint8_t out[kBlockSize]) const;
    void decrypt_block(const uint8_t in[kBlockSize], uint8_t out[kBlockSize]) const;

private:
    Gf8 gf_;
    std::array<uint8_t, 256> sbox_{};
    std::array<uint8_t, 256> inv_sbox_{};
    std::array<uint8_t, 176> round_keys_{};

    void expand_key(const std::array<uint8_t, kKeySize>& key);
    void sub_bytes(uint8_t state[kBlockSize], const std::array<uint8_t, 256>& box) const;
    static void shift_rows(uint8_t state[kBlockSize]);
    static void inv_shift_rows(uint8_t state[kBlockSize]);
    void mix_columns(uint8_t state[kBlockSize]) const;
    void inv_mix_columns(uint8_t state[kBlockSize]) const;
    static void add_round_key(uint8_t state[kBlockSize], const uint8_t* rk);
};
