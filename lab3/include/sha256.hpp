#pragma once

#include <cstddef>
#include <cstdint>

// SHA-256 (FIPS 180-4).
class Sha256 {
public:
    static constexpr int kDigestLen = 32;

    Sha256();

    void update(const uint8_t* data, std::size_t len);
    void final(uint8_t out[kDigestLen]);

    static void hash(const uint8_t* data, std::size_t len, uint8_t out[kDigestLen]);

private:
    uint32_t h_[8];
    uint8_t block_[64];
    std::size_t block_len_;
    uint64_t total_len_;

    void compress(const uint8_t block[64]);
};
