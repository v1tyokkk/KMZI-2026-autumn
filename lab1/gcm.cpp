#include "gcm.hpp"

#include <cstring>

namespace {

void xor_n(uint8_t* dst, const uint8_t* src, std::size_t n) {
    for (std::size_t i = 0; i < n; ++i) {
        dst[i] ^= src[i];
    }
}

}

AesGcm::AesGcm(const std::array<uint8_t, Aes128::kKeySize>& key, const Gf8& gf) : aes_(key, gf) {
    const uint8_t zero[kBlockSize] = {};
    aes_.encrypt_block(zero, h_.data());
}

bool AesGcm::ct_equal(const uint8_t* a, const uint8_t* b, std::size_t n) {
    uint8_t acc = 0;
    for (std::size_t i = 0; i < n; ++i) {
        acc |= static_cast<uint8_t>(a[i] ^ b[i]);
    }
    return acc == 0;
}

void AesGcm::inc32(uint8_t block[kBlockSize]) {
    const uint32_t n = (static_cast<uint32_t>(block[12]) << 24) |
                       (static_cast<uint32_t>(block[13]) << 16) |
                       (static_cast<uint32_t>(block[14]) << 8) |
                       static_cast<uint32_t>(block[15]);
    const uint32_t inc = n + 1u;
    block[12] = static_cast<uint8_t>(inc >> 24);
    block[13] = static_cast<uint8_t>(inc >> 16);
    block[14] = static_cast<uint8_t>(inc >> 8);
    block[15] = static_cast<uint8_t>(inc);
}

void AesGcm::xor_block(uint8_t dst[kBlockSize], const uint8_t src[kBlockSize]) {
    for (int i = 0; i < kBlockSize; ++i) {
        dst[i] ^= src[i];
    }
}

void AesGcm::gf128_mul(uint8_t z[kBlockSize],
                       const uint8_t x[kBlockSize],
                       const uint8_t y[kBlockSize]) {
    uint8_t v[kBlockSize];
    std::memcpy(v, y, kBlockSize);
    std::memset(z, 0, kBlockSize);

    for (int i = 0; i < 128; ++i) {
        const uint8_t bit = static_cast<uint8_t>((x[i >> 3] >> (7 - (i & 7))) & 1);
        const uint8_t mask = static_cast<uint8_t>(-static_cast<int8_t>(bit));
        for (int j = 0; j < kBlockSize; ++j) {
            z[j] ^= static_cast<uint8_t>(v[j] & mask);
        }

        const uint8_t lsb = static_cast<uint8_t>(v[15] & 1);
        for (int j = 15; j > 0; --j) {
            v[j] = static_cast<uint8_t>((v[j] >> 1) | (v[j - 1] << 7));
        }
        v[0] = static_cast<uint8_t>(v[0] >> 1);
        const uint8_t rmask = static_cast<uint8_t>(-static_cast<int8_t>(lsb));
        v[0] ^= static_cast<uint8_t>(0xE1 & rmask);
    }
}

void AesGcm::ghash(const uint8_t* aad, std::size_t aad_len,
                   const uint8_t* ct, std::size_t ct_len,
                   uint8_t out[kBlockSize]) const {
    uint8_t y[kBlockSize] = {};
    uint8_t block[kBlockSize];

    auto absorb = [&](const uint8_t* data, std::size_t len) {
        std::size_t off = 0;
        while (off < len) {
            std::memset(block, 0, kBlockSize);
            const std::size_t n = (len - off < kBlockSize) ? (len - off) : kBlockSize;
            std::memcpy(block, data + off, n);
            xor_block(y, block);
            gf128_mul(block, y, h_.data());
            std::memcpy(y, block, kBlockSize);
            off += n;
        }
    };

    absorb(aad, aad_len);
    absorb(ct, ct_len);

    std::memset(block, 0, kBlockSize);
    const uint64_t aad_bits = static_cast<uint64_t>(aad_len) * 8ull;
    const uint64_t ct_bits = static_cast<uint64_t>(ct_len) * 8ull;
    for (int i = 0; i < 8; ++i) {
        block[i] = static_cast<uint8_t>(aad_bits >> (56 - 8 * i));
        block[8 + i] = static_cast<uint8_t>(ct_bits >> (56 - 8 * i));
    }
    xor_block(y, block);
    gf128_mul(out, y, h_.data());
}

void AesGcm::compute_j0(const uint8_t* iv, std::size_t iv_len, uint8_t j0[kBlockSize]) const {
    if (iv_len == 12) {
        std::memcpy(j0, iv, 12);
        j0[12] = 0;
        j0[13] = 0;
        j0[14] = 0;
        j0[15] = 1;
        return;
    }
    uint8_t y[kBlockSize] = {};
    uint8_t block[kBlockSize];
    std::size_t off = 0;
    while (off < iv_len) {
        std::memset(block, 0, kBlockSize);
        const std::size_t n = (iv_len - off < kBlockSize) ? (iv_len - off) : kBlockSize;
        std::memcpy(block, iv + off, n);
        xor_block(y, block);
        gf128_mul(block, y, h_.data());
        std::memcpy(y, block, kBlockSize);
        off += n;
    }
    std::memset(block, 0, kBlockSize);
    const uint64_t iv_bits = static_cast<uint64_t>(iv_len) * 8ull;
    for (int i = 0; i < 8; ++i) {
        block[8 + i] = static_cast<uint8_t>(iv_bits >> (56 - 8 * i));
    }
    xor_block(y, block);
    gf128_mul(j0, y, h_.data());
}

void AesGcm::ctr_crypt(const uint8_t j0[kBlockSize],
                       const uint8_t* in, std::size_t len, uint8_t* out) const {
    uint8_t ctr[kBlockSize];
    std::memcpy(ctr, j0, kBlockSize);
    inc32(ctr);

    uint8_t ks[kBlockSize];
    std::size_t off = 0;
    while (off < len) {
        aes_.encrypt_block(ctr, ks);
        const std::size_t n = (len - off < kBlockSize) ? (len - off) : kBlockSize;
        std::memcpy(out + off, in + off, n);
        xor_n(out + off, ks, n);
        inc32(ctr);
        off += n;
    }
}

void AesGcm::encrypt(const uint8_t* iv, std::size_t iv_len,
                     const uint8_t* aad, std::size_t aad_len,
                     const uint8_t* pt, std::size_t pt_len,
                     uint8_t* ct,
                     std::array<uint8_t, kTagSize>& tag) const {
    uint8_t j0[kBlockSize];
    compute_j0(iv, iv_len, j0);

    if (pt_len > 0) {
        ctr_crypt(j0, pt, pt_len, ct);
    }

    uint8_t s[kBlockSize];
    ghash(aad, aad_len, ct, pt_len, s);

    uint8_t enc0[kBlockSize];
    aes_.encrypt_block(j0, enc0);
    for (int i = 0; i < kTagSize; ++i) {
        tag[static_cast<std::size_t>(i)] = static_cast<uint8_t>(s[i] ^ enc0[i]);
    }
}

bool AesGcm::decrypt(const uint8_t* iv, std::size_t iv_len,
                     const uint8_t* aad, std::size_t aad_len,
                     const uint8_t* ct, std::size_t ct_len,
                     const uint8_t tag[kTagSize],
                     uint8_t* pt) const {
    uint8_t j0[kBlockSize];
    compute_j0(iv, iv_len, j0);

    uint8_t s[kBlockSize];
    ghash(aad, aad_len, ct, ct_len, s);

    uint8_t enc0[kBlockSize];
    aes_.encrypt_block(j0, enc0);
    uint8_t expected[kTagSize];
    for (int i = 0; i < kTagSize; ++i) {
        expected[i] = static_cast<uint8_t>(s[i] ^ enc0[i]);
    }

    if (ct_len > 0) {
        ctr_crypt(j0, ct, ct_len, pt);
    }

    return ct_equal(expected, tag, kTagSize);
}
