#include "aes128.hpp"

#include <cstring>

Aes128::Aes128(const std::array<uint8_t, kKeySize>& key, const Gf8& gf) : gf_(gf) {
    gf_.generate_sbox(sbox_);
    gf_.generate_inv_sbox(sbox_, inv_sbox_);
    expand_key(key);
}

void Aes128::expand_key(const std::array<uint8_t, kKeySize>& key) {
    std::memcpy(round_keys_.data(), key.data(), kKeySize);

    uint8_t rcon = 0x01;
    std::size_t bytes = kKeySize;
    while (bytes < round_keys_.size()) {
        uint8_t t0 = round_keys_[bytes - 4];
        uint8_t t1 = round_keys_[bytes - 3];
        uint8_t t2 = round_keys_[bytes - 2];
        uint8_t t3 = round_keys_[bytes - 1];

        if (bytes % 16 == 0) {
            const uint8_t rot0 = t1, rot1 = t2, rot2 = t3, rot3 = t0;
            t0 = static_cast<uint8_t>(sbox_[rot0] ^ rcon);
            t1 = sbox_[rot1];
            t2 = sbox_[rot2];
            t3 = sbox_[rot3];
            rcon = gf_.xtime(rcon);
        }

        round_keys_[bytes] = static_cast<uint8_t>(round_keys_[bytes - 16] ^ t0);
        round_keys_[bytes + 1] = static_cast<uint8_t>(round_keys_[bytes - 15] ^ t1);
        round_keys_[bytes + 2] = static_cast<uint8_t>(round_keys_[bytes - 14] ^ t2);
        round_keys_[bytes + 3] = static_cast<uint8_t>(round_keys_[bytes - 13] ^ t3);
        bytes += 4;
    }
}

void Aes128::sub_bytes(uint8_t state[kBlockSize], const std::array<uint8_t, 256>& box) const {
    for (int i = 0; i < kBlockSize; ++i) {
        state[i] = box[state[i]];
    }
}

void Aes128::shift_rows(uint8_t state[kBlockSize]) {
    uint8_t t = state[1];
    state[1] = state[5];
    state[5] = state[9];
    state[9] = state[13];
    state[13] = t;

    t = state[2];
    uint8_t t2 = state[6];
    state[2] = state[10];
    state[6] = state[14];
    state[10] = t;
    state[14] = t2;

    t = state[15];
    state[15] = state[11];
    state[11] = state[7];
    state[7] = state[3];
    state[3] = t;
}

void Aes128::inv_shift_rows(uint8_t state[kBlockSize]) {
    uint8_t t = state[13];
    state[13] = state[9];
    state[9] = state[5];
    state[5] = state[1];
    state[1] = t;

    t = state[2];
    uint8_t t2 = state[6];
    state[2] = state[10];
    state[6] = state[14];
    state[10] = t;
    state[14] = t2;

    t = state[3];
    state[3] = state[7];
    state[7] = state[11];
    state[11] = state[15];
    state[15] = t;
}

void Aes128::mix_columns(uint8_t state[kBlockSize]) const {
    for (int c = 0; c < 4; ++c) {
        const uint8_t a = state[4 * c];
        const uint8_t b = state[4 * c + 1];
        const uint8_t c0 = state[4 * c + 2];
        const uint8_t d = state[4 * c + 3];
        state[4 * c] = static_cast<uint8_t>(gf_.multiply(a, 0x02) ^ gf_.multiply(b, 0x03) ^ c0 ^ d);
        state[4 * c + 1] = static_cast<uint8_t>(a ^ gf_.multiply(b, 0x02) ^ gf_.multiply(c0, 0x03) ^ d);
        state[4 * c + 2] = static_cast<uint8_t>(a ^ b ^ gf_.multiply(c0, 0x02) ^ gf_.multiply(d, 0x03));
        state[4 * c + 3] = static_cast<uint8_t>(gf_.multiply(a, 0x03) ^ b ^ c0 ^ gf_.multiply(d, 0x02));
    }
}

void Aes128::inv_mix_columns(uint8_t state[kBlockSize]) const {
    for (int c = 0; c < 4; ++c) {
        const uint8_t a = state[4 * c];
        const uint8_t b = state[4 * c + 1];
        const uint8_t c0 = state[4 * c + 2];
        const uint8_t d = state[4 * c + 3];
        state[4 * c] = static_cast<uint8_t>(gf_.multiply(a, 0x0E) ^ gf_.multiply(b, 0x0B) ^
                                           gf_.multiply(c0, 0x0D) ^ gf_.multiply(d, 0x09));
        state[4 * c + 1] = static_cast<uint8_t>(gf_.multiply(a, 0x09) ^ gf_.multiply(b, 0x0E) ^
                                               gf_.multiply(c0, 0x0B) ^ gf_.multiply(d, 0x0D));
        state[4 * c + 2] = static_cast<uint8_t>(gf_.multiply(a, 0x0D) ^ gf_.multiply(b, 0x09) ^
                                               gf_.multiply(c0, 0x0E) ^ gf_.multiply(d, 0x0B));
        state[4 * c + 3] = static_cast<uint8_t>(gf_.multiply(a, 0x0B) ^ gf_.multiply(b, 0x0D) ^
                                               gf_.multiply(c0, 0x09) ^ gf_.multiply(d, 0x0E));
    }
}

void Aes128::add_round_key(uint8_t state[kBlockSize], const uint8_t* rk) {
    for (int i = 0; i < kBlockSize; ++i) {
        state[i] ^= rk[i];
    }
}

void Aes128::encrypt_block(const uint8_t in[kBlockSize], uint8_t out[kBlockSize]) const {
    uint8_t state[kBlockSize];
    std::memcpy(state, in, kBlockSize);

    add_round_key(state, round_keys_.data());
    for (int r = 1; r < kRounds; ++r) {
        sub_bytes(state, sbox_);
        shift_rows(state);
        mix_columns(state);
        add_round_key(state, round_keys_.data() + static_cast<std::size_t>(r) * kBlockSize);
    }
    sub_bytes(state, sbox_);
    shift_rows(state);
    add_round_key(state, round_keys_.data() + static_cast<std::size_t>(kRounds) * kBlockSize);

    std::memcpy(out, state, kBlockSize);
}

void Aes128::decrypt_block(const uint8_t in[kBlockSize], uint8_t out[kBlockSize]) const {
    uint8_t state[kBlockSize];
    std::memcpy(state, in, kBlockSize);

    add_round_key(state, round_keys_.data() + static_cast<std::size_t>(kRounds) * kBlockSize);
    inv_shift_rows(state);
    sub_bytes(state, inv_sbox_);
    for (int r = kRounds - 1; r > 0; --r) {
        add_round_key(state, round_keys_.data() + static_cast<std::size_t>(r) * kBlockSize);
        inv_mix_columns(state);
        inv_shift_rows(state);
        sub_bytes(state, inv_sbox_);
    }
    add_round_key(state, round_keys_.data());

    std::memcpy(out, state, kBlockSize);
}
