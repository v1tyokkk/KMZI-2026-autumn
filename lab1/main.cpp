#include "gf8.hpp"
#include "aes128.hpp"
#include "gcm.hpp"

#include <array>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace {

std::string hex_of(const uint8_t* data, std::size_t n) {
    static const char* kDigits = "0123456789abcdef";
    std::string out;
    out.resize(n * 2);
    for (std::size_t i = 0; i < n; ++i) {
        out[2 * i] = kDigits[data[i] >> 4];
        out[2 * i + 1] = kDigits[data[i] & 0x0F];
    }
    return out;
}

std::string hex_of(const std::array<uint8_t, 16>& a) { return hex_of(a.data(), a.size()); }

int g_failed = 0;

void expect_eq(const char* name, const std::string& got, const std::string& want) {
    if (got == want) {
        std::printf("  [OK] %s\n", name);
        return;
    }
    std::printf("  [FAIL] %s\n    got  %s\n    want %s\n", name, got.c_str(), want.c_str());
    ++g_failed;
}

void expect_true(const char* name, bool ok) {
    if (ok) {
        std::printf("  [OK] %s\n", name);
        return;
    }
    std::printf("  [FAIL] %s\n", name);
    ++g_failed;
}

uint8_t schoolbook_mul(uint8_t a, uint8_t b, uint16_t poly) {
    uint8_t r = 0;
    const uint8_t p_byte = static_cast<uint8_t>(poly & 0xFF);
    for (int i = 0; i < 8; ++i) {
        if (b & 1) {
            r ^= a;
        }
        const uint8_t hi = static_cast<uint8_t>(a & 0x80);
        a = static_cast<uint8_t>(a << 1);
        if (hi) {
            a ^= p_byte;
        }
        b = static_cast<uint8_t>(b >> 1);
    }
    return r;
}

void parse_hex(const char* s, std::vector<uint8_t>& out) {
    out.clear();
    std::size_t n = std::strlen(s);
    if (n % 2 != 0) {
        return;
    }
    auto nibble = [](char c) -> int {
        if (c >= '0' && c <= '9') return c - '0';
        if (c >= 'a' && c <= 'f') return c - 'a' + 10;
        if (c >= 'A' && c <= 'F') return c - 'A' + 10;
        return -1;
    };
    out.reserve(n / 2);
    for (std::size_t i = 0; i < n; i += 2) {
        const int hi = nibble(s[i]);
        const int lo = nibble(s[i + 1]);
        if (hi < 0 || lo < 0) {
            out.clear();
            return;
        }
        out.push_back(static_cast<uint8_t>((hi << 4) | lo));
    }
}

void test_gf8() {
    std::printf("\n== Калькулятор GF(2^8), p(x)=0x14D ==\n");
    const Gf8 gf(Gf8::kVariantPoly);

    expect_eq("add(0x57, 0x13)",
              hex_of(std::array<uint8_t, 1>{Gf8::add(0x57, 0x13)}.data(), 1),
              "44");
    expect_eq("multiply Karatsuba 0x57*0x13",
              hex_of(std::array<uint8_t, 1>{gf.multiply(0x57, 0x13)}.data(), 1),
              "bd");
    expect_eq("xtime(0x80) без ветвлений",
              hex_of(std::array<uint8_t, 1>{gf.xtime(0x80)}.data(), 1),
              "4d");
    expect_eq("inverse(0x00) -> 0x00",
              hex_of(std::array<uint8_t, 1>{gf.inverse(0x00)}.data(), 1),
              "00");
    expect_eq("inverse(0x53)",
              hex_of(std::array<uint8_t, 1>{gf.inverse(0x53)}.data(), 1),
              "04");
    expect_eq("0x53 * inverse(0x53)",
              hex_of(std::array<uint8_t, 1>{gf.multiply(0x53, gf.inverse(0x53))}.data(), 1),
              "01");

    int mismatches = 0;
    int inv_fail = 0;
    for (int a = 0; a < 256; ++a) {
        const uint8_t ua = static_cast<uint8_t>(a);
        if (gf.multiply(ua, gf.inverse(ua)) != (a == 0 ? 0 : 1)) {
            ++inv_fail;
        }
        for (int b = 0; b < 256; ++b) {
            const uint8_t ub = static_cast<uint8_t>(b);
            if (gf.multiply(ua, ub) != schoolbook_mul(ua, ub, 0x14D)) {
                ++mismatches;
            }
        }
    }
    expect_true("Карацуба == shift-and-XOR на всех 65536 парах", mismatches == 0);
    expect_true("a * a^{-1} = 1 (и 0*0^{-1}=0)", inv_fail == 0);

    const Gf8 aes_field(0x11B);
    expect_eq("AES-поле: 0x57*0x83 = 0xC1",
              hex_of(std::array<uint8_t, 1>{aes_field.multiply(0x57, 0x83)}.data(), 1),
              "c1");

    std::array<uint8_t, 256> sbox{};
    aes_field.generate_sbox(sbox);
    expect_eq("S-box AES S[0x00]=0x63", hex_of(&sbox[0x00], 1), "63");
    expect_eq("S-box AES S[0x53]=0xED", hex_of(&sbox[0x53], 1), "ed");

    std::array<uint8_t, 256> sbox14{};
    gf.generate_sbox(sbox14);
    expect_eq("S-box варианта 11, первые 16 байт",
              hex_of(sbox14.data(), 16),
              "637c4d5774b179ca4904abfc6eefb72d");
}

void test_aes() {
    std::printf("\n== AES-128 ==\n");
    const std::array<uint8_t, 16> key = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
                                         0x08, 0x09, 0x0A, 0x0B, 0x0C, 0x0D, 0x0E, 0x0F};
    const uint8_t pt[16] = {0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77,
                            0x88, 0x99, 0xAA, 0xBB, 0xCC, 0xDD, 0xEE, 0xFF};

    {
        Aes128 aes(key, Gf8(0x11B));
        uint8_t ct[16] = {};
        uint8_t back[16] = {};
        aes.encrypt_block(pt, ct);
        aes.decrypt_block(ct, back);
        expect_eq("NIST FIPS-197 ciphertext", hex_of(ct, 16),
                  "69c4e0d86a7b0430d8cdb78070b4c55a");
        expect_eq("NIST FIPS-197 roundtrip", hex_of(back, 16), hex_of(pt, 16));
    }

    {
        Aes128 aes(key, Gf8(Gf8::kVariantPoly));
        uint8_t ct[16] = {};
        uint8_t back[16] = {};
        aes.encrypt_block(pt, ct);
        aes.decrypt_block(ct, back);
        expect_eq("AES-128 вариант 11 ciphertext", hex_of(ct, 16),
                  "5a16d072429581b890cb18718121b495");
        expect_eq("AES-128 вариант 11 roundtrip", hex_of(back, 16), hex_of(pt, 16));
    }
}

void test_gcm() {
    std::printf("\n== AES-GCM AEAD ==\n");

    {
        const std::array<uint8_t, 16> key{};
        const uint8_t iv[12] = {};
        AesGcm gcm(key, Gf8(0x11B));
        std::array<uint8_t, 16> tag{};
        gcm.encrypt(iv, 12, nullptr, 0, nullptr, 0, nullptr, tag);
        expect_eq("NIST GCM Test Case 1 tag", hex_of(tag),
                  "58e2fccefa7e3061367f1d57a4e7455a");
    }

    {
        const std::array<uint8_t, 16> key{};
        const uint8_t iv[12] = {};
        const uint8_t pt[16] = {};
        uint8_t ct[16] = {};
        std::array<uint8_t, 16> tag{};
        AesGcm gcm(key, Gf8(0x11B));
        gcm.encrypt(iv, 12, nullptr, 0, pt, 16, ct, tag);
        expect_eq("NIST GCM Test Case 2 ciphertext", hex_of(ct, 16),
                  "0388dace60b6a392f328c2b971b2fe78");
        expect_eq("NIST GCM Test Case 2 tag", hex_of(tag),
                  "ab6e47d42cec13bdf53a67b21257bddf");
    }

    const std::array<uint8_t, 16> key = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
                                         0x08, 0x09, 0x0A, 0x0B, 0x0C, 0x0D, 0x0E, 0x0F};
    const uint8_t iv[12] = {0xCA, 0xFE, 0xBA, 0xBE, 0xFA, 0xCE, 0xDB, 0xAD,
                            0xDE, 0xCA, 0xF8, 0x88};
    const uint8_t aad[4] = {0xFE, 0xED, 0xFA, 0xCE};
    const uint8_t pt[16] = {0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77,
                            0x88, 0x99, 0xAA, 0xBB, 0xCC, 0xDD, 0xEE, 0xFF};
    uint8_t ct[16] = {};
    std::array<uint8_t, 16> tag{};
    AesGcm gcm(key, Gf8(Gf8::kVariantPoly));
    gcm.encrypt(iv, 12, aad, sizeof(aad), pt, sizeof(pt), ct, tag);
    expect_eq("GCM вариант 11 ciphertext", hex_of(ct, 16), "2e7ec25d36941f64c24a248e88dadf3c");
    expect_eq("GCM вариант 11 tag", hex_of(tag), "369d8a24c7c8e3c8f42ee11f297774e7");

    uint8_t back[16] = {};
    const bool ok = gcm.decrypt(iv, 12, aad, sizeof(aad), ct, sizeof(ct), tag.data(), back);
    expect_true("GCM decrypt: тег принят", ok);
    expect_eq("GCM decrypt: plaintext восстановлен", hex_of(back, 16), hex_of(pt, 16));

    uint8_t bad_tag[16];
    std::memcpy(bad_tag, tag.data(), 16);
    bad_tag[0] ^= 0x01;
    uint8_t dummy[16] = {};
    const bool rejected = !gcm.decrypt(iv, 12, aad, sizeof(aad), ct, sizeof(ct), bad_tag, dummy);
    expect_true("GCM decrypt: испорченный тег отвергнут", rejected);

    uint8_t a[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    uint8_t b[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    expect_true("ct_equal: одинаковые теги", AesGcm::ct_equal(a, b, 16));
    b[15] ^= 0xFF;
    expect_true("ct_equal: различие в последнем байте", !AesGcm::ct_equal(a, b, 16));
    b[15] ^= 0xFF;
    b[0] ^= 0xFF;
    expect_true("ct_equal: различие в первом байте (нет early-exit)", !AesGcm::ct_equal(a, b, 16));
}

void demo_variant11() {
    std::printf("\n== Демонстрация варианта 11 ==\n");
    std::printf("Алгоритм: AES-128-GCM\n");
    std::printf("p(x) = x^8 + x^6 + x^3 + x^2 + 1  (0x14D)\n");
    std::printf("f(x) = x^{128} + x^7 + x^2 + x + 1\n");
    std::printf("Constant-time: умножение Карацубы в GF(2^8)\n\n");

    const Gf8 gf(Gf8::kVariantPoly);
    std::printf("Пример Карацубы: 0x57 * 0x13 = 0x%02X\n", gf.multiply(0x57, 0x13));
    std::printf("A = 0x57 = A1||A0 = 0x5||0x7,  B = 0x13 = 0x1||0x3\n");
    std::printf("P0=(A0*B0), P2=(A1*B1), P1=((A0+A1)*(B0+B1)), затем редукция mod p(x)\n");
    std::printf("inverse(0x00)=0x%02X, inverse(0x53)=0x%02X\n", gf.inverse(0x00), gf.inverse(0x53));

    const std::array<uint8_t, 16> key = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
                                         0x08, 0x09, 0x0A, 0x0B, 0x0C, 0x0D, 0x0E, 0x0F};
    const uint8_t iv[12] = {0xCA, 0xFE, 0xBA, 0xBE, 0xFA, 0xCE, 0xDB, 0xAD,
                            0xDE, 0xCA, 0xF8, 0x88};
    const uint8_t aad[4] = {0xFE, 0xED, 0xFA, 0xCE};
    const uint8_t pt[16] = {0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77,
                            0x88, 0x99, 0xAA, 0xBB, 0xCC, 0xDD, 0xEE, 0xFF};

    AesGcm gcm(key, gf);
    uint8_t ct[16] = {};
    std::array<uint8_t, 16> tag{};
    gcm.encrypt(iv, 12, aad, sizeof(aad), pt, sizeof(pt), ct, tag);

    std::printf("\nKey : %s\n", hex_of(key).c_str());
    std::printf("IV  : %s\n", hex_of(iv, 12).c_str());
    std::printf("AAD : %s\n", hex_of(aad, 4).c_str());
    std::printf("PT  : %s\n", hex_of(pt, 16).c_str());
    std::printf("CT  : %s\n", hex_of(ct, 16).c_str());
    std::printf("Tag : %s\n", hex_of(tag).c_str());
}

int run_cli(int argc, char** argv) {
    if (argc < 5) {
        std::printf("Использование: %s encrypt <key16hex> <iv12hex> <pthex> [aadhex]\n", argv[0]);
        return 1;
    }

    std::vector<uint8_t> key, iv, pt, aad;
    parse_hex(argv[2], key);
    parse_hex(argv[3], iv);
    parse_hex(argv[4], pt);
    if (argc >= 6) {
        parse_hex(argv[5], aad);
    }
    if (key.size() != 16 || iv.size() != 12) {
        std::printf("Ключ должен быть 16 байт, IV — 12 байт (hex).\n");
        return 1;
    }

    std::array<uint8_t, 16> key_arr{};
    std::memcpy(key_arr.data(), key.data(), 16);
    AesGcm gcm(key_arr, Gf8(Gf8::kVariantPoly));
    std::vector<uint8_t> ct(pt.size());
    std::array<uint8_t, 16> tag{};
    gcm.encrypt(iv.data(), iv.size(),
                aad.empty() ? nullptr : aad.data(), aad.size(),
                pt.empty() ? nullptr : pt.data(), pt.size(),
                ct.empty() ? nullptr : ct.data(), tag);
    std::printf("CT  %s\n", hex_of(ct.empty() ? nullptr : ct.data(), ct.size()).c_str());
    std::printf("Tag %s\n", hex_of(tag).c_str());
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc > 1 && std::strcmp(argv[1], "encrypt") == 0) {
        return run_cli(argc, argv);
    }

    std::printf("Лабораторная работа №1, вариант 11\n");
    std::printf("AES-128-GCM, p(x)=0x14D, умножение Карацубы в GF(2^8)\n");

    test_gf8();
    test_aes();
    test_gcm();
    demo_variant11();

    std::printf("\n== Итог: %s (%d ошибок) ==\n",
                g_failed == 0 ? "все проверки пройдены" : "есть ошибки", g_failed);
    return g_failed == 0 ? 0 : 1;
}
