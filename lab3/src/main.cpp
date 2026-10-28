#include "ct.hpp"
#include "ecdsa.hpp"
#include "sha256.hpp"

#include <cstdio>
#include <cstring>
#include <string>

namespace {

int g_failed = 0;

std::string hash_hex(const uint8_t* data, std::size_t len) {
    uint8_t digest[Sha256::kDigestLen];
    Sha256::hash(data, len, digest);
    static const char* kDigits = "0123456789abcdef";
    std::string out;
    out.reserve(64);
    for (uint8_t b : digest) {
        out.push_back(kDigits[b >> 4]);
        out.push_back(kDigits[b & 0x0F]);
    }
    return out;
}

bool parse_hex(const char* hex, uint8_t* out, std::size_t len) {
    if (hex == nullptr || std::strlen(hex) != len * 2) {
        return false;
    }
    for (std::size_t i = 0; i < len; ++i) {
        int hi = 0;
        int lo = 0;
        for (int k = 0; k < 2; ++k) {
            const char c = hex[2 * i + k];
            int v;
            if (c >= '0' && c <= '9') {
                v = c - '0';
            } else if (c >= 'a' && c <= 'f') {
                v = c - 'a' + 10;
            } else if (c >= 'A' && c <= 'F') {
                v = c - 'A' + 10;
            } else {
                return false;
            }
            (k == 0 ? hi : lo) = v;
        }
        out[i] = static_cast<uint8_t>((hi << 4) | lo);
    }
    return true;
}

bool parse_bytes32(uint8_t out[32], const char* hex) { return parse_hex(hex, out, 32); }

bool parse_scalar(U256& out, const char* hex) { return u256::from_hex(out, hex); }

U256 sub_one(const U256& v) {
    U256 r = v;
    uint64_t borrow = 1;
    for (int i = 0; i < 4 && borrow != 0; ++i) {
        const uint64_t old = r.w[i];
        r.w[i] = old - borrow;
        borrow = (old < borrow) ? 1u : 0u;
    }
    return r;
}

void expect_true(const char* name, bool ok) {
    if (ok) {
        std::printf("  [OK]   %s\n", name);
    } else {
        std::printf("  [FAIL] %s\n", name);
        ++g_failed;
    }
}

void expect_hex(const char* name, const U256& got, const char* want) {
    const std::string g = u256::to_hex(got);
    if (g == want) {
        std::printf("  [OK]   %s\n", name);
    } else {
        std::printf("  [FAIL] %s\n         got  %s\n         want %s\n", name, g.c_str(), want);
        ++g_failed;
    }
}

void print_point(const char* label, const PointAff& p) {
    if (p.infinity) {
        std::printf("%s\n  infinity\n", label);
        return;
    }
    std::printf("%s\n  %s\n  %s\n", label, u256::to_hex(p.x).c_str(), u256::to_hex(p.y).c_str());
}

int run_demo() {
    std::printf("Лабораторная работа 3, вариант 11\n");
    std::printf("ECDSA, secp256k1, проективные координаты Якоби, Joye's Double-Add\n");

    const Curve& curve = Curve::secp256k1();
    const Mod& fp = curve.field();
    const Mod& fn = curve.order();

    std::printf("\n== SHA-256 и арифметика GF(p) ==\n");
    const uint8_t abc[] = {'a', 'b', 'c'};
    const std::string abc_hash = hash_hex(abc, 3);
    expect_true("SHA-256(abc)",
                abc_hash == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    expect_hex("2 + 3 = 5", fp.from_mont(fp.add(fp.from_small(2), fp.from_small(3))), "0000000000000000000000000000000000000000000000000000000000000005");
    expect_hex("7 * 8 = 56", fp.from_mont(fp.mul(fp.from_small(7), fp.from_small(8))), "0000000000000000000000000000000000000000000000000000000000000038");
    expect_hex("7 * 7^-1 = 1", fp.from_mont(fp.mul(fp.from_small(7), fp.inv(fp.from_small(7)))), "0000000000000000000000000000000000000000000000000000000000000001");

    std::printf("\n== Кривая и скалярное умножение ==\n");
    const PointAff g = ec_to_affine(curve, curve.generator());
    expect_true("G лежит на кривой", curve.on_curve(g));
    U256 two;
    u256::set_word(two, 2);
    const PointAff two_g = ec_to_affine(curve, ct::secure_mul(curve, curve.generator(), two));
    expect_hex("2G.x", two_g.x, "c6047f9441ed7d6d3045406e95c07cd85c778e4b8cef3ca7abac09b95c709ee5");
    expect_hex("2G.y", two_g.y, "1ae168fea63dc339a3c58419466ceaeef7f632653266d0e1236431a950cfe52a");
    expect_true("q * G = O", ec_to_affine(curve, ct::secure_mul(curve, curve.generator(), fn.modulus())).infinity);
    const PointAff neg_g = ec_to_affine(curve, ct::secure_mul(curve, curve.generator(), sub_one(fn.modulus())));
    expect_true("(q - 1) * G = -G", !neg_g.infinity && u256::compare(neg_g.x, g.x) == 0 &&
                                        u256::compare(neg_g.y, fp.from_mont(fp.sub(fp.zero(), fp.to_mont(g.y)))) == 0);

    bool multiples_ok = true;
    PointJac acc = ec_infinity(curve);
    for (uint64_t k = 1; k <= 8; ++k) {
        U256 s;
        u256::set_word(s, k);
        acc = ec_add(curve, acc, curve.generator());
        multiples_ok = multiples_ok && ec_equal_affine(ec_to_affine(curve, acc), ec_to_affine(curve, ct::secure_mul(curve, curve.generator(), s)));
    }
    expect_true("1G..8G совпадают с последовательным сложением", multiples_ok);

    std::printf("\n== ECDH ==\n");
    U256 d1;
    U256 d2;
    u256::set_word(d1, 1);
    u256::set_word(d2, 2);
    PointAff q1;
    PointAff q2;
    U256 s1;
    U256 s2;
    const bool fixed_ecdh = EcProtocols::public_key(curve, d1, q1) && EcProtocols::public_key(curve, d2, q2) &&
                            EcProtocols::ecdh_shared(curve, d1, q2, s1) &&
                            EcProtocols::ecdh_shared(curve, d2, q1, s2);
    expect_true("dA=1, dB=2", fixed_ecdh && u256::compare(s1, s2) == 0 && u256::compare(s1, two_g.x) == 0);

    EcKeyPair alice;
    EcKeyPair bob;
    const bool keygen = EcProtocols::generate_key(curve, alice) && EcProtocols::generate_key(curve, bob);
    const bool both = keygen && EcProtocols::ecdh_shared(curve, alice.d, bob.q, s1) &&
                      EcProtocols::ecdh_shared(curve, bob.d, alice.q, s2);
    expect_true("случайные ключи, общий секрет совпал", both && u256::compare(s1, s2) == 0);
    print_point("Алиса, открытый ключ:", alice.q);
    print_point("Боб, открытый ключ:", bob.q);
    std::printf("Общий секрет (x):\n  %s\n", u256::to_hex(s1).c_str());

    std::printf("\n== ECDSA ==\n");
    const char* message_hex = "4b4d5a492d323032363a206c6162332076617269616e74203131";
    uint8_t message[64];
    const std::size_t message_len = std::strlen(message_hex) / 2;
    expect_true("тестовое сообщение в hex", parse_hex(message_hex, message, message_len));
    uint8_t hash[32];
    Sha256::hash(message, message_len, hash);
    expect_true("Hash(M) по эталону",
                hash_hex(message, message_len) ==
                    "9c1824c56d9512bdfbeaacdfde79a0c4cb051539a4a8a8febb20d9713e3bad93");
    std::printf("Hash(M):\n  ");
    for (uint8_t b : hash) {
        std::printf("%02x", b);
    }
    std::printf("\n");

    U256 r;
    U256 s;
    const bool signed_ok = EcProtocols::sign(curve, alice.d, hash, r, s);
    expect_true("подпись сформирована", signed_ok);
    expect_true("подпись верна", signed_ok && EcProtocols::verify(curve, alice.q, hash, r, s));
    U256 bad_s = s;
    bad_s.w[0] ^= 0xFF;
    expect_true("искажённый байт s отвергнут", signed_ok && !EcProtocols::verify(curve, alice.q, hash, r, bad_s));
    uint8_t bad_hash[32];
    std::memcpy(bad_hash, hash, sizeof(hash));
    bad_hash[0] ^= 1;
    expect_true("искажённый хэш отвергнут", signed_ok && !EcProtocols::verify(curve, alice.q, bad_hash, r, s));
    std::printf("r:\n  %s\ns:\n  %s\n", u256::to_hex(r).c_str(), u256::to_hex(s).c_str());

    U256 fixed_k;
    U256 fixed_r;
    U256 fixed_s;
    u256::set_word(fixed_k, 2);
    u256::zero(fixed_r);
    u256::zero(fixed_s);
    const bool fixed_sign = EcProtocols::sign_with_k(curve, d1, hash, fixed_k, fixed_r, fixed_s);
    expect_hex("ECDSA(d=1, k=2): r", fixed_r,
               "c6047f9441ed7d6d3045406e95c07cd85c778e4b8cef3ca7abac09b95c709ee5");
    expect_hex("ECDSA(d=1, k=2): s", fixed_s,
               "b10e522cd7c148159617f6a73a1d0ece93be51c298cbf2d3336671954d56263c");
    expect_true("фиксированная подпись верна", fixed_sign && EcProtocols::verify(curve, q1, hash, fixed_r, fixed_s));

    std::printf("\n%s, провалов: %d\n", g_failed == 0 ? "Все проверки пройдены" : "Есть ошибки", g_failed);
    return g_failed == 0 ? 0 : 1;
}

void print_usage(const char* argv0) {
    std::printf(
        "Использование:\n"
        "  %s demo\n"
        "  %s pub <d>\n"
        "  %s ecdh <dA> <dB>\n"
        "  %s shared <d> <Qx> <Qy>\n"
        "  %s sign <d> <hash32>\n"
        "  %s verify <Qx> <Qy> <hash32> <r> <s>\n"
        "  %s hash <hex-message>\n"
        "Все значения — сырые hex-строки; hash32 — 64 hex-символа.\n",
        argv0, argv0, argv0, argv0, argv0, argv0, argv0);
}

int run_cli(int argc, char** argv) {
    const Curve& curve = Curve::secp256k1();
    const std::string cmd = argv[1];

    if (cmd == "pub" && argc == 3) {
        U256 d;
        PointAff q;
        if (!parse_scalar(d, argv[2]) || !EcProtocols::public_key(curve, d, q)) {
            return 2;
        }
        std::printf("%s\n%s\n", u256::to_hex(q.x).c_str(), u256::to_hex(q.y).c_str());
        return 0;
    }

    if (cmd == "ecdh" && argc == 4) {
        U256 da;
        U256 db;
        PointAff qa;
        PointAff qb;
        U256 shared_a;
        U256 shared_b;
        if (!parse_scalar(da, argv[2]) || !parse_scalar(db, argv[3]) ||
            !EcProtocols::public_key(curve, da, qa) || !EcProtocols::public_key(curve, db, qb) ||
            !EcProtocols::ecdh_shared(curve, da, qb, shared_a) ||
            !EcProtocols::ecdh_shared(curve, db, qa, shared_b) || u256::compare(shared_a, shared_b) != 0) {
            return 2;
        }
        std::printf("%s\n%s\n%s\n%s\n%s\n", u256::to_hex(qa.x).c_str(), u256::to_hex(qa.y).c_str(),
                    u256::to_hex(qb.x).c_str(), u256::to_hex(qb.y).c_str(),
                    u256::to_hex(shared_a).c_str());
        return 0;
    }

    if (cmd == "shared" && argc == 5) {
        U256 d;
        PointAff peer;
        U256 shared;
        if (!parse_scalar(d, argv[2]) || !parse_scalar(peer.x, argv[3]) ||
            !parse_scalar(peer.y, argv[4])) {
            return 2;
        }
        peer.infinity = false;
        if (!EcProtocols::ecdh_shared(curve, d, peer, shared)) {
            return 2;
        }
        std::printf("%s\n", u256::to_hex(shared).c_str());
        return 0;
    }

    if (cmd == "sign" && argc == 4) {
        U256 d;
        U256 r;
        U256 s;
        uint8_t hash[32];
        if (!parse_scalar(d, argv[2]) || !parse_bytes32(hash, argv[3]) ||
            !EcProtocols::sign(curve, d, hash, r, s)) {
            return 2;
        }
        std::printf("%s\n%s\n", u256::to_hex(r).c_str(), u256::to_hex(s).c_str());
        return 0;
    }

    if (cmd == "verify" && argc == 7) {
        PointAff q;
        U256 r;
        U256 s;
        uint8_t hash[32];
        if (!parse_scalar(q.x, argv[2]) || !parse_scalar(q.y, argv[3]) ||
            !parse_bytes32(hash, argv[4]) || !parse_scalar(r, argv[5]) ||
            !parse_scalar(s, argv[6])) {
            return 2;
        }
        q.infinity = false;
        std::printf("%s\n", EcProtocols::verify(curve, q, hash, r, s) ? "valid" : "invalid");
        return 0;
    }

    if (cmd == "hash" && argc == 3) {
        const std::size_t len = std::strlen(argv[2]);
        if (len % 2 != 0) {
            return 2;
        }
        std::string data(len / 2, '\0');
        if (!parse_hex(argv[2], reinterpret_cast<uint8_t*>(&data[0]), data.size())) {
            return 2;
        }
        std::printf("%s\n", hash_hex(reinterpret_cast<const uint8_t*>(data.data()), data.size()).c_str());
        return 0;
    }

    print_usage(argv[0]);
    return 2;
}

}

int main(int argc, char** argv) {
    if (argc >= 2 && std::strcmp(argv[1], "demo") != 0) {
        return run_cli(argc, argv);
    }
    return run_demo();
}
