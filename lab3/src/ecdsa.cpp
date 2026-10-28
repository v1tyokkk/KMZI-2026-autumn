#include "ecdsa.hpp"

#include "ct.hpp"

#include <fstream>

namespace {

bool scalar_valid(const Mod& fn, const U256& d) {
    return !u256::is_zero(d) && u256::compare(d, fn.modulus()) < 0;
}

}

U256 EcProtocols::hash_scalar(const Mod& fn, const uint8_t hash[32]) {
    U256 h;
    u256::from_bytes(h, hash, 32);
    return fn.reduce(h);
}

bool EcProtocols::random_scalar(const Mod& fn, U256& out) {
    std::ifstream rng("/dev/urandom", std::ios::binary);
    if (!rng) {
        return false;
    }
    uint8_t buf[32];
    for (int attempt = 0; attempt < 128; ++attempt) {
        rng.read(reinterpret_cast<char*>(buf), sizeof(buf));
        if (!rng) {
            return false;
        }
        U256 v;
        u256::from_bytes(v, buf, sizeof(buf));
        if (scalar_valid(fn, v)) {
            out = v;
            return true;
        }
    }
    return false;
}

bool EcProtocols::public_key(const Curve& curve, const U256& d, PointAff& q) {
    if (!scalar_valid(curve.order(), d)) {
        return false;
    }
    const PointJac p = ct::secure_mul(curve, curve.generator(), d);
    q = ec_to_affine(curve, p);
    return !q.infinity;
}

bool EcProtocols::generate_key(const Curve& curve, EcKeyPair& key) {
    if (!random_scalar(curve.order(), key.d)) {
        return false;
    }
    return public_key(curve, key.d, key.q);
}

bool EcProtocols::ecdh_shared(const Curve& curve, const U256& d, const PointAff& peer, U256& shared_x) {
    if (!scalar_valid(curve.order(), d) || peer.infinity || !curve.on_curve(peer)) {
        return false;
    }
    const PointJac p = ec_from_affine(curve, peer.x, peer.y);
    const PointAff shared = ec_to_affine(curve, ct::secure_mul(curve, p, d));
    if (shared.infinity) {
        return false;
    }
    shared_x = shared.x;
    return true;
}

bool EcProtocols::sign_with_k(const Curve& curve, const U256& d, const uint8_t hash[32], const U256& k,
                              U256& r, U256& s) {
    const Mod& fn = curve.order();
    if (!scalar_valid(fn, d) || !scalar_valid(fn, k)) {
        return false;
    }

    const PointAff c = ec_to_affine(curve, ct::secure_mul(curve, curve.generator(), k));
    if (c.infinity) {
        return false;
    }
    r = fn.reduce(c.x);
    if (u256::is_zero(r)) {
        return false;
    }

    const U256 e = hash_scalar(fn, hash);
    // s = k^{-1} * (e + r * d) mod q.
    const U256 km = fn.to_mont(k);
    const U256 dm = fn.to_mont(d);
    const U256 rm = fn.to_mont(r);
    const U256 em = fn.to_mont(e);
    const U256 sum = fn.add(em, fn.mul(rm, dm));
    s = fn.from_mont(fn.mul(fn.inv(km), sum));
    return !u256::is_zero(s);
}

bool EcProtocols::sign(const Curve& curve, const U256& d, const uint8_t hash[32], U256& r, U256& s) {
    for (int attempt = 0; attempt < 64; ++attempt) {
        U256 k;
        if (!random_scalar(curve.order(), k)) {
            return false;
        }
        if (sign_with_k(curve, d, hash, k, r, s)) {
            return true;
        }
    }
    return false;
}

bool EcProtocols::verify(const Curve& curve, const PointAff& q, const uint8_t hash[32], const U256& r,
                         const U256& s) {
    const Mod& fn = curve.order();
    if (!scalar_valid(fn, r) || !scalar_valid(fn, s) || q.infinity || !curve.on_curve(q)) {
        return false;
    }

    const U256 e = hash_scalar(fn, hash);
    const U256 w = fn.from_mont(fn.inv(fn.to_mont(s)));
    const U256 u1 = fn.from_mont(fn.mul(fn.to_mont(e), fn.to_mont(w)));
    const U256 u2 = fn.from_mont(fn.mul(fn.to_mont(r), fn.to_mont(w)));

    const PointJac p1 = ct::secure_mul(curve, curve.generator(), u1);
    const PointJac p2 = ct::secure_mul(curve, ec_from_affine(curve, q.x, q.y), u2);
    const PointAff c = ec_to_affine(curve, ec_add(curve, p1, p2));
    if (c.infinity) {
        return false;
    }
    return u256::compare(fn.reduce(c.x), r) == 0;
}
