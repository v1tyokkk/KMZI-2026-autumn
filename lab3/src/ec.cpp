#include "ec.hpp"

namespace {

// Параметры secp256k1: y^2 = x^3 + 7 mod p.
const char* kP = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F";
const char* kA = "0";
const char* kB = "7";
const char* kGx = "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798";
const char* kGy = "483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8";
const char* kQ = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141";

uint64_t ct_is_zero(const U256& v) {
    uint64_t x = v.w[0] | v.w[1] | v.w[2] | v.w[3];
    x |= x >> 32;
    x |= x >> 16;
    x |= x >> 8;
    x |= x >> 4;
    x |= x >> 2;
    x |= x >> 1;
    return (x & 1u) ^ 1u;
}

PointJac select_point(const PointJac& a, const PointJac& b, uint64_t bit) {
    return PointJac{u256::select(a.x, b.x, bit), u256::select(a.y, b.y, bit),
                    u256::select(a.z, b.z, bit)};
}

}

const Curve& Curve::secp256k1() {
    static const Curve curve;
    return curve;
}

Curve::Curve() : fp_(kP), fn_(kQ) {
    u256::from_hex(a_, kA);
    u256::from_hex(b_, kB);
    a_mont_ = fp_.to_mont(a_);
    b_mont_ = fp_.to_mont(b_);

    U256 gx;
    U256 gy;
    u256::from_hex(gx, kGx);
    u256::from_hex(gy, kGy);
    g_ = ec_from_affine(*this, gx, gy);
}

PointJac ec_infinity(const Curve& curve) {
    PointJac p;
    p.x = curve.field().one_mont();
    p.y = curve.field().one_mont();
    p.z = curve.field().zero();
    return p;
}

bool ec_is_infinity(const PointJac& p) { return u256::is_zero(p.z); }

PointJac ec_from_affine(const Curve& curve, const U256& x, const U256& y) {
    PointJac p;
    p.x = curve.field().to_mont(x);
    p.y = curve.field().to_mont(y);
    p.z = curve.field().one_mont();
    return p;
}

PointAff ec_to_affine(const Curve& curve, const PointJac& p) {
    PointAff out;
    if (ec_is_infinity(p)) {
        u256::zero(out.x);
        u256::zero(out.y);
        out.infinity = true;
        return out;
    }
    const Mod& f = curve.field();
    const U256 zi = f.inv(p.z);
    const U256 zi2 = f.mul(zi, zi);
    const U256 zi3 = f.mul(zi2, zi);
    out.x = f.from_mont(f.mul(p.x, zi2));
    out.y = f.from_mont(f.mul(p.y, zi3));
    out.infinity = false;
    return out;
}

PointJac ec_double(const Curve& curve, const PointJac& p) {
    const Mod& f = curve.field();
    const U256 y2 = f.mul(p.y, p.y);
    const U256 z2 = f.mul(p.z, p.z);
    const U256 z4 = f.mul(z2, z2);

    const U256 s = f.mul(f.mul(p.x, y2), f.from_small(4));
    const U256 m = f.add(f.mul(f.mul(p.x, p.x), f.from_small(3)), f.mul(curve.a_mont(), z4));
    const U256 x3 = f.sub(f.mul(m, m), f.add(s, s));
    const U256 y4 = f.mul(y2, y2);
    const U256 y3 = f.sub(f.mul(m, f.sub(s, x3)), f.mul(y4, f.from_small(8)));
    const U256 z3 = f.mul(f.add(p.y, p.y), p.z);

    PointJac out{x3, y3, z3};
    const uint64_t is_zero = ct_is_zero(p.z) | ct_is_zero(p.y);
    return select_point(out, ec_infinity(curve), is_zero);
}

PointJac ec_add(const Curve& curve, const PointJac& p, const PointJac& q) {
    const Mod& f = curve.field();
    const U256 z1z1 = f.mul(p.z, p.z);
    const U256 z2z2 = f.mul(q.z, q.z);
    const U256 u1 = f.mul(p.x, z2z2);
    const U256 u2 = f.mul(q.x, z1z1);
    const U256 s1 = f.mul(p.y, f.mul(q.z, z2z2));
    const U256 s2 = f.mul(q.y, f.mul(p.z, z1z1));
    const U256 h = f.sub(u2, u1);
    const U256 r = f.sub(s2, s1);

    const U256 h2 = f.mul(h, h);
    const U256 h3 = f.mul(h, h2);
    const U256 u1h2 = f.mul(u1, h2);
    const U256 x3 = f.sub(f.sub(f.mul(r, r), h3), f.add(u1h2, u1h2));
    const U256 y3 = f.sub(f.mul(r, f.sub(u1h2, x3)), f.mul(s1, h3));
    const U256 z3 = f.mul(f.mul(h, p.z), q.z);

    PointJac out{x3, y3, z3};
    const uint64_t h_zero = ct_is_zero(h);
    const uint64_t r_zero = ct_is_zero(r);

    // Особые случаи O, P=Q и P=-Q выбираются масками, без ветвлений.
    out = select_point(out, ec_double(curve, p), h_zero & r_zero);
    out = select_point(out, ec_infinity(curve), h_zero & (r_zero ^ 1u));
    out = select_point(out, q, ct_is_zero(p.z));
    out = select_point(out, p, ct_is_zero(q.z));
    return out;
}

PointJac ec_neg(const Curve& curve, const PointJac& p) {
    PointJac out = p;
    out.y = curve.field().sub(curve.field().zero(), p.y);
    return out;
}

bool ec_equal_affine(const PointAff& a, const PointAff& b) {
    if (a.infinity || b.infinity) {
        return a.infinity == b.infinity;
    }
    return u256::compare(a.x, b.x) == 0 && u256::compare(a.y, b.y) == 0;
}

bool Curve::on_curve(const PointAff& p) const {
    if (p.infinity) {
        return true;
    }
    if (u256::compare(p.x, fp_.modulus()) >= 0 || u256::compare(p.y, fp_.modulus()) >= 0) {
        return false;
    }
    const U256 x = fp_.to_mont(p.x);
    const U256 y = fp_.to_mont(p.y);
    const U256 lhs = fp_.mul(y, y);
    const U256 x2 = fp_.mul(x, x);
    const U256 rhs = fp_.add(fp_.add(fp_.mul(x2, x), fp_.mul(a_mont_, x)), b_mont_);
    return fp_.equal(lhs, rhs);
}
