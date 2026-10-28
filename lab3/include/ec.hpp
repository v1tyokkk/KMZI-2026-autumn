#pragma once

#include "bigint.hpp"

// Точка в проективных координатах Якоби: x = X/Z^2, y = Y/Z^3.
struct PointJac {
    U256 x;
    U256 y;
    U256 z;
};

struct PointAff {
    U256 x;
    U256 y;
    bool infinity;
};

// Модуль 1. Параметры кривой и алгебра точек над GF(p).
class Curve {
public:
    static const Curve& secp256k1();

    const Mod& field() const { return fp_; }
    const Mod& order() const { return fn_; }
    const PointJac& generator() const { return g_; }
    const U256& a_mont() const { return a_mont_; }
    const U256& b_mont() const { return b_mont_; }

    bool on_curve(const PointAff& p) const;

private:
    Curve();

    Mod fp_;
    Mod fn_;
    U256 a_;
    U256 b_;
    U256 a_mont_;
    U256 b_mont_;
    PointJac g_;
};

PointJac ec_infinity(const Curve& curve);
bool ec_is_infinity(const PointJac& p);
PointJac ec_from_affine(const Curve& curve, const U256& x, const U256& y);
PointAff ec_to_affine(const Curve& curve, const PointJac& p);
PointJac ec_add(const Curve& curve, const PointJac& p, const PointJac& q);
PointJac ec_double(const Curve& curve, const PointJac& p);
PointJac ec_neg(const Curve& curve, const PointJac& p);
bool ec_equal_affine(const PointAff& a, const PointAff& b);
