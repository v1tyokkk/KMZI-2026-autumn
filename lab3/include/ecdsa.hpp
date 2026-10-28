#pragma once

#include "ec.hpp"

struct EcKeyPair {
    U256 d;
    PointAff q;
};

// Модуль 2. ECDH и ECDSA над заданной кривой.
class EcProtocols {
public:
    static bool generate_key(const Curve& curve, EcKeyPair& key);
    static bool public_key(const Curve& curve, const U256& d, PointAff& q);

    static bool ecdh_shared(const Curve& curve, const U256& d, const PointAff& peer, U256& shared_x);

    static bool sign(const Curve& curve, const U256& d, const uint8_t hash[32], U256& r, U256& s);
    static bool sign_with_k(const Curve& curve, const U256& d, const uint8_t hash[32], const U256& k,
                            U256& r, U256& s);
    static bool verify(const Curve& curve, const PointAff& q, const uint8_t hash[32], const U256& r,
                       const U256& s);

private:
    static U256 hash_scalar(const Mod& fn, const uint8_t hash[32]);
    static bool random_scalar(const Mod& fn, U256& out);
};
