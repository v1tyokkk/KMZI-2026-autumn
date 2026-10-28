#include "ct.hpp"

namespace ct {

void cswap(PointJac& a, PointJac& b, uint64_t bit) {
    const uint64_t mask = 0 - (bit & 1u);
    U256* pa[3] = {&a.x, &a.y, &a.z};
    U256* pb[3] = {&b.x, &b.y, &b.z};
    for (int k = 0; k < 3; ++k) {
        for (int i = 0; i < 4; ++i) {
            const uint64_t diff = (pa[k]->w[i] ^ pb[k]->w[i]) & mask;
            pa[k]->w[i] ^= diff;
            pb[k]->w[i] ^= diff;
        }
    }
}

PointJac secure_mul(const Curve& curve, const PointJac& p, const U256& d) {
    PointJac r0 = ec_infinity(curve);
    PointJac r1 = p;
    // Биты от старшего к младшему: формула Джойя даёт dP при возврате r0.
    for (int i = 255; i >= 0; --i) {
        const uint64_t bit = u256::bit(d, i) ? 1u : 0u;
        cswap(r0, r1, bit);
        r1 = ec_add(curve, r0, r1);
        r0 = ec_double(curve, r0);
        cswap(r0, r1, bit);
    }
    return r0;
}

}
