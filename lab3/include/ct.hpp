#pragma once

#include "ec.hpp"

// Модуль 3. Constant-time операции со скаляром и точками.
namespace ct {

void cswap(PointJac& a, PointJac& b, uint64_t bit);

// Лестница Джойя: ровно 256 сложений и 256 удвоений, без ветвлений по битам d.
PointJac secure_mul(const Curve& curve, const PointJac& p, const U256& d);

}
