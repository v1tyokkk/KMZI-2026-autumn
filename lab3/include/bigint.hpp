#pragma once

#include <cstdint>
#include <string>

// 256-битное число: 4 64-битных слова, младшее слово первое.
struct U256 {
    uint64_t w[4];
};

namespace u256 {

void zero(U256& v);
void set_word(U256& v, uint64_t x);
bool from_hex(U256& out, const char* hex);
std::string to_hex(const U256& v);
bool from_bytes(U256& out, const uint8_t* be, int len);
void to_bytes(const U256& v, uint8_t* be);
bool is_zero(const U256& v);
int compare(const U256& a, const U256& b);
bool bit(const U256& v, int i);
U256 select(const U256& a, const U256& b, uint64_t bit);

}

// Арифметика по модулю 256-битного простого числа в форме Монтгомери.
class Mod {
public:
    explicit Mod(const char* modulus_hex);

    const U256& modulus() const { return m_; }
    U256 zero() const;
    U256 one_mont() const { return one_; }
    U256 from_small(uint64_t x) const;

    U256 to_mont(const U256& a) const;
    U256 from_mont(const U256& a) const;

    U256 add(const U256& a, const U256& b) const;
    U256 sub(const U256& a, const U256& b) const;
    U256 mul(const U256& a, const U256& b) const;
    U256 inv(const U256& a) const;

    // Приведение обычного числа (не Монтгомери) по модулю.
    U256 reduce(const U256& a) const;

    bool equal(const U256& a, const U256& b) const;
    bool is_zero(const U256& a) const;

private:
    U256 m_;
    U256 r2_;
    U256 one_;
    uint64_t n0_;

    U256 add_raw(const U256& a, const U256& b) const;
    void mont_mul(U256& out, const U256& a, const U256& b) const;
    U256 conditional_sub(const U256& a, uint64_t high) const;
};
