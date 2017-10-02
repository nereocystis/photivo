#ifndef CLAMP_H
#define CLAMP_H

#include <cassert>
#include <functional>

template <class S, class T, class Compare>
constexpr const S &clamp(const T &v, const T &lo, const T &hi, Compare comp) {
  // assert(!comp(hi, lo));
  return comp(v, lo) ? lo : comp(hi, v) ? hi : v;
}

template <class S, class T>
constexpr const S &clamp(const T &v, const T &lo, const T &hi) {
  return clamp(v, lo, hi, std::less<T>());
}

template <class T> constexpr const uint16_t &clamp16(const T &v) {
  return clamp<uint16_t>(static_cast<int32_t>(v), static_cast<int32_t>(0),
                         static_cast<int32_t>(0x10000), std::less<T>());
}

#endif // CLAMP_H
