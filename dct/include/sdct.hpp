#pragma once

#include <array>

namespace sdct {

extern std::array<float, 8> scales;

template<bool scale = true>
void dct(float* data, size_t stride = 1);

template<bool scale = true>
void idct(float* data, size_t stride = 1);

template<bool scale = true>
std::array<float, 8> dct(const std::array<float, 8>& data);

template<bool scale = true>
std::array<float, 8> idct(const std::array<float, 8>& data);

}  // namespace sdct
