#pragma once

#include <array>

namespace fdct {

void dct(float* data);
void idct(float* data);

std::array<float, 8> dct(const std::array<float, 8>& data);
std::array<float, 8> idct(const std::array<float, 8>& data);

}  // namespace fdct
