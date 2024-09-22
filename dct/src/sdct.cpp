#include "sdct.hpp"

namespace sdct {

static constexpr float W1 = 0.980785280;
static constexpr float W2 = 0.923879533;
static constexpr float W3 = 0.831469612;
static constexpr float W4 = 0.707106781;
static constexpr float W5 = 0.555570233;
static constexpr float W6 = 0.382683432;
static constexpr float W7 = 0.195090322;

static constexpr auto W1_inv = 0.25f / W1;
static constexpr auto W2_inv = 0.25f / W2;
static constexpr auto W3_inv = 0.25f / W3;
static constexpr auto W4_inv = 0.25f / W4;
static constexpr auto W5_inv = 0.25f / W5;
static constexpr auto W6_inv = 0.25f / W6;
static constexpr auto W7_inv = 0.25f / W7;

// for external usage
std::array<float, 8> scales = {
    W4_inv, W1_inv, W2_inv, W3_inv, 
    W4_inv, W5_inv, W6_inv, W7_inv
};

static constexpr auto W2_sum_W6 = W2 + W6;
static constexpr auto W2_diff_W6 = W2 - W6;

template<bool scale>
void dct(float* data, size_t stride) {
    // stage 1
    float x0 = data[0 * stride] + data[7 * stride];
    float x1 = data[1 * stride] + data[6 * stride];
    float x2 = data[2 * stride] + data[5 * stride];
    float x3 = data[3 * stride] + data[4 * stride];
    float x4 = data[3 * stride] - data[4 * stride];
    float x5 = data[2 * stride] - data[5 * stride];
    float x6 = data[1 * stride] - data[6 * stride];
    float x7 = data[0 * stride] - data[7 * stride];

    //stage 2
    float y0 = x0 + x3;
    float y1 = x1 + x2;
    float y2 = x1 - x2;
    float y3 = x0 - x3;
    float y4 = x4 + x5;
    float y5 = W4 * (x5 + x6);
    float y6 = x6 + x7;
    float y7 = x7 - y5;
    x5 = y5 + x7;

    //stage 3
    data[0 * stride] = y0 + y1;
    data[4 * stride] = y0 - y1;
    x2 = W4 * (y2 + y3);
    float x8 = W6 * (y4 + y6);
    x4 = x8 + y6 * W2_diff_W6;
    x6 = x8 - y4 * W2_sum_W6;

    //stage 4
    data[2 * stride] = x2 + y3;
    data[6 * stride] = y3 - x2;
    data[1 * stride] = x4 + x5;
    data[3 * stride] = x6 + y7;
    data[5 * stride] = y7 - x6;
    data[7 * stride] = x5 - x4;

    if constexpr (scale) {
        data[0 * stride] *= W4_inv;
        data[1 * stride] *= W1_inv;
        data[2 * stride] *= W2_inv;
        data[3 * stride] *= W3_inv;
        data[4 * stride] *= W4_inv;
        data[5 * stride] *= W5_inv;
        data[6 * stride] *= W6_inv;
        data[7 * stride] *= W7_inv;
    }
}

template<bool scale>
void idct(float* data, size_t stride) {
    // stage 1
    float x4 = data[1 * stride];
    float x3 = data[2 * stride];
    float x2 = data[6 * stride];
    float x7 = data[7 * stride];
    float x0 = data[0 * stride];
    float x1 = data[4 * stride];
    float x5 = data[3 * stride];
    float x6 = data[5 * stride];

    if constexpr (scale) {
        x0 *= W4_inv;
        x1 *= W4_inv;
        x2 *= W6_inv;
        x3 *= W2_inv;
        x4 *= W1_inv;
        x5 *= W3_inv;
        x6 *= W5_inv;
        x7 *= W7_inv;
    }

    float y0 = x0 + x1;
    float y1 = x0 - x1;
    float y2 = x3 - x2;
    float y3 = x2 + x3;
    float y4 = x4 - x7;
    float y5 = x4 + x7;
    float y6 = x5 - x6;
    float y7 = x5 + x6;
    // stage 2
    x0 = y0;
    x1 = y1;
    x2 = W4 * y2;
    x3 = x2 + y3;
    x5 = W4 * (y5 - y7);
    x7 = y5 + y7;
    float y8 = W6 * (y4 + y6);
    x4 = y8 - y6 * W2_sum_W6;
    x6 = y8 + y4 * W2_diff_W6;

    // stage 3
    y0 = x0 + x3;
    y1 = x1 + x2;
    y2 = x1 - x2;
    y3 = x0 - x3;
    y4 = x4;
    y5 = x4 + x5;
    y6 = x5 + x6;
    y7 = x6 + x7;

    // stage 4
    data[0 * stride] = y0 + y7;
    data[1 * stride] = y1 + y6;
    data[2 * stride] = y2 + y5;
    data[3 * stride] = y3 + y4;
    data[4 * stride] = y3 - y4;
    data[5 * stride] = y2 - y5;
    data[6 * stride] = y1 - y6;
    data[7 * stride] = y0 - y7;
}

template<bool scale>
std::array<float, 8> dct(const std::array<float, 8>& data) {
    auto out = data;
    dct<scale>(out.data());
    return out;
}

template<bool scale>
std::array<float, 8> idct(const std::array<float, 8>& data) {
    auto out = data;
    idct<scale>(out.data());
    return out;
}

template void dct<true>(float* data, size_t stride);
template void dct<false>(float* data, size_t stride);

template void idct<true>(float* data, size_t stride);
template void idct<false>(float* data, size_t stride);

template std::array<float, 8> dct<true>(const std::array<float, 8>& data);
template std::array<float, 8> dct<false>(const std::array<float, 8>& data);

template std::array<float, 8> idct<true>(const std::array<float, 8>& data);
template std::array<float, 8> idct<false>(const std::array<float, 8>& data);

}  // namespace sdct
