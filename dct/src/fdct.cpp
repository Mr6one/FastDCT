#include "fdct.hpp"

namespace fdct {

static constexpr float W1 = 0.490392625;
static constexpr float W2 = 0.461939752;
static constexpr float W3 = 0.415734798;
static constexpr float W4 = 0.353553385;
static constexpr float W5 = 0.277785122;
static constexpr float W6 = 0.191341713;
static constexpr float W7 = 0.097545162;

static constexpr auto W1_sum_W7 = W1 + W7;
static constexpr auto W1_diff_W7 = W1 - W7;

static constexpr auto W5_sum_W3 = W3 + W5;
static constexpr auto W5_diff_W3 = W5 - W3;

static constexpr auto W2_sum_W6 = W2 + W6;
static constexpr auto W2_diff_W6 = W2 - W6;
static constexpr auto W4x2 = 2 * W4;

void dct(float* data) {
    // stage 1
    float x0 = data[0] + data[7];
    float x1 = data[1] + data[6];
    float x2 = data[2] + data[5];
    float x3 = data[3] + data[4];
    float x4 = data[3] - data[4];
    float x5 = data[2] - data[5];
    float x6 = data[1] - data[6];
    float x7 = data[0] - data[7];

    //stage 2
    float y0 = x0 + x3;
    float y1 = x1 + x2;
    float y2 = x1 - x2;
    float y3 = x0 - x3;
    float x8 = W7 * (x4 + x7);
    float y4 = x8 + x7 * W1_diff_W7;
    float y7 = x8 - x4 * W1_sum_W7;
    x8 = W3 * (x5 + x6);
    float y5 = x8 + x5 * W5_diff_W3;
    float y6 = x8 - x6 * W5_sum_W3;

    // stage 3
    data[0] = W4 * (y0 + y1);
    data[4] = W4 * (y0 - y1);
    x8 = W6 * (y2 + y3);
    data[6] = x8 - y2 * W2_sum_W6;
    data[2] = x8 + y3 * W2_diff_W6;
    data[1] = y4 + y5;
    data[7] = y6 + y7;

    // stage 4
    x5 = y4 - y5;
    x6 = y6 - y7;
    data[3] = W4x2 * (x5 - x6);
    data[5] = W4x2 * (x5 + x6);
}

void idct(float* data) {
    // stage 1
    float x4 = data[1];
    float x3 = data[2];
    float x2 = data[6];
    float x7 = data[7];
    float x0 = W4 * data[0];
    float x1 = W4 * data[4];
    float x5 = W4x2 * data[3];
    float x6 = W4x2 * data[5];
    
    // stage 2
    float y5 = x5 + x6;
    float y6 = x6 - x5;

    float y0 = x0 + x1;
    float y1 = x0 - x1;
    float x8 = W6 * (x2 + x3);
    float y2 = x8 - x2 * W2_sum_W6;
    float y3 = x8 + x3 * W2_diff_W6;
    float y4 = x4 + y5;
    float y7 = x7 - y6;
    y5 = x4 - y5;
    y6 = y6 + x7;

    // stage 3
    x0 = y0 + y3;
    x1 = y1 + y2;
    x2 = y1 - y2;
    x3 = y0 - y3;
    x8 = W7 * (y4 + y7);
    x4 = x8 - y7 * W1_sum_W7;
    x7 = x8 + y4 * W1_diff_W7;
    x8 = W3 * (y5 + y6);
    x5 = x8 + y5 * W5_diff_W3;
    x6 = x8 - y6 * W5_sum_W3;

    // stage 4
    data[0] = x0 + x7;
    data[1] = x1 + x6;
    data[2] = x2 + x5;
    data[3] = x3 + x4;
    data[4] = x3 - x4;
    data[5] = x2 - x5;
    data[6] = x1 - x6;
    data[7] = x0 - x7;
}

std::array<float, 8> dct(const std::array<float, 8>& data) {
    auto out = data;
    dct(out.data());
    return out;
}

std::array<float, 8> idct(const std::array<float, 8>& data) {
    auto out = data;
    idct(out.data());
    return out;
}

}  // namespace fdct
