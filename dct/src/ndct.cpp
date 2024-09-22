#include "ndct.hpp"

namespace ndct {

static constexpr float W1 = 0.490392625;
static constexpr float W2 = 0.461939752;
static constexpr float W3 = 0.415734798;
static constexpr float W4 = 0.353553385;
static constexpr float W5 = 0.277785122;
static constexpr float W6 = 0.191341713;
static constexpr float W7 = 0.097545162;

static constexpr std::array<float, 64> weights = {
    W4, W4, W4, W4, W4, W4, W4, W4,
    W1, W3, W5, W7, -W7, -W5, -W3, -W1,
    W2, W6, -W6, -W2, -W2, -W6, W6, W2,
    W3, -W7, -W1, -W5, W5, W1, W7, -W3,
    W4, -W4, -W4, W4, W4, -W4, -W4, W4,
    W5, -W1, W7, W3, -W3, -W7, W1, -W5,
    W6, -W2, W2, -W6, -W6, W2, -W2, W6,
    W7, -W5, W3, -W1, W1, -W3, W5, -W7
};

void dct(float* data) {
    float out[8];
    for (size_t i = 0; i < 8; ++i) {
        out[i] = 0;
        for (size_t j = 0; j < 8; ++j) out[i] += weights[i * 8 + j] * data[j];
    }
    for (size_t i = 0; i < 8; ++i) data[i] = out[i];
}

void idct(float* data) {
    float out[8];
    for (size_t i = 0; i < 8; ++i) {
        out[i] = 0;
        for (size_t j = 0; j < 8; ++j) out[i] += weights[i + j * 8] * data[j];
    }
    for (size_t i = 0; i < 8; ++i) data[i] = out[i];
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

}  // namespace ndct
