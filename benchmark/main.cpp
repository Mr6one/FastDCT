#include <iostream>
#include <chrono>
#include <functional>
#include <iomanip>

#include "ndct.hpp"
#include "fdct.hpp"
#include "sdct.hpp"

template<size_t N, typename F, typename ...Args>
void benchmark_function(F&& f, Args&& ...args) {
    auto start = std::chrono::steady_clock::now();
    for (size_t i = 0; i < N; ++i) std::forward<F>(f)(std::forward<Args>(args)...);
    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << std::setprecision(3) << "total: " << duration.count() / 1000.0 << "s" << 
                 "\t iteration: "  << 1e6 * duration.count() / (N + 0.0) << "ns" << std::endl;
}

int main() {
    constexpr std::array<float, 8> samples = {1, 1, 1, 3, 1, 1, 1, 7};
    constexpr size_t iters = 1 << 26;

    using F = std::array<float, 8>(*)(const std::array<float, 8>&);
    std::cout << "naive: \t";
    benchmark_function<iters>(static_cast<F>(ndct::dct), samples);
    std::cout << "fast: \t";
    benchmark_function<iters>(static_cast<F>(fdct::dct), samples);
    std::cout << "scale: \t";
    benchmark_function<iters>(static_cast<F>(sdct::dct<true>), samples);
    std::cout << "scale: \t";
    benchmark_function<iters>(static_cast<F>(sdct::dct<false>), samples);

    return 0;
}
