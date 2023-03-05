//
// Created by Kirill Gavrilov on 21.02.2023.
//

#include <cmath>
#include <iostream>
#include "Statistics.cpp"

float polynomial(float x, const float *a, int n) {
    float result = 0;
    for (int i = n; i >= 0; i--) {
        result = std::fma(result, x, a[i]);
    }
    return result;
}

float kahan_sum(const float *x, int n) {
    float sum = 0;
    float tweak = 0;
    for (int i = 0; i < n; i++) {
        auto tweaked_x = x[i] - tweak;
        auto new_sum = sum + tweaked_x;
        tweak = (new_sum - sum) - tweaked_x;
        sum = new_sum;
    }
    return sum - tweak;
}

void check_kahan() {
    const int check_kahan_n = 10000;
    float check_kahan[check_kahan_n];
    for (int i = 0; i < check_kahan_n; i++) {
        check_kahan[i] = 1.0f / 3;
    }
    float sum = 0;
    for (int i = 0; i < check_kahan_n; i++) {
        sum += check_kahan[i];
    }
    std::cout << kahan_sum(check_kahan, check_kahan_n) << ' ' << sum << '\n';
}

float pairwise_sum_simd(float *x, int n) {
    //
    while (n > 16) {
        while (n % 8 > 0) {
            x[n / 2] += x[n - 1];
            n--;
        }

#pragma omp simd
        for (int i = 0; i < n / 2; i += 8) {
            x[i] += x[n / 2 + i];
            x[i + 1] += x[n / 2 + i + 1];
            x[i + 2] += x[n / 2 + i + 2];
            x[i + 3] += x[n / 2 + i + 3];
            x[i + 4] += x[n / 2 + i + 4];
            x[i + 5] += x[n / 2 + i + 5];
            x[i + 6] += x[n / 2 + i + 6];
            x[i + 7] += x[n / 2 + i + 7];
        }

        n /= 2;
    }

    while (n > 1) {
        x[0] += x[n - 1];
        n--;
    }

    return x[0];
}

float length(const float *x, int n) {
    if (n == 0) {
        return 0.0f;
    }
    float current_sum = 1.0f;
    float current_abs_max = abs(x[0]);
    for (int i = 1; i < n; i++) {
        float next_abs = abs(x[i]);
        if (next_abs > current_abs_max) {
            current_sum *= (current_abs_max / next_abs) * (current_abs_max / next_abs);
            current_abs_max = next_abs;
        }
        current_sum += (x[i] / current_abs_max) * (x[i] / current_abs_max);
    }

    return sqrt(current_sum) * current_abs_max;
}

void check_pairwise_sum_simd() {
    const int check_pairwise_n = 1000000;
    float check_pairwise[check_pairwise_n];
    for (int i = 0; i < check_pairwise_n; i++) {
        check_pairwise[i] = 1.0f / 10;
    }
    float sum = 0;
    for (int i = 0; i < check_pairwise_n; i++) {
        sum += check_pairwise[i];
    }
    std::cout << pairwise_sum_simd(check_pairwise, check_pairwise_n) << ' ' << sum << '\n';
}

void check_statistics() {
    Statistics statistics;
    for (int i = 0; i < 100000; i++) {
        statistics.update((i % 100) / 100.0f);
    }
    statistics.print();
}

void check_length() {
    float vec[5] = {1e12, 1e12, 1e12, 1e12, 1e12};
    std::cout << length(vec, 5) << '\n';
}

int main() {
    check_kahan();
    check_pairwise_sum_simd();
    check_statistics();
    check_pairwise_sum_simd();
    check_length();

    return 0;
}