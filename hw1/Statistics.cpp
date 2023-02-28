//
// Created by Kirill Gavrilov on 21.02.2023.
//

#include <vector>
#include <iostream>

class Statistics {
private:
    std::vector<float> items;
    float current_mean = 0;
    float current_m_n = 0;
public:
    void update(float x) {
        items.push_back(x);
        float old_mean = current_mean;
        current_mean += (x - current_mean) / count();
        current_m_n += (x - current_mean) * (x - old_mean);
    }

    int count() const noexcept {
        return items.size();
    }

    float min() const noexcept {
        if (count() == 0) {
            return 0.0f;
        }
        float current_min = items[0];
        for (int i = 1; i < items.size(); i++) {
            current_min = std::min(current_min, items[i]);
        }
        return current_min;
    }

    float max() const noexcept {
        if (count() == 0) {
            return 0.0f;
        }
        float current_max = items[0];
        for (int i = 1; i < items.size(); i++) {
            current_max = std::max(current_max, items[i]);
        }
        return current_max;
    }

    float sum() const noexcept {
        float sum = 0;
        for (float a: items) {
            sum += a;
        }
        return sum;
    }

    float mean() const noexcept {
        return current_mean;
    }

    float variance() const noexcept {
        return current_m_n / count();
    }

    void print() {
        std::cout << "Statistics" << ", "
                  << "count: " << count() << ", "
                  << "min: " << min() << ", "
                  << "max: " << max() << ", "
                  << "sum: " << sum() << ", "
                  << "mean: " << mean() << ", "
                  << "variance: " << variance() << '\n';
    }
};