#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace sim {
// weights are the column sums of the storage mass matrix. reduce_sum must
// sum across the communicator, including ranks with no local unknowns.
template <class Sum>
double ConservativeBounds(double* c, const double* weights, int size,
                        double lower, double upper, Sum reduce_sum)
{
    double removed = 0.0;
    for (int i = 0; i < size; ++i)
    {
        if (!std::isfinite(c[i]) || !std::isfinite(weights[i]) || weights[i] <= 0)
            removed = INFINITY;
        const double bounded = std::clamp(c[i], lower, upper);
        removed += weights[i] * (c[i] - bounded);
        c[i] = bounded;
    }
    removed = reduce_sum(removed);
    if (!std::isfinite(removed))
        throw std::runtime_error("Non-finite CH concentration or invalid storage weights");
    if (removed == 0.0) return 0.0;

    double capacity = 0.0;
    for (int i = 0; i < size; ++i)
        capacity += weights[i] * (removed > 0 ? upper - c[i] : c[i] - lower);
    capacity = reduce_sum(capacity);
    if (std::abs(removed) > capacity)
        throw std::runtime_error("CH inventory exceeds concentration bounds: particle group has no remaining capacity; reduce timestep or stop at capacity");
    const double fraction = std::abs(removed) / capacity;
    for (int i = 0; i < size; ++i)
        c[i] += fraction * (removed > 0 ? upper - c[i] : lower - c[i]);
    return removed;
}
}
