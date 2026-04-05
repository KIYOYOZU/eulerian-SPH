#pragma once

#include "types.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace eulerian_sph
{
using LoopIndex = long long;

inline bool shouldParallelize(const std::size_t work_items)
{
#ifdef _OPENMP
    return work_items >= 256;
#else
    (void)work_items;
    return false;
#endif
}

namespace parallel
{
template <class Function>
inline void forRange(const std::size_t count, const Function &function)
{
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(shouldParallelize(count))
    for (LoopIndex index = 0; index < static_cast<LoopIndex>(count); ++index)
    {
        function(static_cast<std::size_t>(index));
    }
#else
    for (std::size_t index = 0; index != count; ++index)
    {
        function(index);
    }
#endif
}

template <class Function>
inline Scalar reduceMax(const std::size_t count, const Scalar initial_value, const Function &function)
{
    Scalar result = initial_value;
#ifdef _OPENMP
#pragma omp parallel if(shouldParallelize(count))
    {
        Scalar local_result = result;
#pragma omp for nowait
        for (LoopIndex index = 0; index < static_cast<LoopIndex>(count); ++index)
        {
            local_result = std::max(local_result, function(static_cast<std::size_t>(index)));
        }
#pragma omp critical
        { result = std::max(result, local_result); }
    }
#else
    for (std::size_t index = 0; index != count; ++index)
    {
        result = std::max(result, function(index));
    }
#endif
    return result;
}

template <class Function>
inline Scalar reduceMin(const std::size_t count, const Scalar initial_value, const Function &function)
{
    Scalar result = initial_value;
#ifdef _OPENMP
#pragma omp parallel if(shouldParallelize(count))
    {
        Scalar local_result = result;
#pragma omp for nowait
        for (LoopIndex index = 0; index < static_cast<LoopIndex>(count); ++index)
        {
            local_result = std::min(local_result, function(static_cast<std::size_t>(index)));
        }
#pragma omp critical
        { result = std::min(result, local_result); }
    }
#else
    for (std::size_t index = 0; index != count; ++index)
    {
        result = std::min(result, function(index));
    }
#endif
    return result;
}
} // namespace parallel
} // namespace eulerian_sph
