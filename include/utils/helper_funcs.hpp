#ifndef HELPER_FUNCS_HPP
#define HELPER_FUNCS_HPP

#include "../core/base_types.hpp"
#include <cassert>
#include <math.h>

template <typename T> bool InRange(const T &x, const T &lowerlim, const T &upperlim) {
    return ((x >= lowerlim) && (x <= upperlim));
}

Matrix CholeskyDecomposition(const Matrix &matrix);

#endif
