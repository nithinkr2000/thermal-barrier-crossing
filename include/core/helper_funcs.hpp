#ifndef HELPER_FUNCS_HPP
#define HELPER_FUNCS_HPP

template <typename T>
bool InRange(const T& x, const T& lowerlim, const T& upperlim)
{
    return ((x >= lowerlim) && (x <= upperlim));
}

#endif