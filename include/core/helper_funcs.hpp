#ifndef HELPER_FUNCS_HPP
#define HELPER_FUNCS_HPP

template <typename T>
bool in_range(const T& x, const T& lowerlim, const T& upperlim)
{
    return ((x >= lowerlim) && (x <= upperlim));
}

#endif