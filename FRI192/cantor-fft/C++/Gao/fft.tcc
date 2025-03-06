#include <cstddef>
#include <omp.h>

#include <libff/common/profiling.hpp>
#include <libff/algebra/field_utils/field_utils.hpp>
#include "libiop/algebra/utils.hpp"

#include <Gao/fft.hpp>
#include <libff/algebra/field_utils/algorithms.hpp>

#include <iostream>
#include <vector>
#include "utils/utils.hpp"

namespace gao {

template<typename FieldT>
PreComputedValues_Level1<FieldT> pre_computation_lvl1(const libiop::affine_subspace<FieldT> &domain)
{
    const size_t m = domain.dimension();
    PreComputedValues_Level1<FieldT> values(m);

    std::vector<FieldT> betas2(domain.basis());
    FieldT shift2 = domain.shift();
    for (size_t j = 0; j < m; ++j)
    {
        values.mult_factor_vec.emplace_back(betas2[m-1-j]);

         /* compute deltas used in the reverse process */
        FieldT betainv = values.mult_factor_vec[j].inverse();
        for (size_t i = 0; i < m-1-j; ++i)
        {
            FieldT newbeta = betas2[i] * betainv;
            values.recursed_betas.emplace_back(newbeta);
            betas2[i] = newbeta.squared() - newbeta;
        }

        FieldT newshift = shift2 * betainv;
        values.recursed_shifts.emplace_back(newshift);
        shift2 = newshift.squared() - newshift;
    }
    return values;
}

template<typename FieldT>
PreComputedValues_Level2<FieldT> pre_computation_lvl2(const libiop::affine_subspace<FieldT> &domain)
{
    const size_t m = domain.dimension();
    PreComputedValues_Level2<FieldT> values(m);

    PreComputedValues_Level1<FieldT> values_lvl1 = pre_computation_lvl1<FieldT>(domain);

    size_t recursed_betas_ptr = (m * ( m - 1))/2;
    /* unwind the recursion */
    for (size_t j = 0; j < m; ++j)
    {
        FieldT beta = values_lvl1.mult_factor_vec[j];
        std::vector<FieldT> betai_vec((1ull<<(m-j))-1);

        FieldT betai(1);

        for (size_t k = 0; k < (1ull<<(m-j))-1; k++)
        {
            betai *= beta;
            betai_vec[k] = betai;
        }
        values.betai_vec.emplace_back(betai_vec);

        recursed_betas_ptr -= j;
        /* note that this devolves to empty range for the first loop iteration */
        std::vector<FieldT> popped_betas = std::vector<FieldT>(values_lvl1.recursed_betas.begin()+recursed_betas_ptr,
                                                               values_lvl1.recursed_betas.begin()+recursed_betas_ptr+j);
        const FieldT popped_shift = values_lvl1.recursed_shifts[m-1-j];
        values.sums_vec.emplace_back(libiop::all_subset_sums<FieldT>(popped_betas, popped_shift));
    }

    return values;
}

template<typename FieldT>
PreComputedValues_CO<FieldT> pre_computation_co(const libiop::affine_subspace<FieldT> &domain)
{
    const size_t m = domain.dimension();
    PreComputedValues_CO<FieldT> values(m);
    std::vector<FieldT> betas2(domain.basis());
    FieldT shift2 = domain.shift();

    for (size_t j = 0; j < m; ++j)
    {   
        values.recursed_shifts.emplace_back(shift2);
        shift2 = shift2.squared() - shift2;
    }
    std::vector<FieldT> popped_betas = std::vector<FieldT>(betas2.begin(), betas2.end()-1);
    values.beta_sums_vec = libiop::all_subset_sums<FieldT>(popped_betas, FieldT::zero());
    libiop::bitreverse_vector<FieldT>(values.beta_sums_vec);

    return values;
}


template<typename FieldT>
std::vector<FieldT> additive_FFT(const std::vector<FieldT> &poly_coeffs,
                                 const PreComputedValues_Level2<FieldT> &values)
{
    const size_t m = values.dimension;
    
    std::vector<FieldT> S(poly_coeffs);
    S.resize(1ull<<m, FieldT::zero());
    const size_t n = S.size();
    assert(n == (1ull<<m));

    for (size_t j = 0; j < m; ++j)
    {

        size_t k = 0;
        for (size_t ofs = (1ull<<j); ofs < n; ofs += (1ull<<j))
        {
            FieldT betai = values.betai_vec[j][k++];
            for (size_t p = 0; p < (1ull<<j); ++p)
            {
                S[ofs + p] *= betai;
            }

        }

        /* perform radix conversion */
        for (size_t stride = n/4; stride >= (1ul << j); stride >>= 1)
        {
            for (size_t ofs = 0; ofs < n; ofs += stride*4)
            {
                for (size_t i = 0; i < stride; ++i)
                {
                    S[ofs+2*stride+i] += S[ofs+3*stride+i];
                    S[ofs+1*stride+i] += S[ofs+2*stride+i];
                }
            }
        }
    }

    libiop::bitreverse_vector<FieldT>(S);

    /* unwind the recursion */
    for (size_t j = 0; j < m; ++j)
    {
        size_t stride = 1ull<<j;
        for (size_t ofs = 0; ofs < n; ofs += 2*stride)
        {
            for (size_t i = 0; i < stride; ++i)
            {
                S[ofs+i] += S[ofs+stride+i] * values.sums_vec[j][i];
                S[ofs+stride+i] += S[ofs+i];
            }
        }
    }
    assert(recursed_betas_ptr == 0);

    return S;
}


template<typename FieldT>
std::vector<FieldT> additive_FFT(const std::vector<FieldT> &poly_coeffs,
                                 const PreComputedValues_Level1<FieldT> &values)
{
    const size_t m = values.dimension;
    
    std::vector<FieldT> S(poly_coeffs);
    S.resize(1ull<<m, FieldT::zero());
    const size_t n = S.size();
    assert(n == (1ull<<m));

    std::vector<FieldT> recursed_betas = std::move(values.recursed_betas);
    std::vector<FieldT> recursed_shifts = std::move(values.recursed_shifts);
    size_t recursed_betas_ptr = (m * ( m - 1))/2;

    for (size_t j = 0; j < m; ++j)
    {
        FieldT beta = values.mult_factor_vec[j];
        FieldT betai(1);

        /* twist by beta. TODO: this can often be elided by a careful choice of betas */
        for (size_t ofs = (1ull<<j); ofs < n; ofs += (1ull<<j))
        {
            betai *= beta;
            for (size_t p = 0; p < (1ull<<j); ++p)
            {
                S[ofs + p] *= betai;
            }
        }

        /* perform radix conversion */
        for (size_t stride = n/4; stride >= (1ul << j); stride >>= 1)
        {
            for (size_t ofs = 0; ofs < n; ofs += stride*4)
            {
                for (size_t i = 0; i < stride; ++i)
                {
                    S[ofs+2*stride+i] += S[ofs+3*stride+i];
                    S[ofs+1*stride+i] += S[ofs+2*stride+i];
                }
            }
        }
    }

    libiop::bitreverse_vector<FieldT>(S);

    /* unwind the recursion */
    for (size_t j = 0; j < m; ++j)
    {
        recursed_betas_ptr -= j;
        /* note that this devolves to empty range for the first loop iteration */
        std::vector<FieldT> popped_betas = std::vector<FieldT>(recursed_betas.begin()+recursed_betas_ptr,
                                                               recursed_betas.begin()+recursed_betas_ptr+j);
        const FieldT popped_shift = recursed_shifts[m-1-j];
        std::vector<FieldT> sums = libiop::all_subset_sums<FieldT>(popped_betas, popped_shift);

        size_t stride = 1ull<<j;
        for (size_t ofs = 0; ofs < n; ofs += 2*stride)
        {
            for (size_t i = 0; i < stride; ++i)
            {
                S[ofs+i] += S[ofs+stride+i] * sums[i];
                S[ofs+stride+i] += S[ofs+i];
            }
        }
    }
    assert(recursed_betas_ptr == 0);

    return S;
}

template<typename FieldT>
std::vector<FieldT> additive_FFT_CO(const std::vector<FieldT> &poly_coeffs,
                                 const PreComputedValues_CO<FieldT> &values)
{
    const size_t m = values.dimension;
    
    std::vector<FieldT> S(poly_coeffs);
    S.resize(1ull<<m, FieldT::zero());

    const size_t n = S.size();
    assert(n == (1ull<<m));

    for (size_t j = 0; j < m; ++j)
    {
        /* perform radix conversion */
        for (size_t stride = n/4; stride >= (1ul << j); stride >>= 1)
        {
            for (size_t ofs = 0; ofs < n; ofs += stride*4)
            {
                for (size_t i = 0; i < stride; ++i)
                {
                    S[ofs+2*stride+i] += S[ofs+3*stride+i];
                    S[ofs+1*stride+i] += S[ofs+2*stride+i];
                }
            }
        }
    }

    // libiop::bitreverse_vector<FieldT>(S);

    // /* unwind the recursion */
    // for (size_t j = 0; j < m; ++j)
    // {
    //     size_t stride = 1ull<<j;
    //     size_t cnt = m-j-1;
    //     const FieldT popped_shift = values.recursed_shifts[m-1-j];

    //     for (size_t ofs = 0; ofs < n; ofs += 2*stride)
    //     {   
    //         for (size_t i = 0; i < stride; ++i)
    //         {
    //             S[ofs+i] += S[ofs+stride+i] * (values.beta_sums_vec[i<<(cnt)] + popped_shift);
    //             S[ofs+stride+i] += S[ofs+i];
    //         }
    //     }
    // }
    // return S;

    for (int j = m-1; j>=0; --j)
    {
        size_t step = 1ull << (j);
        size_t cnt = 0;
        const FieldT popped_shift = values.recursed_shifts[j];
        for (size_t ofs = 0; ofs < n; ofs += 2*step)
        {   
            const FieldT mult_factor = values.beta_sums_vec[cnt++] + popped_shift;
            for (size_t i = 0; i < step; ++i)
            {
                S[ofs+i] += S[ofs+step+i] * (mult_factor);
                S[ofs+step+i] += S[ofs+i];
            }       
        }
    }
    libiop::bitreverse_vector<FieldT>(S);
    return S;
}


template<typename FieldT>
std::vector<FieldT> additive_FFT_CO(const std::vector<FieldT> &poly_coeffs,
                                 const libiop::affine_subspace<FieldT> &domain)
{
    std::vector<FieldT> S(poly_coeffs);
    S.resize(domain.num_elements(), FieldT::zero());

    const size_t n = S.size();
    const size_t m = domain.dimension();
    assert(n == (1ull<<m));

    std::vector<FieldT> recursed_shifts(m, FieldT(0));

    std::vector<FieldT> betas2(domain.basis());
    FieldT shift2 = domain.shift();
    for (size_t j = 0; j < m; ++j)
    {

        /* perform radix conversion */
        for (size_t stride = n/4; stride >= (1ul << j); stride >>= 1)
        {
            for (size_t ofs = 0; ofs < n; ofs += stride*4)
            {
                for (size_t i = 0; i < stride; ++i)
                {
                    S[ofs+2*stride+i] += S[ofs+3*stride+i];
                    S[ofs+1*stride+i] += S[ofs+2*stride+i];
                }
            }
        }
        recursed_shifts[j] = shift2;
        shift2 = shift2.squared() - shift2;
    }

    libiop::bitreverse_vector<FieldT>(S);

    /* unwind the recursion */
    for (size_t j = 0; j < m; ++j)
    {
        /* note that this devolves to empty range for the first loop iteration */
        std::vector<FieldT> popped_betas = std::vector<FieldT>(betas2.end()-1 - j,
                                                               betas2.end()-1);
        const FieldT popped_shift = recursed_shifts[m-1-j];
        std::vector<FieldT> sums = libiop::all_subset_sums<FieldT>(popped_betas, popped_shift);

        size_t stride = 1ull<<j;
        for (size_t ofs = 0; ofs < n; ofs += 2*stride)
        {
            for (size_t i = 0; i < stride; ++i)
            {
                S[ofs+i] += S[ofs+stride+i] * sums[i];
                S[ofs+stride+i] += S[ofs+i];
            }
        }
    }

    return S;
}

template<typename FieldT>
std::vector<FieldT> additive_IFFT(const std::vector<FieldT> &evals,
                                  const libiop::affine_subspace<FieldT> &domain)
{
    const size_t n = evals.size();
    const size_t m = domain.dimension();
    assert(n == (1ull<<m));

    std::vector<FieldT> S(evals);
    std::vector<FieldT> recursed_twists(m, FieldT(0));

    std::vector<FieldT> betas2(domain.basis());
    FieldT shift2 = domain.shift();
    for (size_t j = 0; j < m; ++j)
    {
        const FieldT beta = betas2[m-1-j];
        const FieldT betainv = beta.inverse();
        recursed_twists[j] = betainv;

        std::vector<FieldT> newbetas(m-1-j, FieldT(0));

        for (size_t i = 0; i < m-1-j; ++i)
        {
            FieldT newbeta = betas2[i] * betainv;
            newbetas[i] = newbeta;
            betas2[i] = newbeta.squared() - newbeta;
        }

        FieldT newshift = shift2 * betainv;
        shift2 = newshift.squared() - newshift;

        const std::vector<FieldT> sums = libiop::all_subset_sums<FieldT>(newbetas, newshift);

        const size_t half = 1ull<<(m-1-j);
        for (size_t ofs = 0; ofs < n; ofs += 2*half)
        {
            for (size_t p = 0; p < half; ++p)
            {
                S[ofs + half + p] += S[ofs + p];
                S[ofs + p] += S[ofs + half + p] * sums[p];
            }
        }
    }

    libiop::bitreverse_vector<FieldT>(S);

    for (size_t j = 0; j < m; ++j)
    {
        size_t N = 4ull<<(m-1-j);
        /* perform radix combinations */
        while (N <= n)
        {
            const size_t quarter = N/4;
            for (size_t ofs = 0; ofs < n; ofs += N)
            {
                for (size_t i = 0; i < quarter; ++i)
                {
                    S[ofs+1*quarter+i] += S[ofs+2*quarter+i];
                    S[ofs+2*quarter+i] += S[ofs+3*quarter+i];
                }
            }
            N *= 2;
        }

        /* twist by \beta^{-1} */
        const FieldT betainv = recursed_twists[m-1-j];
        FieldT betainvi(1);
        for (size_t ofs = 0; ofs < n; ofs += (1ull<<(m-1-j)))
        {
            for (size_t p = 0; p < (1ull<<(m-1-j)); ++p)
            {
                S[ofs + p] *= betainvi;
            }
            betainvi *= betainv;
        }
    }

    return S;
}

} // namespace gao

