#include <benchmark/benchmark.h>
#include <cstdlib>
#include <libff/algebra/fields/binary/gf256.hpp>
#include "libiop/algebra/subspace.hpp"
#include "libiop/algebra/utils.hpp"
#include "libiop/fft.hpp"
#include "Cantor/fft.hpp"
#include "Gao/fft.hpp"

// Benchmark for libiop::naive_FFT -----------------------------------------------------------------------------------------------------------------------
static void BM_libiop_naive_fft(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    libiop::field_subset<FieldT> domain = libiop::field_subset<FieldT>(libiop::affine_subspace<FieldT>::random_affine_subspace(m));
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(poly_coeffs);
        benchmark::DoNotOptimize(domain);
        benchmark::DoNotOptimize(result = libiop::naive_FFT<FieldT>(poly_coeffs, domain));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

// Benchmark for libiop::additive_FFT -----------------------------------------------------------------------------------------------------------------------
static void BM_libiop_additive_fft(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    libiop::field_subset<FieldT> domain = libiop::field_subset<FieldT>(libiop::affine_subspace<FieldT>::random_affine_subspace(m));
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(poly_coeffs);
        benchmark::DoNotOptimize(domain);
        benchmark::DoNotOptimize(result = libiop::additive_FFT<FieldT>(poly_coeffs, domain.subspace()));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

// Benchmark for gao::additive_FFT (lvl1) -----------------------------------------------------------------------------------------------------------------------
static void BM_gao_additive_fft_lvl1(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    libiop::field_subset<FieldT> domain = libiop::field_subset<FieldT>(libiop::affine_subspace<FieldT>::random_affine_subspace(m));
    gao::PreComputedValues_Level1<FieldT> values_lvl1 = gao::pre_computation_lvl1<FieldT>(domain.subspace());
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(poly_coeffs);
        benchmark::DoNotOptimize(domain);
        benchmark::DoNotOptimize(result = gao::additive_FFT<FieldT>(poly_coeffs, values_lvl1));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

// Benchmark for gao::additive_FFT (lvl2) -----------------------------------------------------------------------------------------------------------------------
static void BM_gao_additive_fft_lvl2(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    libiop::field_subset<FieldT> domain = libiop::field_subset<FieldT>(libiop::affine_subspace<FieldT>::random_affine_subspace(m));
    gao::PreComputedValues_Level2<FieldT> values_lvl2 = gao::pre_computation_lvl2<FieldT>(domain.subspace());
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(poly_coeffs);
        benchmark::DoNotOptimize(domain);
        benchmark::DoNotOptimize(result = gao::additive_FFT<FieldT>(poly_coeffs, values_lvl2));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

// Benchmark for gao::additive_FFT_CO -----------------------------------------------------------------------------------------------------------------------
static void BM_gao_additive_fft_co(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    std::vector<FieldT> basis(cantor_basis<FieldT>(m));
    std::reverse(basis.begin(), basis.end());

    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    libiop::field_subset<FieldT> domain{libiop::affine_subspace<FieldT>(basis, FieldT::random_element())};
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(poly_coeffs);
        benchmark::DoNotOptimize(domain);
        benchmark::DoNotOptimize(result = gao::additive_FFT_CO<FieldT>(poly_coeffs, domain.subspace()));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

// Benchmark for gao::additive_FFT_CO (lvl2) -----------------------------------------------------------------------------------------------------------------------
static void BM_gao_additive_fft_co_lvl2(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    std::vector<FieldT> basis(cantor_basis<FieldT>(m));
    std::reverse(basis.begin(), basis.end());

    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    libiop::field_subset<FieldT> domain{libiop::affine_subspace<FieldT>(basis, FieldT::random_element())};
    gao::PreComputedValues_CO<FieldT> values = gao::pre_computation_co<FieldT>(domain.subspace());
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(poly_coeffs);
        benchmark::DoNotOptimize(domain);
        benchmark::DoNotOptimize(result = gao::additive_FFT_CO<FieldT>(poly_coeffs, values));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

// Benchmark for cantor::additive_FFT -----------------------------------------------------------------------------------------------------------------------
static void BM_cantor_additive_fft(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    std::vector<FieldT> basis(cantor_basis<FieldT>(m));

    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    libiop::field_subset<FieldT> domain{libiop::affine_subspace<FieldT>(basis, FieldT::random_element())};
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(poly_coeffs);
        benchmark::DoNotOptimize(domain);
        benchmark::DoNotOptimize(result = cantor::additive_FFT<FieldT>(poly_coeffs, domain.subspace()));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

// Benchmark for cantor::additive_FFT (precmp)-----------------------------------------------------------------------------------------------------------------------
static void BM_cantor_additive_fft_precmp(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    std::vector<FieldT> basis(cantor_basis<FieldT>(m));

    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    libiop::field_subset<FieldT> domain{libiop::affine_subspace<FieldT>(basis, FieldT::random_element())};
    cantor::PreComputedValues<FieldT> values = cantor::pre_computation(domain.subspace());
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(poly_coeffs);
        benchmark::DoNotOptimize(domain);
        benchmark::DoNotOptimize(result = cantor::additive_FFT<FieldT>(poly_coeffs, values));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

const int MIN_RANGE = std::stoi(std::getenv("BM_MIN_RANGE"));
const int MAX_RANGE = std::stoi(std::getenv("BM_MAX_RANGE"));
const int STEP = std::stoi(std::getenv("BM_STEP"));

BENCHMARK(BM_libiop_additive_fft)->DenseRange(MIN_RANGE, MAX_RANGE, STEP)->Unit(benchmark::kMicrosecond)->ReportAggregatesOnly(true);
BENCHMARK(BM_gao_additive_fft_lvl1)->DenseRange(MIN_RANGE, MAX_RANGE, STEP)->Unit(benchmark::kMicrosecond)->ReportAggregatesOnly(true);
BENCHMARK(BM_gao_additive_fft_lvl2)->DenseRange(MIN_RANGE, MAX_RANGE, STEP)->Unit(benchmark::kMicrosecond)->ReportAggregatesOnly(true);
BENCHMARK(BM_gao_additive_fft_co)->DenseRange(MIN_RANGE, MAX_RANGE, STEP)->Unit(benchmark::kMicrosecond)->ReportAggregatesOnly(true);
BENCHMARK(BM_gao_additive_fft_co_lvl2)->DenseRange(MIN_RANGE, MAX_RANGE, STEP)->Unit(benchmark::kMicrosecond)->ReportAggregatesOnly(true);
BENCHMARK(BM_cantor_additive_fft)->DenseRange(MIN_RANGE, MAX_RANGE, STEP)->Unit(benchmark::kMicrosecond)->ReportAggregatesOnly(true);
BENCHMARK(BM_cantor_additive_fft_precmp)->DenseRange(MIN_RANGE, MAX_RANGE, STEP)->Unit(benchmark::kMicrosecond)->ReportAggregatesOnly(true);
BENCHMARK(BM_libiop_naive_fft)->DenseRange(MIN_RANGE, 10, STEP)->Unit(benchmark::kMicrosecond)->ReportAggregatesOnly(true);


BENCHMARK_MAIN();
