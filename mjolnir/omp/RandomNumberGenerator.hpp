#ifndef MJOLNIR_OMP_RANDOM_NUMBER_GENERATOR
#define MJOLNIR_OMP_RANDOM_NUMBER_GENERATOR
#include <mjolnir/util/logger.hpp>
#include <mjolnir/util/aligned_allocator.hpp>
#include <mjolnir/util/aligned_storage.hpp>
#include <mjolnir/core/RandomNumberGenerator.hpp>
#include <mjolnir/omp/OpenMPSimulatorTraits.hpp>

namespace mjolnir
{

template<typename realT, template<typename, typename> class boundaryT>
class RandomNumberGenerator<OpenMPSimulatorTraits<realT, boundaryT>>
{
  public:
    using traits_type     = OpenMPSimulatorTraits<realT, boundaryT>;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;

    static constexpr std::size_t cache_alignment = 64;

    using rng_type = aligned_storage<std::mt19937, cache_alignment>;
    using nrm_type = aligned_storage<std::normal_distribution<real_type>,
                                     cache_alignment>;

  public:
    explicit RandomNumberGenerator(const std::uint32_t seed)
        : seed_(seed),
          rngs_(omp_get_max_threads()),
          nrms_(omp_get_max_threads(),
                nrm_type(std::normal_distribution<real_type>(0.0, 1.0)))
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        std::seed_seq              sseq{seed};
        std::vector<std::uint32_t> buf(omp_get_max_threads());
        sseq.generate(buf.begin(), buf.end());
        for(std::size_t i=0, e=omp_get_max_threads(); i<e; ++i)
        {
            MJOLNIR_LOG_INFO("the ", i, "-th RNG is seeded by ", buf.at(i));
            rngs_.at(i).value.seed(buf.at(i));
        }
    }
    ~RandomNumberGenerator() = default;

    explicit RandomNumberGenerator(const std::string& internal_state)
        : rngs_(omp_get_max_threads()), nrms_(omp_get_max_threads())
    {
        std::istringstream iss(internal_state);
        iss >> this->seed_;
        for(auto& rng : this->rngs_)
        {
            iss >> rng.value;
        }
        for(auto& nrm : this->nrms_)
        {
            iss >> nrm.value;
        }
        if(iss.fail())
        {
            throw_exception<std::runtime_error>("[error] mjolnir::"
                "RandomNumberGenerator<OMP>: parse error in ", internal_state);
        }
    }
    std::string internal_state() const
    {
        std::ostringstream oss;
        oss << this->seed_ << ' ';
        for(const auto& rng : this->rngs_)
        {
            oss << rng.value << ' ';
        }
        for(const auto& nrm : this->nrms_)
        {
            oss << nrm.value << ' ';
        }
        return oss.str();
    }

    real_type uniform_real01()
    {
        auto& rng = this->rngs_.at(omp_get_thread_num()).value;
        return std::generate_canonical<
            real_type, std::numeric_limits<real_type>::digits>(rng);
    }
    real_type uniform_real(const real_type min, const real_type max)
    {
        return this->uniform_real01() * (max - min) + min;
    }

    real_type gaussian()
    {
        const std::size_t thread_id = omp_get_thread_num();
        auto& rng = this->rngs_.at(thread_id).value;
        auto& nrm = this->nrms_.at(thread_id).value;
        return nrm(rng);
    }
    real_type gaussian(const real_type mean, const real_type stddev)
    {
        const std::size_t thread_id = omp_get_thread_num();
        auto& rng = this->rngs_.at(thread_id).value;
        auto& nrm = this->nrms_.at(thread_id).value;
        return nrm(rng) * stddev + mean;
    }

    bool operator==(const RandomNumberGenerator<traits_type>& other) const
    {
        return this->rngs_.size() == other.rngs_.size() && std::equal(
                this->rngs_.begin(), this->rngs_.end(), other.rngs_.begin(),
                [](const rng_type& lhs, const rng_type& rhs) noexcept {
                    return lhs.value == rhs.value;
                }) && this->nrms_.size() == other.nrms_.size() && std::equal(
                this->nrms_.begin(), this->nrms_.end(), other.nrms_.begin(),
                [](const nrm_type& lhs, const nrm_type& rhs) noexcept {
                    return lhs.value == rhs.value;
                });
    }
    bool operator!=(const RandomNumberGenerator<traits_type>& other) const
    {
        return !(*this == other);
    }

  private:
    std::uint32_t seed_;

    std::vector<rng_type, aligned_allocator<rng_type>> rngs_;
    std::vector<nrm_type, aligned_allocator<nrm_type>> nrms_;
};

template<typename realT, template<typename, typename> class boundaryT>
constexpr std::size_t
RandomNumberGenerator<OpenMPSimulatorTraits<realT, boundaryT>>::cache_alignment;

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class RandomNumberGenerator<OpenMPSimulatorTraits<double, UnlimitedBoundary>>;
extern template class RandomNumberGenerator<OpenMPSimulatorTraits<float,  UnlimitedBoundary>>;
extern template class RandomNumberGenerator<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class RandomNumberGenerator<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif

} // mjolnir
#endif /*MJOLNIR_CORE_RANDOM_NUMBER_GENERATOR*/
