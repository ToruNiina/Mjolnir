#define BOOST_TEST_MODULE "test_cell_list"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/util/empty.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/UnlimitedGridCellList.hpp>
#include <mjolnir/core/PeriodicGridCellList.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/Topology.hpp>
#include <random>

template<typename T>
struct dummy_potential
{
    using real_type      = T;
    using parameter_type = mjolnir::empty_t;

    using topology_type        = mjolnir::Topology;
    using molecule_id_type     = typename topology_type::molecule_id_type;
    using connection_kind_type = typename topology_type::connection_kind_type;

    explicit dummy_potential(const real_type cutoff): cutoff_(cutoff) {}

    real_type max_cutoff_length() const noexcept {return this->cutoff_;}

    parameter_type prepare_params(std::size_t, std::size_t) const noexcept
    {
        return parameter_type{};
    }

    bool is_ignored_molecule(std::size_t, std::size_t) const {return false;}

    std::vector<std::pair<connection_kind_type, std::size_t>> ignore_within() const
    {
        return std::vector<std::pair<connection_kind_type, std::size_t>>{};
    }

    std::string name() const {return "dummy potential";}

    real_type cutoff_;
};

BOOST_AUTO_TEST_CASE(test_CellList_UnlimitedBoundary)
{
    mjolnir::LoggerManager::set_default_logger("test_cell_list.log");
    using traits_type     = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type       = typename traits_type::real_type;
    using boundary_type   = typename traits_type::boundary_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using parameter_type  = typename dummy_potential<real_type>::parameter_type;

    constexpr std::size_t N = 1000;
    constexpr double      L = 10.0;
    constexpr double cutoff = 2.0;
    constexpr double margin = 0.25; // threshold = 2.0 * (1+0.25) = 2.5
    constexpr double threshold = cutoff * (1.0 + margin);

    const auto distribute_particle = [](std::mt19937& mt, double l) -> coordinate_type
    {
        return coordinate_type(
            l * std::generate_canonical<real_type, std::numeric_limits<real_type>::digits>(mt),
            l * std::generate_canonical<real_type, std::numeric_limits<real_type>::digits>(mt),
            l * std::generate_canonical<real_type, std::numeric_limits<real_type>::digits>(mt)
        );
    };

    dummy_potential<real_type> pot(cutoff);

    mjolnir::System<traits_type> sys(N, boundary_type{});

    std::mt19937 mt(123456789);
    for(std::size_t i=0; i < N; ++i)
    {
        sys.at(i).mass     = 1.0;
        sys.at(i).position = distribute_particle(mt, L);
    }
    sys.topology().construct_molecules();

    mjolnir::UnlimitedGridCellList<traits_type, parameter_type> vlist(margin);

    using neighbor_type = typename decltype(vlist)::neighbor_type;

    BOOST_TEST(!vlist.valid());

    vlist.initialize(sys, pot);
    vlist.make(sys, pot);
    BOOST_TEST(vlist.valid());

    for(std::size_t i=0; i<N; ++i)
    {
        for(std::size_t j=i+1; j<N; ++j)
        {
            const auto partners = vlist.partners(i);
            if(std::find_if(partners.begin(), partners.end(),
                [=](const neighbor_type& elem) -> bool {return elem.index == j;}
                        ) == partners.end())
            {
                // should be enough distant (>= threshold)
                const auto dist = mjolnir::math::length(sys.adjust_direction(
                            sys.position(j) - sys.position(i)));
                BOOST_TEST(dist >= threshold);
            }
            else
            {
                // should be enough close (< threshold)
                const auto dist = mjolnir::math::length(sys.adjust_direction(
                            sys.position(j) - sys.position(i)));
                BOOST_TEST(dist < threshold);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(test_CellList_PeriodicBoundary)
{
    mjolnir::LoggerManager::set_default_logger("test_cell_list.log");
    using traits_type     = mjolnir::SimulatorTraits<double, mjolnir::CuboidalPeriodicBoundary>;
    using real_type       = typename traits_type::real_type;
    using boundary_type   = typename traits_type::boundary_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using parameter_type  = typename dummy_potential<real_type>::parameter_type;

    constexpr std::size_t N = 1000;
    constexpr double      L = 10.0;
    constexpr double cutoff = 1.0;
    constexpr double margin = 0.5;
    constexpr double threshold = cutoff * (1.0 + margin);

    const auto distribute_particle = [](std::mt19937& mt, double l) -> coordinate_type
    {
        return coordinate_type(
            l * std::generate_canonical<real_type, std::numeric_limits<real_type>::digits>(mt),
            l * std::generate_canonical<real_type, std::numeric_limits<real_type>::digits>(mt),
            l * std::generate_canonical<real_type, std::numeric_limits<real_type>::digits>(mt)
        );
    };

    dummy_potential<real_type> pot(cutoff);

    mjolnir::System<traits_type> sys(N, boundary_type(coordinate_type(0.0, 0.0, 0.0), coordinate_type(L, L, L)));

    std::mt19937 mt(123456789);
    for(std::size_t i=0; i < N; ++i)
    {
        sys.at(i).mass     = 1.0;
        sys.at(i).position = distribute_particle(mt, L);
    }
    sys.topology().construct_molecules();

    mjolnir::PeriodicGridCellList<traits_type, parameter_type> vlist(margin);
    using neighbor_type = typename decltype(vlist)::neighbor_type;

    BOOST_TEST(!vlist.valid());

    vlist.initialize(sys, pot);
    vlist.make(sys, pot);
    BOOST_TEST(vlist.valid());

    for(std::size_t i=0; i<N; ++i)
    {
        for(std::size_t j=i+1; j<N; ++j)
        {
            const auto partners = vlist.partners(i);
            if(std::find_if(partners.begin(), partners.end(),
                    [=](const neighbor_type& elem){return elem.index == j;}
                        ) == partners.end())
            {
                // should be enough distant (>= threshold)
                const auto dist = mjolnir::math::length(sys.adjust_direction(
                            sys.position(j) - sys.position(i)));
                BOOST_TEST(dist >= threshold);
            }
            else
            {
                // should be enough close (< threshold)
                const auto dist = mjolnir::math::length(sys.adjust_direction(
                            sys.position(j) - sys.position(i)));
                BOOST_TEST(dist < threshold);
            }
        }
    }
}
