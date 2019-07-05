#define BOOST_TEST_MODULE "test_read_uniform_lennard_jones_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif
#include <mjolnir/input/read_global_potential.hpp>

#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <tuple>

using test_types = std::tuple<double, float>;

constexpr inline float  tolerance_value(float)  noexcept {return 1e-4;}
constexpr inline double tolerance_value(double) noexcept {return 1e-8;}

template<typename Real>
decltype(boost::test_tools::tolerance(std::declval<Real>()))
tolerance() {return boost::test_tools::tolerance(tolerance_value(Real()));}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_uniform_lennard_jones_noenv, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_uniform_lennard_jones.log");
    using real_type = T;

    // a dummy system for testing `initialize` method
    using traits_type   = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    using boundary_type = typename traits_type::boundary_type;
    mjolnir::System<traits_type> sys(10, boundary_type{});
    {
        using namespace toml::literals;
        const auto v = u8R"(
            interaction = "Pair"
            potential   = "UniformLennardJones"
            spatial_partition.type  = "CellList"
            ignore.molecule         = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            sigma   = 2.0
            epsilon = 1.5
        )"_toml;

        auto g = mjolnir::read_uniform_lennard_jones_potential<real_type>(v);

        const auto ignore_within = g.ignore_within();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(g.ignore_within().size() == 2u);
        BOOST_TEST(within.at("bond")    == 3ul);
        BOOST_TEST(within.at("contact") == 1ul);
        BOOST_TEST(g.sigma()   == real_type(2.0), tolerance<real_type>());
        BOOST_TEST(g.epsilon() == real_type(1.5), tolerance<real_type>());
        BOOST_TEST(g.participants().empty());

        g.initialize(sys);
        BOOST_TEST(g.participants().size() == 10u);
        BOOST_TEST(g.participants().at(0) == 0u);
        BOOST_TEST(g.participants().at(1) == 1u);
        BOOST_TEST(g.participants().at(2) == 2u);
        BOOST_TEST(g.participants().at(3) == 3u);
        BOOST_TEST(g.participants().at(4) == 4u);
        BOOST_TEST(g.participants().at(5) == 5u);
        BOOST_TEST(g.participants().at(6) == 6u);
        BOOST_TEST(g.participants().at(7) == 7u);
        BOOST_TEST(g.participants().at(8) == 8u);
        BOOST_TEST(g.participants().at(9) == 9u);
    }
    {
        using namespace toml::literals;
        const auto v = u8R"(
            interaction = "Pair"
            potential   = "UniformLennardJones"
            spatial_partition.type  = "CellList"
            ignore.molecule         = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            "σ" = 2.0
            "ε" = 1.5
        )"_toml;

        auto g = mjolnir::read_uniform_lennard_jones_potential<real_type>(v);

        const auto ignore_within = g.ignore_within();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(g.ignore_within().size() == 2u);
        BOOST_TEST(within.at("bond")    == 3ul);
        BOOST_TEST(within.at("contact") == 1ul);
        BOOST_TEST(g.sigma()   == real_type(2.0), tolerance<real_type>());
        BOOST_TEST(g.epsilon() == real_type(1.5), tolerance<real_type>());
        BOOST_TEST(g.participants().empty());

        g.initialize(sys);
        BOOST_TEST(g.participants().size() == 10u);
        BOOST_TEST(g.participants().at(0) == 0u);
        BOOST_TEST(g.participants().at(1) == 1u);
        BOOST_TEST(g.participants().at(2) == 2u);
        BOOST_TEST(g.participants().at(3) == 3u);
        BOOST_TEST(g.participants().at(4) == 4u);
        BOOST_TEST(g.participants().at(5) == 5u);
        BOOST_TEST(g.participants().at(6) == 6u);
        BOOST_TEST(g.participants().at(7) == 7u);
        BOOST_TEST(g.participants().at(8) == 8u);
        BOOST_TEST(g.participants().at(9) == 9u);
    }
    {
        using namespace toml::literals;
        const auto v = u8R"(
            interaction = "Pair"
            potential   = "UniformLennardJones"
            spatial_partition.type  = "CellList"
            ignore.molecule         = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            "σ" = 2.0
            "ε" = 1.5
            parameters = [
                {index = 1},
                {index = 2},
                {index = 3},
                {index = 4},
                {index = 5},
            ]
        )"_toml;

        auto g = mjolnir::read_uniform_lennard_jones_potential<real_type>(v);

        const auto ignore_within = g.ignore_within();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(g.ignore_within().size() == 2u);
        BOOST_TEST(within.at("bond")    == 3ul);
        BOOST_TEST(within.at("contact") == 1ul);
        BOOST_TEST(g.sigma()   == real_type(2.0), tolerance<real_type>());
        BOOST_TEST(g.epsilon() == real_type(1.5), tolerance<real_type>());
        BOOST_TEST(g.participants().size() == 5u);

        g.initialize(sys);
        BOOST_TEST(g.participants().size() == 5u);
        BOOST_TEST(g.participants().at(0)  == 1u);
        BOOST_TEST(g.participants().at(1)  == 2u);
        BOOST_TEST(g.participants().at(2)  == 3u);
        BOOST_TEST(g.participants().at(3)  == 4u);
        BOOST_TEST(g.participants().at(4)  == 5u);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_uniform_lennard_jones_env, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_uniform_lennard_jones.log");
    using real_type = T;

    // a dummy system for testing `initialize` method
    using traits_type   = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    using boundary_type = typename traits_type::boundary_type;
    mjolnir::System<traits_type> sys(10, boundary_type{});

    {
        using namespace toml::literals;
        const auto v = u8R"(
            interaction = "Pair"
            potential   = "UniformLennardJones"
            spatial_partition.type  = "CellList"
            ignore.molecule         = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            "σ" = 2.0
            "ε" = 1.5
            env.three = 3
            parameters = [
                {index = 1},
                {index = 2},
                {index = "three"},
                {index = 4},
                {index = 5},
            ]
        )"_toml;

        auto g = mjolnir::read_uniform_lennard_jones_potential<real_type>(v);

        const auto ignore_within = g.ignore_within();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(g.ignore_within().size() == 2u);
        BOOST_TEST(within.at("bond")    == 3ul);
        BOOST_TEST(within.at("contact") == 1ul);
        BOOST_TEST(g.sigma()   == real_type(2.0), tolerance<real_type>());
        BOOST_TEST(g.epsilon() == real_type(1.5), tolerance<real_type>());
        BOOST_TEST(g.participants().size() == 5u);

        g.initialize(sys);
        BOOST_TEST(g.participants().size() == 5u);
        BOOST_TEST(g.participants().at(0)  == 1u);
        BOOST_TEST(g.participants().at(1)  == 2u);
        BOOST_TEST(g.participants().at(2)  == 3u);
        BOOST_TEST(g.participants().at(3)  == 4u);
        BOOST_TEST(g.participants().at(4)  == 5u);
    }
}