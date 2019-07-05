#define BOOST_TEST_MODULE "test_read_lennard_jones_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/input/read_global_potential.hpp>
#include <tuple>

using test_types = std::tuple<double, float>;

constexpr inline float  tolerance_value(float)  noexcept {return 1e-4;}
constexpr inline double tolerance_value(double) noexcept {return 1e-8;}

template<typename Real>
decltype(boost::test_tools::tolerance(std::declval<Real>()))
tolerance() {return boost::test_tools::tolerance(tolerance_value(Real()));}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_lennard_jones_noenv, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_lennard_jones.log");

    using real_type = T;
    {
        using namespace toml::literals;
        const auto v = u8R"(
            interaction = "Pair"
            potential   = "LennardJones"
            spatial_partition.type  = "Naive"
            ignore.molecule         = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            parameters = [
                {index =   0, sigma =   2.0, epsilon = 1.5},
                {index =   1, sigma =   2.0, "ε"     = 1.5},
                {index =   2, "σ"   =   2.0, epsilon = 1.5},
                {index =   3, "σ"   =   2.0, "ε"     = 1.5},
                {index =   5, sigma =   5.0, epsilon = 0.5},
                {index =   7, sigma =   7.0, epsilon = 0.7},
                {index = 100, sigma = 100.0, epsilon = 0.1},
            ]
        )"_toml;

        const auto g = mjolnir::read_lennard_jones_potential<real_type>(v);

        const auto ignore_within = g.ignore_within();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(g.ignore_within().size() == 2u);
        BOOST_TEST(within.at("bond")    == 3ul);
        BOOST_TEST(within.at("contact") == 1ul);

        BOOST_TEST(g.participants().size() ==   7u);
        BOOST_TEST(g.participants().at(0)  ==   0u);
        BOOST_TEST(g.participants().at(1)  ==   1u);
        BOOST_TEST(g.participants().at(2)  ==   2u);
        BOOST_TEST(g.participants().at(3)  ==   3u);
        BOOST_TEST(g.participants().at(4)  ==   5u);
        BOOST_TEST(g.participants().at(5)  ==   7u);
        BOOST_TEST(g.participants().at(6)  == 100u);

        BOOST_TEST(g.parameters().at(  0).first  == real_type(  2.0), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(  1).first  == real_type(  2.0), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(  2).first  == real_type(  2.0), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(  3).first  == real_type(  2.0), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(  5).first  == real_type(  5.0), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(  7).first  == real_type(  7.0), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(100).first  == real_type(100.0), tolerance<real_type>());

        BOOST_TEST(g.parameters().at(  0).second == real_type(1.5), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(  1).second == real_type(1.5), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(  2).second == real_type(1.5), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(  3).second == real_type(1.5), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(  5).second == real_type(0.5), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(  7).second == real_type(0.7), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(100).second == real_type(0.1), tolerance<real_type>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_lennard_jones_env, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_lennard_jones.log");

    using real_type = T;
    {
        using namespace toml::literals;
        const auto v = u8R"(
            interaction = "Pair"
            potential   = "LennardJones"
            spatial_partition.type  = "Naive"
            ignore.molecule         = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            env.100   = 100.0
            env.small = 0.1
            env.half  = 0.5
            parameters = [
                {index =   0, sigma =   2.0, epsilon = 1.5},
                {index =   1, sigma =   2.0, "ε"     = 1.5},
                {index =   2, "σ"   =   2.0, epsilon = 1.5},
                {index =   3, "σ"   =   2.0, "ε"     = 1.5},
                {index =   5, sigma =   5.0, epsilon = "half"},
                {index =   7, sigma =   7.0, epsilon = 0.7},
                {index = 100, sigma = "100", epsilon = "small"},
            ]
        )"_toml;

        const auto g = mjolnir::read_lennard_jones_potential<real_type>(v);

        const auto ignore_within = g.ignore_within();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(g.ignore_within().size() == 2u);
        BOOST_TEST(within.at("bond")    == 3ul);
        BOOST_TEST(within.at("contact") == 1ul);

        BOOST_TEST(g.participants().size() ==   7u);
        BOOST_TEST(g.participants().at(0)  ==   0u);
        BOOST_TEST(g.participants().at(1)  ==   1u);
        BOOST_TEST(g.participants().at(2)  ==   2u);
        BOOST_TEST(g.participants().at(3)  ==   3u);
        BOOST_TEST(g.participants().at(4)  ==   5u);
        BOOST_TEST(g.participants().at(5)  ==   7u);
        BOOST_TEST(g.participants().at(6)  == 100u);

        BOOST_TEST(g.parameters().at(  0).first  == real_type(  2.0), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(  1).first  == real_type(  2.0), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(  2).first  == real_type(  2.0), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(  3).first  == real_type(  2.0), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(  5).first  == real_type(  5.0), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(  7).first  == real_type(  7.0), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(100).first  == real_type(100.0), tolerance<real_type>());

        BOOST_TEST(g.parameters().at(  0).second == real_type(1.5), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(  1).second == real_type(1.5), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(  2).second == real_type(1.5), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(  3).second == real_type(1.5), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(  5).second == real_type(0.5), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(  7).second == real_type(0.7), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(100).second == real_type(0.1), tolerance<real_type>());
    }
}