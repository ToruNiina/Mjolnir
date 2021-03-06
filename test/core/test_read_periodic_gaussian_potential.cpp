#define BOOST_TEST_MODULE "test_read_periodic_gaussian_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/input/read_local_potential.hpp>

BOOST_AUTO_TEST_CASE(read_periodic_gaussian_double)
{
    mjolnir::LoggerManager::set_default_logger("test_read_periodic_gaussian.log");
    using real_type = double;
    using namespace toml::literals;

    constexpr real_type tol = 1e-8;
    {
        const toml::value env; // empty env
        const auto v = u8R"(
            indices = [1, 2]
            k       = 3.14
            sigma   = 0.577
            v0      = 2.71
        )"_toml;

        const auto g = mjolnir::read_periodic_gaussian_potential<real_type>(v, env);
        BOOST_TEST(g.k()     == 3.14,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.sigma() == 0.577, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.v0()    == 2.71,  boost::test_tools::tolerance(tol));
    }

    {
        const toml::value env; // empty env
        const auto v = u8R"(
            indices = [1, 2]
            k       = 3.14
            "σ"     = 0.577
            v0      = 2.71
        )"_toml;

        const auto g = mjolnir::read_periodic_gaussian_potential<real_type>(v, env);
        BOOST_TEST(g.k()     == 3.14,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.sigma() == 0.577, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.v0()    == 2.71,  boost::test_tools::tolerance(tol));
    }

    {
        const toml::value env = u8R"(
            indices = [1, 2]
            k       = 3.14
            sigma   = 0.577
            v0      = 2.71
        )"_toml;

        const auto v = u8R"(
            indices = "indices"
            k       = "k"
            sigma   = "sigma"
            v0      = "v0"
        )"_toml;

        const auto g = mjolnir::read_periodic_gaussian_potential<real_type>(v, env);
        BOOST_TEST(g.k()     == 3.14,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.sigma() == 0.577, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.v0()    == 2.71,  boost::test_tools::tolerance(tol));
    }
}

BOOST_AUTO_TEST_CASE(read_periodic_gaussian_float)
{
    mjolnir::LoggerManager::set_default_logger("test_read_periodic_gaussian.log");
    using real_type = float;
using namespace toml::literals;
    constexpr real_type tol = 1e-4;

    {
        const toml::value env; // empty env
        const auto v = u8R"(
            indices = [1, 2]
            k       = 3.14
            sigma   = 0.577
            v0      = 2.71
        )"_toml;
        const auto g = mjolnir::read_periodic_gaussian_potential<real_type>(v, env);
        BOOST_TEST(g.k()     == 3.14f,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.sigma() == 0.577f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.v0()    == 2.71f,  boost::test_tools::tolerance(tol));
    }

    {
        const toml::value env; // empty env
        const auto v = u8R"(
            indices = [1, 2]
            k       = 3.14
            "σ"     = 0.577
            v0      = 2.71
        )"_toml;

        const auto g = mjolnir::read_periodic_gaussian_potential<real_type>(v, env);
        BOOST_TEST(g.k()     == 3.14f,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.sigma() == 0.577f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.v0()    == 2.71f,  boost::test_tools::tolerance(tol));
    }

    {
        const toml::value env = u8R"(
            indices = [1, 2]
            k       = 3.14
            sigma   = 0.577
            v0      = 2.71
        )"_toml;

        const auto v = u8R"(
            indices = "indices"
            k       = "k"
            sigma   = "sigma"
            v0      = "v0"
        )"_toml;

        const auto g = mjolnir::read_periodic_gaussian_potential<real_type>(v, env);
        BOOST_TEST(g.k()     == 3.14f,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.sigma() == 0.577f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.v0()    == 2.71f,  boost::test_tools::tolerance(tol));
    }
}

// ---------------------------------------------------------------------------
// read_local_potential

BOOST_AUTO_TEST_CASE(read_local_potential_periodic_gaussian_double)
{
    mjolnir::LoggerManager::set_default_logger("test_read_periodic_gaussian.log");

    using real_type = double;
    constexpr real_type tol = 1e-8;
    {
        using namespace toml::literals;
        const auto v = u8R"(
            parameters = [
                {indices = [1, 2], k = 3.14, sigma = 0.577, v0 = 2.71}
            ]
        )"_toml;

        const auto g = mjolnir::read_local_potential<2,
              mjolnir::PeriodicGaussianPotential<real_type>>(v);

        const std::array<std::size_t, 2> ref_idx{{1, 2}};

        BOOST_TEST(g.size() == 1u);
        BOOST_TEST(g.at(0).first == ref_idx);
        BOOST_TEST(g.at(0).second.k()     == 3.14,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.at(0).second.sigma() == 0.577, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.at(0).second.v0()    == 2.71,  boost::test_tools::tolerance(tol));
    }
    {
        using namespace toml::literals;
        const auto v = u8R"(
            env.pi = 3.14
            parameters = [
                {indices = [1, 2], k = "pi", sigma = 0.577, v0 = 2.71}
            ]
        )"_toml;

        const auto g = mjolnir::read_local_potential<2,
              mjolnir::PeriodicGaussianPotential<real_type>>(v);

        const std::array<std::size_t, 2> ref_idx{{1, 2}};

        BOOST_TEST(g.size() == 1u);
        BOOST_TEST(g.at(0).first == ref_idx);
        BOOST_TEST(g.at(0).second.k()     == 3.14,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.at(0).second.sigma() == 0.577, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.at(0).second.v0()    == 2.71,  boost::test_tools::tolerance(tol));
    }

}

BOOST_AUTO_TEST_CASE(read_local_potential_periodic_gaussian_float)
{
    mjolnir::LoggerManager::set_default_logger("test_read_periodic_gaussian.log");

    using real_type = float;
    constexpr real_type tol = 1e-4;
    {
        using namespace toml::literals;
        const auto v = u8R"(
            parameters = [
                {indices = [1, 2], k = 3.14, sigma = 0.577, v0 = 2.71}
            ]
        )"_toml;

        const auto g = mjolnir::read_local_potential<2,
              mjolnir::PeriodicGaussianPotential<real_type>>(v);

        const std::array<std::size_t, 2> ref_idx{{1, 2}};

        BOOST_TEST(g.size() == 1u);
        BOOST_TEST(g.at(0).first == ref_idx);
        BOOST_TEST(g.at(0).second.k()     == 3.14f,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.at(0).second.sigma() == 0.577f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.at(0).second.v0()    == 2.71f,  boost::test_tools::tolerance(tol));
    }
}
