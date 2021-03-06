#define BOOST_TEST_MODULE "test_repulsive_gocontact_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/check_potential.hpp>
#include <mjolnir/forcefield/local/GoContactRepulsivePotential.hpp>
#include <iomanip>

BOOST_AUTO_TEST_CASE(GoContactRepulsive_double)
{
    using real_type = double;
    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-6;
    constexpr real_type tol = 1e-6;

    const real_type e  = 1.0;
    const real_type r0 = 5.0;

    mjolnir::GoContactRepulsivePotential<real_type> pot(e, r0);

    {
        const real_type x_min = 0.8 * r0;
        const real_type x_max = pot.cutoff();
        mjolnir::test::check_potential(pot, x_min, x_max, tol, h, N);
    }
    {
        const real_type x_min = pot.cutoff() + 1.01;
        const real_type x_max = 2 * pot.cutoff();
        const real_type dx = (x_max - x_min) / N;
        for(std::size_t i=0; i<N; ++i)
        {
            const real_type x = x_min + i * dx;
            BOOST_TEST(pot.potential(x)  == 0.0, boost::test_tools::tolerance(h));
            BOOST_TEST(pot.derivative(x) == 0.0, boost::test_tools::tolerance(h));
        }
    }
}

BOOST_AUTO_TEST_CASE(GoContactRepulsive_float)
{
    using real_type = float;
    constexpr std::size_t N = 100;
    constexpr real_type   h = 1e-2;
    constexpr real_type tol = 1e-2;

    const real_type e  = 1.0;
    const real_type r0 = 5.0;

    mjolnir::GoContactRepulsivePotential<real_type> pot(e, r0);
    {
        const real_type x_min = 0.8f * r0;
        const real_type x_max = pot.cutoff();

        mjolnir::test::check_potential(pot, x_min, x_max, tol, h, N);
    }
    {
        const real_type x_min = pot.cutoff() + 1.01;
        const real_type x_max = 2 * pot.cutoff();
        const real_type dx = (x_max - x_min) / N;
        for(std::size_t i=0; i<N; ++i)
        {
            const real_type x = x_min + i * dx;
            BOOST_TEST(pot.potential(x)  == 0.0, boost::test_tools::tolerance(h));
            BOOST_TEST(pot.derivative(x) == 0.0, boost::test_tools::tolerance(h));
        }
    }
}
