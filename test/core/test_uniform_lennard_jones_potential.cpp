#define BOOST_TEST_MODULE "test_uniform_lennard_jones_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/check_potential.hpp>
#include <mjolnir/forcefield/global/UniformLennardJonesPotential.hpp>

BOOST_AUTO_TEST_CASE(UniformLennardJones_double)
{
    using real_type = double;
    using potential_type = mjolnir::UniformLennardJonesPotential<real_type>;
    using parameter_type = potential_type::parameter_type;

    constexpr std::size_t N = 10000;
    constexpr real_type   h = 1e-6;
    constexpr real_type tol = 1e-6;

    const real_type sigma   = 3.0;
    const real_type epsilon = 1.0;

    potential_type potential(/*cutoff = */ 2.5, sigma, epsilon);

    const real_type x_min = 0.8 * sigma;
    const real_type x_max = potential.cutoff_ratio() * sigma;

    mjolnir::test::check_potential(potential, parameter_type{},
                                   x_min, x_max, tol, h, N);
}

BOOST_AUTO_TEST_CASE(UniformLennardJones_float)
{
    using real_type = float;
    using potential_type = mjolnir::UniformLennardJonesPotential<real_type>;
    using parameter_type = potential_type::parameter_type;

    constexpr std::size_t N = 1000;
    constexpr real_type   h = 0.002f;
    constexpr real_type tol = 0.005f;

    const real_type sigma   = 3.0f;
    const real_type epsilon = 1.0f;

    potential_type potential(/*cutoff = */ 2.5f, sigma, epsilon);

    const real_type x_min = 0.8f * sigma;
    const real_type x_max = potential.cutoff_ratio() * sigma;

    mjolnir::test::check_potential(potential, parameter_type{},
                                   x_min, x_max, tol, h, N);
}
