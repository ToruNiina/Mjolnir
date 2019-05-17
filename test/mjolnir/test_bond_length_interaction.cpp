#define BOOST_TEST_MODULE "test_bond_length_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/interaction/local/BondLengthInteraction.hpp>
#include <mjolnir/potential/local/HarmonicPotential.hpp>
#include <mjolnir/util/make_unique.hpp>

BOOST_AUTO_TEST_CASE(BondLength_calc_force)
{
    using traits_type      = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type        = traits_type::real_type;
    using coord_type       = traits_type::coordinate_type;
    using boundary_type    = traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using harmonic_type    = mjolnir::HarmonicPotential<real_type>;
    using bond_length_type = mjolnir::BondLengthInteraction<traits_type, harmonic_type>;

    constexpr real_type tol = 1e-8;

    auto normalize = [](const coord_type& v){return v / mjolnir::math::length(v);};

    const real_type k(100.);
    const real_type native(2.0);

    harmonic_type    potential(k, native);
    bond_length_type interaction("none", {{ {{0,1}}, potential}});

    system_type sys(2, boundary_type{});

    sys.at(0).mass = 1.0;
    sys.at(1).mass = 1.0;
    sys.at(0).rmass = 1.0;
    sys.at(1).rmass = 1.0;

    sys.at(0).position = coord_type(0,0,0);
    sys.at(1).position = coord_type(0,0,0);
    sys.at(0).velocity = coord_type(0,0,0);
    sys.at(1).velocity = coord_type(0,0,0);
    sys.at(0).force    = coord_type(0,0,0);
    sys.at(1).force    = coord_type(0,0,0);

    sys.at(0).name  = "X";
    sys.at(1).name  = "X";
    sys.at(0).group = "NONE";
    sys.at(1).group = "NONE";

    const real_type dr = 1e-3;
    real_type dist = 1e0;
    for(int i = 0; i < 2000; ++i)
    {
        sys[0].position = coord_type(0,0,0);
        sys[1].position = coord_type(0,0,0);
        sys[0].force    = coord_type(0,0,0);
        sys[1].force    = coord_type(0,0,0);
        sys[1].position[0] = dist;

        const real_type deriv = potential.derivative(dist);
        const real_type coef  = std::abs(deriv);

        interaction.calc_force(sys);

        const real_type force_strength1 = mjolnir::math::length(sys[0].force);
        const real_type force_strength2 = mjolnir::math::length(sys[1].force);


        // direction
        if(i == 1000) // most stable point
        {
            BOOST_TEST(force_strength1 == 0.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(force_strength2 == 0.0, boost::test_tools::tolerance(tol));
        }
        else if(i < 1000) // repulsive
        {
            BOOST_TEST(coef == force_strength1, boost::test_tools::tolerance(tol));
            BOOST_TEST(coef == force_strength2, boost::test_tools::tolerance(tol));

            const real_type dir1 = mjolnir::math::dot_product(
                normalize(sys[0].force), normalize(sys[0].position - sys[1].position));
            const real_type dir2 = mjolnir::math::dot_product(
                normalize(sys[1].force), normalize(sys[1].position - sys[0].position));

            BOOST_TEST(dir1 == 1.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(dir2 == 1.0, boost::test_tools::tolerance(tol));
        }
        else if(i > 1000) // attractive
        {
            BOOST_TEST(coef == force_strength1, boost::test_tools::tolerance(tol));
            BOOST_TEST(coef == force_strength2, boost::test_tools::tolerance(tol));

            const real_type dir1 = mjolnir::math::dot_product(
                normalize(sys[0].force), normalize(sys[1].position - sys[0].position));
            const real_type dir2 = mjolnir::math::dot_product(
                normalize(sys[1].force), normalize(sys[0].position - sys[1].position));

            BOOST_TEST(dir1 == 1e0, boost::test_tools::tolerance(tol));
            BOOST_TEST(dir2 == 1e0, boost::test_tools::tolerance(tol));
        }
        BOOST_TEST(mjolnir::math::length(sys[0].force + sys[1].force) == 0.0,
                   boost::test_tools::tolerance(tol));

        dist += dr;
    }
}
