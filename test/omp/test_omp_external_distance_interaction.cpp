#define BOOST_TEST_MODULE "test_omp_external_distance_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/math/math.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/AxisAlignedPlane.hpp>
#include <mjolnir/forcefield/external/LennardJonesWallPotential.hpp>
#include <mjolnir/omp/System.hpp>
#include <mjolnir/omp/RandomNumberGenerator.hpp>
#include <mjolnir/omp/UnlimitedGridCellList.hpp>
#include <mjolnir/omp/ExternalDistanceInteraction.hpp>
#include <mjolnir/util/make_unique.hpp>

BOOST_AUTO_TEST_CASE(omp_ExternalDistacne_calc_force)
{
    constexpr double tol = 1e-8;
    mjolnir::LoggerManager::set_default_logger("test_omp_external_distance_interaction.log");

    using traits_type      = mjolnir::OpenMPSimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type        = typename traits_type::real_type;
    using coordinate_type  = typename traits_type::coordinate_type;
    using boundary_type    = typename traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using topology_type    = mjolnir::Topology;
    using potential_type   = mjolnir::LennardJonesWallPotential<real_type>;
    using parameter_type   = typename potential_type::parameter_type;
    using shape_type       = mjolnir::AxisAlignedPlane<traits_type, mjolnir::PositiveZDirection<traits_type>>;
    using interaction_type = mjolnir::ExternalDistanceInteraction<traits_type, potential_type, shape_type>;
    using rng_type         = mjolnir::RandomNumberGenerator<traits_type>;

    using sequencial_system_type      = mjolnir::System<
        mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>>;
    using sequencial_shape_type       = mjolnir::AxisAlignedPlane<
        mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>,
        mjolnir::PositiveZDirection<traits_type>>;
    using sequencial_interaction_type = mjolnir::ExternalDistanceInteraction<
        mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>,
        potential_type, sequencial_shape_type>;

    const int max_number_of_threads = omp_get_max_threads();
    BOOST_TEST_WARN(max_number_of_threads > 2);
    BOOST_TEST_MESSAGE("maximum number of threads = " << max_number_of_threads);


    const std::size_t N_particle = 64;
    std::vector<std::pair<std::size_t, parameter_type>> parameters(N_particle);
    for(std::size_t i=0; i<N_particle; ++i)
    {
        parameters.emplace_back(i, parameter_type(1.0, 1.0));
    }
    const potential_type potential(2.5, parameters);

    for(int num_thread=1; num_thread<=max_number_of_threads; ++num_thread)
    {
        omp_set_num_threads(num_thread);
        BOOST_TEST_MESSAGE("maximum number of threads = " << omp_get_max_threads());


        rng_type    rng(123456789);
        system_type sys(N_particle, boundary_type{});
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            const auto i_x = i % 4;
            const auto i_y = i / 4;
            const auto i_z = i / 16;

            sys.mass(i)     = 1.0;
            sys.position(i) = mjolnir::math::make_coordinate<coordinate_type>(i_x*2.0, i_y*2.0, i_z*2.0);
            sys.velocity(i) = mjolnir::math::make_coordinate<coordinate_type>(0, 0, 0);
            sys.force(i)    = mjolnir::math::make_coordinate<coordinate_type>(0, 0, 0);
            sys.name(i)     = "X";
            sys.group(i)    = "TEST";
        }

        topology_type topol(N_particle);

        // add perturbation
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            mjolnir::math::X(sys.position(i)) += rng.uniform_real(-0.1, 0.1);
            mjolnir::math::Y(sys.position(i)) += rng.uniform_real(-0.1, 0.1);
            mjolnir::math::Z(sys.position(i)) += rng.uniform_real(-0.1, 0.1);
        }

        // init sequential one with the same coordinates
        sequencial_system_type seq_sys(N_particle, boundary_type{});
        assert(sys.size() == seq_sys.size());
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            seq_sys.mass(i)     = sys.mass(i);
            seq_sys.position(i) = sys.position(i);
            seq_sys.velocity(i) = sys.velocity(i);
            seq_sys.force(i)    = sys.force(i);
            seq_sys.name(i)     = sys.name(i);
            seq_sys.group(i)    = sys.group(i);
        }

        shape_type            xyplane    (0.0);
        sequencial_shape_type seq_xyplane(0.0);

        topol.construct_molecules();

        xyplane    .initialize(sys,     potential);
        seq_xyplane.initialize(seq_sys, potential);

        interaction_type            interaction(
                std::move(xyplane), potential_type(potential));
        sequencial_interaction_type seq_interaction(
                std::move(seq_xyplane), potential_type(potential));

        interaction    .initialize(sys);
        seq_interaction.initialize(seq_sys);

        // calculate forces with openmp
        interaction.calc_force(sys);
        sys.postprocess_forces();

        // calculate forces without openmp
        seq_interaction.calc_force(seq_sys);

        // check the values are the same
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            BOOST_TEST(mjolnir::math::X(seq_sys.force(i)) == mjolnir::math::X(sys.force(i)),
                       boost::test_tools::tolerance(tol));
            BOOST_TEST(mjolnir::math::Y(seq_sys.force(i)) == mjolnir::math::Y(sys.force(i)),
                       boost::test_tools::tolerance(tol));
            BOOST_TEST(mjolnir::math::Z(seq_sys.force(i)) == mjolnir::math::Z(sys.force(i)),
                       boost::test_tools::tolerance(tol));
        }
        BOOST_TEST(interaction.calc_energy(sys) == seq_interaction.calc_energy(seq_sys),
                   boost::test_tools::tolerance(tol));
    }
}

BOOST_AUTO_TEST_CASE(omp_ExternalDistance_calc_force_and_energy)
{
    mjolnir::LoggerManager::set_default_logger(
            "test_omp_external_distance_interaction.log");

    using traits_type      = mjolnir::OpenMPSimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type        = traits_type::real_type;
    using coordinate_type  = traits_type::coordinate_type;
    using boundary_type    = traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using potential_type   = mjolnir::LennardJonesWallPotential<real_type>;
    using shape_type       = mjolnir::AxisAlignedPlane<traits_type, mjolnir::PositiveZDirection<traits_type>>;
    using interaction_type = mjolnir::ExternalDistanceInteraction<traits_type, potential_type, shape_type>;

    interaction_type interaction(shape_type(0.0, 0.5),
        potential_type(/* cutoff = */2.0, {
            {0, {1.0, 1.0}}, {1, {1.0, 1.0}}
        }));

    system_type sys(2, boundary_type{});

    sys.at(0).mass  = 1.0;
    sys.at(1).mass  = 1.0;
    sys.at(0).rmass = 1.0;
    sys.at(1).rmass = 1.0;

    sys.at(0).position = coordinate_type( 0.0, 0.0, 1.0);
    sys.at(1).position = coordinate_type( 0.0, 0.0, 1.0);
    sys.at(0).velocity = coordinate_type( 0.0, 0.0, 0.0);
    sys.at(1).velocity = coordinate_type( 0.0, 0.0, 0.0);
    sys.at(0).force    = coordinate_type( 0.0, 0.0, 0.0);
    sys.at(1).force    = coordinate_type( 0.0, 0.0, 0.0);

    sys.at(0).name  = "X";
    sys.at(1).name  = "X";
    sys.at(0).group = "NONE";
    sys.at(1).group = "NONE";

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);
    std::normal_distribution<real_type> gauss(0.0, 1.0);

    for(int i = 0; i < 10000; ++i)
    {
        sys.at(0).position = coordinate_type(0.0, 0.0, 1.0);
        sys.at(1).position = coordinate_type(0.0, 0.0, 1.0);

        // move particles a bit, randomly. and reset forces.
        for(std::size_t idx=0; idx<sys.size(); ++idx)
        {
            sys.position(idx) += coordinate_type(0.01 * uni(mt), 0.01 * uni(mt), 0.01 * uni(mt));
            sys.force(idx)     = coordinate_type(0.0, 0.0, 0.0);
        }

        system_type ref_sys = sys;

        constexpr real_type tol = 1e-4;

        const auto energy = interaction.calc_force_and_energy(sys);
        const auto ref_energy = interaction.calc_energy(ref_sys);
        interaction.calc_force(ref_sys);
        BOOST_TEST(ref_energy == energy, boost::test_tools::tolerance(tol));

        for(std::size_t idx=0; idx<sys.size(); ++idx)
        {
            BOOST_TEST(mjolnir::math::X(sys.force(idx)) == mjolnir::math::X(ref_sys.force(idx)), boost::test_tools::tolerance(tol));
            BOOST_TEST(mjolnir::math::Y(sys.force(idx)) == mjolnir::math::Y(ref_sys.force(idx)), boost::test_tools::tolerance(tol));
            BOOST_TEST(mjolnir::math::Z(sys.force(idx)) == mjolnir::math::Z(ref_sys.force(idx)), boost::test_tools::tolerance(tol));
        }
    }
}
