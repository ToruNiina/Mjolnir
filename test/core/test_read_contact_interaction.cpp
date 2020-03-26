#define BOOST_TEST_MODULE "test_read_contact_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/input/read_local_interaction.hpp>
#include <mjolnir/input/read_local_potential.hpp>

BOOST_AUTO_TEST_CASE(read_contact_go_contact)
{
    mjolnir::LoggerManager::set_default_logger("test_read_contact_interaction.log");

    using traits_type = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type   = traits_type::real_type;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction = "Contact"
            potential   = "GoContact"
            topology    = "contact"
            parameters  = [] # empty
        )"_toml;

        const auto base = mjolnir::read_local_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<mjolnir::ContactInteraction<
            traits_type, mjolnir::GoContactPotential<real_type>>*
            >(base.get()); // check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));
    }
}

BOOST_AUTO_TEST_CASE(read_contact_repl_go_contact)
{
    mjolnir::LoggerManager::set_default_logger("test_read_contact_interaction.log");

    using traits_type = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type   = traits_type::real_type;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction = "Contact"
            potential   = "RepulsiveGoContact"
            topology    = "contact"
            parameters  = [] # empty
        )"_toml;

        const auto base = mjolnir::read_local_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<mjolnir::ContactInteraction<
            traits_type, mjolnir::GoContactRepulsivePotential<real_type>>*
            >(base.get()); // check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));
    }
}

BOOST_AUTO_TEST_CASE(read_contact_attr_go_contact)
{
    mjolnir::LoggerManager::set_default_logger("test_read_contact_interaction.log");

    using traits_type = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type   = traits_type::real_type;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction = "Contact"
            potential   = "AttractiveGoContact"
            topology    = "contact"
            parameters  = [] # empty
        )"_toml;

        const auto base = mjolnir::read_local_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<mjolnir::ContactInteraction<
            traits_type, mjolnir::GoContactAttractivePotential<real_type>>*
            >(base.get()); // check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));
    }
}

BOOST_AUTO_TEST_CASE(read_contact_gaussian)
{
    mjolnir::LoggerManager::set_default_logger("test_read_contact_interaction.log");

    using traits_type = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type   = traits_type::real_type;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction = "Contact"
            potential   = "Gaussian"
            topology    = "none"
            parameters  = [] # empty
        )"_toml;

        const auto base = mjolnir::read_local_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<mjolnir::ContactInteraction<
            traits_type, mjolnir::GaussianPotential<real_type>>*
            >(base.get()); // check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));
    }
}

BOOST_AUTO_TEST_CASE(read_contact_multiple_go_contact)
{
    mjolnir::LoggerManager::set_default_logger("test_read_contact_interaction.log");

    using traits_type = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type   = traits_type::real_type;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction = "Contact"
            potential   = "GoContact"
            topology    = "contact"
            parameters  = [
            {indices = [ 0,      1    ], k = 1.0, v0 = 5.0},
            {indices = [[0, 1], [2, 3]], k = 1.0, v0 = 5.0},
            ]
        )"_toml;

        const auto base = mjolnir::read_local_interaction<traits_type>(v);
        BOOST_TEST_REQUIRE(static_cast<bool>(base));

        const auto derv = dynamic_cast<mjolnir::ContactInteraction<
            traits_type, mjolnir::GoContactPotential<real_type>>*
            >(base.get()); // check the expected type is contained
        BOOST_TEST_REQUIRE(static_cast<bool>(derv));

        BOOST_TEST(derv->potentials().size() == 5);
        BOOST_TEST(derv->potentials().at(0).first.at(0) == 0);
        BOOST_TEST(derv->potentials().at(0).first.at(1) == 1);

        BOOST_TEST(derv->potentials().at(1).first.at(0) == 0);
        BOOST_TEST(derv->potentials().at(1).first.at(1) == 2);

        BOOST_TEST(derv->potentials().at(2).first.at(0) == 0);
        BOOST_TEST(derv->potentials().at(2).first.at(1) == 3);

        BOOST_TEST(derv->potentials().at(3).first.at(0) == 1);
        BOOST_TEST(derv->potentials().at(3).first.at(1) == 2);

        BOOST_TEST(derv->potentials().at(4).first.at(0) == 1);
        BOOST_TEST(derv->potentials().at(4).first.at(1) == 3);
    }
}
