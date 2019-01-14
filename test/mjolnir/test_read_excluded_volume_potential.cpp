#define BOOST_TEST_MODULE "test_read_excluded_volume_potential"

#include <boost/test/included/unit_test.hpp>
#include <mjolnir/input/read_global_potential.hpp>

BOOST_AUTO_TEST_CASE(read_excluded_volume_double)
{
    mjolnir::LoggerManager::set_default_logger("test_read_excluded_volume.log");

    using real_type = double;
    constexpr real_type tol = 1e-8;
    {
        const toml::value v = toml::table{
            {"interaction",       toml::value("Pair")},
            {"potential",         toml::value("ExcludedVolume")},
            {"spatial_partition", toml::value(toml::table{
                        {"type", toml::value("Nothing")}
            })},
            {"epsilon",           toml::value(3.14)},
            {"ignore",            toml::value(toml::table{
                {"molecule",         toml::value("Nothing")},
                {"particles_within", toml::table{{"bond", 3}, {"contact", 1}}},
            })},
            {"parameters",        toml::value(toml::array{
                toml::table{{"index", 0}, {"radius", 2.0}},
                toml::table{{"index", 1}, {"radius", 2.0}}
            })}
        };
        const auto g = mjolnir::read_excluded_volume_potential<real_type>(v);

        const auto ignore_within = g.ignore_within();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(g.ignore_within().size() == 2);
        BOOST_TEST(within.at("bond")    == 3ul);
        BOOST_TEST(within.at("contact") == 1ul);
        BOOST_TEST(g.parameters().size() == 2);
        BOOST_TEST(g.parameters().at(0)  == 2.0,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.parameters().at(1)  == 2.0,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.epsilon()           == 3.14, boost::test_tools::tolerance(tol));
    }
}

BOOST_AUTO_TEST_CASE(read_excluded_volume_float)
{
    mjolnir::LoggerManager::set_default_logger("test_read_excluded_volume.log");
    using real_type = float;
    constexpr real_type tol = 1e-4;

    {
        const toml::value v = toml::table{
            {"interaction",       toml::value("Pair")},
            {"potential",         toml::value("ExcludedVolume")},
            {"spatial_partition", toml::value(toml::table{
                        {"type", toml::value("Nothing")}
            })},
            {"epsilon",           toml::value(3.14)},
            {"ignore",            toml::value(toml::table{
                {"molecule",         toml::value("Nothing")},
                {"particles_within", toml::table{{"bond", 3}, {"contact", 1}}},
            })},
            {"parameters",        toml::value(toml::array{
                toml::table{{"index", 0}, {"radius", 2.0}},
                toml::table{{"index", 1}, {"radius", 2.0}}
            })}
        };
        const auto g = mjolnir::read_excluded_volume_potential<real_type>(v);

        const auto ignore_within = g.ignore_within();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(g.ignore_within().size() == 2);
        BOOST_TEST(within.at("bond")    == 3ul);
        BOOST_TEST(within.at("contact") == 1ul);
        BOOST_TEST(g.parameters().size() == 2);
        BOOST_TEST(g.parameters().at(0)  == 2.0f,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.parameters().at(1)  == 2.0f,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.epsilon()           == 3.14f, boost::test_tools::tolerance(tol));
    }
}