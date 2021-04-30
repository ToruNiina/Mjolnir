#ifndef MJOLNIR_UTILITY_COLOR_HPP
#define MJOLNIR_UTILITY_COLOR_HPP
#include <mjolnir/util/io.hpp>

// This utility manipulators output some ANSI escape codes.
// On Windows (not supported by Mjolnir currently), it does not check
// the stream is connected to a tty.

namespace mjolnir
{
namespace io
{
inline std::ostream& bold(std::ostream& os)
{
    if(detail::isatty(os))
    {
        os << "\x1b[01m";
    }
    return os;
}
inline std::ostream& red(std::ostream& os)
{
    if(detail::isatty(os))
    {
        os << "\x1b[31m";
    }
    return os;
}
inline std::ostream& green(std::ostream& os)
{
    if(detail::isatty(os))
    {
        os << "\x1b[32m";
    }
    return os;
}
inline std::ostream& yellow(std::ostream& os)
{
    if(detail::isatty(os))
    {
        os << "\x1b[33m";
    }
    return os;
}
inline std::ostream& blue(std::ostream& os)
{
    if(detail::isatty(os))
    {
        os << "\x1b[34m";
    }
    return os;
}
inline std::ostream& magenta(std::ostream& os)
{
    if(detail::isatty(os))
    {
        os << "\x1b[35m";
    }
    return os;
}
inline std::ostream& cyan(std::ostream& os)
{
    if(detail::isatty(os))
    {
        os << "\x1b[36m";
    }
    return os;
}
inline std::ostream& white(std::ostream& os)
{
    if(detail::isatty(os))
    {
        os << "\x1b[37m";
    }
    return os;
}
inline std::ostream& nocolor(std::ostream& os)
{
    if(detail::isatty(os))
    {
        os << "\x1b[0m";
    }
    return os;
}

} // io
} // mjolnir
#endif// MJOLNIR_UTILITY_COLOR_HPP
