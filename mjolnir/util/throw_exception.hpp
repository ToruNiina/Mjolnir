#ifndef MJOLNIR_UTIL_THROW_EXCEPTION_HPP
#define MJOLNIR_UTIL_THROW_EXCEPTION_HPP
#include <sstream>
#include <string>

namespace mjolnir
{

template<typename Exception, typename ...Args>
[[noreturn]] void throw_exception(Args&& ... args)
{
    std::ostringstream oss;
    (oss << ... << args);
    throw Exception(oss.str());
}

} // mjolnir
#endif// MJOLNIR_UTIL_THROW_EXCEPTION_HPP
