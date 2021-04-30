#ifndef MJOLNIR_UTIL_LOGGER_HPP
#define MJOLNIR_UTIL_LOGGER_HPP
#include <mjolnir/util/macro.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <mjolnir/util/color.hpp>
#include <mjolnir/util/range.hpp>
#include <array>
#include <vector>
#include <map>
#include <unordered_map>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>

namespace mjolnir
{

/* these macros are defined.                                                 *
 * - MJOLNIR_LOG_DEBUG  -- for debugging. it makes simulation too slow.      *
 * - MJOLNIR_LOG_INFO   -- output in-detail info unless it slows down.       *
 * - MJOLNIR_LOG_NOTICE -- write current progress and status to console.     *
 * - MJOLNIR_LOG_WARN   -- may not be wrong, but undesirable stuff happens.  *
 * - MJOLNIR_LOG_ERROR  -- something wrong.                                  *
 *                                                                           *
 * - MJOLNIR_LOG_SCOPE    -- helps logging by indentation.                   *
 * - MJOLNIR_LOG_FUNCTION -- helps logging by indentation.                   *
 *
 * - MJOLNIR_GET_DEFAULT_LOGGER -- get a default logger.                     */

namespace logger_detail
{
// since there is no standard way to format containers, the standard library
// does not provide output operators for containers. To make logging easier,
// here output operators are defined for some of the containers.
// To avoid fixing the output format, these are defined in the special namespace
// `logger_detail`. loggers first imports this namespace and output using these
// operators.

template<typename Iterator>
std::ostream& operator<<(std::ostream& os, const range<Iterator>& rg)
{
    os << '[';
    for(const auto& v : rg){os << v << ", ";}
    os << ']';
    return os;
}

template<typename T, std::size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T, N>& ar)
{
    os << '[';
    for(const auto& v : ar){os << v << ", ";}
    os << ']';
    return os;
}

template<typename T, typename Alloc>
std::ostream& operator<<(std::ostream& os, const std::vector<T, Alloc>& vc)
{
    os << '[';
    for(const auto& item : vc){os << item << ", ";}
    os << ']';
    return os;
}

template<typename Key, typename Value, typename Comp, typename Alloc>
std::ostream& operator<<(std::ostream& os, const std::map<Key, Value, Comp, Alloc>& mp)
{
    os << '{';
    for(const auto& kv : mp){os << kv.first << '=' << kv.second << ", ";}
    os << '}';
    return os;
}

template<typename Key, typename Value, typename Hash, typename Pred, typename Alloc>
std::ostream& operator<<(std::ostream& os, const std::unordered_map<Key, Value, Hash, Pred, Alloc>& mp)
{
    os << '{';
    for(const auto& kv : mp){os << kv.first << '=' << kv.second << ", ";}
    os << '}';
    return os;
}

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& os, const std::pair<T1, T2>& pr)
{
    os << '[' << pr.first << ", " << pr.second << ']';
    return os;
}

template<typename ... Ts, std::size_t ... Is>
void tuple_output_helper(std::ostream& os, const std::tuple<Ts...>& t,
                         std::index_sequence<Is...>)
{
    (..., (os << (Is == 0 ? "" : ", ") << std::get<Is>(t)));
    return;
}

template<typename ... Ts>
std::ostream& operator<<(std::ostream& os, const std::tuple<Ts...>& t)
{
    os << '[';
    tuple_output_helper(os, t, std::make_index_sequence<sizeof...(Ts)>{});
    os << ']';
    return os;
}

} // logger_detail

class Logger
{
  public:

    static constexpr inline std::size_t indent_size = 2;

    enum class Level : std::uint8_t
    {
        None   = 0, // internal use only
        Debug  = 1,
        Info   = 2,
        Notice = 3,
        Warn   = 4,
        Error  = 5
    };

  public:

    explicit Logger(const std::string& fname)
        : indentation_needed_(true), indent_(0), fname_(fname)
    {
        // clear the contents and check the file can be opened
        std::ofstream ofs(this->fname_);
        if(!ofs.good())
        {
            std::cerr << "logger: file open error: " << fname_ << std::endl;
            std::exit(1);
        }
    }
    ~Logger() = default;

    void indent()   noexcept {indent_ += 1; return;}
    void unindent() noexcept {indent_ -= 1; return;}

    template<typename ... Ts>
    void log(Level level, Ts&& ... args)
    {
        std::fstream ofs(this->fname_, std::ios_base::out | std::ios_base::app);
        if(!ofs.good())
        {
            throw_exception<std::runtime_error>(
                "Logger: file open error: ", this->fname_);
        }
        if(this->indentation_needed_)
        {
            ofs << std::string(indent_size * indent_, ' ');
        }

        if(level == Level::Notice || level == Level::Warn || level == Level::Error)
        {
            // message is also printed to stderr
            std::cerr << "-- ";
            if(level == Level::Warn)
            {
                std::cerr << '[' << io::bold << io::yellow << "warning" << io::nocolor << "] ";
            }
            else if(level == Level::Error)
            {
                std::cerr << '[' << io::bold << io::red << "error" << io::nocolor << "] ";
            }
            output_message(std::cerr, io::bold, std::forward<Ts>(args)..., io::nocolor);
            std::cerr << std::endl; // output immediately (flush)
        }

        output_message(ofs, this->stringize(level), std::forward<Ts>(args)...);
        ofs << std::endl; // output immediately (flush)

        this->indentation_needed_ = true;
        return;
    }

    template<typename ... Ts>
    void log_no_lf(Level level, Ts&& ... args)
    {
        std::fstream ofs(this->fname_, std::ios_base::out | std::ios_base::app);
        if(!ofs.good())
        {
            throw_exception<std::runtime_error>(
                "Logger: file open error: ", this->fname_);
        }
        if(this->indentation_needed_)
        {
            ofs << std::string(indent_size * indent_, ' ');
        }

        if(level == Level::Notice || level == Level::Warn || level == Level::Error)
        {
            // message is also printed to stderr
            std::cerr << "-- ";
            if(level == Level::Warn)
            {
                std::cerr << '[' << io::bold << io::yellow << "warning" << io::nocolor << "] ";
            }
            else if(level == Level::Error)
            {
                std::cerr << '[' << io::bold << io::red << "error" << io::nocolor << "] ";
            }
            output_message(std::cerr, std::forward<Ts>(args)...);
            std::cerr << std::flush; // no Line Feed
        }

        output_message(ofs, this->stringize(level), std::forward<Ts>(args)...);
        std::cerr << std::flush;

        this->indentation_needed_ = false;
        return;
    }

  private:

    std::string stringize(const Level lv)
    {
        switch(lv)
        {
            case Level::None  : return "";
            case Level::Debug : return "";
            case Level::Info  : return "";
            case Level::Notice: return "";
            case Level::Warn  : return "[warning] ";
            case Level::Error : return "[error] ";
            default: return "UNKNOWN LOG LEVEL: ";
        }
    }

    template<typename ... Args>
    static void output_message(std::ostream& os, Args&& ... args)
    {
        using namespace logger_detail; // import container output operators
        (os << ... << args);
        return;
    }

  private:

    bool indentation_needed_;
    std::size_t indent_;
    std::string fname_;
};

class LoggerManager
{
  public:
    using logger_type    = Logger;
    using resource_type  = std::unique_ptr<logger_type>;
    using container_type = std::map<std::string, resource_type>;

  public:

    static void set_default_logger(const std::string& fname)
    {
        if(default_ == fname)
        {
            std::cerr << "[warning] for Developers: Default Logger("
                      << fname << ") is set twice." << std::endl;
            return;
        }

        default_ = fname;
        if(loggers_.count(fname) != 0)
        {
            std::cerr << "[warning] for Developers: Logger(" << fname << ") is "
                      << "already set. from now, it becomes the default logger."
                      << std::endl;
            return;
        }
        loggers_.emplace(fname, std::make_unique<logger_type>(fname));
        return;
    }

    static logger_type& get_default_logger()
    {
        if(default_.empty())
        {
            throw_exception<std::out_of_range>("mjolnir::LoggerManager: "
                "default logger is not set yet.");
        }
        if(loggers_.count(default_) == 0)
        {
            throw_exception<std::out_of_range>("mjolnir::LoggerManager: "
                "default logger (", default_, ") does not exist");
        }
        return *(loggers_.at(default_));
    }

    static logger_type& get_logger(const std::string& name)
    {
        if(loggers_.count(name) == 0)
        {
            loggers_.emplace(name, std::make_unique<logger_type>(name));
        }
        return *(loggers_.at(name));
    }

  private:

    static inline std::string    default_;
    static inline container_type loggers_;
};

class Scope
{
  public:

    Scope(Logger& trc, const std::string& name, const std::string& loc)
      : start_(std::chrono::system_clock::now()), logger_(trc),
        name_(name), location_(loc)
    {
#if !defined(MJOLNIR_DEBUG)
        // remove template typenames (e.g. [with T = int]) from the name
        // because it often become too long to read
        const auto offset = name_.find('[');
        if(offset != std::string::npos)
        {
            name_.erase(name_.begin() + offset, name_.end());
        }
#endif
        logger_.log(Logger::Level::None, this->name_, " {");
        logger_.log(Logger::Level::None, "--> ", this->location_, ':');
        logger_.indent();
    }
    ~Scope()
    {
        logger_.unindent();
        logger_.log(Logger::Level::None, "} ", this->format_duration(
                    std::chrono::system_clock::now() - this->start_));
    }

    std::string const& name()     const noexcept {return name_;}
    std::string const& location() const noexcept {return location_;}

  private:

    std::string format_duration(const std::chrono::system_clock::duration& dur)
    {
        const auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(dur);
        if(ns.count() < 1000)
        {
            return std::to_string(ns.count()) + " [ns]";
        }
        else if(ns.count() < 1000000) // 10^6, less than milliseconds
        {
            return std::to_string(ns.count() * 1e-3) + " [us]";
        }
        else if(ns.count() < 1000000000) // 10^9, less than seconds
        {
            return std::to_string(ns.count() * 1e-6) + " [ms]";
        }
        else
        {
            return std::to_string(ns.count() * 1e-9) + " [sec]";
        }
    }

  private:
    std::chrono::system_clock::time_point start_;
    Logger& logger_;
    std::string name_;
    std::string location_;
};

// set name of the default log file.
#define MJOLNIR_SET_DEFAULT_LOGGER(name)  LoggerManager::set_default_logger(name)

// get logger
#define MJOLNIR_GET_DEFAULT_LOGGER() auto& l_o_g_g_e_r_ = LoggerManager::get_default_logger()

// normal log
#define MJOLNIR_LOG_INFO(...)   l_o_g_g_e_r_.log(Logger::Level::Info,   __VA_ARGS__)
#define MJOLNIR_LOG_NOTICE(...) l_o_g_g_e_r_.log(Logger::Level::Notice, __VA_ARGS__)
#define MJOLNIR_LOG_WARN(...)   l_o_g_g_e_r_.log(Logger::Level::Warn,   __VA_ARGS__)
#define MJOLNIR_LOG_ERROR(...)  l_o_g_g_e_r_.log(Logger::Level::Error,  __VA_ARGS__)

// no linefeed at the end of line.
#define MJOLNIR_LOG_INFO_NO_LF(...)   l_o_g_g_e_r_.log_no_lf(Logger::Level::Info,   __VA_ARGS__)
#define MJOLNIR_LOG_NOTICE_NO_LF(...) l_o_g_g_e_r_.log_no_lf(Logger::Level::Notice, __VA_ARGS__)
#define MJOLNIR_LOG_WARN_NO_LF(...)   l_o_g_g_e_r_.log_no_lf(Logger::Level::Warn,   __VA_ARGS__)
#define MJOLNIR_LOG_ERROR_NO_LF(...)  l_o_g_g_e_r_.log_no_lf(Logger::Level::Error,  __VA_ARGS__)

// write current scope to log file
#define MJOLNIR_LOG_SCOPE(name) Scope s_c_o_p_e_##__LINE__ (l_o_g_g_e_r_, MJOLNIR_STRINGIZE(name), __FILE__ ":" MJOLNIR_STRINGIZE(__LINE__))
#define MJOLNIR_LOG_FUNCTION()  Scope s_c_o_p_e_##__LINE__ (l_o_g_g_e_r_, MJOLNIR_FUNC_NAME,       __FILE__ ":" MJOLNIR_STRINGIZE(__LINE__))

// loggers that are only enabled when MJOLNIR_DEBUG is defined
#ifdef MJOLNIR_DEBUG
#  define MJOLNIR_GET_DEFAULT_LOGGER_DEBUG() auto& l_o_g_g_e_r_ = LoggerManager::get_default_logger()
#  define MJOLNIR_LOG_DEBUG(...)             l_o_g_g_e_r_.log(Logger::Level::Debug, __VA_ARGS__)
#  define MJOLNIR_LOG_DEBUG_NO_LF(...)       l_o_g_g_e_r_.log_no_lf(Logger::Level::Debug, __VA_ARGS__)
#  define MJOLNIR_LOG_SCOPE_DEBUG(name)      Scope s_c_o_p_e_##__LINE__ (l_o_g_g_e_r_, MJOLNIR_STRINGIZE(name), __FILE__ ":" MJOLNIR_STRINGIZE(__LINE__))
#  define MJOLNIR_LOG_FUNCTION_DEBUG()       Scope s_c_o_p_e_##__LINE__ (l_o_g_g_e_r_, MJOLNIR_FUNC_NAME,       __FILE__ ":" MJOLNIR_STRINGIZE(__LINE__))
#else
#  define MJOLNIR_GET_DEFAULT_LOGGER_DEBUG() /**/
#  define MJOLNIR_LOG_DEBUG(...)             /**/
#  define MJOLNIR_LOG_DEBUG_NO_LF(...)       /**/
#  define MJOLNIR_LOG_SCOPE_DEBUG(name)      /**/
#  define MJOLNIR_LOG_FUNCTION_DEBUG()       /**/
#endif // MJOLNIR_DEBUG

} // mjolnir
#endif /* MJOLNIR_UTIL_LOGGER */
