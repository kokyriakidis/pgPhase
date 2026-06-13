#ifndef PGPHASE_ARG_PARSE_HPP
#define PGPHASE_ARG_PARSE_HPP

// Validating wrappers around std::stoi / std::stod / std::stoll for CLI option
// parsing. The bare std functions throw std::invalid_argument with no context
// when given a non-numeric token (e.g. when a flag that expects a value is
// followed by another flag), which aborts the program. These helpers throw
// std::runtime_error with the offending value and option name so the top-level
// handler can print an actionable message and exit cleanly.

#include <cstddef>
#include <stdexcept>
#include <string>

namespace pgphase_collect {

/// Parse an integer option value, naming @p opt on failure.
inline int parse_int_arg(const char* value, const char* opt) {
    try {
        std::size_t consumed = 0;
        const int result = std::stoi(value, &consumed);
        if (consumed != std::string(value).size())
            throw std::invalid_argument(value);
        return result;
    } catch (const std::exception&) {
        throw std::runtime_error(std::string("invalid integer value \"") + value +
                                 "\" for option " + opt);
    }
}

/// Parse a long-long option value (e.g. chunk sizes), naming @p opt on failure.
inline long long parse_ll_arg(const char* value, const char* opt) {
    try {
        std::size_t consumed = 0;
        const long long result = std::stoll(value, &consumed);
        if (consumed != std::string(value).size())
            throw std::invalid_argument(value);
        return result;
    } catch (const std::exception&) {
        throw std::runtime_error(std::string("invalid integer value \"") + value +
                                 "\" for option " + opt);
    }
}

/// Parse a floating-point option value, naming @p opt on failure.
inline double parse_double_arg(const char* value, const char* opt) {
    try {
        std::size_t consumed = 0;
        const double result = std::stod(value, &consumed);
        if (consumed != std::string(value).size())
            throw std::invalid_argument(value);
        return result;
    } catch (const std::exception&) {
        throw std::runtime_error(std::string("invalid numeric value \"") + value +
                                 "\" for option " + opt);
    }
}

} // namespace pgphase_collect

#endif
