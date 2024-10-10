#pragma once

#include <errut/booltype.h>
#include <cstdlib>
#include <stdexcept>
#include <limits>
#include <sstream>
#include <iomanip>

namespace grale 
{

inline errut::bool_t getenv(const std::string &key, std::string &value)
{
	if (!std::getenv(key.c_str()))
		return "No environment variable '" + key + "' is found";
	value = std::string(std::getenv(key.c_str()));
	return true;
}

inline errut::bool_t getenv(const std::string &key, int &value,
							int minValue = (std::numeric_limits<int>::min)(), // wrapping in parentheses because of conflict with windows macros
							int maxValue = (std::numeric_limits<int>::max)())
{
	std::string strValue;
	auto r = getenv(key, strValue);
	if (!r)
		return r;

	auto getErrStr = [key, strValue](const std::string &reason)
	{
		return "Can't convert environment variable " + key + " contents '" + strValue + "' to integer: " + reason;
	};

	try
	{
		value = std::stoi(strValue);
	}
	catch(const std::invalid_argument &e)
	{
		return getErrStr("Doesn't seem an integer value");
	}
	catch(const std::out_of_range &e)
	{
		return getErrStr("Value out of range");
	}
	catch( ... )
	{
		return getErrStr("Unknown error");
	}

	if (value < minValue)
		return getErrStr("Value smaller than minimum value " + std::to_string(minValue));
	if (value > maxValue)
		return getErrStr("Value larger than maximum value " + std::to_string(maxValue));

	return true;
}

inline std::string float_to_string(float f)
{
	std::stringstream ss;
	ss << std::setprecision(std::numeric_limits<float>::max_digits10) << f;
	return ss.str();
}

} // namespace grale
