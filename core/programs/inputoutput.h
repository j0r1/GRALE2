#ifndef INPUTOUTPUT_H

#define INPUTOUTPUT_H

#define __STDC_FORMAT_MACROS // Need this for PRId64
#include <inttypes.h>
#include <errut/booltype.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <string>
#include <vector>

extern int stdInFileDescriptor;
extern int stdOutFileDescriptor;

errut::bool_t ReadLine(int fd, int msecTimeout, std::string &s);
errut::bool_t ReadLineStdin(int msecTimeout, std::string &s);
errut::bool_t ReadBytes(int fd, std::vector<uint8_t> &bytes);
errut::bool_t ReadBytes(int fd, uint8_t *pBytes, size_t numBytes);
errut::bool_t ReadBytesStdin(std::vector<uint8_t> &bytes);
errut::bool_t ReadBytesStdin(uint8_t *pBytes, size_t numBytes);
errut::bool_t WriteBytes(int fd, const std::vector<uint8_t> &bytes);
errut::bool_t WriteBytes(int fd, const void *pBytes, size_t len);
errut::bool_t WriteBytesStdout(const std::vector<uint8_t> &bytes);
errut::bool_t WriteBytesStdout(const void *pBytes, size_t len);
errut::bool_t WriteLineStdout(const std::string &line);

void SplitLine(const std::string &line, std::vector<std::string> &args, const std::string &separatorChars = " \t",
	       const std::string &quoteChars = "\"'", const std::string &commentStartChars = "#", 
	       bool ignoreZeroLengthFields = true);

bool parseAsInt(const std::string &str, int &number);
bool parseAsInt(const std::string &str, int64_t &number);
bool parseAsDouble(const std::string &str, double &number);
std::string doubleToString(double x);
std::string intToString(int x);
std::string intToString(int64_t x);
std::string stringToString(const std::string &str);
std::string trim(const std::string &str, const std::string &trimChars = " \t\r\n");
bool startsWith(const std::string &str, const std::string &prefix);

std::string strprintf_cstr(const char *format, ...);

#define strprintf(format, ...) strprintf_cstr(stringToString(format).c_str(), __VA_ARGS__ )

//inline std::string strprintf(const std::string &format, ...)
inline std::string strprintf_cstr(const char *format, ...)
{
	const int MAXBUFLEN = 8192;
	char buf[MAXBUFLEN+1];
	va_list ap;

	va_start(ap, format);
#ifndef WIN32
	vsnprintf(buf, MAXBUFLEN, format, ap);
#else
	vsnprintf_s(buf, MAXBUFLEN, _TRUNCATE, format, ap);
#endif // WIN32
	va_end(ap);

	buf[MAXBUFLEN] = 0;
	return std::string(buf);
}

inline std::string doubleToString(double x)
{
	std::string s = strprintf("%.15g", x);
	if (s == "1.#INF")
		return "inf";
	if (s == "-1.#INF")
		return "-inf";
	return s;
}

inline std::string intToString(int x)
{
	return strprintf("%d", x);
}

inline std::string intToString(int64_t x)
{
	return strprintf("%" PRId64, x);
}

inline std::string stringToString(const std::string &str)
{
	return str;
}

#endif // INPUTOUTPUT_H
