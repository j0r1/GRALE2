#ifndef GRALE_LOG_H

#define GRALE_LOG_H

#include "graleconfig.h"
#include <string>
#include <fstream>

namespace grale
{

#ifdef GRALECONFIG_SUPPORT_LOGGING
class Log
{
public:
	enum Level { NONE = -1, ERR, WRN, INF, DBG };

	Log();
	~Log();

	void init(const std::string &execName);
	void operator()(Level lvl, const std::string &s);
private:
	void openFile(const std::string &fileName);
	void openStdOut();
	void openStdErr();
	void setLogLevel(Level lvl)									{ m_outputLevel = lvl; }

	std::ostream *m_pStream;
	std::fstream m_fileStream;
	Level m_outputLevel;
};
#else
class Log
{
public:
	enum Level { NONE = -1, ERR, WRN, INF, DBG };

	Log() { }
	~Log() { }

	void init(const std::string &execName) { }
	void operator()(Level lvl, const std::string &s) { }
};
#endif // GRALECONFIG_SUPPORT_LOGGING

extern Log LOG;

} // end namespace

#endif // GRALE_LOG_H
