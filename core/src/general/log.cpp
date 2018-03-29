#include "log.h"
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <sstream>

#ifdef GRALECONFIG_SUPPORT_LOGGING
#ifdef _WIN32
#include <process.h>

int GetPID()
{
	return (int)_getpid();
}

#else
#include <unistd.h>
#include <libgen.h>
#include <string.h>

int GetPID()
{
	return getpid();
}

std::string GetBaseName(const std::string &s)
{
	char tmp[4096];
	strncpy(tmp, s.c_str(), 4095);
	tmp[4095] = 0;
	return std::string(basename(tmp));
}

#endif // _WIN32
#endif // GRALECONFIG_SUPPORT_LOGGING

using namespace std;

namespace grale
{

#ifdef GRALECONFIG_SUPPORT_LOGGING

string getDateTimeString()
{
	time_t t = time(0);
	struct tm *lt = localtime(&t);
	char str[1024];

	sprintf(str, "[%02d:%02d:%02d %04d/%02d/%02d]", lt->tm_hour, lt->tm_min, lt->tm_sec, 
			                                        1900+lt->tm_year, lt->tm_mon, lt->tm_mday);
	return string(str);
}

string getDateTimeStringForFileName()
{
	time_t t = time(0);
	struct tm *lt = localtime(&t);
	char str[1024];

	sprintf(str, "%04d%02d%02d%02d%02d%02d", 1900+lt->tm_year, lt->tm_mon, lt->tm_mday, lt->tm_hour, lt->tm_min, lt->tm_sec);
	return string(str);
}

Log::Log()
{
	m_pStream = nullptr;
	m_outputLevel = NONE;
}

void Log::init(const string &execName)
{
	char *pTmp = nullptr;
	if ((pTmp = getenv("GRALELOG_LOGNAME")) != nullptr)
	{
		string logName(pTmp);

		if (logName == "stdout")
			openStdOut();
		else if (logName == "stderr")
			openStdErr();
		else
		{
			if (logName == "auto" or logName.length() == 0)
			{
				stringstream ss;

				ss << "/tmp/gralelog_" << GetBaseName(execName) << "_" << getDateTimeStringForFileName() << "_" << GetPID() << ".log";
				logName = ss.str();
			}

			openFile(logName);
		}

		if ((pTmp = getenv("GRALELOG_LEVEL")) != nullptr)
		{
			string level(pTmp);

			if (level == "err" || level == "error" || level == "ERR" || level == "ERROR")
				setLogLevel(ERR);
			else if (level == "wrn" || level == "warn" || level == "WRN" || level == "WARN")
				setLogLevel(WRN);
			else if (level == "inf" || level == "info" || level == "INF" || level == "INFO")
				setLogLevel(INF);
			else if (level == "dbg" || level == "debug" || level == "DBG" || level == "DEBUG")
				setLogLevel(DBG);
			else
			{
				cerr << "Unknown log level specified: '" << level << "'" << endl;
				abort();
			}
		}
		else
		{
			setLogLevel(DBG);
			(*this)(DBG, "No log level specified, setting to DEBUG");
		}

		(*this)(DBG, "Using log: " + logName);
	}
}

Log::~Log()
{
}

void Log::openFile(const string &fileName)
{
	if (m_fileStream.is_open())
		m_fileStream.close();

	m_fileStream.open(fileName, std::fstream::out);
	if (!m_fileStream.is_open())
	{
		cerr << "Unable to open logfile '" << fileName << "'" << endl;
		abort();
	}
	m_pStream = &m_fileStream;
}

void Log::openStdOut()
{
	if (m_fileStream.is_open())
		m_fileStream.close();

	m_pStream = &cout;
}

void Log::openStdErr()
{
	if (m_fileStream.is_open())
		m_fileStream.close();

	m_pStream = &cerr;
}

void Log::operator()(Level lvl, const string &s)
{
	if (lvl > m_outputLevel)
		return;

	if (!m_pStream)
	{
		cerr << "No logging stream was opened" << endl;
		abort();
	}


	string lvlStr = "?????";
	if (lvl == DBG) lvlStr = "DEBUG";
	else if (lvl == INF) lvlStr = "INFO ";
	else if (lvl == WRN) lvlStr = "WARN ";
	else if (lvl == ERR) lvlStr = "ERROR";

	(*m_pStream) << getDateTimeString() << " " << lvlStr << ": " << s << endl;
	m_pStream->flush();
}

#endif // GRALECONFIG_SUPPORT_LOGGING

Log LOG;

} // end namespace

