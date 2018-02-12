#ifdef WIN32
#include <windows.h>
#endif // WIN32
#include "inputoutput.h"
#ifdef WIN32
#include <io.h>
#include <fcntl.h>

typedef int ssize_t;
#else
#include <unistd.h>
#include <poll.h>
#endif // WIN32
#include <stdio.h>
#include <assert.h>
#include <signal.h>
#include <string>
#include <chrono>
#include <iostream>

using namespace std;
using namespace chrono;
using namespace errut;

// This is nasty but quick workaround
int stdInFileDescriptor = 0;
int stdOutFileDescriptor = 1;

bool_t WriteLineStdout(const string &line)
{
	string l = line + "\n";
	return WriteBytesStdout(l.c_str(), l.length());
}

bool_t WriteBytes(int fd, const void *pBytes, size_t len)
{
	if (len == 0)
		return "Can't write null bytes";

	ssize_t num = write(fd, pBytes, len);
	if (num != (ssize_t)len)
		return strprintf("Incomplete write: %d", num);

	return true;
}

bool_t WriteBytes(int fd, const vector<uint8_t> &bytes)
{
	return WriteBytes(fd, &(bytes[0]), bytes.size());
}

bool_t WriteBytesStdout(const vector<uint8_t> &bytes)
{
	return WriteBytes(stdOutFileDescriptor, bytes);
}

bool_t WriteBytesStdout(const void *pBytes, size_t len)
{
	return WriteBytes(stdOutFileDescriptor, pBytes, len);
}

bool_t ReadBytes(int fd, uint8_t *pBytes, size_t numBytes)
{
	int numLeft = (int)numBytes;
	if (numLeft == 0)
		return "Number of bytes should be at least one";

	uint8_t *pArray = pBytes;
	while (numLeft > 0)
	{
		ssize_t num = read(fd, pArray, numLeft);
		if (num < 0)
			return strprintf("Error in 'read': %d", (int)num);

		if (num == 0)
		{
#ifndef WIN32
			usleep(10000); // avoid a tight loop
#else
			Sleep(10);
#endif // !WIN32
		}
		else
		{
			pArray += num;
			numLeft -= (int)num;
		}
	}
	return true;
}

bool_t ReadBytes(int fd, vector<uint8_t> &bytes)
{
	if (bytes.size() == 0)
		return "Number of bytes should be at least one";
	return ReadBytes(fd, &bytes[0], bytes.size());
}

bool_t ReadBytesStdin(vector<uint8_t> &bytes)
{
	return ReadBytes(stdInFileDescriptor, bytes);
}

bool_t ReadBytesStdin(uint8_t *pBytes, size_t numBytes)
{
	return ReadBytes(stdInFileDescriptor, pBytes, numBytes);
}

bool_t ReadLineStdin(int msecTimeout, string &s)
{
	return ReadLine(stdInFileDescriptor, msecTimeout, s);
}

#ifndef WIN32
bool_t ReadLine(int fd, int msecTimeout, string &s)
{
	string currentLine;
	char inputChar[2] = {0, 0};
	bool done = false;

	// This just reads one character at a time, inefficient but easy
	// Should not be a problem if not much data is read at once

	auto startTime = high_resolution_clock::now();

	while (!done)
	{
		high_resolution_clock::duration diff = high_resolution_clock::now() - startTime;
		if (duration_cast<milliseconds>(diff).count() > msecTimeout)
			return strprintf("Timeout of %d milliseconds is exceeded", msecTimeout);
	
		struct pollfd p;

		p.fd = fd;
		p.events = POLLIN;
		p.revents = 0;

		int status = poll(&p, 1, 10);
		if (status < 0 && errno != EINTR)
		{
			perror("poll:");
			fflush(stderr);
			return strprintf("Error in 'poll', status is %d", status);
		}

		if (status > 0) // something to read
		{
			int r = read(fd, inputChar, 1);
			if (r == 0)
				return "Nothing to read, stdin closed";

			if (r < 0)
			{
				if (errno != EINTR)
				{
					perror("read:");
					fflush(stderr);
					return strprintf("Error in 'read', status is %d", r);
				}
				// Ignore EINTR
			}
			else
			{
				if (inputChar[0] == '\n')
				{
					done = true;
				}
				else
					currentLine += string(inputChar);
			}
		}
	}

	s = currentLine;
	return true;
}
#else
bool_t ReadLine(int fd, int msecTimeout, string &s)
{
	if (fd != 0)
		return "Only stdin input is supported in ReadLine";

	string currentLine;
	HANDLE inputHandle = GetStdHandle(STD_INPUT_HANDLE);
	char inputChar[2] = { 0, 0 };

	bool done = false;
	DWORD avail = 0;
	BOOL success = PeekNamedPipe(inputHandle, NULL, 0, NULL, &avail, NULL);
	if (!success) // assume stdin is redirected, not console input
		return "Only input using a pipe is supported in the windows version";

	auto startTime = high_resolution_clock::now();

	// This just reads one character at a time, inefficient but easy
	// Should not be a problem if not much data is read at once
	while (!done)
	{
		high_resolution_clock::duration diff = high_resolution_clock::now() - startTime;
		if (duration_cast<milliseconds>(diff).count() > msecTimeout)
			return strprintf("Timeout of %d milliseconds is exceeded", msecTimeout);

		bool charRead = false;

		success = PeekNamedPipe(inputHandle, NULL, 0, NULL, &avail, NULL);
		if (!success)
			return "Nothing to read, stdin closed";

		if (avail > 0)
		{
			DWORD numRead = 0;
			BOOL s = ReadFile(inputHandle, (void *)inputChar, 1, &numRead, NULL);
			if (s && numRead == 1)
			{
				// cerr << "success = " << success << " s = " << s << "numRead = " << numRead << endl;
				charRead = true;
			}
		}

		if (!charRead)
			Sleep(100);
		else
		{
			if (inputChar[0] == '\n')
			{
				done = true;
			}
			else
			{
				if (inputChar[0] != '\r') // ignore this on windows
					currentLine.append(inputChar);
			}
		}
	}

	s = currentLine;
	return true;
}
#endif // !WIN32

bool HasCharacter(const string &charList, char c)
{
	for (size_t i = 0 ; i < charList.length() ; i++)
	{
		if (c == charList[i])
			return true;
	}
	return false;
}

void SplitLine(const string &line, vector<string> &args, const string &separatorChars,
	       const string &quoteChars, const string &commentStartChars, bool ignoreZeroLengthFields)
{
	vector<string> arguments;
	size_t startPos = 0;

	while (startPos < line.length() && HasCharacter(separatorChars, line[startPos]))
	{
		startPos++;

		if (!ignoreZeroLengthFields)
			arguments.push_back("");
	}

	string curString("");
	bool done = false;

	if (startPos >= line.length())
	{
		if (!ignoreZeroLengthFields)
			arguments.push_back("");
	}

	while (startPos < line.length() && !done)
	{
		size_t endPos = startPos;
		bool endFound = false;
		bool gotSeparator = false;

		while (!endFound && endPos < line.length())
		{
			if (HasCharacter(separatorChars, line[endPos]) || HasCharacter(commentStartChars, line[endPos]))
			{
				curString += line.substr(startPos, endPos-startPos);
				endFound = true;

				if (HasCharacter(separatorChars, line[endPos]))
				{
					gotSeparator = true;
					endPos++;
				}
			}
			else if (HasCharacter(quoteChars, line[endPos]))
			{
				curString += line.substr(startPos, endPos-startPos);

				char quoteStartChar = line[endPos];

				endPos += 1;
				startPos = endPos;

				while (endPos < line.length() && line[endPos] != quoteStartChar)
					endPos++;

				curString += line.substr(startPos, endPos-startPos);

				if (endPos < line.length())
					endPos++;

				startPos = endPos;
			}
			else
				endPos++;
		}

		if (!endFound)
		{
			if (endPos-startPos > 0)
				curString += line.substr(startPos, endPos-startPos);
		}

		if (curString.length() > 0 || !ignoreZeroLengthFields)
			arguments.push_back(curString);

		if (endPos < line.length() && HasCharacter(commentStartChars, line[endPos]))
			done = true;
		else
		{
			startPos = endPos;
			curString = string("");


			while (startPos < line.length() && HasCharacter(separatorChars, line[startPos]))
			{
				gotSeparator = true;
				startPos++;

				if (!ignoreZeroLengthFields)
					arguments.push_back("");
			}
			
			if (gotSeparator)
			{
				if (startPos >= line.length())
				{
					if (!ignoreZeroLengthFields)
						arguments.push_back("");
				}
			}
		}
	}

	args = arguments;
}

string trim(const string &str, const string &trimChars)
{
	if (str.length() == 0)
		return "";

	bool foundStart = false;
	size_t startIdx = 0;

	while (startIdx < str.length() && !foundStart)
	{
		char c = str[startIdx];

		if (!HasCharacter(trimChars, c))
			foundStart = true;
		else
			startIdx++;
	}

	if (!foundStart || startIdx == str.length()) // trimmed everything
		return "";

	bool foundEnd = false;
	int endIdx = (int)(str.length()) - 1;

	while (endIdx >= 0 && !foundEnd)
	{
		char c = str[endIdx];

		if (!HasCharacter(trimChars, c))
			foundEnd = true;
		else
			endIdx--;
	}

	assert(foundEnd);

	int len = endIdx+1 - startIdx;
	assert(len > 0);

	return str.substr(startIdx, len);
}

bool parseAsInt(const string &str, int &value)
{
	string valueStr = trim(str);
	if (valueStr.length() == 0)
		return false;

	const char *nptr = valueStr.c_str();
	char *endptr;
	
	long int v = strtol(nptr,&endptr,10); // base 10
	
	if (*nptr != '\0')
	{
		if (*endptr != '\0')
		{
			return false;
		}
	}

	value = (int)v;

	if ((long)value != v)
	{
		return false;
	}

	return true;
}

bool parseAsInt(const string &str, int64_t &value)
{
	string valueStr = trim(str);
	if (valueStr.length() == 0)
		return false;

	const char *nptr = valueStr.c_str();
	char *endptr;
	
	long int v = strtol(nptr,&endptr,10); // base 10
	
	if (*nptr != '\0')
	{
		if (*endptr != '\0')
		{
			return false;
		}
	}

	value = v;

	return true;
}

bool parseAsDouble(const string &str, double &value)
{
	string valueStr = trim(str);
	if (valueStr.length() == 0)
		return false;

	if (valueStr == "inf" || valueStr == "+inf")
	{
		value = numeric_limits<double>::infinity();
		return true;
	}
	if (valueStr == "-inf")
	{
		value = -numeric_limits<double>::infinity();
		return true;
	}

	const char *nptr;
	char *endptr;
	
	nptr = valueStr.c_str();
	value = strtod(nptr, &endptr);

	if (*nptr != '\0')
	{
		if (*endptr != '\0')
		{
			return false;
		}
	}
	
	return true;
}

bool startsWith(const std::string &str, const std::string &prefix)
{
	if (prefix.length() > str.length())
		return false;
	if (std::equal(prefix.begin(), prefix.end(), str.begin()))
		return true;
	return false;
}

