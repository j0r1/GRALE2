#include "pernodecounter.h"

// TODO: Error version if windows
// TODO: make an implementation using shared memory instead of files

#ifndef WIN32

#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>

using namespace std;

namespace grale
{

PerNodeCounter::PerNodeCounter(const string &shmFileName)
{
	m_file = -1;
	m_count = 0;
	m_fileName = shmFileName;
}

PerNodeCounter::~PerNodeCounter()
{
	if (m_file < 0)
		return;

	bool success = readWrite(false); // This stores in m_count
	close(m_file);
	if (success && m_count == 1) 
	{
		// Means we reduced it to zero, try to unlink the file
		unlink(m_fileName.c_str());
	}
}

int PerNodeCounter::getCount()
{
	if (m_file < 0)
	{
		if (!openFile())
			return -1;

		if (!readWrite(true))
			return -1;
	}
	return m_count;
}

bool PerNodeCounter::openFile()
{
	string prefix("/dev/shm/");
	if (m_fileName.compare(0, prefix.size(), prefix) != 0)
	{
		setErrorString("Filename must start with /dev/shm");
		return false;
	}

	m_file = open(m_fileName.c_str(), O_CREAT|O_RDWR, S_IRUSR|S_IWUSR);
	if (m_file < 0)
	{
		setErrorString("Unable to open file " + m_fileName + ": " + string(strerror(errno)));
		return false;
	}

	return true;
}

bool PerNodeCounter::readWrite(bool increaseOrDecrease)
{
	if (m_file < 0)
	{
		setErrorString("No file open");
		return false;
	}

	struct flock lck;
	memset(&lck, 0, sizeof(struct flock));
	lck.l_whence = SEEK_SET;
	lck.l_start = 0;
	lck.l_len = 2;

	lck.l_type = F_WRLCK;

	// TODO: wanted to use F_OFD_SETLKW, but this isn't in the anaconda headers??
	if (fcntl(m_file, F_SETLKW, &lck) < 0)
	{
		setErrorString("Error creating lock on file");
		return false;
	}

	auto releaseLock = [this,&lck]()
	{
		lck.l_type = F_UNLCK;
		// TODO: wanted to use F_OFD_SETLK, but this isn't in the anaconda headers??
		fcntl(m_file, F_SETLK, &lck);
	};

	lseek(m_file, 0, SEEK_SET);
	auto status = read(m_file, &m_count, 2);
	if (status < 0)
	{
		releaseLock();
		setErrorString("Can't read from file: " + string(strerror(errno)));
		return false;
	}

	if (status < 2) // Not enough data yet
		m_count = 0;

	if (m_count == 0xffff && increaseOrDecrease)
	{
		releaseLock();
		setErrorString("Overflow");
		return false;
	}
	if (m_count == 0 && !increaseOrDecrease)
	{
		releaseLock();
		setErrorString("Underflow");
		return false;
	}

	uint16_t newCount = m_count;
	if (increaseOrDecrease)
		newCount++;
	else
		newCount--;

	lseek(m_file, 0, SEEK_SET);
	if (write(m_file, &newCount, 2) != 2)
	{
		releaseLock();
		setErrorString("Unable to write new value");
		return false;
	}

	releaseLock();
	return true;
}

} // end namespace

#else

namespace grale
{

PerNodeCounter::PerNodeCounter(const string &shmFileName)
{
}

PerNodeCounter::~PerNodeCounter()
{
}

int PerNodeCounter::getCount()
{
	return 0;
}

bool PerNodeCounter::openFile()
{
	setErrorString("Not available on windows");
	return false;
}

bool PerNodeCounter::readWrite(bool increaseOrDecrease)
{
	setErrorString("Not available on windows");
	return false;
}

} // end namespace
#endif // !WIN32

