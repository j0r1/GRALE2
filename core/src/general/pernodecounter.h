#pragma once

#include "graleconfig.h"
#include <errut/errorbase.h>
#include <stdint.h>

namespace grale
{

class GRALE_IMPORTEXPORT PerNodeCounter : public errut::ErrorBase
{
public:
	PerNodeCounter(const std::string &shmFileName);
	~PerNodeCounter();

	int getCount();
private:
	bool openFile();
	bool readWrite(bool increaseOrDecrease);

	std::string m_fileName;
	int m_file;
	uint16_t m_count;
};

} // namespace grale

