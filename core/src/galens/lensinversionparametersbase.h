#pragma once

#include "graleconfig.h"
#include <serut/serializationinterface.h>
#include <errut/errorbase.h>

namespace grale
{

class GRALE_IMPORTEXPORT LensInversionParametersBase : public errut::ErrorBase
{
public:
	virtual bool write(serut::SerializationInterface &si) const = 0;
	virtual bool read(serut::SerializationInterface &si) = 0;
};

}