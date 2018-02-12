/*

  This file is a part of GRALE, a library to facilitate the simulation
  and inversion of gravitational lenses.

  Copyright (C) 2008-2012 Jori Liesenborgs

  Contact: jori.liesenborgs@gmail.com
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
  
*/

#include "graleconfig.h" 
#include "configurationparameters.h"

#include "debugnew.h"

using namespace std;

namespace grale
{

#define TYPEDPARAMETER_TYPE_NONE 0
#define TYPEDPARAMETER_TYPE_BOOLEAN 1
#define TYPEDPARAMETER_TYPE_INTEGER 2
#define TYPEDPARAMETER_TYPE_REAL 3
#define TYPEDPARAMETER_TYPE_STRING 4

bool TypedParameter::read(serut::SerializationInterface &si)
{
	int32_t t;
	
	if (!si.readInt32(&t))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	clear();

	switch(t)
	{
	case TYPEDPARAMETER_TYPE_NONE:
		break;
	case TYPEDPARAMETER_TYPE_BOOLEAN:
	case TYPEDPARAMETER_TYPE_INTEGER:
		{
			int32_t v;

			if (!si.readInt32(&v))
			{
				setErrorString(si.getErrorString());
				return false;
			}

			if (t == TYPEDPARAMETER_TYPE_BOOLEAN)
			{
				m_type = Boolean;
				m_boolValue = (v != 0);
			}
			else
			{
				m_type = Integer;
				m_intValue = v;
			}
			break;
		}
	case TYPEDPARAMETER_TYPE_REAL:
		{
			m_type = Real;
			if (!si.readDouble(&m_realValue))
			{
				setErrorString(si.getErrorString());
				return false;
			}
			break;
		}
	case TYPEDPARAMETER_TYPE_STRING:
		{
			m_type = String;
			if (!si.readString(m_strValue))
			{
				setErrorString(si.getErrorString());
				return false;
			}
			break;
		}
	default:
		{
			setErrorString("Unknown value type");
			return false;
		}
	}
	return true;
}

bool TypedParameter::write(serut::SerializationInterface &si) const
{
	bool r;

	if (m_type == None)
		r = si.writeInt32(TYPEDPARAMETER_TYPE_NONE);
	else if (m_type == Boolean)
	{
		int32_t x = (m_boolValue)?1:0;
		r = si.writeInt32(TYPEDPARAMETER_TYPE_BOOLEAN) && si.writeInt32(x);
	}
	else if (m_type == Integer)
		r = si.writeInt32(TYPEDPARAMETER_TYPE_INTEGER) && si.writeInt32((int32_t)m_intValue);
	else if (m_type == Real)
		r = si.writeInt32(TYPEDPARAMETER_TYPE_REAL) && si.writeDouble(m_realValue);
	else if (m_type == String)
		r = si.writeInt32(TYPEDPARAMETER_TYPE_STRING) && si.writeString(m_strValue);
	else
	{
		setErrorString("Internal error: unknown TypedParameter::m_type value");
		return false;
	}

	if (!r)
	{
		setErrorString(si.getErrorString());
		return false;
	}

	return true;
}

void TypedParameter::copyFrom(const TypedParameter &src)
{
	m_type = src.m_type;
	m_boolValue = src.m_boolValue;
	m_intValue = src.m_intValue;
	m_realValue = src.m_realValue;
	m_strValue = src.m_strValue;
}

ConfigurationParameters::ConfigurationParameters()
{ 
	m_pParameters = new map<string, ParamWithMarker>();
}

ConfigurationParameters::ConfigurationParameters(const ConfigurationParameters &cfg) 
{
	m_pParameters = new map<string, ParamWithMarker>();
	copyFrom(cfg);
}

ConfigurationParameters::~ConfigurationParameters()							
{ 
	delete m_pParameters;
}

#define CONFIGURATIONPARAMETERSID 0x43464750

bool ConfigurationParameters::read(serut::SerializationInterface &si)
{
	int32_t id;
	
	if (!si.readInt32(&id))
	{
		setErrorString(std::string("Error reading configuration parameters ID: ") + si.getErrorString());
		return false;
	}
	if (id != CONFIGURATIONPARAMETERSID)
	{
		setErrorString("Read invalid configuration parameters ID");
		return false;
	}

	int32_t num;

	if (!si.readInt32(&num))
	{
		setErrorString(std::string("Error reading number of parameters: ") + si.getErrorString());
		return false;
	}

	m_pParameters->clear();
	
	for (int32_t i = 0 ; i < num ; i++)
	{
		ParamWithMarker param;
		string key;

		if (!si.readString(key))
		{
			setErrorString("Error reading parameter key: " + si.getErrorString());
			return false;
		}

		if (!param.read(si))
		{
			setErrorString("Error reading parameter value for key '" + key + "':" + param.getErrorString());
			return false;
		}

		(*m_pParameters)[key] = param;
	}

	return true;
}

bool ConfigurationParameters::write(serut::SerializationInterface &si) const
{
	if (!si.writeInt32(CONFIGURATIONPARAMETERSID))
	{
		setErrorString(std::string("Error writing configuration parameters ID: ") + si.getErrorString());
		return false;
	}

	if (!si.writeInt32(m_pParameters->size()))
	{
		setErrorString(std::string("Error writing number of parameters: ") + si.getErrorString());
		return false;
	}

	for (auto it = m_pParameters->begin() ; it != m_pParameters->end() ; it++)
	{
		if (!si.writeString(it->first))
		{
			setErrorString("Error writing parameter key name: " + si.getErrorString());
			return false;
		}
		if (!it->second.write(si))
		{
			setErrorString("Error writing parameter: " + it->second.getErrorString());
			return false;
		}
	}
	
	return true;
}

const ConfigurationParameters &ConfigurationParameters::operator=(const ConfigurationParameters &cfg)
{
	copyFrom(cfg);
	return *this; 
}

void ConfigurationParameters::copyFrom(const ConfigurationParameters &cfg)
{
	*m_pParameters = *(cfg.m_pParameters);
}

void ConfigurationParameters::clearParameters()
{
	m_pParameters->clear();
}

void ConfigurationParameters::setParameter(const std::string &key, bool v)
{
	(*m_pParameters)[key] = ParamWithMarker(v);
}

void ConfigurationParameters::setParameter(const std::string &key, int v)
{
	(*m_pParameters)[key] = ParamWithMarker(v);
}

void ConfigurationParameters::setParameter(const std::string &key, double v)
{
	(*m_pParameters)[key] = ParamWithMarker(v);
}

void ConfigurationParameters::setParameter(const std::string &key, const std::string &v)
{
	(*m_pParameters)[key] = ParamWithMarker(v);
}

bool ConfigurationParameters::hasParameter(const std::string &key) const
{
	if (m_pParameters->find(key) != m_pParameters->end())
		return true;
	return false;
}

bool ConfigurationParameters::getParameter(const std::string &key, TypedParameter &dst) const
{
	auto it = m_pParameters->find(key);
	if (it == m_pParameters->end())
	{
		setErrorString("Specified key '" + key + "' was not found in the parameters");
		return false;
	}

	it->second.setMarker();
	dst = it->second;
	return true;
}

bool ConfigurationParameters::getParameter(const std::string &key, bool &v) const
{
	TypedParameter param;
	if (!getParameter(key, param))
		return false; // error already set

	if (!param.isBoolean())
	{
		setErrorString("Parameter for specified key '" + key + "' is not a boolean");
		return false;
	}
	v = param.getBooleanValue();
	return true;
}

bool ConfigurationParameters::getParameter(const std::string &key, int &v) const
{
	TypedParameter param;
	if (!getParameter(key, param))
		return false; // error already set

	if (!param.isInteger())
	{
		setErrorString("Parameter for specified key '" + key + "' is not an integer");
		return false;
	}
	v = param.getIntegerValue();
	return true;
}

bool ConfigurationParameters::getParameter(const std::string &key, double &v) const
{
	TypedParameter param;
	if (!getParameter(key, param))
		return false; // error already set

	if (!param.isReal())
	{
		setErrorString("Parameter for specified key '" + key + "' is not a real value");
		return false;
	}
	v = param.getRealValue();
	return true;
}

bool ConfigurationParameters::getParameter(const std::string &key, std::string &v) const
{
	TypedParameter param;
	if (!getParameter(key, param))
		return false; // error already set

	if (!param.isString())
	{
		setErrorString("Parameter for specified key '" + key + "' is not a string");
		return false;
	}
	v = param.getStringValue();
	return true;
}

size_t ConfigurationParameters::getNumberOfParameters() const
{
	return m_pParameters->size();
}

void ConfigurationParameters::clearRetrievalMarkers()
{
	for (auto &p : *m_pParameters)
		p.second.setMarker(false);
}

void ConfigurationParameters::getUnretrievedKeys(std::vector<std::string> &keys)
{
	keys.clear();
	for (auto &p : *m_pParameters)
	{
		if (!p.second.isMarkerSet())
			keys.push_back(p.first);
	}
}

void ConfigurationParameters::getAllKeys(std::vector<std::string> &keys)
{
	keys.clear();
	for (auto &p : *m_pParameters)
		keys.push_back(p.first);
}

} // end namespace

