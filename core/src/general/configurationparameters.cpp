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
#include <assert.h>
#include <iostream>

using namespace std;

namespace grale
{

#define TYPEDPARAMETER_TYPE_NONE 0
#define TYPEDPARAMETER_TYPE_BOOLEAN 1
#define TYPEDPARAMETER_TYPE_INTEGER 2
#define TYPEDPARAMETER_TYPE_REAL 3
#define TYPEDPARAMETER_TYPE_STRING 4

void TypedParameter::dump(const string &prefix) const
{
	cerr << prefix << "m_type    = " << m_type << endl;
	cerr << prefix << "m_isArray = " << ((m_isArray)?"true":"false") << endl;

	auto dumpValues = [&prefix](const string &name, const auto &v, auto convert)
	{
		cerr << prefix << name << " = ";
		cerr << "["; 
		for (auto x: v)
			cerr << convert(x) << ",";
		cerr << "]" << endl;
	};

	dumpValues("m_boolValues", m_boolValues, [](bool v) { return (v)?1:0; });
	dumpValues("m_intValues", m_intValues, [](int v) { return v; });
	dumpValues("m_realValues", m_realValues, [](double v) { return v; });
	dumpValues("m_strValues", m_strValues, [](const string v) { return v; });
}

bool TypedParameter::read(serut::SerializationInterface &si)
{
	int32_t t, num, isArr;
	
	if (!si.readInt32(&t) || !si.readInt32(&num) || !si.readInt32(&isArr))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	if (num < 0)
	{
		setErrorString("Read invalid number of entries (" + to_string(num) + ")");
		return false;
	}

	assert(isArr == 0 || isArr == 1);

	clear();
	m_isArray = (isArr != 0)?true:false;

	if (t == TYPEDPARAMETER_TYPE_NONE)
	{
		if (num != 0)
		{
			setErrorString("Invalid number of entries (" + to_string(num) + ") for type None");
			return false;
		}
	}
	else if (t == TYPEDPARAMETER_TYPE_BOOLEAN || t == TYPEDPARAMETER_TYPE_INTEGER)
	{
		vector<int32_t> values(num);
		if (!si.readInt32s(values))
		{
			setErrorString(si.getErrorString());
			return false;
		}

		auto transfer = [&values](auto &dst, auto convert)
		{
			dst.clear();
			for (auto v : values)
				dst.push_back(convert(v));
		};

		if (t == TYPEDPARAMETER_TYPE_BOOLEAN)
		{
			m_type = Boolean;
			transfer(m_boolValues, [](int32_t v) -> bool { return v != 0; });
		}
		else // integer
		{
			m_type = Integer;
			transfer(m_intValues, [](int32_t v) -> int { return v; });
		}
	}
	else if (t == TYPEDPARAMETER_TYPE_REAL)
	{
		m_type = Real;
		m_realValues.resize(num);
		if (!si.readDoubles(m_realValues))
		{
			setErrorString(si.getErrorString());
			return false;
		}
	}
	else if (t == TYPEDPARAMETER_TYPE_STRING)
	{
		m_type = String;
		for (int i = 0 ; i < num ; i++)
		{
			string s;
			if (!si.readString(s))
			{
				setErrorString(si.getErrorString());
				return false;
			}
			m_strValues.push_back(s);
		}
	}
	else
	{
		setErrorString("Unknown value type");
		return false;		
	}
	return true;
}

bool TypedParameter::write(serut::SerializationInterface &si) const
{
	int32_t isArr = (m_isArray)?1:0;

	if (m_type == None)
	{
		assert(!m_isArray);
		if (!si.writeInt32s({ TYPEDPARAMETER_TYPE_NONE, 0, isArr }))
		{
			setErrorString(si.getErrorString());
			return false;
		}
	}
	else if (m_type == Boolean || m_type == Integer)
	{
		int32_t t;
		vector<int32_t> values;

		auto transfer = [&values](const auto &src, auto convert)
		{
			for (auto x : src)
				values.push_back(convert(x));
		};

		if (m_type == Boolean)
		{
			t = TYPEDPARAMETER_TYPE_BOOLEAN;
			transfer(m_boolValues, [](bool v) -> int32_t { return (v)?1:0; });
		}
		else
		{
			t = TYPEDPARAMETER_TYPE_INTEGER;
			transfer(m_intValues, [](int v) -> int32_t { return v; });
		}

		if (!si.writeInt32s({ t, (int32_t)values.size(), isArr}) || !si.writeInt32s(values))
		{
			setErrorString(si.getErrorString());
			return false;
		}
	}
	else if (m_type == Real)
	{
		if (!si.writeInt32s({TYPEDPARAMETER_TYPE_REAL, (int32_t)m_realValues.size(), isArr}) || 
		    !si.writeDoubles(m_realValues))
		{
			setErrorString(si.getErrorString());
			return false;
		}
	}
	else if (m_type == String)
	{
		if (!si.writeInt32s({TYPEDPARAMETER_TYPE_STRING , (int32_t)m_strValues.size(), isArr}))
		{
			setErrorString(si.getErrorString());
			return false;
		}
		for (const auto &s : m_strValues)
		{
			if (!si.writeString(s))
			{
				setErrorString(si.getErrorString());
				return false;
			}
		}
	}
	else
	{
		setErrorString("Internal error: unknown TypedParameter::m_type value");
		return false;
	}

	return true;
}

void TypedParameter::copyFrom(const TypedParameter &src)
{
	m_isArray = src.m_isArray;
	m_type = src.m_type;
	m_boolValues = src.m_boolValues;
	m_intValues = src.m_intValues;
	m_realValues = src.m_realValues;
	m_strValues = src.m_strValues;
}

ConfigurationParameters::ConfigurationParameters()
{ 
	m_pParameters = make_unique<map<string, ParamWithMarker>>();
}

ConfigurationParameters::ConfigurationParameters(const ConfigurationParameters &cfg) 
{
	m_pParameters = make_unique<map<string, ParamWithMarker>>();
	copyFrom(cfg);
}

ConfigurationParameters::~ConfigurationParameters()							
{ 
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

void ConfigurationParameters::setParameterEmpty(const std::string &key)
{
	(*m_pParameters)[key] = ParamWithMarker();
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

void ConfigurationParameters::setParameter(const std::string &key, const vector<bool> &v)
{
	(*m_pParameters)[key] = ParamWithMarker(v);
}

void ConfigurationParameters::setParameter(const std::string &key, const vector<int> &v)
{
	(*m_pParameters)[key] = ParamWithMarker(v);
}

void ConfigurationParameters::setParameter(const std::string &key, const vector<double> &v)
{
	(*m_pParameters)[key] = ParamWithMarker(v);
}

void ConfigurationParameters::setParameter(const std::string &key, const vector<string> &v)
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

const TypedParameter *ConfigurationParameters::getParameter(const std::string &key) const
{
	auto it = m_pParameters->find(key);
	if (it == m_pParameters->end())
	{
		setErrorString("Specified key '" + key + "' was not found in the parameters");
		return nullptr;
	}
	it->second.setMarker();
	return &(it->second);
}

template<typename X, typename F1, typename F2> 
bool ConfigurationParameters::getParameterTemplate(const string &key, X &v,
												   F1 typeCheck, const string &typeName, F2 getValue) const
{
	const TypedParameter *pParam = getParameter(key);
	if (!pParam)
		return false;
	if (!typeCheck(pParam))
	{
		setErrorString("Parameter for specified key '" + key + "' is not " + typeName);
		return false;
	}
	v = getValue(pParam);
	return true;
}

bool ConfigurationParameters::getParameter(const std::string &key, bool &v) const
{
	return getParameterTemplate(key, v, [](auto p) { return (p->isBoolean() && !p->isArray()); },
	                            "a boolean", [](auto p) { return p->getBooleanValue(); });
}

bool ConfigurationParameters::getParameter(const std::string &key, vector<bool> &v) const
{
	return getParameterTemplate(key, v, [](auto p) { return (p->isBoolean() && p->isArray()); },
	                            "a boolean array", [](auto p) { return p->getBooleanValues(); });
}


bool ConfigurationParameters::getParameter(const std::string &key, int &v) const
{
	return getParameterTemplate(key, v, [](auto p) { return (p->isInteger() && !p->isArray()); },
							"an integer", [](auto p) { return p->getIntegerValue(); });
}

bool ConfigurationParameters::getParameter(const std::string &key, vector<int> &v) const
{
	return getParameterTemplate(key, v, [](auto p) { return (p->isInteger() && p->isArray()); },
							"an integer array", [](auto p) { return p->getIntegerValues(); });
}

bool ConfigurationParameters::getParameter(const std::string &key, double &v) const
{
	return getParameterTemplate(key, v, [](auto p) { return (p->isReal() && !p->isArray()); },
							"a real value", [](auto p) { return p->getRealValue(); });
}

bool ConfigurationParameters::getParameter(const std::string &key, vector<double> &v) const
{
	return getParameterTemplate(key, v, [](auto p) { return (p->isReal() && p->isArray()); },
							"a real value array", [](auto p) { return p->getRealValues(); });
}

bool ConfigurationParameters::getParameter(const std::string &key, std::string &v) const
{
	return getParameterTemplate(key, v, [](auto p) { return (p->isString() && !p->isArray()); },
							"a string", [](auto p) { return p->getStringValue(); });
}

bool ConfigurationParameters::getParameter(const std::string &key, vector<string> &v) const
{
	return getParameterTemplate(key, v, [](auto p) { return (p->isString() && p->isArray()); },
							"a string array", [](auto p) { return p->getStringValues(); });
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

void ConfigurationParameters::getUnretrievedKeys(std::vector<std::string> &keys) const
{
	keys.clear();
	for (auto &p : *m_pParameters)
	{
		if (!p.second.isMarkerSet())
			keys.push_back(p.first);
	}
}

void ConfigurationParameters::getAllKeys(std::vector<std::string> &keys) const
{
	keys.clear();
	for (auto &p : *m_pParameters)
		keys.push_back(p.first);
}

void ConfigurationParameters::dump() const
{
	for (auto &p : *m_pParameters)
	{
		cerr << p.first << " =" << endl;
		p.second.dump("    ");
	}
}

} // end namespace

