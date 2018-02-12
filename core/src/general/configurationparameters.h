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

#ifndef GRALE_CONFIGURATIONPARAMETERS_H

#define GRALE_CONFIGURATIONPARAMETERS_H

#include "graleconfig.h"
#include <errut/errorbase.h>
#include <serut/serializationinterface.h>
#include <map>

namespace grale
{

class GRALE_IMPORTEXPORT TypedParameter : public errut::ErrorBase
{
public:
	enum Type { None, Boolean, String, Integer, Real };

	TypedParameter(const TypedParameter &src)							{ copyFrom(src); }
	TypedParameter()													{ clear(); }
	TypedParameter(bool v)												{ clear(); m_type = Boolean; m_boolValue = v; }
	TypedParameter(int v)												{ clear(); m_type = Integer; m_intValue = v; }
	TypedParameter(double v)											{ clear(); m_type = Real; m_realValue = v; }
	TypedParameter(const std::string &v)								{ clear(); m_type = String; m_strValue = v; }
	~TypedParameter()													{ }

	Type getType() const												{ return m_type; }
	bool isEmpty() const												{ return m_type == None; }
	bool isBoolean() const												{ return m_type == Boolean; }
	bool isInteger() const												{ return m_type == Integer; }
	bool isReal() const													{ return m_type == Real; }
	bool isString() const												{ return m_type == String; }

	bool getBooleanValue() const										{ if (!isBoolean()) return false; return m_boolValue; }
	int getIntegerValue() const											{ if (!isInteger()) return 0; return m_intValue; }
	double getRealValue() const											{ if (!isReal()) return 0; return m_realValue; }
	std::string getStringValue() const									{ if (!isString()) return ""; return m_strValue; }

	const TypedParameter &operator=(const TypedParameter &src)			{ copyFrom(src); return *this; }

	bool read(serut::SerializationInterface &si);
	bool write(serut::SerializationInterface &si) const;
private:
	void copyFrom(const TypedParameter &src);
	void clear()														{ m_type = None; m_boolValue = false; m_intValue = 0; m_realValue = 0; m_strValue = ""; }

	Type m_type;
	bool m_boolValue;
	int m_intValue;
	double m_realValue;
	std::string m_strValue;
};

class GRALE_IMPORTEXPORT ConfigurationParameters : public errut::ErrorBase
{
public:
	ConfigurationParameters();
	ConfigurationParameters(const ConfigurationParameters &cfg);
	~ConfigurationParameters();

	bool read(serut::SerializationInterface &si);
	bool write(serut::SerializationInterface &si) const;

	const ConfigurationParameters &operator=(const ConfigurationParameters &cfg);

	size_t getNumberOfParameters() const;
	void getAllKeys(std::vector<std::string> &keys);

	void clearParameters();
	void setParameter(const std::string &key, bool v);
	void setParameter(const std::string &key, int v);
	void setParameter(const std::string &key, double v);
	void setParameter(const std::string &key, const std::string &v);

	bool hasParameter(const std::string &key) const;
	bool getParameter(const std::string &key, TypedParameter &dst) const;
	bool getParameter(const std::string &key, bool &v) const;
	bool getParameter(const std::string &key, int &v) const;
	bool getParameter(const std::string &key, double &v) const;
	bool getParameter(const std::string &key, std::string &v) const;

	void clearRetrievalMarkers();
	void getUnretrievedKeys(std::vector<std::string> &keys);
private:
	void copyFrom(const ConfigurationParameters &cfg);

	class ParamWithMarker : public TypedParameter
	{
	public:
		ParamWithMarker()													{ m_marker = false; }
		ParamWithMarker(const ParamWithMarker &src) : TypedParameter(src)	{ m_marker = src.m_marker; }
		ParamWithMarker(bool v) : TypedParameter(v)							{ m_marker = false; }
		ParamWithMarker(int v) : TypedParameter(v)							{ m_marker = false; }
		ParamWithMarker(double v) : TypedParameter(v)						{ m_marker = false; }
		ParamWithMarker(const std::string &v) : TypedParameter(v)			{ m_marker = false; }

		const ParamWithMarker &operator=(const ParamWithMarker &src)		{ TypedParameter::operator=(src); m_marker = src.m_marker; return *this; }

		bool isMarkerSet() const											{ return m_marker; }
		void setMarker(bool v = true)										{ m_marker = v; }
	private:	
		bool m_marker;
	};

	std::map<std::string, ParamWithMarker> *m_pParameters;
};

} // end namespace

#endif // GRALE_CONFIGURATIONPARAMETERS_H

