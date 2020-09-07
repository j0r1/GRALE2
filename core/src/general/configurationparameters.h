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

	template<typename X> void commonConst(Type t, std::vector<X> &vec, const X &v)
	{
		clear();
		m_type = t;
		vec.resize(1, v);
		m_isArray = false;
	}

	template<typename X> void commonConst(Type t, std::vector<X> &vec, const std::vector<X> &src)
	{
		clear();
		m_type = t;
		vec = src;
		m_isArray = true;
	}

	TypedParameter(const TypedParameter &src)							{ copyFrom(src); }
	TypedParameter()													{ clear(); }
	TypedParameter(bool v) 												{ commonConst(Boolean, m_boolValues, v); }
	TypedParameter(int v)												{ commonConst(Integer, m_intValues, v); }
	TypedParameter(double v)											{ commonConst(Real, m_realValues, v); }
	TypedParameter(const std::string &v)								{ commonConst(String, m_strValues, v); }
	TypedParameter(const char *) = delete;

	TypedParameter(const std::vector<bool> &v) 							{ commonConst(Boolean, m_boolValues, v); }
	TypedParameter(const std::vector<int> &v)							{ commonConst(Integer, m_intValues, v); }
	TypedParameter(const std::vector<double> &v)						{ commonConst(Real, m_realValues, v); }
	TypedParameter(const std::vector<std::string> &v)					{ commonConst(String, m_strValues, v); }

	~TypedParameter()													{ }

	void dump(const std::string &prefix = std::string()) const;

	Type getType() const												{ return m_type; }
	bool isArray() const												{ return m_isArray; }

	bool isEmpty() const												{ return m_type == None; }
	bool isBoolean() const												{ return m_type == Boolean; }
	bool isInteger() const												{ return m_type == Integer; }
	bool isReal() const													{ return m_type == Real; }
	bool isString() const												{ return m_type == String; }
	
	int getNumberOfEntries() const
	{ 
		if (isBoolean())
			return (int)m_boolValues.size();
		if (isInteger())
			return (int)m_intValues.size();
		if (isReal())
			return (int)m_realValues.size();
		if (isString())
			return (int)m_strValues.size();
		return 0;
	}

	template<typename X> X getValue(bool check, int i, X bad, const std::vector<X> &v) const
	{
		if (!check)
			return bad;
		if (i < 0 || i >= (int)v.size())
			return bad;
		return v[i];
	}

	bool getBooleanValue(int i = 0) const								{ return getValue(isBoolean(), i, false, m_boolValues); }
	int getIntegerValue(int i = 0) const								{ return getValue(isInteger(), i, 0, m_intValues); }
	double getRealValue(int i = 0) const								{ return getValue(isReal(), i, 0.0, m_realValues); }
	std::string getStringValue(int i = 0) const							{ return getValue(isString(), i, std::string(""), m_strValues); }

	const std::vector<bool> &getBooleanValues() const					{ return m_boolValues; }
	const std::vector<int> &getIntegerValues() const					{ return m_intValues; }
	const std::vector<double> &getRealValues() const					{ return m_realValues; }
	const std::vector<std::string> &getStringValues() const				{ return m_strValues; }

	const TypedParameter &operator=(const TypedParameter &src)			{ copyFrom(src); return *this; }

	bool read(serut::SerializationInterface &si);
	bool write(serut::SerializationInterface &si) const;
private:
	void copyFrom(const TypedParameter &src);
	void clear()
	{
		m_isArray = false;
		m_type = None;
		m_boolValues.resize(0);
		m_intValues.resize(0);
		m_realValues.resize(0);
		m_strValues.resize(0);
	}

	Type m_type;
	std::vector<bool> m_boolValues;
	std::vector<int> m_intValues;
	std::vector<double> m_realValues;
	std::vector<std::string> m_strValues;
	bool m_isArray;
};

class GRALE_IMPORTEXPORT ConfigurationParameters : public errut::ErrorBase
{
public:
	ConfigurationParameters();
	ConfigurationParameters(const ConfigurationParameters &cfg);
	~ConfigurationParameters();

	void dump() const;

	bool read(serut::SerializationInterface &si);
	bool write(serut::SerializationInterface &si) const;

	const ConfigurationParameters &operator=(const ConfigurationParameters &cfg);

	size_t getNumberOfParameters() const;
	void getAllKeys(std::vector<std::string> &keys) const;

	void clearParameters();

	void setParameterEmpty(const std::string &key);
	void setParameter(const std::string &key, bool v);
	void setParameter(const std::string &key, int v);
	void setParameter(const std::string &key, double v);
	void setParameter(const std::string &key, const std::string &v);
	void setParameter(const std::string &key, const char *) = delete;

	void setParameter(const std::string &key, const std::vector<bool> &v);
	void setParameter(const std::string &key, const std::vector<int> &v);
	void setParameter(const std::string &key, const std::vector<double> &v);
	void setParameter(const std::string &key, const std::vector<std::string> &v);

	bool hasParameter(const std::string &key) const;
	const TypedParameter *getParameter(const std::string &key) const;
	bool getParameter(const std::string &key, TypedParameter &dst) const;

	bool getParameter(const std::string &key, bool &v) const;
	bool getParameter(const std::string &key, int &v) const;
	bool getParameter(const std::string &key, double &v) const;
	bool getParameter(const std::string &key, std::string &v) const;
	bool getParameter(const std::string &key, const char *) const = delete;

	bool getParameter(const std::string &key, std::vector<bool> &v) const;
	bool getParameter(const std::string &key, std::vector<int> &v) const;
	bool getParameter(const std::string &key, std::vector<double> &v) const;
	bool getParameter(const std::string &key, std::vector<std::string> &v) const;

	void clearRetrievalMarkers();
	void getUnretrievedKeys(std::vector<std::string> &keys) const;
private:
	template<typename X, typename F1, typename F2> 
	bool getParameterTemplate(const std::string &key, X &v,
							  F1 typeCheck, const std::string &typeName, F2 getValue) const;

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
		ParamWithMarker(const std::vector<bool> &v) : TypedParameter(v)		{ m_marker = false; }
		ParamWithMarker(const std::vector<int> &v) : TypedParameter(v)		{ m_marker = false; }
		ParamWithMarker(const std::vector<double> &v) : TypedParameter(v)	{ m_marker = false; }
		ParamWithMarker(const std::vector<std::string> &v) : TypedParameter(v) { m_marker = false; }

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

