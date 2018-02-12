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
#include "imagesdataextended.h"
#include "configurationparameters.h"

#include "debugnew.h"

using namespace std;

namespace grale
{

ImagesDataExtended::ImagesDataExtended()
{ 
	m_pExtraParam = new ConfigurationParameters();
	m_Ds = 0;
	m_Dds = 0; 
}

ImagesDataExtended::ImagesDataExtended(const ImagesDataExtended &dat) : ImagesData(dat)
{
	m_pExtraParam = new ConfigurationParameters();
	copyFrom(dat);
}

ImagesDataExtended::ImagesDataExtended(const ImagesData &dat) : ImagesData(dat)
{
	m_pExtraParam = new ConfigurationParameters();
	m_Ds = 0;
	m_Dds = 0; 
}

ImagesDataExtended::ImagesDataExtended(double Ds, double Dds)
{
	m_pExtraParam = new ConfigurationParameters();
	m_Ds = Ds;
	m_Dds = Dds;
}

ImagesDataExtended::~ImagesDataExtended()							
{ 
	delete m_pExtraParam;
}

#define IMAGESDATAEXTENDEDID 0x41544451

bool ImagesDataExtended::read(serut::SerializationInterface &si)
{
	int32_t id;
	
	if (!si.readInt32(&id))
	{
		setErrorString(std::string("Error reading extended images data ID: ") + si.getErrorString());
		return false;
	}
	if (id != IMAGESDATAEXTENDEDID)
	{
		setErrorString("Read invalid extended images data ID");
		return false;
	}

	if (!si.readDouble(&m_Ds))
	{
		setErrorString(std::string("Error reading D_s: ") + si.getErrorString());
		return false;
	}
	
	if (!si.readDouble(&m_Dds))
	{
		setErrorString(std::string("Error reading D_ds: ") + si.getErrorString());
		return false;
	}

	if (!m_pExtraParam->read(si))
	{
		setErrorString("Error reading extra parameters: " + m_pExtraParam->getErrorString());
		return false;
	}

	return ImagesData::read(si);
}

bool ImagesDataExtended::write(serut::SerializationInterface &si) const
{
	if (!si.writeInt32(IMAGESDATAEXTENDEDID))
	{
		setErrorString(std::string("Error writing extended images data ID: ") + si.getErrorString());
		return false;
	}

	if (!si.writeDouble(m_Ds))
	{
		setErrorString(std::string("Error writing D_s: ") + si.getErrorString());
		return false;
	}

	if (!si.writeDouble(m_Dds))
	{
		setErrorString(std::string("Error writing D_ds: ") + si.getErrorString());
		return false;
	}

	if (!m_pExtraParam->write(si))
	{
		setErrorString("Error writing extra parameters: " + m_pExtraParam->getErrorString());
		return false;
	}
	
	return ImagesData::write(si);
}

const ImagesDataExtended &ImagesDataExtended::operator=(const ImagesDataExtended &dat)
{
	copyFrom(dat);
	ImagesData::operator=(dat);
	return *this; 
}

const ImagesDataExtended &ImagesDataExtended::operator=(const ImagesData &dat)
{
	m_Ds = 0;
	m_Dds = 0; 
	m_pExtraParam->clearParameters();
	ImagesData::operator=(dat);
	return *this; 
}

void ImagesDataExtended::copyFrom(const ImagesDataExtended &dat)
{
	m_Ds = dat.m_Ds;
	m_Dds = dat.m_Dds;
	*m_pExtraParam = *(dat.m_pExtraParam);
}

void ImagesDataExtended::clearExtraParameters()
{
	m_pExtraParam->clearParameters();
}

size_t ImagesDataExtended::getNumberOfExtraParameters() const
{
	return m_pExtraParam->getNumberOfParameters();
}

void ImagesDataExtended::getAllExtraParameterKeys(std::vector<std::string> &keys)
{
	m_pExtraParam->getAllKeys(keys);
}

void ImagesDataExtended::setExtraParameter(const std::string &key, bool v)
{
	m_pExtraParam->setParameter(key, v);
}

void ImagesDataExtended::setExtraParameter(const std::string &key, int v)
{
	m_pExtraParam->setParameter(key, v);
}

void ImagesDataExtended::setExtraParameter(const std::string &key, double v)
{
	m_pExtraParam->setParameter(key, v);
}

void ImagesDataExtended::setExtraParameter(const std::string &key, const std::string &v)
{
	m_pExtraParam->setParameter(key, v);
}

bool ImagesDataExtended::hasExtraParameter(const std::string &key) const
{
	return m_pExtraParam->hasParameter(key);
}

bool ImagesDataExtended::getExtraParameter(const std::string &key, TypedParameter &dst) const
{
	if (!m_pExtraParam->getParameter(key, dst))
	{
		setErrorString(m_pExtraParam->getErrorString());
		return false;
	}
	return true;
}

bool ImagesDataExtended::getExtraParameter(const std::string &key, bool &v) const
{
	if (!m_pExtraParam->getParameter(key, v))
	{
		setErrorString(m_pExtraParam->getErrorString());
		return false;
	}
	return true;
}

bool ImagesDataExtended::getExtraParameter(const std::string &key, int &v) const
{
	if (!m_pExtraParam->getParameter(key, v))
	{
		setErrorString(m_pExtraParam->getErrorString());
		return false;
	}
	return true;
}

bool ImagesDataExtended::getExtraParameter(const std::string &key, double &v) const
{
	if (!m_pExtraParam->getParameter(key, v))
	{
		setErrorString(m_pExtraParam->getErrorString());
		return false;
	}
	return true;
}

bool ImagesDataExtended::getExtraParameter(const std::string &key, std::string &v) const
{
	if (!m_pExtraParam->getParameter(key, v))
	{
		setErrorString(m_pExtraParam->getErrorString());
		return false;
	}
	return true;
}

void ImagesDataExtended::clearRetrievalMarkers()
{
	m_pExtraParam->clearRetrievalMarkers();
}

void ImagesDataExtended::getUnretrievedKeys(std::vector<std::string> &keys)
{
	m_pExtraParam->getUnretrievedKeys(keys);
}

} // end namespace

