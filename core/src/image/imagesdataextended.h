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

#ifndef GRALE_IMAGESDATAEXTENDED_H

#define GRALE_IMAGESDATAEXTENDED_H

#include "graleconfig.h"
#include "imagesdata.h"
#include <memory>

namespace grale
{

class TypedParameter;
class ConfigurationParameters;

class GRALE_IMPORTEXPORT ImagesDataExtended : public ImagesData
{
public:
	ImagesDataExtended();
	ImagesDataExtended(const ImagesDataExtended &dat);
	ImagesDataExtended(const ImagesData &dat);
	ImagesDataExtended(double Ds, double Dds);
	~ImagesDataExtended();

	double getDs()	const												{ return m_Ds; }
	double getDds() const												{ return m_Dds; }

	void setDs(double v)												{ m_Ds = v; }
	void setDds(double v)												{ m_Dds = v; }

	bool read(serut::SerializationInterface &si);
	bool write(serut::SerializationInterface &si) const;

	const ImagesDataExtended &operator=(const ImagesDataExtended &dat);
	const ImagesDataExtended &operator=(const ImagesData &dat);

	size_t getNumberOfExtraParameters() const;
	void getAllExtraParameterKeys(std::vector<std::string> &keys);

	void clearExtraParameters();
	void setExtraParameter(const std::string &key, bool v);
	void setExtraParameter(const std::string &key, int v);
	void setExtraParameter(const std::string &key, double v);
	void setExtraParameter(const std::string &key, const std::string &v);

	bool hasExtraParameter(const std::string &key) const;
	bool getExtraParameter(const std::string &key, TypedParameter &dst) const;
	bool getExtraParameter(const std::string &key, bool &v) const;
	bool getExtraParameter(const std::string &key, int &v) const;
	bool getExtraParameter(const std::string &key, double &v) const;
	bool getExtraParameter(const std::string &key, std::string &v) const;

	void clearRetrievalMarkers();
	void getUnretrievedKeys(std::vector<std::string> &keys);
private:
	void copyFrom(const ImagesDataExtended &dat);

	double m_Ds, m_Dds;
	std::unique_ptr<ConfigurationParameters> m_pExtraParam;
};

} // end namespace

#endif // GRALE_IMAGESDATAEXTENDED_H

