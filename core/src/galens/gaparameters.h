/*

  This file is a part of GRALE, a library to facilitate the simulation
  and inversion of gravitational lenses.

  Copyright (C) 2008-2017 Jori Liesenborgs

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

#pragma once

#include "graleconfig.h"
#include "eaparameters.h"

namespace grale
{

class GRALE_IMPORTEXPORT GAParameters : public EAParameters
{
public:
	GAParameters(double selectionPressure = 2.5, bool useElitism = true, 
	             bool alwaysIncludeBest = true, double crossoverRate = 0.9);
	~GAParameters();

	double getSelectionPressure() const													{ return m_selectionPressure; }
	bool getUseElitism() const															{ return m_useElitism; }
	bool getAlwaysIncludeBest() const													{ return m_alwaysIncludeBest; }
	double getCrossOverRate() const														{ return m_crossOverRate; }
protected:
	errut::bool_t readInternal(serut::SerializationInterface &si) override;
	errut::bool_t writeInternal(serut::SerializationInterface &si) const override;
private:
	double m_selectionPressure;
	bool m_useElitism;
	bool m_alwaysIncludeBest;
	double m_crossOverRate;
};

} // end namespace

