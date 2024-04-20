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

#include "gaparameters.h"

using namespace errut;

namespace grale
{

GAParameters::GAParameters(double selectionPressure, bool useElitism, bool alwaysIncludeBest, double crossOverRate,
		                   double smallMutationSize)
	: EAParameters(EAParameters::GA)
{
	m_selectionPressure = selectionPressure;
	m_useElitism = useElitism;
	m_alwaysIncludeBest = alwaysIncludeBest;
	m_crossOverRate = crossOverRate;
	m_smallMutSize = smallMutationSize;
}

GAParameters::~GAParameters()
{
}

#define GAPARAMETERSFLAG_BIT_USEELITISM			1
#define GAPARAMETERSFLAG_BIT_ALWAYSINCLUDEBEST	2

bool_t GAParameters::readInternal(serut::SerializationInterface &si)
{
	int32_t flags = 0;
	if (!si.readDouble(&m_selectionPressure) || 
	    !si.readInt32(&flags) ||
		!si.readDouble(&m_crossOverRate) ||
		!si.readDouble(&m_smallMutSize))
		return si.getErrorString();

	m_useElitism = (flags&GAPARAMETERSFLAG_BIT_USEELITISM)?true:false;
	m_alwaysIncludeBest = (flags&GAPARAMETERSFLAG_BIT_ALWAYSINCLUDEBEST)?true:false;
	return true;
}

bool_t GAParameters::writeInternal(serut::SerializationInterface &si) const
{
	int32_t flags = 0;

	if (m_useElitism)
		flags |= GAPARAMETERSFLAG_BIT_USEELITISM;
	if (m_alwaysIncludeBest)
		flags |= GAPARAMETERSFLAG_BIT_ALWAYSINCLUDEBEST;

	if (!si.writeDouble(m_selectionPressure) ||
	    !si.writeInt32(flags) ||
		!si.writeDouble(m_crossOverRate) ||
		!si.writeDouble(m_smallMutSize))
		return si.getErrorString();

	return true;
}

} // end namespace

