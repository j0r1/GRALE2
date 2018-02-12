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
#include "galensmodule.h"

#include "debugnew.h"

namespace grale
{

GALensModule::GALensModule()
{
	registerFunctionName("CreateFitnessObject");
}

GALensModule::~GALensModule()
{
}

LensFitnessObject *GALensModule::createFitnessObject() const
{
	if (!isOpen())
	{
		setErrorString("No module has been opened");
		return 0;
	}
	
	return m_createFitness();
}

void GALensModule::onFoundFunctionName(const std::string &name, void *function)
{
	GAModule::onFoundFunctionName(name, function);

	if (name == std::string("CreateFitnessObject"))
		m_createFitness = (CreateFitnessObjectFunction)function;
}
	
} // end namespace


