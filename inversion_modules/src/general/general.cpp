/*

  This file is a part of the GRALE Lens inversion modules, which determine
  the fitness of a lens model for the genetic algorithm used by the
  GRALE library.

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

#include <grale/lensinversiongafactorysingleplanecpu.h>
#include "general_template.h"

extern "C"
{
	GAMODS_EXPORT mogal::GAFactory *CreateFactoryInstance()
	{
		return new grale::LensInversionGAFactory_General<grale::LensInversionGAFactorySinglePlaneCPU>();
	}

	GAMODS_EXPORT grale::LensFitnessObject *CreateFitnessObject()
	{
		return new grale::LensFitnessGeneral();
	}
}

