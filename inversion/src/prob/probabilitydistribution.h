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

/**
 * \file probabilitydistribution.h
 */

#ifndef GRALE_PROBABILITYDISTRIBUTION_H

#define GRALE_PROBABILITYDISTRIBUTION_H

#include "graleconfig.h"
#include <mogal/randomnumbergenerator.h>

namespace grale
{

/** Base class for generating random numbers according to a specific probability distribution. */
class GRALE_IMPORTEXPORT ProbabilityDistribution
{
protected:
	ProbabilityDistribution(const mogal::RandomNumberGenerator *pRndGen);
public:
	virtual ~ProbabilityDistribution();
	virtual double pickNumber() const = 0;
protected:
	double pickRandomNumber() const							{ return m_pRndGen->pickRandomNumber(); }
private:
	const mogal::RandomNumberGenerator *m_pRndGen;
};

} // end namespace

#endif // GRALE_PROBABILITYDISTRIBUTION_H

