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

#ifndef GRALE_RANDOMNUMBERGENERATOR_H

#define GRALE_RANDOMNUMBERGENERATOR_H

#include "graleconfig.h"
#include <eatk/randomnumbergenerator.h>
#include <stdint.h>
#include <gsl/gsl_rng.h>
#include <stdlib.h>

namespace grale
{

class GRALE_IMPORTEXPORT RandomNumberGenerator : public eatk::RandomNumberGenerator
{
public:
	RandomNumberGenerator();
	~RandomNumberGenerator();

    double getRandomDouble() override;
    float getRandomFloat() override { return (float)getRandomDouble(); }
	uint32_t getRandomUint32() override;

	uint32_t getSeed() const { return m_seed; }
private:
	static uint32_t pickSeed();

	mutable gsl_rng *m_pRng;
	uint32_t m_seed;
};

} // end namespace

#endif // GRALE_RANDOMNUMBERGENERATOR_H

