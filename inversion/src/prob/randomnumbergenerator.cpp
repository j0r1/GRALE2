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
#include "randomnumbergenerator.h"
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <random>
#include <cassert>

namespace grale
{

RandomNumberGenerator::RandomNumberGenerator()
{
	m_id = std::this_thread::get_id();

	uint32_t seed = pickSeed();
	char *debugSeed;

	if ((debugSeed = getenv("GRALE_DEBUG_SEED")) != 0)
		seed = (uint32_t)strtol(debugSeed, 0, 10);

	//std::cerr << "USING GSL RNG SEED " << seed << std::endl;
	m_pRng = gsl_rng_alloc(gsl_rng_env_setup());
	gsl_rng_set(m_pRng, seed);

	m_seed = seed;
}

RandomNumberGenerator::~RandomNumberGenerator()
{
	gsl_rng_free(m_pRng);
}

double RandomNumberGenerator::getRandomDouble()
{
	assert(std::this_thread::get_id() == m_id);
	return gsl_rng_uniform(m_pRng);
}

uint32_t RandomNumberGenerator::getRandomUint32()
{
	assert(std::this_thread::get_id() == m_id);
	return (uint32_t)gsl_rng_uniform_int(m_pRng, (unsigned long int)0xffffffff);
}

uint32_t RandomNumberGenerator::pickSeed()
{
	std::random_device rd;
	uint32_t seed = rd();

	return seed;
}

} // end namespace

