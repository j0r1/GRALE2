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

#ifndef GRALE_MULTIFITNESSHISTORY_H

#define GRALE_MULTIFITNESSHISTORY_H

#include "graleconfig.h"
#include "fitnesshistory.h"

namespace grale
{

class GRALE_IMPORTEXPORT MultiFitnessHistory : public errut::ErrorBase
{
public:
	MultiFitnessHistory(int numFitnessComponents, int historySize, float threshold);
	~MultiFitnessHistory();

	void reset(int historySize, float threshold);
	void reset(float threshold);

	void processValue(int fitnessComponent, float x)					{ m_fitnessHistories[fitnessComponent]->processValue(x); }
	bool advance();
	bool isFinished();
	void printDebugInfo();
public:
	std::vector<FitnessHistory *> m_fitnessHistories;
};

} // end namespace

#endif // GRALE_MULTIFITNESSHISTORY_H

