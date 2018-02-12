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

#ifndef GRALE_FITNESSHISTORY_H

#define GRALE_FITNESSHISTORY_H

#include "graleconfig.h"
#include <errut/errorbase.h>
#include <vector>

namespace grale
{

class GRALE_IMPORTEXPORT FitnessHistory : public errut::ErrorBase
{
public:
	FitnessHistory(int historySize, float threshold)					{ reset(historySize, threshold); }
	~FitnessHistory()									{ }

	void reset(int historySize, float threshold);
	void reset(float threshold)								{ reset(m_fitnessHistory.size(), threshold); }
	void processValue(float x);
	bool advance();
	bool isFinished();
	void getDebugInfo(float *pCurValue, float *pRefValue, float *pConvergence, float *pThresHold);
private:
	std::vector<float> m_fitnessHistory;
	int m_currentPosition;
	bool m_justAdvanced;
	float m_threshold;
};

} // end namespace

#endif // GRALE_FITNESSHISTORY_H
