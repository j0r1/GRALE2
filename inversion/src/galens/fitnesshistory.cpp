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

#include "fitnesshistory.h"
#include "mathfunctions.h"

#include "debugnew.h"

namespace grale
{

void FitnessHistory::reset(int historySize, float threshold)
{
	m_currentPosition = 0;
	m_fitnessHistory.resize(historySize);
	m_justAdvanced = true;
	m_threshold = threshold;
}

void FitnessHistory::processValue(float x)
{
	if (m_justAdvanced)
	{
		m_justAdvanced = false;
		m_fitnessHistory[m_currentPosition%m_fitnessHistory.size()] = x;
	}
	else
	{
		int pos = m_currentPosition%m_fitnessHistory.size();

		if (m_fitnessHistory[pos] > x)
			m_fitnessHistory[pos] = x;
	}
}

bool FitnessHistory::advance()
{
	if (m_justAdvanced)
	{
		setErrorString("No value has been processed yet");
		return false;
	}

	m_currentPosition++;
	m_justAdvanced = true;
	return true;
}

bool FitnessHistory::isFinished()
{
	int pos = m_currentPosition;

	if (m_justAdvanced)
		pos--;

	if (pos < m_fitnessHistory.size())
		return false;

	int realPos = pos%m_fitnessHistory.size();
	int compPos = (pos+1)%m_fitnessHistory.size();

	float x0 = m_fitnessHistory[compPos];
	float x1 = m_fitnessHistory[realPos];

	if (x0 == x1)
		return true;

	if (ABS((x1-x0)/x0) < m_threshold)
		return true;
	return false;
}

void FitnessHistory::getDebugInfo(float *pCurValue, float *pRefValue, float *pConvergence, float *pThresHold)
{
	int pos = m_currentPosition;

	if (m_justAdvanced)
		pos--;

	int realPos = pos%m_fitnessHistory.size();

       	*pCurValue = m_fitnessHistory[realPos];

	if (pos < m_fitnessHistory.size())
	{
		*pRefValue = -1;
		*pConvergence = -1;
	}
	else
	{
		int compPos = (pos+1)%m_fitnessHistory.size();
 		*pRefValue = m_fitnessHistory[compPos];

		if ((*pCurValue) == (*pRefValue))
			*pConvergence = 1;
		else
			*pConvergence = ABS(((*pCurValue)-(*pRefValue))/(*pRefValue));
	}

	*pThresHold = m_threshold;
}

} // end namespace
