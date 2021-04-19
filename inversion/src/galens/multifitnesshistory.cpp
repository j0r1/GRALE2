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
#include "multifitnesshistory.h"
#include <iostream>
#include <sstream>

#include "debugnew.h"

using namespace std;

namespace grale
{

MultiFitnessHistory::MultiFitnessHistory(int numFitnessComponents, int historySize, float threshold)
{
	for (int i = 0 ; i < numFitnessComponents ; i++)
		m_fitnessHistories.push_back(make_unique<FitnessHistory>(historySize, threshold));
}

MultiFitnessHistory::~MultiFitnessHistory()
{
}

void MultiFitnessHistory::reset(int historySize, float threshold)
{
	for (int i = 0 ; i < m_fitnessHistories.size() ; i++)
		m_fitnessHistories[i]->reset(historySize, threshold);
}

void MultiFitnessHistory::reset(float threshold)
{
	for (int i = 0 ; i < m_fitnessHistories.size() ; i++)
		m_fitnessHistories[i]->reset(threshold);
}

bool MultiFitnessHistory::advance()
{
	for (int i = 0 ; i < m_fitnessHistories.size() ; i++)
	{
		if (!m_fitnessHistories[i]->advance())
		{
			setErrorString(m_fitnessHistories[i]->getErrorString());
			return false;
		}
	}
	return true;
}

bool MultiFitnessHistory::isFinished()
{
	for (int i = 0 ; i < m_fitnessHistories.size() ; i++)
	{
		if (!m_fitnessHistories[i]->isFinished())
			return false;
	}
	return true;
}

std::string MultiFitnessHistory::getDebugInfo()
{
	std::stringstream ss;

	for (int i = 0 ; i < m_fitnessHistories.size() ; i++)
	{
		float curValue, refValue, convergence, thres;

		m_fitnessHistories[i]->getDebugInfo(&curValue, &refValue, &convergence, &thres);
		ss << curValue << "/" << refValue << "(" << convergence << "/" << thres << "/";
		if (m_fitnessHistories[i]->isFinished())
			ss << "true) ";
		else
			ss << "false) ";
	}
	return ss.str();
}

void MultiFitnessHistory::printDebugInfo()
{
	std::cerr << getDebugInfo() << std::endl;
}

} // end namespace
