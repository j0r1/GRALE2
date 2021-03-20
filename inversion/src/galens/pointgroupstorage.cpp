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

#include "pointgroupstorage.h"
#include <grale/imagesdataextended.h>

namespace grale
{

PointGroupStorage::PointGroupStorage()
{
}

PointGroupStorage::~PointGroupStorage()
{
	clear();
}

void PointGroupStorage::init(const std::list<ImagesDataExtended *> &images)
{
	clear();
	
	for (auto it = images.begin() ; it != images.end() ; it++)
		add(*it);
}

void PointGroupStorage::add(const ImagesDataExtended *pImg)
{
	int src = m_groupInfo.size();
	m_groupInfo.resize(src+1);

	if (!pImg)
		return;

	int numgroups = pImg->getNumberOfGroups();
	
	m_groupInfo[src].resize(numgroups);

	for (int group = 0 ; group < numgroups ; group++)
	{
		int np = pImg->getNumberOfGroupPoints(group);

		m_groupInfo[src][group].resize(np*2);

		for (int p = 0 ; p < np ; p++)
		{
			int img, ptidx;
			
			pImg->getGroupPointIndices(group, p, &img, &ptidx);
			m_groupInfo[src][group][p*2] = img;
			m_groupInfo[src][group][p*2+1] = ptidx;
		}
	}
}

void PointGroupStorage::clear()
{
	m_groupInfo.clear();
}

} // end namespace

