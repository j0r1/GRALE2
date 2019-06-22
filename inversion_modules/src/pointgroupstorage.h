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

#ifndef GRALE_POINTGROUPSTORAGE_H

#define GRALE_POINTGROUPSTORAGE_H

#include <vector>
#include <list>

namespace grale
{

class ImagesDataExtended;

class PointGroupStorage
{
public:
	PointGroupStorage();
	~PointGroupStorage();

	void init(const std::list<ImagesDataExtended *> &images);
	void add(const ImagesDataExtended *pImg);
	void add(const ImagesDataExtended &img)										{ return add(&img); } 
	void clear();

	int getNumberOfSources() const												{ return m_groupInfo.size(); }
	int getNumberOfGroups(int source) const										{ return m_groupInfo[source].size(); }
	int getNumberOfGroupPoints(int source, int group) const						{ return m_groupInfo[source][group].size()/2; }
	void getGroupPointIndices(int source, int group, int idx, int *imageNumber, 
	                          int *pointNumber) const							{ *imageNumber = m_groupInfo[source][group][idx*2]; *pointNumber = m_groupInfo[source][group][idx*2+1]; }
private:
	std::vector<std::vector<std::vector<int> > > m_groupInfo;
};

} // end namespace

#endif // GRALE_POINTGROUPSTORAGE_H

