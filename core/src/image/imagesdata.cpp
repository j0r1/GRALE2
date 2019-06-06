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
#include "imagesdata.h"
#include <serut/fileserializer.h>
#include <map>

#include "debugnew.h"

namespace grale
{

ImagesData::ImagesData()
{
}

void ImagesData::copyFrom(const ImagesData &data)
{
	m_images = data.m_images;
	m_intensities = data.m_intensities;
	m_shearComponent1s = data.m_shearComponent1s;
	m_shearComponent2s = data.m_shearComponent2s;
	m_shearWeights = data.m_shearWeights;
	m_groupPoints = data.m_groupPoints;
	m_timeDelayInfo = data.m_timeDelayInfo;
	m_triangulations = data.m_triangulations;
	m_topRight = data.m_topRight;
	m_bottomLeft = data.m_bottomLeft;
}

ImagesData::~ImagesData()
{
	clear();
}

void ImagesData::clear()
{
	m_images.clear();
	m_intensities.clear();
	m_shearComponent1s.clear();
	m_shearComponent2s.clear();
	m_shearWeights.clear();
	m_groupPoints.clear();
	m_triangulations.clear();
	m_timeDelayInfo.clear();
	m_bottomLeft = Vector2D<double>(0,0);
	m_topRight = Vector2D<double>(0,0);
}

ImagesData *ImagesData::createCopy() const
{
	ImagesData *imgdata = new ImagesData();

	imgdata->m_images = m_images;
	imgdata->m_intensities = m_intensities;
	imgdata->m_shearComponent2s = m_shearComponent2s;
	imgdata->m_shearComponent1s = m_shearComponent1s;
	imgdata->m_shearWeights = m_shearWeights;
	imgdata->m_groupPoints = m_groupPoints;
	imgdata->m_triangulations = m_triangulations;
	imgdata->m_topRight = m_topRight;
	imgdata->m_bottomLeft = m_bottomLeft;

	return imgdata;
}

bool ImagesData::create(int numImages, bool intensities, bool shearInfo)
{
	if (numImages <= 0)
	{
		setErrorString("The number of images must be greater than zero");
		return false;
	}
		
	clear();
	
	m_images.resize(numImages);
	if (intensities)
		m_intensities.resize(numImages);
	if (shearInfo)
	{
		m_shearComponent1s.resize(numImages);
		m_shearComponent2s.resize(numImages);
		m_shearWeights.resize(numImages);
	}
	return true;
}

int ImagesData::addImage()
{
	int newNumImages = m_images.size()+1;

	m_images.resize(newNumImages);
	if (m_intensities.size() > 0)
		m_intensities.resize(newNumImages);
	if (m_shearComponent1s.size() > 0)
	{
		m_shearComponent1s.resize(newNumImages);
		m_shearComponent2s.resize(newNumImages);
		m_shearWeights.resize(newNumImages);
	}
	if (m_triangulations.size() > 0)
		m_triangulations.resize(newNumImages);

	return (newNumImages-1); // index of the new image
}
	
int ImagesData::addPoint(int imageNumber, Vector2D<double> point)
{
	if (imageNumber < 0 || imageNumber >= m_images.size())
	{
		setErrorString("Invalid image number");
		return -1;
	}
	
	if (m_intensities.size() > 0)
	{
		setErrorString("Intensity information must be specified for this images data set");
		return -1;
	}
	if (m_shearComponent1s.size() > 0)
	{
		setErrorString("Shear information must be specified for this images data set");
		return -1;
	}
	
	int pointIndex = m_images[imageNumber].size();
	
	m_images[imageNumber].resize(pointIndex+1);
	m_images[imageNumber][pointIndex] = point;
	
	findExtremes();
	
	return pointIndex;
}

int ImagesData::addPoint(int imageNumber, Vector2D<double> point, double intensity)
{
	if (imageNumber < 0 || imageNumber >= m_images.size())
	{
		setErrorString("Invalid image number");
		return false;
	}

	if (m_intensities.size() == 0)
	{
		setErrorString("No intensity information can be specified for this images data set");
		return false;
	}
	if (m_shearComponent1s.size() > 0)
	{
		setErrorString("Shear information must be specified for this images data set");
		return -1;
	}

	int pointIndex = m_images[imageNumber].size();
	
	m_images[imageNumber].resize(pointIndex+1);
	m_images[imageNumber][pointIndex] = point;
	m_intensities[imageNumber].resize(pointIndex+1);
	m_intensities[imageNumber][pointIndex] = intensity;
	
	findExtremes();

	return pointIndex;
}

int ImagesData::addPoint(int imageNumber, Vector2D<double> point, Vector2D<double> shearComponents, double shearWeight)
{
	if (imageNumber < 0 || imageNumber >= m_images.size())
	{
		setErrorString("Invalid image number");
		return false;
	}

	if (m_intensities.size() > 0)
	{
		setErrorString("Intensity information must be specified for this images data set");
		return false;
	}
	if (m_shearComponent1s.size() == 0)
	{
		setErrorString("No shear information can be specified for this images data set");
		return -1;
	}

	int pointIndex = m_images[imageNumber].size();
	
	m_images[imageNumber].resize(pointIndex+1);
	m_images[imageNumber][pointIndex] = point;

	double shearComponent1 = shearComponents.getX();
	double shearComponent2 = shearComponents.getY();
	m_shearComponent1s[imageNumber].resize(pointIndex+1);
	m_shearComponent1s[imageNumber][pointIndex] = shearComponent1;
	m_shearComponent2s[imageNumber].resize(pointIndex+1);
	m_shearComponent2s[imageNumber][pointIndex] = shearComponent2;
	m_shearWeights[imageNumber].resize(pointIndex+1);
	m_shearWeights[imageNumber][pointIndex] = shearWeight;
	
	findExtremes();
	return pointIndex;
}

int ImagesData::addPoint(int imageNumber, Vector2D<double> point, double intensity, Vector2D<double> shearComponents, double shearWeight)
{
	if (imageNumber < 0 || imageNumber >= m_images.size())
	{
		setErrorString("Invalid image number");
		return false;
	}

	if (m_intensities.size() == 0)
	{
		setErrorString("No intensity information can be specified for this images data set");
		return false;
	}
	if (m_shearComponent1s.size() == 0)
	{
		setErrorString("No shear information can be specified for this images data set");
		return -1;
	}

	int pointIndex = m_images[imageNumber].size();
	
	m_images[imageNumber].resize(pointIndex+1);
	m_images[imageNumber][pointIndex] = point;
	m_intensities[imageNumber].resize(pointIndex+1);
	m_intensities[imageNumber][pointIndex] = intensity;

	double shearComponent1 = shearComponents.getX();
	double shearComponent2 = shearComponents.getY();
	m_shearComponent1s[imageNumber].resize(pointIndex+1);
	m_shearComponent1s[imageNumber][pointIndex] = shearComponent1;
	m_shearComponent2s[imageNumber].resize(pointIndex+1);
	m_shearComponent2s[imageNumber][pointIndex] = shearComponent2;
	m_shearWeights[imageNumber].resize(pointIndex+1);
	m_shearWeights[imageNumber][pointIndex] = shearWeight;
	
	findExtremes();
	return pointIndex;
}

int ImagesData::addGroup()
{
	if (m_images.size() == 0)
	{
		setErrorString("No images exist yet");
		return -1;
	}

	int groupIndex = m_groupPoints.size();

	m_groupPoints.resize(groupIndex+1);
	
	return groupIndex;
}

bool ImagesData::addGroupPoint(int groupNumber, int imageIndex, int pointIndex)
{
	if (groupNumber < 0 || groupNumber >= m_groupPoints.size())
	{
		setErrorString("Invalid group number");
		return false;
	}

	if (imageIndex < 0 || imageIndex >= m_images.size())
	{
		setErrorString("Invalid image index");
		return false;
	}

	if (pointIndex < 0 || pointIndex >= m_images[imageIndex].size())
	{
		setErrorString("Invalid point index");
		return false;
	}
	
	int pos = m_groupPoints[groupNumber].size();
	
	m_groupPoints[groupNumber].resize(pos+2);
	m_groupPoints[groupNumber][pos] = imageIndex;
	m_groupPoints[groupNumber][pos+1] = pointIndex;
	
	return true;
}
	
bool ImagesData::addTimeDelayInfo(int imageIndex, int pointIndex, double timeDelay)
{
	if (imageIndex < 0 || imageIndex >= m_images.size())
	{
		setErrorString("Invalid image index");
		return false;
	}

	if (pointIndex < 0 || pointIndex >= m_images[imageIndex].size())
	{
		setErrorString("Invalid point index");
		return false;
	}

	m_timeDelayInfo.push_back(TimeDelayPoint(imageIndex, pointIndex, timeDelay));

	return true;
}

bool ImagesData::addTriangle(int imageIndex, int index1, int index2, int index3)
{
	if (m_triangulations.size() == 0)
		m_triangulations.resize(m_images.size());

	if (imageIndex < 0 || imageIndex >= m_images.size())
	{
		setErrorString("Invalid image index");
		return false;
	}

	if (index1 < 0 || index1 >= m_images[imageIndex].size())
	{
		setErrorString("Invalid point index 1");
		return false;
	}
	if (index2 < 0 || index2 >= m_images[imageIndex].size())
	{
		setErrorString("Invalid point index 2");
		return false;
	}
	if (index3 < 0 || index3 >= m_images[imageIndex].size())
	{
		setErrorString("Invalid point index 3");
		return false;
	}

	if (index1 == index2 || index2 == index3 || index1 == index3)
	{
		setErrorString("Two point indices are equal");
		return false;
	}

	m_triangulations[imageIndex].push_back(TriangleIndices(index1, index2, index3));

	return true;
}

#define IMAGESDATAID 0x41544449

bool ImagesData::load(const std::string &fname)
{
	serut::FileSerializer fser;
	
	if (!fser.open(fname, serut::FileSerializer::ReadOnly))
	{
		setErrorString(std::string("Can't load images data: ") + fser.getErrorString());
		return false;
	}
	if (!read(fser))
		return false;
	return true;
}

bool ImagesData::save(const std::string &fname) const
{
	serut::FileSerializer fser;

	if (!fser.open(fname, serut::FileSerializer::WriteOnly))
	{
		setErrorString(std::string("Can't save images data: ") + fser.getErrorString());
		return false;
	}
	if (!write(fser))
		return false;
	return true;
}

#define IMAGESDATA_FLAG_INTENSITIES	1
#define IMAGESDATA_FLAG_TRIANGULATIONS	2
#define IMAGESDATA_FLAG_TIMEDELAY	4 
#define IMAGESDATA_FLAG_SHEAR		8
#define IMAGESDATA_FLAG_SHEARWEIGHTS	16

bool ImagesData::read(serut::SerializationInterface &si)
{
	int32_t id, numimgs, flag;
	bool gotintens, gottriang, gotTimeDelay, gotShear, gotShearWeights;
	std::vector<std::vector<Vector2D<double> > > imgs;
	std::vector<std::vector<double> > intens, shearComponent1s, shearComponent2s, shearWeights;
	std::vector<TimeDelayPoint> timeDelays;
	std::vector<int32_t> numpoints;
	
	if (!si.readInt32(&id))
	{
		setErrorString(std::string("Error reading images data ID: ") + si.getErrorString());
		return false;
	}
	if (id != IMAGESDATAID)
	{
		setErrorString("Read invalid images data ID");
		return false;
	}
	if (!si.readInt32(&flag))
	{
		setErrorString(std::string("Error reading feature flag: ") + si.getErrorString());
		return false;
	}

	gotintens = ((flag&IMAGESDATA_FLAG_INTENSITIES) == 0)?false:true;
	gottriang = ((flag&IMAGESDATA_FLAG_TRIANGULATIONS) == 0)?false:true;
	gotTimeDelay = ((flag&IMAGESDATA_FLAG_TIMEDELAY) == 0)?false:true;
	gotShear = ((flag&IMAGESDATA_FLAG_SHEAR) == 0)?false:true;
	gotShearWeights = ((flag&IMAGESDATA_FLAG_SHEARWEIGHTS) == 0)?false:true;

	if (gotShearWeights && !gotShear)
	{
		setErrorString("Error in format: expecting shear weights, but not shear");
		return false;
	}
	
	if (!si.readInt32(&numimgs))
	{
		setErrorString(std::string("Error reading number of images: ") + si.getErrorString());
		return false;
	}

	imgs.resize(numimgs);
	numpoints.resize(numimgs);
	if (gotintens)
		intens.resize(numimgs);
	if (gotShear)
	{
		shearComponent1s.resize(numimgs);
		shearComponent2s.resize(numimgs);
		shearWeights.resize(numimgs);
	}
	
	if (!si.readInt32s(numpoints))
	{
		setErrorString(std::string("Error reading images: ") + si.getErrorString());
		return false;
	}

	for (int i = 0 ; i < numimgs ; i++)
	{
		int num = numpoints[i];

		imgs[i].resize(num);
		if (gotintens)
			intens[i].resize(num);
		if (gotShear)
		{
			shearComponent1s[i].resize(num);
			shearComponent2s[i].resize(num);
			shearWeights[i].resize(num);
		}
	}

	for (int i = 0 ; i < numimgs ; i++)
	{
		int num = numpoints[i];
		bool err = false;

		for (int j = 0 ; !err && j < imgs[i].size() ; j++)
		{
			if (!si.readDoubles(imgs[i][j].getComponents(), 2))
				err = true;
		}

		if (!err)
		{
			if (gotintens)
			{
				if (!si.readDoubles(intens[i]))
					err = true;
			}
			if (gotShear && !err)
			{
				if (!si.readDoubles(shearComponent1s[i]))
					err = true;
				else if (!si.readDoubles(shearComponent2s[i]))
					err = true;
				for (auto &w : shearWeights[i])
					w = 1;
			}
			if (gotShearWeights && !err)
			{
				if (!si.readDoubles(shearWeights[i]))
					err = true;
			}
		}

		if (err)
		{
			setErrorString(std::string("Error reading images: ") + si.getErrorString());
			return false;
		}
	}

	int32_t groups;

	if (!si.readInt32(&groups))
	{
		setErrorString(std::string("Error reading number of groups: ") + si.getErrorString());
		return false;
	}

	std::vector<int32_t> numgrppts;
	std::vector<std::vector<int32_t> > grppts;

	if (groups > 0)
	{
		numgrppts.resize(groups);
		grppts.resize(groups);

		bool err = false;

		if (!si.readInt32s(numgrppts))
			err = true;

		for (int i = 0 ; !err && i < groups ; i++)
		{
			int num = numgrppts[i]*2;

			grppts[i].resize(num);
			if (!si.readInt32s(grppts[i]))
				err = true;
		}

		if (err)
		{
			setErrorString(std::string("Error reading group info: ") + si.getErrorString());
			return false;
		}
	}

	std::vector<std::vector<TriangleIndices> > triangulations;
	
	if (gottriang)
	{
		triangulations.resize(imgs.size());

		for (int i = 0 ; i < imgs.size() ; i++)
		{
			int32_t num;

			if (!si.readInt32(&num))
			{
				setErrorString(std::string("Error reading number of triangles: ") + si.getErrorString());
				return false;
			}

			for (int32_t j = 0 ; j < num ; j++)
			{
				int32_t indices[3];

				if (!si.readInt32s(indices, 3))
				{
					setErrorString(std::string("Error reading triangle indices: ") + si.getErrorString());
					return false;
				}

				triangulations[i].push_back(TriangleIndices(indices[0], indices[1], indices[2]));
			}
		}
	}

	if (gotTimeDelay)
	{
		int32_t numDelays;

		if (!si.readInt32(&numDelays))
		{
			setErrorString(std::string("Error reading number of time delay points: ") + si.getErrorString());
			return false;
		}

		for (int32_t i = 0 ; i < numDelays ; i++)
		{
			int32_t indices[2];
			double delay;

			if (!si.readInt32s(indices, 2))
			{
				setErrorString(std::string("Error reading time delay point indices: ") + si.getErrorString());
				return false;
			}
			if (!si.readDouble(&delay))
			{
				setErrorString(std::string("Error reading time delay: ") + si.getErrorString());
				return false;
			}
			timeDelays.push_back(TimeDelayPoint(indices[0], indices[1], delay));
		}
	}

	clear(); // clear old data

	// save new data

	m_images = imgs;
	m_intensities = intens;
	m_groupPoints = grppts;
	m_triangulations = triangulations;
	m_timeDelayInfo = timeDelays;
	m_shearComponent1s = shearComponent1s;
	m_shearComponent2s = shearComponent2s;
	m_shearWeights = shearWeights;

	findExtremes();

	return true;
}

bool ImagesData::write(serut::SerializationInterface &si) const
{
	int flag = 0;

	if (m_intensities.size() != 0)
		flag |= IMAGESDATA_FLAG_INTENSITIES;
	if (m_triangulations.size() != 0)
		flag |= IMAGESDATA_FLAG_TRIANGULATIONS;
	if (m_timeDelayInfo.size() != 0)
		flag |= IMAGESDATA_FLAG_TIMEDELAY;
	if (m_shearComponent1s.size() != 0)
	{
		flag |= IMAGESDATA_FLAG_SHEAR;
		flag |= IMAGESDATA_FLAG_SHEARWEIGHTS;
	}

	if (!si.writeInt32(IMAGESDATAID))
	{
		setErrorString(std::string("Error writing images data ID: ") + si.getErrorString());
		return false;
	}
	if (!si.writeInt32(flag))
	{
		setErrorString(std::string("Error writing feature flag: ") + si.getErrorString());
		return false;
	}
	if (!si.writeInt32(m_images.size()))
	{
		setErrorString(std::string("Error writing number of images: ") + si.getErrorString());
		return false;
	}

	std::vector<int32_t> numImagePoints(m_images.size());

	for (int i = 0 ; i < numImagePoints.size() ; i++)
		numImagePoints[i] = m_images[i].size();

	if (!si.writeInt32s(numImagePoints))
	{
		setErrorString(std::string("Error writing images: ") + si.getErrorString());
		return false;
	}

	for (int i = 0 ; i < m_images.size() ; i++)
	{
		for (int j = 0 ; j < m_images[i].size() ; j++)
		{
			if (!si.writeDoubles(m_images[i][j].getComponents(), 2))
			{
				setErrorString(std::string("Error writing images: ") + si.getErrorString());
				return false;
			}
		}
		if (m_intensities.size() != 0)
		{
			if (!si.writeDoubles(m_intensities[i]))
			{
				setErrorString(std::string("Error writing images: ") + si.getErrorString());
				return false;
			}
		}
		if (m_shearComponent1s.size() != 0)
		{
			if (!si.writeDoubles(m_shearComponent1s[i]))
			{
				setErrorString(std::string("Error writing images: ") + si.getErrorString());
				return false;
			}
			if (!si.writeDoubles(m_shearComponent2s[i]))
			{
				setErrorString(std::string("Error writing images: ") + si.getErrorString());
				return false;
			}
			if (!si.writeDoubles(m_shearWeights[i]))
			{
				setErrorString("Error writing images (shear weights): " + si.getErrorString());
				return false;
			}
		}
	}

	if (!si.writeInt32(m_groupPoints.size()))
	{
		setErrorString(std::string("Error writing number of groups: ") + si.getErrorString());
		return false;
	}

	if (m_groupPoints.size() > 0)
	{
		std::vector<int32_t> numgrouppoints(m_groupPoints.size());

		for (int i = 0 ; i < m_groupPoints.size() ; i++)
			numgrouppoints[i] = m_groupPoints[i].size()/2;

		if (!si.writeInt32s(numgrouppoints))
		{
			setErrorString(std::string("Error writing number of group points: ") + si.getErrorString());
			return false;
		}

		for (int i = 0 ; i < m_groupPoints.size() ; i++)
		{
			if (!si.writeInt32s(m_groupPoints[i]))
			{
				setErrorString(std::string("Error writing group points: ") + si.getErrorString());
				return false;
			}
		}
	}

	if (m_triangulations.size() != 0)
	{
		for (int i = 0 ; i < m_images.size() ; i++) // this is the amount of triangulations we should have
		{
			int32_t indices[3];
			
			// write size

			if (!si.writeInt32(m_triangulations[i].size()))
			{
				setErrorString(std::string("Error writing triangulation size: ") + si.getErrorString());
				return false;
			}

			for (auto it = m_triangulations[i].begin() ; it != m_triangulations[i].end() ; it++)
			{
				indices[0] = (*it).getIndex(0);
				indices[1] = (*it).getIndex(1);
				indices[2] = (*it).getIndex(2);

				if (!si.writeInt32s(indices, 3))
				{
					setErrorString(std::string("Error writing triangulation indices: ") + si.getErrorString());
					return false;
				}
			}
		}
	}

	if (m_timeDelayInfo.size() != 0)
	{
		if (!si.writeInt32(m_timeDelayInfo.size()))
		{
			setErrorString(std::string("Error writing number of time delay points: ") + si.getErrorString());
			return false;
		}

		for (int i = 0 ; i < m_timeDelayInfo.size() ; i++)
		{
			int32_t indices[2];

			indices[0] = m_timeDelayInfo[i].getImageIndex();
			indices[1] = m_timeDelayInfo[i].getPointIndex();

			if (!si.writeInt32s(indices, 2))
			{
				setErrorString(std::string("Error writing time delay point indices: ") + si.getErrorString());
				return false;
			}

			if (!si.writeDouble(m_timeDelayInfo[i].getTimeDelay()))
			{
				setErrorString(std::string("Error writing time delay: ") + si.getErrorString());
				return false;
			}
		}
	}

	return true;
}
	
void ImagesData::findExtremes()
{
	double minx,miny,maxx,maxy;
	bool extset = false;

	minx = 0;
	miny = 0;
	maxx = 0;
	maxy = 0;

	for (int img = 0 ; img < m_images.size() ; img++)
	{
		for (int imgpoint = 0 ; imgpoint < m_images[img].size() ; imgpoint++)
		{
			Vector2D<double> pos = getImagePointPosition(img, imgpoint);

			if (!extset)
			{
				extset = true;
				minx = pos.getX();
				maxx = minx;
				miny = pos.getY();
				maxy = miny;
			}
			else
			{
				if (pos.getX() < minx) minx = pos.getX();
				if (pos.getY() < miny) miny = pos.getY();
				if (pos.getX() > maxx) maxx = pos.getX();
				if (pos.getY() > maxy) maxy = pos.getY();
			}
		}
	}

	m_topRight = Vector2D<double>(maxx,maxy);
	m_bottomLeft = Vector2D<double>(minx,miny);
}

void ImagesData::centerOnPosition(double a, double d)
{
	for (size_t img = 0 ; img < m_images.size() ; img++)
	{
		for (size_t imgpoint = 0 ; imgpoint < m_images[img].size() ; imgpoint++)
		{
			Vector2D<double> pos = m_images[img][imgpoint];
			
			double A = pos.getX();
			double D = pos.getY();
			
			double x2 = COS(d)*COS(D)*COS(A-a)+SIN(D)*SIN(d);
			double y2 = COS(D)*SIN(A-a);
			double z2 = -COS(D)*SIN(d)*COS(A-a)+SIN(D)*COS(d);

			double D2 = ASIN(z2);

			double xproj = x2/COS(D2);
			double yproj = y2/COS(D2);

			double A2 = ACOS(xproj);
			if (yproj < 0)
				A2 = -A2;
			
			m_images[img][imgpoint] = Vector2D<double>(A2,D2);
		}
	}
}

void rotZ(double &x, double &y, double &z, double phi)
{
	double x1 = x*COS(phi) + y*SIN(phi);
	double y1 =-x*SIN(phi) + y*COS(phi);
	x = x1;
	y = y1;
}

void rotY(double &x, double &y, double &z, double phi)
{
	double x1 = x*COS(phi) + z*SIN(phi);
	double z1 =-x*SIN(phi) + z*COS(phi);
	x = x1;
	z = z1;
}

void ImagesData::uncenterOnPosition(double a, double d)
{
	for (size_t img = 0 ; img < m_images.size() ; img++)
	{
		for (size_t imgpoint = 0 ; imgpoint < m_images[img].size() ; imgpoint++)
		{
			Vector2D<double> pos = m_images[img][imgpoint];
			
			double A = pos.getX();
			double D = pos.getY();
			
			double x2 = COS(D)*COS(A);
			double y2 = COS(D)*SIN(A);
			double z2 = SIN(D);

			rotY(x2, y2, z2, -d);
			rotZ(x2, y2, z2, -a);

			double D2 = ASIN(z2);
			double xproj = x2/COS(D2);
			double yproj = y2/COS(D2);

			double A2 = ACOS(xproj);
			if (yproj < 0)
				A2 = -A2;
			
			m_images[img][imgpoint] = Vector2D<double>(A2,D2);
		}
	}
}

void ImagesData::subtractIntensity(double v)
{
	if (!hasIntensities())
		return;
	
	for (size_t img = 0 ; img < m_images.size() ; img++)
		for (size_t imgpoint = 0 ; imgpoint < m_images[img].size() ; imgpoint++)
			m_intensities[img][imgpoint] -= v;
}

void ImagesData::clearTriangulation()
{
	if (m_triangulations.size() == 0)
		return;
	
	m_triangulations.clear();
}

bool ImagesData::getTriangles(int image, std::vector<TriangleIndices> &triangles) const
{
	if (m_triangulations.size() == 0)
	{
		setErrorString("No triangulations have been created");
		return false;
	}

	if (image < 0 || image >= m_triangulations.size())
	{
		setErrorString("Invalid image index was specified");
		return false;
	}

	triangles = m_triangulations[image];
	return true;
}

} // end namespace
