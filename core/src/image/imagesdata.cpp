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
#include <iostream>
#include <string>
#include <sstream>

#include "debugnew.h"

using namespace std;

namespace grale
{

ImagesData::ImagesData()
{
}

void ImagesData::copyFrom(const ImagesData &data)
{
	m_properties = data.m_properties;
	m_images = data.m_images;
	m_imagePointProperties = data.m_imagePointProperties;
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
	m_properties.clear();
	m_images.clear();
	m_imagePointProperties.clear();
	m_groupPoints.clear();
	m_triangulations.clear();
	m_timeDelayInfo.clear();
	m_bottomLeft = Vector2D<double>(0,0);
	m_topRight = Vector2D<double>(0,0);
}

ImagesData *ImagesData::createCopy() const
{
	ImagesData *imgdata = new ImagesData();

	imgdata->m_properties = m_properties;
	imgdata->m_images = m_images;
	imgdata->m_imagePointProperties = m_imagePointProperties;
	imgdata->m_groupPoints = m_groupPoints;
	imgdata->m_triangulations = m_triangulations;
	imgdata->m_topRight = m_topRight;
	imgdata->m_bottomLeft = m_bottomLeft;

	return imgdata;
}

bool ImagesData::create(int numImages, const vector<PropertyName> &properties)
{
	if (numImages < 0)
	{
		setErrorString("The number of images must be greater than zero");
		return false;
	}

	for (auto n : properties)
	{
		if (n < 0 || n >= MaxProperty)
		{
			setErrorString("Invalid property identifier " + to_string((int)n));
			return false;
		}
	}

	clear();
	
	m_images.resize(numImages);

	m_properties.resize(MaxProperty, false);
	m_imagePointProperties.resize(MaxProperty);
	for (auto n : properties)
	{
		m_properties[n] = true;
		m_imagePointProperties[n].resize(numImages);
	}
	return true;
}

bool ImagesData::create(int numImages, bool intensities, bool shearInfo)
{
	vector<PropertyName> props;
	if (intensities)
		props.push_back(Intensity);
	if (shearInfo)
		props.insert(props.end(), { ShearComponent1, ShearComponent2, Weight });
	return create(numImages, props);
}

int ImagesData::addImage()
{
	if (m_properties.size() != MaxProperty)
	{
		setErrorString("Must create the instance first");
		return false;
	}
	
	int newNumImages = m_images.size()+1;

	m_images.resize(newNumImages);
	for (int i = 0 ; i < m_properties.size() ; i++)
	{
		if (m_properties[i])
			m_imagePointProperties[i].resize(newNumImages);
	}

	if (m_triangulations.size() > 0)
		m_triangulations.resize(newNumImages);

	return (newNumImages-1); // index of the new image
}

int ImagesData::addPoint(int imageNumber, Vector2Dd point, const vector<pair<PropertyName, double>> &properties)
{
	auto error = [this](const string &s) -> int
	{
		setErrorString(s);
		return -1;
	};

	if (imageNumber < 0 || imageNumber >= m_images.size())
		return error("Invalid image number");

	bool specifiedProps[MaxProperty] = { false }; // False should be 0, so that should be fine

	for (auto pv : properties)
	{
		const PropertyName n = pv.first;
		const double value = pv.second;

		if (n < 0 || n >= MaxProperty)
			return error("Invalid property " + to_string((int)n));

		if (specifiedProps[n])
			return error("Property " + to_string((int)n) + " is specified more than once");

		specifiedProps[n] = true;
	}

	auto listProps = [](auto &v) -> string
	{
		stringstream ss;
		ss << "[ ";
		for (int i = 0 ; i < MaxProperty ; i++)
		{
			bool x = v[i];
			ss << ((x)?"true ":"false ");
		}
		ss << "]";
		return ss.str();
	};

	assert(MaxProperty == m_properties.size());
	for (int i = 0 ; i < MaxProperty ; i++)
	{
		if (specifiedProps[i] != m_properties[i])
			return error("Not all required properties are specified, differs from creation time (" + listProps(specifiedProps) + " != " + listProps(m_properties) + ")");
	}
	
	int pointIndex = m_images[imageNumber].size();
	m_images[imageNumber].push_back(point);
	
	for (auto pv : properties)
	{
		const PropertyName n = pv.first;
		const double value = pv.second;
	
		assert(m_imagePointProperties[n].size() == m_images.size());
		assert(m_imagePointProperties[n][imageNumber].size() == pointIndex);
		m_imagePointProperties[n][imageNumber].push_back(value);
	}
		
	findExtremes();
	return pointIndex;
}

int ImagesData::addPoint(int imageNumber, Vector2D<double> point)
{
	return addPoint(imageNumber, point, {});
}

int ImagesData::addPoint(int imageNumber, Vector2D<double> point, double intensity)
{
	return addPoint(imageNumber, point, { { Intensity, intensity }});
}

int ImagesData::addPoint(int imageNumber, Vector2D<double> point, Vector2D<double> shearComponents, double shearWeight)
{
	return addPoint(imageNumber, point, { 
		{ ShearComponent1, shearComponents.getX() },
		{ ShearComponent2, shearComponents.getY() },
		{ Weight, shearWeight }
	});
}

int ImagesData::addPoint(int imageNumber, Vector2D<double> point, double intensity, Vector2D<double> shearComponents, double shearWeight)
{
	return addPoint(imageNumber, point, {
		{ Intensity, intensity },
		{ ShearComponent1, shearComponents.getX() },
		{ ShearComponent2, shearComponents.getY() },
		{ Weight, shearWeight }
	});
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

#define IMAGESDATAID_OLD 0x41544449
#define IMAGESDATAID 0x4154444a

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

bool ImagesData::read(serut::SerializationInterface &si)
{
	int32_t id;
	if (!si.readInt32(&id))
	{
		setErrorString(std::string("Error reading images data ID: ") + si.getErrorString());
		return false;
	}

	if (id == IMAGESDATAID_OLD)
		return readOld(si);

	if (id != IMAGESDATAID)
	{
		setErrorString("Read invalid images data ID");
		return false;
	}

	// Clear the data so that we can fill it in directly
	clear();

	auto error = [this](const string &errMsg) -> bool
	{
		setErrorString(errMsg);
		clear(); // Clear it again so that we don't get an inconsistent state
		return false;
	};

	int32_t flag;
	if (!si.readInt32(&flag))
		return error("Can't read flag: " + si.getErrorString());
	if (flag != 0)
		return error("Can't handle flag " + to_string(flag) + ", currently only 0 is known");

	int32_t numProps;
	if (!si.readInt32(&numProps))
		return error("Can't read number of properties: " + si.getErrorString());
	if (numProps < 0)
		return error("Don't understand number of properties " + to_string(numProps));
	
	m_properties.resize(numProps, false);

	// TODO: check if a property is _used_ that's beyond the known properties

	int32_t numImages;
	if (!si.readInt32(&numImages))
		return error("Can't read number of images: " + si.getErrorString());
	if (numImages < 0)
		return error("Invalid number of images " + to_string(numImages));
	m_images.resize(numImages);

	for (auto &imgPts : m_images)
	{
		int32_t numPts;
		if (!si.readInt32(&numPts))
			return error("Can't read number of image points: " + si.getErrorString());
		if (numPts < 0)
			return error("Invalid number of image points " + to_string(numPts));
		
		imgPts.resize(numPts);
		for (auto &pt : imgPts)
		{
			if (!si.readDoubles(pt.getComponents(), 2))
				return error("Can't read image point coordinates: " + si.getErrorString());
		}
	}

	m_imagePointProperties.resize(numProps);
	for (int i = 0 ; i < numProps ; i++)
	{
		auto &prop = m_imagePointProperties[i];
		int32_t numPropImgs;

		if (!si.readInt32(&numPropImgs))
			return error("Can't read number of images for property: " + si.getErrorString());
		
		if (numPropImgs == -1) // m_properties[i] already set to false;
			continue;

		m_properties[i] = true;
		if (numPropImgs != numImages)
			return error("Inconsistent number of images, expecting " + to_string(numImages) + " but got " + to_string(numPropImgs));
		
		prop.resize(numImages);
		for (size_t j = 0 ; j < prop.size() ; j++)
		{
			int32_t numPropImgPts;
			if (!si.readInt32(&numPropImgPts))
				return error("Can't read number of image points according to property: " + si.getErrorString());

			const int32_t expectedPoints = (int32_t)m_images[j].size();
			if (numPropImgPts != expectedPoints)
				return error("Inconsistent number of image points, expecting " + to_string(expectedPoints) + " but got " + to_string(numPropImgPts));

			auto &propPts = prop[j];
			propPts.resize(numPropImgPts);
			if (!si.readDoubles(propPts))
				return error("Error reading property values: " + si.getErrorString());
		}
	}

	int32_t numGroups;
	if (!si.readInt32(&numGroups))
		return error("Error reading number of point groups: " + si.getErrorString());
	if (numGroups < 0)
		return error("Invalid number of groups " + to_string(numGroups));
	
	m_groupPoints.resize(numGroups);
	for (auto &grp : m_groupPoints)
	{
		int32_t numGrpPts;
		if (!si.readInt32(&numGrpPts))
			return error("Can't read number of group points: " + si.getErrorString());
		if (numGrpPts < 0 || numGrpPts%2 != 0)
			return error("Invalid number of points, should be positive and even ((imgidx, pointidx) pairs), but is " + to_string(numGrpPts));
		grp.resize(numGrpPts);

		if (!si.readInt32s(grp))
			return error("Can't read point group info: " + si.getErrorString());
		
		// Check that the image,point pairs exist
		for (size_t i = 0 ; i < grp.size() ; i += 2)
		{
			int imgIdx = grp[i];
			int ptIdx = grp[i+1];

			if (imgIdx < 0 || imgIdx >= m_images.size())
				return error("Invalid image index " + to_string(imgIdx));
			if (ptIdx < 0 || ptIdx >= m_images[imgIdx].size())
				return error("Invalid point index " + to_string(ptIdx) + " in image idx " + to_string(imgIdx));
		}
	}

	int32_t numTriangulations;
	if (!si.readInt32(&numTriangulations))
		return error("Can't read number of triangulations: " + si.getErrorString());
	if (!(numTriangulations == 0 || numTriangulations == m_images.size()))
		return error("Invalid number of triangulations " + to_string(numTriangulations));

	if (numTriangulations != 0)
	{
		m_triangulations.resize(numTriangulations);
		for (int imgIdx = 0 ; imgIdx < numTriangulations ; imgIdx++)
		{
			auto &triang = m_triangulations[imgIdx];
			int32_t numTriangles;
			if (!si.readInt32(&numTriangles))
				return error("Can't read number of triangles: " + si.getErrorString());
			if (numTriangles < 0)
				return error("Invalid number of triangles " + to_string(numTriangles));

			for (int32_t i = 0 ; i < numTriangles ; i++)
			{
				int32_t indices[3];
				if (!si.readInt32s(indices, 3))
					return error("Can't read triangle indices");

				for (int j = 0 ; j < 3 ; j++)
				{
					int ptIdx = indices[j];
					if (ptIdx < 0 || ptIdx >= m_images[imgIdx].size())
						return error("Invalid point index " + to_string(ptIdx) + " in triangulation of image idx " + to_string(imgIdx));
				}

				triang.push_back(TriangleIndices(indices[0], indices[1], indices[2]));
			}
		}
	}

	int32_t numTds;
	if (!si.readInt32(&numTds))
		return error("Can't read number of time delays: " + si.getErrorString());
	if (numTds < 0)
		return error("Read invalid number of time delays " + to_string(numTds));
	
	for (int32_t i = 0 ; i < numTds ; i++)
	{
		int32_t indices[2];
		if (!si.readInt32s(indices, 2))
			return error("Can't read point indices for time delay: " + si.getErrorString());
		
		const int32_t imgIdx = indices[0];
		if (imgIdx < 0 || imgIdx >= m_images.size())
			return error("Read invalid image index for time delay " + to_string(imgIdx));

		const int32_t ptIdx = indices[1];
		if (ptIdx < 0 || ptIdx >= m_images[imgIdx].size())
			return error("Read invalid point index for time delay " + to_string(ptIdx) + " in image idx " + to_string(imgIdx));

		double delay;
		if (!si.readDouble(&delay))
			return error("Error reading time delay: " + si.getErrorString());
		
		m_timeDelayInfo.push_back(TimeDelayPoint(imgIdx, ptIdx, delay));
	}

	findExtremes();
	return true;
}

#define IMAGESDATAOLD_FLAG_INTENSITIES	1
#define IMAGESDATAOLD_FLAG_TRIANGULATIONS	2
#define IMAGESDATAOLD_FLAG_TIMEDELAY	4 
#define IMAGESDATAOLD_FLAG_SHEAR		8
#define IMAGESDATAOLD_FLAG_SHEARWEIGHTS	16

bool ImagesData::readOld(serut::SerializationInterface &si)
{
	int32_t numimgs, flag;
	bool gotintens, gottriang, gotTimeDelay, gotShear, gotShearWeights;
	std::vector<std::vector<Vector2D<double> > > imgs;
	std::vector<std::vector<double> > intens, shearComponent1s, shearComponent2s, shearWeights;
	std::vector<TimeDelayPoint> timeDelays;
	std::vector<int32_t> numpoints;
	
	if (!si.readInt32(&flag))
	{
		setErrorString(std::string("Error reading feature flag: ") + si.getErrorString());
		return false;
	}

	gotintens = ((flag&IMAGESDATAOLD_FLAG_INTENSITIES) == 0)?false:true;
	gottriang = ((flag&IMAGESDATAOLD_FLAG_TRIANGULATIONS) == 0)?false:true;
	gotTimeDelay = ((flag&IMAGESDATAOLD_FLAG_TIMEDELAY) == 0)?false:true;
	gotShear = ((flag&IMAGESDATAOLD_FLAG_SHEAR) == 0)?false:true;
	gotShearWeights = ((flag&IMAGESDATAOLD_FLAG_SHEARWEIGHTS) == 0)?false:true;

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

	m_properties.resize(MaxProperty, false);
	m_images = imgs;
	if (gotintens)
	{
		m_properties[Intensity] = true;
		m_imagePointProperties[Intensity] = intens;
	}

	m_groupPoints = grppts;
	m_triangulations = triangulations;
	m_timeDelayInfo = timeDelays;

	if (gotShear)
	{
		m_properties[ShearComponent1] = true;
		m_properties[ShearComponent2] = true;
		m_properties[Weight] = true;

		m_imagePointProperties[ShearComponent1] = shearComponent1s;
		m_imagePointProperties[ShearComponent2] = shearComponent2s;
		m_imagePointProperties[Weight] = shearWeights;
	}

	findExtremes();
	return true;
}

bool ImagesData::write(serut::SerializationInterface &si) const
{
	if (!si.writeInt32(IMAGESDATAID))
	{
		setErrorString(std::string("Error writing images data ID: ") + si.getErrorString());
		return false;
	}

	int32_t flag = 0; // For now just zero, perhaps this is useful later
	if (!si.writeInt32(flag))
	{
		setErrorString(std::string("Error writing feature flag: ") + si.getErrorString());
		return false;
	}

	const int32_t numProps = m_properties.size();
	if (!si.writeInt32(numProps))
	{
		setErrorString("Can't write image property flags: " + si.getErrorString());
		return false;
	}

	if (!si.writeInt32(m_images.size()))
	{
		setErrorString("Error writing number of images: " + si.getErrorString());
		return false;
	}

	for (auto &imgPts : m_images)
	{
		const int32_t numPts = imgPts.size();
		if (!si.writeInt32(numPts))
		{
			setErrorString("Can't write number of image points: " + si.getErrorString());
			return false;
		}

		for (auto &pt : imgPts)
		{
			if (!si.writeDoubles(pt.getComponents(), 2))
			{
				setErrorString("Error writing image point: " + si.getErrorString());
				return false;
			}
		}
	}

	if (m_imagePointProperties.size() != numProps)
	{
		setErrorString("Internal error: mismatch in property count");
		return false;
	}

	for (int i = 0 ; i < numProps ; i++)
	{
		auto &prop = m_imagePointProperties[i];
		int32_t numPropImgs = -1; // this property is not enabled
		if (m_properties[i])
		{
			numPropImgs = prop.size();
		
			if (numPropImgs != m_images.size())
			{
				setErrorString("Internal error: mismatch between number of images and properties");
				return false;
			}
		}

		if (!si.writeInt32(numPropImgs))
		{
			setErrorString("Error writing number of images according to property: " + si.getErrorString());
			return false;
		}

		for (size_t j = 0 ; j < prop.size() ; j++)
		{
			auto &propPts = prop[j];
			if (propPts.size() != m_images[j].size())
			{
				setErrorString("Internal error: number of image points does not match that of property");
				return false;
			}

			if (!si.writeInt32(propPts.size()) || !si.writeDoubles(propPts))
			{
				setErrorString("Unable to write property values: " + si.getErrorString());
				return false;
			}
		}
	}

	if (!si.writeInt32(m_groupPoints.size()))
	{
		setErrorString("Error writing number of groups: " + si.getErrorString());
		return false;
	}

	for (auto &grp : m_groupPoints)
	{
		if (!si.writeInt32(grp.size()) || !si.writeInt32s(grp))
		{
			setErrorString("Unable to write point group info: " + si.getErrorString());
			return false;
		}
	}

	if (!si.writeInt32(m_triangulations.size()))
	{
		setErrorString("Error writing number of triangulations: " + si.getErrorString());
		return false;
	}

	if (m_triangulations.size() != 0)
	{
		for (int i = 0 ; i < m_images.size() ; i++) // this is the amount of triangulations we should have
		{
			int32_t indices[3];
			
			// write size

			if (!si.writeInt32(m_triangulations[i].size()))
			{
				setErrorString("Error writing triangulation size: " + si.getErrorString());
				return false;
			}

			for (auto it = m_triangulations[i].begin() ; it != m_triangulations[i].end() ; it++)
			{
				indices[0] = (*it).getIndex(0);
				indices[1] = (*it).getIndex(1);
				indices[2] = (*it).getIndex(2);

				if (!si.writeInt32s(indices, 3))
				{
					setErrorString("Error writing triangulation indices: " + si.getErrorString());
					return false;
				}
			}
		}
	}

	if (!si.writeInt32(m_timeDelayInfo.size()))
	{
		setErrorString("Error writing number of time delay points: " + si.getErrorString());
		return false;
	}

	if (m_timeDelayInfo.size() != 0)
	{
		for (int i = 0 ; i < m_timeDelayInfo.size() ; i++)
		{
			int32_t indices[2];

			indices[0] = m_timeDelayInfo[i].getImageIndex();
			indices[1] = m_timeDelayInfo[i].getPointIndex();

			if (!si.writeInt32s(indices, 2))
			{
				setErrorString("Error writing time delay point indices: " + si.getErrorString());
				return false;
			}

			if (!si.writeDouble(m_timeDelayInfo[i].getTimeDelay()))
			{
				setErrorString("Error writing time delay: " + si.getErrorString());
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
			m_imagePointProperties[Intensity][img][imgpoint] -= v;
}

void ImagesData::clearTriangulation()
{
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

bool ImagesData::operator==(const ImagesData &d) const
{
	// auto printVec = [](const auto &v, auto f, const string &extra = string("  "))
	// {
	// 	cerr << extra << "[ ";
	// 	for (auto e : v)
	// 		cerr << f(e) << ",";
	// 	cerr << "]" << endl;
	// };

	// auto ident = [](auto x) { return x; };

	if (m_properties != d.m_properties)
	{
		// cerr << "Properties mismatch" << endl;
		// printVec(m_properties, ident);
		// printVec(d.m_properties, ident);
		return false;
	}
	if (m_images != d.m_images)
		return false;
	if (m_imagePointProperties != d.m_imagePointProperties)
		return false;
	if (m_groupPoints != d.m_groupPoints)
		return false;
	if (m_triangulations != d.m_triangulations)
		return false;
	if (m_timeDelayInfo != d.m_timeDelayInfo)
		return false;
	if (!(m_bottomLeft == d.m_bottomLeft))
		return false;
	if (!(m_topRight == d.m_topRight))
		return false;
	return true;
}

} // end namespace
