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

#ifndef GRALE_IMAGESDATA_H

#define GRALE_IMAGESDATA_H

#include "graleconfig.h"
#include "vector2d.h"
#include "triangle2d.h"
#include "triangleindices.h"
#include <serut/serializationinterface.h>
#include <sys/types.h>
#include <stdint.h>
#include <string>
#include <vector>

namespace grale
{

class GRALE_IMPORTEXPORT ImagesData : public errut::ErrorBase
{
public:
	enum PropertyName
	{
		Intensity = 0,
		ShearComponent1,
		ShearComponent2,
		ShearWeight,
		DistanceFraction,
		ShearComponent1Uncertainty,
		ShearComponent2Uncertainty,
		MaxProperty
	};

	ImagesData();
	ImagesData(const ImagesData &data)						{ copyFrom(data); }
	virtual ~ImagesData();
	ImagesData *createCopy() const;

	bool create(int numImages, bool intensities, bool shearInfo);
	bool create(int numImages, const std::vector<PropertyName> &properties);

	int addImage();
	int addPoint(int imageNumber, Vector2Dd point, const std::vector<std::pair<PropertyName, double>> &properties);

	// TODO: remove these later
	int addPoint(int imageNumber, Vector2D<double> point); // returns index of the point, -1 on error
	int addPoint(int imageNumber, Vector2D<double> point, double intensity); // returns index of the point, -1 on error
	int addPoint(int imageNumber, Vector2D<double> point, Vector2D<double> shearComponents, double shearWeight = 1);
	int addPoint(int imageNumber, Vector2D<double> point, double intensity, Vector2D<double> shearComponents, double shearWeight = 1);

	int addGroup(); // returns the index of the group, -1 on error
	bool addGroupPoint(int groupNumber, int imageIndex, int pointIndex);
	bool addTimeDelayInfo(int imageIndex, int pointIndex, double timeDelay); // timeDelay in days
	bool addTriangle(int imageNumber, int index1, int index2, int index3);
	
	bool load(const std::string &fname);
	bool save(const std::string &fname) const;
	bool read(serut::SerializationInterface &si);
	bool write(serut::SerializationInterface &si) const;

	bool hasProperty(PropertyName n) const				{ if (n < 0 || n >= m_properties.size()) return false; return m_properties[n]; }
	// TODO: remove these, for now they're set for backwards compatibility
	bool hasIntensities() const							{ return hasProperty(Intensity); }
	bool hasShearInfo() const							{ return (hasProperty(ShearComponent1) && hasProperty(ShearComponent2) && hasProperty(ShearWeight)); }

	bool hasTimeDelays() const							{ if (m_timeDelayInfo.empty()) return false; return true; }

	int getNumberOfImages() const 							{ return m_images.size(); }
	int getNumberOfImagePoints(int i) const						{ return m_images[i].size(); }
	Vector2D<double> getImagePointPosition(int image, int point) const		{ return m_images[image][point]; }
	double getImagePointProperty(PropertyName n, int image, int point) const { return m_imagePointProperties[n][image][point]; }

	// TODO: remove these
	double getImagePointIntensity(int image, int point) const			{ return getImagePointProperty(Intensity, image, point); }
	double getShearComponent1(int image, int point) const				{ return getImagePointProperty(ShearComponent1, image, point); }
	double getShearComponent2(int image, int point) const				{ return getImagePointProperty(ShearComponent2, image, point); }
	double getShearWeight(int image, int point) const					{ return getImagePointProperty(ShearWeight, image, point); }

	void setImagePointPosition(int image, int point, Vector2D<double> position) { m_images[image][point] = position; findExtremes(); }

	bool hasTriangulation()	const							{ if (m_triangulations.size() != 0) return true; return false; }
	bool getTriangles(int image, std::vector<TriangleIndices> &triangles) const;
	void clearTriangulation();

	int getNumberOfGroups() const							{ return m_groupPoints.size(); }
	int getNumberOfGroupPoints(int group) const					{ return m_groupPoints[group].size()/2; }
	void getGroupPointIndices(int group, int pointnr, int *img, int *point) const	{ *img = m_groupPoints[group][pointnr*2]; *point = m_groupPoints[group][pointnr*2+1]; }
	void clearGroupInfo()								{ m_groupPoints.clear(); }

	int getNumberOfTimeDelays() const						{ return m_timeDelayInfo.size(); }
	void getTimeDelay(int index, int *pImg, int *pPoint, double *pDelay) const	{ *pImg = m_timeDelayInfo[index].getImageIndex(); *pPoint = m_timeDelayInfo[index].getPointIndex(); *pDelay = m_timeDelayInfo[index].getTimeDelay(); }
	void clearTimeDelayInfo() 							{ m_timeDelayInfo.clear(); }

	Vector2D<double> getTopRightCorner() const					{ return m_topRight; }
	Vector2D<double> getBottomLeftCorner() const					{ return m_bottomLeft; }

	void centerOnPosition(double ra, double dec);
	void uncenterOnPosition(double ra, double dec);
	void subtractIntensity(double v);

	const ImagesData &operator=(const ImagesData &dat)				{ clear(); copyFrom(dat); return *this; }
	bool operator==(const ImagesData &d) const;
private:
	bool readOld(serut::SerializationInterface &si);

	void clear();
	void findExtremes();
	void copyFrom(const ImagesData &data);

	class TimeDelayPoint
	{
	public:
		TimeDelayPoint(int imageIndex, int pointIndex, double timeDelay)
		{
			m_imageIndex = imageIndex;
			m_pointIndex = pointIndex;
			m_timeDelay = timeDelay;
		}

		int getImageIndex() const						{ return m_imageIndex; }
		int getPointIndex() const						{ return m_pointIndex; }
		double getTimeDelay() const						{ return m_timeDelay; }
		bool operator==(const TimeDelayPoint &src) const
		{ 
			return (m_imageIndex == src.m_imageIndex && 
			        m_pointIndex == src.m_pointIndex && 
					m_timeDelay == src.m_timeDelay);
		}
	private:
		int m_imageIndex;
		int m_pointIndex;
		double m_timeDelay;
	};

	std::vector<bool> m_properties;
	std::vector<std::vector<Vector2D<double> > > m_images;
	std::vector<std::vector<std::vector<double>>> m_imagePointProperties;
	std::vector<std::vector<int32_t> > m_groupPoints;
	std::vector<std::vector<TriangleIndices> > m_triangulations;
	std::vector<TimeDelayPoint> m_timeDelayInfo;
	Vector2D<double> m_bottomLeft;
	Vector2D<double> m_topRight;
};

} // end namespace

#endif // GRALE_IMAGESDATA_H

