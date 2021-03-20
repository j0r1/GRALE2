#include "fitnessutil.h"
#include "pointgroupstorage.h"
#include "fitnesscomponent.h"
#include "projectedimagesinterface.h"
#include <assert.h>
#include <limits>

using namespace std;

namespace grale
{

float getScaleFactor_PointImages(const ProjectedImagesInterface &interface,
		                         const std::vector<int> &sourceIndices)
{
	vector<int> sourceGroup(sourceIndices.size(), 0);
	vector<float> scaleFactors(1, 0);

	ScaleFactorWorkspace ws;
	ws.m_floats.resize(4);

	getScaleFactors_PointImages(interface, sourceIndices, sourceGroup, scaleFactors, ws, MinMax);
	float scale = scaleFactors[0];
	assert(scale != 0);
	return scale;
}

void getScaleFactors_PointImages(const ProjectedImagesInterface &interface,
		                         const vector<int> &sourceIndices,
								 const vector<int> &sourceGroup,
								 vector<float> &scaleFactors, // one for each group
								 ScaleFactorWorkspace &ws,
								 PointImageScaleType scaleType)
{
	float defaultScale = (interface.getAngularScale()*interface.getAngularScale())/(ANGLE_ARCSEC*ANGLE_ARCSEC);
	for (auto &s : scaleFactors)
		s = defaultScale;

	assert(scaleFactors.size() > 0);
	assert(sourceGroup.size() == sourceIndices.size());

	if (sourceIndices.size() > 1) // We're assuming that different source indices are present
	{
		if (scaleType == MinMax)
		{
			vector<float> &workSpace = ws.m_floats;

			workSpace.resize(scaleFactors.size()*4);
			for (int i = 0 ; i < scaleFactors.size() ; i++)
			{
				float minx = numeric_limits<float>::max();
				float maxx = -numeric_limits<float>::max();
				float miny = numeric_limits<float>::max();
				float maxy = -numeric_limits<float>::max();

				workSpace[i*4+0] = minx;
				workSpace[i*4+1] = maxx;
				workSpace[i*4+2] = miny;
				workSpace[i*4+3] = maxy;
			}

			for (int sIdx = 0 ; sIdx < sourceIndices.size() ; sIdx++)
			{
				const int s = sourceIndices[sIdx];
				const int group = sourceGroup[sIdx];

				assert(group >= 0 && group < scaleFactors.size());
				assert(s < interface.getNumberOfSources());
				const int numImages = interface.getNumberOfImages(s);

				if (numImages > 1)
				{
					for (int i = 0 ; i < numImages ; i++)
					{
						Vector2D<float> beta = interface.getBetas(s, i)[0];
						float x = beta.getX();
						float y = beta.getY();
						
						float &minx = workSpace[group*4+0];
						float &maxx = workSpace[group*4+1];
						float &miny = workSpace[group*4+2];
						float &maxy = workSpace[group*4+3];

						maxx = MAX(maxx, x);
						minx = MIN(minx, x);
						maxy = MAX(maxy, y);
						miny = MIN(miny, y);
					}
				}
			}

			for (int group = 0 ; group < scaleFactors.size() ; group++)
			{
				float minx = workSpace[group*4+0];
				float maxx = workSpace[group*4+1];
				float miny = workSpace[group*4+2];
				float maxy = workSpace[group*4+3];

				float dx = maxx-minx;
				float dy = maxy-miny;

				if (dx > 0 && dy > 0)
				{
					float scale = 1.0f/(dx*dx + dy*dy);
					assert(!isnan(scale));

					scaleFactors[group] = scale;
				}
			}
		}
		else // MAD
		{
			assert(scaleType == MAD);

			vector<vector<float>> &vecWorkspace = ws.m_vecFloat;
			const int numGroups = scaleFactors.size();

			vecWorkspace.resize(numGroups*2); // one for x components, one for y

			for (int group = 0 ; group < scaleFactors.size() ; group++)
			{
				auto &xPoints = vecWorkspace[group*2+0];
				auto &yPoints = vecWorkspace[group*2+1];
				xPoints.resize(0);
				yPoints.resize(0);
			}

			for (int sIdx = 0 ; sIdx < sourceIndices.size() ; sIdx++)
			{
				const int s = sourceIndices[sIdx];
				const int group = sourceGroup[sIdx];

				assert(group >= 0 && group < scaleFactors.size());
				assert(s < interface.getNumberOfSources());
				const int numImages = interface.getNumberOfImages(s);

				auto &xPoints = vecWorkspace[group*2+0];
				auto &yPoints = vecWorkspace[group*2+1];

				if (numImages > 1)
				{
					for (int i = 0 ; i < numImages ; i++)
					{
						Vector2D<float> beta = interface.getBetas(s, i)[0];
						float x = beta.getX();
						float y = beta.getY();
						
						xPoints.push_back(x);
						yPoints.push_back(y);
					}
				}
			}

			for (int group = 0 ; group < scaleFactors.size() ; group++)
			{
				auto &xPoints = vecWorkspace[group*2+0];
				auto &yPoints = vecWorkspace[group*2+1];
				size_t len = xPoints.size();
				size_t len2 = len/2;

				float xMed = 0, yMed = 0;

				for (int k = 0 ; k < 2 ; k++)
				{
					sort(xPoints.begin(), xPoints.end());
					sort(yPoints.begin(), yPoints.end());

					if (len % 2 == 1) // odd
					{
						xMed = xPoints[len2];
						yMed = yPoints[len2];
					}
					else // even
					{
						xMed = 0.5f*(xPoints[len2] + xPoints[len2-1]);
						yMed = 0.5f*(yPoints[len2] + yPoints[len2-1]);
					}

					if (k == 0)
					{
						for (size_t i = 0 ; i < len ; i++)
						{
							xPoints[i] = std::abs(xPoints[i]-xMed);
							yPoints[i] = std::abs(yPoints[i]-yMed);
						}
					}
				}

				// Actually the square of the scale is used in the other routine
				float scaleSquared = 1.0f/(xMed*xMed + yMed*yMed);
				assert(!isnan(scaleSquared));

				scaleFactors[group] = scaleSquared;
			}
		}
	}
}

float calculateOverlapFitness_PointImages(const ProjectedImagesInterface &interface, 
		                                  const std::vector<int> &sourceIndices,
										  float scale)
{
	vector<int> sourceGroup(sourceIndices.size(), 0);
	vector<float> scaleFactors(1, scale);

	return calculateOverlapFitness_PointImages(interface, sourceIndices, sourceGroup, scaleFactors);
}

float calculateOverlapFitness_PointImages(const ProjectedImagesInterface &interface, 
		                                  const vector<int> &sourceIndices,
										  const vector<int> &sourceGroups,
										  const vector<float> &scaleFactors)
{
	float fitness = 0;
	int sourceCount = 0;

	assert(sourceGroups.size() == sourceIndices.size());
	assert(scaleFactors.size() > 0);

	for (int sIdx = 0 ; sIdx < sourceIndices.size() ; sIdx++)
	{
		const int s = sourceIndices[sIdx];
		const int group = sourceGroups[sIdx];

		assert(group >= 0 && group < scaleFactors.size());
		assert(s < interface.getNumberOfSources());
		const int numImages = interface.getNumberOfImages(s);
		const float scale = scaleFactors[group];

		if (numImages > 1)
		{
			float sourceFitness = 0;
			
			for (int i = 0 ; i < numImages ; i++)
			{
				Vector2D<float> beta1 = interface.getBetas(s, i)[0];

				for (int j = i+1 ; j < numImages ; j++)
				{
					Vector2D<float> beta2 = interface.getBetas(s, j)[0];
					Vector2D<float> diff = beta2-beta1;

					sourceFitness += diff.getLengthSquared();
				}
			}

			int numSprings = (numImages*(numImages-1))/2;
			
			sourceFitness /= (float)numSprings;
			fitness += sourceFitness * scale;
			sourceCount++; // Only count the sources that are actually used
		}
	}

	if (sourceCount > 0)
		fitness /= (float)sourceCount;

	assert(!isnan(fitness));
	return fitness;
}

float calculateOverlapFitness_PointGroups(const PointGroupStorage &pointGroups,
		                                  const ProjectedImagesInterface &interface,
										  const std::vector<int> &sourceIndices,
										  PointGroupRMSType t,
										  // I added this so that the same calculation
										  // can be used in the bayesian fitness
										  vector<Vector2Df> *pAllThetaDiffs
										  )
{
	assert(sourceIndices.size() == pointGroups.getNumberOfSources());
	
	int numContribs = 0;
	float fitness = 0;

	vector<Vector2Df> betas;
	if (pAllThetaDiffs)
		pAllThetaDiffs->clear();

	for (int sIdx = 0 ; sIdx < sourceIndices.size() ; sIdx++)
	{
		const int src = sourceIndices[sIdx];

		const int groups = pointGroups.getNumberOfGroups(sIdx);
		for (int g = 0 ; g < groups ; g++)
		{
			int numPoints = pointGroups.getNumberOfGroupPoints(sIdx, g);
			if (numPoints <= 1) // can't do anything with only one point
				continue;

			float groupRms = 0;

			auto getGroupPointBetaAndDerivs = [src, sIdx, g, &pointGroups, &interface](int idx,  
					                                       float &axx, float &ayy, float &axy,
														   bool withDerivs = true) -> Vector2Df
			{
				int imgNum = 0, ptNum = 0;
				pointGroups.getGroupPointIndices(sIdx, g, idx, &imgNum, &ptNum);

				assert(imgNum >= 0 && imgNum < interface.getNumberOfImages(src));
				assert(ptNum >= 0 && ptNum < interface.getNumberOfImagePoints(src, imgNum));

				if (withDerivs)
				{
					axx = interface.getDerivativesXX(src, imgNum)[ptNum];
					ayy = interface.getDerivativesYY(src, imgNum)[ptNum];
					axy = interface.getDerivativesXY(src, imgNum)[ptNum];
				}

				return interface.getBetas(src, imgNum)[ptNum];
			};

			auto getGroupPointBeta = [getGroupPointBetaAndDerivs](int idx) -> Vector2Df
			{
				float dummy;
				return getGroupPointBetaAndDerivs(idx, dummy, dummy, dummy, false);
			};

			auto getBetas_All = [numPoints, getGroupPointBeta](vector<Vector2Df> &betas)
			{
				betas.clear();
				for (int p = 0 ; p < numPoints ; p++)
					betas.push_back(getGroupPointBeta(p));
			};

			auto getBetas_Average = [numPoints, getGroupPointBeta](vector<Vector2Df> &betas)
			{
				betas.clear();
				Vector2Df avgBeta { 0, 0 };
				for (int p = 0 ; p < numPoints ; p++)
					avgBeta += getGroupPointBeta(p);
				avgBeta /= numPoints;
				betas.push_back(avgBeta);
			};

			if (t == AllBetas)
				getBetas_All(betas);
			else
			{
				assert(t == AverageBeta);
				getBetas_Average(betas);
			}

			int ptCount = 0;
			for (int p1 = 0 ; p1 < numPoints ; p1++)
			{
				// The point theta0 corresponds to a source plane position beta0
				float axx0, ayy0, axy0;
				Vector2Df beta0 = getGroupPointBetaAndDerivs(p1, axx0, ayy0, axy0);

				float det = (1.0f-axx0)*(1.0f-ayy0)-axy0*axy0;
				const float epsilon = 1e-10; // to avoid division by zero
				if (ABS(det) < epsilon)
					det = copysign(epsilon, det);
				const float invDet = 1.0f/det;

				// inverse of magnification matrix
				const float AI[2][2] = { { invDet*(1.0f-ayy0), invDet*axy0 },
					                     { invDet*axy0, (1.0f-axx0)*invDet } };

				for (auto dBeta : betas)
				{
					ptCount++;

					dBeta -= beta0;

					// The beta corresponding to this point (if close to beta0),
					// would correspond to a theta that differs by this much from
					// theta0
					Vector2Df dTheta { AI[0][0]*dBeta.getX() + AI[0][1]*dBeta.getY(),
						               AI[1][0]*dBeta.getX() + AI[1][1]*dBeta.getY() };

					// TODO: optionally use something else here?
					groupRms += dTheta.getLengthSquared();

					if (pAllThetaDiffs)
						pAllThetaDiffs->push_back(dTheta);
				}
			}

			groupRms /= ptCount;

			fitness += groupRms;
			numContribs++;
		}
	}

	if (numContribs > 0)
		fitness = SQRT(fitness/numContribs);
	
	// Let's express this in arcsec, to keep things understandable
	float s = (float)(interface.getAngularScale()/ANGLE_ARCSEC);
	return fitness*s;
}

float calculateOverlapFitness_Extended(const PointGroupStorage &pointGroups, const ProjectedImagesInterface &iface,
		                               const vector<int> &sourceIndices,
									   const std::vector<bool> &rectFlags,
									   const std::vector<bool> &groupFlags,
									   FitnessComponentCache *pCache
									   )
{
	float posfitness = 0;
	int sourceCount = 0;

	float avgImgScale = 0;
	int imgScaleCount = 0;

	assert(pointGroups.getNumberOfSources() == sourceIndices.size());
	assert(rectFlags.size() == sourceIndices.size());
	assert(groupFlags.size() == sourceIndices.size());

	vector<bool> processedSource(sourceIndices.size(), false);

	// In the first iteration, we'll handle the extended images
	// In the second iteration we'll process the point images in a similar
	// manner, but with the source scale based on the extended images
	for (int it = 0 ; it < 2 ; it++)
	{
		for (int sIdx = 0 ; sIdx < sourceIndices.size() ; sIdx++)
		{
			if (processedSource[sIdx])
				continue;

			const int s = sourceIndices[sIdx];
			//std::cerr << "source = " << s << endl;

			assert(s < iface.getNumberOfSources());
			int numimages = iface.getNumberOfImages(s);
			if (numimages < 2) // Sources with a single image may be present, they can be used in null space
				continue;

			float sourcefitness = 0;
			vector<Vector2D<float> > bottomleftpoints(numimages);
			vector<Vector2D<float> > toprightpoints(numimages);
			float scale = 0;

			bool pointImage = false;
			for (int i = 0 ; i < numimages ; i++)
			{
				int numpoints = iface.getNumberOfImagePoints(s,i);
				if (numpoints == 1)
					pointImage = true;

				const Vector2D<float> *betas = iface.getBetas(s,i);

				float maxx = betas[0].getX();
				float minx = betas[0].getX();
				float maxy = betas[0].getY();
				float miny = betas[0].getY();

				//cerr << "x = " << betas[0].getX() << " " << "y = " << betas[0].getY() << endl;
				
				for (int p = 1 ; p < numpoints ; p++)
				{
					maxx = MAX(maxx,betas[p].getX());
					maxy = MAX(maxy,betas[p].getY());
					minx = MIN(minx,betas[p].getX());
					miny = MIN(miny,betas[p].getY());
					
					//cerr << "x = " << betas[p].getX() << " " << "y = " << betas[p].getY() << endl;
				}

				bottomleftpoints[i] = Vector2D<float>(minx,miny);
				toprightpoints[i] = Vector2D<float>(maxx,maxy);
				float imgscale = SQRT((maxx-minx)*(maxx-minx)+(maxy-miny)*(maxy-miny));

				scale += imgscale;

			}

			scale /= (float)numimages;
			//cerr << "scale for " << s << " is " << scale << endl;

			if (it == 0)
			{
				if (pointImage) // handle point images when the scale has been set
				{
					//cerr << "is point image, skipping for now" << endl;
					continue;
				}

				avgImgScale += scale;
				imgScaleCount++;
			}
			else
			{
				assert(!pointImage); // in the second iteraction we should only process pointimages
				scale = avgImgScale;
				//cerr << "forcing scale to " << scale << endl;

				if (scale == 0)
					return std::numeric_limits<float>::quiet_NaN();
			}

			// Cache the scale for use by the caustic penalty
			if (pCache)
				pCache->setEstimatedSourceScale(s, scale);

			int numsprings = 0;

			// Surrounding square
			if (rectFlags[sIdx])
			{
				for (int i = 0 ; i < numimages ; i++)
				{
					Vector2D<float> corneri1(bottomleftpoints[i].getX(),toprightpoints[i].getY());
					Vector2D<float> corneri2(toprightpoints[i].getX(),toprightpoints[i].getY());
					Vector2D<float> corneri3(toprightpoints[i].getX(),bottomleftpoints[i].getY());
					Vector2D<float> corneri4(bottomleftpoints[i].getX(),bottomleftpoints[i].getY());
						
					for (int j = i+1 ; j < numimages ; j++)
					{
						Vector2D<float> cornerj1(bottomleftpoints[j].getX(),toprightpoints[j].getY());
						Vector2D<float> cornerj2(toprightpoints[j].getX(),toprightpoints[j].getY());
						Vector2D<float> cornerj3(toprightpoints[j].getX(),bottomleftpoints[j].getY());
						Vector2D<float> cornerj4(bottomleftpoints[j].getX(),bottomleftpoints[j].getY());
						
						Vector2D<float> diff1 = (corneri1-cornerj1)/scale;
						Vector2D<float> diff2 = (corneri2-cornerj2)/scale;
						Vector2D<float> diff3 = (corneri3-cornerj3)/scale;
						Vector2D<float> diff4 = (corneri4-cornerj4)/scale;

						//cerr << "between " << i << " and " << j << ": " << diff1.getLengthSquared() << " "
						//	<< diff2.getLengthSquared() << " "
						//	<< diff3.getLengthSquared() << " "
						//	<< diff4.getLengthSquared() << endl;
						sourcefitness += diff1.getLengthSquared()+diff2.getLengthSquared()+diff3.getLengthSquared()+diff4.getLengthSquared();
					}
				}

				numsprings += 4*(numimages*(numimages-1))/2;
			}
		
			// Point groups (if any)
			if (groupFlags[sIdx])
			{
				assert(sIdx < pointGroups.getNumberOfSources());
				int ng = pointGroups.getNumberOfGroups(sIdx);

				for (int g = 0 ; g < ng ; g++)
				{
					int np = pointGroups.getNumberOfGroupPoints(sIdx,g);
					
					for (int p = 0 ; p < np ; p++)
					{
						int img,point;

						pointGroups.getGroupPointIndices(sIdx,g,p,&img,&point);

						assert(s >= 0 && s < iface.getNumberOfSources());
						assert(img >= 0 && img < iface.getNumberOfImages(s));
						assert(point >= 0 && point < iface.getNumberOfImagePoints(s, img));
						Vector2D<float> point1 = iface.getBetas(s,img)[point];
						
						for (int q = p+1 ; q < np ; q++)
						{
							pointGroups.getGroupPointIndices(sIdx,g,q,&img,&point);
							assert(img >= 0 && img < iface.getNumberOfImages(s));
							assert(point >= 0 && point < iface.getNumberOfImagePoints(s, img));
							Vector2D<float> point2 = iface.getBetas(s,img)[point];

							Vector2D<float> diff = (point2-point1)/scale;

							sourcefitness += diff.getLengthSquared();
						}
					}
					
					numsprings += (np*(np-1))/2;
				}
			}

			if (numsprings != 0)
				sourcefitness /= (float)numsprings;

			//cerr << "sourcefitness: " << sourcefitness << endl;
			posfitness += sourcefitness;
			sourceCount++;

			processedSource[sIdx] = true;

			if (it == 0)
				avgImgScale /= imgScaleCount;

			//std::cerr << "sourcefitness = " << sourcefitness << std::endl;	
		}
	}

	//std::cerr << std::endl;

	if (sourceCount > 0)
		posfitness /= sourceCount;
	
	assert(!isnan(posfitness));
	return posfitness;
}

// Create one point group for the point images
shared_ptr<ImagesDataExtended> addGroupsToPointImages(const ImagesDataExtended &imgDat, string &errStr)
{
	shared_ptr<ImagesDataExtended> Bad;
	shared_ptr<ImagesDataExtended> grpDat0 { new ImagesDataExtended(imgDat) };
	ImagesData &grpDat = *grpDat0.get();
	bool hasPt = false;

	if (grpDat.getNumberOfGroups() != 0)
	{
		errStr = "Image already contains point groups";
		return Bad;
	}

	int newGrp = grpDat.addGroup();

	for (int img = 0 ; img < grpDat.getNumberOfImages() ; img++)
	{
		int numPts = grpDat.getNumberOfImagePoints(img);
		if (numPts != 1)
		{
			errStr = "Image id " + to_string(img) + " has " + to_string(numPts) + " points (should be exactly 1)";
			return Bad;
		}

		grpDat.addGroupPoint(newGrp, img, 0);
		hasPt = true;
	}

	if (!hasPt)
	{
		errStr = "No points were present, no groups could be added";
		return Bad;
	}
	return grpDat0;
}

} // end namespace
