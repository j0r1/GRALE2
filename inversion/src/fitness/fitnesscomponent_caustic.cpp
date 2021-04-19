#include "fitnesscomponent_caustic.h"
#include "fitnessutil.h"
#include "imagesdataextended.h"

using namespace std;

namespace grale
{

// FitnessComponent_CausticPenalty

FitnessComponent_CausticPenalty::FitnessComponent_CausticPenalty(const std::shared_ptr<FitnessComponentCache> &pCache)
	: FitnessComponent("causticpenalty", pCache)
{
	addRecognizedTypeName("causticgrid");
}

FitnessComponent_CausticPenalty::~FitnessComponent_CausticPenalty()
{
}

bool FitnessComponent_CausticPenalty::inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence)
{
	string typeName;

	imgDat.getExtraParameter("type", typeName);
	if (typeName == "extendedimages") 
	{
		int numImg = imgDat.getNumberOfImages();
		if (numImg < 2)
		{
			setErrorString("Images data set must contain at least two images per source");
			return false;
		}

		for (int i = 0 ; i < numImg ; i++)
		{
			int numPoints = imgDat.getNumberOfImagePoints(i);

			if (numPoints < 1)
			{
				setErrorString("An empty image is present");
				return false;
			}
		}
		
		// Just record some information at this point

		needCalcDeflections = true;

		m_lastImgIdx = idx;
		m_lastImgDds = imgDat.getDds();
		m_lastImgDs = imgDat.getDs();
		m_lastImgNumImgs = numImg;

		return true;
	}

	if (typeName != "causticgrid") // ignore the rest
		return true;

	if (m_lastImgIdx < 0)
	{
		setErrorString("No previous extended image information that this caustic grid can refer to");
		return false;
	}

	// TODO: use another check? Dds doesn't make much sense in a multi-plane scenario
	if (m_lastImgDds != imgDat.getDds() || m_lastImgDs != imgDat.getDs())
	{
		setErrorString("Caustic grid distances are not the same as the ones from the last recorded extended images");
		return false;
	}

	if (!imgDat.hasTriangulation())
	{
		setErrorString("Caustic grid doesn't contain a triangulation");
		return false;
	}

	needCalcDeflections = true;
	needCalcDeflDeriv = true;
	needCalcInverseMag = true;

	int numImages = imgDat.getNumberOfImages();
	int sourcePos = m_sourceIndices.size();

	m_lineSegments.resize(sourcePos+1);
	m_lineSegmentFlags.resize(sourcePos+1);
	m_lineSegmentIntersections.resize(sourcePos+1);
	m_critTriangles.resize(sourcePos+1);

	m_lineSegments[sourcePos].resize(numImages);
	m_lineSegmentFlags[sourcePos].resize(numImages);
	m_lineSegmentIntersections[sourcePos].resize(numImages);
	m_critTriangles[sourcePos].resize(numImages);

	for (int i = 0 ; i < numImages ; i++)
	{
		vector<TriangleIndices> t;

		if (!imgDat.getTriangles(i, t))
		{
			setErrorString("Unable to obtain triangulation data from the caustic grid: " + imgDat.getErrorString());
			return false;
		}

		m_critTriangles[sourcePos][i].resize(t.size());

		int j = 0;
		vector<TriangleIndices>::const_iterator it;

		for (it = t.begin(); it != t.end() ; it++, j++)
		{
			TriangleIndices triangle = *it;

			int lineIndex1 = storeTriangleLine(sourcePos, i, triangle.getIndex(0), triangle.getIndex(1));
			int lineIndex2 = storeTriangleLine(sourcePos, i, triangle.getIndex(1), triangle.getIndex(2));
			int lineIndex3 = storeTriangleLine(sourcePos, i, triangle.getIndex(2), triangle.getIndex(0));
		
			m_critTriangles[sourcePos][i][j] = TriangleIndices(lineIndex1, lineIndex2, lineIndex3);
		}

		m_lineSegmentFlags[sourcePos][i].resize(m_lineSegments[sourcePos][i].size());
		m_lineSegmentIntersections[sourcePos][i].resize(m_lineSegments[sourcePos][i].size());
	}

	m_sourceIndices.push_back(m_lastImgIdx);
	m_critIndices.push_back(idx);

	addImagesDataIndex(m_lastImgIdx);
	addImagesDataIndex(idx);
	addToUsedImagesCount(m_lastImgNumImgs);
	//Note: disabled this again, so we can use the same image for caustics and null space
	//      for the same image
	//m_lastImgIdx = -1; // reset the referenced image
	return true;
}

int FitnessComponent_CausticPenalty::storeTriangleLine(int sourcePos, int imgIndex, int point1, int point2)
{
	int p1 = point1;
	int p2 = point2;

	if (p1 > p2)
	{
		p1 = point2;
		p2 = point1;
	}

	for (int i = 0 ; i < m_lineSegments[sourcePos][imgIndex].size() ; i++)
	{
		if (m_lineSegments[sourcePos][imgIndex][i].first == p1 && m_lineSegments[sourcePos][imgIndex][i].second == p2)
			return i;
	}

	int index = m_lineSegments[sourcePos][imgIndex].size();

	m_lineSegments[sourcePos][imgIndex].push_back(std::pair<int,int>(p1, p2));

	return index;
}

bool FitnessComponent_CausticPenalty::calculateFitness(const ProjectedImagesInterface &iface, float &fitness)
{
	fitness = calculateCausticPenaltyFitness(iface, m_sourceIndices, m_critIndices,
			                                 m_lineSegments, m_critTriangles, m_lineSegmentFlags,
											 m_lineSegmentIntersections, getCache());
	return true;
}

} // end namespace
