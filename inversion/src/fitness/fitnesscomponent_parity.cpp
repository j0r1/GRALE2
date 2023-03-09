#include "fitnesscomponent_parity.h"
#include "projectedimagesinterface.h"

using namespace std;

namespace grale
{

FitnessComponent_ParityPenalty::FitnessComponent_ParityPenalty(const std::shared_ptr<FitnessComponentCache> &pCache)
	: FitnessComponent("paritypenalty", pCache)
{
	addRecognizedTypeName("paritygroups");
}

FitnessComponent_ParityPenalty::~FitnessComponent_ParityPenalty()
{
}

bool FitnessComponent_ParityPenalty::inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence)
{
	string typeName;

	imgDat.getExtraParameter("type", typeName);
	if (typeName != "paritygroups")
		return true; // Just ignore

	if (imgDat.getNumberOfGroups() < 1 || imgDat.getNumberOfGroups() > 2)
	{
		setErrorString("Need one or two point groups per images data set, to describe same or opposite parities");
		return false;
	}

	m_pointGroups.add(imgDat);
	addImagesDataIndex(idx);

	needCalcDeflDeriv = true;
	needCalcInverseMag = true;

	return true;
}

bool FitnessComponent_ParityPenalty::calculateFitness(const ProjectedImagesInterface &iface, float &fitness)
{
	auto &srcs = getUsedImagesDataIndices();
	assert(m_pointGroups.getNumberOfSources() == srcs.size());

	auto getInvMag = [this,&iface](int srcIdx, int srcIdxGroups, int group, int groupPoint) -> float
	{
		int im1, pt1;
		assert(srcIdxGroups >= 0 && srcIdxGroups < m_pointGroups.getNumberOfSources());
		m_pointGroups.getGroupPointIndices(srcIdxGroups, group, groupPoint, &im1, &pt1);

		assert(srcIdx >= 0 && srcIdx < iface.getNumberOfSources());
		assert(im1 >= 0 && im1 < iface.getNumberOfImages(srcIdx));
		assert(pt1 >= 0 && pt1 < iface.getNumberOfImagePoints(srcIdx, im1));
		return iface.getInverseMagnifications(srcIdx, im1)[pt1];
	};

	// TODO: for now, we're not doing anything to compensate for differnt amount
	//       of points per image or per source
	float total = 0;
	for (size_t i = 0 ; i < srcs.size() ; i++)
	{
		int srcIdx = srcs[i];
		int srcIdxGroups = i;
		int numGroups = m_pointGroups.getNumberOfGroups(srcIdxGroups);

		assert(numGroups > 0 && numGroups < 3);

		// For each individual group, make sure that points have same parity
		for (int g = 0 ; g < numGroups ; g++)
		{
			int numPts = m_pointGroups.getNumberOfGroupPoints(srcIdxGroups, g);
			// TODO: here we're checking every pair of points, is there something better?
			for (int p = 0 ; p < numPts-1 ; p++)
			{
				float mu1 = getInvMag(srcIdx, srcIdxGroups, g, p);

				for (int q = 0 ; q < numPts ; q++)
				{
					float mu2 = getInvMag(srcIdx, srcIdxGroups, g, q);

					if (mu1 * mu2 <= 0) // different sign, penalize
						total++;
				}
			}
		}

		if (numGroups == 2) // If there are two groups, check that they have different parity
		{
			int g1 = 0, g2 = 0;
			int g1Size = m_pointGroups.getNumberOfGroupPoints(srcIdxGroups, g1);
			int g2Size = m_pointGroups.getNumberOfGroupPoints(srcIdxGroups, g2);

			for (int p = 0 ; p < g1Size ; p++)
			{
				float mu1 = getInvMag(srcIdx, srcIdxGroups, g1, p);

				for (int q = 0 ; q < g2Size ; q++)
				{
					float mu2 = getInvMag(srcIdx, srcIdxGroups, g2, q);

					if (mu1 * mu2 >= 0) // same sign, penalize
						total++;
				}
			}
		}
	}

	fitness = total;
	return true;
}

}

