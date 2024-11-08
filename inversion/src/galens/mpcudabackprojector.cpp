#include "mpcudabackprojector.h"
#include "multiplanecuda.h"
#include "cosmology.h"

using namespace std;

namespace grale
{

MPCUDABackProjector::MPCUDABackProjector()
{
	m_pMPCU = nullptr;
}

MPCUDABackProjector::~MPCUDABackProjector()
{
}

bool MPCUDABackProjector::init(const string &libraryPath, int deviceIndex,
		const Cosmology &cosmology,
		const vector<float> &lensRedshifts,
		const vector<vector<PlummerLensInfo>> &lenses, 
		const vector<float> &sourceRedshifts,
		const vector<ImagesDataExtended *> &images)
{
	storeOriginalData(images);

	m_angularScale = getAngularScale(images);
	m_thetas.clear();

	// Store the image points, expressed in the angular scale
	for (auto &img : images)
	{
		const int numImg = img->getNumberOfImages();
		if (numImg == 0)
		{
			setErrorString("One of the images data sets does not contain any images");
			return false;
		}

		m_thetas.push_back(vector<Vector2Df>());
		vector<Vector2Df> &curThetas = m_thetas.back();

		for (int i = 0 ; i < numImg ; i++)
		{
			const int numPts = img->getNumberOfImagePoints(i);
			if (numPts == 0)
			{
				setErrorString("A specific image does not contain any points");
				return false;
			}

			for (int p = 0 ; p < numPts ; p++)
			{
				Vector2Dd pt = img->getImagePointPosition(i, p);
				pt /= m_angularScale;

				curThetas.push_back(Vector2Df((float)pt.getX(), (float)pt.getY()));
			}
		}
	}

	// Store the plummer parameters, with this angular unit
	vector<vector<MultiPlaneCUDA::PlummerInfo>> fixedPlummerParams;
	for (auto &plummers : lenses)
	{
		fixedPlummerParams.push_back(vector<MultiPlaneCUDA::PlummerInfo>());
		auto &curPlummers = fixedPlummerParams.back();

		for (auto p : plummers)
		{
			float scaledWidth = (float)(p.getAngularWidth()/m_angularScale);
			Vector2Dd scaledPos = p.getAngularPosition();
			scaledPos /= m_angularScale;

			curPlummers.push_back({ Vector2Df((float)scaledPos.getX(), (float)scaledPos.getY()),
									scaledWidth, p.getMass() });
		}
	}

	auto pMPCU = make_unique<MultiPlaneCUDA>();
	if (!pMPCU->init(libraryPath, deviceIndex, m_angularScale, cosmology.getH(),
					   cosmology.getOmegaM(), cosmology.getOmegaR(),
					   cosmology.getOmegaV(), cosmology.getW(),
					   lensRedshifts, fixedPlummerParams, sourceRedshifts, m_thetas))
	{
		setErrorString("Unable to initialize multi-plane CUDA based calculator: " + pMPCU->getErrorString());
		return false;
	}

	m_pMPCU = move(pMPCU);

	return true;
}

double MPCUDABackProjector::getAngularScale(const vector<ImagesDataExtended*> &images)
{
	double dx = 0, dy = 0;

	for (auto &img : images)
	{
		Vector2Dd tr = img->getTopRightCorner();
		Vector2Dd bl = img->getBottomLeftCorner();
		dx = MAX(dx, ABS(tr.getX()));
		dx = MAX(dx, ABS(bl.getX()));
		dy = MAX(dy, ABS(tr.getY()));
		dy = MAX(dy, ABS(bl.getY()));
	}
	return SQRT(dx*dx + dy*dy)/10.0;
}

bool MPCUDABackProjector::calculateSourcePositions(const vector<std::vector<float>> &massFactors,
												   const vector<float> &sheetDensities)
{
	if (!m_pMPCU)
	{
		setErrorString("Not initialized");
		return false;
	}
	
	if (!m_pMPCU->calculateSourcePositions(massFactors, sheetDensities))
	{
		setErrorString(m_pMPCU->getErrorString());
		return false;
	}

	return true;
}

const Vector2D<float> *MPCUDABackProjector::getBetas(int sourceNumber) const
{
	assert(m_pMPCU);
	assert(sourceNumber < m_thetas.size());
	const vector<Vector2Df> &srcPos = *(m_pMPCU->getSourcePositions(sourceNumber));
	return srcPos.data();
}

const Vector2D<float> *MPCUDABackProjector::getBetas(int sourceNumber, int imageNumber) const
{
	assert(m_pMPCU);
	assert(sourceNumber < m_thetas.size() && sourceNumber < m_offsets.size());
	assert(imageNumber < m_offsets[sourceNumber].size());

	int offset = m_offsets[sourceNumber][imageNumber];
	assert(offset < m_thetas[sourceNumber].size());
	
	const vector<Vector2Df> &srcPos = *(m_pMPCU->getSourcePositions(sourceNumber));
	return &(srcPos[offset]);
}

} // end namespace
