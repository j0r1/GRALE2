#include "lensinversionparametersmultiplanegpu.h"

using namespace std;

namespace grale
{

LensInversionParametersMultiPlaneGPU::LensInversionParametersMultiPlaneGPU()
	: m_massEstimate(0),
	  m_useSheets(false),
	  m_allowNeg(false),
	  m_deviceIdx(0)
{

}

LensInversionParametersMultiPlaneGPU::LensInversionParametersMultiPlaneGPU(const Cosmology &cosmology,
	const std::vector<double> &lensRedshifts,
	const std::vector<std::vector<std::shared_ptr<LensInversionBasisLensInfo>>> &basisLenses,
	const std::vector<std::shared_ptr<GravitationalLens>> &baseLensesPerPlane,
	const std::vector<std::shared_ptr<ImagesDataExtended>> &sourceImages,
	double massEstimate,
	bool useMassSheets,
	const ConfigurationParameters *pFitnessObjectParameters,
	bool allowNegativeWeights,
	const ScaleSearchParameters &massScaleSearchParams,
	int deviceIndex)
{
	m_cosmology = cosmology;
	m_lensRedshifts = lensRedshifts;

	auto createCopyImgData = [](shared_ptr<ImagesDataExtended> &img)
	{
		return make_shared<ImagesDataExtended>(*img);
	};
	auto createCopyBasisLens = [](shared_ptr<LensInversionBasisLensInfo> &bl)
	{
		shared_ptr<GravitationalLens> lensCopy(bl->m_pLens->createCopy());
		return make_shared<LensInversionBasisLensInfo>(lensCopy, bl->m_center, bl->m_relevantLensingMass);
	};
	auto createCopyLens = [](shared_ptr<GravitationalLens> &l)
	{
		shared_ptr<GravitationalLens> lensCopy;
		if (l.get())
			lensCopy = l->createCopy();
		return lensCopy;
	};

	auto deepCopyVector = [](const auto &srcVec, auto createCopy)
	{
		remove_const_t<remove_reference_t<decltype(srcVec)>> dstVec;

		for (auto x : srcVec)
		{
			auto copy = createCopy(x);
			dstVec.push_back(copy);
		}
		return dstVec;
	};

	for (auto &lp : basisLenses)
		m_basisLenses.push_back(deepCopyVector(lp, createCopyBasisLens));

	m_baseLensesPerPlane = deepCopyVector(baseLensesPerPlane, createCopyLens);

	m_images = deepCopyVector(sourceImages, createCopyImgData);

	m_massEstimate = massEstimate;
	m_useSheets = useMassSheets;
	if (pFitnessObjectParameters)
		m_fitnessObjectParams = make_shared<ConfigurationParameters>(*pFitnessObjectParameters);
	m_allowNeg = allowNegativeWeights;
	m_scaleSearchParams = massScaleSearchParams;
	m_deviceIdx = deviceIndex;
}


LensInversionParametersMultiPlaneGPU::~LensInversionParametersMultiPlaneGPU()
{

}

bool LensInversionParametersMultiPlaneGPU::write(serut::SerializationInterface &si) const
{
	if (m_lensRedshifts.size() != m_basisLenses.size())
	{
		setErrorString("Invalid configuration: the number of lens redshifts (" + to_string(m_lensRedshifts.size()) +
		               ") does not equal the number of lens planes with basis lenses (" + to_string(m_basisLenses.size()) + ")");
		return false;
	}
	if (!si.writeInt32((int32_t)m_lensRedshifts.size()) || !si.writeDoubles(m_lensRedshifts))
	{
		setErrorString("Unable to write the lens redshifts: " + si.getErrorString());
		return false;
	}
	if (!m_cosmology.write(si))
	{
		setErrorString("Can't write cosmology: " + m_cosmology.getErrorString());
		return false;
	}

	for (auto &plane : m_basisLenses)
	{
		if (!si.writeInt32((int32_t)plane.size()))
		{
			setErrorString("Can't write number of basis lenses in lens plane: " + si.getErrorString());
			return false;
		}
		for (auto &bl : plane)
		{
			assert(bl->m_pLens.get());
			double values[3] = { bl->m_center.getX(), bl->m_center.getY(), bl->m_relevantLensingMass };
			if (!bl->m_pLens->write(si) || !si.writeDoubles(values, 3))
			{
				setErrorString("Unable to write basis lens or its properties");
				return false;
			}
		}
	}

	if (!si.writeInt32((int32_t)m_baseLensesPerPlane.size()))
	{
		setErrorString("Can't write number of base lenses for the lens planes: " + si.getErrorString());
		return false;
	}
	for (auto &l : m_baseLensesPerPlane)
	{
		int32_t haveLensInt = (l.get() == nullptr)?0:1;
		if (!si.writeInt32(haveLensInt))
		{
			setErrorString("Can't write base lens flag: " + si.getErrorString());
			return false;
		}
		if (l.get())
		{
			if (!l->write(si))
			{
				setErrorString("Can't write base lens for a plane: " + l->getErrorString());
				return false;
			}
		}
	}

	if (!si.writeInt32((int32_t)m_images.size()))
	{
		setErrorString("Can't write number of images: " + si.getErrorString());
		return false;
	}
	for (auto &img : m_images)
	{
		assert(img.get());
		if (!img->write(si))
		{
			setErrorString("Can't write images data set: " + img->getErrorString());
			return false;
		}
	}

	int32_t intProps[3] = { 
		(m_useSheets)?1:0,
		(m_fitnessObjectParams.get())?1:0,
		(m_allowNeg)?1:0
	};

	if (!si.writeDouble(m_massEstimate) || !si.writeInt32s(intProps, 3))
	{
		setErrorString("Unable to write mass estimate or properties: " + si.getErrorString());
		return false;
	}

	if (m_fitnessObjectParams.get())
	{
		if (!m_fitnessObjectParams->write(si))
		{
			setErrorString("Unable to write fitness object parameters: " + m_fitnessObjectParams->getErrorString());
			return false;
		}
	}

	if (!m_scaleSearchParams.write(si))
	{
		setErrorString("Unable to write the mass scale search parameters: " + si.getErrorString());
		return false;
	}

	if (!si.writeInt32(m_deviceIdx))
	{
		setErrorString("Unable to write device index: " + si.getErrorString());
		return false;
	}
	return true;
}

bool LensInversionParametersMultiPlaneGPU::read(serut::SerializationInterface &si)
{
	int32_t numPlanes;
	if (!si.readInt32(&numPlanes))
	{
		setErrorString("Unable to read the number of lens planes: " + si.getErrorString());
		return false;
	}
	m_lensRedshifts.resize(numPlanes);
	if (!si.readDoubles(m_lensRedshifts))
	{
		setErrorString("Unable to read different lens plane redshifts: " + si.getErrorString());
		return false;
	}

	if (!m_cosmology.read(si))
	{
		setErrorString("Can't read cosmological model: " + m_cosmology.getErrorString());
		return false;
	}

	m_basisLenses.clear();
	for (int32_t i = 0 ; i < numPlanes ; i++)
	{
		m_basisLenses.push_back(vector<shared_ptr<LensInversionBasisLensInfo>>());
		auto &plane = m_basisLenses.back();

		int32_t numBasisLenses;
		if (!si.readInt32(&numBasisLenses))
		{
			setErrorString("Can't read number of basis lenses in a lens plane: " + si.getErrorString());
			return false;
		}
		
		for (int32_t j = 0 ; j < numBasisLenses ; j++)
		{
			string errStr;
			unique_ptr<GravitationalLens> pLens;

			if (!GravitationalLens::read(si, pLens, errStr))
			{
				setErrorString("Can't read a basis lens: " + errStr);
				return false;
			}

			shared_ptr<GravitationalLens> basisLens(move(pLens));
			double values[3];

			if (!si.readDoubles(values, 3))
			{
				setErrorString("Can't read basis lens properties: " + si.getErrorString());
				return false;
			}
			plane.push_back(make_shared<LensInversionBasisLensInfo>(basisLens, Vector2Dd(values[0], values[1]), values[2]));
		}
	}

	int32_t numBaseLensesPerPlane;
	if (!si.readInt32(&numBaseLensesPerPlane))
	{
		setErrorString("Can't read number of base lenses per plane: " + si.getErrorString());
		return false;
	}
	if (numBaseLensesPerPlane < 0 || numBaseLensesPerPlane > numPlanes) // sanity check
	{
		setErrorString("Invalid value for number of base lenses per plane: " + to_string(numBaseLensesPerPlane));
		return false;
	}

	m_baseLensesPerPlane.resize(numBaseLensesPerPlane);
	for (auto &l : m_baseLensesPerPlane)
	{
		int32_t haveBaseLens;
		if (!si.readInt32(&haveBaseLens))
		{
			setErrorString("Can't read base lens flag: " + si.getErrorString());
			return false;
		}
		if (haveBaseLens)
		{
			unique_ptr<GravitationalLens> lens;
			string errStr;
			if (!GravitationalLens::read(si, lens, errStr))
			{
				setErrorString("Can't read a base lens: " + errStr);
				return false;
			}
			l = move(lens);
		}
	}

	int32_t numImages;
	if (!si.readInt32(&numImages))
	{
		setErrorString("Can't read number of images: " + si.getErrorString());
		return false;
	}

	m_images.clear();
	for (int32_t i = 0 ; i < numImages ; i++)
	{
		shared_ptr<ImagesDataExtended> img = make_shared<ImagesDataExtended>();
		if (!img->read(si))
		{
			setErrorString("Unable to read images data set: " + img->getErrorString());
			return false;
		}
		m_images.push_back(img);
	}

	int32_t intProps[3];
	if (!si.readDouble(&m_massEstimate) || !si.readInt32s(intProps, 3))
	{
		setErrorString("Can't read mass estimate or properties: " + si.getErrorString());
		return false;
	}

	m_useSheets = (intProps[0] == 0)?false:true;
	bool haveFitnessObjectParams = (intProps[1] == 0)?false:true;
	m_allowNeg = (intProps[2] == 0)?false:true;

	if (haveFitnessObjectParams)
	{
		m_fitnessObjectParams = make_unique<ConfigurationParameters>();
		if (!m_fitnessObjectParams->read(si))
		{
			setErrorString("Unable to read fitness object parameters: " + m_fitnessObjectParams->getErrorString());
			return false;
		}
	}
	else
		m_fitnessObjectParams.reset();

	if (!m_scaleSearchParams.read(si))
	{
		setErrorString("Unable to read mass scale search parameters: " + m_scaleSearchParams.getErrorString());
		return false;
	}

	int32_t devIdx;
	if (!si.readInt32(&devIdx))
	{
		setErrorString("Unable to read device index: " + si.getErrorString());
		return false;
	}
	m_deviceIdx = devIdx;

	return true;
}

}
