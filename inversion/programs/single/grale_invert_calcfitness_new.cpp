#include "inversioncommunicatornewga.h"
#include "inputoutput.h"
#include "lensfitnessobject.h"
#include "imagesbackprojector.h"
#include "precalculatedbackprojector.h"
#include "lensinversionparameterssingleplanecpu.h"
#include "gravitationallens.h"
#include "imagesdataextended.h"
#include "configurationparameters.h"
#include <errut/booltype.h>
#include <memory>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;
using namespace grale;
using namespace errut;

class CalcFitnessCommunicator : public InversionCommunicator
{
public:
	CalcFitnessCommunicator() { }
	~CalcFitnessCommunicator() { }
protected:
	string getVersionInfo() const override { return "Inversion module usage information"; }

	bool_t runModule(const std::string &lensFitnessObjectType,
					 std::unique_ptr<grale::LensFitnessObject> fitnessObject,
					 const std::string &calculatorType);
};

class GravitationalLensWrapper : public errut::ErrorBase
{
public:
	GravitationalLensWrapper() { }
	~GravitationalLensWrapper() { }

	bool read(serut::SerializationInterface &s)
	{
		unique_ptr<GravitationalLens> pLens;
		string errorString;

		if (!GravitationalLens::read(s, pLens, errorString))
		{
			setErrorString(errorString);
			return false;
		}

		m_lens = move(pLens);
		return true;
	}

	const shared_ptr<GravitationalLens> get() const { return m_lens; }
private:
	shared_ptr<GravitationalLens> m_lens;
};

bool_t CalcFitnessCommunicator::runModule(const std::string &lensFitnessObjectType, 
					 std::unique_ptr<grale::LensFitnessObject> fitnessObject,
					 const std::string &calculatorType)
{
	bool_t r;
	string type;
	if (!(r= readLineWithPrefix("TYPE", type, 10000)))
			return "Error reading fitness info/calculation type: " + r.getErrorString();

	vector<uint8_t> factoryParamBytes;
	if (!(r = readLineAndBytesWithPrefix("GAFACTORYPARAMS", factoryParamBytes, 10000)))
		return "Error reading GA factory parameters: " + r.getErrorString();

	LensInversionParametersSinglePlaneCPU factoryParams;
	if (!(r = loadFromBytes(factoryParams, factoryParamBytes)))
		return "Unable to load GA factory parameters from received data: " + r.getErrorString();

	// We're just abusing the factory parameters to get the images, check
	// that we're using some dummy settings to be sure that this is the purpose
	if (factoryParams.getBasisLenses().size() != 1 ||
		factoryParams.getBaseLens() != nullptr)
		return "Expecting dummy settings for grid and base lens";

	auto &f = fitnessObject;
	if (!f.get())
		return "Unable to create lens fitness object for module: " + lensFitnessObjectType;

	list<ImagesDataExtended *> dummyShortList;
	// list can be modified
	list<ImagesDataExtended *> images;
	for (auto img : factoryParams.getImages())
		images.push_back(img.get());

	double z_d = factoryParams.getZ_d();
	if (!f->init(z_d, images, dummyShortList, factoryParams.getFitnessObjectParameters()))
		return "Unable to initialize fitness object: " + f->getErrorString();
		  
	vector<string> unusedKeys;
	factoryParams.getFitnessObjectParameters()->getUnretrievedKeys(unusedKeys);
	if (unusedKeys.size() > 0)
		return "Key '" + unusedKeys[0] + "' in configuration parameters is not used"; 

	string fitnessDesc = f->getFitnessComponentsDescription();
	WriteLineStdout("FITDESC:" + fitnessDesc);

	if (type == "fitnessdescription")
	{
		// Nothing to do
	}
	else if (type == "precalculated" || type == "lens")
	{
		unique_ptr<ProjectedImagesInterface> iface;
		GravitationalLensWrapper lensWrapper;

		if (type == "lens")
		{
			vector<uint8_t> lensParams;
			if (!(r = readLineAndBytesWithPrefix("LENS", lensParams, 10000)))
				return "Error reading GA factory parameters: " + r.getErrorString();

			if (!(r = loadFromBytes(lensWrapper, lensParams)))
				return "Error loading lens: " + lensWrapper.getErrorString();

			iface = make_unique<ImagesBackProjector>(lensWrapper.get(), images, z_d);
		}
		else // precalculated
		{
			vector<ImagesData *> imagesVector;
			for (auto img : factoryParams.getImages())
				imagesVector.push_back(img.get());

			vector<ImagesData *> bpImagesVector;
			vector<shared_ptr<ImagesData>> bpImagesVectorSharedPtr;

			for (size_t i = 0 ; i < imagesVector.size() ; i++)
			{
				vector<uint8_t> imgDataBytes;
				if (!(r = readLineAndBytesWithPrefix("IMGDATA", imgDataBytes, 10000)))
					return "Error reading backprojected image data bytes: " + r.getErrorString();

				shared_ptr<ImagesData> imgDat = make_shared<ImagesData>();
				if (!(r = loadFromBytes(*(imgDat.get()), imgDataBytes)))
					return "Couldn't load images data set from bytes: " + r.getErrorString();
				bpImagesVectorSharedPtr.push_back(imgDat);
				bpImagesVector.push_back(imgDat.get());
			}

			auto pBp = make_unique<PreCalculatedBackProjector>();
			if (!pBp->init(imagesVector, bpImagesVector))
				return "Couldn't initialize pre-calculated backprojector interface: " + pBp->getErrorString();

			iface = move(pBp);
		}

		vector<float> fitnessComp(f->getNumberOfFitnessComponents());
		if (!f->calculateOverallFitness(*iface.get(), &fitnessComp[0])) // TODO: do this in background thread and write PING/PONG stuff?
			return "Error calculating fitness: " + f->getErrorString();

		stringstream ss;
		ss << fitnessComp[0];
		for (size_t i = 1 ; i < fitnessComp.size() ; i++)
			ss << "," << fitnessComp[i];

		WriteLineStdout("FITNESS:" + ss.str());
	}
	else
		return "Invalid type '" + type + "' for fitness info/calculation";

	WriteLineStdout("DONE");

	return true;
}

int main(int argc, char *argv[])
{
	CalcFitnessCommunicator comm;

	bool_t r = comm.run();
	if (!r)
	{
		cerr << "ERROR: " << r.getErrorString() << endl;
		return -1;
	}

	return 0;
}
