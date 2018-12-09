#include "inversioncommunicator.h"
#include "inputoutput.h"
#include "galensmodule.h"
#include "lensfitnessobject.h"
#include "imagesbackprojector.h"
#include "gridlensinversiongafactoryparams.h"
#include "gravitationallens.h"
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

	bool_t runModule(const string &moduleDir, const string &moduleFile, GALensModule *pModule);
};

bool_t CalcFitnessCommunicator::runModule(const string &moduleDir, const string &moduleFile, GALensModule *pModule)
{
	bool_t r;
	int popSize = 0;

	vector<uint8_t> factoryParamBytes;
	if (!(r = readLineAndBytesWithPrefix("GAFACTORYPARAMS", factoryParamBytes, 10000)))
		return "Error reading GA factory parameters: " + r.getErrorString();

	GridLensInversionGAFactoryParams factoryParams;
	if (!(r = loadFromBytes(factoryParams, factoryParamBytes)))
		return "Unable to load GA factory parameters from received data: " + r.getErrorString();

	// We're just abusing the factory parameters to get the images, check
	// that we're using some dummy settings to be sure that this is the purpose
	if (factoryParams.getMaximumNumberOfGenerations() != 1 ||
		factoryParams.getBasisLenses().size() != 1)
		return "Expecting dummy settings for max number of generations and grid";

	unique_ptr<LensFitnessObject> f(pModule->createFitnessObject());
	if (!f.get())
		return "Unable to create lens fitness object for module: " + moduleFile;

	list<ImagesDataExtended *> dummyShortList;
	// list can be modified
	list<ImagesDataExtended *> images;
	for (auto img : factoryParams.getImages())
		images.push_back(img.get());

	double z_d = factoryParams.getZ_d();
	if (!f->init(z_d, images, dummyShortList, factoryParams.getFitnessObjectParameters()))
		return "Unable to initialize fitness object: " + f->getErrorString();

	string fitnessDesc = f->getFitnessComponentsDescription();
	WriteLineStdout("FITDESC:" + fitnessDesc);

	// if a lens is set (we're abusing the base lens in the parameters for this), we'll
	// calculate the fitness as well
	const GravitationalLens *pLens = factoryParams.getBaseLens();
	if (pLens)
	{
		// We'll abuse the base lens as the lens for which we're calculating the fitness

		unique_ptr<GravitationalLens> lens(pLens->createCopy());
		ImagesBackProjector bp(*(lens.get()), images, z_d, false);
		f->postInit(images, dummyShortList, bp.getAngularScale());

		vector<float> fitnessComp(f->getNumberOfFitnessComponents());
		f->calculateOverallFitness(bp, &fitnessComp[0]); // TODO: do this in background thread and write PING/PONG stuff?

		stringstream ss;
		ss << fitnessComp[0];
		for (size_t i = 1 ; i < fitnessComp.size() ; i++)
			ss << "," << fitnessComp[i];

		WriteLineStdout("FITNESS:" + ss.str());
	}
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
