#include "lensinversionparametersmultiplanegpu.h"
#include "cosmology.h"
#include "constants.h"
#include "plummerlens.h"
#include "randomnumbergenerator.h"
#include "lensinversiongafactorymultiplanegpu.h"
#include "lensfitnessgeneral.h"
#include "lensgaindividual.h"
#include <eatk/multithreadedpopulationfitnesscalculation.h>
#include <eatk/population.h>
#include <iostream>

using namespace std;
using namespace grale;
using namespace errut;

Cosmology cosm;
auto pRng = make_shared<RandomNumberGenerator>();
RandomNumberGenerator &rng = *pRng;

shared_ptr<GravitationalLens> plummer(double zd, double mass, double width)
{
	auto lens = make_shared<PlummerLens>();
	double Dd = cosm.getAngularDiameterDistance(zd);
	PlummerLensParams params { mass, width };
	if (!lens->init(Dd, &params))
		throw runtime_error("Can't init Plummer lens: " + lens->getErrorString());
	return lens;
}

shared_ptr<LensInversionBasisLensInfo> positionedPlummer(double zd, double mass, double width, Vector2Dd pos)
{
	return make_shared<LensInversionBasisLensInfo>(plummer(zd, mass, width), pos, mass);
}

shared_ptr<ImagesDataExtended> createRandomImages()
{
	auto img = make_shared<ImagesDataExtended>();
	size_t num = rng.getRandomUint32()%8 + 2;
	//size_t num = 2;
	if (!img->create(num, {}))
		throw runtime_error("Can't create image");
	
	for (size_t i = 0 ; i < num ; i++)
	{
		double x = ((rng.getRandomDouble()-0.5)*2)* (10*ANGLE_ARCSEC);
		double y = ((rng.getRandomDouble()-0.5)*2)* (10*ANGLE_ARCSEC);
		if (img->addPoint(i, {x, y}) < 0)
			throw runtime_error("Can't add point: " + img->getErrorString());
	}

	double z = rng.getRandomDouble()*6.0;
	img->setExtraParameter("z", z);
	img->setExtraParameter("type", string("pointimages"));
	img->save("test.imgdata");
	return img;
}

int main(void)
{	
	vector<double> zds { 0.5, 0.8, 1.5 };
	vector<vector<shared_ptr<LensInversionBasisLensInfo>>> basisLenses {
		// Lenses for first plane
		{
			positionedPlummer(zds[0], 1e13*MASS_SOLAR, 1*ANGLE_ARCSEC, { 1*ANGLE_ARCSEC, 0}),
			positionedPlummer(zds[0], 1e13*MASS_SOLAR, 1.5*ANGLE_ARCSEC, { -1*ANGLE_ARCSEC, 0})
		},
		{
			positionedPlummer(zds[1], 2e13*MASS_SOLAR, 2*ANGLE_ARCSEC, { 0, 1*ANGLE_ARCSEC }),
			positionedPlummer(zds[1], 2e13*MASS_SOLAR, 2.5*ANGLE_ARCSEC, { 0, -1*ANGLE_ARCSEC })
		},
		{
			positionedPlummer(zds[1], 3e13*MASS_SOLAR, 3*ANGLE_ARCSEC, { 1*ANGLE_ARCSEC, 1*ANGLE_ARCSEC }),
			positionedPlummer(zds[1], 3e13*MASS_SOLAR, 3.5*ANGLE_ARCSEC, { -1*ANGLE_ARCSEC, -1*ANGLE_ARCSEC })
		}
	};
	vector<shared_ptr<ImagesDataExtended>> images(rng.getRandomUint32()%15 + 5); // 20 sources
	for (auto &i : images)
		i = createRandomImages();

	auto fitObj = make_unique<LensFitnessGeneral>();
	auto confParams = fitObj->getDefaultParametersInstance();
	ScaleSearchParameters searchParams(false);
	LensInversionParametersMultiPlaneGPU invParams(cosm, zds, basisLenses, images, 
	                                               1e14*MASS_SOLAR, false, confParams.get(),
												   false, searchParams, 0);

	bool_t r;
	size_t numThreads = 4;
	vector<shared_ptr<eatk::GenomeFitnessCalculation>> calculators(numThreads);

	size_t numBasisFunction, numSheets, numObj;
	bool allowNeg;

	for (auto &gc : calculators)
	{
		auto genomeCalculator = make_shared<LensInversionGAFactoryMultiPlaneGPU>(make_unique<LensFitnessGeneral>());
		if (!(r = genomeCalculator->init(invParams)))
			throw runtime_error("Can't init multi-plane genome calculator: " + r.getErrorString());
	
		numBasisFunction = genomeCalculator->getNumberOfBasisFunctions();
		numSheets = genomeCalculator->getNumberOfSheets();
		numObj = genomeCalculator->getNumberOfObjectives();
		allowNeg = genomeCalculator->allowNegativeValues();

		gc = genomeCalculator;
	}

	cout << "#basisfunctions = " << numBasisFunction << endl;
	cout << "#sheets = " << numSheets << endl;
	cout << "#obj = " << numObj << endl;
	cout << "allowneg = " << ((allowNeg)?"true":"false") << endl;

	// Create a population

	LensGAIndividualCreation indivCreation(pRng, numBasisFunction, numSheets, allowNeg, numObj);
	size_t numPops = 1;
	size_t popSize = 16;

	auto createPop = [&indivCreation](size_t popSize) -> shared_ptr<eatk::Population>
	{
		auto pop = make_shared<eatk::Population>();
		for (size_t i = 0 ; i < popSize ; i++)
		{
			auto ind = indivCreation.createReferenceIndividual();
			ind->genome() = indivCreation.createInitializedGenome();
			ind->fitness() = indivCreation.createEmptyFitness();
			pop->append(ind);
		}
		return pop;
	};

	vector<shared_ptr<eatk::Population>> populations(numPops);
	for (auto &pop : populations)
		pop = createPop(popSize);
	
	// Calculate the fitness of its members

	eatk::MultiThreadedPopulationFitnessCalculation fitCalc;

	if (!(r = fitCalc.initThreadPool(calculators)))
		throw runtime_error("Can't init thread pool: " + r.getErrorString());

	const int numCalcs = 4;
	for (int i = 0 ; i < numCalcs ; i++)
	{
		cout << "Calculating " << (i+1) << "/" << numCalcs << endl; 
	
		if (!(r = fitCalc.calculatePopulationFitness(populations)))
			throw runtime_error("Can't calculate fitnesses: " + r.getErrorString());
		for (auto &pop : populations)
			for (auto &i : pop->individuals())
				i->fitness()->setCalculated(false);
	}

	fitCalc.destroyThreadPool();
	return 0;
}

