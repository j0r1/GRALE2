#include "lensfitnessobject.h"
#include "backprojectmatrix.h"
#include "imagesbackprojector.h"
#include "precalculatedbackprojector.h"
#include "randomnumbergenerator.h"
#include "plummerlens.h"
#include "multipleplummerlens.h"
#include "constants.h"
#include "imagesdataextended.h"
#include "compositelens.h"
#include "masssheetlens.h"
#include "lensfitnessobject.h"
#include "lensinversionparameterssingleplanecpu.h"
#include <vector>
#include <iostream>
#include <memory>
#include <stdexcept>

using namespace grale;
using namespace std;

void process(int N, const Vector2Df *v1, double s1, const Vector2Df *v2, double s2, double &maxDiff)
{
	for (int i = 0 ; i < N ; i++)
	{
		double b1x = v1[i].getX()*s1;
		double b1y = v1[i].getY()*s1;
		double b2x = v2[i].getX()*s2;
		double b2y = v2[i].getY()*s2;
		//cout << b1x << " " << b2x << " " << b1x-b2x << endl;
		//cout << b1y << " " << b2y << " " << b1y-b2y << endl;
		double dx = std::abs(b1x-b2x);
		double dy = std::abs(b1y-b2y);
		if (dx > maxDiff)
			maxDiff = dx;
		if (dy > maxDiff)
			maxDiff = dy;
	}
}

void process(int N, const float *v1, double s1, const float *v2, double s2, double &maxDiff)
{
	for (int i = 0 ; i < N ; i++)
	{
		double x1 = v1[i]*s1;
		double x2 = v2[i]*s2;
		//cout << x1 << " " << x2 << " " << x1-x2 << " " << x1/x2 << endl;
		double d = std::abs(x1-x2);
		if (d > maxDiff)
			maxDiff = d;
	}
}

int main(int argc, char *argv[])
{
	RandomNumberGenerator rng;
	DeflectionMatrix matrix;
	vector<ImagesDataExtended *> images;
	vector<bool> ones;
#if 1
	const int NSOURCES = 10;
	const int NLENSES = 10;
#else
	const int NSOURCES = 1;
	const int NLENSES = 1;
#endif

	for (int i = 0 ; i < NSOURCES ; i++) // images
	{
		double frac = (rng.getRandomDouble()*0.4+0.5);
		ImagesDataExtended *pImg = new ImagesDataExtended(1.0, frac);

		int Nimg = 3;
		if (!pImg->create(Nimg, true, true))
			throw runtime_error("Couldn't create image");

		pImg->setExtraParameter("type", string("extendedimages"));

		for (int j = 0 ; j < Nimg ; j++)
		{
			for (int k = 0 ; k < 5 ; k++)
			{
				double x = (rng.getRandomDouble()*2.0-1.0)*30*ANGLE_ARCSEC;
				double y = (rng.getRandomDouble()*2.0-1.0)*30*ANGLE_ARCSEC;
				double intens = rng.getRandomDouble();
				double shear1 = rng.getRandomDouble();
				double shear2 = rng.getRandomDouble();
				Vector2Dd pt(x, y);
				Vector2Dd shear(shear1, shear2);

				int idx = pImg->addPoint(j, pt, intens, shear);

				double td = rng.getRandomDouble()*50;
				pImg->addTimeDelayInfo(j, idx, td);
			}
		}
		images.push_back(pImg);
		ones.push_back(true);
	}

	if (!matrix.startInit())
		throw runtime_error("Couldn't init matrix");

	double D_d = 900*DIST_MPC;
	double z_d = 0.7;

	GravitationalLens *pBaseLens = nullptr;
	if (rng.getRandomDouble() < 0.5)
	{
		double x = (rng.getRandomDouble()*2.0-1.0)*25*ANGLE_ARCSEC;
		double y = (rng.getRandomDouble()*2.0-1.0)*25*ANGLE_ARCSEC;
		double w = (rng.getRandomDouble()*0.25+0.75)*20*ANGLE_ARCSEC;
		double mass = (rng.getRandomDouble()*0.9+0.1)*MASS_SOLAR*1e11;
		MultiplePlummerLensParams params({PlummerLensInfo(mass, w, {x, y})});

		pBaseLens = new MultiplePlummerLens();
		if (!pBaseLens->init(D_d, &params))
			throw runtime_error("Couldn't init base lens");
		cerr << "Using base lens" << endl;
	}
	else
		cerr << "Not using base lens" << endl;

	bool useMassSheet = (rng.getRandomDouble() < 0.5)?true:false;
	shared_ptr<GravitationalLens> sheetLens;
	if (useMassSheet)
	{
		cerr << "Using mass sheet" << endl;
		sheetLens = LensInversionParametersSinglePlaneCPU::createDefaultSheetLens(LensInversionParametersSinglePlaneCPU::Genome, D_d); 
	}
	else
		cerr << "Not using mass sheet" << endl;

	BackProjectMatrix bpMatrix;
	
	if (!bpMatrix.startInit(z_d, D_d, &matrix, images, ones, ones, ones, ones, pBaseLens, sheetLens.get()))
		throw runtime_error("Couldn't init BackProjectMatrixNew");

	vector<pair<shared_ptr<GravitationalLens>, Vector2Dd>> basisLenses;
	for (int i = 0 ; i < NLENSES ; i++) // lenses
	{
		double x = (rng.getRandomDouble()*2.0-1.0)*15*ANGLE_ARCSEC;
		double y = (rng.getRandomDouble()*2.0-1.0)*15*ANGLE_ARCSEC;
		Vector2Dd ctr(x,y);

		double w = (rng.getRandomDouble()*0.9+0.1)*7*ANGLE_ARCSEC;
		double mass = (rng.getRandomDouble()*0.9+0.1)*MASS_SOLAR*1e13;
		PlummerLensParams params(mass, w);
		shared_ptr<GravitationalLens> pLens(new PlummerLens());
		if (!pLens->init(D_d, &params))
			throw runtime_error("Couldn't init lens");

		basisLenses.push_back({pLens, ctr});
	}

	if (!matrix.endInit(basisLenses))
		throw runtime_error("Couldn't end matrix init");

	if (!bpMatrix.endInit())
		throw runtime_error("Couldn't end BackProjectMatrixNew");

	vector<float> weights(basisLenses.size());
	for (auto &w : weights)
		w = (float)(rng.getRandomDouble()*0.5+0.5);

	if (!matrix.calculateBasisMatrixProducts(weights, true, true, true, true))
		throw runtime_error("Couldn't calculate matrix product");

	bpMatrix.storeDeflectionMatrixResults();
	double scaleFactor = rng.getRandomDouble()+0.5;
	double sheetFactor = rng.getRandomDouble()+0.3;
	bpMatrix.calculate(scaleFactor, sheetFactor);

	// Build corresponding composite lens
	CompositeLensParams compParams;
	for (int i = 0 ; i < weights.size() ; i++)
		compParams.addLens(weights[i]*scaleFactor, basisLenses[i].second, 0, *(basisLenses[i].first));

	if (pBaseLens)
		compParams.addLens(1.0, Vector2Dd(0,0), 0, *pBaseLens);

	if (useMassSheet)
		compParams.addLens(sheetFactor, Vector2Dd(0,0), 0, *sheetLens.get());

	auto compLens = make_shared<CompositeLens>();
	if (!compLens->init(D_d, &compParams))
		throw runtime_error("Couldn't init composite lens");

	list<ImagesDataExtended *> imageList;
	for (auto i : images)
		imageList.push_back(i);

	ImagesBackProjector bpSlow(compLens, imageList, z_d);

	double angScaleMatrix = bpMatrix.getAngularScale()/ANGLE_ARCSEC;
	cout << "angScaleMatrix = " << angScaleMatrix << endl; 
	double angScaleSlow = bpSlow.getAngularScale()/ANGLE_ARCSEC;
	cout << "angScaleSlow = " << angScaleSlow << endl; 

	vector<ImagesData *> imgs;
	vector<ImagesData *> srcs;
	for (size_t s = 0 ; s < images.size() ; s++)
	{
		ImagesData *pImg = images[s];
		ImagesData *pSrc = new ImagesData();
		imgs.push_back(pImg);
		srcs.push_back(pSrc);

		pSrc->create(pImg->getNumberOfImages(), false, false);
		for (int i = 0 ; i < pImg->getNumberOfImages() ; i++)
		{
			for (int p = 0 ; p < pImg->getNumberOfImagePoints(i) ; p++)
			{
				Vector2Df b = bpSlow.getBetas(s, i)[p];
				pSrc->addPoint(i, Vector2Dd(b.getX(), b.getY())*angScaleSlow);
			}
		}
	}

	PreCalculatedBackProjector preBp;
	if (!preBp.init(imgs, srcs))
		throw runtime_error("Couldn't init PreCalculatedBackProjector: " + preBp.getErrorString());

	double angScalePrecalc = preBp.getAngularScale()/ANGLE_ARCSEC;
	cout << "angScalePrecalc = " << angScalePrecalc << endl; 

	if (argc != 2)
	{
		cerr << "Specify a lens module" << endl;
		return -1;
	}

	unique_ptr<LensFitnessObject> pLFO = LensFitnessObjectRegistry::instance().createFitnessObject(argv[1]);
	auto pConfParams = pLFO->getDefaultParametersInstance();

	list<ImagesDataExtended *> imagesList;
	for (auto i : images)
		imagesList.push_back(i);
	list<ImagesDataExtended *> dummy;

	if (!pLFO->init(z_d, imagesList, dummy, pConfParams.get()))
	{
		cerr << "Unable to initialize lens fitness object: " << pLFO->getErrorString() << endl;
		return -1;
	}

	float fitnessMatrix, fitnessSlow, fitnessPrecalc;
	if (!pLFO->calculateMassScaleFitness(bpMatrix, fitnessMatrix))
		throw runtime_error("Can't calculate fitness based on matrix");
	if (!pLFO->calculateMassScaleFitness(bpSlow, fitnessSlow))
		throw runtime_error("Can't calculate fitness based on slow");
	if (!pLFO->calculateMassScaleFitness(bpMatrix, fitnessPrecalc))
		throw runtime_error("Can't calculate fitness based on precalc bp");

	cerr << "fitnessMatrix: " << fitnessMatrix << endl;
	cerr << "fitnessSlow: " << fitnessSlow << endl;
	cerr << "fitnessPrecalc:" << fitnessPrecalc << endl;
	return 0;

}
