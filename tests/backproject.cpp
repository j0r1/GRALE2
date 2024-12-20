#include "backprojectmatrix.h"
#include "imagesbackprojector.h"
#include "randomnumbergenerator.h"
#include "plummerlens.h"
#include "multipleplummerlens.h"
#include "constants.h"
#include "imagesdataextended.h"
#include "compositelens.h"
#include "masssheetlens.h"
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

int main(void)
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
		MassSheetLensParams params(D_d, 1, 1);
		sheetLens = shared_ptr<GravitationalLens>(new MassSheetLens());
		if (!sheetLens->init(D_d, &params))
			throw runtime_error("Unable to initialize sheet lens: " + sheetLens->getErrorString());
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
	{
		MassSheetLensParams params(((SPEED_C*SPEED_C)/(4.0*CONST_PI*CONST_G*D_d)));
		MassSheetLens l;
		if (!l.init(D_d, &params))
			throw runtime_error("Couldn't init mass sheet lens");
		compParams.addLens(sheetFactor, Vector2Dd(0,0), 0, l);
	}
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

	int nSrc = bpMatrix.getNumberOfSources();
	double maxDiff = 0;
	for (int s = 0 ; s < nSrc ; s++)
	{
		int npts = bpMatrix.getNumberOfImagePoints(s);
		process(npts, bpMatrix.getBetas(s), angScaleMatrix, bpSlow.getBetas(s), angScaleSlow, maxDiff);
		process(npts, bpMatrix.getAlphas(s), angScaleMatrix, bpSlow.getAlphas(s), angScaleSlow, maxDiff);
	}
	cout << "max diff alpha/beta = " << maxDiff << " arcsec" << endl;

	maxDiff = 0;
	double pScaleMatrix = bpMatrix.getAngularScale()*bpMatrix.getAngularScale()/(ANGLE_ARCSEC*ANGLE_ARCSEC);
	double pScaleSlow = bpSlow.getAngularScale()*bpSlow.getAngularScale()/(ANGLE_ARCSEC*ANGLE_ARCSEC);
	for (int s = 0 ; s < nSrc ; s++)
	{
		int npts = bpMatrix.getNumberOfImagePoints(s);
		process(npts, bpMatrix.getLensPotential(s), pScaleMatrix, bpSlow.getLensPotential(s), pScaleSlow, maxDiff);
	}
	cout << "max diff potential = " << maxDiff << " arcsec^2" << endl;

	maxDiff = 0;
	for (int s = 0 ; s < nSrc ; s++)
	{
		int npts = bpMatrix.getNumberOfImagePoints(s);
		process(npts, bpMatrix.getSecondDerivativesXXX(s), 1.0/angScaleMatrix, 
		        bpSlow.getSecondDerivativesXXX(s), 1.0/angScaleSlow, maxDiff);
		process(npts, bpMatrix.getSecondDerivativesYYY(s), 1.0/angScaleMatrix, 
		        bpSlow.getSecondDerivativesYYY(s), 1.0/angScaleSlow, maxDiff);
		process(npts, bpMatrix.getSecondDerivativesXXY(s), 1.0/angScaleMatrix, 
		        bpSlow.getSecondDerivativesXXY(s), 1.0/angScaleSlow, maxDiff);
		process(npts, bpMatrix.getSecondDerivativesYYX(s), 1.0/angScaleMatrix, 
		        bpSlow.getSecondDerivativesYYX(s), 1.0/angScaleSlow, maxDiff);
	}

	cout << "max diff a_xxx,a_yyy,a_xxy,a_yyx = " << maxDiff << " arcsec^-1" << endl;

	vector<float> tdMatrix, tdSlow;
	auto getTimeDelays = [](ProjectedImagesInterface &iface, vector<float> &v)
	{
		int nsrc = iface.getNumberOfSources();
		for (int s = 0 ; s < nsrc ; s++)
		{
			int nimg = iface.getNumberOfImages(s);
			for (int i = 0 ; i < nimg ; i++)
			{
				int npts = iface.getNumberOfImagePoints(s, i);
				for (int p = 0 ; p < npts ; p++)
				{
					double t = iface.getTimeDelay(s, i, p, {0,0});
					v.push_back(t);
				}
			}
		}
	};
	
	getTimeDelays(bpMatrix, tdMatrix);
	getTimeDelays(bpSlow, tdSlow);

	maxDiff = 0;
	process(tdMatrix.size(), tdMatrix.data(), 1.0/365.0, tdSlow.data(), 1.0/365.0, maxDiff);
	cout << "max diff = " << maxDiff << " years" << endl;

	auto getRelativeTimeDelaysTDFormula = [D_d,z_d](ProjectedImagesInterface &iface, vector<float> &v)
	{
		int nsrc = iface.getNumberOfSources();
		for (int s = 0 ; s < nsrc ; s++)
		{
			int nimg = iface.getNumberOfImages(s);
			Vector2Df beta0 = iface.getBetas(s, 0)[0];
			double t0 = iface.getTimeDelay(s, 0, 0, beta0);
			for (int i = 0 ; i < nimg ; i++)
			{
				int npts = iface.getNumberOfImagePoints(s, i);
				for (int p = 0 ; p < npts ; p++)
				{
					Vector2Df beta = iface.getBetas(s, i)[p];
					double t = iface.getTimeDelay(s, i, p, beta);
					v.push_back(t-t0);
				}
			}
		}
	};

	auto getRelativeTimeDelaysDeflectionAngle = [D_d,z_d](ProjectedImagesInterface &iface, vector<float> &v)
	{
		int nsrc = iface.getNumberOfSources();
		for (int s = 0 ; s < nsrc ; s++)
		{
			int nimg = iface.getNumberOfImages(s);
			Vector2Df alpha0 = iface.getAlphas(s, 0)[0];
			float phi0 = iface.getLensPotential(s, 0)[0];
			double X = (D_d*(1.0+z_d)/SPEED_C)/iface.getDistanceFraction(s) * iface.getAngularScale()*iface.getAngularScale();
			X /= (60*60*24); // convert seconds to days

			for (int i = 0 ; i < nimg ; i++)
			{
				int npts = iface.getNumberOfImagePoints(s, i);
				for (int p = 0 ; p < npts ; p++)
				{
					Vector2Df alpha1 = iface.getAlphas(s, i)[p];
					float phi1 = iface.getLensPotential(s, i)[p];
					double phiDiff = (phi1-phi0);
					double alphaDiff = 0.5*(alpha1.getLengthSquared() - alpha0.getLengthSquared());

					double td = (alphaDiff - phiDiff)*X;
					v.push_back(td);
				}
			}
		}
	};

	maxDiff = 0;
	vector<float> tdRel1, tdRel2;
	getRelativeTimeDelaysTDFormula(bpSlow, tdRel1);
	getRelativeTimeDelaysDeflectionAngle(bpSlow, tdRel2);
	process(tdRel1.size(), tdRel1.data(), 1.0/365, tdRel2.data(), 1.0/365, maxDiff);
	cout << "max diff = " << maxDiff << " days" << endl;

	maxDiff = 0;
	tdRel1.clear();
	tdRel2.clear();
	getRelativeTimeDelaysTDFormula(bpMatrix, tdRel1);
	getRelativeTimeDelaysDeflectionAngle(bpMatrix, tdRel2);
	process(tdRel1.size(), tdRel1.data(), 1.0/365, tdRel2.data(), 1.0/365, maxDiff);
	cout << "max diff = " << maxDiff << " days" << endl;

	return 0;
}
