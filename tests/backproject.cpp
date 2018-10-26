#include "backprojectmatrixnew.h"
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

int main(void)
{
	RandomNumberGenerator rng;
	DeflectionMatrix matrix;
	vector<ImagesDataExtended *> images;
	vector<bool> ones;

	for (int i = 0 ; i < 10 ; i++) // images
	{
		double frac = (rng.pickRandomNumber()*0.4+0.5);
		ImagesDataExtended *pImg = new ImagesDataExtended(1.0, frac);

		int Nimg = 3;
		if (!pImg->create(Nimg, true, true))
			cerr << "Couldn't create image" << endl;

		for (int j = 0 ; j < Nimg ; j++)
		{
			for (int k = 0 ; k < 5 ; k++)
			{
				double x = (rng.pickRandomNumber()*2.0-1.0)*30*ANGLE_ARCSEC;
				double y = (rng.pickRandomNumber()*2.0-1.0)*30*ANGLE_ARCSEC;
				double intens = rng.pickRandomNumber();
				double shear1 = rng.pickRandomNumber();
				double shear2 = rng.pickRandomNumber();
				Vector2Dd pt(x, y);

				int idx = pImg->addPoint(j, pt, intens, shear1, shear2);

				double td = rng.pickRandomNumber()*50;
				pImg->addTimeDelayInfo(j, idx, td);
			}
		}
		images.push_back(pImg);
		ones.push_back(true);
	}

	if (!matrix.startInit())
		cerr << "Couldn't init matrix" << endl;

	double D_d = 900*DIST_MPC;
	double z_d = 0.7;

	GravitationalLens *pBaseLens = nullptr;
	if (rng.pickRandomNumber() < 0.5)
	{
		double x = (rng.pickRandomNumber()*2.0-1.0)*25*ANGLE_ARCSEC;
		double y = (rng.pickRandomNumber()*2.0-1.0)*25*ANGLE_ARCSEC;
		double w = (rng.pickRandomNumber()*0.25+0.75)*20*ANGLE_ARCSEC;
		double mass = (rng.pickRandomNumber()*0.9+0.1)*MASS_SOLAR*1e11;
		MultiplePlummerLensParams params({PlummerLensInfo(mass, w, {x, y})});

		pBaseLens = new MultiplePlummerLens();
		if (!pBaseLens->init(D_d, &params))
			cerr << "Couldn't init base lens" << endl;
		cerr << "Using base lens" << endl;
	}
	else
		cerr << "Not using base lens" << endl;

	bool useMassSheet = (rng.pickRandomNumber() < 0.5)?true:false;
	if (useMassSheet)
		cerr << "Using mass sheet" << endl;
	else
		cerr << "Not using mass sheet" << endl;

	BackProjectMatrixNew bpMatrix;
	
	if (!bpMatrix.startInit(z_d, D_d, &matrix, images, ones, ones, ones, pBaseLens, true, true, true, useMassSheet))
		cerr << "Couldn't init BackProjectMatrixNew" << endl;

	vector<pair<GravitationalLens *, Vector2Dd>> basisLenses;
	for (int i = 0 ; i < 10 ; i++) // lenses
	{
		double x = (rng.pickRandomNumber()*2.0-1.0)*15*ANGLE_ARCSEC;
		double y = (rng.pickRandomNumber()*2.0-1.0)*15*ANGLE_ARCSEC;
		Vector2Dd ctr(x,y);

		double w = (rng.pickRandomNumber()*0.9+0.1)*7*ANGLE_ARCSEC;
		double mass = (rng.pickRandomNumber()*0.9+0.1)*MASS_SOLAR*1e13;
		PlummerLensParams params(mass, w);
		GravitationalLens *pLens = new PlummerLens();
		if (!pLens->init(D_d, &params))
			cerr << "Couldn't init lens" << endl;

		basisLenses.push_back({pLens, ctr});
	}

	if (!matrix.endInit(basisLenses))
		cerr << "Couldn't end matrix init" << endl;

	if (!bpMatrix.endInit())
		cerr << "Couldn't end BackProjectMatrixNew" << endl;

	vector<float> weights(basisLenses.size());
	for (auto &w : weights)
		w = (float)(rng.pickRandomNumber()*0.5+0.5);

	if (!matrix.calculateBasisMatrixProducts(weights, true, true, true))
		cerr << "Couldn't calculate matrix product" << endl;

	bpMatrix.storeDeflectionMatrixResults();
	double scaleFactor = rng.pickRandomNumber()+0.5;
	double sheetFactor = rng.pickRandomNumber()+0.3;
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
			cerr << "Couldn't init mass sheet lens" << endl;
		compParams.addLens(sheetFactor, Vector2Dd(0,0), 0, l);
	}
	CompositeLens compLens;
	if (!compLens.init(D_d, &compParams))
		cerr << "Couldn't init composite lens" << endl;

	list<ImagesDataExtended *> imageList;
	for (auto i : images)
		imageList.push_back(i);

	ImagesBackProjector bpSlow(compLens, imageList, z_d, false);

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

	cout << "max diff = " << maxDiff << " arcsec" << endl;

	return 0;
}
