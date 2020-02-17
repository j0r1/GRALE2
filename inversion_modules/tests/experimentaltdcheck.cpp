/*
 * This is a test file for an experimental time delay fitness measure, to
 * make sure that it reaches a near-zero value when the correct lens is used.
 *
 * First, for a certain lens model, a random source position is used, the
 * corresponding images are calculated and the time delays for these images
 * are determined. This is then used as input for fake lenses, which are
 * simply rescaled (in mass density value) versions of the original lens,
 * so for a scale factor of 1, the correct lens is again used and the time
 * delay fitness values should be (nearly) zero at that point. 
 *
 * The output should be saved to 'dat' and the experimentaltdcheck_plot.py
 * file can be used to create a plot.
 */ 
#include "fitnessutil.h"
#include <grale/gravitationallens.h>
#include <grale/lensplane.h>
#include <grale/imageplane.h>
#include <grale/constants.h>
#include <grale/randomnumbergenerator.h>
#include <grale/imagesdataextended.h>
#include <grale/imagesbackprojector.h>
#include <grale/compositelens.h>
#include <iostream>
#include <limits>

using namespace std;
using namespace grale;

int main(void)
{
	double zd = 0.4;
	GravitationalLens *pLens = nullptr;
	string errstr;

	if (!GravitationalLens::load("reallens_nosheet.lensdata", &pLens, errstr))
		cerr << errstr << endl;

	double Dd = pLens->getLensDistance();
	LensPlane lp;
	if (!lp.init(pLens, Vector2Dd(-120,-120)*ANGLE_ARCSEC, Vector2Dd(120,120)*ANGLE_ARCSEC, 1024, 1024))
		cerr << lp.getErrorString() << endl;

	RandomNumberGenerator rng;
	double DdsOverDs = rng.pickRandomNumber()*0.4+0.5;

	ImagePlane ip;
	if (!ip.init(&lp, 1.0, DdsOverDs))
		cerr << ip.getErrorString() << endl;

	double betax = (rng.pickRandomNumber()*2.0-1.0)*10*ANGLE_ARCSEC;
	double betay = (rng.pickRandomNumber()*2.0-1.0)*10*ANGLE_ARCSEC;
	Vector2Dd beta { betax, betay };
	vector<Vector2Dd> thetas;

	if (!ip.traceBeta(beta, thetas))
	{
		cerr << "Couldn't find theta positions" << endl;
		return -1;
	}
	if (thetas.size() < 2)
	{
		cerr << "Only found " << thetas.size() << " images" << endl;
		return -1;
	}
	cerr << "Found " << thetas.size() << " images" << endl;

	vector<double> timeDelays;
	double minTd = numeric_limits<double>::max();
	for (int i = 0 ; i < thetas.size() ; i++)
	{
		double td = 0;
		pLens->getTimeDelay(zd, 1.0, DdsOverDs, thetas[i], beta, &td);
		td /= 24*60*60;
		if (td < minTd)
			minTd = td;
		timeDelays.push_back(td);
	}

	for (auto &td : timeDelays)
		td -= minTd;

	ImagesDataExtended imgDat(1.0, DdsOverDs);
	imgDat.create(thetas.size(), false, false);

	for (int i = 0 ; i < thetas.size() ; i++)
	{
		int pt = imgDat.addPoint(i, thetas[i]);
		double td = timeDelays[i];
		cerr << "td = " << td << " days" << endl;
		imgDat.addTimeDelayInfo(i, pt, td);
	}

	ImagesBackProjector bpReal(*pLens, { &imgDat }, zd, false);
	vector<pair<int,int>> refPts;

	cerr << "calculateTimeDelayFitness: " << calculateTimeDelayFitness(bpReal, { 0 }) << std::endl;
	cerr << "calculateTimeDelayFitnessExperimental2: " << calculateTimeDelayFitnessExperimental2(bpReal, { 0 }) << std::endl;
	cerr << "calculateTimeDelayFitness_Relative: " << calculateTimeDelayFitness_Relative(bpReal, { 0 }, refPts) << std::endl;
	cerr << endl;

	cout << "# frac normal experimental experimental2 normal_rel" << endl;
	for (int i = 10 ; i < 200 ; i++)
	{
		double frac = (double)i/100.0;
		cout << frac << " ";

		CompositeLensParams params;
		CompositeLens lens;
		params.addLens(frac, { 0, 0 }, 0, *pLens);

		if (!lens.init(Dd, &params))
			cerr << "Couldn't init composite lens" << endl;

		ImagesBackProjector bp(lens, { &imgDat }, zd, false);
		cout << calculateTimeDelayFitness(bp, { 0 }) << " ";
		cout << calculateTimeDelayFitnessExperimental2(bp, { 0 }) << " ";
		cout << calculateTimeDelayFitness_Relative(bp, { 0 }, refPts) << std::endl;
	}

	return 0;
}
