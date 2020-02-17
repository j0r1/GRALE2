#include "fitnessutil.h"
#include <grale/gravitationallens.h>
#include <grale/lensplane.h>
#include <grale/imageplane.h>
#include <grale/constants.h>
#include <grale/randomnumbergenerator.h>
#include <grale/imagesdataextended.h>
#include <grale/imagesbackprojector.h>
#include <grale/compositelens.h>
#include <grale/masssheetlens.h>
#include <iostream>
#include <limits>
#include <memory>

using namespace std;
using namespace grale;

int main(void)
{
	double zd = 0.4;
	GravitationalLens *pLoadedLens = nullptr;
	string errstr;

	if (!GravitationalLens::load("reallens_nosheet.lensdata", &pLoadedLens, errstr))
		cerr << errstr << endl;
	
	unique_ptr<GravitationalLens> loadedLens(pLoadedLens);
	double Dd = pLoadedLens->getLensDistance();

	MassSheetLens sheetLens;
	MassSheetLensParams sheetLensParams{ 1.0 };
	if (!sheetLens.init(Dd, &sheetLensParams))
		cerr << "Can't init sheet lens: " << sheetLens.getErrorString() << endl;

	CompositeLens realLens;
	CompositeLensParams realLensParams;

	double loadedLensFrac = 1.0;
	double sheetLensFrac = 3.0;
	realLensParams.addLens(loadedLensFrac, {0, 0}, 0, *pLoadedLens);
	realLensParams.addLens(sheetLensFrac, {0, 0}, 0, sheetLens);

	if (!realLens.init(Dd, &realLensParams))
		cerr << "Can't init real lens: " << realLens.getErrorString() << endl;

	LensPlane lp;
	if (!lp.init(&realLens, Vector2Dd(-120,-120)*ANGLE_ARCSEC, Vector2Dd(120,120)*ANGLE_ARCSEC, 1024, 1024))
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
		realLens.getTimeDelay(zd, 1.0, DdsOverDs, thetas[i], beta, &td);
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

	ImagesBackProjector bpReal(realLens, { &imgDat }, zd, false);
	vector<pair<int,int>> refPoints;

	cerr << "calculateTimeDelayFitness: " << calculateTimeDelayFitness(bpReal, { 0 }) << std::endl;
	cerr << "calculateTimeDelayFitnessExperimental2: " << calculateTimeDelayFitnessExperimental2(bpReal, { 0 }) << std::endl;
	cerr << endl;

	cout << "# frac normal experimental experimental2" << endl;
	for (int i = 10 ; i < 200 ; i++)
	{
		double frac = (double)i/100.0;
		cout << frac << " ";

		CompositeLensParams params;
		CompositeLens lens;
		params.addLens(loadedLensFrac, { 0, 0 }, 0, *loadedLens.get());
		params.addLens(frac*sheetLensFrac, { 0, 0 }, 0, sheetLens);

		if (!lens.init(Dd, &params))
			cerr << "Couldn't init composite lens" << endl;

		ImagesBackProjector bp(lens, { &imgDat }, zd, false);
		cout << calculateTimeDelayFitness(bp, { 0 }) << " ";
		cout << calculateTimeDelayFitnessExperimental2(bp, { 0 }) << std::endl;
	}

	return 0;
}
