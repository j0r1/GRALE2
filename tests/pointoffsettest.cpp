#include "imagesbackprojector.h"
#include "nsielens.h"
#include "positionrandomizationbpwrapper.h"
#include "randomnumbergenerator.h"
#include "imagesdataextended.h"
#include <iostream>

using namespace std;
using namespace grale;
using namespace errut;

int main(void)
{
	shared_ptr<eatk::RandomNumberGenerator> rng = make_shared<RandomNumberGenerator>();
	size_t numSrcs = rng->getRandomUint32() % 10 + 5;
	vector<bool> haveRndOffsets;
	size_t rndCount = 0;
	for (size_t i = 0 ; i < numSrcs ; i++)
	{
		if (rng->getRandomDouble() > 0.5)
		{
			haveRndOffsets.push_back(true);
			rndCount++;
		}
		else
			haveRndOffsets.push_back(false);
	}

	cout << "Using " << rndCount << " sources of " << numSrcs << " with random offsets" << endl;
	
	vector<shared_ptr<ImagesDataExtended>> baseImages;
	vector<shared_ptr<ImagesDataExtended>> imagesWithOffsets;
	vector<Vector2Dd> offsetsDouble;

	const double LargeSize = 40;
	const double SmallSize = 0.1;

	for (size_t s = 0 ; s < numSrcs ; s++)
	{
		size_t numImgs = rng->getRandomUint32() % 5 + 2;
		double distFrac = rng->getRandomDouble(0.7, 0.9);
		shared_ptr<ImagesDataExtended> imgDat = make_shared<ImagesDataExtended>(1.0, distFrac);
		shared_ptr<ImagesDataExtended> imgDatWithOff = make_shared<ImagesDataExtended>(1.0, distFrac);

		baseImages.push_back(imgDat);
		imagesWithOffsets.push_back(imgDatWithOff);
		
		if (!imgDat->create(numImgs, {}))
			throw runtime_error("Can't create img: " + imgDat->getErrorString());

		if (!imgDatWithOff->create(numImgs, {}))
			throw runtime_error("Can't create img: " + imgDat->getErrorString());

		for (size_t i = 0 ; i < numImgs ; i++)
		{
			size_t numPts = rng->getRandomUint32() % 5 + 3;
			for (size_t p = 0 ; p < numPts ; p++)
			{
				double xPos = rng->getRandomDouble(-LargeSize, LargeSize)*ANGLE_ARCSEC;
				double yPos = rng->getRandomDouble(-LargeSize, LargeSize)*ANGLE_ARCSEC;

				imgDat->addPoint(i, { xPos, yPos });
				if (!haveRndOffsets[s])
					imgDatWithOff->addPoint(i, { xPos, yPos });
				else
				{
					double xOff = rng->getRandomDouble(-SmallSize, SmallSize)*ANGLE_ARCSEC;
					double yOff = rng->getRandomDouble(-SmallSize, SmallSize)*ANGLE_ARCSEC;

					imgDatWithOff->addPoint(i, { xPos+xOff, yPos+yOff });
					offsetsDouble.push_back({ xOff, yOff });
					//offsetsDouble.push_back({ 0, 0 });
				}
			}
		}
	}

	cout << "Number of small offsets: " << offsetsDouble.size() << endl;
	
	shared_ptr<GravitationalLens> lens = make_shared<NSIELens>();
	NSIELensParams params(1000000, 0.75, 1*ANGLE_ARCSEC);
	if (!lens->init(1000*DIST_MPC, &params))
		throw runtime_error("Can't init lens: " + lens->getErrorString());
	double zd = 0.4;

	auto toList = [](auto &vecofSharedPtr)
	{
		list<ImagesDataExtended*> l;
		for (auto &x : vecofSharedPtr)
			l.push_back(x.get());
		return l;
	};

	shared_ptr<ImagesBackProjector> bpBase = make_shared<ImagesBackProjector>(lens, toList(baseImages), zd);
	shared_ptr<ImagesBackProjector> bpWithOffsets = make_shared<ImagesBackProjector>(lens, toList(imagesWithOffsets), zd);

	//cerr << bpBase.getAngularScale() << " " << bpWithOffsets.getAngularScale() << endl;

	shared_ptr<PositionRandomizationBackprojectWrapper> rndBp = make_shared<PositionRandomizationBackprojectWrapper>(bpBase);
	bool_t r;

	size_t numExpectedOffsets;
	if (!(r = rndBp->initRandomization(haveRndOffsets, numExpectedOffsets)))
		throw runtime_error(r.getErrorString());

	cout << "Expecting " << numExpectedOffsets << " offsets" << endl;
	cout << "Have " << offsetsDouble.size() << " offsets" << endl;

	vector<Vector2Df> offsetsFloat;
	for (auto x : offsetsDouble)
	{
		x /= bpBase->getAngularScale();
		offsetsFloat.push_back({ (float)x.getX(), (float)x.getY()});
	}

	if (!(r = rndBp->setRandomOffsets(offsetsFloat)))
		throw runtime_error(r.getErrorString());

	auto getRescaledVector = [](auto v, auto &bp)
	{
		double scale = bp->getAngularScale();
		Vector2Dd x(v.getX(), v.getY());
		x *= scale;
		x /= ANGLE_ARCSEC;
		return x;
	};

	auto getRescaledPhi = [](auto v, auto &bp)
	{
		double scale = bp->getAngularScale();
		scale *= scale;
		double x = v;
		x *= scale;
		x /= ANGLE_ARCSEC*ANGLE_ARCSEC;
		return x;
	};
	
	double maxAlpha = 0, maxBeta = 0, maxTheta = 0, maxPsi = 0;

	auto checkMaxValue = [](auto &maxVar, auto value1, auto value2)
	{
		auto diff = std::abs(value1-value2);
		if (diff > maxVar)
			maxVar = diff;
	};

	auto checkMaxVector = [&checkMaxValue](auto &maxVar, auto value1, auto value2)
	{
		checkMaxValue(maxVar, value1.getX(), value2.getX());
		checkMaxValue(maxVar, value1.getY(), value2.getY());
	};

	for (size_t s = 0 ; s < rndBp->getNumberOfSources() ; s++)
	{
		if (haveRndOffsets[s])
			cout << "OFFSETS!" << endl;
		else
			cout << "No offsets, should be same" << endl;

		for (size_t i = 0 ; i < rndBp->getNumberOfImages(s) ; i++)
		{
			for (size_t p = 0 ; p < rndBp->getNumberOfImagePoints(s, i) ; p++)
			{
				auto theta1 = getRescaledVector(rndBp->getThetas(s, i)[p], rndBp);
				auto theta2 = getRescaledVector(bpWithOffsets->getThetas(s, i)[p], bpWithOffsets);
				auto alpha1 = getRescaledVector(rndBp->getAlphas(s, i)[p], rndBp);
				auto alpha2 = getRescaledVector(bpWithOffsets->getAlphas(s, i)[p], bpWithOffsets);
				auto beta1 = getRescaledVector(rndBp->getBetas(s, i)[p], rndBp);
				auto beta2 = getRescaledVector(bpWithOffsets->getBetas(s, i)[p], bpWithOffsets);
				auto phi1 = getRescaledPhi(rndBp->getLensPotential(s, i)[p], rndBp);
				auto phi2 = getRescaledPhi(bpWithOffsets->getLensPotential(s, i)[p], bpWithOffsets);

				cout << "thetaX " << theta1.getX() << " " << theta2.getX() << " " << std::abs(theta1.getX() - theta2.getX()) << endl;
				cout << "thetaY " << theta1.getY() << " " << theta2.getY() << " " << std::abs(theta1.getY() - theta2.getY()) << endl;
				cout << "alphaX " << alpha1.getX() << " " << alpha2.getX() << " " << std::abs(alpha1.getX() - alpha2.getX()) << endl;
				cout << "alphaY " << alpha1.getY() << " " << alpha2.getY() << " " << std::abs(alpha1.getY() - alpha2.getY()) << endl;
				cout << "betaX " << beta1.getX() << " " << beta2.getX() << " " << std::abs(beta1.getX() - beta2.getX()) << endl;
				cout << "betaY " << beta1.getY() << " " << beta2.getY() << " " << std::abs(beta1.getY() - beta2.getY()) << endl;
				cout << "phi " << phi1 << " " << phi2 << " " << std::abs(phi1 - phi2) << endl;

				checkMaxVector(maxTheta, theta1, theta2);
				checkMaxVector(maxAlpha, alpha1, alpha2);
				checkMaxVector(maxBeta, beta1, beta2);
				checkMaxValue(maxPsi, phi1, phi2);
			}
		}
	}

	cout << endl;
	cout << "maxTheta: " << maxTheta << endl;
	cout << "maxAlpha: " << maxAlpha << endl;
	cout << "maxBeta: " << maxBeta << endl;
	cout << "maxPsi: " << maxPsi << endl;
	return 0;
}
