#include "randomnumbergenerator.h"
#include "imagesdata.h"
#include <errut/booltype.h>
#include <serut/vectorserializer.h>
#include <iostream>
#include <fstream>
#include <set>

using namespace grale;
using namespace errut;
using namespace std;
using namespace serut;

int pickInt(RandomNumberGenerator &rnd, const int m)
{
    return (int)(rnd.pickRandomNumber()*m);
}

bool_t addPointsToImages(RandomNumberGenerator &rnd, ImagesData &imgDat, const vector<ImagesData::PropertyName> &props)
{
    const int numImgs = imgDat.getNumberOfImages();
    vector<pair<ImagesData::PropertyName, double>> propValues;

    for (auto n : props)
        propValues.push_back({n, 123});
    
    const int maxPts = 50;

    for (int i = 0 ; i < numImgs ; i++)
    {
        const int numPts = pickInt(rnd, maxPts);
        for (int p = 0 ; p < numPts ; p++)
        {
            // Fill in new property values
            for (auto &pv : propValues)
                pv.second = rnd.pickRandomNumber();

            const double x = rnd.pickRandomNumber();
            const double y = rnd.pickRandomNumber();
            if (imgDat.addPoint(i, { x, y }, propValues) < 0)
                return "Error adding point: " + imgDat.getErrorString();
        }
    }
    return true;
}

bool_t addPointGroups(RandomNumberGenerator &rnd, ImagesData &imgDat)
{
    const int maxGroups = 10;
    const int maxGroupPts = 20;
    const int numGroups = pickInt(rnd, maxGroups);
    const int numImgs = imgDat.getNumberOfImages();
    if (numImgs == 0)
        return true;

    cerr << "# groups: " << numGroups << endl;

    for (int g = 0 ; g < numGroups ; g++)
    {
        if (imgDat.addGroup() < 0)
            return "Error adding group: " + imgDat.getErrorString();

        const int numGroupPts = pickInt(rnd, maxGroupPts);
        for (int gp = 0 ; gp < numGroupPts ; gp++)
        {
            const int imgIdx = pickInt(rnd, numImgs);
            const int numImgPts = imgDat.getNumberOfImagePoints(imgIdx);
            if (numImgPts < 0)
                return "Error getting number of image points: " + imgDat.getErrorString();
    
            if (numImgPts == 0)
                continue;

            const int ptIdx = pickInt(rnd, numImgPts);
            if (!imgDat.addGroupPoint(g, imgIdx, ptIdx))
                return "Error adding group point: " + imgDat.getErrorString();
        }
    }

    return true;
}

bool_t addTriangulations(RandomNumberGenerator &rnd, ImagesData &imgDat)
{
    const int numImgs = imgDat.getNumberOfImages();
    if (numImgs == 0)
        return true;
    
    const int maxTriang = 100;
    const int numTriang = pickInt(rnd, maxTriang);
    cerr << "# triangles: " << numTriang << endl;

    int count = 0;
    for (int t = 0 ; t < numTriang ; t++)
    {
        const int imgIdx = pickInt(rnd, numImgs);
        const int numPts = imgDat.getNumberOfImagePoints(imgIdx);
        if (numPts < 0)
            return "Can't get number of image points: " + imgDat.getErrorString();
        if (numPts == 0)
            continue;

        const int idx1 = pickInt(rnd, numPts);
        const int idx2 = pickInt(rnd, numPts);
        const int idx3 = pickInt(rnd, numPts);

        set<int> s { idx1, idx2, idx3 };
        if (s.size() != 3)
            continue;

        if (!imgDat.addTriangle(imgIdx, idx1, idx2, idx3))
            return "Can't add triangle: " + imgDat.getErrorString();
        
        count++;
    }
    cerr << "# triangles added: " << count << endl;
    return true;
}

bool_t addTimeDelays(RandomNumberGenerator &rnd, ImagesData &imgDat)
{
    const int numImgs = imgDat.getNumberOfImages();
    if (numImgs == 0)
        return true;

    const int maxTds = 100;
    const int numTds = pickInt(rnd, maxTds);
    cerr << "# tds: " << numTds << endl;

    for (int i = 0 ; i < numTds ; i++)
    {
        const int imgIdx = pickInt(rnd, numImgs);
        const int numPts = imgDat.getNumberOfImagePoints(imgIdx);
        if (numPts < 0)
            return "Error getting number of image points: " + imgDat.getErrorString();
        
        if (numPts == 0)
            continue;
        
        const int ptIdx = pickInt(rnd, numPts);
        const double delay = rnd.pickRandomNumber();

        if (!imgDat.addTimeDelayInfo(imgIdx, ptIdx, delay))
            return "Can't add time delay: " + imgDat.getErrorString();
    }
    return true;
}

vector<ImagesData::PropertyName> pickProperties(RandomNumberGenerator &rnd)
{
    int numProps = pickInt(rnd, (int)ImagesData::MaxProperty);
    set<ImagesData::PropertyName> propSet;
    for (int i = 0 ; i < numProps ; i++)
        propSet.insert((ImagesData::PropertyName)(pickInt(rnd, (int)ImagesData::MaxProperty)));
    return vector<ImagesData::PropertyName>(propSet.begin(), propSet.end());
}

bool_t initImageData(ImagesData &imgDat)
{
    RandomNumberGenerator rnd;
    vector<ImagesData::PropertyName> props = pickProperties(rnd);
    cerr << "# properties: " << props.size() << endl;

    const int maxNumImgs = 10;
    const int numImages = pickInt(rnd, maxNumImgs);

    if (!imgDat.create(numImages, props))
        return "Error creating imagesdata instance: " + imgDat.getErrorString();
    cerr << "# initial images: " << numImages << endl;
    
    auto add = [&rnd, &imgDat, &props]() -> bool_t
    {
        bool_t r;
        if (!(r = addPointsToImages(rnd, imgDat, props)) ||
           (!(r = addPointGroups(rnd, imgDat))) ||
           (!(r = addTriangulations(rnd, imgDat))) ||
           (!(r = addTimeDelays(rnd, imgDat)))
        )
            return r;
        return true;
    };

    bool_t r = add();
    if (!r)
        return r;

    // TODO: add images
    const int numAddImages = pickInt(rnd, maxNumImgs);
    for (int i = 0 ; i < numAddImages ; i++)
    {
        if (imgDat.addImage() < 0)
            return "Error adding image: " + imgDat.getErrorString();
    }
    cerr << "# added images: " << numAddImages << endl;

    r = add();
    if (!(r))
        return r;

    return true;
}

bool_t runTest(const ImagesData &imgDat)
{
    // Check copy constructor
    {    
        ImagesData imgDatCopy(imgDat);
        if (!(imgDatCopy == imgDat))
            return "Error using copy constructor";
    }

    // Check assignment
    {
        ImagesData imgDatCopy;

        imgDatCopy = imgDat;
        if (!(imgDatCopy == imgDat))
            return "Error using assignment operator";
    }

    // Check serialization
    {
        ImagesData imgDatCopy;
        VectorSerializer vSer;

        if (!imgDat.write(vSer))
            return "Error in serialization test: can't write: " + imgDat.getErrorString();
        if (!imgDatCopy.read(vSer))
            return "Error in serialization test: can't read: " + imgDatCopy.getErrorString();
        
        if (!(imgDatCopy == imgDat))
            return "Error using serialization";
        
        cerr << "# serialization bytes: " << vSer.getBuffer().size() << endl;
        // ofstream f("test.dat");
        // f.write((char*)vSer.getBufferPointer(), vSer.getBufferSize());

    }

    return true;
}

int main(int argc, char *argv[])
{
	if (argc == 1)
	{
		const int numTests = 10000;
		for (int i = 0 ; i < numTests ; i++)
		{
			ImagesData imgDat;
			bool_t r = initImageData(imgDat);
			if (!r)
			{
				cerr << "Error: " << r.getErrorString() << endl;
				return -1;
			}

			r = runTest(imgDat);
			if (!r)
			{
				cerr << "Error: " << r.getErrorString() << endl;
				return -1;
			}
		}
		cout << "Ran " << numTests << " tests" << endl;
	}
	else
	{
		for (int i = 1 ; i < argc ; i++)
		{
			ImagesData imgDat;
			if (!imgDat.load(argv[i]))
			{
				cerr << "Can't load " << argv[i] << ": " << imgDat.getErrorString() << endl;
				return -1;
			}

			bool_t r = runTest(imgDat);
			if (!r)
			{
				cerr << "Error: " << r.getErrorString() << endl;
				return -1;
			}
		}
	}
    return 0;
}
