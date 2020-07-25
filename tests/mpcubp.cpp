#include "vector2d.h"
#include "plummerlensinfo.h"
#include "cosmology.h"
#include "mpcudabackprojector.h"
#include "imagesdata.h"
#include <vector>
#include <iostream>
#include <sstream>

using namespace std;
using namespace grale;
#include "mpcubp_generated.h"

bool approxSame(double x1, double x2)
{
    if (std::abs(x1-x2)/std::abs(x1) > 1e-5)
        return false;
    if (std::abs(x1-x2)/std::abs(x2) > 1e-5)
        return false;
    return true;        
}

bool approxSame(Vector2Dd v1, Vector2Dd v2)
{
    return approxSame(v1.getX(), v2.getX()) && approxSame(v1.getY(), v2.getY());
}

string S(Vector2Dd v, double unit = ANGLE_ARCSEC)
{
    v /= unit;
    stringstream ss;
    ss << "(" << v.getX() << ", " << v.getY() << ")";
    return ss.str();
}

string S(Vector2Df v, double unit = ANGLE_ARCSEC)
{
    return S(Vector2Dd(v.getX(), v.getY()), unit);
}

Vector2Dd V(Vector2Df x)
{
    return Vector2Dd(x.getX(), x.getY());
}

shared_ptr<ImagesData> toImage(const vector<pair<Vector2Dd, Vector2Dd>> &thetasAndBetas)
{
    auto imgDat = make_shared<ImagesData>();
    for (int img = 0 ; img < thetasAndBetas.size() ; img++)
    {
        imgDat->addImage();
        imgDat->addPoint(img, thetasAndBetas[img].first);
        // cout << "Added " << S(thetasAndBetas[img].first) << endl;
    }

    return imgDat;
}

int main(int argc, char const *argv[])
{
    vector<float> redshifts;
    vector<vector<PlummerLensInfo>> plummers;

    for (auto &zAndPlummers : redshiftsAndPlummers)
    {
        redshifts.push_back((float)zAndPlummers.first);
        plummers.push_back(zAndPlummers.second);
        
        cout << "Lensplane redshift: " << zAndPlummers.first << endl;
        cout << "  #plummers = " << zAndPlummers.second.size() << endl;
    }
    // cout << "massFactors.size() == " << massFactors.size() << endl;
    // for (auto &m : massFactors)
    //     cout << "  lenses: " << m.size() << endl;
    
    vector<float> sourceRedshifts;
    vector<shared_ptr<ImagesData>> images;
    
    for (auto &zsAndThetaBeta : thetaBetaMappings)
    {
        sourceRedshifts.push_back(zsAndThetaBeta.first);
        images.push_back(toImage(zsAndThetaBeta.second));
    }

    MPCUDABackProjector bp;
    string libKey = "MPCUDALIBRARY";
    if (!getenv(libKey.c_str()))
    {
        cerr << "Please set " << libKey << " to the full path of the library" << endl;
        return -1;
    }
    string libPath(getenv(libKey.c_str()));

    if (!bp.init(libPath, cosm, redshifts, plummers, sourceRedshifts, images))
    {
        cerr << "Unable to initialize MPCUDABackProjector: " << bp.getErrorString() << endl;
        return -1;
    }

    if (!bp.calculateSourcePositions(massFactors, sheetValues))
    {
        cerr << "Error calculating source positions: " << bp.getErrorString() << endl;
        return -1;
    }

    double angUnit = bp.getAngularScale();
    cout << "Angular unit is " << angUnit/ANGLE_ARCSEC << " arcsec" << endl;
    
    for (int s = 0 ; s < bp.getNumberOfSources() ; s++)
    {
        cout << "Source: " << s << " z = " << sourceRedshifts[s]
             << ", #images = " << bp.getNumberOfImages(s) << endl;
        bool ok = true;

        for (int i = 0 ; i < bp.getNumberOfImages(s) ; i++)
        {
            int numPoints = bp.getNumberOfImagePoints(s, i);
            assert(numPoints == 1);

            const Vector2Df *pThetas = bp.getThetas(s, i);
            const Vector2Df *pBetas = bp.getBetas(s, i);
            
            auto t = V(pThetas[0]), b = V(pBetas[0]);
            auto bReal = thetaBetaMappings[s].second[i].second;

            t *= angUnit;
            b *= angUnit;

            if (!approxSame(b, bReal))            
            {
                cout << "  ! " << S(t) << " -> "
                         << S(b) << " != "
                         << S(bReal) << endl;
                ok = false;
            }
        }
        if (ok)
            cout << "  all ok" << endl;
    }

    return 0;
}
