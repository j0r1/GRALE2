#include "imagesdata.h"
#include <iostream>

using namespace std;
using namespace grale;

void dump(const string &fn)
{
    ImagesData imgDat;

    if (!imgDat.load(fn))
        cout << "ERROR: can't load " << fn << ": " << imgDat.getErrorString() << endl;
    cout << "FN: " << fn << endl;
    cout << "  Number of images: " << imgDat.getNumberOfImages() << endl;
    for (int i = 0 ; i < imgDat.getNumberOfImages() ; i++)
        cout << "  Image " << i << ": " << imgDat.getNumberOfImagePoints(i) << " points" << endl;
}

int main(int argc, char const *argv[])
{
    for (int i = 1 ; i < argc ; i++)
        dump(argv[i]);
    return 0;
}
