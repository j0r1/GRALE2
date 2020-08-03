#include "multiplanecuda.h"
#include "constants.h"
#include <iostream>

using namespace grale;
using namespace std;

int main(int argc, char *argv[])
{
	if (argc != 2)
	{
		cerr << "Specify path to library" << endl;
		return -1;
	}
	MultiPlaneCUDA mpcu;
	vector<float> sourceRedshifts = { 1.25f, 2.0f };

	if (!mpcu.init(argv[1], 0, // device id
	            ANGLE_ARCSEC,
				0.7f, 0.3f, 0.f, 0.7f, -1.f, // Cosmological model
				{ 0.5f, 1.5f }, // lens redshifts
				{ { { { 0.f, 0.f }, 5.f, 1e14*MASS_SOLAR } },
				  { { { 1.0f, -1.0f}, 3.0f, 5e13*MASS_SOLAR } } }, // plummer parameters
				sourceRedshifts, // source redshifts
				{ { { 10.f, 8.f } },
				  { { -12.0f, 2.0f } } } // theta's
				))
	{
		cerr << "Unable to initialize MultiPlaneCUDA: " << mpcu.getErrorString() << endl;
		return -1;
	}

	cout << "Initialized" << endl;
	
	vector<float> factors { 0.5, 1, 2.0 };
	for (auto f : factors)
	{
		cout << "Using factor " << f << endl;
		if (!mpcu.calculateSourcePositions({ { f }, { f } }))
		{
			cerr << "Can't calculate source positions: " << mpcu.getErrorString() << endl;
			return -1;
		}
		for (size_t sIdx = 0 ; sIdx < sourceRedshifts.size() ; sIdx++)
		{
			cout << "Positions for source " << sIdx << endl;
			auto pPositions = mpcu.getSourcePositions(sIdx);
			if (!pPositions)
			{
				cerr << "Can't get source positions: " << mpcu.getErrorString() << endl;
				return -1;
			}

			auto positions = *pPositions;
			for (auto p : positions)
				cout << " " << p.getX() << " " << p.getY() << endl;
		}
	}

	return 0;
}
