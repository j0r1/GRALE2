#include "randomnumbergenerator.h"
#include <iostream>

using namespace grale;
using namespace std;

int main(void)
{
	RandomNumberGenerator rng;

	cout << "This rng should work in both optimized and debug mode" << endl;
	cout << rng.getRandomDouble() << endl;
#ifndef NDEBUG
	cout << "Debug mode, assertion should be triggered!" << endl;
#else
	cout << "Optimized mode, no assertion should be triggered for next random number" << endl; // Even if use is not the way it's intended
#endif

	thread t([&rng](){
			cout << rng.getRandomDouble() << endl;
	});

	t.join();

	return 0;
}
