#include "pernodecounter.h"
#include <iostream>

using namespace std;

#ifndef WIN32 // The PerNodeCounter class only works on non-windows anyway
#include <unistd.h>

using namespace grale;

int main(int argc, char const *argv[])
{
	if (argc != 3)
	{
		cerr << "Specify file name and delay (seconds)" << endl;
		return -1;
	}

	PerNodeCounter counter(argv[1]);
	int delay = atoi(argv[2]);
	if (delay < 0 || delay > 1000)
	{
		cerr << "Invalid delay " << delay << endl;
		return -1;
	}

	int num = counter.getCount();
	if (num < 0)
	{
		cerr << "Error: " << counter.getErrorString() << endl;
		return -1;
	}
	cout << "Count: " << num << endl;
	cout << "Sleeping " << delay << " seconds" << endl;
	sleep(delay);
	return 0;
}

#else

int main(void)
{
	cerr << "Not available on windows" << endl;
	return -1;
}

#endif // !WIN32