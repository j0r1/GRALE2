#include "inversioncommunicator.h"
#include "inputoutput.h"
#include "galensmodule.h"
#include "lensfitnessobject.h"
#include <errut/booltype.h>
#include <memory>
#include <iostream>
#include <string>

using namespace std;
using namespace grale;
using namespace errut;

class UsageCommunicator : public InversionCommunicator
{
public:
	UsageCommunicator() { }
	~UsageCommunicator() { }
protected:
	string getVersionInfo() const override { return "Inversion module usage information"; }

	bool_t runModule(const string &moduleDir, const string &moduleFile, GALensModule *pModule)
	{
		unique_ptr<LensFitnessObject> f(pModule->createFitnessObject());
		if (!f.get())
			return "Unable to create lens fitness object for module: " + moduleFile;

		string usage = f->getUsage();
		vector<char> txt(usage.length());
		for (int i = 0 ; i < txt.size() ; i++)
			txt[i] = usage[i];

		WriteLineStdout(strprintf("USAGE:%d", (int)txt.size()));
		WriteBytesStdout(&txt[0], txt.size());

		return true;
	}
};

int main(int argc, char *argv[])
{
	UsageCommunicator comm;

	bool_t r = comm.run();
	if (!r)
	{
		cerr << "ERROR: " << r.getErrorString() << endl;
		return -1;
	}

	return 0;
}
