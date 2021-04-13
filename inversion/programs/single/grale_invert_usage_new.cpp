#include "inversioncommunicatornewga.h"
#include "inputoutput.h"
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

	bool_t runModule(const std::string &lensFitnessObjectType, 
	                 std::unique_ptr<grale::LensFitnessObject> fitnessObject,
					 const std::string &calculatorType) override
    {
		auto &f = fitnessObject;
		if (!f.get())
			return "Unable to create lens fitness object for module: " + lensFitnessObjectType;

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
