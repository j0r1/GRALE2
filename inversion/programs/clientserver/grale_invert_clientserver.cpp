#include <enut/ipv4address.h>
#include <errut/booltype.h>
#include "log.h"
#include "inputoutput.h"
#include "gaparameters.h"
#include "lensinversiongafactoryparamssingleplanecpu.h"
#include "inversioncommunicator.h"
#include <serut/memoryserializer.h>
#include <mogal/geneticalgorithm.h>

#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <chrono>

using namespace std;
using namespace grale;
using namespace serut;
using namespace mogal;
using namespace nut;
using namespace errut;

#include "commonga.h"

class CSCommunicator : public InversionCommunicator
{
public:
	CSCommunicator(IPv4Address serverAddress, uint16_t serverPort)
	{
		m_serverAddress = serverAddress;
		m_serverPort = serverPort;
	}
	~CSCommunicator() { }
protected:
	string getVersionInfo() const override { return "Client-server method"; }

	bool_t runGA(int popSize, GAFactory &factory, GeneticAlgorithmParams &params,
	             const std::string &moduleDir, const std::string &moduleFile,
	             const std::vector<uint8_t> &factoryParamBytes) override
	{
		if (!m_ga.run(m_serverAddress, m_serverPort, moduleFile, popSize, factory, &params))
			return "Error running GA: " + m_ga.getErrorString();

		return true;
	}

	const mogal::GeneticAlgorithm *getGeneticAlgorithm() const override { return &m_ga; }
private:
	CommonGA<GeneticAlgorithm> m_ga;
	IPv4Address m_serverAddress;
	uint16_t m_serverPort;
};

int main0(int argc, char *argv[])
{
	if (argc != 3)
	{
		cerr << "ERROR: need a server address and server port" << endl;
		return -1;
	}

	int ip1 = 0, ip2 = 0, ip3 = 0, ip4 = 0;
	if (sscanf(argv[1], "%d.%d.%d.%d", &ip1, &ip2, &ip3, &ip4) != 4)
	{
		cerr << "ERROR: first argument does not seem like an IPv4 address" << endl;
		return -1;
	}

	int portNumber;
	if (!parseAsInt(argv[2], portNumber) || portNumber < 1 || portNumber > 65535)
	{
		cerr << "ERROR: can't interpret second argument as port number" << endl;
		return -1;
	}

	CSCommunicator comm(IPv4Address((uint8_t)ip1,(uint8_t)ip2,(uint8_t)ip3,(uint8_t)ip4), (uint16_t)portNumber);

	bool_t r = comm.run();
	if (!r)
	{
		cerr << "ERROR: " << r.getErrorString() << endl;
		return -1;
	}

	return 0;
}

int main(int argc, char *argv[])
{
	LOG.init(argv[0]);
#ifdef _WIN32
	WSADATA dat;
	WSAStartup(MAKEWORD(2,2),&dat);
#endif // _WIN32
	int r = main0(argc, argv);
#ifdef _WIN32
	WSACleanup();
#endif // _WIN32
	cerr << "MAIN: REACHED END" << endl;
	exit(r);
	return r;
}
