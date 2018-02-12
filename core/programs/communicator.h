#ifndef COMMUNICATOR_H

#define COMMUNICATOR_H

#include <errut/booltype.h>
#include <stdint.h>
#include <string>
#include <vector>

namespace grale 
{
	class GravitationalLens;
}

class Communicator
{
public:
	typedef errut::bool_t bool_t;

	Communicator();
	virtual ~Communicator();

	bool_t render();
protected:
	virtual std::string getRenderType() const = 0;
	virtual std::string getVersionInfo() const = 0;
	
	void setStatus(const std::string &s);
	void setProgress(double current, double target);
private:
	virtual bool_t renderGrid(const std::vector<uint8_t> &lensData, grale::GravitationalLens *pLens, 
	                          double x0, double y0, double dX, double dY, int numX, int numY, 
							  std::vector<double> &renderPoints) = 0;
	virtual bool_t renderPointVector(const std::vector<uint8_t> &lensData, grale::GravitationalLens *pLens, 
	                          const std::vector<double> &inputXY, std::vector<double> &renderPoints) = 0;
};

#endif // COMMUNICATOR_H
