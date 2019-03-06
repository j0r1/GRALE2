/*

  This file is a part of GRALE, a library to facilitate the simulation
  and inversion of gravitational lenses.

  Copyright (C) 2008-2012 Jori Liesenborgs

  Contact: jori.liesenborgs@gmail.com
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
  
*/

#include "graleconfig.h"
#include "gravitationallens.h"
#include "pointmasslens.h"
#include "sislens.h"
#include "gausslens.h"
#include "plummerlens.h"
#include "multipleplummerlens.h"
#include "nsielens.h"
#include "nsislens.h"
#include "sielens.h"
#include "squarelens.h"
#include "multiplesquarelens.h"
#include "multiplegausslens.h"
#include "masssheetlens.h"
#include "massdisklens.h"
#include "compositelens.h"
#include "profilelens.h"
#include "polynomialmassprofilelens.h"
#include "multiplewendlandlens.h"
#include "deflectiongridlens.h"
#include "nfwlens.h"
#include "ellipticnfwlens.h"
#include "sersiclens.h"
#include "ellipticsersiclens.h"
#include "piemdlens.h"
#include "pimdlens.h"
#include "alphapotlens.h"
#include "constants.h"
#include <serut/fileserializer.h>
#include <serut/memoryserializer.h>
#include <serut/dummyserializer.h>
#include <sys/types.h>
#include <stdint.h>
#include <iostream>

#include "debugnew.h"

using namespace std;

namespace grale
{

GravitationalLens::GravitationalLens(GravitationalLens::LensType t)
{
	m_init = false;
	m_lensType = t; 
	m_pParameters = 0;
	m_Dd = 0;
}

GravitationalLens::~GravitationalLens()
{
	if (m_pParameters != 0)
		delete m_pParameters;
}

bool GravitationalLens::init(double D_d, const GravitationalLensParams *pLensParams)
{
	if (m_init)
	{
		setErrorString("Already initialized");
		return false;
	}
	
	m_Dd = D_d;
	
	if (pLensParams != 0)
	{
		m_pParameters = pLensParams->createCopy();
		if (m_pParameters == 0)
		{
			setErrorString("Can't create copy of the lens parameters");
			return false;
		}
	}
	else
		m_pParameters = 0;
	
	if (!processParameters(pLensParams))
	{
		if (m_pParameters != 0)
		{
			delete m_pParameters;
			m_pParameters = 0;
		}
		return false;
	}

	m_derivDistScale = -1.0;
	m_init = true;
	return true;
}

bool GravitationalLens::traceTheta(double D_s, double D_ds, Vector2D<double> theta, Vector2D<double> *pBeta) const
{
	if (!m_init)
	{
		setErrorString("Lens is not initialized");
		return false;
	}
	
	Vector2D<double> alphaVector;
	
	if (!getAlphaVector(theta, &alphaVector))
		return false;
	
	*pBeta = theta - (D_ds/D_s)*alphaVector;
	return true;
}

// If no explicit code is provided, a numerical approximation will be used.
bool GravitationalLens::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
	Vector2D<double> alpha, alphadx, alphady;
	double dtheta = m_derivDistScale;
	
	if (dtheta < 0)
	{
		setErrorString("Using derivative numerical approximation code but typical scale was not set!");
		return false;
	}
	
	if (!getAlphaVector(theta, &alpha))
		return false;
	
	if (!getAlphaVector(theta+Vector2D<double>(dtheta,0.0), &alphadx))
		alphadx = alpha;
	
	if (!getAlphaVector(theta+Vector2D<double>(0.0,dtheta), &alphady))
		alphady = alpha;
	
	axx = (alphadx.getX()-alpha.getX())/dtheta;
	ayy = (alphady.getY()-alpha.getY())/dtheta;

	// We'll use the average of two calculation which should yield the same value
	axy = ((alphady.getX()-alpha.getX())/dtheta + (alphadx.getY()-alpha.getY())/dtheta)/2.0;
	
	return true;
}

double GravitationalLens::getInverseMagnification(double D_s, double D_ds, Vector2D<double> theta) const
{
	double factor = D_ds/D_s;
	double a11, a22, a12;
	
	if (!getAlphaVectorDerivatives(theta, a11, a22, a12))
	{
		std::cerr << "ERROR: " << getErrorString() << std::endl;
		return 1.0; // to avoid a crash
	}

	a11 *= factor;
	a22 *= factor;
	a12 *= factor;
	
	return ((1.0-a11)*(1.0-a22)-a12*a12);
}

#define GRAVLENSIDENTIFIER2						0x57415247
#define LENSNUMBER_POINTMASS						0
#define LENSNUMBER_SIS							1
#define LENSNUMBER_GAUSSIAN						2
#define LENSNUMBER_PLUMMER						3
#define LENSNUMBER_MULTIPLEPLUMMERS					4
#define LENSNUMBER_NSIE							5
#define LENSNUMBER_NSIS							6
#define LENSNUMBER_SIE							7
#define LENSNUMBER_SQUARE						8
#define LENSNUMBER_MULTIPLESQUARES					9
#define LENSNUMBER_MULTIPLEGAUSSIANS					10
#define LENSNUMBER_MASSSHEET						11
#define LENSNUMBER_COMPOSITE						12
#define LENSNUMBER_MASSDISK						13
#define LENSNUMBER_PROFILE						14
#define LENSNUMBER_POLYNOMIALMASSPROFILE				15
#define LENSNUMBER_MULTIPLEWENDLANDLENS					16
#define LENSNUMBER_DEFLECTIONGRIDLENS					17
#define LENSNUMBER_NFW							18
#define LENSNUMBER_ELLNFW						19
#define LENSNUMBER_SERSIC						20
#define LENSNUMBER_ELLSERSIC						21
#define LENSNUMBER_PIEMD							22
#define LENSNUMBER_PIMD								23
#define LENSNUMBER_ALPHAPOT								24
//#define LENSNUMBER_TESTLENS						234	

bool GravitationalLens::write(serut::SerializationInterface &si) const
{
	if (!m_init)
	{
		setErrorString("Lens was not initialized");
		return false;
	}
	
	int32_t id, lensNumber, gotParams;

	id = GRAVLENSIDENTIFIER2;
	switch(m_lensType)
	{
	case Pointmass:
		lensNumber = LENSNUMBER_POINTMASS;
		break;	
	case SIS:
		lensNumber = LENSNUMBER_SIS;
		break;
	case Gaussian:
		lensNumber = LENSNUMBER_GAUSSIAN;
		break;
	case Plummer:
		lensNumber = LENSNUMBER_PLUMMER;
		break;
	case MultiplePlummers:
		lensNumber = LENSNUMBER_MULTIPLEPLUMMERS;
		break;
	case NSIE:
		lensNumber = LENSNUMBER_NSIE;
		break;
	case NSIS:
		lensNumber = LENSNUMBER_NSIS;
		break;
	case SIE:
		lensNumber = LENSNUMBER_SIE;
		break;
	case Square:
		lensNumber = LENSNUMBER_SQUARE;
		break;
	case MultipleSquares:
		lensNumber = LENSNUMBER_MULTIPLESQUARES;
		break;
	case MultipleGaussians:
		lensNumber = LENSNUMBER_MULTIPLEGAUSSIANS;
		break;
	case MassSheet:
		lensNumber = LENSNUMBER_MASSSHEET;
		break;
	case Composite:
		lensNumber = LENSNUMBER_COMPOSITE;
		break;
	case MassDisk:
		lensNumber = LENSNUMBER_MASSDISK;
		break;
	case Profile:
		lensNumber = LENSNUMBER_PROFILE;
		break;
	case PolynomialMassProfile:
		lensNumber = LENSNUMBER_POLYNOMIALMASSPROFILE;
		break;
	case MultipleWendland:
		lensNumber = LENSNUMBER_MULTIPLEWENDLANDLENS;
		break;
	case DeflectionGrid:
		lensNumber = LENSNUMBER_DEFLECTIONGRIDLENS;
		break;
	case NFW:
		lensNumber = LENSNUMBER_NFW;
		break;
	case EllipticNFW:
		lensNumber = LENSNUMBER_ELLNFW;
		break;
	case Sersic:
		lensNumber = LENSNUMBER_SERSIC;
		break;
	case EllipticSersic:
		lensNumber = LENSNUMBER_ELLSERSIC;
		break;
	case PIEMD:
		lensNumber = LENSNUMBER_PIEMD;
		break;
	case PIMD:
		lensNumber = LENSNUMBER_PIMD;
		break;
	case AlphaPot:
		lensNumber = LENSNUMBER_ALPHAPOT;
		break;
	default:
		setErrorString("Lens type not recognized");
		return false;
	}

	if (m_pParameters == 0)
		gotParams = 0;
	else
		gotParams = 1;

	if (!si.writeInt32(id))
	{
		setErrorString("Error writing ID");
		return false;
	}
	if (!si.writeInt32(lensNumber))
	{
		setErrorString("Error writing lens type");
		return false;
	}
	if (!si.writeDouble(m_Dd))
	{
		setErrorString("Error writing lens distance");
		return false;
	}
	if (!si.writeDouble(m_derivDistScale))
	{
		setErrorString("Error writing derivative distance");
		return false;
	}		
	if (!si.writeInt32(gotParams))
	{
		setErrorString("Error writing lens parameters");
		return false;
	}
	if (m_pParameters)
	{
		if (!m_pParameters->write(si))
		{
			setErrorString("Error writing lens parameters");
			return false;
		}
	}
	return true;
}

bool GravitationalLens::read(serut::SerializationInterface &si, GravitationalLens **pLens,std::string &errorString)
{
	GravitationalLens *pTmpLens;
	GravitationalLensParams *pParams = 0;
	int32_t id, gotParams, lensNumber;
	double Dd, derivDistScale;
	
	if (!si.readInt32(&id))
	{
		errorString = std::string("Error reading ID");
		return false;
	}
	if (id != GRAVLENSIDENTIFIER2)
	{
		errorString = std::string("An incorrect ID was read");
		return false;
	}

	if (!si.readInt32(&lensNumber))
	{
		errorString = std::string("Error reading lens type");
		return false;
	}
	if (!si.readDouble(&Dd))
	{
		errorString = std::string("Error reading lens distance");
		return false;
	}

	if (!si.readDouble(&derivDistScale))
	{
		errorString = std::string("Error reading derivative distance");
		return false;
	}

	if (!si.readInt32(&gotParams))
	{
		errorString = std::string("Error reading lens parameters");
		return false;
	}
	
	switch (lensNumber)
	{
	case LENSNUMBER_POINTMASS:
		pTmpLens = new PointmassLens();
		pParams = new PointmassLensParams();
		break;
	case LENSNUMBER_SIS:
		pTmpLens = new SISLens();
		pParams = new SISLensParams();
		break;
	case LENSNUMBER_GAUSSIAN:
		pTmpLens = new GaussLens();
		pParams = new GaussLensParams();
		break;
	case LENSNUMBER_PLUMMER:
		pTmpLens = new PlummerLens();
		pParams = new PlummerLensParams();
		break;
	case LENSNUMBER_MULTIPLEPLUMMERS:
		pTmpLens = new MultiplePlummerLens();
		pParams = new MultiplePlummerLensParams();
		break;
	case LENSNUMBER_NSIE:
		pTmpLens = new NSIELens();
		pParams = new NSIELensParams();
		break;
	case LENSNUMBER_NSIS:
		pTmpLens = new NSISLens();
		pParams = new NSISLensParams();
		break;
	case LENSNUMBER_SIE:
		pTmpLens = new SIELens();
		pParams = new SIELensParams();
		break;
	case LENSNUMBER_SQUARE:
		pTmpLens = new SquareLens();
		pParams = new SquareLensParams();
		break;
	case LENSNUMBER_MULTIPLESQUARES:
		pTmpLens = new MultipleSquareLens();
		pParams = new MultipleSquareLensParams();
		break;
	case LENSNUMBER_MULTIPLEGAUSSIANS:
		pTmpLens = new MultipleGaussLens();
		pParams = new MultipleGaussLensParams();
		break;
	case LENSNUMBER_MASSSHEET:
		pTmpLens = new MassSheetLens();
		pParams = new MassSheetLensParams();
		break;
	case LENSNUMBER_COMPOSITE:
		pTmpLens = new CompositeLens();
		pParams = new CompositeLensParams();
		break;
	case LENSNUMBER_MASSDISK:
		pTmpLens = new MassDiskLens();
		pParams = new MassDiskLensParams();
		break;
	case LENSNUMBER_PROFILE:
		pTmpLens = new ProfileLens();
		pParams = new ProfileLensParams();
		break;
	case LENSNUMBER_POLYNOMIALMASSPROFILE:
		pTmpLens = new PolynomialMassProfileLens();
		pParams = new PolynomialMassProfileLensParams();
		break;
	case LENSNUMBER_MULTIPLEWENDLANDLENS:
		pTmpLens = new MultipleWendlandLens();
		pParams = new MultipleWendlandLensParams();
		break;
	case LENSNUMBER_DEFLECTIONGRIDLENS:
		pTmpLens = new DeflectionGridLens();
		pParams = new DeflectionGridLensParams();
		break;
	case LENSNUMBER_NFW:
		pTmpLens = new NFWLens();
		pParams = new NFWLensParams();
		break;
	case LENSNUMBER_ELLNFW:
		pTmpLens = new EllipticNFWLens();
		pParams = new EllipticNFWLensParams();
		break;
	case LENSNUMBER_SERSIC:
		pTmpLens = new SersicLens();
		pParams = new SersicLensParams();
		break;
	case LENSNUMBER_ELLSERSIC:
		pTmpLens = new EllipticSersicLens();
		pParams = new EllipticSersicLensParams();
		break;
	case LENSNUMBER_PIEMD:
		pTmpLens = new PIEMDLens();
		pParams = new PIEMDLensParams();
		break;
	case LENSNUMBER_PIMD:
		pTmpLens = new PIMDLens();
		pParams = new PIMDLensParams();
		break;
	case LENSNUMBER_ALPHAPOT:
		pTmpLens = new AlphaPotLens();
		pParams = new AlphaPotLensParams();
		break;
	default:
		errorString = std::string("Can't recognize lens type");
		return false;
	}

	bool initStatus = true;
	
	if (gotParams)
	{
		if (!pParams->read(si))
		{
			delete pTmpLens;
			delete pParams;
			errorString = std::string("Error reading lens parameters");
			return false;
		}
		initStatus = pTmpLens->init(Dd, pParams);
	}
	else
		initStatus = pTmpLens->init(Dd, 0);
	
	if (pParams)
		delete pParams;

	if (!initStatus)
	{
		errorString = std::string(std::string("Couldn't initialize the new lens: ") + pTmpLens->getErrorString());
		delete pTmpLens;
		return false;
	}
	//cerr << "derivDistScale = " << derivDistScale << endl;
	pTmpLens->m_derivDistScale = derivDistScale;
	*pLens = pTmpLens;
	return true;
}

bool GravitationalLens::load(const std::string &fileName, GravitationalLens **pLens, std::string &errorString)
{
	serut::FileSerializer fs;

	if (!fs.open(fileName, serut::FileSerializer::ReadOnly))
	{
		errorString = std::string(std::string("Error opening lens file: ") + fs.getErrorString());
		return false;
	}
	
	return GravitationalLens::read(fs, pLens, errorString);
}

bool GravitationalLens::save(const std::string &fileName) const
{
	serut::FileSerializer fs;

	if (!fs.open(fileName, serut::FileSerializer::WriteOnly))
	{
		setErrorString(std::string("Error opening lens file: ") + fs.getErrorString());
		return false;
	}
	return write(fs);
}
	
bool GravitationalLens::getTimeDelay(double z_d, double D_s, double D_ds, Vector2D<double> theta, 
		                     Vector2D<double> beta, double *pTimeVal) const
{
	double potentialValue = 0;

	if (!getProjectedPotential(D_s, D_ds, theta, &potentialValue))
		return false;

	Vector2D<double> diff = theta - beta;
	*pTimeVal = ((m_Dd*(1.0+z_d))/SPEED_C)*(D_s/D_ds)*(0.5*diff.getLengthSquared()-potentialValue);
	
	return true;
}

bool GravitationalLens::getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const
{
	setErrorString("No projected potential is implemented for this lens");
	return false;
}

GravitationalLens *GravitationalLens::createCopy() const
{
	serut::DummySerializer dumSer;

	if (!write(dumSer))
	{
		setErrorString(std::string("Couldn't count the number of bytes to allocate: ") + getErrorString());
		return 0;
	}

	size_t length = dumSer.getBytesWritten();
	std::vector<uint8_t> buffer(length);
	serut::MemorySerializer memSer(&(buffer[0]), length, &(buffer[0]), length);

	if (!write(memSer))
	{
		setErrorString(std::string("Couldn't write the lens data into a new memory buffer: ") + getErrorString());
		return 0;
	}

	GravitationalLens *pLens = 0;
	std::string errStr;

	if (!GravitationalLens::read(memSer, &pLens, errStr))
	{
		setErrorString(std::string("Couldn't recreate the lens from the data in the memory buffer: ") + errStr);
		return 0;
	}

	return pLens;
}

bool GravitationalLens::getSuggestedScales(double *pDeflectionScale, double *pPotentialScale) const
{
	setErrorString("Scale suggestion not implemented for this lens");
	return false;
}

bool GravitationalLens::getCLParameterCounts(int *pNumIntParams, int *pNumFloatParams) const
{
	setErrorString("No OpenCL implementation is provided");
	return false;
}

bool GravitationalLens::getCLParameters(double deflectionScale, double potentialScale, int *pIntParams, float *pFloatParams) const
{
	setErrorString("No OpenCL implementation is provided");
	return false;
}

std::string GravitationalLens::getCLProgram(std::string &subRoutineName) const
{
	subRoutineName = "NOT_IMPLEMENTED";

	return std::string("");
}

std::string GravitationalLens::getCLLensQuantitiesStructure() const
{
	std::string str;

	str += "typedef struct\n";
	str += "{\n";
	str += " 	float alphaX;\n";
	str += " 	float alphaY;\n";
	str += "	float potential;\n";
	str += " 	float axx;\n";
	str += "	float ayy;\n";
	str += " 	float axy;\n";
	str += "} LensQuantities;\n";

	return str;
}

std::string GravitationalLens::getCLLensProgram(std::string &subRoutineName) const
{
	std::string prog = getCLLensQuantitiesStructure();
	
	prog += getCLProgram(subRoutineName);

	return prog;
}

} // end namespace

