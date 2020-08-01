#include "multiplanecontainer.h"
#include <limits>

using namespace std;

namespace grale
{

MultiPlaneContainerParams::MultiPlaneContainerParams()
{
}

MultiPlaneContainerParams::~MultiPlaneContainerParams()
{
}

bool MultiPlaneContainerParams::add(shared_ptr<GravitationalLens> lens, double z)
{
    if (z < 0)
    {
        setErrorString("Redshift is negative");
        return false;
    }
    if (!lens.get())
    {
        setErrorString("No lens specified");
        return false;
    }

    m_lensesAndRedshifts.push_back({lens, z});
    return true;
}

bool MultiPlaneContainerParams::add(GravitationalLens *pLens, double z)
{
    if (z < 0)
    {
        setErrorString("Redshift is negative");
        return false;
    }
    if (!pLens)
    {
        setErrorString("No lens specified");
        return false;
    }

    shared_ptr<GravitationalLens> copy(pLens->createCopy());
    if (!copy.get())
    {
        setErrorString("Couldn't create copy of lens");
        return false;
    }

    m_lensesAndRedshifts.push_back({ copy, z });
    return true;
}

GravitationalLensParams *MultiPlaneContainerParams::createCopy() const
{
    MultiPlaneContainerParams *pCopy = new MultiPlaneContainerParams();
    for (auto lensZ : m_lensesAndRedshifts)
    {
        // The .get() will cause a copy to be made
        if (!pCopy->add(lensZ.first.get(), lensZ.second))
        {
            setErrorString("Error copying sublens: " + pCopy->getErrorString());
            delete pCopy;
            return nullptr;
        }
    }
    return pCopy;
}

bool MultiPlaneContainerParams::write(serut::SerializationInterface &si) const
{
    int32_t num = (int32_t)m_lensesAndRedshifts.size();
    if (!si.writeInt32(num))
    {
        setErrorString("Couldn't write number of contained lenses: " + si.getErrorString());
        return false;
    }

    for (auto lensZ : m_lensesAndRedshifts)
    {
        auto lens = lensZ.first;
        double z = lensZ.second;
        if (!lens->write(si))
        {
            setErrorString("Couldn't write a contained lens: " + lens->getErrorString());
            return false;
        }

        if (!si.writeDouble(z))
        {
            setErrorString("Couldn't write a contained lens's redshift: " + si.getErrorString());
            return false;
        }
    }
    return true;
}

bool MultiPlaneContainerParams::read(serut::SerializationInterface &si)
{
    int32_t num;
    if (!si.readInt32(&num))
    {
        setErrorString("Couldn't read number of contained lenses: " + si.getErrorString());
        return false;
    }

    m_lensesAndRedshifts.clear();
    for (int32_t i = 0 ; i < num ; i++)
    {
        GravitationalLens *pLens = nullptr;
        string errStr;

        if (!GravitationalLens::read(si, &pLens, errStr))
        {
            setErrorString("Unable to read a contained lens: " + errStr);
            return false;
        }

        shared_ptr<GravitationalLens> lens(pLens);
        double z;

        if (!si.readDouble(&z))
        {
            setErrorString("Unable to read contained lens's readshift: " + si.getErrorString());
            return false;
        }

        m_lensesAndRedshifts.push_back({ lens, z });
    }    
    return true;
}

MultiPlaneContainer::MultiPlaneContainer() : GravitationalLens(MPContainer)
{

}

MultiPlaneContainer::~MultiPlaneContainer()
{

}

string MultiPlaneContainer::s_onlyContainerError = "This 'lens' is only a container for multiple lenses at different redshifts";

bool MultiPlaneContainer::processParameters(const GravitationalLensParams *pLensParams)
{
    if (getLensDistance() != 0)
    {
        setErrorString("For this container, the lens distance must be set to 0");
        return false;
    }

    const MultiPlaneContainerParams *pParams = dynamic_cast<const MultiPlaneContainerParams *>(pLensParams);
    if (!pParams)
    {
        setErrorString("No lens parameters specified, or parameters of wrong type");
        return false;
    }

    int numLenses = pParams->getNumberOfLenses();
    if (numLenses == 0)
    {
        setErrorString("No lenses specified in parameters");
        return false;
    }

    for (int i = 0 ; i < numLenses ; i++)
    {
        double z = pParams->getRedshift(i);
        if (z < 0)
        {
            setErrorString("Invalid contained lens redshift (" + to_string(z) + ")");
            return false;
        }
    }

    return true;
}

bool MultiPlaneContainer::getAlphaVector(Vector2D<double> theta, Vector2D<double> *pAlpha) const
{
    *pAlpha = Vector2Dd { numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN() };
    setErrorString(s_onlyContainerError);
    return false;
}

double MultiPlaneContainer::getSurfaceMassDensity(Vector2D<double> theta) const
{
    setErrorString(s_onlyContainerError);
    return numeric_limits<double>::quiet_NaN();
}

bool MultiPlaneContainer::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
    axx = numeric_limits<double>::quiet_NaN();
    ayy = numeric_limits<double>::quiet_NaN();
    axy = numeric_limits<double>::quiet_NaN();
    setErrorString(s_onlyContainerError);
    return false;
}

bool MultiPlaneContainer::getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const
{
    *pPotentialValue = numeric_limits<double>::quiet_NaN();
    setErrorString(s_onlyContainerError);
    return false;
}

} // end namespace