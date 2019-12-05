from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr
from libcpp cimport bool

cimport grale.vector2d as vector2d
cimport grale.serut as serut
cimport grale.errut as errut

ctypedef vector2d.Vector2Dd Vector2Dd

cdef extern from "grale/plummerlensinfo.h" namespace "grale":
    cdef cppclass PlummerLensInfo:
        PlummerLensInfo(double mass, double angularWidth, Vector2Dd angularPosition)
        double getMass()
        double getAngularWidth()
        Vector2Dd getAngularPosition()

cdef extern from "grale/squarelensinfo.h" namespace "grale":
    cdef cppclass SquareLensInfo:
        SquareLensInfo(double mass, double angularWidth, Vector2Dd angularPosition)
        double getMass()
        double getAngularWidth()
        Vector2Dd getAngularPosition()

cdef extern from "grale/gausslensinfo.h" namespace "grale":
    cdef cppclass GaussLensInfo:
        GaussLensInfo(double mass, double angularWidth, Vector2Dd angularPosition)
        double getMass()
        double getAngularWidth()
        Vector2Dd getAngularPosition()

cdef extern from "grale/gravitationallens.h" namespace "grale::GravitationalLens":
    cdef enum LensType:
        Gaussian,
        MultiplePlummers,
        Plummer,
        Pointmass,
        SIS,
        NSIE,
        NSIS,
        SIE,
        Square,
        MultipleSquares,
        MultipleGaussians,
        MassSheet,
        Composite,
        MassDisk,
        Profile,
        PolynomialMassProfile,
        MultipleWendland,
        DeflectionGrid,
        NFW,
        EllipticNFW,
        Sersic,
        EllipticSersic,
        PIEMD,
        PIMD,
        AlphaPot,
        Harmonic,
        PotentialGrid,
        CircularPieces

cdef extern from "grale/gravitationallens.h" namespace "grale":

    cdef cppclass GravitationalLensParams(errut.ErrorBase):
        pass

    cdef cppclass GravitationalLens(errut.ErrorBase):
        LensType getLensType()
        GravitationalLens *createCopy() 
        bool init(double Dd, GravitationalLensParams *pParams)
        bool traceTheta(double Ds, double Dds, Vector2Dd theta, Vector2Dd *beta)
        bool getAlphaVector(Vector2Dd theta, Vector2Dd *alpha)
        bool getAlphaVectorDerivatives(Vector2Dd theta, double &axx, double &ayy, double &axy)
        double getInverseMagnification(double Ds, double Dds, Vector2Dd theta)
        double getSurfaceMassDensity(Vector2Dd theta)
        double getLensDistance()
        void setLensDistance(double Dd)
        bool getTimeDelay(double zd, double Ds, double Dds, Vector2Dd theta, Vector2Dd beta, double *pTimeDelay)
        bool getProjectedPotential(double Ds, double Dds, Vector2Dd theta, double *pPotentialValue)
        void setDerivativeAngularDistanceScale(double distanceScale)
        @staticmethod
        bool load(string fileName, GravitationalLens **pLens, string &errorString)
        @staticmethod
        bool read(serut.SerializationInterface &si, GravitationalLens **pLens, string &errorString)
        bool save(string fileName)
        bool write(serut.SerializationInterface &si)
        const GravitationalLensParams *getLensParameters() const

cdef extern from "grale/symmetriclens.h" namespace "grale":

    cdef cppclass SymmetricLens(GravitationalLens):
        double getMassInside(double thetaLength) const
        double getProfileSurfaceMassDensity(double thetaLength) const
        @staticmethod
        SymmetricLens *cast(GravitationalLens *pLens)

cdef extern from "grale/gausslens.h" namespace "grale":

    cdef cppclass GaussLens(GravitationalLens): 
        pass
    cdef cppclass GaussLensParams(GravitationalLensParams):
        GaussLensParams(double mass, double angularWidth)
        double getMass()
        double getAngularWidth()

ctypedef const GaussLensParams* GaussLensParamsPtrConst

cdef extern from "grale/multipleplummerlens.h" namespace "grale":

    cdef cppclass MultiplePlummerLens(GravitationalLens): 
        pass
    cdef cppclass MultiplePlummerLensParams(GravitationalLensParams):
        MultiplePlummerLensParams(vector[PlummerLensInfo] &lensInfo)
        const vector[PlummerLensInfo] &getLensInfo()

ctypedef const MultiplePlummerLensParams* MultiplePlummerLensParamsPtrConst

cdef extern from "grale/plummerlens.h" namespace "grale":

    cdef cppclass PlummerLens(GravitationalLens): 
        pass
    cdef cppclass PlummerLensParams(GravitationalLensParams):
        PlummerLensParams(double mass, double angularWidth)
        double getLensMass()
        double getAngularWidth()

ctypedef const PlummerLensParams* PlummerLensParamsPtrConst

cdef extern from "grale/pointmasslens.h" namespace "grale":

    cdef cppclass PointmassLens(GravitationalLens): 
        pass
    cdef cppclass PointmassLensParams(GravitationalLensParams):
        PointmassLensParams(double mass)
        double getLensMass()

ctypedef const PointmassLensParams* PointmassLensParamsPtrConst

cdef extern from "grale/sislens.h" namespace "grale":

    cdef cppclass SISLens(GravitationalLens):
        pass
    cdef cppclass SISLensParams(GravitationalLensParams):
        SISLensParams(double velocityDispersion)
        double getVelocityDispersion()

ctypedef const SISLensParams* SISLensParamsPtrConst

cdef extern from "grale/nsielens.h" namespace "grale":

    cdef cppclass NSIELens(GravitationalLens):
        pass
    cdef cppclass NSIELensParams(GravitationalLensParams):
        NSIELensParams(double velocityDispersion, double ellipticity, double angularCoreRadius)
        double getVelocityDispersion()
        double getEllipticity()
        double getAngularCoreRadius()

ctypedef const NSIELensParams* NSIELensParamsPtrConst

cdef extern from "grale/nsislens.h" namespace "grale":

    cdef cppclass NSISLens(GravitationalLens):
        pass
    cdef cppclass NSISLensParams(GravitationalLensParams):
        NSISLensParams(double velocityDispersion, double angularCoreRadius)
        double getVelocityDispersion()
        double getAngularCoreRadius()

ctypedef const NSISLensParams* NSISLensParamsPtrConst

cdef extern from "grale/sielens.h" namespace "grale":

    cdef cppclass SIELens(GravitationalLens):
        pass
    cdef cppclass SIELensParams(GravitationalLensParams):
        SIELensParams(double velocityDispersion, double ellipticity)
        double getVelocityDispersion()
        double getEllipticity()

ctypedef const SIELensParams* SIELensParamsPtrConst

cdef extern from "grale/squarelens.h" namespace "grale":

    cdef cppclass SquareLens(GravitationalLens):
        pass
    cdef cppclass SquareLensParams(GravitationalLensParams):
        SquareLensParams(double mass, double width)
        double getLensMass()
        double getAngularWidth()

ctypedef const SquareLensParams* SquareLensParamsPtrConst

cdef extern from "grale/multiplesquarelens.h" namespace "grale":

    cdef cppclass MultipleSquareLens(GravitationalLens):
        pass
    cdef cppclass MultipleSquareLensParams(GravitationalLensParams):
        MultipleSquareLensParams(vector[SquareLensInfo] &lensInfo)
        const vector[SquareLensInfo] &getLensInfo()

ctypedef const MultipleSquareLensParams* MultipleSquareLensParamsPtrConst

cdef extern from "grale/multiplegausslens.h" namespace "grale":

    cdef cppclass MultipleGaussLens(GravitationalLens):
        pass
    cdef cppclass MultipleGaussLensParams(GravitationalLensParams):
        MultipleGaussLensParams(vector[GaussLensInfo] &lensInfo)
        const vector[GaussLensInfo] &getLensInfo()

ctypedef const MultipleGaussLensParams* MultipleGaussLensParamsPtrConst

cdef extern from "grale/masssheetlens.h" namespace "grale":

    cdef cppclass MassSheetLens(GravitationalLens):
        pass
    cdef cppclass MassSheetLensParams(GravitationalLensParams):
        MassSheetLensParams(double density)
        MassSheetLensParams(double Dd, double Ds, double Dds)
        double getDensity()

ctypedef const MassSheetLensParams* MassSheetLensParamsPtrConst

cdef extern from "grale/compositelens.h" namespace "grale":

    cdef cppclass CompositeLens(GravitationalLens):
        int getNumberOfSubLenses() const
        const GravitationalLens *getSubLens(int i) const
        Vector2Dd getSubLensPosition(int i) const
        double getSubLensAngleOriginal(int i) const
        double getSubLensFactor(int i) const
    cdef cppclass CompositeLensParams(GravitationalLensParams):
        CompositeLensParams()
        bool addLens(double factor, Vector2Dd position, double angle, const GravitationalLens &lens)

ctypedef CompositeLens *CompositeLensPtr 

cdef extern from "grale/massdisklens.h" namespace "grale":

    cdef cppclass MassDiskLens(GravitationalLens):
        pass
    cdef cppclass MassDiskLensParams(GravitationalLensParams):
        MassDiskLensParams(double density, double angularRadius)
        MassDiskLensParams(double Dd, double Ds, double Dds, double angularRadius)
        double getDensity()
        double getAngularRadius()

ctypedef const MassDiskLensParams* MassDiskLensParamsPtrConst

cdef extern from "grale/polynomialmassprofilelens.h" namespace "grale":

    cdef cppclass PolynomialPart:
        PolynomialPart(double xOffset, double yOffset, double xScale, double yScale, double xEnd, const vector[double] &coefficients)
        double getXOffset()
        double getYOffset()
        double getXScale()
        double getYScale()
        double getEndPosition()
        const vector[double] &getCoefficients()

    cdef cppclass PolynomialMassProfileLens(GravitationalLens):
        pass

    cdef cppclass PolynomialMassProfileLensParams(GravitationalLensParams):
        PolynomialMassProfileLensParams()
        void addPolynomialPart(const PolynomialPart &p)
        const vector[PolynomialPart] &getPolynomialParts()

    cdef cppclass PolynomialZeroMassLensParams(GravitationalLensParams):
        PolynomialZeroMassLensParams(double Dd, double densityScale, double angularRadius, double zeroPoint)

    cdef cppclass PolynomialTimeDelayLensParams(GravitationalLensParams):
        PolynomialTimeDelayLensParams(double theta1, double theta2, double timeDiff, double zLens)

ctypedef const PolynomialMassProfileLensParams* PolynomialMassProfileLensParamsPtrConst

cdef extern from "grale/multiplewendlandlens.h" namespace "grale":

    cdef cppclass WendlandLensInfo:
        WendlandLensInfo()
        WendlandLensInfo(double heightFactor, double angularScale, Vector2Dd angularPosition)
        double getHeightFactor()
        double getAngularScale()
        Vector2Dd getAngularPosition()

    cdef cppclass MultipleWendlandLensParams(GravitationalLensParams):
        MultipleWendlandLensParams()
        MultipleWendlandLensParams(vector[WendlandLensInfo] &phiXInfo, vector[WendlandLensInfo] &phiYInfo)
        const vector[WendlandLensInfo] &getPhiXInfo()
        const vector[WendlandLensInfo] &getPhiYInfo()
        bool matchDeflections(vector[Vector2Dd] &deflectionPoints, vector[Vector2Dd] &deflectionAngles, double angularScale)

    cdef cppclass MultipleWendlandLens(GravitationalLens):
        pass

ctypedef const MultipleWendlandLensParams* MultipleWendlandLensParamsPtrConst

cdef extern from "grale/deflectiongridlens.h" namespace "grale":

    cdef cppclass DeflectionGridLensParams(GravitationalLensParams):
        DeflectionGridLensParams(vector[double] &alphaX, vector[double] &alphaY, int width, int height,
                                 Vector2Dd bottomLeft, Vector2Dd topRight)
        const vector[double] &getAlphaX()
        const vector[double] &getAlphaY()
        int getWidth()
        int getHeight()
        Vector2Dd getBottomLeft()
        Vector2Dd getTopRight()

    cdef cppclass DeflectionGridLens(GravitationalLens):
        pass

ctypedef const DeflectionGridLensParams* DeflectionGridLensParamsPtrConst

cdef extern from "grale/nfwlens.h" namespace "grale":

    cdef cppclass NFWLensParams(GravitationalLensParams):
        NFWLensParams(double rho_s, double theta_s)
        double get3DDensityScale()
        double getAngularRadiusScale()

    cdef cppclass NFWLens(GravitationalLens):
        pass

ctypedef const NFWLensParams* NFWLensParamsPtrConst

cdef extern from "grale/ellipticnfwlens.h" namespace "grale":

    cdef cppclass EllipticNFWLensParams(GravitationalLensParams):
        EllipticNFWLensParams(double rho_s, double theta_s, double q)
        double get3DDensityScale()
        double getAngularRadiusScale()
        double getEllipticity()

    cdef cppclass EllipticNFWLens(GravitationalLens):
        pass

ctypedef const EllipticNFWLensParams* EllipticNFWLensParamsPtrConst

cdef extern from "grale/sersiclens.h" namespace "grale":

    cdef cppclass SersicLensParams(GravitationalLensParams):
        SersicLensParams(double centralDensity, double angularScale, double index)
        double getCentralDensity()
        double getAngularScale()
        double getSersicIndex()

    cdef cppclass SersicLens(GravitationalLens):
        pass

ctypedef const SersicLensParams* SersicLensParamsPtrConst

cdef extern from "grale/ellipticsersiclens.h" namespace "grale":

    cdef cppclass EllipticSersicLensParams(GravitationalLensParams):
        EllipticSersicLensParams(double centralDensity, double angularScale, double index, double q)
        double getCentralDensity()
        double getAngularScale()
        double getSersicIndex()
        double getEllipticity()

    cdef cppclass EllipticSersicLens(GravitationalLens):
        pass

ctypedef const EllipticSersicLensParams* EllipticSersicLensParamsPtrConst

cdef extern from "grale/profilelens.h" namespace "grale":

    cdef cppclass ProfileLensParams(GravitationalLensParams):
        ProfileLensParams(double endRadius, vector[double] &profile)
        double getEndRadius()
        const vector[double] getProfile()

    cdef cppclass ProfileLens(GravitationalLens):
        pass

ctypedef const ProfileLensParams* ProfileLensParamsPtrConst

cdef extern from "grale/piemdlens.h" namespace "grale":

    cdef cppclass PIEMDLensParams(GravitationalLensParams):
        PIEMDLensParams(double sigma0, double coreRadius, double scaleRadius, double epsilon)
        double getCentralDensity()
        double getCoreRadius()
        double getScaleRadius()
        double getEpsilon()

    cdef cppclass PIEMDLens(GravitationalLens):
        pass

ctypedef const PIEMDLensParams* PIEMDLensParamsPtrConst

cdef extern from "grale/pimdlens.h" namespace "grale":

    cdef cppclass PIMDLensParams(GravitationalLensParams):
        PIMDLensParams(double sigma0, double coreRadius, double scaleRadius)
        double getCentralDensity()
        double getCoreRadius()
        double getScaleRadius()

    cdef cppclass PIMDLens(GravitationalLens):
        pass

ctypedef const PIMDLensParams* PIMDLensParamsPtrConst

cdef extern from "grale/alphapotlens.h" namespace "grale":

    cdef cppclass AlphaPotLensParams(GravitationalLensParams):
        AlphaPotLensParams(double b, double s, double q, double K2, double alpha)
        double getB()
        double getS()
        double getQ()
        double getK2()
        double getAlpha()

    cdef cppclass AlphaPotLens(GravitationalLens):
        pass

ctypedef const AlphaPotLensParams* AlphaPotLensParamsPtrConst

cdef extern from "grale/harmoniclens.h" namespace "grale":

    cdef cppclass HarmonicLensParams(GravitationalLensParams):
        HarmonicLensParams(double sigma0, double k, double l, double phiX, double phiY)
        double getDensityScale()
        double getK()
        double getL()
        double getPhiX()
        double getPhiY()

    cdef cppclass HarmonicLens(GravitationalLens):
        pass

ctypedef const HarmonicLensParams* HarmonicLensParamsPtrConst

cdef extern from "grale/potentialgridlens.h" namespace "grale":

    cdef cppclass PotentialGridLensParams(GravitationalLensParams):
        PotentialGridLensParams()

    cdef cppclass PotentialGridLens(GravitationalLens):
        pass

ctypedef const PotentialGridLensParams* PotentialGridLensParamsPtrConst

cdef extern from "grale/circularpieceslens.h" namespace "grale":

    cdef cppclass CircularPieceInfo:
        CircularPieceInfo(const shared_ptr[GravitationalLens] &lens,
			                 double startRadius, double endRadius,
							 double potentialScale, double potentialOffset)

        const shared_ptr[GravitationalLens] &getLens()
        double getStartRadius()
        double getEndRadius()
        double getPotentialScale()
        double getPotentialOffset()

    cdef cppclass CircularPiecesLensParams(GravitationalLensParams):
        CircularPiecesLensParams(const vector[CircularPieceInfo] &pieces, const vector[double] interpolationFunctionCoeffs)
        const vector[CircularPieceInfo] &getPiecesInfo()
        const vector[double] &getInterpolationFunctionCoefficients()

    cdef cppclass CircularPiecesLens(GravitationalLens):
        pass

ctypedef const CircularPiecesLensParams* CircularPiecesLensParamsPtrConst

