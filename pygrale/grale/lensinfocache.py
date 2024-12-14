"""This module provides the :class:`LensInfoCache <grale.lensinfocache.LensInfoCache>`
class, which can be used to cache the calculations performed in a 
:class:`LensInfo <grale.plotutil.LensInfo>` object. This can be done either
just in memory, or written to disk as well. In the latter case, if at a
later time the cache does not exist in memory anymore but can be loaded from disk,
this can save time recalculating e.g. the mass map."""

import pickle
import hashlib
import os
from . import plotutil
from . import privutil
from . import lenses

class LensInfoCacheException(Exception):
    """An exception that's generated if something goes wrong in a function
    provided by this module."""
    pass

class LensInfoCache(object):
    """The main class to keep track of :class:`LensInfo <grale.plotutil.LensInfo>` 
    objects, for different lens models (at same redshift), but for certain
    sizes. This can come in handy when several inversion results need to
    be compared."""

    def __init__(self, z_d, rangeInfos = {}, cosmology = "default"):
        """Initialize the instance, where a lens redshift of `z_d` will
        be assumed for all stored :class:`LensInfo <grale.plotutil.LensInfo>`
        objects. Together with `cosmology`, this will be passed on to
        the constructor of each :class:`LensInfo <grale.plotutil.LensInfo>`
        instance.

        Further information about what to calculate is stored in
        `rangeInfos`. This is a dictionary with keys that can be chosen
        freely, the values are passed on to a :class:`LensInfo <grale.plotutil.LensInfo>`
        constructor to specify for example the size of the region.
        This way, for example, you could prepare the cache to hold `LensInfo` 
        instances for a smaller, strong lensing region, and a wider, weak
        lensing region:

        .. code-block:: python

            { "strong": { "size": 120*ANGLE_ARCSEC },
              "weak": { "size": 400*ANGLE_ARCSEC, "zs": 2.5 } }

        When e.g. :func:`getLensInfo(someModel, "strong") <grale.lensinfocache.LensInfoCache.getLensInfo>`
        is called for the ``"strong"`` parameters, either a new :class:`LensInfo <grale.plotutil.LensInfo>`
        is constructed, or a cached one (from memory or from disk) is returned,
        depending on whether or not the lens model `someModel` 
        (specified as a :class:`lens object <grale.lenses.GravitationalLens>` 
        or a file name) has been seen before.
        """

        self.liCache = { }
        self.cacheFileNameFormat = "licache_{}_.pickle"
        self.zd = z_d
        self.rangeInfos = rangeInfos
        
        cosmology = privutil.initCosmology(cosmology)
        if not cosmology:
            raise LensInfoCacheException("No cosmological model was specified")
        
        self.cosm = cosmology

    def _getCacheFileName(self, hsh):
        return self.cacheFileNameFormat.format(hsh)

    def clear(self):
        """Clears the caches that are stored in memory, but does **not** remove
        anything stored on disk."""
        self.liCache = { }

    def saveCaches(self):
        """Save all caches that are currently in memory to disk."""
        for n in self.liCache:
            # TODO: use a feedbackobject
            cfn = self._getCacheFileName(n)
            print("Saving", self.liCache[n]["filename"], "as", cfn)
            pickle.dump(self.liCache[n], open(cfn,"wb"))

    def getLensInfoEntry(self, fileNameOrLens):
        """Returns the complete entry for the lens model in `fileNameOrLens`,
        which can either be the name of a file, or a :class:`lens model <grale.lenses.GravitationalLens>`
        in memory. This contains the information about all currently known
        :class:`LensInfo <grale.plotutil.LensInfo>` objects, with sizes as
        specified in the constructor, for the lens model. If it doesn't exist
        in the in-memory or on-disk cache yet, new :class:`LensInfo <grale.plotutil.LensInfo>` 
        objects are created.
        """
        if type(fileNameOrLens) == str:
            data = open(fileNameOrLens, "rb").read()
        else:
            data = fileNameOrLens.toBytes()
            
        s = hashlib.sha256()
        s.update(data)
        h = s.hexdigest()
        if h in self.liCache:
            return self.liCache[h]

        cfn = self._getCacheFileName(h)

        if os.path.exists(cfn):
            self.liCache[h] = pickle.load(open(cfn, "rb"))
            assert h == self.liCache[h]["hash"], "Hash of loaded data doesn't match expected value"
            return self.liCache[h]

        if type(fileNameOrLens) == str:
            lens = lenses.GravitationalLens.load(fileNameOrLens)
            fileName = fileNameOrLens
        else:
            lens = fileNameOrLens
            fileName = "(lens model)"

        d = { "li": {}, "filename": fileName, "hash": h, "lens": lens}
        for rangeName in self.rangeInfos:
            kwargs = self.rangeInfos[rangeName]
            d["li"][rangeName] = plotutil.LensInfo(lens, **kwargs, zd=self.zd, cosmology=self.cosm)
                    
        self.liCache[h] = d
        return d
    
    def getLensInfo(self, fileNameOrLens, rangeName):
        """Returns the :class:`LensInfo <grale.plotutil.LensInfo>` object for
        the model specified by `fileNameOrLens`, for the size specified by
        `rangeName`. This last parameter should have been specified in the
        constructor."""
        d = self.getLensInfoEntry(fileNameOrLens)
        return d["li"][rangeName]
