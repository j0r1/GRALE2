"""TODO"""

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
    """TODO"""
    def __init__(self, z_d, rangeInfos = {}, cosmology = "default"):
        """TODO"""
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
        """TODO"""
        self.liCache = { }

    def saveCaches(self):
        """TODO"""
        for n in self.liCache:
            # TODO: use a feedbackobject
            print("Saving", n, "as", self.liCache[n]["filename"])
            cfn = self._getCacheFileName(n)
            pickle.dump(self.liCache[n], open(cfn,"wb"))

    def getLensInfoEntry(self, fileNameOrLens):       
        """TODO"""
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
        """TODO"""
        d = self.getLensInfoEntry(fileNameOrLens)
        return d["li"][rangeName]
