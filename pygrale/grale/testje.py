from  grale.inversion import *

data = { "string": "123", "int": 123, "float": 123.0, "boolean": True }
c = ConfigurationParameters(data)
d = c.asDict()
print(d)

b = open("dbg_gafactparams.dat", "rb").read()
p = GridLensInversionParameters.fromBytes(b)

b2 = p.toBytes()
open("tmp.dat", "wb").write(b2)
