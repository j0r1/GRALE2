import grale.inversionparams as ip
import pprint
import numpy as np

tests = [
    { 'a': 1,
      'b': 2.0,
      'c': True,
      'd': 'Hello',
      'e': None
    },
    { 'A': np.array([2,3,5,7,11,13], dtype=np.int32),
      'B': np.array([1.0,2.0,-3.0]),
      'C': np.array([False,True,True,False]),
      'D': ['Hello', 'Hi'],
      'E': []
    },
]

for test in tests:
    d = ip.ConfigurationParameters(test)
    pprint.pprint(d.asDict())

