import grale.renderers as renderers
import grale.feedback as feedback

# Allowed string values are 'threads' and 'mpi'
renderers.setDefaultMassRenderer("threads")

# Allowed string values are 'threads', 'mpi' and 'opencl'
# (but not all lens models have an OpenCL implementation available)
renderers.setDefaultLensPlaneRenderer("threads")

# Allowed string values are 'none', 'stdout' and 'notebook'
feedback.setDefaultFeedback("notebook")

