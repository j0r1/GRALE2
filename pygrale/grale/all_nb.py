"""Similar to the :mod:`'all' <grale.all>` module, but sets the default
:mod:`feedback <grale.feedback>` mechanism to ``notebook``.
"""
from .all import *

feedback.setDefaultFeedback("notebook")
print("Set feedback style to 'notebook'")
