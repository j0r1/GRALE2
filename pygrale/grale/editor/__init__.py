import os
import sys

editorDir = os.path.join(os.path.dirname(__file__))
if not editorDir in sys.path:
    sys.path.append(editorDir)
