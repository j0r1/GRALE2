import grale.cosmology as cosmology
from grale.constants import DIST_MPC
import numpy as np

for i in range(20):
    h = np.random.uniform(0.6, 0.8)
    Wm = np.random.uniform(0.25, 0.35)
    Wr = np.random.uniform(0, 0.1)
    if np.random.random() < 0.5: # flat
        Wv = 1.0-Wr-Wm
    else:
        Wv = np.random.uniform(0.6, 0.8)
    
    w = np.random.uniform(-1.1, -0.9)
    
    cosm = cosmology.Cosmology(h, Wm, Wr, Wv, w)
    
    z1 = 0 if np.random.random() < 0.5 else np.random.uniform(0, 1.5)
    z2 = np.random.uniform(0.1, 4.5)
    resultmpc = cosm.getAngularDiameterDistance(z1, z2) / DIST_MPC
    print(f'{{ {{ "h", {h} }}, ', end='')
    print(f'{{ "Wm", {Wm} }}, ', end='')
    print(f'{{ "Wr", {Wr} }}, ', end='')
    print(f'{{ "Wv", {Wv} }}, ', end='')
    print(f'{{ "w", {w} }}, ', end='')
    print(f'{{ "z1", {z1} }}, ', end='')
    print(f'{{ "z2", {z2} }}, ', end='')
    print(f'{{ "resultmpc", {resultmpc} }} }},')
