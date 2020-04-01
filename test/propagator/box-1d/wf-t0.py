import numpy as np
from numpy import exp, sin, pi

from qprop.param.file import ParameterFile
param = ParameterFile("in.param")
Nx, dx = (param[key] for key in ("Nx", "dx"))

Nx_tot = 1 + Nx + 1
xarr = 0. + dx * np.arange(Nx_tot)


# Superposition of eigenstates
L = xarr[-1] - xarr[0]
wf_E1 = sin(pi/L*xarr).astype(np.complex)
wf_E2 = sin(2*pi/L*xarr).astype(np.complex)
wf_E3 = sin(3*pi/L*xarr).astype(np.complex)
wf_t0 = wf_E1 + 0.3 * wf_E2 - 0.22 * wf_E3
#wf_t0 = wf_E1

# Gaussian wave packet
xmin, xmax = xarr[[0,-1]]
xmid = 0.5 * (xmin + xmax)
kx = 0.
#wf_t0 = exp(-(xarr-xmid)**2).astype(np.complex)
#wf_t0 *= exp(1.j*kx*xarr)


wf_t0[[0,-1]] = 0.
from tdse.propagator.box1d import Wavefunction_Uniform_1D_Box 
wf_obj = Wavefunction_Uniform_1D_Box(Nx, dx, x0=0.0)
wf_obj.normalize(wf_t0, dx)

wf_t0.tofile("wf-t0.bin")
