import numpy as np

# Load data from files
debug_dxp_arr = np.fromfile("./debug-dxp-arr.bin", dtype=np.float)
debug_impl_dxp_arr = np.fromfile("./debug-impl-eq-dxp-arr.bin", dtype=np.float)
debug_wf_dxp_arr = np.fromfile("./debug-wf-dxp-arr.bin", dtype=np.complex)
debug_velo_dxp_arr = np.fromfile("./debug-velocity-dxp-arr.bin", dtype=np.float)

# Load parameters
from qprop.param.file import ParameterFile
param = ParameterFile("in.param")
dt = param["dt"]
dx = param["dx"]

dx_grid = (debug_dxp_arr[-1] - debug_dxp_arr[0]) / 3
assert abs(dx-dx_grid) < 1e-11
x_grid_arr = debug_dxp_arr[0] + dx_grid * np.arange(4)

import matplotlib.pyplot as plt
fig, ax = plt.subplots()
l_impl_eq, = ax.plot(debug_dxp_arr, debug_impl_dxp_arr, '.-k')
ax.plot(ax.get_xlim(), [0,0], color='gray')
ax.plot([0,0], ax.get_ylim(), color='gray')

ax.vlines(x_grid_arr, *ax.get_ylim())

axwf = ax.twinx()
l_wf_dxp_arr, = axwf.plot(debug_dxp_arr, debug_wf_dxp_arr.real)
axwf.plot(debug_dxp_arr, np.abs(debug_wf_dxp_arr))

#axv = ax.twinx()
#l_velo, = axv.plot(debug_dxp_arr, debug_velo_dxp_arr, linestyle='--')
v_dt_arr = dt * debug_velo_dxp_arr
l_v_dt, = ax.plot(debug_dxp_arr, v_dt_arr, linestyle='--')

ax.legend(
        [l_impl_eq, l_wf_dxp_arr, l_v_dt], 
        ["impl-eq", "wf-dxp-arr", "velocity * dt"]
        )

plt.show()
fig.savefig("debug-impl-dxp-arr-test-1.png")

