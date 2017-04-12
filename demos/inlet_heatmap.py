
from __future__ import (absolute_import, division, print_function)


from netCDF4 import Dataset
import numpy as np
import pysgrid

import matplotlib.pyplot as plt

import cartopy.crs as ccrs
# from cartopy.io import shapereader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

# rotation is still ugly...
from pysgrid.processing_2d import rotate_vectors, vector_sum

# url = ('http://geoport.whoi.edu/thredds/dodsC/clay/usgs/users/jcwarner/Projects/Sandy/triple_nest/00_dir_NYB05.ncml') # noqa
url2 = ('http://geoport-dev.whoi.edu/thredds/dodsC/clay/usgs/users/zdefne/run076/his/00_dir_roms_display.ncml') # noqa

nc = Dataset(url2)
sgrid = pysgrid.load_grid(nc)
sgrid  # We need a better __repr__ and __str__ !!!
lons, lats = np.mgrid[-74.38:-74.26:600j, 39.45:39.56:600j]

points = np.stack((lons, lats), axis=-1)

print(points.shape)
time_idx = 0
v_idx = 0

# interp_u = sgrid.interpolate_var_to_points(
#     points, sgrid.u, slice=[time_idx, v_idx])
# interp_v = sgrid.interpolate_var_to_points(
#     points, sgrid.v, slice=[time_idx, v_idx])
# sgrid.interpolate_var_to_points(
#     points[19:21, 241:243], sgrid.u, slices=[time_idx, v_idx])
# interp_u = sgrid.interpolate_var_to_points(
#     points, sgrid.u, slices=[time_idx, v_idx])
interp_u = sgrid.interpolate_var_to_points(
    points, sgrid.u[time_idx, v_idx], slices=None)
interp_v = sgrid.interpolate_var_to_points(
    points, sgrid.v, slices=[time_idx, v_idx])

ind = sgrid.locate_faces(points)
ang_ind = ind + [1, 1]
angles = sgrid.angles[:][ang_ind[:, 0], ang_ind[:, 1]]
u_rot, v_rot = rotate_vectors(interp_u, interp_v, angles)
u_rot = u_rot.reshape(600, -1)
v_rot = v_rot.reshape(600, -1)

uv_vector_sum = vector_sum(u_rot, v_rot)


def make_map(projection=ccrs.PlateCarree(), figsize=(20, 20)):
    fig, ax = plt.subplots(figsize=figsize,
                           subplot_kw=dict(projection=projection))
    gl = ax.gridlines(draw_labels=True)
    gl.xlabels_top = gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    return fig, ax


mscale = 1
vscale = 10
scale = 0.03
lon_data = lons
lat_data = lats

fig, ax = make_map()

kw = dict(scale=1.0 / scale, pivot='middle', width=0.003, color='black')
# q = plt.quiver(lon_data[::vscale, ::vscale], lat_data[::vscale, ::vscale],
# u_rot[::vscale, ::vscale], v_rot[::vscale, ::vscale], zorder=2, **kw)

cs = plt.pcolormesh(lon_data[::mscale, ::mscale],
                    lat_data[::mscale, ::mscale],
                    uv_vector_sum[::mscale, ::mscale], zorder=1, cmap=plt.cm.rainbow)
ax.coastlines('10m')
# plt.subplots(figsize=(20,20))
# cs = plt.imshow(uv_vector_sum, cmap=plt.cm.rainbow)
plt.show()
