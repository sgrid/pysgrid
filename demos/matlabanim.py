#!/usr/bin/env python

"""
An animated image
"""

from __future__ import (absolute_import, division, print_function)

import netCDF4 as nc4
import numpy as np
import pysgrid
from datetime import timedelta
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FFMpegWriter
from pysgrid.processing_2d import rotate_vectors, vector_sum

import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

"""
~~HOW TO USE THIS~~

The four lines below are the primary 'controls'.
`url`: filename or URL of dataset
`lons,lats` mesh of points that cover the specified bounds
`maxslice` The number of slices from 0 to turn into an animation
`fps` The number of frames in between each data slice.

This is a WIP. Reader beware.
"""

# url = ('C:\Users\Jay.Hennen\Documents\Code\pygnome\py_gnome\scripts\script_curv_field\TBOFS.nc')
# lons, lats = np.mgrid[-82.8:-82.5:600j, 27.5:27.75:600j]
url = ('http://geoport-dev.whoi.edu/thredds/dodsC/clay/usgs/users/zdefne/run076/his/00_dir_roms_display.ncml') # noqa
lons, lats = np.mgrid[-74.38:-74.26:600j, 39.45:39.56:600j]
# lons, lats = np.mgrid[-82.8:-82.5:600j, 27.5:27.75:600j]
maxslice = 3
fps = 10


def interpolated_velocities(grid, points, ind, timeobj, tindex, u, v, depth=-1.):
    '''
    Finds velocities at the points at the time specified, interpolating in 2D
    over the u and v grids to do so.
    :param time: The time in the simulation
    :param points: a numpy array of points that you want to find interpolated velocities for
    :param indices: Numpy array of indices of the points, if already known.
    :return: interpolated velocities at the specified points
    '''
    t_alphas = timeobj.ialphas(tindex)
    t_index = int(np.floor(tindex))
    mem = True
    _hash = grid._hash_of_pts(points)

    u0 = grid.interpolate_var_to_points(points,
                                        grid.u,
                                        slices=[t_index, -1],
                                        memo=mem,
                                        _hash=_hash)
    u1 = grid.interpolate_var_to_points(points,
                                        grid.u,
                                        slices=[t_index + 1, -1],
                                        memo=mem, _hash=_hash)

    v0 = grid.interpolate_var_to_points(points,
                                        grid.v,
                                        slices=[t_index, -1],
                                        memo=mem,
                                        _hash=_hash)
    v1 = grid.interpolate_var_to_points(points,
                                        grid.v,
                                        slices=[t_index + 1, -1],
                                        memo=mem,
                                        _hash=_hash)

    u_vels = u0 + (u1 - u0) * t_alphas
    v_vels = v0 + (v1 - v0) * t_alphas

    return u_vels, v_vels


class Time(object):

    def __init__(self, data, base_dt_str=None):
        """

        :param data: A netCDF, biggus, or dask source for time data
        :return:
        """
        self.time = nc4.num2date(data[:], units=data.units)

    def ialphas(self, index):
        '''
        given a floating point index between 0 and max index,
        give interpolation alphas for that time
        '''
        i0 = np.floor(index)
        i1 = np.ceil(index)
        frac = index - i0
        return frac
        t0 = self.time[i0]
        t1 = self.time[i1]
        if i0 == i1:
            return t0
        else:
            return t0 * frac + t1 * (1 - frac)

    def time_str(self, index):
        i0 = np.floor(index)
        i1 = np.ceil(index)
        frac = index - i0
        t0 = self.time[i0]
        t1 = self.time[i1]
        time = t0 + timedelta(seconds=(t1 - t0).total_seconds() * frac)
        return time.strftime('%c')


def f(time):
    '''
    time: float index
    '''
    vels = interpolated_velocities(sgrid,
                                   points,
                                   timeobj,
                                   time,
                                   sgrid.u,
                                   sgrid.v,
                                   # u_alphas,
                                   # v_alphas,
                                   # u_ind,
                                   # v_ind,
                                   )

    u_rot = vels[:, 0]
    v_rot = vels[:, 1]
    u_rot, v_rot = rotate_vectors(u_rot, v_rot, angles)
    u_rot = u_rot.reshape(600, -1)
    v_rot = v_rot.reshape(600, -1)

    uv_vector_sum = vector_sum(u_rot, v_rot)
    return uv_vector_sum


def make_map(projection=ccrs.PlateCarree(), figsize=(9, 9)):
    fig, ax = plt.subplots(figsize=figsize,
                           subplot_kw=dict(projection=projection))
    gl = ax.gridlines(draw_labels=True)
    gl.xlabels_top = gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    return fig, ax


nc = nc4.Dataset(url)
timeobj = sgrid = None
if ('ocean_time' in nc.variables.keys()):
    timeobj = Time(nc['ocean_time'])
else:
    timeobj = Time(nc['time'])


if 'grid' in nc.variables.keys():
    sgrid = pysgrid.load_grid(nc)
else:
    sgrid = pysgrid.SGrid(node_lon=nc['lon_psi'],
                          node_lat=nc['lat_psi'],
                          edge1_lon=nc['lon_u'],
                          edge1_lat=nc['lat_u'],
                          edge2_lon=nc['lon_v'],
                          edge2_lat=nc['lat_v'],
                          )
    sgrid.u = pysgrid.variables.SGridVariable(data=nc['u'])
    sgrid.v = pysgrid.variables.SGridVariable(data=nc['v'])
    sgrid.angles = pysgrid.variables.SGridVariable(data=nc['angle'])


points = np.stack((lons, lats), axis=-1).reshape(-1, 2)

ind = sgrid.locate_faces(points)

ang_ind = ind + [1, 1]
angles = sgrid.angles[:][ang_ind[:, 0], ang_ind[:, 1]]

# ims is a list of lists, each row is a list of artists to draw in the
# current frame; here we are just animating one artist, the image, in
# each frame
fig, ax = make_map()
print(fig)
print(ax)
index = 0
ax.coastlines('10m')

t = np.linspace(0, maxslice, maxslice * fps)
cs = qv = tl = None
time_str = timeobj.time_str(0)
tl = ax.text(0, 1, time_str, bbox=dict(
    facecolor='white', alpha=0.8), transform=ax.transAxes)


def gen_map(k):
    global t, index, cs, qv, tl, timeobj
    tindex = t[index]
    if cs is not None:
        cs.remove()
        qv.remove()
    time_str = timeobj.time_str(tindex)
    tl.set_text(time_str)
    mscale = 1
    vscale = 15
    scale = 0.04
    lon_data = lons
    lat_data = lats

    print(tindex)
    print(time_str)
    u_rot, v_rot = interpolated_velocities(sgrid, points, ind, timeobj, tindex, sgrid.u, sgrid.v)
    u_rot, v_rot = rotate_vectors(u_rot, v_rot, angles)
    u_rot = u_rot.reshape(600, -1)
    v_rot = v_rot.reshape(600, -1)

    uv_vector_sum = vector_sum(u_rot, v_rot)

    kw = dict(scale=1.0 / scale, pivot='middle', width=0.003, color='black')
    cs = plt.pcolormesh(lon_data[::mscale, ::mscale],
                        lat_data[::mscale, ::mscale],
                        uv_vector_sum[::mscale, ::mscale], zorder=1, cmap=plt.cm.rainbow)
    qv = plt.quiver(lon_data[::vscale, ::vscale], lat_data[::vscale, ::vscale],
                    u_rot[::vscale, ::vscale], v_rot[::vscale, ::vscale], zorder=2, **kw)
    index += 1
    return cs, qv, tl


print('creating animation')
ani = animation.FuncAnimation(fig,
                              gen_map,
                              frames=maxslice * fps - 1,
                              interval=100,
                              blit=True,
                              repeat=False)

writer = FFMpegWriter(fps=fps, bitrate=1500)
# plt.show()

print('saving')
ani.save('currents_movie.mp4', writer=writer)
print('done')
