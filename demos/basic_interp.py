import pysgrid
import numpy as np
import matplotlib.pyplot as plt

node_lon = np.array(([1, 3, 5], [1, 3, 5], [1, 3, 5]))
node_lat = np.array(([1, 1, 1], [3, 3, 3], [5, 5, 5]))
edge2_lon = np.array(([0, 2, 4, 6], [0, 2, 4, 6], [0, 2, 4, 6]))
edge2_lat = np.array(([1, 1, 1, 1], [3, 3, 3, 3], [5, 5, 5, 5]))
edge1_lon = np.array(([1, 3, 5], [1, 3, 5], [1, 3, 5], [1, 3, 5]))
edge1_lat = np.array(([0, 0, 0], [2, 2, 2], [4, 4, 4], [6, 6, 6]))
center_lon = np.array(([0, 2, 4, 6], [0, 2, 4, 6], [0, 2, 4, 6], [0, 2, 4, 6]))
center_lat = np.array(([0, 0, 0, 0], [2, 2, 2, 2], [4, 4, 4, 4], [6, 6, 6, 6]))

sgrid = pysgrid.SGrid(node_lon=node_lon,
                      node_lat=node_lat,
                      edge1_lon=edge1_lon,
                      edge1_lat=edge1_lat,
                      edge2_lon=edge2_lon,
                      edge2_lat=edge2_lat,
                      center_lon=center_lon,
                      center_lat=center_lat)

c_var = np.array(([0, 0, 0, 0], [0, 1, 2, 0], [0, 2, 1, 0], [0, 0, 0, 0]))
e2_var = np.array(([1, 0, 0, 1], [0, 1, 2, 0], [0, 0, 0, 0]))
e1_var = np.array(([1, 1, 0], [0, 1, 0], [0, 2, 0], [1, 1, 0]))
n_var = np.array(([0, 1, 0], [1, 0, 1], [0, 1, 0]))

ptsx, ptsy = np.mgrid[0:6:600j, 0:6:600j]
pts = np.stack((ptsx, ptsy), axis=-1)


interp_c = sgrid.interpolate_var_to_points(pts, c_var).reshape(600, 600)
interp_e1 = sgrid.interpolate_var_to_points(pts, e1_var).reshape(600, 600).T
interp_e2 = sgrid.interpolate_var_to_points(pts, e2_var).reshape(600, 600).T
interp_n = sgrid.interpolate_var_to_points(pts, n_var).reshape(600, 600)

plt.subplot(221)
plt.imshow(interp_c, extent=(0, 6, 0, 6), origin='lower')
plt.vlines(center_lon, center_lat[0], center_lat[-1])
plt.hlines(center_lon, center_lat[0], center_lat[-1])
plt.title('rho grid interpolation')

plt.subplot(222)
plt.imshow(interp_e1, extent=(0, 6, 0, 6), origin='lower')
plt.vlines(edge2_lon, center_lat[0], center_lat[-1])
plt.hlines(center_lon, edge1_lat[0], edge1_lat[-1])
plt.title('u grid interpolation')

plt.subplot(223)
plt.imshow(interp_e2, extent=(0, 6, 0, 6), origin='lower')
plt.vlines(center_lon, node_lat[0], node_lat[-1])
plt.hlines(edge2_lon, center_lat[0], center_lat[-1])
plt.title('v grid interpolation')

plt.subplot(224)
plt.imshow(interp_n, extent=(0, 6, 0, 6), origin='lower')
plt.vlines(node_lon, node_lat[0], node_lat[-1])
plt.hlines(node_lon, node_lat[0], node_lat[-1])
plt.title('psi grid interpolation')
plt.show()
