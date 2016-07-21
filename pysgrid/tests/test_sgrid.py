'''
Created on Apr 7, 2015

@author: ayan
'''

from __future__ import (absolute_import, division, print_function)

import os

import pytest
import numpy as np

from ..sgrid import SGrid, load_grid
from .write_nc_test_files import (deltares_sgrid,
                                  deltares_sgrid_no_optional_attr,
                                  non_compliant_sgrid,
                                  roms_sgrid,
                                  wrf_sgrid_2d)


CURRENT_DIR = os.path.dirname(__file__)
TEST_FILES = os.path.join(CURRENT_DIR, 'files')
write_path = os.path.join(CURRENT_DIR, 'test_sgrid_write.nc')


def test_exception_raised(non_compliant_sgrid):
    with pytest.raises(ValueError):
        load_grid(non_compliant_sgrid)


# TestSGridCreate.
def test_load_from_file(roms_sgrid):
    fname = 'tmp_sgrid_roms.nc'  # FIXME: name handling in the fixture.
    sg_obj = load_grid(fname)
    assert isinstance(sg_obj, SGrid)


def test_load_from_dataset(roms_sgrid):
    sg_obj = load_grid(roms_sgrid)
    assert isinstance(sg_obj, SGrid)


@pytest.fixture
def sgrid_roms_sgrid(roms_sgrid):
    return load_grid(roms_sgrid)


def test_center_lon(sgrid_roms_sgrid):
    center_lon = sgrid_roms_sgrid.center_lon
    assert center_lon.shape == (4, 4)


def test_center_lat(sgrid_roms_sgrid):
    center_lat = sgrid_roms_sgrid.center_lat
    assert center_lat.shape == (4, 4)


def test_variables(sgrid_roms_sgrid):
    dataset_vars = sgrid_roms_sgrid.variables
    expected_vars = [u's_rho',
                     u's_w',
                     u'time',
                     u'xi_rho',
                     u'eta_rho',
                     u'xi_psi',
                     u'eta_psi',
                     u'xi_u',
                     u'eta_u',
                     u'xi_v',
                     u'eta_v',
                     u'grid',
                     u'u',
                     u'v',
                     u'fake_u',
                     u'lon_rho',
                     u'lat_rho',
                     u'lon_psi',
                     u'lat_psi',
                     u'lat_u',
                     u'lon_u',
                     u'lat_v',
                     u'lon_v',
                     u'salt',
                     u'zeta']
    assert len(dataset_vars) == len(expected_vars)
    assert dataset_vars == expected_vars


def test_grid_variables(sgrid_roms_sgrid):
    dataset_grid_variables = sgrid_roms_sgrid.grid_variables
    expected_grid_variables = [u'u', u'v', u'fake_u', u'salt']
    assert len(dataset_grid_variables) == len(expected_grid_variables)
    assert set(dataset_grid_variables) == set(expected_grid_variables)


def test_non_grid_variables(sgrid_roms_sgrid):
    dataset_non_grid_variables = sgrid_roms_sgrid.non_grid_variables
    expected_non_grid_variables = [u's_rho',
                                   u's_w',
                                   u'time',
                                   u'xi_rho',
                                   u'eta_rho',
                                   u'xi_psi',
                                   u'eta_psi',
                                   u'xi_u',
                                   u'eta_u',
                                   u'xi_v',
                                   u'eta_v',
                                   u'grid',
                                   u'lon_rho',
                                   u'lat_rho',
                                   u'lon_psi',
                                   u'lat_psi',
                                   u'lat_u',
                                   u'lon_u',
                                   u'lat_v',
                                   u'lon_v',
                                   u'zeta']
    assert len(dataset_non_grid_variables) == len(expected_non_grid_variables)
    assert set(dataset_non_grid_variables) == set(expected_non_grid_variables)


def test_variable_slicing(sgrid_roms_sgrid):
    u_center_slices = sgrid_roms_sgrid.u.center_slicing
    v_center_slices = sgrid_roms_sgrid.v.center_slicing
    u_center_expected = (np.s_[:], np.s_[:], np.s_[1:-1], np.s_[:])
    v_center_expected = (np.s_[:], np.s_[:], np.s_[:], np.s_[1:-1])
    assert u_center_slices == u_center_expected
    assert v_center_slices == v_center_expected


def test_grid_variable_average_axes(sgrid_roms_sgrid):
    uc_axis = sgrid_roms_sgrid.u.center_axis
    uc_axis_expected = 1
    un_axis = sgrid_roms_sgrid.u.node_axis
    un_axis_expected = 0
    lon_rho_c_axis = sgrid_roms_sgrid.lon_rho.center_axis
    lon_rho_n_axis = sgrid_roms_sgrid.lon_rho.node_axis
    assert uc_axis == uc_axis_expected
    assert un_axis == un_axis_expected
    assert lon_rho_c_axis is None
    assert lon_rho_n_axis is None


def test_optional_grid_attrs(sgrid_roms_sgrid):
    face_coordinates = sgrid_roms_sgrid.face_coordinates
    node_coordinates = sgrid_roms_sgrid.node_coordinates
    edge1_coordinates = sgrid_roms_sgrid.edge1_coordinates
    edge2_coordinates = sgrid_roms_sgrid.edge2_coordinates
    fc_expected = ('lon_rho', 'lat_rho')
    nc_expected = ('lon_psi', 'lat_psi')
    e1c_expected = ('lon_u', 'lat_u')
    e2c_expected = ('lon_v', 'lat_v')
    assert face_coordinates == fc_expected
    assert node_coordinates == nc_expected
    assert edge1_coordinates == e1c_expected
    assert edge2_coordinates == e2c_expected


def test_write_sgrid_to_netcdf(sgrid_roms_sgrid):
    sgrid_roms_sgrid.save_as_netcdf(write_path)
    result = load_grid(write_path)
    os.remove(write_path)
    assert isinstance(result, SGrid)  # TODO: Add more "round-trip" tests.


# TestSGridNoCoordinates
"""
Test to make sure that if no coordinates (e.g. face, edge1, etc)
are specified, those coordinates can be inferred from the dataset.

A file is representing a delft3d dataset is used for this test.

"""


@pytest.fixture
def sgrid_deltares(deltares_sgrid_no_optional_attr):
    return load_grid(deltares_sgrid_no_optional_attr)


def test_face_coordinate_inference(sgrid_deltares):
    face_coordinates = sgrid_deltares.face_coordinates
    expected_face_coordinates = (u'XZ', u'YZ')
    assert face_coordinates == expected_face_coordinates


def test_center_lon_deltares_no_coord(sgrid_deltares):
    center_lon = sgrid_deltares.center_lon
    assert center_lon.shape == (4, 4)


def test_center_lat_deltares_no_coord(sgrid_deltares):
    center_lat = sgrid_deltares.center_lat
    assert center_lat.shape == (4, 4)


def test_node_lon(sgrid_deltares):
    node_lon = sgrid_deltares.node_lon
    assert node_lon.shape == (4, 4)


def test_node_lat(sgrid_deltares):
    node_lat = sgrid_deltares.node_lat
    assert node_lat.shape == (4, 4)


def test_grid_angles(sgrid_deltares):
    angles = sgrid_deltares.angles
    angles_shape = (4, 4)
    assert angles.shape == angles_shape


# TestSGridWRFDataset
"""
Test a representative WRF file.

"""


@pytest.fixture
def sgrid_wrf(wrf_sgrid_2d):
    return load_grid(wrf_sgrid_2d)


def test_topology_dimension(sgrid_wrf):
    topology_dim = sgrid_wrf.topology_dimension
    expected_dim = 2
    assert topology_dim == expected_dim


def test_variable_slicing_wrf(sgrid_wrf):
    u_slice = sgrid_wrf.U.center_slicing
    u_expected = (np.s_[:], np.s_[:], np.s_[:], np.s_[:])
    v_slice = sgrid_wrf.V.center_slicing
    v_expected = (np.s_[:], np.s_[:], np.s_[:], np.s_[:])
    assert u_slice == u_expected
    assert v_slice == v_expected


def test_variable_average_axes(sgrid_wrf):
    u_avg_axis = sgrid_wrf.U.center_axis
    u_axis_expected = 1
    v_avg_axis = sgrid_wrf.V.center_axis
    v_axis_expected = 0
    assert u_avg_axis == u_axis_expected
    assert v_avg_axis == v_axis_expected


# TestSGridDelft3dDataset
"""
Test using a representative delft3d file.

"""


@pytest.fixture
def sgrid_obj_deltares(deltares_sgrid):
    return load_grid(deltares_sgrid)


def test_center_lon_deltares(sgrid_obj_deltares):
    center_lon = sgrid_obj_deltares.center_lon
    assert center_lon.shape == (4, 4)


def test_center_lat_deltares(sgrid_obj_deltares):
    center_lat = sgrid_obj_deltares.center_lat
    assert center_lat.shape == (4, 4)


def test_topology_dimension_deltares(sgrid_obj_deltares):
    topology_dim = sgrid_obj_deltares.topology_dimension
    expected_dim = 2
    assert topology_dim == expected_dim


def test_variable_slice(sgrid_obj_deltares):
    u_center_slices = sgrid_obj_deltares.U1.center_slicing
    v_center_slices = sgrid_obj_deltares.V1.center_slicing
    u_center_expected = (np.s_[:], np.s_[:], np.s_[:], np.s_[1:])
    v_center_expected = (np.s_[:], np.s_[:], np.s_[1:], np.s_[:])
    xz_center_slices = sgrid_obj_deltares.XZ.center_slicing
    xcor_center_slices = sgrid_obj_deltares.XCOR.center_slicing
    xz_center_expected = (np.s_[1:], np.s_[1:])
    xcor_center_expected = (np.s_[:], np.s_[:])
    assert u_center_slices == u_center_expected
    assert v_center_slices == v_center_expected
    assert xz_center_slices == xz_center_expected
    assert xcor_center_slices == xcor_center_expected


def test_averaging_axes(sgrid_obj_deltares):
    u1c_axis = sgrid_obj_deltares.U1.center_axis
    u1c_expected = 0
    v1n_axis = sgrid_obj_deltares.V1.node_axis
    v1n_expected = 0
    latitude_c_axis = sgrid_obj_deltares.latitude.center_axis
    latitude_n_axis = sgrid_obj_deltares.latitude.node_axis
    assert u1c_axis == u1c_expected
    assert v1n_axis == v1n_expected
    assert latitude_c_axis is None
    assert latitude_n_axis is None


def test_grid_optional_attrs(sgrid_obj_deltares):
    face_coordinates = sgrid_obj_deltares.face_coordinates
    node_coordinates = sgrid_obj_deltares.node_coordinates
    edge1_coordinates = sgrid_obj_deltares.edge1_coordinates
    edge2_coordinates = sgrid_obj_deltares.edge2_coordinates
    fc_expected = ('XZ', 'YZ')
    nc_expected = ('XCOR', 'YCOR')
    assert face_coordinates == fc_expected
    assert node_coordinates == nc_expected
    assert edge1_coordinates is None
    assert edge2_coordinates is None


def test_grid_variables_deltares(sgrid_obj_deltares):
    grid_variables = sgrid_obj_deltares.grid_variables
    expected_grid_variables = [u'U1', u'V1', u'FAKE_U1', u'W', u'FAKE_W']
    assert set(grid_variables) == set(expected_grid_variables)


def test_angles(sgrid_obj_deltares):
    angles = sgrid_obj_deltares.angles
    expected_shape = (4, 4)
    assert angles.shape == expected_shape


def test_no_3d_attributes(sgrid_obj_deltares):
    assert not hasattr(sgrid_obj_deltares, 'volume_padding')
    assert not hasattr(sgrid_obj_deltares, 'volume_dimensions')
    assert not hasattr(sgrid_obj_deltares, 'volume_coordinates')
    assert not hasattr(sgrid_obj_deltares, 'face1_padding')
    assert not hasattr(sgrid_obj_deltares, 'face1_coordinates')
    assert not hasattr(sgrid_obj_deltares, 'face1_dimensions')
    assert not hasattr(sgrid_obj_deltares, 'face2_padding')
    assert not hasattr(sgrid_obj_deltares, 'face2_coordinates')
    assert not hasattr(sgrid_obj_deltares, 'face2_dimensions')
    assert not hasattr(sgrid_obj_deltares, 'face3_padding')
    assert not hasattr(sgrid_obj_deltares, 'face3_coordinates')
    assert not hasattr(sgrid_obj_deltares, 'edge3_padding')
    assert not hasattr(sgrid_obj_deltares, 'edge3_coordinates')
    assert not hasattr(sgrid_obj_deltares, 'edge3_dimensions')


def test_2d_attributes(sgrid_obj_deltares):
    assert hasattr(sgrid_obj_deltares, 'face_padding')
    assert hasattr(sgrid_obj_deltares, 'face_coordinates')
    assert hasattr(sgrid_obj_deltares, 'face_dimensions')
    assert hasattr(sgrid_obj_deltares, 'vertical_padding')
    assert hasattr(sgrid_obj_deltares, 'vertical_dimensions')


# TestSGridSaveNoNodeCoordinates
"""
Test that SGrid.save_as_netcdf is saving content
when there are no nodes or node coordinates specified.

This scenario will typically occur with WRF datasets.

"""


def round_trip_wrf(wrf_sgrid_2d):
    sgrid_target = os.path.join(TEST_FILES, 'tmp_sgrid.nc')
    sg_obj = load_grid(wrf_sgrid_2d)
    sg_obj.save_as_netcdf(sgrid_target)
    os.remove(sgrid_target)


@pytest.fixture
def sgrid_wrf_sgrid_2d(wrf_sgrid_2d):
    return load_grid(wrf_sgrid_2d)


def test_sgrid(sgrid_wrf_sgrid_2d):
    assert isinstance(sgrid_wrf_sgrid_2d, SGrid)


def test_nodes(sgrid_wrf_sgrid_2d):
    node_lon = sgrid_wrf_sgrid_2d.node_lon
    node_lat = sgrid_wrf_sgrid_2d.node_lat
    assert node_lon is None
    assert node_lat is None


def test_node_coordinates(sgrid_wrf_sgrid_2d):
    node_coordinates = sgrid_wrf_sgrid_2d.node_coordinates
    assert node_coordinates is None


def test_node_dimesnions(sgrid_wrf_sgrid_2d):
    node_dims = sgrid_wrf_sgrid_2d.node_dimensions
    expected = 'west_east_stag south_north_stag'
    assert node_dims == expected


# TestSGridSaveNodeCoordinates
"""
Test that SGrid.save_as_netcdf is saving
content correctly.

There maybe a better way to do this using
mocks, but this will do for now.

"""


def round_trip_deltares(deltares_sgrid):
    sgrid_target = os.path.join(TEST_FILES, 'tmp_sgrid.nc')
    sg_obj = load_grid(deltares_sgrid)
    sg_obj.save_as_netcdf(sgrid_target)
    os.remoce(sgrid_target)


@pytest.fixture
def sgrid_deltares_sgrid(deltares_sgrid):
    return load_grid(deltares_sgrid)


def test_save_as_netcdf(sgrid_deltares_sgrid):
    """
    Test that the attributes in the
    saved netCDF file are as expected.

    """
    target_dims = sgrid_deltares_sgrid.dimensions
    expected_target_dims = [(u'MMAXZ', 4),
                            (u'NMAXZ', 4),
                            (u'MMAX', 4),
                            (u'NMAX', 4),
                            (u'KMAX', 2),
                            (u'KMAX1', 3),
                            (u'time', 2)
                            ]
    target_vars = sgrid_deltares_sgrid.variables
    expected_target_vars = [u'XZ',
                            u'YZ',
                            u'XCOR',
                            u'YCOR',
                            u'grid',
                            u'U1',
                            u'FAKE_U1',
                            u'V1',
                            u'W',
                            u'FAKE_W',
                            u'time',
                            u'latitude',
                            u'longitude',
                            u'grid_latitude',
                            u'grid_longitude'
                            ]
    target_grid_vars = sgrid_deltares_sgrid.grid_variables
    expected_target_grid_vars = [u'U1',
                                 u'FAKE_U1',
                                 u'V1',
                                 u'W',
                                 u'FAKE_W']
    target_face_coordinates = sgrid_deltares_sgrid.face_coordinates
    expected_target_face_coordinates = (u'XZ', u'YZ')
    assert isinstance(sgrid_deltares_sgrid, SGrid)
    assert len(target_dims) == len(expected_target_dims)
    assert set(target_dims) == set(expected_target_dims)
    assert len(target_vars) == len(expected_target_vars)
    assert set(target_vars) == set(expected_target_vars)
    assert len(target_grid_vars) == len(expected_target_grid_vars)
    assert set(target_grid_vars) == set(expected_target_grid_vars)
    assert target_face_coordinates == expected_target_face_coordinates


def test_saved_sgrid_attributes(sgrid_deltares_sgrid):
    """
    Test that calculated/inferred attributes
    are as expected from the saved filed.

    """
    u1_var = sgrid_deltares_sgrid.U1
    u1_var_center_avg_axis = u1_var.center_axis
    expected_u1_center_axis = 0
    u1_vector_axis = u1_var.vector_axis
    expected_u1_vector_axis = 'X'
    original_angles = sgrid_deltares_sgrid.angles
    saved_angles = sgrid_deltares_sgrid.angles
    assert u1_var_center_avg_axis == expected_u1_center_axis
    assert u1_vector_axis == expected_u1_vector_axis
    np.testing.assert_almost_equal(original_angles, saved_angles, decimal=3)  # noqa
