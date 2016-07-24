"""
Created on Apr 7, 2015

@author: ayan

"""

from __future__ import (absolute_import, division, print_function)

import pytest
import numpy as np

from ..sgrid import SGrid, load_grid

from .write_nc_test_files import wrf_sgrid_2d


"""
Test SGrid WRF Dataset.

"""


@pytest.fixture
def sgrid(wrf_sgrid_2d):
    return load_grid(wrf_sgrid_2d)


def test_topology_dimension(sgrid):
    topology_dim = sgrid.topology_dimension
    expected_dim = 2
    assert topology_dim == expected_dim


def test_variable_slicing_wrf(sgrid):
    u_slice = sgrid.U.center_slicing
    u_expected = (np.s_[:], np.s_[:], np.s_[:], np.s_[:])
    v_slice = sgrid.V.center_slicing
    v_expected = (np.s_[:], np.s_[:], np.s_[:], np.s_[:])
    assert u_slice == u_expected
    assert v_slice == v_expected


def test_variable_average_axes(sgrid):
    u_avg_axis = sgrid.U.center_axis
    u_axis_expected = 1
    v_avg_axis = sgrid.V.center_axis
    v_axis_expected = 0
    assert u_avg_axis == u_axis_expected
    assert v_avg_axis == v_axis_expected


"""
Test SGrid Save No-Node Coordinates.

This scenario will typically occur with WRF datasets.

"""


def test_roundtrip(wrf_sgrid_2d, tmpdir):
    """
    TODO: add more "round-trip" tests.

    """
    fname = tmpdir.mkdir('files').join('wrf_roundtrip.nc')
    sg_obj = load_grid(wrf_sgrid_2d)
    sg_obj.save_as_netcdf(fname)


@pytest.fixture
def sgrid_sgrid_2d(wrf_sgrid_2d):
    return load_grid(wrf_sgrid_2d)


def test_sgrid(sgrid_sgrid_2d):
    assert isinstance(sgrid_sgrid_2d, SGrid)


def test_nodes(sgrid_sgrid_2d):
    node_lon = sgrid_sgrid_2d.node_lon
    node_lat = sgrid_sgrid_2d.node_lat
    assert node_lon is None
    assert node_lat is None


def test_node_coordinates(sgrid_sgrid_2d):
    node_coordinates = sgrid_sgrid_2d.node_coordinates
    assert node_coordinates is None


def test_node_dimensions(sgrid_sgrid_2d):
    node_dims = sgrid_sgrid_2d.node_dimensions
    expected = 'west_east_stag south_north_stag'
    assert node_dims == expected
