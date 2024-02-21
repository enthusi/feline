import struct
from unittest.mock import patch

import numpy as np
import sys
import os
import mpdaf


from pytest_mock import mocker

sys.path.insert(0, '/home/runner/work/feline/src/postprocessing')
from feline.src.postprocessing.combination import load_data
from feline.src.postprocessing.combination import create_masking_plot
from feline.src.postprocessing.combination import write_to_file
from feline.src.postprocessing.combination import process_cube_data


def test_load_data_returns_correct_sum():
    file = 'test_file.fits'
    ext = 1
    result = load_data(file, ext)
    assert isinstance(result, mpdaf.obj.Image)
    assert str(result) == "<Image(shape=(10, 10), unit='1e-20 erg / (Angstrom s cm2)', dtype='float64')>"   # replace expected_sum with the expected value


@patch('mpdaf.obj.Image.write')
def test_create_masking_plot(mock_write, mocker, tmpdir):
    mocker.patch('sys.argv', ['combination.py', 'test_file.fits'])
    file = tmpdir.mkdir('data').join('image00.fits')
    mocker.patch('os.path.join', return_value='data/raw/test_file.fits')
    mock_write().write('image00.fits')
    create_masking_plot()
    assert mock_write.called


def test_write_to_file_writes_correct_data(tmp_path):
    file_path = tmp_path / 'test_file.dat'
    data = 123.456
    write_to_file(file_path, data)
    with open(file_path, 'rb') as fin:
        result = struct.unpack('f', fin.read())
    assert round(result[0], 3) == data


def test_process_cube_data_writes_correct_data(tmp_path):
    cube_data = np.random.rand(10, 10, 10)
    dy, dx = cube_data.shape[1:]
    file_path = tmp_path / 'test_file.dat'
    process_cube_data(cube_data, dy, dx, file_path)
    with open(file_path, 'rb') as fin:
        result = struct.unpack('f' * cube_data.size, fin.read())
    assert np.allclose(np.sort(result), np.sort(cube_data.flatten()), rtol=0.01)
