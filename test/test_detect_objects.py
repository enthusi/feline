# test/test_detect_objects.py
import numpy as np
import sys
import unittest.mock
import pytest
from feline.src.afterprocess.detect_objects import gauss2d
from feline.src.afterprocess.detect_objects import twoD_Gaussian
from feline.src.afterprocess.detect_objects import world_to_pix

def test_world_to_pix():
    # Create a WCS object (you should replace this with your actual coord)
    header = {'NAXIS': 2, 'CRVAL1': 0, 'CRVAL2': 0, 'CRPIX1': 1, 'CRPIX2': 1, 'CTYPE1': 'RA---TAN', 'CTYPE2': 'DEC--TAN'}
    coord = WCS(header)

    # Set up sample radius
    rad = [10.0, 20.0]

    # Call the function
    result = world_to_pix(coord, rad)

    # Check the result against expected values
    # Replace these with the actual expected values based on your use case
    expected_x = 1.0
    expected_y = 1.0

    assert result[0] == pytest.approx(expected_x)
    assert result[1] == pytest.approx(expected_y)
    
def test_gauss2d():
    # Define parameters for testing
    xy = (1.0, 2.0)
    amp = 3.0
    x0 = 1.5
    y0 = 2.5
    a = 0.1
    b = 0.2
    c = 0.3

    # Call the gauss2d method with the defined parameters
    result = gauss2d(xy, amp, x0, y0, a, b, c)

    # Calculate the expected result
    inner = a * (xy[0] - x0) ** 2 + 2 * b * (xy[0] - x0) ** 2 * (xy[1] - y0) ** 2 + c * (xy[1] - y0) ** 2
    expected_result = amp * np.exp(-inner)

    # Check the result against the expected values
    assert result == expected_result

def test_twoD_Gaussian():        
    # Define parameters for testing
    amplitude = 1.0
    xo = 0.0
    yo = 0.0
    sigma_x = 1.0
    sigma_y = 1.0
    theta = 0.0
    offset = 0.0

    # Create a grid of coordinates for testing
    x, y = np.meshgrid(np.linspace(-5, 5, 11), np.linspace(-5, 5, 11))
    coordinates = np.vstack((x.ravel(), y.ravel()))

    # Call the twoD_Gaussian method with the defined parameters
    result = twoD_Gaussian(coordinates, amplitude, xo, yo, sigma_x, sigma_y, theta, offset)

    # Calculate the expected result
    a = (np.cos(theta) ** 2) / (2 * sigma_x ** 2) + (np.sin(theta) ** 2) / (2 * sigma_y ** 2)
    b = -(np.sin(2 * theta)) / (4 * sigma_x ** 2) + (np.sin(2 * theta)) / (4 * sigma_y ** 2)
    c = (np.sin(theta) ** 2) / (2 * sigma_x ** 2) + (np.cos(theta) ** 2) / (2 * sigma_y ** 2)
    expected_result = offset + amplitude * np.exp(- (a * ((x - xo) ** 2) + 2 * b * (x - xo) * (y - yo) + c * ((y - yo) ** 2))).ravel()

    # Check the result against the expected values
    assert np.array_equal(result, expected_result)
