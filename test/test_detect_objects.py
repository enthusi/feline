# test/test_detect_objects.py
import numpy as np
import sys
import unittest.mock
import pytest

def test_gauss2d():
    with unittest.mock.patch.object(sys, 'argv', ["detect_objects.py", "/home/runner/work/feline/feline/data/processed/dummy.fits"]):
        # Now, import the module that uses sys.argv
        from feline.src.afterprocess.detect_objects import gauss2d
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
    with unittest.mock.patch.object(sys, 'argv', ["detect_objects.py", "/home/runner/work/feline/feline/data/processed/dummy.fits"]):
        # Now, import the module that uses sys.argv
        from feline.src.afterprocess.detect_objects import twoD_Gaussian
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
