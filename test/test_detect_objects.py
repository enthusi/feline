# test/test_detect_objects.py
import numpy as np
import sys
from unittest.mock import Mock

sys.path.insert(0, '/home/runner/work/feline/src/postprocessing')
from feline.src.postprocessing.detect_objects import gauss2d
from feline.src.postprocessing.detect_objects import twoD_Gaussian
from feline.src.postprocessing.detect_objects import onclick


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
	expected_result = offset + amplitude * np.exp(
		- (a * ((x - xo) ** 2) + 2 * b * (x - xo) * (y - yo) + c * ((y - yo) ** 2))).ravel()

	# Check the result against the expected values
	assert np.array_equal(result, expected_result)


def test_gauss2d_zero_parameters():
	# Define parameters for testing
	xy = (0.0, 0.0)
	amp = 0.0
	x0 = 0.0
	y0 = 0.0
	a = 0.0
	b = 0.0
	c = 0.0

	# Call the gauss2d method with the defined parameters
	result = gauss2d(xy, amp, x0, y0, a, b, c)

	# Check the result against the expected values
	assert result == 0.0


def test_gauss2d_negative_parameters():
	# Define parameters for testing
	xy = (-1.0, -1.0)
	amp = -1.0
	x0 = -1.0
	y0 = -1.0
	a = -1.0
	b = -1.0
	c = -1.0

	# Call the gauss2d method with the defined parameters
	result = gauss2d(xy, amp, x0, y0, a, b, c)

	# Check the result against the expected values
	assert result == amp


def test_onclick_prints_correctly(capfd):
	# Create a mock event
	event = Mock()

	# Set the attributes of the mock event that onclick uses
	event.xdata = 1.0
	event.ydata = 2.0
	event.button = 1
	event.x = 10
	event.y = 20

	# Call the onclick method with the mock event
	onclick(event)

	# Capture the output
	out, err = capfd.readouterr()

	# Check the output
	assert out == "#button=1, x=10, y=20, xdata=1.000000, ydata=2.000000\n"


def test_twoD_Gaussian_with_zero_amplitude_returns_offset():
	# Define parameters for testing
	amplitude = 0.0
	xo = 0.0
	yo = 0.0
	sigma_x = 1.0
	sigma_y = 1.0
	theta = 0.0
	offset = 1.0

	# Create a grid of coordinates for testing
	x, y = np.meshgrid(np.linspace(-5, 5, 11), np.linspace(-5, 5, 11))
	coordinates = np.vstack((x.ravel(), y.ravel()))

	# Call the twoD_Gaussian method with the defined parameters
	result = twoD_Gaussian(coordinates, amplitude, xo, yo, sigma_x, sigma_y, theta, offset)

	# Check the result against the expected values
	assert np.allclose(result, np.full_like(coordinates, offset), atol=1e-8)
