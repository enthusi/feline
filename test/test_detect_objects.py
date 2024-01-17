# test/test_detect_objects.py
import numpy as np
import astropy.wcs
import pytest
from afterprocess.detect_objects import world_to_pix

# Mock Astropy's WCS class for testing
class MockWCS:
    def wcs_world2pix(self, radarray, dummy):
        # Implement a simple mock behavior
        return np.array([[10.0, 20.0]])

def test_world_to_pix_with_dummy_file(monkeypatch):
    # Create an instance of the MockWCS class
    mock_coord = MockWCS()

    # Define a radius for testing
    radius = [30.0, 40.0]

    # Use monkeypatch to set sys.argv[1] to "dummy.fits" during the test
    monkeypatch.setattr("sys.argv", ["script_name.py", "dummy.fits"])

    # Call the world_to_pix method with the modified sys.argv
    result = world_to_pix(mock_coord, radius)

    # Check the result against the expected values
    expected_result = (10.0, 20.0)
    assert result == expected_result
