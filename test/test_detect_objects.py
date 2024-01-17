import numpy as np
import astropy.wcs
from feline.src.afterprocess.detect_objects import world_to_pix

class MockWCS:
    def wcs_world2pix(self, radarray, dummy):
        # Implement a simple mock behavior
        return np.array([[10.0, 20.0]])

def test_world_to_pix():
    # Create an instance of the MockWCS class
    mock_coord = MockWCS()

    # Define a radius for testing
    radius = [30.0, 40.0]

    # Call the world_to_pix method with the mock data
    result = world_to_pix(mock_coord, radius)

    # Check the result against the expected values
    expected_result = (10.0, 20.0)
    assert result == expected_result
