import numpy as np
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from afterprocess.detect_objects import world_to_pix
import afterprocess.project_path_config


def test_world_to_pix():
    # Assuming you have a Coord object with wcs_world2pix method
    class MockCoord:
        def wcs_world2pix(self, radarray, _):
            # Mock implementation for testing
            return np.array([[10, 20]])

    coord = MockCoord()
    rad = [1, 1]

    result = world_to_pix(coord, rad)

    # Assert the expected result based on the mock implementation
    assert result == (10, 20)
