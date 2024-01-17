import numpy as np
from feline.src.afterprocess.detect_objects import world_to_pix


try:
    import matplotlib.pyplot as plt
    plt.switch_backend('agg')  # Use a non-interactive backend
except ImportError:
    pass
    
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
