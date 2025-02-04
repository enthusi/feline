Verifying Code Integrity After Changes
======================================

To ensure that modifications do not affect expected results, we provide a verification process using SHA-256 hashing.
Instead of storing the raw data array and the reference catalog, we store only their hashes.

After running the software with the test data cube ``DATACUBE_UDF-10.fits``, users can verify the output by navigating to the ``test`` directory and executing:

.. code-block:: bash

    cd test
    python test_verify_hash.py

This script generates a hash of the raw data array and the produced catalog and compares it to the stored reference hashes,
confirming that the output remains consistent.
