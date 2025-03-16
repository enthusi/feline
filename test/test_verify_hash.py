import hashlib
import typing
import numpy as np

# Hardcoded expected SHA-256 hash
EXPECTED_HASH_SORTED_CATALOG = "21d8cd338f8dca4fbc7a8b0c3015a41ace5c984c59a1854158997a68f5925790"
EXPECTED_HASH_RAW_ARRAY = "4edb1ab22c28737dcc3c44a76e7e15693d749f7e32bd012cd903016e08e700bd"

# Path to the file to check
FILE_PATH_SORTED_CATALOG = "../data/runtime_files/sorted_catalog.txt"
FILE_PATH_RAW_ARRAY = "../data/runtime_files/feline_float32_array.raw"
FILE_PATH_ROUNDED_ARRAY = "../data/runtime_files/rounded_feline_float32_array.raw"


def calculate_sha256(file_path: str) -> typing.Union[str, None]:
    """Calculate the SHA-256 hash of a file."""
    sha256 = hashlib.sha256()
    try:
        with open(file_path, 'rb') as f:
            for chunk in iter(lambda: f.read(4096), b""):
                sha256.update(chunk)
        return sha256.hexdigest()
    except FileNotFoundError:
        print("Error: File not found.")
        return None


def show_result(expected_hash: str, file_path: str) -> None:
    """Compare hashes and show result."""
    print(f"expected hash:                {expected_hash}")
    file_hash = calculate_sha256(file_path)
    print(f"hash value:                   {file_hash}")
    if file_hash is None:
        print("❌Failed: File does not exist.")
    elif file_hash == expected_hash:
        print("✅Passed: File hash matches expected value.")
    else:
        print("❌Failed: File hash does not match expected value.")


def round_raw_array(file_path: str) -> None:
    """Read and round the raw array to the first decimal place."""
    np.round(np.fromfile(file_path, dtype=np.float32), 1).astype(np.float32).tofile(FILE_PATH_ROUNDED_ARRAY)


if __name__ == "__main__":
    print("Hash Test Raw Array")
    round_raw_array(FILE_PATH_RAW_ARRAY)
    show_result(EXPECTED_HASH_RAW_ARRAY, FILE_PATH_ROUNDED_ARRAY)
    print()
    print("Hash Test Sorted Catalog Array")
    show_result(EXPECTED_HASH_SORTED_CATALOG, FILE_PATH_SORTED_CATALOG)
