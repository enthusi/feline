import hashlib

# Hardcoded expected SHA-256 hash
EXPECTED_HASH_SORTED_CATALOG = "8478f7da42d44e9d99bff025d76e07b1fd8fdae814327621d5f99621322cbd59"
EXPECTED_HASH_RAW_ARRAY = "ef43be7379c84430a48095d9ebe4aa4043673d43443e524753301c87f702de44"

# Path to the file to check
FILE_PATH_SORTED_CATALOG = "../data/runtime_files/sorted_catalog.txt"
FILE_PATH_RAW_ARRAY = "../data/runtime_files/feline_float32_array.raw"


def calculate_sha256(file_path):
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


def show_result(expected_hash, file_path):
    print(f"expected hash raw array:      {expected_hash}")
    file_hash = calculate_sha256(file_path)
    print(f"hash value:                   {file_hash}")
    if file_hash is None:
        print("❌Failed: File does not exist.")
    elif file_hash == expected_hash:
        print("✅Passed: File hash matches expected value.")
    else:
        print("❌Failed: File hash does not match expected value.")


def main():
    show_result(EXPECTED_HASH_RAW_ARRAY, FILE_PATH_RAW_ARRAY)
    print()
    show_result(EXPECTED_HASH_SORTED_CATALOG, FILE_PATH_SORTED_CATALOG)


if __name__ == "__main__":
    main()
