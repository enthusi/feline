import hashlib

# Hardcoded expected SHA-256 hash
EXPECTED_HASH = "528df34011402095eee3833c35abe34c14a5687bcf4e047b193654ac7425d3a5"

# Path to the file to check
FILE_PATH = "../data/runtime_files/sorted_catalog.txt"

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

def main():
    print(f"expected hash: {EXPECTED_HASH}")
    file_hash = calculate_sha256(FILE_PATH)
    print(f"hash value:    {file_hash}")
    if file_hash is None:
        print("❌Failed: File does not exist.")
    elif file_hash == EXPECTED_HASH:
        print("✅Passed: File hash matches expected value.")
    else:
        print("❌Failed: File hash does not match expected value.")

if __name__ == "__main__":
    main()
