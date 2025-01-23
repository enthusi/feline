import hashlib


def hash_raw_feline_output(file_path: str, algorithm: str ="sha256", chunk_size: int =1048576) -> str:
	hasher = hashlib.new(algorithm)
	with open(file_path, "rb") as f:
		# 1 MB chunks
		while chunk := f.read(chunk_size):
			hasher.update(chunk)
	return hasher.hexdigest()


if __name__ == "__main__":
	print(hash_raw_feline_output("../data/runtime_files/float32_array_omp4.raw"))
