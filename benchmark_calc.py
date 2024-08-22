import re

def parse_time(text):
    """Extracts and converts time from a text string into seconds."""
    match = re.search(r'Run time: ([\d\.]+) (seconds|ns)', text)
    if not match:
        raise ValueError("Run time not found in text")
    time, unit = match.groups()
    time = float(time)
    if unit == 'ns':
        time /= 1e9  # Convert nanoseconds to seconds
    return time

def extract_times(file_path):
    """Reads the file and extracts all run times."""
    with open(file_path, 'r') as file:
        text = file.read()
    # Find all run times in the text
    times = re.findall(r'Run time: ([\d\.]+) (seconds|ns)', text)
    if not times:
        raise ValueError("No run times found in file")
    # Convert all times to seconds
    times_in_seconds = [parse_time(f"Run time: {time} {unit}") for time, unit in times]
    return times_in_seconds

def get_cores(file_path):
    """Extracts the number of CPU cores from the file."""
    with open(file_path, 'r') as file:
        text = file.read()
    match = re.search(r'CPU Cores: (\d+)', text)
    if not match:
        raise ValueError("Number of CPU cores not found in file")
    return int(match.group(1))

def average(lst):
    """Returns the average of a list of numbers."""
    return sum(lst) / len(lst)

def main(single_core_file, multi_core_file):
    # Extract and calculate average times
    single_core_times = extract_times(single_core_file)
    multi_core_times = extract_times(multi_core_file)
    average_single_core_time = average(single_core_times)
    average_multi_core_time = average(multi_core_times)
    
    # Extract the number of cores and calculate speedup and efficiency
    num_cores = get_cores(multi_core_file)
    speedup = average_single_core_time / average_multi_core_time
    efficiency = speedup / num_cores
    
    # Print results
    print(f"Average single-core time: {average_single_core_time:.6f} seconds")
    print(f"Average multi-core time: {average_multi_core_time:.6f} seconds")
    print(f"Average speedup: {speedup:.6f}")
    print(f"Average efficiency: {efficiency:.6f}")

# Example usage
if __name__ == "__main__":
    single_core_file = 'benchmarks/benchmark_results_11_04_35.txt'
    multi_core_file = 'benchmarks/benchmark_results_10_35_54.txt'
    main(single_core_file, multi_core_file)

