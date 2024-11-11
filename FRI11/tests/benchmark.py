import subprocess
import re
import statistics

def run_benchmark(command, iterations=10):
    times = []

    # Regex pattern to extract the time from your program's output
    time_pattern = re.compile(r"(\d+\.\d+)\s*seconds")

    for i in range(iterations):
        print(f"Running iteration {i + 1}/{iterations}...")

        # Run your program and capture its output
        #try:
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error in iteration {i + 1}: {result.stderr}")
            continue

        # Extract the time from the program's output
        match = time_pattern.search(result.stdout)
        if match:
            elapsed_time = float(match.group(1))
            times.append(elapsed_time)
            print(f"Iteration {i + 1} completed in {elapsed_time:.6f} seconds")
        else:
            print(f"Time not found in iteration {i + 1} output: {result.stdout}")

        # except Exception as e:
        #     print(f"Exception occurred in iteration {i + 1}: {e}")

    # Calculate statistics
    if times:
        average_time = statistics.mean(times)
        min_time = min(times)
        max_time = max(times)

        print("\nBenchmark Results:")
        print(f"Total iterations: {iterations}")
        print(f"Average time: {average_time:.6f} seconds")
        print(f"Min time: {min_time:.6f} seconds")
        print(f"Max time: {max_time:.6f} seconds")
    else:
        print("No valid times collected.")

if __name__ == "__main__":
    command = "./fri"
    run_benchmark(command)