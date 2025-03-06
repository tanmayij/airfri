import subprocess
import json
import matplotlib.pyplot as plt
import os
import csv

# Build and run the benchmark (assuming the build directory is './build' and cmake is configured)
def run_benchmark(min_range, max_range, step, benchmark_repetitions, output_file="benchmark_output.json"):
    os.environ['BM_MIN_RANGE'] = str(min_range)
    os.environ['BM_MAX_RANGE'] = str(max_range)
    os.environ['BM_STEP'] = str(step)
    # Command to run the benchmark and output results to JSON format
    benchmark_command = [
        "./build/run_benchmark",  # Replace with your actual executable name
        f"--benchmark_repetitions={benchmark_repetitions}",
        f"--benchmark_out={output_file}"
    ]
    
    # Running the benchmark
    subprocess.run(benchmark_command, check=True)
    print(f"Benchmark results saved to {output_file}")

# Function to read and parse the JSON output
def parse_benchmark_output(output_file="benchmark_output.json"):
    with open(output_file, 'r') as f:
        data = json.load(f)
    
    results = {}
    
    for benchmark in data['benchmarks']:
        if (benchmark['aggregate_name'] == 'mean'):
            algo_name = benchmark['run_name'].split('/')[0]
            range_val = int(benchmark['run_name'].split('/')[1])

            
            if algo_name not in results:
                results[algo_name] = {'range': [], 'mean': []}
            
            results[algo_name]['range'].append(range_val)
            results[algo_name]['mean'].append(benchmark['cpu_time'])  # Convert to microseconds

    
    return results

def write_to_csv_per_algorithm(results, path):
    for algo_name, values in results.items():
        output_csv = path + f"{algo_name}.csv"
        with open(output_csv, mode='w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Range', 'Mean CPU Time (microseconds)'])
            for r, mean in zip(values['range'], values['mean']):
                writer.writerow([r, mean])
        print(f"Data for {algo_name} written to {output_csv}")


# Function to plot the benchmark results
def plot_benchmark_results(results):
    plt.figure(figsize=(10, 6))
    markers = ['o', '^', 'v', '<', '>', 'x', 'p', '*', 's', 'D']

    for i, (algo, data) in enumerate(results.items()):
        marker = markers[i % len(markers)]
        plt.plot(data['range'], data['mean'], label=algo, marker=marker, markersize=3, linewidth=1)
    
    plt.xlabel('Range (log2 size)')
    plt.ylabel('Time (microseconds)')
    plt.yscale('log')
    plt.title('Benchmark Comparison of Additive FFT Algorithms (Logarithmic Scale)')
    plt.legend()
    plt.grid(True, which="both", ls="--")
    plt.show()

# Main script execution
if __name__ == "__main__":
    output_root = "build/"
    output_file =  output_root + "benchmark_output.json" #"../data_backup/benchmark_output.json"
    min_range = 5
    max_range = 20
    step = 5
    benchmark_repetitions = 100
    # Run the benchmark if the output file does not exist
    if not os.path.exists(output_file):
        run_benchmark(min_range, max_range, step, benchmark_repetitions, output_file)
    
    # Parse the benchmark output
    benchmark_results = parse_benchmark_output(output_file)
    write_to_csv_per_algorithm(benchmark_results, output_root)
    
    # Plot the results
    plot_benchmark_results(benchmark_results)
