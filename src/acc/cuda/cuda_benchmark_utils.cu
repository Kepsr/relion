

#include "src/acc/cuda/cuda_benchmark_utils.h"

// Non-concurrent benchmarking tools (only for Linux)
#include <vector>
#include <algorithm>
#include <time.h>
#include <string>
#include <signal.h>

#include "src/macros.h"
#include "src/error.h"

void relion_timer::cuda_cpu_tic(const std::string& id) {
	const auto& ids = cuda_cpu_benchmark_identifiers;
	if (std::find(ids.begin(), ids.end(), id) == ids.end()) {
		cuda_cpu_benchmark_identifiers.push_back(id);
		cuda_cpu_benchmark_start_times.push_back(clock());
	} else {
		printf("DEBUG_ERROR: Provided identifier '%s' already exists in call to cuda_cpu_tic.\n", id.c_str());
		CRITICAL(ERRCTIC);
	}
}

void relion_timer::cuda_cpu_toc(const std::string& id) {
	const auto& ids = cuda_cpu_benchmark_identifiers;
	const auto it = std::find(ids.begin(), ids.end(), id);
	if (it == ids.end()) {
		printf("DEBUG_ERROR: Provided identifier '%s' not found in call to cuda_cpu_toc.\n", id.c_str());
		// exit(EXIT_FAILURE);
		return;
	}
	const size_t i = it - ids.end();
	clock_t t = clock() - cuda_cpu_benchmark_start_times[i];
	cuda_cpu_benchmark_identifiers.erase(cuda_cpu_benchmark_identifiers.begin() + i);
	cuda_cpu_benchmark_start_times.erase(cuda_cpu_benchmark_start_times.begin() + i);
	fprintf(cuda_cpu_benchmark_fPtr, "%06.2f ms ......", (float) t / CLOCKS_PER_SEC * 1000.0);
	for (const auto& x: cuda_cpu_benchmark_identifiers)
		fprintf(cuda_cpu_benchmark_fPtr, "......");
	fprintf(cuda_cpu_benchmark_fPtr, " %s\n", id.c_str());
	// printf("%s \t %.2f ms\n", id.c_str(), (float) t / CLOCKS_PER_SEC * 1000.0);
}

void relion_timer::cuda_gpu_tic(const std::string& id) {
	const auto& ids = cuda_gpu_benchmark_identifiers;
	if (std::find(ids.begin(), ids.end(), id) == ids.end()) {
		cudaEvent_t start, stop;
		cudaEventCreate(&start);
		cudaEventCreate(&stop);
		cudaEventRecord(start, 0);
		cuda_gpu_benchmark_identifiers.push_back(id);
		cuda_gpu_benchmark_start_times.push_back(start);
		cuda_gpu_benchmark_stop_times.push_back(stop);
	} else {
		printf("DEBUG_ERROR: Provided identifier '%s' already exists in call to cuda_gpu_tic.\n",
				id.c_str());
		CRITICAL(ERRGTIC);
	}
}

void relion_timer::cuda_gpu_toc(const std::string& id) {
	const auto& ids = cuda_gpu_benchmark_identifiers;
	const auto it = std::find(ids.begin(), ids.end(), id);
	if (it == ids.end()) {
		printf("DEBUG_ERROR: Provided identifier '%s' not found in call to cuda_gpu_tac.\n",
				id.c_str());
		CRITICAL(ERRGTOC);
	} else {
		const size_t i = it - ids.end();
		cudaEventRecord(cuda_gpu_benchmark_stop_times[i], 0);
		cudaEventSynchronize(cuda_gpu_benchmark_stop_times[i]);
	}
}

void relion_timer::cuda_gpu_printtictoc() {
	if (cuda_gpu_benchmark_identifiers.empty()) {
		printf("DEBUG_ERROR: There were no identifiers found in the list, on call to cuda_gpu_toc.\n");
		CRITICAL(ERRTPC);
	} else {
		for (int i = 0; i < cuda_gpu_benchmark_identifiers.size(); i++) {
			float time;
			cudaEventElapsedTime(&time,
					cuda_gpu_benchmark_start_times[i],
					cuda_gpu_benchmark_stop_times[i]);
			cudaEventDestroy(cuda_gpu_benchmark_start_times[i]);
			cudaEventDestroy(cuda_gpu_benchmark_stop_times[i]);
			fprintf(cuda_gpu_benchmark_fPtr,"%.2f ms \t %s\n",
					time, cuda_gpu_benchmark_identifiers[i].c_str());
		}
		cuda_gpu_benchmark_identifiers.clear();
		cuda_gpu_benchmark_start_times.clear();
		cuda_gpu_benchmark_stop_times.clear();
	}
}
