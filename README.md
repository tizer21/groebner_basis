## Testing
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTING=ON
cmake --build build --parallel

ctest --test-dir build --output-on-failure



## Benchmarking
cmake -S . -B build-bench \
  -DCMAKE_BUILD_TYPE=Release \
  -DGROEBNER_BASIS_BUILD_BENCHMARKS=ON \
  -DBUILD_TESTING=OFF
cmake --build build-bench --config Release --parallel

./build-bench/bench/groebner_basis_benchmarks \
  --benchmark_repetitions=5 \
  --benchmark_report_aggregates_only=true
