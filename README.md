# TSS-PV Benchmark
1. Run the venv:
   ``
   source ./venv/bin/activate
   ``
   - Note 1: deactivate conda if necessary: `conda deactivate`
   - Note 2: `Dalek-Curve25519` is wrapped using `Pyo3` and `Maturin` as `Curve25519_Python` package in `./venv/lib/python3.10/site-packages/curve25519_python`

2. Now it is ready to run the benchmark: `python3 ./benchmark/benchmark.py`
3. Historical versions can be found in `./history`
