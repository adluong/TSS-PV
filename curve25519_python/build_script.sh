#!/bin/bash

# Activate the virtual environment
# source ../venv/bin/activate

cargo clean

# Build the project with RUSTFLAGS for SIMD optimization
RUSTFLAGS='--cfg curve25519_dalek_backend="simd" -C target_feature=+avx2' maturin build --release

# Find the generated wheel file (assuming one .whl file is produced in target/wheels)
WHEEL_FILE=$(ls ./target/wheels/*.whl | head -n 1)

# Install the wheel if found
if [ -f "$WHEEL_FILE" ]; then
    pip install "$WHEEL_FILE"
else
    echo "Error: No wheel file found in ./target/wheels"
    exit 1
fi
