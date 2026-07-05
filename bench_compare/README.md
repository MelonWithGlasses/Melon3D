# Reproducing the Melon3D vs box3d comparison

`box3d_scenes.cpp` mirrors the scenes of `../bench/bench.cpp` on top of the
[box3d](https://github.com/erincatto/box3d) API, so both engines simulate
identical geometry with identical timesteps and step counts.

```sh
# 1. get and build box3d (tested at commit 1bec63c)
git clone https://github.com/erincatto/box3d extern/box3d
cmake -S extern/box3d -B extern/box3d/build -DCMAKE_BUILD_TYPE=Release \
      -DBOX3D_SAMPLES=OFF -DBOX3D_UNIT_TESTS=OFF -DBOX3D_BENCHMARKS=OFF \
      -DCMAKE_C_FLAGS="-O3 -mavx2 -mfma" -DCMAKE_CXX_FLAGS="-O3 -mavx2 -mfma"
cmake --build extern/box3d/build -j

# 2. build the harness
g++ -O3 -mavx2 -mfma -std=c++17 -I extern/box3d/include box3d_scenes.cpp \
    -L extern/box3d/build/src -lbox3d -lm -o box3d_scenes

# 3. run both
./box3d_scenes 4          # box3d, 4 substeps
../build/melon3d_bench    # Melon3D
```

Methodology notes (please keep these in mind when quoting numbers):

- box3d runs **single-threaded** here because this harness does not wire up
  its task-system callbacks. box3d fully supports multithreading when a
  task scheduler is provided; the single-thread configuration simply
  reflects its out-of-the-box behavior in this harness.
- Melon3D numbers use its built-in thread pool (results are bit-identical
  for any worker count).
- Both engines are timed the same way: wall time around the full step call
  (Melon3D's internal profiler reports exactly that). The churn scene also
  includes body create/destroy time in both harnesses.
- Every scene uses the same 60 Hz timestep, 4 substeps, identical body
  counts, sizes, densities, friction values, and the same PRNG sequence
  where randomness is involved.
