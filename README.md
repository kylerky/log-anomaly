# Log Anomaly

This is intended for online/offline detection of anomalies using logs.

Log Anomaly makes use of systemd to read logs in the system, so make sure you are running it on a systemd-enabled system. In addition, in order to compile it from source, you will need to have Eigen3 library installed.

## how to compile

make sure these are installed:

- cmake (for generating build files)

- a build system (for example, make)

- systemd (should be run on a system that is booted by systemd)

- Eigen3

- openmp (optional)

- a c++ compiler that supports C++17

then run cmake like this:

```
cd /path/to/build/dir
cmake /path/to/source
make
```

that's it.

BTW, no install scripts.

A few executables will be generated. Two lie in the 'detection' directory under the build path, which are detector and viewer respectively. Others lie in 'graph' and 'journal'. They are for testing.

## References

1. Elisabeth Baseman, Sean Blanchard, Zongze Li and Song Fu, [Relational Synthesis of Text and Numeric Data for Anomaly Detection on Computing System Logs](https://ieeexplore.ieee.org/document/7838262)

2. Xiaoyu Wang and Osamu Hasegawa, [Adaptive density estimation based on self-organizing incremental neural network using Gaussian process](https://ieeexplore.ieee.org/abstract/document/7966401)

3. Yoshihiro Nakamura and Osamu Hasegawa, [Nonparametric Density Estimation Based on Self-Organizing Incremental Neural Network for Large Noisy Data](https://ieeexplore.ieee.org/document/7389432)
