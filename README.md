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
