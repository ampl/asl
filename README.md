# asl

This repository contains the up-to-date version of ASL (AMPL Solver Library) as maintained by David Gay. It supports ```cmake``` multiplatform builds.
There are a few cmake switches that define which optional modules will be built:

- *BUILD_EXAMPLES* builds the ASL examples
- *BUILD_F2c* build the f2c (Fortran To C) library; implied in case the examples are being compiled
- *BUILD_MT_LIBS* builds the asl-mt and asl2-mt libraries, multithreaded, compiled using OpenMP switches

## Linux systems
To build the static library under x86 Unix/Linux systems and the examples, simply do the following: 

```
mkdir build
cd build
cmake .. -DBUILD_EXAMPLES=1
make .
```

This by default builds 64 bits versions of the libraries, to build 32 bits builds, define the variable `ARCH` when calling cmake:

```
cmake .. -DARCH=32
```

## Windows systems
To build the library on Windows (this assumes Visual Studio 2019 as a generator); note also that after the third step you'll have a file `ASL.sln` that can 
be opened from Visual Studio to continue the build from there.

```
md build
cd build
cmake .. 
cmake --build .
```

To build other flavours of the library, replace the ```cmake ..``` command with:

* **VS2019 32 bits:** `cmake .. -A Win32`
* **VS2017 64 bits:** `cmake .. -G "Visual Studio 15 2017 Win64"
* **VS2017 32 bits:** `cmake .. -G "Visual Studio 15 2017"


## License

BSD-3

***
Copyright © 2023 AMPL Optimization inc. All rights reserved.

