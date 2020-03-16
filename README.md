# asl

This repository contains the a branched version of ASL (AMPL Solver Library) that has two aims:
1. Minimize code duplication by having one copy of the files shared betweetn ASL and ASL2
2. Have the solver library thread-safe, so that multiple instances can be used concurrently

## Build instructions
### Linux systems
To build the static library under x86 Unix/Linux systems, simply do the following: 

```
mkdir build
cd build
cmake .. 
make .
```

This by default builds 64 bits versions of the libraries, to build 32 bits builds, define the variable `ARCH` when calling cmake:

```
cmake .. -DARCH=32
```

### Windows systems
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


