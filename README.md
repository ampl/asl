# asl

This repository contains the up-to-date version of ASL (AMPL Solver Library) as maintained by David Gay. It supports ```cmake``` multiplatform builds.

## Linux systems
To build the static library under x86 Unix/Linux systems, simply to the following (note the definition of ARCH to 64 to build a 64 bits version of the library. Set it to 32 for 32 bits).

```
mkdir build
cd build
cmake .. -DARCH=64
make .
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


