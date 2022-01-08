# MHSim

MHSim: A Simulation Framework for Memristor-based Heterogeneous Computing Architectures

## Table of Contents

- [Introduction](#Introduction)
- [Setup](#Setup)
- [Usage](#usage)
- [License](#origianl-license--copyright-of-zsim)

## Introduction

MHSim is implemented with ZSim and NeuroSim. MHSim is used to evaluate the performance and energy of memristor-based accelerators (MBAs) for general-purpose applications written in C/C++.

The NeuroSim here is reorganized to support floating-point general matrix-matrix multiplication (GEMM). It can simulate the MVM results considering the non-ideal properties of memristor devices and circuits. The estimated latency and energy of performing MVM operations are also supported. We use ZSim to simulate CPUs and the memory hierarchy and exploit NeuroSim to simulate the MBA.


## Setup

**1. External Dependencies**

You need to install packages below before compiling MHSim.
- gcc (>=4.6)
- Intel Pin Toolkit (>2.8 and <3.0) (You can skip this package since the [pin](https://github.com/burymyname/pin-2.14) (2.14) has been included in this project)
- [libconfig](http://www.hyperrealm.com/libconfig)
- [libhdf5](https://github.com/HDFGroup/hdf5) (v1.8.4 path 1 or higher)
- [libelf](https://github.com/WolfgangSt/libelf)


**2. Compiling and Installation**

Update the environment script [env.sh](env.sh) according to your machine configuration.

```sh
MHSim=`pwd`
NEUROSIMPATH=$MHSim/NeuroSim/ #root path of NeuroSim
PINPATH=$MHSim/zsim/pin-2.14/  #path of pin_tool
LIBCONFIG= #path of libconfig
HDF5= #path of hdf5
LIBELF= #path of libelf
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HDF5/lib:$LIBCONFIG:$LIBELF/lib:/lib
INCLUDE=$INCLUDE:$HDF5/include:$LIBCONFIG:/include
LIBRARY_PATH=$LIBRARY_PATH:$HDF5/lib
export PINPATH NEUROSIMPATH HDF5 LD_LIBRARY_PATH LIBRARY_PATH
```

Compile MHSim with the below instructions:
```sh
$ cd MHSim
$ source env.sh
$ cd zsim
$ scons -j16
```

## Usage

After compiling the MHSim, you can evaluate any applications with: 

```sh
$ ./build/opt/zsim tests/mba.cfg
# Assume that you have specified the command in mba.cfg.
```

The simulation results can be found in zsim.out. We note that the applications without cblas_sgemm functions do not work in the MBA since the MBA mainly accelerates the MVM operations.

You can also integrate the MHSim with other application to simulate the computation error of memristor crossbar arrays. You can run the instructions below to compile a .so file:

```sh
$ cd MHSim/NeuroSim
$ make sharedobj
```
These instructions will generate a libmhsim.so file. You can link this library with other applications such as Caffe to replace the cblas_sgemm function. The main.cpp file shows examples of replacing the cblas_sgemm function and estimating the latency and energy of performing the analog MVMs.

## Origianl License & Copyright of zsim

zsim is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 2.

zsim was originally written by Daniel Sanchez at Stanford University, and per Stanford University policy, the copyright of this original code remains with Stanford (specifically, the Board of Trustees of Leland Stanford Junior University). Since then, zsim has been substantially modified and enhanced at MIT by Daniel Sanchez, Nathan Beckmann, and Harshad Kasture. zsim also incorporates contributions on main memory performance models from Krishna Malladi, Makoto Takami, and Kenta Yasufuku.

zsim was also modified and enhanced while Daniel Sanchez was an intern at Google. Google graciously agreed to share these modifications under a GPLv2 license. This code is (C) 2011 Google Inc. Files containing code developed at Google have a different license header with the correct copyright attribution.

Additionally, if you use this software in your research, we request that you reference the zsim paper ("ZSim: Fast and Accurate Microarchitectural Simulation of Thousand-Core Systems", Sanchez and Kozyrakis, ISCA-40, June 2013) as the source of the simulator in any publications that use this software, and that you send us a citation of your work.
