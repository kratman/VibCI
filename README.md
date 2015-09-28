[//]: # (Mixture of GitHub markdown and HTML. HTML is needed for formatting.)

***
<div align=center> <h2>
LOVCI: Ladder Operator Vibrational Configuration Interaction
</h2> </div>

<div align=center> <h4> By: Eric G. Kratz </h4> </div>

***

<h4>
NOTICE: This repository is a work in progress. Currently only the harmonic
calculations are functioning properly.
</h4>

### Introduction

LOVCI is a simple program to calculate anharmonic vibrational spectra.

### Installation

Currently, the binary and user's manual are not included in the repository.
The Makefile can be used to generate both files. Since LOVCI is designed to
be simple, only a small number of packages are required to compile the code.
An approximate list of packages is given below.
```
 LOVCI binary: OpenMP, Eigen3
```

To install LOVCI, clone the git repository
```
user:$ mkdir VibCI
user:$ git clone https://github.com/kratman/VibCI.git ./VibCI/
```

or unpack the zipped source code
```
user:$ mkdir VibCI
user:$ cd VibCI/
user:$ unzip VibCI-master.zip
user:$ mv VibCI-master/* .
user:$ rmdir VibCI-master
```

The source code can be compiled with the Makefile provided with LOVCI.
On Ubuntu boxes, the Makefile should function without modifications. However,
with other operating systems it may be necessary to change the path to the
Eigen3 package.
```
Default: -I/usr/include/eigen3/
```

The Makefile can produce both the documentation and the binary.
```
user:$ make install
```

### Running a calculation



### Theory: Vibrational Configuration Interaction


