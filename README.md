[//]: # (Mixture of GitHub markdown and HTML. HTML is needed for formatting.)

***
<div align=center> <h2>
LOVCI: Ladder Operator Vibrational Configuration Interaction
</h2> </div>

<div align=center> <h4> By: Eric G. Kratz </h4> </div>

***

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

Only a single input file is required to run the LOVCI program.

#### Direct product basis sets

Basis: Product
Broadening: [type] [width] [resolution] [cutoff]
Modes: [Nm]
 [id] [freq] [quanta] [intensity]
 ...
Spectator_modes: [Ns]
 [id] [freq] [intensity]
 ...
Force_constants: [Nfc]
 [power] [modes] [value]
 ...

#### Progression basis sets

Basis: Progression
Prog_mode: [pmode] [pquanta]
Mixed_modes: [Np]
 [modes]
Broadening: [type] [width] [resolution] [cutoff]
Modes: [Nm]
 [id] [freq] [quanta] [intensity]
 ...
Spectator_modes: [Ns]
 [id] [freq] [intensity]
 ...
Force_constants: [Nfc]
 [power] [modes] [value]
 ...

#### Keywords

Basis: There are two options for the basis keyword (product,progression).
A product basis is a large complete basis. Progression basis sets are
small and assume that there are 1D modes combined with a Franck-Condon
progression with a low frequency mode.

Prog_mode: The prog_mode keyword needs two values (pmode,pquanta). The value
of pmode is the integer ID of the low-frequency mode used in the progression.
The value of pquanta is the number of quanta to include in the progressions.

Mixed_modes: The mixed_modes keyword takes the number of modes that have
a progression (Np) followed by a list of integer IDs for the modes.

Broadening: 


Modes: 


Spectator_modes: 


Force_constants: 


### Theory: Vibrational Configuration Interaction


