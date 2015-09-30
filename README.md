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

### Performing VCI calculations

Only a single input file is required to run the LOVCI program. Example input
can be found in the tests and doc directories.

#### Direct product basis sets

Basis: Product <br>
Broadening: [type] [width] [resolution] [cutoff] <br>
Modes: [Nm] <br>
 [id] [freq] [quanta] [intensity] <br>
 ... <br>
Spectator_modes: [Ns] <br>
 [id] [freq] [intensity] <br>
 ... <br>
Force_constants: [Nfc] <br>
 [power] [modes] [value] <br>
 ...

#### Progression basis sets

Basis: Progression <br>
Prog_mode: [pmode] [pquanta] <br>
Mixed_modes: [Np] <br>
 [modes] <br>
Broadening: [type] [width] [resolution] [cutoff] <br>
Modes: [Nm] <br>
 [id] [freq] [quanta] [intensity] <br>
 ... <br>
Spectator_modes: [Ns] <br>
 [id] [freq] [intensity] <br>
 ... <br>
Force_constants: [Nfc] <br>
 [power] [modes] [value] <br>
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

Broadening: There are two types of line shapes in LOVCI (Lorentzian,Gaussian).
Both types of line shapes require a line width, the resolution for
printing the spectrum, and a maximum frequency for the final spectrum
(cutoff).

Modes: The keyword needs the total number of vibrational modes (Nm), followed
by the properties of each mode. Every mode should be labeled with an integer
ID, a frequency (cm^-1), the maximum number of quanta in the basis set, and
the intensity. The modes are labeled from 0 to (Nm-1).

Spectator_modes: The input for the (Ns) spectator modes is nearly the same as
the input for active modes. The spectators are labeled from (Nm) to (Nm+Ns+1).
Since spectator modes are not included in the basis set, there is no need to
specify the number of quanta.

Force_constants: The (Nfc) anharmonic force constants are defined by the
total power of the modes (cubic=3,quartic=4,etc), the list of integer IDs
for the modes, and the value of the force constant. The list of modes can
be given in any order and the force constants can include any positive
power. The values of the force constants (cm^-1) are the same as those
given in the output of common QM packages (e.g. Gaussian, GAMESS).

### Theory: Vibrational Configuration Interaction


