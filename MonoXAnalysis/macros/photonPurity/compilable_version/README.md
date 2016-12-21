# photonPurityCode
code to run photon purity for monojet analysis

There are two ways to use it at the moment:
1) using root prompt as for any macro
2) using as a C++ executable after compiling with make command (there is a Makefile in this directory)


1) [prompt]$ root -l -q 'MakeFitPurityNewId.cc+("MEDIUM")'

where "MEDIUM" is the photon ID. Livia suggested using MEDIUM for now (this is one of the IDs accepted by the macro at the moment)

2) [prompt]$ make
   [prompt]$ ./tmp/MakeFitPurityNewId MEDIUM

The Makefile will compile each .C source in the current directory and use them when linking. Therefore, unless you change this behaviour, make sure there are no other sources you don't want to be linked to executable.

Objects files (.o) will be store in tmp/objects directory created by Makefile in the current one. The executable is in tmp/.



==========================
--> WARNING
==========================
At the moment, the code compiles but it crashes during execution. Work in progress ...