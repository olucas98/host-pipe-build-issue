# host pipe build issue


I created a version of the k-mers counting app with most of the computation removed, leaving just the pipeline structure. Both the full design and this simplified version can be compiled for emulation and run without issue, only having problems when trying to build for hardware.

Attempting to build produces this error: "host_pipes.hpp:369: Compiler Error: Compiler was unable to determine any resolvable address spaces. Please check if an appropriate board resource was specified."
 
Included files:

•	two_bit_v19_hostpipe.cpp – source code, basic pipeline structure is Phase 1: Host -> 1A -> 1B -> 1C, then Phase 2: Host -> 2A -> 2B -> 2C -> Host (3 total host pipes)

•	build_fpga_hw_HP_new.sh – script to build hardware on DevCloud, submitted with “! qsub -l nodes=1:fpga_compile:ppn=2 -l walltime=24:00:00 -d . build_fpga_hw_HP_new.sh”

•	build_fpga_hw_HP_new.sh.e2287893– error file produced when running the build script, this example is from building the smaller version that’s attached, but full app produced similar results

•	fasta.c & fasta.h – helper files needed to compile the project
