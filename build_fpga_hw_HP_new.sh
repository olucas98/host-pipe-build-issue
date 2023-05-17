#!/bin/bash
export PATH=/glob/intel-python/python3/bin/:/glob/intel-python/python2/bin/:${PATH}
source /glob/development-tools/versions/oneapi/2023.1/oneapi/setvars.sh --force
rm -Rf bin_HP
cp -pR bin-original bin
icpx -fsycl -fintelfpga ./src/two_bit_v19_hostpipe_basic.cpp -Xshardware -Xsenable-unequal-tc-fusion -I/glob/development-tools/versions/oneapi/2023.1/oneapi/compiler/2023.1.0/linux/lib/oclfpga/include -I/glob/development-tools/versions/oneapi/2023.1/oneapi/compiler/2023.1.0/linux/lib/oclfpga/include/sycl/ext/intel/prototype/ -Xsboard=/opt/intel/oneapi/intel_s10sx_pac:pac_s10 -Xsprofile -o ./bin_HP/kmers_s10.fpga
