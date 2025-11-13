# OpenMMDihedralPlugin
An OpenMM implementation of the HPS-SS dihedral potential by Rizuan et al.


This is an openmm implementation of the HPS-SS dihedral potential by Rizuan et al. [DOI: 10.1021/acs.jcim.2c00450]  
The original code can be found here :  https://github.com/azamat-rizuan/HPS-SS-model  

## Installation notes
A conda python environment with openmm installed is required.

Tested with a Linux x86_64 machine and  


   condas gxx compiler x86_64-conda-linux-gnu-gcc [conda install -c conda-forge gxx_linux-64]  
   openmm8.1.1 and openmm8.2  

#setup compiler  
setenv CXX '`which x86_64-conda-linux-gnu-g++`'  
setenv CC '`which x86_64-conda-linux-gnu-gcc`'  

#Install  
`git clone https://github.com/karamanoslab/OpenMMDihedralPlugin`  
`mkdir build`  
`cd build`  
`cmake3 ../ -DCMAKE_CXX_FLAGS="-std=c++11 -D_GLIBCXX_USE_CXX11_ABI=1" -DCMAKE_CXX_STANDARD=11 -DCMAKE_INSTALL_PREFIX=./install`  

`cmake3 --build . `

After installation set your LD_LIBRARY_PATH

`setenv LD_LIBRARY_PATH "path to hpss_plugin"/build:${LD_LIBRARY_PATH}"`  

This builds the C++ plugin and a swig generated python wheel which is then installed with pip (> v20 is required).

Test analytical forces and regenerate the energy function as in Rizuan et al Figure 1 by running the two python scripts in tests
