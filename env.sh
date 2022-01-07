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
