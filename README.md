# echo_state_network
An echo state network (a recurrent neural network) written in fortran for the gnu fortran compiler. Currently the network is trained reproduce chaotic motion about a Lorenz attractor.

Program makes use of LAPACK function "dgeev" to normalise the spectral radius of the ESN. 
