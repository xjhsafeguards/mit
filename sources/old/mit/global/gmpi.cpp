#include "gmpi.h"

auto GMPI = new Gmpi;

void Gmpi::init(int argc, char** argv){
#ifdef __MPI
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &NCOR);
    MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
#endif
}

void Gmpi::end(){
#ifdef __MPI
    MPI_Finalize();
#endif
}

int Gmpi::ncore(){
    return NCOR;
}

int Gmpi::rank(){
    return RANK;
}
