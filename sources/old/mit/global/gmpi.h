#ifndef GMPI_H
#define GMPI_H

class Gmpi{
    int NCOR = 1;
    int RANK = 0;
public:
    Gmpi(){}
    
    void init(int argc, char** argv);
    void end();
    
    int ncore();
    int rank();
    
};

extern Gmpi *GMPI;

#endif
