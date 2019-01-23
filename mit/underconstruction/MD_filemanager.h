#ifndef MD_FILEMANAGER_H
#define MD_FILEMANAGER_H

#include <fstream>
#include <string>
#include <vector>
#include <memory>

#include "Cell.h"

class MD_filemanager{
public:
    virtual void read_infile() = 0;
    virtual std::shared_ptr<Cell> read_next() = 0;
};

class QE_MD_npt:public MD_filemanager{
public:
    QE_MD(std::string in_cel="data.cel", std::string in_pos="data.pos", std::vector<int> in_n = std::vector<int>({0}), std::vector<std::string> in_name = std::vector<std::string>({""})): cell_file(make_shared<ifstream>(in_cel)), pos_file(make_shared<ifstream>(in_pos)),atom_n(in_n),atom_name(in_name){}
private:
    std::shared_ptr<ifstream> cell_file;
    //std::shared_ptr<ifstream> wan_file;
    std::shared_ptr<ifstream> pos_file;
    
    std::vector<int> atom_n;
    std::vector<std::string> atom_name;
    
}
#endif
