# MIT README #

Molecule Investigation Toolkit

## Developer ##
* Jianhang Xu: [jianhang.xu@temple.edu](mailto:jianhang.xu@temple.edu)

## Small code usage
* 

## Install ##

## Classes in mit ##

### Global
#### Glog:
* Glog(){}
* void init();
* void end();
//set
* ofstream& setprecision(size_t p=10);
* ofstream& log_setprecision(size_t p=10);
* ostream& out_setprecision(size_t p=10);
//write
* Glog* write(const string& s);
* Glog* write_line(const string& s);
* ofstream& stream();
* ofstream& log_stream();
#### Gmpi:
* Gmpi(){}
* void init(int argc, char** argv);
* void end();
//read
* int ncore();
* int rank();

### IO
#### Input:
* Input(){}
* void init();
* void iarg(int argc,char** argv);

### Data
#### Dlist:
in progress
