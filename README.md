# MIT README #

Molecule Investigation Toolkit

## Developer ##
* Jianhang Xu: [jianhang.xu@temple.edu](mailto:jianhang.xu@temple.edu)

## Small code usage

### File convert
#### Qe_deepmd

#### Ipi_qe
Read in Ipi output file pos.xyz write to qe format .pos and .cel file
* options
> ./ipi_qe.x -p output_prefix -n input_name -t time_step -angs -f

### Read properties
#### Qe_volume
Read volume and density in QE output .cel write to given filename 
* options
> ./qe_volume.x -n input_name -o output_name -m mass_of_system -angs -f
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