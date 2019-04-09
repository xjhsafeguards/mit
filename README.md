# MIT README #

Molecule Investigation Toolkit

## Developer ##
* Jianhang Xu: [jianhang.xu@temple.edu](mailto:jianhang.xu@temple.edu)

## Class set up for mit ##
[Readme for mit folder](doc/mit.md)

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

#### Ipi_temperature
Read temperature in ipi data.out and write temperature and average to given filename
* options
> ./ipi_temperature.x -n input_name -o output_name -c temperature_col -t time_col -s step_col

#### select_thd
Read in temperature.txt, average_hbonds.txt and volume.txt and write the deviation from average to given filename
* options
> ./select_thd.x -t temperature_name -h hbond_name -d density_name -o output_name 

## Install ##

## Classes in mit/mit ##



## Classes in mit cmake version ##

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
