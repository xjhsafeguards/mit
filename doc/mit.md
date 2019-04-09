# MIT folder #

This setup used Eigen for data structures of cells and provide some functions to analyze the system

## Class Cell functions ##

#### constructor ####
* Cell()
#### operators ####
* bool operator==(const Cell& cel) const 
#### I/O ####
* virtual void read(std::istream&)
* virtual void skip(std::istream&) const 
* virtual void read(std::istream&,std::istream&)
* virtual void skip(std::istream&,std::istream&) const
* virtual void read_cell(std::istream&)
* virtual void skip_cell(std::istream&) const
* virtual void read_pos(std::istream&)
* virtual void skip_pos(std::istream&) const 

* virtual void set_box(double a, double b, double c)

* virtual std::ostream& write_cell(std::ostream& os) const 
* virtual std::ostream& write_positions(std::ostream& os) const 
* virtual std::ostream& write_position(std::ostream& os, int i) const
#### interaction with other class ####
* virtual Atom atom(int index) const
* virtual Atom atom(int index) const
* virtual Atom_group atoms() const
* virtual Atom_group atoms(std::string in_type) const
* virtual Atom_group atoms(int start_index,int end_index) const
* virtual Box  box() const
#### print values ####
* virtual double volume() const
* virtual double snapshot() const 
* virtual double time() const
* virtual std::string type(int) const
* virtual double distance(int i,int j) const
* virtual double angle(int i,int j, int k) const

## Class Atom functions ##
#### constructor
* Atom(const Cell& in_cel,int in_index):cel(in_cel),index(in_index)
* Atom(const Atom&)
* Atom(Atom&&)

#### print values ####
Using functions from Cell class
* double distance(int j) const
* double distance(const Atom& a1) const
* double angle(int j, int k) const
* double angle(const Atom& a1, const Atom& a2) const
* std::string type() const

## Class Atom_group ## 
#### constructor
* Atom_group(const Cell& in_cel,std::vector<int> in_indexs=std::vector<int>())
#### operations 
* void add(const Atom& a)
* Atom_iterator begin()
* Atom_iterator end()
### class Atom_iterator ###
#### constructor
* Atom_iterator(const Atom_group& in_ag,std::vector<int>::iterator in_it)
* Atom operator*() 
* Atom_iterator operator++() 
* bool operator!=(const Atom_iterator& in_AI)



