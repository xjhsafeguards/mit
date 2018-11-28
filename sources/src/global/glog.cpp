#include "glog.h"

auto GLOG = new Glog;

void Glog::init(){
    start_time = clock();
    ofs_log.open("Core" + to_string(GMPI->rank()) + ".log");
    assert(ofs_log.good());
    ofs_log << " ------------------------------ " << endl;
    ofs_log << "     Initiate MIT Version" << MIT_VERSION_MAJOR << "." << MIT_VERSION_MINOR << "   " << endl;
    ofs_log << " ------------------------------ " << endl;
    time_t time_start = std::time(NULL);
    ofs_log << "Start Time: " << ctime(&time_start) << endl;
}

void Glog::end(){
    if(GMPI->rank() == 0)
        cout << "Job Done!" << endl;
    assert(ofs_log.good());
    ofs_log << " ------------------------------ " << endl;
    ofs_log << "            Finished            " << endl;
    ofs_log << " ------------------------------ " << endl;
    time_t time_end = std::time(NULL);
    ofs_log << "End Time: " << ctime(&time_end) << endl;
    end_time = clock();
    ofs_log << "Totle Running Time : " << static_cast<double>(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;
}

ofstream& Glog::setprecision(size_t p){
    assert(ofs_log.good());
    ofs_log.precision(p);
    return ofs_log;
}

ofstream& Glog::log_setprecision(size_t p){
    assert(ofs_log.good());
    ofs_log.precision(p);
    return ofs_log;
}

ostream& Glog::out_setprecision(size_t p){
    std::cout.precision(p);
    return std::cout;
}

Glog* Glog::write(const string& s){
    assert(ofs_log.good());
    ofs_log << s;
    return this;
}

Glog* Glog::write_line(const string& s){
    assert(ofs_log.good());
    ofs_log << s << endl;
    return this;
}

ofstream& Glog::stream(){
    assert(ofs_log.good());
    return ofs_log;
}

ofstream& Glog::log_stream(){
    assert(ofs_log.good());
    return ofs_log;
}

ostream& Glog::out_stream(){
    return std::cout;
}
