#ifndef UTILITY_TIME_H
#define UTILITY_TIME_H

#include <ctime>
#include <map>
#include <string>
#include <iostream>
#include <iomanip>

/*
 The Timer class
  with function reset(),start(),stop(),get_total()
 The Timers class
  with function reset(str),start(str),stop(str)
  print(os),summerize(os) // show all timer
  show_time(os,note) // current time
 */

namespace Utility_time {

    class Timer{
        bool is_on = false;
        clock_t begin_t;
        clock_t total = 0;
    public:
        // timer controls
        void reset();
        void start();
        void stop();
        // return the total timer count in second
        double get_total() const;
    };
    
    class Timers: public std::map<std::string,Timer>{
        friend std::ostream& operator<<(std::ostream&,const Timers&);
        
    public:
        typedef std::map<std::string,Timer> super;
        using super::super;
        
        // timer controls
        void reset(std::string);
        void start(std::string);
        void stop(std::string);
        // summerize the result to given ostream
        void print(std::ostream&) const;
        void summerize(std::ostream&) const;
        // show current time
        void show_time(std::ostream&,std::string) const;
    };
    
    
    inline void Timer::reset(){
        is_on = false;
        total = 0;
    }
    inline void Timer::start(){
        if(!is_on){
            is_on = true;
            begin_t = clock();
        }
    }
    inline void Timer::stop(){
        if(is_on){
            is_on = false;
            total += clock() - begin_t;
        }
    }
    inline double Timer::get_total() const{
        return (static_cast<double>(total)/CLOCKS_PER_SEC);
    }
    
    
    inline void Timers::reset(std::string name){
        this->at(name).reset();
    }
    inline void Timers::start(std::string name){
        (*this)[name].start();
    }
    inline void Timers::stop(std::string name){
        (*this)[name].stop();
    }
    inline void Timers::print(std::ostream& os) const{
        os << "----------------------------------" << "\n";
        os << "Timer summerize:\n";
        os << std::setw(15) << "Name" << std::setw(15) << "Time(s)"<< std::endl;
        for( const auto& v : *this ){
            os << std::setw(15) << v.first << std::setw(15) << v.second.get_total() << std::endl;
        }
    }
    inline void Timers::summerize(std::ostream& os) const{
        print(os);
        show_time(os,"End Time");
    }
    inline void Timers::show_time(std::ostream& os,std::string in_string = "Current Time") const{
        time_t end_t=std::time(NULL);
        os << in_string << ": " << ctime(&end_t) << std::endl;
    }
}

extern Utility_time::Timers GTIMER;

using namespace Utility_time;

#endif
