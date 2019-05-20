#include "Utility_time.h"

std::ostream& operator<<(std::ostream& os,const Utility_time::Timers& in_T){
    in_T.print(os);
    return os;
}

Utility_time::Timers GTIMER;
