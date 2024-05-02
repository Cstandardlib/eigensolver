#include "lobpcg.h"
#include <iostream>
int main(){
    auto tp_start = get_current_time(); // tp for time point, different from duration
    for(int i=0; i<100; ++i){std::cout << "1" ;}
    std::cout << std::endl;
    auto tp_end = get_current_time();
    std::chrono::duration<double> t_total;
    t_total += tp_end - tp_start;
    std::cout << t_total.count() << "seconds" << std::endl;

}