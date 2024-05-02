#include <iostream>
#include <chrono>


/**
 * wrap for std::chrono::high_resolution_clock::now();
 */
std::chrono::time_point<std::chrono::high_resolution_clock> get_current_time() {
    return std::chrono::high_resolution_clock::now();
}
int main() {
    // 开始计时
    auto start = std::chrono::high_resolution_clock::now();

    // 在这里放置您要计时的代码
    // ...
    std::cout << "Starting"<< std::endl;

    // 结束计时
    auto end = std::chrono::high_resolution_clock::now();

    // 计算执行时间
    std::chrono::duration<double> duration = end - start;
    std::cout << "time: " << duration.count() << " seconds" << std::endl;

    auto start2 = get_current_time();

    // 在这里放置您要计时的代码
    // ...
    std::cout << "2nd Starting"<< std::endl;

    // 结束计时
    auto end2 = get_current_time();

    // 计算执行时间
    std::chrono::duration<double> duration2 = end2 - start2;
    std::cout << "2nd time: " << duration2.count() << " seconds" << std::endl;
    auto start3 = get_current_time();
    for(int i=0; i<1000; ++i);
    auto end3 = get_current_time();
    auto duration3 = std::chrono::duration_cast<std::chrono::milliseconds>(end3 - start3);
    std::cout << "3rd time: " << duration3.count() << " seconds" << std::endl;

    return 0;
}