#ifndef GLASSFIRE_UTIL_H
#define GLASSFIRE_UTIL_H
#include <iomanip>
#include <sstream>

namespace glassfire{

template<typename T>
string fmt_string(T input, bool with_plus=false, size_t precision=10){
    std::ostringstream streamObj;
    streamObj << std::fixed;
    streamObj << std::setprecision(precision);
    if (input >= 0 && with_plus) streamObj << "+";
    streamObj << input;
    return streamObj.str();
}

}
#endif