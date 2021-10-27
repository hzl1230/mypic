#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
// #include "../PIC_MCC/ESPIC/src/espic_info.h"
typedef double Real;

std::string unknown_cmd_info(const std::string& cmd, 
                             const std::string& file)
{
  std::ostringstream oss;
  oss << "Unknown command \"" << cmd << "\" in [" << file << "]";
  return oss.str();
}
void espic_error(const std::string& msg)
{
  std::cerr << "--> Error: " << msg << ". <--" << std::endl;
  exit(-1);
}
const Real& print(Real* arr, int size) ;

int main(int argc, char** argv)
{   
    Real arr[3]{3,4,5};
    std::vector<std::string> word{"1","2","3"};
    std::vector<Real> vec1;

    std::vector<Real> result;
    std::string cmd{"property"}, infile{"csection.in"} ;
    // result.resize(word.size());
    // espic_error(unknown_cmd_info(cmd, infile));
    transform(word.begin(), word.end(), back_inserter(result), [](auto& str){ return (Real)atof(str.c_str()); });

    for_each(result.begin(), result.end(),[](Real& str){std::cout << str << std::endl;});
    size_t inn(result.size());
    std::cout << inn << std::endl;
    vec1.resize(inn);
    for(int i = 0; i < 4; i++){ 
        vec1[i] = i;
        std::cout << vec1[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}

const Real& print(const Real* & src, Real* target, int size)
{
    for (int i = 0; i < size; i++) {
        std::cout << *(src+i) << " ";
        target[i] = *(src+i);
    }
    std::cout << std::endl;
}

