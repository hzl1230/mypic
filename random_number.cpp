#include <iostream>
#include <random>
#include <ctime>
#include <vector>
#include <iterator>
#include <algorithm>
#include <string>
#include <fstream>
using namespace std;
void fprint(vector<double>, ofstream& );
int main(int argc, char *argv[])
{ 
    uniform_int_distribution<> u(0,9);
    default_random_engine e((size_t)time(NULL));
    uniform_real_distribution<double> r(0,1); 
    double DEinv = 0.1;
    const string str = "Hello";
    vector<double> vec(10);
    for (size_t i = 0; i < 10; i++){
        vec[i] = r(e);
        cout << vec[i] << " ";
    }
    cout << endl;
    // for_each(vec.begin(), vec.end(), [=](auto& x){ x *= DEinv; });
    // for_each(vec.begin(), vec.end(), [=](auto x){ cout << x << " "; });
    ofstream ff("a.dat");
    fprint(vec, ff);
    
    return 0;
}

void fprint(vector<double> vec, ofstream& ff){
    for(auto item : vec)
        ff << item << " ";
    ff << endl;
}