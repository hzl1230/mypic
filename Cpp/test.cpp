#include <iostream>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <array>
#include <algorithm>
#include <numeric>
#include <string>
#include <vector>
#include <map>

using namespace std;
// int sum(int ){

// }

int main(int argc, char* argv[]){
    // vector<int> a(5), b(5), c(5), d(5);
    enum AAA {a, b, c};
    array<int,10UL> arr;
    const int TYPE1 = 0, TYPE2 = 1;
    const int TYPETOT = 2;
    map<AAA, vector<int> > mt;
    array<int,10UL> res[TYPETOT];
    // fill(arr.begin(),arr.end(),0.0);
    // generate(a.begin(), a.end()+1, [i=1]()mutable{ return (i++); });
    // int ld = 5;
    mt[AAA::a].resize(10);
    // mt[AAA::b]
    generate(arr.begin(), arr.end(), [i=1]()mutable{ return (i++); });
    transform(arr.begin(), arr.end(), mt[AAA::a].begin(), [](int x) ->int{if(x>5) return x; else return 5; });
    for(auto iter = mt.begin(); iter != mt.end(); iter++){
        iter->second.resize(5);
        cout << iter->first << " " << iter->second.size() << endl;
    }
    // for_each(mt[AAA::a].begin(),mt[AAA::a].end(),[&](int x){ cout << x << ","; });
    // for(int i=0; i<10; i++) { cout << mt[AAA::a][i] << endl; }
    return 0;
}