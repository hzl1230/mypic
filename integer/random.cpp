#include <random>
#include <vector>
#include <iostream>
#include <algorithm>
// #include <numeric>
#include <iterator>

using std::string;
using std::vector;
typedef double Real;

int main(int argc, char **argv)
{
    vector<Real> arr{1.1, 2.3, 3.5, 4.7, 5.2, 6.1, 7.0}, result;
    vector<Real>::const_iterator beg = arr.cbegin(), end = arr.cend();
    // size_t num = end-beg;
    // std::cout << num << std::endl;
    vector<string> words{"one", "two", "three", "four", "five", "six"};
    srand((unsigned int)time(NULL));
    std::random_device rd; // 产生一个真随机数作为随机数生成器的种子
    std::mt19937 g(rd());
    std::uniform_real_distribution<> RDr(0., 10.);
    std::normal_distribution<> RMB(0., sqrt(2));
    vector<Real> x, y, vx, vy, z(6, 0.), vz(6, 0.);
    for (int i = 0; i < 6; ++i)
    {
        x.emplace_back(RDr(g));
        y.push_back(RDr(g));
        vx.push_back(RMB(g));
        vy.push_back(RMB(g));
    }
    std::copy(x.begin(), x.end(), std::ostream_iterator<Real>(std::cout, ","));
    std::cout << std::endl;
    std::copy(y.begin(), y.end(), std::ostream_iterator<Real>(std::cout, ","));
    std::cout << std::endl;
    std::copy(z.begin(), z.end(), std::ostream_iterator<Real>(std::cout, ","));
    std::cout << std::endl;
    std::copy(vx.begin(), vx.end(), std::ostream_iterator<Real>(std::cout, ","));
    std::cout << std::endl;
    std::copy(vy.begin(), vy.end(), std::ostream_iterator<Real>(std::cout, ","));
    std::cout << std::endl;
    std::copy(vz.begin(), vz.end(), std::ostream_iterator<Real>(std::cout, ","));
    std::cout << std::endl;

    // std::shuffle(arr.begin(), arr.end(), g);
    // std::shuffle(words.begin(), words.end(), g);
    // std::copy(arr.begin(), arr.end(), std::ostream_iterator<Real>(std::cout, " "));
    // std::copy(words.begin(), words.end(), std::ostream_iterator<string>(std::cout, " "));

    // std::cout << std::endl;
    // std::copy(arr.begin(), arr.begin()+3, back_inserter(result));
    // std::copy(result.begin(), result.end(), std::ostream_iterator<Real>(std::cout, " "));
    // std::cout << std::endl;
}