#include <iostream>
#include <fstream>
#include <string>

using std::cout;
using std::endl;
using std::string;
using std::ifstream;
using std::streamsize;
using std::ios;

int main()
{
    
    ifstream in("/panasas/scratch/kmarcus2/out");
    string t;

    getline(in, t);
    cout << t << endl;

    getline(in, t);
    cout << t << endl;
}
