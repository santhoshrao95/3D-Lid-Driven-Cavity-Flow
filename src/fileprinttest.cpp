#include <iostream>
#include<fstream>
using namespace std;

int main()
{
    for (size_t i = 0; i < 10; i++)
    {
        std::ofstream output("whereiam", std::ios::app);
        output<<"I have entered count: " << i << std::endl;
        output.close();
    }
    return 0;
}