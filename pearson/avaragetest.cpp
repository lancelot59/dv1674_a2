#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


double getnumber(const std::string& input)
{
    std::istringstream stream(input);
    double number;
    stream >> number;
    if(stream >> number)
        return number;

}

int main()
{
    double sum;
    for(int i = 1; i < 11; i++)
    {
        std::ifstream file("outputfiles/1024gprof-" + std::to_string(i) + ".txt");
        if(!file.is_open())
            return 1;

        std::string line;
        int line_count = 0;
        while(line_count < 10 && std::getline(file,line))
        {
            sum += getnumber(line);
        }
        file.close();
    }
    std::cout << "avarage runtime: " << sum/10 << std::endl;
    return 0;
}