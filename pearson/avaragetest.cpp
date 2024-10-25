#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>


double getnumber(const std::string& input)
{
    double number;
    std::string temp = input.substr(12,4);
    std::istringstream stream(temp);
    //stream >> number;
    //std::cout << "temp: " << temp << " stream >> number: " << number << std::endl;
    stream >> number;
    return number;

}

int main()
{
    double sum;
    std::vector<double> median;
    for(int i = 2; i < 13; i++)
    {
        std::ifstream file("outputfiles/1024gprof-" + std::to_string(i) + ".txt");
        if(!file.is_open())
            return 1;

        std::string line;
        int line_count = 0;
        while(line_count < 10 && std::getline(file,line))
        {
            if(line_count == 9)
            {
                sum += getnumber(line);
                median.push_back(getnumber(line));
                std::cout << median.back() << std::endl;
            }
            line_count++;
        }
        file.close();
    }
    double temp = median[0];
    int count = 0, i = 0, j = 0;
    while(true)
    {
        if(temp <= median[i])
            count++;
        if(i == median.size())
            if(count == 4)
                break;
            else
            {
                i = 0;
                count = 0;
                temp = median[j++];
            }
        i++;
    }
    std::cout << "avarage runtime: " << sum/10 << std::endl;
    std::cout << "median runtime: " << temp << std::endl;
    
    return 0;
}