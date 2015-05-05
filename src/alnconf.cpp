/*
 Copyright (C) 2015 Héctor Condori Alagón.

 This file is part of ALN, the massive Smith-Waterman pairwise sequence aligner.

 ALN is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 ALN is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with ALN.  If not, see <http://www.gnu.org/licenses/>.
 */

//#define __CL_ENABLE_EXCEPTIONS

#include <iostream>
#include <vector>
#include <CL/cl.hpp>

std::string RED_BEGIN = "\033[1;31m";
std::string RED_END = "\033[0m";

std::string enabled_string =  "(ENABLED)  ";
std::string disabled_string = "(DISABLED) ";

std::string colored_string(std::string str)
{
  return RED_BEGIN + str + RED_END;
}

int main(int argc, char* argv[])
{
  std::vector<cl::Device> all_devices;
  std::vector<cl::Device> enabled_devices;

  std::vector<cl::Platform> platforms;
  cl::Platform::get(&platforms);

  std::string str;
  std::string command;

  while(true)
  {
    std::cout << colored_string("Platforms available:") << std::endl;
    for(int i = 0; i < platforms.size(); i++)
    {
      platforms[i].getInfo(CL_PLATFORM_NAME, &str);
      std::cout << " " << i + 1 << ": " <<str << std::endl;
      std::vector<cl::Device> devices;
      platforms[i].getDevices(CL_DEVICE_TYPE_ALL, &devices);
      for(int j = 0; j < devices.size(); j++)
      {
        devices[j].getInfo(CL_DEVICE_NAME, &str);
        std::cout << "  " << disabled_string << i + 1 << "." << j+1 << ": " << str;
        devices[j].getInfo(CL_DEVICE_VENDOR, &str);
        std::cout << ", " << str << std::endl;
      }
    }
    std::cout << "command: ";
    std::cin >> command;
  }
}
