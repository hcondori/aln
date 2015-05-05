#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>

#include "helper.hpp"

#include "../sw_f32.blob"

std::string config_filename = ".aln.conf";

float* substitution_matrix()
{
  float* sm = new float[128 * 128];// = (float*) malloc (sizeof(float) * 128 * 128);

  for (int i = 0; i < 128; ++i)
    {
      for (int j = 0; j < 128; ++j)
      {
        if (i == j)
          sm[128 * i + j] = 5;
        else
          sm[128 * i + j] = -4;
      }
    }
  sm[0] = -4;
  return sm;
}

std::vector<cl::Device> all_ocl2_devices()
{
  std::ifstream infile(getenv("HOME") + std::string("/") + config_filename);
  std::string line;
  std::vector<std::string> device_names;
  while (std::getline(infile, line))
  {
    device_names.push_back(line);
    std::cout << "Leyendo" << std::endl;
  }
  
  std::vector<cl::Platform> platforms;
  std::vector<cl::Device> devices;
  
  cl::Platform::get(&platforms);  
  std::cout << "Platform number found: " << platforms.size() << std::endl;
  
  std::string platform_str;
  for(auto &platform : platforms)
  {
    std::vector<cl::Device> temp;
    platform.getDevices(CL_DEVICE_TYPE_ALL, &temp);
    platform.getInfo(CL_PLATFORM_NAME, &platform_str);
    for(auto &device : temp)
    {
      std::string device_str;
      device.getInfo(CL_DEVICE_NAME, &device_str);
      if(std::find(device_names.begin(),
                  device_names.end(), (platform_str + ": " + device_str)) != device_names.end())
      {
        devices.push_back(device);
        std::cout << "Device added: " << device_str << std::endl;
      }
    }
  }
  std::cout << "HIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII" << std::endl;
  return devices;
}

cl::Program getProgram(cl::Context& context, const char* filename)
{
  //std::ifstream sourceFile(filename);
  std::string sourceCode((const char*)sw_f32_cl);
  /*std::string sourceCode(std::istreambuf_iterator<char>(sourceFile), 
                         (std::istreambuf_iterator<char>()));*/
  cl::Program::Sources source(1, std::make_pair(sourceCode.c_str(), sourceCode.length() + 1));
  cl::Program program = cl::Program(context, source);
  return program;
}
