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
#include <thread>

#include <boost/program_options.hpp>

#ifdef OCL_SUPPORT
  #ifdef __APPLE__
    #include "cl.hpp"
  #else
    #include <CL/cl.hpp>
  #endif
#endif

//#include "helper.hpp"
#include "cpudispatcher.hpp"
#include "matrix.hpp"


namespace po = boost::program_options;

int main(int argc, char* argv[])
{
  std::string filename1;
  std::string filename2;
  std::string outfilename;
  po::options_description desc("Allowed options");
  desc.add_options()
      ("help,h", "produce help message")
      ("file1,a", po::value<std::string>()->required(), "first sequence file")
      ("file2,b", po::value<std::string>()->required(), "second sequence file")
      ("out,o", po::value<std::string>(), "output sequence file")
      ("gapopen,O", po::value<float>(), "gap open penalty")
      ("gapextend,E", po::value<float>(), "gap extend penalty")
      ("match,M", po::value<float>(), "match reward")
      ("mismatch,m", po::value<float>(), "mismatch penalty")
      ("fraction,f", po::value<float>(), "table percentage to analize")
      ("datatype,d", po::value<std::string>(), "datatype to perform the alignment")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);    

  if (vm.count("help")) {
      std::cout << desc << "\n";
      return 1;
  }
  filename1 = vm["file1"].as<std::string>();
  filename2 = vm["file2"].as<std::string>();
  outfilename = vm["out"].as<std::string>();
  
  if (vm.count("help")) {
      std::cout << desc << "\n";
      return 1;
  }
  
  float gapopen = 10.0f;
  float gapextend = 0.5f;
  float match = 5.0f;
  float mismatch = -4.0f;
  float fraction = 1.0f;
  std::string datatype = "float";
  
  if (vm.count("gapopen"))
    gapopen = vm["gapopen"].as<float>();
  
  if (vm.count("gapextend"))
    gapextend = vm["gapextend"].as<float>();
  
  if (vm.count("datatype"))
    datatype = vm["datatype"].as<std::string>();
  
  bool with_sm = true;
  bool do_banded = false;
  
  if(vm.count("match") && vm.count("mismatch"))
  {
    with_sm = false;
    match = vm["match"].as<float>();
    mismatch  = vm["mismatch"].as<float>();
  }
  
  if(vm.count("fraction"))
  {
    do_banded = true;
    fraction = vm["fraction"].as<float>();
  }
  
  
  float* sm = substitution_matrix();
  
  aln::CPUDispatcher dispatcher(filename1, filename2, outfilename);
  
  if(datatype == "float")
  {
    if(with_sm)
      dispatcher.run(gapopen, gapextend, sm, fraction);
    else
      dispatcher.run(gapopen, gapextend, match, mismatch, fraction);
  }
  else
  {
    if(with_sm)
      dispatcher.run(gapopen, gapextend, sm, fraction);
    else
      dispatcher.run((int16_t)gapopen, (int16_t)gapextend, (int16_t)match, (int16_t)mismatch, fraction);
  }
    

  std::cout << "Execution successful." << std::endl;
  delete[] sm;
}
