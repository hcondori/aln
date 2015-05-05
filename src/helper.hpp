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

#ifndef HELPER_HPP_
#define HELPER_HPP_

#include <vector>

#ifdef OCL_SUPPORT
  #ifdef __APPLE__
    #include "cl.hpp"
  #else
    #include <CL/cl.hpp>
  #endif
#endif



float* substitution_matrix();

std::vector<cl::Device> all_ocl2_devices();

/*
 * Returns de devices enabled to be used.
 * This configuration is saved in a config file
 */
void enabled_devices(cl::Platform platform);

cl::Program getProgram(cl::Context& context, const char* filename);

#endif /* HELPER_HPP_ */