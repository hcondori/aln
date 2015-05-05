include_directories(${ALN_SOURCE_DIR}/include)

#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -fabi-version=6 -fopenmp")

FIND_PACKAGE( Boost 1.54 COMPONENTS program_options REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

add_executable(aln aln.cpp
                   cpudispatcher.cpp
                   fastareader.cpp
                   #sequence.cpp
                   #sequencebuffer.cpp
                   sw.cpp
                   avx2_sw.cpp
                   avx2_swi.cpp
                   backtrack.cpp
                   #helper.cpp
              )
target_link_libraries(aln OpenCL ${Boost_LIBRARIES} )

#add_executable(alnconf alnconf.cpp)
#target_link_libraries(alnconf ncurses OpenCL)
