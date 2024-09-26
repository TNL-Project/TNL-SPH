cmake_minimum_required(VERSION 3.28)

# Set default CUDA architecture
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
   set(CMAKE_CUDA_ARCHITECTURES "native")
endif()

project(TNL-SPH
   LANGUAGES CXX CUDA
)

# Options to build
option( TNL-SPH_USE_SYSTEM_GTEST "Use GTest installed in the local system and do not download the latest version" ON )
option( TNL-SPH_BUILD_TESTS "Build tests" ON )
option( TNL_BUILD_BENCHMARKS "Compile the 'src/Benchmarks' directory" OFF )

# Require C++17, and disable compiler-specific extensions (if possible)
foreach(lang CXX CUDA)
    set(CMAKE_${lang}_STANDARD 17)
    set(CMAKE_${lang}_STANDARD_REQUIRED ON)
    set(CMAKE_${lang}_EXTENSIONS OFF)
endforeach()

# Set build flags for CXX
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -Wall -Werror=vla")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS_DEBUG}")

# Set build flags for CUDA
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Wall")
set(CMAKE_CUDA_FLAGS_DEBUG "-g")
set(CMAKE_CUDA_FLAGS_RELEASE "-O3 -DNDEBUG --expt-relaxed-constexpr --expt-extended-lambda -use_fast_math")
set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "${CMAKE_CUDA_FLAGS_RELEASE} ${CMAKE_CUDA_FLAGS_DEBUG}")
if(CMAKE_CUDA_COMPILER_ID STREQUAL "NVIDIA")
    set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "${CMAKE_CUDA_FLAGS_RELWITHDEBINFO} --generate-line-info")
endif()

# Warn about redundant semicolons
if( CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
    CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR
    CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang"
)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra-semi")
    if( CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR
        CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang"
    )
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra-semi-stmt")
    endif()
endif()

# Build tests using gtest
if( TNL-SPH_BUILD_TESTS )
   enable_testing()

   # repeat failed tests to tolerate sporadic failures
   set( CMAKE_CTEST_ARGUMENTS ${CMAKE_CTEST_ARGUMENTS} --repeat until-pass:2 )

   # let CTest write test results in the JUnit XML format
   set( CMAKE_CTEST_ARGUMENTS ${CMAKE_CTEST_ARGUMENTS} --output-junit "${CMAKE_BINARY_DIR}/tests-report.xml" )

   if( TNL-SPH_USE_SYSTEM_GTEST OR TNL-SPH_OFFLINE_BUILD )
      # find GoogleTest installed in the local system
      find_package( GTest REQUIRED )
   else()
      # fetch and build GoogleTest from source
      include( FetchGoogleTest )
   endif()
   set( CXX_TESTS_FLAGS )
   set( CUDA_TESTS_FLAGS )
   set( TESTS_LIBRARIES GTest::gtest_main )
   set( TESTS_LINKER_FLAGS "" )

   if( TNL_BUILD_COVERAGE AND CMAKE_BUILD_TYPE STREQUAL "Debug" )
      # set compiler flags needed for code coverage
      set( CXX_TESTS_FLAGS ${CXX_TESTS_FLAGS} --coverage )
      set( CUDA_TESTS_FLAGS ${CUDA_TESTS_FLAGS} -Xcompiler --coverage )
      set( TESTS_LINKER_FLAGS ${TESTS_LINKER_FLAGS} --coverage )
   endif()

   # TODO: Add mpi tests.
endif()





# Set the project's include path
include_directories(include)

# Find or fetch dependencies
include(FetchContent)
FetchContent_Declare(TNL
    GIT_REPOSITORY https://gitlab.com/tnl-project/tnl.git
    GIT_TAG main
    EXCLUDE_FROM_ALL
)
set(TNL_EXPORT_INTERFACE_TARGETS ON)  # needed for adding to the interface target
FetchContent_MakeAvailable(TNL)

# Set the install path for Python modules
# FIXME: for some reason relying on find_package from libs/pytnl fails (dependent target 'Python3::Module' is not defined)
find_package(Python 3 COMPONENTS Interpreter Development REQUIRED)

# TODO:
# Require python vtk

# Add project subdirectories
add_subdirectory(include)
add_subdirectory(examples)
add_subdirectory(src)
