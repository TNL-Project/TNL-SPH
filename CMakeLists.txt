cmake_minimum_required(VERSION 3.28)

# Set default CUDA architecture
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
   set(CMAKE_CUDA_ARCHITECTURES "native")
endif()

project(TNL-SPH
   LANGUAGES CXX CUDA
)

# Declare all CMake options for the project
option(TNL_SPH_EXPORT_INTERFACE_TARGETS "Instruct CMake to generate rules to export the interface target" ${PROJECT_IS_TOP_LEVEL})
option(TNL_SPH_ENABLE_MPI "Enable MPI support, if available" ON)

# Require C++20, and disable compiler-specific extensions (if possible)
foreach(lang CXX CUDA)
    set(CMAKE_${lang}_STANDARD 20)
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
set(CMAKE_CUDA_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "${CMAKE_CUDA_FLAGS_RELEASE} ${CMAKE_CUDA_FLAGS_DEBUG}")

# Enforce (more or less) warning-free builds
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wno-error=deprecated -Wno-error=deprecated-declarations -Wno-error=uninitialized")

if(CMAKE_CUDA_COMPILER_ID STREQUAL "NVIDIA")
    set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "${CMAKE_CUDA_FLAGS_RELWITHDEBINFO} --generate-line-info")
    if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC")
        # Enforce (more or less) warning-free builds for host code
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler -Werror -Xcompiler -Wno-error=deprecated -Xcompiler -Wno-error=deprecated-declarations")
    endif()
    # Force colorized output for nvcc's host compiler (assuming it is GCC...)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler -fdiagnostics-color")
endif()

if(CMAKE_CUDA_COMPILER_ID STREQUAL "Clang")
    # Enforce (more or less) warning-free builds
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Werror -Wno-error=deprecated -Wno-error=deprecated-declarations -Wno-error=unknown-cuda-version")
    # workaround for Clang 15
    # https://github.com/llvm/llvm-project/issues/58491
    set(CMAKE_CUDA_FLAGS_DEBUG "-g -Xarch_device -g0")
endif()

# Disable GCC's infamous warnings (they produce many false positives)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-maybe-uninitialized -Wno-stringop-overflow")
    if(CMAKE_CUDA_COMPILER_ID STREQUAL "NVIDIA")
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler -Wno-maybe-uninitialized -Wno-stringop-overflow")
    endif()
endif()

# Warn about redundant semicolons
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
   CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR
   CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" OR
   CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra-semi")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR
       CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" OR
       CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra-semi-stmt")
    endif()
endif()

set(TNL_SPH_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/include" CACHE INTERNAL "Directories where TNL-SPH headers are located")

# Create the exported targets
add_library(TNL_SPH INTERFACE)
# Create aliases to match exported targets
add_library(TNL::TNL_SPH ALIAS TNL_SPH)

# Set the project's include path
include_directories(include)

# Find or fetch dependencies
include(FetchContent)
FetchContent_Declare(TNL
    GIT_REPOSITORY https://gitlab.com/tnl-project/tnl.git
    GIT_TAG TH/particles
    EXCLUDE_FROM_ALL
)
set(TNL_EXPORT_INTERFACE_TARGETS ON)  # needed for adding to the interface target
FetchContent_MakeAvailable(TNL)

# FIXME: We should also require python with vtk module

# Download the gencase from DualSPHysics
option(ENABLE_GENCASE_DOWNLOAD "Enable downloading GenCase_linux64" FALSE)
if(ENABLE_GENCAGE_DOWNLOAD)
   set(DOWNLOAD_URL "https://github.com/DualSPHysics/DualSPHysics/raw/master/bin/linux/GenCase_linux64")
   set(DESTINATION_DIR "${CMAKE_SOURCE_DIR}/src/tools/3rdparty")
   set(DESTINATION_FILE "${DESTINATION_DIR}/GenCase_linux64")

   # Create the destination directory if it doesn't exist
   file(MAKE_DIRECTORY ${DESTINATION_DIR})
   # Download the file during the configuration phase
   file(DOWNLOAD
       ${DOWNLOAD_URL}
       ${DESTINATION_FILE}
       SHOW_PROGRESS
       STATUS DOWNLOAD_STATUS
   )
   # Set executable permissions for the downloaded file
   file(CHMOD ${DESTINATION_FILE} PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
else()
    message(STATUS "GenCase download is disabled (ENABLE_GENCAGE_DOWNLOAD is FALSE)")
endif()

# Add project subdirectories
add_subdirectory(include)
add_subdirectory(examples)
add_subdirectory(src)
