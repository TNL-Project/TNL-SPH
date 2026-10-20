# Unified test suite root. The regression tier is pure Python and needs no CMake wiring.

# The unit tier's Python runner is registered with CTest under unit/. Tests are
# optional: without a Python interpreter they are simply not registered - the
# configure must never fail on a missing interpreter, only report it.
find_package( Python3 COMPONENTS Interpreter )
if( NOT Python3_Interpreter_FOUND )
   message( STATUS "Python3 interpreter not found - the 'unit-tier' CTest test will not be registered." )
endif()

add_subdirectory( unit )
add_subdirectory( infrastructure )
