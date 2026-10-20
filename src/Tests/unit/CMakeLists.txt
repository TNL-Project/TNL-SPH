# Unit tier: C++ tests live under cpp/ (todo 2); the Python unit collector is
# registered with CTest below (todo 30).
add_subdirectory( cpp )

# Python unit collector runner. `--here` discovers and executes every unit
# check group located in this directory, so the test tracks the tier without
# further CMake maintenance. It is invoked as a module (`-m unit.runner`):
# runner.py uses package-relative imports (`from . import collector`), so
# executing it as a plain script cannot work; the interpreter path also keeps
# this independent of shebangs and executable bits. PYTHONPATH makes the
# `unit` package importable from any working directory. Python3 is located at
# the parent src/Tests level; without an interpreter the test is simply not
# registered (see the STATUS message there).
if( Python3_Interpreter_FOUND )
   add_test( NAME unit-tier
             COMMAND ${Python3_EXECUTABLE} -m unit.runner --here )
   set_tests_properties( unit-tier PROPERTIES
      LABELS "unit"
      ENVIRONMENT "PYTHONPATH=${CMAKE_SOURCE_DIR}/src/Tests" )
endif()
