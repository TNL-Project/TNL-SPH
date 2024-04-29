#include <gtest/gtest.h>
#include "GtestPrintToOverrides.h"

#if defined( HAVE_MPI )
   #include <TNL/MPI/ScopedInitializer.h>
   #include <TNL/MPI/Wrappers.h>

   #include <sstream>

class MinimalistBufferedPrinter : public ::testing::EmptyTestEventListener
{
private:
   std::stringstream sout;

public:
   // Called before a test starts.
   void
   OnTestStart( const ::testing::TestInfo& test_info ) override
   {
      sout << test_info.test_case_name() << "." << test_info.name() << " Start." << std::endl;
   }

   // Called after a failed assertion or a SUCCEED() invocation.
   void
   OnTestPartResult( const ::testing::TestPartResult& test_part_result ) override
   {
      sout << ( test_part_result.failed() ? "====Failure=== " : "===Success=== " ) << test_part_result.file_name() << " "
           << test_part_result.line_number() << std::endl
           << test_part_result.summary() << std::endl;
   }

   // Called after a test ends.
   void
   OnTestEnd( const ::testing::TestInfo& test_info ) override
   {
      const int rank = TNL::MPI::GetRank();
      sout << test_info.test_case_name() << "." << test_info.name() << " End." << std::endl;
      std::cout << rank << ":" << std::endl << sout.str() << std::endl;
      sout.str( std::string() );
      sout.clear();
   }
};
#endif

int
main( int argc, char* argv[] )
{
   ::testing::InitGoogleTest( &argc, argv );

#ifdef HAVE_MPI
   TNL::MPI::ScopedInitializer mpi( argc, argv );

   if( TNL::MPI::GetSize() > 1 ) {
      ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();

      delete listeners.Release( listeners.default_result_printer() );
      listeners.Append( new MinimalistBufferedPrinter );
   }
#endif

   return RUN_ALL_TESTS();
}

