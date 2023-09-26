#include <cfloat> //FLT_MAX

namespace TNL {
namespace ParticleSystem {

template< typename ParticlesConfig, typename Reader >
class ReadParticles
{
public:
   using GlobalIndexType = typename ParticlesConfig::GlobalIndexType;

   ReadParticles( const std::string& inputFileName, const GlobalIndexType& numberOfParticles, const GlobalIndexType numberOfAllocatedParticles )
   : reader( inputFileName ), numberOfParticles( numberOfParticles ), numberOfAllocatedParticles( numberOfAllocatedParticles )
   {
      reader.detectParticleSystem();
   }

   template< typename PointArray >
   void readParticles( PointArray& particles )
   {
      using ParticleSystemToReadData = typename ParticleSystem::Particles< ParticlesConfig, Devices::Host >;
      ParticleSystemToReadData particlesToRead( numberOfParticles, numberOfParticles, 0. );
      reader.template loadParticle< ParticleSystemToReadData >( particlesToRead );

      PointArray pointsLoaded( numberOfParticles );
      //std::cout << "ParticlesLoaded size: " << particlesToRead.getPoints() << std::endl;
      pointsLoaded = particlesToRead.getPoints();
      pointsLoaded.resize( numberOfAllocatedParticles, FLT_MAX );
      particles = pointsLoaded;
   }

   template< typename Array, typename Type >
   void readParticleVariable( Array& array, const std::string& name )
   {
      //Array arrayLoaded( array.getSize() );
      Array arrayLoaded(  reader.getNumberOfPoints()  );
      arrayLoaded = std::get< std::vector< Type > >( reader.readPointData( name ) );
      //it would be nice to have type from Array::ValueType, but I need float for vector anyway.
      arrayLoaded.resize( numberOfAllocatedParticles, FLT_MAX );
      array = arrayLoaded;
   }

   template< typename Array, typename Type >
   void readParticleVariable2D( Array& array, const std::string& name )
   {
      //Array arrayLoaded( array.getSize() );
      Array arrayLoaded(  reader.getNumberOfPoints()  );
      std::vector< Type > temporary = std::get< std::vector< Type > >( reader.readPointData( name ) );
      std::cout << " Loaded " << name << " size: " << temporary.size() << std::endl;


      using HostArray = typename Array::
         template Self< std::remove_const_t< typename Array::ValueType >, Devices::Host, typename Array::IndexType >;

      //HostArray hostArray( array.getSize() );
      HostArray hostArray( reader.getNumberOfPoints() );

      std::vector< Type > hostBuffer;
      for( int i = 1; i < temporary.size() + 1; i ++ )
      {
         if( ( i % 3 == 0 ) )
          continue;

         hostBuffer.push_back( temporary[ i - 1 ] );
      }

      int counter = 0;
      for( int i = 0; i < hostBuffer.size(); i++ )
      {
         //std::cout << hostBuffer[ i ] << " ";
         if( i % 2 == 0 )
            continue;

         typename HostArray::ValueType aux = { hostBuffer[ i - 1 ], hostBuffer[ i ] };
         hostArray[ counter ] = aux;
         counter++;
         //std::cout << aux << std::endl;
      }

      array = hostArray;
   }

protected:

   Reader reader;
   GlobalIndexType numberOfParticles;
   GlobalIndexType numberOfAllocatedParticles;

};

} //ParticleSystem
} //TNL

