namespace TNL {
namespace ParticleSystem {

template< typename ParticlesConfig, typename Reader >
class ReadParticles
{
public:

   ReadParticles( const std::string& inputFileName )
   : reader( inputFileName )
   {
      reader.detectParticleSystem();
   }

   template< typename PointArray >
   void readParticles( PointArray& particles )
   {
      using ParticleSystemToReadData = typename ParticleSystem::Particles< ParticlesConfig, Devices::Host >;
      ParticleSystemToReadData particlesToRead( ParticlesConfig::numberOfParticles, ParticlesConfig::searchRadius );
      reader.template loadParticle< ParticleSystemToReadData >( particlesToRead );

      PointArray pointsLoaded( ParticlesConfig::numberOfParticles );
      pointsLoaded = particlesToRead.getPoints();
      particles = pointsLoaded;
   }

   template< typename Array, typename Type >
   void readParticleVariable( Array& array, const std::string& name )
   {
      Array arrayLoaded( array.getSize() );
      arrayLoaded = std::get< std::vector< Type > >( reader.readPointData( name ) );
      //it would be nice to have type from Array::ValueType, but I need float for vector anyway.
      array = arrayLoaded;
   }

   template< typename Array, typename Type >
   void readParticleVariable2D( Array& array, const std::string& name )
   {
      Array arrayLoaded( array.getSize() );
      std::vector< Type > temporary = std::get< std::vector< Type > >( reader.readPointData( name ) );
      std::cout << " Loaded " << name << " size: " << temporary.size() << std::endl;


      using HostArray = typename Array::
         template Self< std::remove_const_t< typename Array::ValueType >, Devices::Host, typename Array::IndexType >;

      HostArray hostArray( array.getSize() );

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

};

} //ParticleSystem
} //TNL

