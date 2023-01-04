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

protected:

   Reader reader;

};

} //ParticleSystem
} //TNL

