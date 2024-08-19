namespace TNL {
namespace SPH {
namespace features {

template< typename ObjectPointer, typename ModelParams, typename RealType >
void checkWeaklyCompressibleAssumption( ObjectPointer& object,
                                        ModelParams& modelParams,
                                        const RealType& dt,
                                        TNL::Logger & logger )
{
   using DeviceType = typename ModelParams::DeviceType;

   const auto rhoView = object->getVariables()->rho.getConstView();
   const auto drhoView = object->getVariables()->drho.getConstView();

   auto fetch = [=] __cuda_callable__ ( int i )
   {
      const RealType drhoRhoFrac_i = TNL::fabs( drho[ i ] * dt / rho[ i ] );
      return ( drhoRhoFrac_i < 0.01f ) ? ( 1 ) : ( 0 );
   };

   const RealType violatedCount = Algorithms::reduce< DeviceType >( 0,
                                                                    object->getNumberOfParticles(),
                                                                    fetch,
                                                                    TNL::Sum() );
   if( violatedCount > 0 ){
      logger.writeParameter( "WCSPH error: ", "WC assumption violated" );
      logger.writeParameter( "Number of particles violating the assumption:", violatedCount );
   }
}

} // features
} // SPH
} // TNL

