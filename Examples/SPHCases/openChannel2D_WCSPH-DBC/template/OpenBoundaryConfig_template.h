namespace TNL {
namespace ParticleSystem {
namespace SPH {

/**
 * PARAMETERS OF OPEN BOUNDARY PATCH
 *
 * This class is used to store core parameters for inlet boundary patch
 * i.e. inlet or outlet. The values are used only for initialization.
 *
 * It is necessary to enter:
 *
 * - orientation_x - x component of normal buffer vector
 * - orientation_y - y component of normal buffer vector
 * - velocity_x - initial x component of open boundary patch velocity
 * - velocity_y - initial x component of open boundary patch velocity
 * - position_x - referential position of open boundary buffer TODO: Move to centre.
 * - inlet_density - referential position of open boundary buffer TODO: Move to centre.
 * - bufferWidth_x - width of buffer - dependent on number of layers
 * - bufferWidth_y - width of buffer - dependent on number of layers
 * - bufferEdge - DEPRECATED
 *
 */
struct InletBuffer
{
   static constexpr float orientation_x = placeholderOBP1Orientation_xf;
   static constexpr float orientation_y = placeholderOBP1Orientation_yf;
   static constexpr float velocity_x = placeholderOBP1Velocity_xf;
   static constexpr float velocity_y = placeholderOBP1Velocity_yf;
   static constexpr float position_x = placeholderOBP1Position_xf;
   static constexpr float position_y = placeholderOBP1Position_yf;
   static constexpr float inlet_density = placeholderOBP1Densityf;
   static constexpr float bufferWidth_x = placeholderOBP1Width_xf; //ie 4 layers
   static constexpr float bufferWidth_y = placeholderOBP1Width_yf; //ie 4 layers
   static constexpr float bufferEdge = placeholderOBP1BufferEdgef; //ie 4 layers
};

/**
 * PARAMETERS OF OPEN BOUNDARY PATCH
 *
 * This class is used to store core parameters for inlet boundary patch
 * i.e. inlet or outlet. The values are used only for initialization.
 *
 * It is necessary to enter:
 *
 * - orientation_x - x component of normal buffer vector
 * - orientation_y - y component of normal buffer vector
 * - velocity_x - initial x component of open boundary patch velocity
 * - velocity_y - initial x component of open boundary patch velocity
 * - position_x - referential position of open boundary buffer TODO: Move to centre.
 * - inlet_density - referential position of open boundary buffer TODO: Move to centre.
 * - bufferWidth_x - width of buffer - dependent on number of layers
 * - bufferWidth_y - width of buffer - dependent on number of layers
 * - bufferEdge - DEPRECATED
 *
 */
struct OutletBuffer
{
   static constexpr float orientation_x = placeholderOBP2Orientation_xf;
   static constexpr float orientation_y = placeholderOBP2Orientation_yf;
   static constexpr float velocity_x = placeholderOBP2Velocity_xf;
   static constexpr float velocity_y = placeholderOBP2Velocity_yf;
   static constexpr float position_x = placeholderOBP2Position_xf;
   static constexpr float position_y = placeholderOBP2Position_yf;
   static constexpr float inlet_density = placeholderOBP2Densityf;
   static constexpr float bufferWidth_x = placeholderOBP2Width_xf; //ie 4 layers
   static constexpr float bufferWidth_y = placeholderOBP2Width_yf; //ie 4 layers
   static constexpr float bufferEdge = placeholderOBP2BufferEdgef; //ie 4 layers
};

} //SPH
} //namespace ParticleSystem
} //namespace TNL

