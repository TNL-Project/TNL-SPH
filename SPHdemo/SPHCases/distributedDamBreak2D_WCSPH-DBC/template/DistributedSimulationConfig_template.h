#pragma once

namespace TNL {
namespace ParticleSystem {
namespace Config {

class SubdomainConfig_G1
{
   public:
   static constexpr int particleIdxStart = placeholderParticleIdxStartG1;
   static constexpr int particleIdxRealStart = placeholderParticleIdxRealStartG1;

   static constexpr int particleIdxEnd = placeholderParticleIdxEndG1;
   static constexpr int particleIdxRealEnd = placeholderParticleIdxRealEndG1;

   static constexpr int gridIdxOverlapStar = placeholderGridIdxOverlapStartG1;
   static constexpr int gridIdxStart = placeholderGridIdxStartG1;

   static constexpr int gridIdxOverlapEnd = placeholderGridIdxOverlapEndG1;
   static constexpr int gridIdxEnd = placeholderGridIdxEndG1;
};

class SubdomainConfig_G2
{
   public:
   static constexpr int particleIdxStart = placeholderParticleIdxStartG2;
   static constexpr int particleIdxRealStart = placeholderParticleIdxRealStartG2;

   static constexpr int particleIdxEnd = placeholderParticleIdxEndG2;
   static constexpr int particleIdxRealEnd = placeholderParticleIdxRealEndG2;

   static constexpr int gridIdxOverlapStar = placeholderGridIdxOverlapStartG2;
   static constexpr int gridIdxStart = placeholderGridIdxStartG2;

   static constexpr int gridIdxOverlapEnd = placeholderGridIdxOverlapEndG2;
   static constexpr int gridIdxEnd = placeholderGridIdxEndG2;
};


} //namespace ParticleSystemConfig
} //namespace ParticleSystem
} //namespace TNL

