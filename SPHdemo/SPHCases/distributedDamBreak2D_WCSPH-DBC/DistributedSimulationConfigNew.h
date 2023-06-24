#pragma once

namespace TNL {
namespace ParticleSystem {
namespace Config {

class SubdomainConfig_G1
{
   public:
   int particleIdxStart = 0;
   int particleIdxRealStart = 0;

   int particleIdxEnd = 22049;
   int particleIdxRealEnd = 22499;

   int gridIdxOverlapStar = 0;
   int gridIdxStart = 0;

   int gridIdxOverlapEnd = 54;
   int gridIdxEnd = 53;
};

class SubdomainConfig_G2
{
   public:
   int particleIdxStart = 0;
   int particleIdxRealStart = 0;

   int particleIdxEnd = 22049;
   int particleIdxRealEnd = 22499;

   int gridIdxOverlapStar = 53;
   int gridIdxStart = 54;

   int gridIdxOverlapEnd = 290;
   int gridIdxEnd = 290;
};


} //namespace ParticleSystemConfig
} //namespace ParticleSystem
} //namespace TNL

