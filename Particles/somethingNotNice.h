//: //GlobalIndexType localBeginX = 0;
//: //GlobalIndexType localBeginY = 0;
//: //GlobalIndexType interiorEndX = 7;
//: //GlobalIndexType interiorEndY = 7;
//:
//: GlobalIndexType localBeginX = 0;
//: GlobalIndexType localBeginY = 0;
//: GlobalIndexType interiorEndX = particles.grid->getInteriorEnd()[0];
//: GlobalIndexType interiorEndY = particles.grid->getInteriorEnd()[1];
//:
//: if( centralCell.isBoundary() == false )
//: {
//:   const LocalIndexType mp = neighborEntities.template getEntityIndex< -1, 1 >();
//:   NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mp);
//:   const LocalIndexType zp = neighborEntities.template getEntityIndex< 0, 1 >();
//:   NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zp);
//:   const LocalIndexType pp = neighborEntities.template getEntityIndex< 1, 1 >();
//:   NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, pp);
//:   const LocalIndexType mz = neighborEntities.template getEntityIndex< -1, 0 >();
//:   NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mz);
//:   const LocalIndexType zz = neighborEntities.template getEntityIndex< 0, 0 >();
//:   NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zz);
//:   const LocalIndexType pz = neighborEntities.template getEntityIndex< 1, 0 >();
//:   NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, pz);
//:   const LocalIndexType mm = neighborEntities.template getEntityIndex< -1, -1 >();
//:   NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mm);
//:   const LocalIndexType zm = neighborEntities.template getEntityIndex< 0, -1 >();
//:   NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zm);
//:   const LocalIndexType pm = neighborEntities.template getEntityIndex< 1, -1 >();
//:   NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, pm);
//:   //printf("[%d, %d, %d, %d, %d, %d, %d, %d, %d]\n", mp, zp, pp, mz, zz, pz, mm, zm, pm);
//: }
//: else
//: {
//:   if(centralCell.getCoordinates()[0] == localBeginX)
//:   {
//:     if(centralCell.getCoordinates()[1] == localBeginY)
//:     {
//:       //printf( " /0,0/ " );
//:       const LocalIndexType zp = neighborEntities.template getEntityIndex< 0, 1 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zp);
//:       const LocalIndexType pp = neighborEntities.template getEntityIndex< 1, 1 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, pp);
//:       const LocalIndexType zz = neighborEntities.template getEntityIndex< 0, 0 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zz);
//:       const LocalIndexType pz = neighborEntities.template getEntityIndex< 1, 0 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, pz);
//:       //printf("[%d, %d, %d, %d]\n", zp, pp, zz, pz);
//:     }
//:     // else if(centralCell.getCoordinates()[1] == 7)
//:     else if(centralCell.getCoordinates()[1] == interiorEndY)
//:     {
//:       //printf(" /0,9/ ");
//:       const LocalIndexType zz = neighborEntities.template getEntityIndex< 0, 0 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zz);
//:       const LocalIndexType zp = neighborEntities.template getEntityIndex< 1, 0 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zp);
//:       const LocalIndexType zm = neighborEntities.template getEntityIndex< 0, -1 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zm);
//:       const LocalIndexType mm = neighborEntities.template getEntityIndex< 1, -1 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mm);
//:       //printf("[%d, %d, %d, %d]\n", zz, zp, zm, mm);
//:     }
//:     else
//:     {
//:       const LocalIndexType zp = neighborEntities.template getEntityIndex< 0, 1 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zp);
//:       const LocalIndexType pp = neighborEntities.template getEntityIndex< 1, 1 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, pp);
//:       const LocalIndexType zz = neighborEntities.template getEntityIndex< 0, 0 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zz);
//:       const LocalIndexType pz = neighborEntities.template getEntityIndex< 1, 0 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, pz);
//:       const LocalIndexType zm = neighborEntities.template getEntityIndex< 0, -1 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zm);
//:       const LocalIndexType pm = neighborEntities.template getEntityIndex< 1, -1 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, pm);
//:       //printf("[%d, %d, %d, %d, %d, %d]\n", zp, pp, zz, pz, zm, pm);
//:     }
//:   }
//:   else if(centralCell.getCoordinates()[0] == interiorEndX)
//:   {
//:     if(centralCell.getCoordinates()[1] == localBeginY)
//:     {
//:       //printf(" /9,0/ ");
//:       const LocalIndexType mp = neighborEntities.template getEntityIndex< -1, 1 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mp);
//:       const LocalIndexType zp = neighborEntities.template getEntityIndex< 0, 1 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zp);
//:       const LocalIndexType mz = neighborEntities.template getEntityIndex< -1, 0 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mz);
//:       const LocalIndexType zz = neighborEntities.template getEntityIndex< 0, 0 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zz);
//:       //printf("[%d, %d, %d, %d]\n", mp, zp, mz, zz);
//:     }
//:     else if(centralCell.getCoordinates()[1] == interiorEndY)
//:     {
//:       //printf(" /9,9/ ");
//:       const LocalIndexType mz = neighborEntities.template getEntityIndex< -1, 0 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mz);
//:       const LocalIndexType zz = neighborEntities.template getEntityIndex< 0, 0 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zz);
//:       const LocalIndexType mm = neighborEntities.template getEntityIndex< -1, -1 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mm);
//:       const LocalIndexType zm = neighborEntities.template getEntityIndex< 0, -1 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zm);
//:       //printf("[%d, %d, %d, %d]\n", mz, zz, mm, zm);
//:     }
//:     else
//:     {
//:       const LocalIndexType mp = neighborEntities.template getEntityIndex< -1, 1 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mp);
//:       const LocalIndexType zp = neighborEntities.template getEntityIndex< 0, 1 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zp);
//:       const LocalIndexType mz = neighborEntities.template getEntityIndex< -1, 0 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mz);
//:       const LocalIndexType zz = neighborEntities.template getEntityIndex< 0, 0 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zz);
//:       const LocalIndexType mm = neighborEntities.template getEntityIndex< -1, -1 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mm);
//:       const LocalIndexType zm = neighborEntities.template getEntityIndex< 0, -1 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zm);
//:       //printf("[%d, %d, %d, %d, %d, %d]\n", mp, zp, mz, zz, mm, zm);
//:     }
//:   }
//:   else
//:   {
//:     if(centralCell.getCoordinates()[1] == localBeginY)
//:     {
//:       const LocalIndexType mp = neighborEntities.template getEntityIndex< -1, 1 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mp);
//:       const LocalIndexType zp = neighborEntities.template getEntityIndex< 0, 1 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zp);
//:       const LocalIndexType pp = neighborEntities.template getEntityIndex< 1, 1 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, pp);
//:       const LocalIndexType mz = neighborEntities.template getEntityIndex< -1, 0 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mz);
//:       const LocalIndexType zz = neighborEntities.template getEntityIndex< 0, 0 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zz);
//:       const LocalIndexType pz = neighborEntities.template getEntityIndex< 1, 0 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, pz);
//:       //printf("[%d, %d, %d, %d, %d, %d]\n", mp, zp, pp, mz, zz, pz);
//:     }
//:     else
//:     {
//:       const LocalIndexType mz = neighborEntities.template getEntityIndex< -1, 0 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mz);
//:       const LocalIndexType zz = neighborEntities.template getEntityIndex< 0, 0 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zz);
//:       const LocalIndexType pz = neighborEntities.template getEntityIndex< 1, 0 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, pz);
//:       const LocalIndexType mm = neighborEntities.template getEntityIndex< -1, -1 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mm);
//:       const LocalIndexType zm = neighborEntities.template getEntityIndex< 0, -1 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zm);
//:       const LocalIndexType pm = neighborEntities.template getEntityIndex< 1, -1 >();
//:       NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, pm);
//:       //printf("[%d, %d, %d, %d, %d, %d]\n", mz, zz, pz, mm, zm, pm);
//:     }
//:   }
//: }

//GlobalIndexType localBeginX = 0;
//GlobalIndexType localBeginY = 0;
//GlobalIndexType interiorEndX = 7;
//GlobalIndexType interiorEndY = 7;

GlobalIndexType localBeginX = 0;
GlobalIndexType localBeginY = 0;
GlobalIndexType interiorEndX = particles.grid->getInteriorEnd()[0];
GlobalIndexType interiorEndY = particles.grid->getInteriorEnd()[1];

if( centralCell.isBoundary() == false )
{
	const LocalIndexType mp = centralCell.template getNeighbourEntity< 2, -1, 1 >().getIndex();
  NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mp);

	const LocalIndexType zp = centralCell.template getNeighbourEntity< 2, 0, 1 >().getIndex();
  NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zp);

	const LocalIndexType pp = centralCell.template getNeighbourEntity< 2, 1, 1 >().getIndex();
  NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, pp);

	const LocalIndexType mz = centralCell.template getNeighbourEntity< 2, -1, 0 >().getIndex();
  NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mz);

  const LocalIndexType zz = centralCell.template getNeighbourEntity< 2, 0, 0 >().getIndex();
  NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zz);

  const LocalIndexType pz = centralCell.template getNeighbourEntity< 2, 1, 0 >().getIndex();
  NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, pz);

  const LocalIndexType mm = centralCell.template getNeighbourEntity< 2, -1, -1 >().getIndex();
  NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mm);

  const LocalIndexType zm = centralCell.template getNeighbourEntity< 2, 0, -1 >().getIndex();
  NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zm);

  const LocalIndexType pm = centralCell.template getNeighbourEntity< 2, 1, -1 >().getIndex();
  NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, pm);

  //printf("[%d, %d, %d, %d, %d, %d, %d, %d, %d]\n", mp, zp, pp, mz, zz, pz, mm, zm, pm);
}
else
{
  if(centralCell.getCoordinates()[0] == localBeginX)
  {
    if(centralCell.getCoordinates()[1] == localBeginY)
    {
      //printf( " /0,0/ " );
      const LocalIndexType zp = centralCell.template getNeighbourEntity< 2, 0, 1 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zp);

      const LocalIndexType pp = centralCell.template getNeighbourEntity< 2, 1, 1 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, pp);

      const LocalIndexType zz = centralCell.template getNeighbourEntity< 2, 0, 0 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zz);

      const LocalIndexType pz = centralCell.template getNeighbourEntity< 2, 1, 0 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, pz);

      //printf("[%d, %d, %d, %d]\n", zp, pp, zz, pz);
    }
    // else if(centralCell.getCoordinates()[1] == 7)
    else if(centralCell.getCoordinates()[1] == interiorEndY)
    {
      //printf(" /0,9/ ");
      const LocalIndexType zz = centralCell.template getNeighbourEntity< 2, 0, 0 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zz);

      const LocalIndexType zp = centralCell.template getNeighbourEntity< 2, 1, 0 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zp);

      const LocalIndexType zm = centralCell.template getNeighbourEntity< 2, 0, -1 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zm);

      const LocalIndexType mm = centralCell.template getNeighbourEntity< 2, 1, -1 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mm);

      //printf("[%d, %d, %d, %d]\n", zz, zp, zm, mm);
    }
    else
    {
      const LocalIndexType zp = centralCell.template getNeighbourEntity< 2, 0, 1 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zp);

      const LocalIndexType pp = centralCell.template getNeighbourEntity< 2, 1, 1 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, pp);

      const LocalIndexType zz = centralCell.template getNeighbourEntity< 2, 0, 0 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zz);

      const LocalIndexType pz = centralCell.template getNeighbourEntity< 2, 1, 0 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, pz);

      const LocalIndexType zm = centralCell.template getNeighbourEntity< 2, 0, -1 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zm);

      const LocalIndexType pm = centralCell.template getNeighbourEntity< 2, 1, -1 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, pm);

      //printf("[%d, %d, %d, %d, %d, %d]\n", zp, pp, zz, pz, zm, pm);
    }
  }
  else if(centralCell.getCoordinates()[0] == interiorEndX)
  {
    if(centralCell.getCoordinates()[1] == localBeginY)
    {
      //printf(" /9,0/ ");
      const LocalIndexType mp = centralCell.template getNeighbourEntity< 2, -1, 1 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mp);

      const LocalIndexType zp = centralCell.template getNeighbourEntity< 2, 0, 1 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zp);

      const LocalIndexType mz = centralCell.template getNeighbourEntity< 2, -1, 0 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mz);

      const LocalIndexType zz = centralCell.template getNeighbourEntity< 2, 0, 0 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zz);

      //printf("[%d, %d, %d, %d]\n", mp, zp, mz, zz);
    }
    else if(centralCell.getCoordinates()[1] == interiorEndY)
    {
      //printf(" /9,9/ ");
      const LocalIndexType mz = centralCell.template getNeighbourEntity< 2, -1, 0 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mz);

      const LocalIndexType zz = centralCell.template getNeighbourEntity< 2, 0, 0 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zz);

      const LocalIndexType mm = centralCell.template getNeighbourEntity< 2, -1, -1 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mm);

      const LocalIndexType zm = centralCell.template getNeighbourEntity< 2, 0, -1 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zm);

      //printf("[%d, %d, %d, %d]\n", mz, zz, mm, zm);
    }
    else
    {
      const LocalIndexType mp = centralCell.template getNeighbourEntity< 2, -1, 1 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mp);

      const LocalIndexType zp = centralCell.template getNeighbourEntity< 2, 0, 1 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zp);

      const LocalIndexType mz = centralCell.template getNeighbourEntity< 2, -1, 0 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mz);

      const LocalIndexType zz = centralCell.template getNeighbourEntity< 2, 0, 0 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zz);

      const LocalIndexType mm = centralCell.template getNeighbourEntity< 2, -1, -1 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mm);

      const LocalIndexType zm = centralCell.template getNeighbourEntity< 2, 0, -1 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zm);

      //printf("[%d, %d, %d, %d, %d, %d]\n", mp, zp, mz, zz, mm, zm);
    }
  }
  else
  {
    if(centralCell.getCoordinates()[1] == localBeginY)
    {
      const LocalIndexType mp = centralCell.template getNeighbourEntity< 2, -1, 1 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mp);

      const LocalIndexType zp = centralCell.template getNeighbourEntity< 2, 0, 1 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zp);

      const LocalIndexType pp = centralCell.template getNeighbourEntity< 2, 1, 1 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, pp);

      const LocalIndexType mz = centralCell.template getNeighbourEntity< 2, -1, 0 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mz);

      const LocalIndexType zz = centralCell.template getNeighbourEntity< 2, 0, 0 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zz);

      const LocalIndexType pz = centralCell.template getNeighbourEntity< 2, 1, 0 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, pz);

      //printf("[%d, %d, %d, %d, %d, %d]\n", mp, zp, pp, mz, zz, pz);
    }
    else
    {
      const LocalIndexType mz = centralCell.template getNeighbourEntity< 2, -1, 0 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mz);

      const LocalIndexType zz = centralCell.template getNeighbourEntity< 2, 0, 0 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zz);

      const LocalIndexType pz = centralCell.template getNeighbourEntity< 2, 1, 0 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, pz);

      const LocalIndexType mm = centralCell.template getNeighbourEntity< 2, -1, -1 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mm);

      const LocalIndexType zm = centralCell.template getNeighbourEntity< 2, 0, -1 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zm);

      const LocalIndexType pm = centralCell.template getNeighbourEntity< 2, 1, -1 >().getIndex();
      NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, pm);

      //printf("[%d, %d, %d, %d, %d, %d]\n", mz, zz, pz, mm, zm, pm);
    }
  }
}

