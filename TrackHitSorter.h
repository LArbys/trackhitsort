#ifndef __TRACKHITSORTER_H__
#define __TRACKHITSORTER_H__

#include <vector>

// larlite
#include "DataFormat/vertex.h"
#include "DataFormat/track.h"
#include "DataFormat/hit.h"

namespace thsort {

  class HitOrder {
  public:
    HitOrder ( const larlite::hit* phit_, float s_, float r_, float d_ ) : phit(phit_), s(s_), r(r_), d(d_) {};
    ~HitOrder() {};
    
    const larlite::hit* phit; // pointer to hit
    float s; // distance to start along track
    float r; // distance from track seg to hit
    float d; // distance from point to vertex in 3d	

    bool operator< ( const HitOrder& rh ) const {
      if ( s < rh.s ) return true;
      return false;
    };
    
  };
  
  class TrackHitSorter {

  public:
    TrackHitSorter(){};
    ~TrackHitSorter(){};

    void buildSortedHitList( const larlite::vertex& vtx, const larlite::track& track, const std::vector<larlite::hit>& hit_v,
			     const float max_radius, std::vector<int>& hitmask_v );
    void dump() const;
    
    std::vector<HitOrder> ordered[3]; // per plane

    
  };


}


#endif
