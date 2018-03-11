#include "TrackHitSorter.h"

#include "LArUtil/Geometry.h"
#include "LArUtil/LArProperties.h"

#include "Geo2D/Core/Geo2D.h"
#include "Geo2D/Core/LineSegment.h"

namespace thsort {

  void TrackHitSorter::buildSortedHitList( const larlite::track& track, const std::vector<larlite::hit>& hit_v, const float max_radius, std::vector<int>& hitmask_v ) {

    // geo utility
    const larutil::Geometry* geo = larutil::Geometry::GetME();
    const float driftv = larutil::LArProperties::GetME()->DriftVelocity();
    const float cm_per_tick = driftv*0.5;
    
    // convert track into line segments
    std::vector< geo2d::LineSegment<float> > seg_v[3]; // segment per plane
    std::vector< float > segdist_v[3]; // distance to the segment

    int numpts = track.NumberTrajectoryPoints();
    int ipt = 0;
    float dist_s = 0;
    while ( ipt+1<numpts ) {
      
      const TVector3& here = track.LocationAtPoint( ipt );
      const TVector3& next = track.LocationAtPoint( ipt+1 );

      //float tick1 = here.X()/cm_per_tick + 3200;
      //float tick2 = next.X()/cm_per_tick + 3200;

      std::cout << "[ipt" << ipt << "] ";
      for (int p=0; p<3; p++) {
	float wire1 = geo->NearestWire( here, p );
	float wire2 = geo->NearestWire( next, p );
	geo2d::LineSegment<float> ls( wire1*0.3, here.X(), wire2*0.3, next.X() ); // tick and wire in cm units. we want the distance to be in cm.
	std::cout << "p" << p << "=[ (" << wire1*0.3 << "," << here.X() << ") (" << wire2*0.3 << "," << next.X() << ") ] ";
	if ( geo2d::length2(ls)>0 ) {
	  seg_v[p].emplace_back( std::move(ls) );
	  segdist_v[p].push_back( dist_s );
	  dist_s += sqrt(geo2d::length2(ls));
	}
      }
      std::cout << std::endl;
      ipt++;
    }
    
    // convert hit positions to same coordinate system
    std::vector< geo2d::Vector<float> > hitco_v[3];
    std::vector< const larlite::hit* > hitp_v[3];
    std::vector< int > hitidx_v[3];
    int ihit=-1;
    for ( auto const& hit : hit_v ) {
      ihit++;
      if ( hitmask_v[ihit]==0 ) continue; // masked
      
      float peakx = (2400+hit.PeakTime() - 3200)*cm_per_tick; // x position
      int plane = (int)hit.WireID().Plane;
      int wire  = (int)hit.WireID().Wire;
      
      //std::cout << "[hit " << ihit << "] x=" << peakx << " p=" << plane << " w=" << wire*0.3 << std::endl;

      geo2d::Vector<float> pt( wire*0.3, peakx );
      hitco_v[plane].emplace_back( std::move(pt) );
      hitp_v[plane].push_back( &hit );
      hitidx_v[plane].push_back( ihit );
    }

    std::cout << "converted hits: p0=" << hitco_v[0].size() << " p1=" << hitco_v[1].size() << " p2=" << hitco_v[2].size() << std::endl;

    // ok, now we associate
    for (int p=0; p<3; p++) {
      // dumb N^2 loop...
      int ihit=-1;
      for ( auto& pt : hitco_v[p] ) {
	ihit++;

	int iseg=-1;
	for ( auto& seg : seg_v[p] ) {
	  iseg++;
	  
	  // check if part of seg
	  geo2d::Vector<float> ab = seg.pt2 - seg.pt1;
	  float s = ab.ddot( pt-seg.pt1 )/geo2d::length2(ab);
	  if ( s<0 || s > 1.0 || std::isnan(s) )
	    continue;

	  
	  geo2d::Vector<float> pt1; // point on line
	  pt1 = seg.pt1 + s*ab;
	  float r = geo2d::dist(pt,pt1);

	  // std::cout << "(" << ihit << "," << iseg << "): "
	  // 	    << " hit[" << pt.x << "," << pt.y << "] "
	  // 	    << " seg[(" << seg.pt1.x << "," << seg.pt1.y << ")->(" << seg.pt2.x << "," << seg.pt2.y << ")] "
	  // 	    << " s= "<< s << " r=" << r << " len(ab)=" << geo2d::length2(ab) << std::endl;
	  
	  
	  if ( r > max_radius ) {
	    continue;
	  }
	  HitOrder ho( hitp_v[p].at(ihit), s+segdist_v[p][iseg], r );
	  ordered[p].emplace_back( std::move(ho) );
	  hitmask_v[ hitidx_v[p].at(ihit) ] = 0; // mask out
	  break;
	}//end of loop over track segment in plane
      }//end of loop over hits in plane
    }//end of loop over plane hits

    // sort hits
    for (int p=0; p<3; p++) {
      std::sort( ordered[p].begin(), ordered[p].end() ); 
    }
  }
  
  void TrackHitSorter::dump() const {
    std::cout << "=========================================================" << std::endl;
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    const larutil::Geometry* geo = larutil::Geometry::GetME();
    const float driftv = larutil::LArProperties::GetME()->DriftVelocity();
    const float cm_per_tick = driftv*0.5;
    
    for (int p=0; p<3; p++) {
      std::cout << "-------------------------------------------" << std::endl;
      std::cout << "Hits on Plane " << p << std::endl;
      for (auto const& ho : ordered[p] ) {
	const larlite::hit* phit = ho.phit;
	std::cout << "  (" << ho.s << "," << ho.r << ") x=" << (2400+phit->PeakTime() - 3200)*cm_per_tick << " w=" << phit->WireID().Wire << std::endl;
      }
    }
    std::cout << "=========================================================" << std::endl;    
  }
  
}
