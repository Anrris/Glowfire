//
// Created by Tai, Yuan-yen on 11/17/17.
// All rights reserved.
//

#ifndef CLASSIFIER_BASE_H
#define CLASSIFIER_BASE_H

#include <cmath>                               /* STL Library*/
#include <string>                               /* STL Library*/
#include <sstream>                              /* STL Library*/
#include <iostream>                             /* STL Library*/
#include <vector>                               /* STL Library*/
#include <list>                                 /* STL Library*/
#include <unordered_set>                        /* STL Library*/
#include <memory>                               /* STL Library*/
#include <limits>                               /* STL Library*/
#include <exception>                            /* STL Library*/
#include <functional>                           /* STL Library*/

#include <boost/geometry.hpp>                   /* Boost Library: For RTree*/
#include <boost/geometry/geometries/point.hpp>  /* Boost Library */
#include <boost/geometry/geometries/box.hpp>    /* Boost Library */
#include <boost/geometry/index/rtree.hpp>       /* Boost Library */

#include <Eigen/Dense>                          /* Eigen Library: For matrix operation*/

#include "glassfire_util.h"
#include "glassfire_cluster_model.h"
#include "glassfire_copier.h"

namespace glassfire{


namespace bg  = boost::geometry;
namespace bgm = boost::geometry::model;
namespace bgi = boost::geometry::index;

/* Brief:
 * Feature              : The basic feature type.
 *
 * RtreePoint           : The element for the Rtree of data points.
 * RtreeBox             : The Rtree search box of a certain region.
 * RtreeValue           : A pair of value for the Rtree and related index.
 * Rtree                : The Rtree.
 *
 * static_for           : Looping through the template static variable.
 *
 * RtreeFeature         : The parent of RtreeFeature is RtreePoint.
 * RtreeFeature_s       : A container of RtreeFeatures
 *
 * Centroid             : The centroid of a cluster.
 * CentroidList         : A collection of the centroid.
 * CentroidListIterator : The iterator of the CentroidList.
 * CentroidRtree        : The Rtree of a Centroid.
 * */


template<typename AxisType, size_t Dimension, typename FeatureInfo>
class GlassfireType
{
public:
    typedef std::vector<AxisType>                Feature;
    typedef std::vector<Feature>                 Feature_s;

    typedef bgm::point<AxisType, Dimension, bg::cs::cartesian>  RtreePoint;
    typedef bgm::box<RtreePoint>                                RtreeBox;
    typedef std::pair<RtreePoint, FeatureInfo>                  RtreeValue;
    typedef bgi::rtree<RtreeValue, bgi::linear<16>>             Rtree;

    typedef Eigen::Matrix<AxisType, Dimension, Dimension>       Matrix;
    typedef Eigen::Matrix<AxisType, 1, Dimension>               RowVector;
    typedef Eigen::Matrix<AxisType, Dimension, 1>               ColVector;

    //===================================
    //--- Implementation of RtreeFeature
    //-----------------------------------
    class RtreeFeature: public RtreePoint
    {
        //static_for  m_static_for;
        point_copier<AxisType, Dimension>    copier;
        Feature     mFeature;
    public:
        RtreeFeature(){}
        RtreeFeature(const Feature & feature): mFeature(feature) {
            setFeature(feature);
        }
        void setFeature(const Feature & feature){
            mFeature = feature;
            copier.copy_to(mFeature, *this);
        }
        auto getFeature() -> const Feature & {
            return mFeature;
        }
        auto at(size_t index) -> const AxisType & {
            return mFeature[index];
        }
        auto createBox(AxisType centroid_distance) -> RtreeBox {
            RtreePoint pt1, pt2;
            copier.create_box_corner(centroid_distance/2, mFeature, pt1, pt2);
            return RtreeBox(pt1, pt2);
        }
        auto getPoint() -> const RtreePoint & { return *this; }
    };

    typedef std::vector<RtreeFeature>    RtreeFeature_s;
    //-----------------------------------
    //--- Implementation of RtreeFeature
    //===================================

    typedef ClusterModel<AxisType, FeatureInfo>    ClusterModel;

    class Centroid;

    typedef typename Centroid::CentroidList         CentroidList;
    typedef typename Centroid::CentroidListPtr      CentroidListPtr;
    typedef typename Centroid::CentroidListIterator CentroidListIterator;

    typedef typename Centroid::CentroidRtreeValue   CentroidRtreeValue;
    typedef typename Centroid::CentroidRtree        CentroidRtree;
    typedef std::shared_ptr<CentroidRtree>               CentroidRtreePtr;

    class ScorerSet; 
};

}

#include "glassfire_centroid.h"
#include "glassfire_scorer.h"

#endif //CLASSIFIER_BASE_H
