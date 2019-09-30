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

namespace glassfire{


using namespace std;
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
    typedef vector<AxisType>                Feature;
    typedef vector<Feature>                 Feature_s;

    typedef bgm::point<AxisType, Dimension, bg::cs::cartesian>  RtreePoint;
    typedef bgm::box<RtreePoint>                                RtreeBox;
    typedef pair<RtreePoint, FeatureInfo>                       RtreeValue;
    typedef bgi::rtree<RtreeValue, bgi::linear<16>>             Rtree;

    typedef Eigen::Matrix<AxisType, Dimension, Dimension>       Matrix;
    typedef Eigen::Matrix<AxisType, 1, Dimension>               RowVector;
    typedef Eigen::Matrix<AxisType, Dimension, 1>               ColVector;

    //=================================
    //--- Implementation of static_for
    //---------------------------------
    template<size_t x, size_t to>
    struct t_static_for
    {
        AxisType m_val_min = numeric_limits<AxisType>::max();
        AxisType m_val_max = numeric_limits<AxisType>::min();

        void checkDimension(size_t dim){
            if(dim != Dimension){
                throw std::length_error("Input data does not match required dimension.");
            }
        }

        void copy_to(Feature &coord, RtreePoint &pt) {
            checkDimension(coord.size());

            m_val_min = min(m_val_min, coord.at(x));
            m_val_max = max(m_val_max, coord.at(x));

            pt.set<x>(coord.at(x));
            t_static_for<x+1,to>().copy_to(coord, pt);
        }

        void create_box_corner(AxisType dist, Feature &coord, RtreePoint &pt1, RtreePoint &pt2) {
            checkDimension(coord.size());

            pt1.set<x>(coord.at(x)-dist);
            pt2.set<x>(coord.at(x)+dist);
            t_static_for<x+1,to>().create_box_corner(dist, coord, pt1, pt2);
        }
    };

    template<size_t to>
    struct t_static_for<to, to>
    {
        void copy_to(Feature &arr, RtreePoint &pt) {}
        void create_box_corner(AxisType dist, Feature &coord, RtreePoint &pt1, RtreePoint &pt2) {}
    };

    typedef t_static_for<0, Dimension>    static_for;
    //---------------------------------
    //--- Implementation of static_for
    //=================================

    //===================================
    //--- Implementation of RtreeFeature
    //-----------------------------------
    class RtreeFeature: public RtreePoint
    {
        static_for  m_static_for;
        Feature     mFeature;
    public:
        RtreeFeature(const Feature & feature): mFeature(feature) {
            setFeature(feature);
        }
        void setFeature(const Feature & feature){
            mFeature = feature;
            m_static_for.copy_to(mFeature, *this);
        }
        auto getFeature() -> const Feature & {
            return mFeature;
        }
        auto at(size_t index) -> const AxisType & {
            return mFeature[index];
        }
        auto createBox(AxisType centroid_distance) -> RtreeBox {
            RtreePoint pt1, pt2;
            m_static_for.create_box_corner(centroid_distance/2, mFeature, pt1, pt2);
            return RtreeBox(pt1, pt2);
        }
        auto getPoint() -> const RtreePoint & { return *this; }
    };

    typedef vector<RtreeFeature>    RtreeFeature_s;
    //-----------------------------------
    //--- Implementation of RtreeFeature
    //===================================

    class ClusterModel;

    class Centroid;

    typedef typename Centroid::CentroidList         CentroidList;
    typedef typename Centroid::CentroidListPtr      CentroidListPtr;
    typedef typename Centroid::CentroidListIterator CentroidListIterator;

    typedef typename Centroid::CentroidRtreeValue   CentroidRtreeValue;
    typedef typename Centroid::CentroidRtree        CentroidRtree;
    typedef shared_ptr<CentroidRtree>               CentroidRtreePtr;

    class ScorerSet; 
};

}

#include "glassfire_cluster_model.h"
#include "glassfire_centroid.h"
#include "glassfire_scorer.h"

#endif //CLASSIFIER_BASE_H
