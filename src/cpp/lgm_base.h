//
// Created by Tai, Yuan-yen on 11/17/17.
// All rights reserved.
//

#ifndef LOCALGAUSSIANMODEL_LGM_BASE_H
#define LOCALGAUSSIANMODEL_LGM_BASE_H

#include <string>                               /* STL Library*/
#include <sstream>                              /* STL Library*/
#include <iostream>                             /* STL Library*/
#include <vector>                               /* STL Library*/
#include <list>                                 /* STL Library*/
#include <unordered_set>                        /* STL Library*/
#include <memory>                               /* STL Library*/
#include <limits>                               /* STL Library*/
#include <exception>                            /* STL Library*/

#include <boost/geometry.hpp>                   /* Boost Library: For RTree*/
#include <boost/geometry/geometries/point.hpp>  /* Boost Library */
#include <boost/geometry/geometries/box.hpp>    /* Boost Library */
#include <boost/geometry/index/rtree.hpp>       /* Boost Library */

#include <Eigen/Dense>                          /* Eigen Library: For matrix operation*/

namespace LGM{

    using namespace std;
    namespace bg  = boost::geometry;
    namespace bgm = boost::geometry::model;
    namespace bgi = boost::geometry::index;

    template<typename AxisType>
    class LgmBase
    {
    public:
        typedef shared_ptr<LgmBase<AxisType>>   Ptr;
        typedef vector<AxisType>                Feature;
        typedef vector<Feature>                 Feature_s;

        virtual void append_feature(Feature) = 0;
        virtual auto run_cluster(AxisType _centroid_distance, AxisType ratio_of_minimum_diff = 0.01) -> size_t = 0;
        virtual auto calc_score(const Feature & feature) -> map<AxisType, size_t, std::greater<AxisType>> = 0;
    };

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
    template<typename AxisType, size_t Dimension>
    class LgmType
    {
    public:
        typedef typename LgmBase<AxisType>::Feature                 Feature;
        typedef typename LgmBase<AxisType>::Feature_s               Feature_s;

        typedef bgm::point<AxisType, Dimension, bg::cs::cartesian>  RtreePoint;
        typedef bgm::box<RtreePoint>                                RtreeBox;
        typedef pair<RtreePoint, size_t>                            RtreeValue;
        typedef bgi::rtree<RtreeValue, bgi::linear<16>>             Rtree;

        typedef Eigen::Matrix<AxisType, Dimension, Dimension>       LgmMatrix;
        typedef Eigen::Matrix<AxisType, 1, Dimension>               LgmRowVector;
        typedef Eigen::Matrix<AxisType, Dimension, 1>               LgmColVector;

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

        //===================================
        //--- Implementation of Centroid ----
        //-----------------------------------
        class Centroid: public RtreeFeature
        {
        public:
            typedef list<Centroid>                                  CentroidList;
            typedef typename CentroidList::iterator                 CentroidListIterator;
            typedef pair<RtreePoint, CentroidListIterator>          CentroidRtreeValue;
            typedef bgi::rtree<CentroidRtreeValue, bgi::linear<16>> CentroidRtree;

        private:
            Rtree &             mRtree_ref;
            RtreeFeature_s &    mRtreeFeature_s_ref;
            AxisType            centroid_distance;

            AxisType        m_diff = numeric_limits<AxisType>::max();
            AxisType        m_minimum_difference;

            LgmMatrix       mCmat;
            LgmMatrix       mInvCmat;
            AxisType        mCmatDet;

            size_t          m_occupation_count = 0;

            const AxisType  _pi_ = 3.141592653589793;
        public:
            Centroid(
                Rtree &             _Rtree_ref,
                RtreeFeature_s&     _RtreeFeature_s_ref,
                Feature             _Mean,
                AxisType            _centroid_distance
            ):
            RtreeFeature        (_Mean),
            mRtree_ref          (_Rtree_ref),
            mRtreeFeature_s_ref (_RtreeFeature_s_ref),
            centroid_distance   (_centroid_distance)
            {
                // Zerolize mCmat
                for (int idx_r = 0; idx_r < mCmat.rows(); idx_r++)
                for (int idx_c = 0; idx_c < mCmat.cols(); idx_c++)
                    if(idx_r != idx_c) mCmat(idx_r, idx_c) = 0;
            }

            auto distance_to(RtreeFeature & other) -> AxisType {
                AxisType total_del_retval = 0;
                for(size_t i=0; i<this->getFeature().size(); i++){
                    AxisType del = this->at(i) - other.at(i);
                    total_del_retval += del * del;
                }
                return sqrt(total_del_retval);
            }
            auto distance_to(Feature & feature) -> AxisType {
                AxisType total_del_retval = 0;
                for(size_t i=0; i<feature.size(); i++){
                    AxisType del = this->at(i) - feature[i];
                    total_del_retval += del * del;
                }
                return sqrt(total_del_retval);
            }

            auto updateMean(AxisType _minimum_difference) -> AxisType {
                m_minimum_difference = _minimum_difference;

                AxisType max_diff = 0;

                vector<RtreeValue> result_s;
                mRtree_ref.query(
                        bgi::intersects(this->createBox(centroid_distance)), back_inserter(result_s)
                );

                m_occupation_count = result_s.size();

                Feature old_mean = this->getFeature();
                Feature new_mean(Dimension, 0);

                // Calculate mean from all points within region.
                for(auto & vpair: result_s){
                    size_t index = vpair.second;
                    for(size_t i=0 ; i<new_mean.size() ; i++){
                        new_mean[i] += mRtreeFeature_s_ref[index].at(i) / result_s.size();
                    }
                }

                m_diff = distance_to(new_mean);

                this->setFeature(new_mean);

                return m_diff;
            }

            auto count() -> size_t {return m_occupation_count;}

            auto needUpdateMean() -> bool {
                return m_diff > m_minimum_difference;
            }

            auto printMean() -> string {
                stringstream ss;
                for(size_t i=0 ; i<this->getFeature().size(); i++){
                    ss << this->at(i) << " ";
                }
                return ss.str();
            }

            auto updateCovariantMatrix(CentroidRtree & centroidRtreeRef) -> AxisType {
                // TODO: using nearest neighbor distance for centroid mCmat update.
                vector<CentroidRtreeValue> nearest_result;
                centroidRtreeRef.query(bgi::nearest((RtreePoint)*this, 2), std::back_inserter(nearest_result));

                AxisType nearest_neighbor_distance = this->distance_to(
                        *nearest_result.front().second
                );

                vector<RtreeValue> result_s;
                mRtree_ref.query(
                        bgi::intersects(this->createBox(nearest_neighbor_distance)), back_inserter(result_s)
                );

                LgmMatrix tmpCmat;
                bool is_fresh_start = true;

                auto reCalcCovMat = [&]() {
                    // Zerolize tmpCmat
                    for (int idx_r = 0; idx_r < tmpCmat.rows(); idx_r++)
                    for (int idx_c = 0; idx_c < tmpCmat.cols(); idx_c++)
                            tmpCmat(idx_r, idx_c) = 0;

                    AxisType p_denumerator = 0;

                    list<AxisType> weightList;
                    for (auto &it : result_s) {
                        auto & feature = mRtreeFeature_s_ref[it.second].getFeature();
                        if(is_fresh_start){
                            weightList.push_back(1.0);
                        }
                        else{
                            weightList.push_back(scoreOfFeature(feature));
                        }
                        p_denumerator += weightList.back();
                    }
                    for (auto &it : result_s) {
                        auto & feature = mRtreeFeature_s_ref[it.second].getFeature();
                        auto & rt_feature = mRtreeFeature_s_ref[it.second];
                        if (this->distance_to(rt_feature) > centroid_distance ) continue;

                        auto colMeanVec = feature_sub_mean_to_colVector(feature);
                        auto weight = weightList.front()/p_denumerator;
                        weightList.pop_front();

                        for (int idx_r = 0; idx_r < tmpCmat.rows(); idx_r++)
                        for (int idx_c = 0; idx_c < tmpCmat.rows(); idx_c++) {
                            tmpCmat(idx_r, idx_c) += weight
                                * colMeanVec(idx_r, 0)
                                * colMeanVec(idx_c, 0);
                        }
                    }

                };


                AxisType max_diff = numeric_limits<AxisType>::max();
                do{
                    reCalcCovMat();

                    LgmMatrix diffCmat = mCmat - tmpCmat;

                    if( is_fresh_start ){
                        mCmat = tmpCmat;
                    }
                    else{
                        mCmat = 0.5 * (mCmat + tmpCmat);
                    }
                    mInvCmat = tmpCmat.inverse();
                    mCmatDet = tmpCmat.determinant();

                    if( is_fresh_start ){
                        is_fresh_start = false;
                        continue;
                    }

                    max_diff = 0;
                    for (int idx_r = 0; idx_r < tmpCmat.rows(); idx_r++)
                    for (int idx_c = 0; idx_c < tmpCmat.cols(); idx_c++) {
                        max_diff = max(max_diff, abs(diffCmat(idx_r, idx_c)) );
                    }
                } while (max_diff > m_minimum_difference);

                return max_diff;
            }

            auto scoreOfFeature(const Feature &feature) -> AxisType {
                auto rowVec = feature_sub_mean_to_rowVector(feature);
                auto colVec = feature_sub_mean_to_colVector(feature);
                AxisType mahalanDistance = (rowVec * mInvCmat * colVec)(0, 0);
                AxisType result = (AxisType)exp( -0.5* mahalanDistance) / sqrt( 2 * _pi_ * mCmatDet );
                return result;
            }

            auto getOccupationCount() -> size_t { return m_occupation_count; }

            auto getCovariantMatrix() -> const LgmMatrix & { return mCmat; }

        private:
            auto feature_sub_mean_to_rowVector(const Feature &feature) -> LgmRowVector {
                LgmRowVector rowVector;
                for(size_t i=0; i<feature.size(); i++){
                    rowVector(0,i) = feature[i] - this->at(i);
                }
                return rowVector;
            }
            auto feature_sub_mean_to_colVector(const Feature &feature) -> LgmColVector {
                LgmColVector colVector;
                for(size_t i=0; i<feature.size(); i++){
                    colVector(i,0) = feature[i] - this->at(i);
                }
                return colVector;
            }
        };
        //-----------------------------------
        //--- Implementation of Centroid ----
        //===================================

        typedef typename Centroid::CentroidList         CentroidList;
        typedef shared_ptr<CentroidList>                CentroidListPtr;
        typedef typename Centroid::CentroidListIterator CentroidListIterator;
        typedef typename Centroid::CentroidRtreeValue   CentroidRtreeValue;
        typedef typename Centroid::CentroidRtree        CentroidRtree;

        //===================================
        //--- Implementation of Scorer ----
        //-----------------------------------
        class Scorer{
            CentroidListPtr mCentroidListPtr;
        public:
            Scorer(){}
            Scorer(CentroidListPtr centroidListPtr):
                mCentroidListPtr(centroidListPtr)
            {}
            auto calc_score(const Feature & feature) -> map<AxisType, size_t, std::greater<AxisType>>
            {
                // This method will return a descending order map
                // Where this map contain a [key: value] pair of [calculated-score: cluster-id]
                // The cluster-id are all positive integers

                auto retval = map<AxisType, size_t, std::greater<AxisType>>();

                size_t centroid_index = 0;
                for(auto & iter: *mCentroidListPtr){
                    auto score = iter.scoreOfFeature(feature);
                    retval[score] = centroid_index;
                    centroid_index++;
                }
                return retval;
            }

            auto get_centroids() -> vector<Centroid> {
                vector<Centroid> retval;
                for(const auto & iter: *mCentroidListPtr){
                    retval.push_back(iter);
                }
                return retval;
            }

            auto cluser_count() -> size_t { return mCentroidListPtr->size(); }
        };
        //-----------------------------------
        //--- Implementation of Scorer ----
        //===================================
    };
}

#endif //LOCALGAUSSIANMODEL_LGM_BASE_H
