//
// Created by Tai, Yuan-yen on 11/11/17.
// All rights reserved.
//

#ifndef LOCAL_GAUSSIAN_MODEL_H
#define LOCAL_GAUSSIAN_MODEL_H

#include "glassfire_base.h"

namespace glassfire 
{

//===================================
//--- Implementation of Classifier
//-----------------------------------
template<typename AxisType, size_t Dimension, typename FeatureInfo>
class Classifier
{
public:
    typedef GlassfireType<AxisType,Dimension,FeatureInfo>   GlassfireType;

    typedef typename GlassfireType::Feature                 Feature;
    typedef typename GlassfireType::Feature_s               Feature_s;

    typedef typename GlassfireType::CentroidRtreeValue      CentroidRtreeValue;

    typedef typename GlassfireType::Rtree                   Rtree;
    typedef typename GlassfireType::RtreeFeature            RtreeFeature;
    typedef typename GlassfireType::RtreeFeature_s          RtreeFeature_s;

    typedef typename GlassfireType::Centroid                Centroid;
    typedef typename GlassfireType::CentroidList            CentroidList;
    typedef typename GlassfireType::CentroidListPtr         CentroidListPtr;
    typedef typename GlassfireType::CentroidListIterator    CentroidListIterator;
    typedef typename GlassfireType::CentroidRtree           CentroidRtree;
    typedef typename GlassfireType::CentroidRtreePtr        CentroidRtreePtr;

    typedef typename GlassfireType::ScorerSet                  ScorerSet;

    typedef unordered_set<size_t>                           CentroidHashSet;
    typedef typename GlassfireType::ClusterModel          ClusterModel;

private:
    Rtree               mRtreeRoot;
    RtreeFeature_s      mRtreeFeature_s;
    CentroidListPtr     mCentroidListPtr;

    CentroidRtreePtr    mCentroidRtreePtr;

public:
    Classifier():
        mRtreeRoot(),
        mRtreeFeature_s()
    { }

    void append_feature(Feature feature, FeatureInfo featureInfo) {
        mRtreeFeature_s.push_back(RtreeFeature(feature));
        mRtreeRoot.insert({mRtreeFeature_s.back(), featureInfo});
    }

    auto run_cluster(AxisType centroid_distance, uint16_t minimal_count = 1, AxisType ratio_of_minimum_diff = 0.01) -> ScorerSet {
        mCentroidListPtr = make_shared<CentroidList>();

        const auto minimum_diff = centroid_distance * ratio_of_minimum_diff;

        // --------------------------------------------------------
        // First part: Generate centroid grid from a given distance
        // --------------------------------------------------------
        CentroidHashSet    centroidHashSet;
        std::hash<string>   stringToHash;
        auto createCentroidFromRtreeFeature = [&](RtreeFeature & rt_feature){
            auto centroidKeyStr = string();
            Feature centroidFeature;

            // Partition the space into small box and assign it to the centroid.
            for(const auto & elem: rt_feature.getFeature()){
                long nKey = 0;
                if(elem >= 0){
                    // Handling positive value
                    nKey = elem / centroid_distance;
                }
                else {
                    // Handling negative value
                    nKey = elem / centroid_distance - 1;
                }
                AxisType centroid_axis = nKey*centroid_distance + 0.5*centroid_distance;
                centroidFeature.push_back(centroid_axis);
                centroidKeyStr += to_string(nKey)+":";
            }
            centroidKeyStr += +":"+to_string(centroid_distance);

            size_t centroidHash = stringToHash(centroidKeyStr);
            if(centroidHashSet.find(centroidHash) == centroidHashSet.end()){
                mCentroidListPtr->push_front(
                        {mRtreeRoot,
                         mRtreeFeature_s,
                         centroidFeature,
                         centroid_distance,
                         centroidKeyStr}
                );

                centroidHashSet.insert(centroidHash);
            }
        };

        for(auto & iter: mRtreeFeature_s){
            createCentroidFromRtreeFeature(iter);
        }

        // ----------------------------------------------------
        // Second part: Optimize centroid
        // ----------------------------------------------------
        CentroidRtree   centroidRtreeRoot;
        auto cleanUpCentroidCollision = [&](bool needCleanTooFewCountCentroid){
            // Clean up centroid collision.
            mCentroidListPtr->sort([](Centroid & L, Centroid & R){return L.count() < R.count();});

            centroidRtreeRoot = CentroidRtree();
            for(auto c_iter = mCentroidListPtr->begin(); c_iter != mCentroidListPtr->end(); c_iter++){
                centroidRtreeRoot.insert({c_iter->getPoint(), c_iter});
            }

            for(auto c_iter = mCentroidListPtr->begin(); c_iter != mCentroidListPtr->end();){
                // If centroid obtain too few data points.
                // Immediately delete.
                if( needCleanTooFewCountCentroid && c_iter->count() <= minimal_count){
                    auto tmp_iter = c_iter;
                    c_iter++;
                    centroidRtreeRoot.remove({tmp_iter->getPoint(), tmp_iter});
                    mCentroidListPtr->erase(tmp_iter);
                    continue;
                }

                vector<CentroidRtreeValue> result_s;
                centroidRtreeRoot.query(
                    bgi::intersects( c_iter->createBox(centroid_distance) ),
                    back_inserter(result_s)
                );

                bool c_iter_advanced = false;
                for(auto rs_iter: result_s){
                    auto cmp = rs_iter.second;
                    if(c_iter != rs_iter.second)
                    if(c_iter->count() <= cmp->count() && c_iter->distance_to(*cmp) < centroid_distance){
                        auto tmp_iter = c_iter;
                        c_iter++;
                        c_iter_advanced = true;
                        centroidRtreeRoot.remove({tmp_iter->getPoint(), tmp_iter});
                        mCentroidListPtr->erase(tmp_iter);
                    }
                }
                if(!c_iter_advanced) c_iter++;
            }
        };

        auto optimizeCentroidPtr = [&](){
            size_t count_needUpdate = 0;
            for(auto & centroid: *mCentroidListPtr){
                if(centroid.needUpdateMean()){
                    count_needUpdate++;
                    centroid.updateMean(minimum_diff);
                }
            }

            return count_needUpdate > 0;
        };

        bool needUpdate;
        size_t iteration=0;
        do{
            needUpdate = optimizeCentroidPtr();
            cleanUpCentroidCollision( iteration > 2);
            iteration++;
        }
        while(needUpdate);

        // ----------------------------------------------------
        // Third part: Calculate covariance matrix
        // ----------------------------------------------------

        for(auto & centroid: *mCentroidListPtr){
            centroid.count_in_range_feature_s(centroidRtreeRoot);
        }
        
        for(auto & centroid: *mCentroidListPtr){
            centroid.updateCovariantMatrix(centroidRtreeRoot, mCentroidListPtr);
        }

        // Build the Centroid Rtree(Ptr) from the Centroid List(Ptr) 
        mCentroidRtreePtr = make_shared<CentroidRtree>();
        for(auto c_iter = mCentroidListPtr->begin(); c_iter != mCentroidListPtr->end(); c_iter++){
            mCentroidRtreePtr->insert({c_iter->getPoint(), c_iter});
        }

        //return mCentroidListPtr->size();
        return ScorerSet(mCentroidListPtr, mCentroidRtreePtr);
    }
};
//-----------------------------------
//--- Implementation of Classifier
//===================================

};

#endif //LOCAL_GAUSSIAN_MODEL_H
