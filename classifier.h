//
// Created by Tai, Yuan-yen on 11/11/17.
// All rights reserved.
//

#ifndef LOCAL_GAUSSIAN_MODEL_H
#define LOCAL_GAUSSIAN_MODEL_H

#include "classifier_base.h"

namespace glowfire 
{
    //===================================
    //--- Implementation of Classifier
    //-----------------------------------
    template<typename AxisType, size_t Dimension>
    class Classifier: public GlowfireBase<AxisType>
    {
    public:
        typedef AxisType                                        VarType;
        typedef typename GlowfireBase<AxisType>::Feature             Feature;
        typedef typename GlowfireBase<AxisType>::Feature_s           Feature_s;

        typedef GlowfireType<AxisType,Dimension>                     GlowfireTypeCollect;
        typedef typename GlowfireTypeCollect::CentroidRtreeValue     CentroidRtreeValue;

        typedef typename GlowfireTypeCollect::Rtree                  Rtree;
        typedef typename GlowfireTypeCollect::RtreeFeature           RtreeFeature;
        typedef typename GlowfireTypeCollect::RtreeFeature_s         RtreeFeature_s;

        typedef typename GlowfireTypeCollect::Centroid               Centroid;
        typedef typename GlowfireTypeCollect::CentroidList           CentroidList;
        typedef typename GlowfireTypeCollect::CentroidListPtr        CentroidListPtr;
        typedef typename GlowfireTypeCollect::CentroidListIterator   CentroidListIterator;
        typedef typename GlowfireTypeCollect::CentroidRtree          CentroidRtree;

        typedef typename GlowfireTypeCollect::Scorer                 Scorer;

        typedef unordered_set<string>                           CentroidStringSet;

    private:
        Rtree               mRtreeRoot;
        RtreeFeature_s      mRtreeFeature_s;
        CentroidListPtr     mCentroidListPtr;

    public:
        Classifier():
            mRtreeRoot(),
            mRtreeFeature_s()
        { }

        void append_feature(Feature feature) {
            mRtreeFeature_s.push_back(RtreeFeature(feature));
            mRtreeRoot.insert({mRtreeFeature_s.back(), mRtreeRoot.size()});
        }

        auto run_cluster(AxisType centroid_distance, uint16_t minimal_count = 1, AxisType ratio_of_minimum_diff = 0.01) -> size_t {
            mCentroidListPtr = make_shared<CentroidList>();

            const auto minimum_diff = centroid_distance * ratio_of_minimum_diff;

            // --------------------------------------------------------
            // First part: Generate centroid grid from a given distance
            // --------------------------------------------------------
            CentroidStringSet    centroidStringSet;
            auto createCentroidFromRtreeFeature = [&](RtreeFeature & rt_feature){
                // TODO: Use hash id to replace the string as look up key.
                auto centroidKey = string();
                Feature centroidFeature;

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
                    centroidKey += to_string(nKey)+" ";
                }

                if(centroidStringSet.find(centroidKey) == centroidStringSet.end()){
                    mCentroidListPtr->push_front(
                            {mRtreeRoot,
                             mRtreeFeature_s,
                             centroidFeature,
                             centroid_distance}
                    );

                    centroidStringSet.insert(centroidKey);
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
                centroid.updateCovariantMatrix(centroidRtreeRoot);
            }

            return mCentroidListPtr->size();
        }

        auto calc_score(const Feature & feature) -> map<AxisType, size_t, std::greater<AxisType>> {
            // This method will return a descending order map
            // Where this map contain a [key: value] pair of [calculated-score: cluster-id]
            // The cluster-id are all positive integers

            auto scorer = Scorer(mCentroidListPtr);
            return scorer.calc_score(feature);
        }

        auto create_scorer() -> Scorer { return Scorer(mCentroidListPtr); }
    };
    //-----------------------------------
    //--- Implementation of Classifier
    //===================================
};

#endif //LOCAL_GAUSSIAN_MODEL_H
