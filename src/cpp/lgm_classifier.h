//
// Created by Tai, Yuan-yen on 11/11/17.
// All rights reserved.
//

#ifndef LOCAL_GAUSSIAN_MODEL_H
#define LOCAL_GAUSSIAN_MODEL_H


#include "lgm_base.h"

namespace LGM
{
    //===================================
    //--- Implementation of LgmClassifier
    //-----------------------------------
    template<typename AxisType, size_t Dimension>
    class LgmClassifier: public LgmBase<AxisType>
    {
    public:
        typedef typename LgmBase<AxisType>::Feature             Feature;
        typedef typename LgmBase<AxisType>::Feature_s           Feature_s;

        typedef LgmType<AxisType,Dimension>                     LgmTypeCollect;
        typedef typename LgmTypeCollect::CentroidRtreeValue     CentroidRtreeValue;

        typedef typename LgmTypeCollect::Rtree                  Rtree;
        typedef typename LgmTypeCollect::RtreeFeature           RtreeFeature;
        typedef typename LgmTypeCollect::RtreeFeature_s         RtreeFeature_s;

        typedef typename LgmTypeCollect::Centroid               Centroid;
        typedef typename LgmTypeCollect::CentroidList           CentroidList;
        typedef typename LgmTypeCollect::CentroidListIterator   CentroidListIterator;
        typedef typename LgmTypeCollect::CentroidRtree          CentroidRtree;

        typedef unordered_set<string>                           CentroidStringSet;

    private:
        Rtree               mRtreeRoot;
        RtreeFeature_s   mRtreeRtreeFeature_s;
        CentroidList        mCentroidList;

        AxisType            centroid_distance;

    public:
        LgmClassifier() = default;

        void append_feature(Feature feature) {
            mRtreeRtreeFeature_s.push_back(RtreeFeature(feature));
            mRtreeRoot.insert({mRtreeRtreeFeature_s.back(), mRtreeRoot.size()});
        }

        auto run_cluster(AxisType _centroid_distance, AxisType ratio_of_minimum_diff = 0.01) -> size_t {
            centroid_distance = _centroid_distance;
            mCentroidList.clear();

            const auto minimum_diff = centroid_distance * ratio_of_minimum_diff;

            // --------------------------------------------------------
            // First part: Generate centroid grid from a given distance
            // --------------------------------------------------------
            CentroidStringSet    centroidStringSet;
            auto createCentroidFromRtreeFeature = [&](RtreeFeature & rt_feature){
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
                    mCentroidList.push_front(
                            {mRtreeRoot,
                             mRtreeRtreeFeature_s,
                             centroidFeature,
                             centroid_distance
                            }
                    );

                    centroidStringSet.insert(centroidKey);
                }
            };

            for(auto & iter: mRtreeRtreeFeature_s){
                createCentroidFromRtreeFeature(iter);
            }

            // ----------------------------------------------------
            // Second part: Optimize centroid
            // ----------------------------------------------------
            auto cleanUpCentroidCollision = [&](){
                mCentroidList.sort([](Centroid & L, Centroid & R){return L.count() > R.count();});

                CentroidRtree   centroidRtreeRoot;

                for(auto c_iter = mCentroidList.begin(); c_iter != mCentroidList.end(); c_iter++){
                    centroidRtreeRoot.insert({c_iter->getPoint(), c_iter});
                }
                for(auto c_iter = mCentroidList.begin(); c_iter != mCentroidList.end(); c_iter++){
                    vector<CentroidRtreeValue> result_s;
                    centroidRtreeRoot.query(
                        bgi::intersects(
                            c_iter->createBox(centroid_distance*2)
                        ),
                        back_inserter(result_s)
                    );

                    for(auto rs_iter: result_s){
                        auto cmp = rs_iter.second;
                        if(c_iter->count() >= cmp->count() && c_iter != rs_iter.second){
                            centroidRtreeRoot.remove({cmp->getPoint(), cmp});
                            mCentroidList.erase(cmp);
                        }
                    }
                }

            };
            auto optimizeCentroidMean = [&](){
                bool needUpdate = true;
                double diff = 0;
                for(auto & centroid: mCentroidList){
                    diff = centroid.updateMean(minimum_diff);
                    needUpdate = needUpdate && centroid.needUpdateMean();
                }

                cleanUpCentroidCollision();

                return needUpdate;
            };

            while(optimizeCentroidMean()){ }

            // ----------------------------------------------------
            // Third part: Calculate covariance matrix
            // ----------------------------------------------------
            for(auto & centroid: mCentroidList){
                centroid.updateCovariantMatrix();
            }

            return mCentroidList.size();
        }

        auto calc_score(const Feature & feature) -> map<AxisType, size_t, std::greater<AxisType>> {
            // This method will return an descending order map
            // Where this map contain a [key: value] pair of [calculated-score: cluster-id]
            // The cluster-id are all positive integers

            auto retval = map<AxisType, size_t, std::greater<AxisType>>();

            size_t centroid_index = 0;
            for(auto & iter: mCentroidList){
                auto score = iter.scoreOfFeature(feature);
                retval[score] = centroid_index;
                centroid_index++;
            }
            return retval;
        }
    };
    //-----------------------------------
    //--- Implementation of LgmClassifier
    //===================================
};

#endif //LOCAL_GAUSSIAN_MODEL_H
