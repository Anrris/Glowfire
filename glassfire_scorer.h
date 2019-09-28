
#ifndef SCORER_H
#define SCORER_H

#include "glassfire_base.h"

namespace glassfire{
    //===================================
    //--- Implementation of Scorer --
    //-----------------------------------
    template<typename AxisType, size_t Dimension, typename FeatureInfo>
    class GlassfireType<AxisType, Dimension, FeatureInfo>::Scorer{ 
            CentroidListPtr mCentroidListPtr;
            CentroidRtreePtr mCentroidRtreePtr;
        public:
        Scorer(){}
        Scorer(CentroidListPtr centroidListPtr, CentroidRtreePtr centroidRtreePtr):
            mCentroidListPtr(centroidListPtr),
            mCentroidRtreePtr(centroidRtreePtr)
        {}

        auto calc_scores(const Feature & feature) -> map<AxisType, string, std::greater<AxisType>> {
            // This method will return a descending order map
            // Where this map contain a [key: value] pair of [calculated-score: cluster-id]
            // The cluster-id are all positive integers

            auto retval = map<AxisType, string, std::greater<AxisType>>();

            size_t centroid_index = 0;
            for(auto & iter: *mCentroidListPtr){
                centroid_index++;

                auto score = iter.scoreOfFeature(feature);
                retval[score] = iter.getKeyStr();
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

        auto query_centroids(const Feature & feature, AxisType box_distance) -> pair<AxisType, ClusterModel> {

            auto rtreeFeature = RtreeFeature(feature);

            vector<CentroidRtreeValue> result_s;
            mCentroidRtreePtr->query(
                    bgi::intersects(rtreeFeature.createBox(box_distance)), back_inserter(result_s)
            );

            auto centroidMap = map<AxisType, ClusterModel, std::greater<AxisType>>();
            for(auto & crv_iter: result_s){
                auto score = crv_iter.second->scoreOfFeature(feature);
                centroidMap.insert({ score, crv_iter.second->get_model() });
            }

            return *centroidMap.begin();
        }

        auto cluser_count() -> size_t { return mCentroidListPtr->size(); }
    };
    //-----------------------------------
    //--- Implementation of Scorer --
    //===================================
}

#endif //CLASSIFIER_BASE_H