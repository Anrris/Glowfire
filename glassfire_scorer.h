
#ifndef SCORER_H
#define SCORER_H

#include "glassfire_base.h"
#include "glassfire_cluster_model.h"

namespace glassfire{

    template<typename AxisType, typename FeatureInfo>
    class ScorerSetBase{
    public:
        typedef std::vector<AxisType> Feature;

        virtual ~ScorerSetBase(){}
        virtual auto calc_scores(const Feature & feature) -> std::map<AxisType, std::string, std::greater<AxisType>>  =0;
        virtual auto get_model_set() -> std::vector<ClusterModel<AxisType, FeatureInfo>> =0;
        virtual auto query_centroids(const Feature & feature, AxisType box_distance) -> std::pair<AxisType, glassfire::ClusterModel<AxisType, FeatureInfo>> = 0;

    };
    //===================================
    //--- Implementation of ScorerSet --
    //-----------------------------------
    template<typename AxisType, size_t Dimension, typename FeatureInfo>
    class GlassfireType<AxisType, Dimension, FeatureInfo>::ScorerSet: public ScorerSetBase<AxisType, FeatureInfo>{ 
            CentroidListPtr mCentroidListPtr;
            CentroidRtreePtr mCentroidRtreePtr;
        public:
        ScorerSet(){}
        ScorerSet(CentroidListPtr centroidListPtr, CentroidRtreePtr centroidRtreePtr):
            mCentroidListPtr(centroidListPtr),
            mCentroidRtreePtr(centroidRtreePtr)
        {}

        auto calc_scores(const Feature & feature) -> std::map<AxisType, std::string, std::greater<AxisType>> {
            // This method will return a descending order map
            // Where this map contain a [key: value] pair of [calculated-score: cluster-id]
            // The cluster-id are all positive integers

            auto retval = std::map<AxisType, std::string, std::greater<AxisType>>();

            size_t centroid_index = 0;
            for(auto & iter: *mCentroidListPtr){
                centroid_index++;

                auto score = iter.scoreOfFeature(feature);
                retval[score] = iter.getKeyStr();
            }
            return retval;
        }

        auto get_centroids() -> std::vector<Centroid> {
            std::vector<Centroid> retval;
            for(const auto & iter: *mCentroidListPtr){
                retval.push_back(iter);
            }
            return retval;
        }

        auto get_model_set() -> std::vector<ClusterModel> {
            std::vector<ClusterModel> retval;
            for(auto & iter: *mCentroidListPtr){
                retval.push_back(iter.get_model());
            }
            return retval;
        }

        auto query_centroids(const Feature & feature, AxisType box_distance) -> std::pair<AxisType, glassfire::ClusterModel<AxisType, FeatureInfo>> {

            auto rtreeFeature = RtreeFeature(feature);

            std::vector<CentroidRtreeValue> result_s;
            mCentroidRtreePtr->query(
                    bgi::intersects(rtreeFeature.createBox(box_distance)), back_inserter(result_s)
            );

            auto centroidMap = std::map<AxisType, ClusterModel, std::greater<AxisType>>();
            for(auto & crv_iter: result_s){
                auto score = crv_iter.second->scoreOfFeature(feature);
                centroidMap.insert({ score, crv_iter.second->get_model() });
            }

            return *centroidMap.begin();
        }

        auto cluser_count() -> size_t { return mCentroidListPtr->size(); }
    };
    //-----------------------------------
    //--- Implementation of ScorerSet --
    //===================================
}

#endif //CLASSIFIER_BASE_H