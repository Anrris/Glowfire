#ifndef SCORER_H
#define SCORER_H

#include <tuple>

#include "glassfire_base.h"
#include "glassfire_cluster_model.h"

namespace glassfire{

    template<typename AxisType, typename FeatureInfo>
    class ScorerSetBase{
    public:
        typedef std::vector<AxisType> Feature;

        virtual ~ScorerSetBase(){}
        virtual auto calc_scores(const Feature & feature)
                        -> std::map<AxisType, std::string, std::greater<AxisType>>  =0;
        virtual auto get_model_set(AxisType regularize)
                        -> std::vector<ClusterModel<AxisType, FeatureInfo>> =0;

        virtual auto query(const Feature & feature, AxisType regularize, int nearest_count = -1)
                        -> std::tuple<bool, AxisType, glassfire::ClusterModel<AxisType, FeatureInfo>, std::string> = 0;
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

        auto get_model_set(AxisType regularize) -> std::vector<ClusterModel> {
            std::vector<ClusterModel> result;
            for(auto & iter: *mCentroidListPtr){
                auto model = iter.get_model(regularize);
                result.push_back(model);
            }
            return result;
        }

        auto query(
            const Feature & feature, AxisType regularize, int nearest_count = -1
            ) -> std::tuple<bool, AxisType, glassfire::ClusterModel<AxisType, FeatureInfo>, std::string> {

            auto rtreeFeature = RtreeFeature(feature);

            if (nearest_count == -1){
                nearest_count = 2 * Dimension * ceil(sqrt((AxisType)Dimension));
            }
            std::vector<CentroidRtreeValue> result_s;
            mCentroidRtreePtr->query(
                    bgi::nearest((RtreePoint)rtreeFeature, nearest_count), back_inserter(result_s)
            );

            auto centroidMap = std::map<AxisType, ClusterModel, std::greater<AxisType>>();
            for(auto & crv_iter: result_s){
                auto model = crv_iter.second->get_model(regularize);
                //auto score = crv_iter.second->scoreOfFeature(feature);
                auto score = model.eval(feature);
                centroidMap.insert({ score, model });
            }
            if(centroidMap.size()==0){
                return std::make_tuple(false, -1, glassfire::ClusterModel<AxisType,FeatureInfo>(), "Cluster not in range!");
            }

            return std::make_tuple(true, centroidMap.begin()->first, centroidMap.begin()->second, "Cluster has been found.");
        }

        auto cluser_count() -> size_t { return mCentroidListPtr->size(); }
    };
    //-----------------------------------
    //--- Implementation of ScorerSet --
    //===================================
}

#endif //CLASSIFIER_BASE_H