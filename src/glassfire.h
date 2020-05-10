//
// Created by Tai, Yuan-yen on 11/11/17.
// All rights reserved.
//

#ifndef LOCAL_GAUSSIAN_MODEL_H
#define LOCAL_GAUSSIAN_MODEL_H

#include "glassfire_base.h"
#include "glassfire_util.h"

namespace glassfire 
{

template<typename AxisType, typename FeatureInfo>
class ClassifierBase{
    typedef vector<AxisType> Feature;
public:
    typedef std::vector<std::tuple<size_t, AxisType, FeatureInfo, std::string, Feature>> DataQueryReturnType;

    virtual ~ClassifierBase(){}
    virtual void append_feature(Feature feature, FeatureInfo featureInfo) =0 ;
    virtual auto run_cluster(
                    AxisType centroid_distance,
                    uint16_t minimal_count = 1,
                    AxisType ratio_of_minimum_diff = 0.01,
                    bool     using_centroid_distance = true
                    ) -> std::shared_ptr<ScorerSetBase<AxisType, FeatureInfo>> =0 ;

    virtual auto query_data(ClusterModel<AxisType, FeatureInfo> & cm, AxisType box_size) -> DataQueryReturnType =0 ;
};

//===================================
//--- Implementation of Classifier
//-----------------------------------
template<typename AxisType, size_t Dimension, typename FeatureInfo>
class Classifier: public ClassifierBase<AxisType, FeatureInfo>
{
public:
    typedef GlassfireType<AxisType,Dimension,FeatureInfo>   GlassfireType;

    typedef typename ClassifierBase<AxisType, FeatureInfo>::DataQueryReturnType DataQueryReturnType;

    typedef typename GlassfireType::Feature                 Feature;
    typedef typename GlassfireType::Feature_s               Feature_s;

    typedef typename GlassfireType::CentroidRtreeValue      CentroidRtreeValue;

    typedef typename GlassfireType::Rtree                   Rtree;
    typedef typename GlassfireType::RtreeValue              RtreeValue;
    typedef typename GlassfireType::RtreePoint              RtreePoint;
    typedef typename GlassfireType::RtreeFeature            RtreeFeature;
    typedef typename GlassfireType::RtreeFeature_s          RtreeFeature_s;

    typedef typename GlassfireType::Centroid                Centroid;
    typedef typename GlassfireType::CentroidList            CentroidList;
    typedef typename GlassfireType::CentroidListPtr         CentroidListPtr;
    typedef typename GlassfireType::CentroidListIterator    CentroidListIterator;
    typedef typename GlassfireType::CentroidRtree           CentroidRtree;
    typedef typename GlassfireType::CentroidRtreePtr        CentroidRtreePtr;

    typedef typename GlassfireType::ScorerSet               ScorerSet;

    typedef std::unordered_set<size_t>                      CentroidHashSet;

private:
    //using ClusterModel = ClusterModel<AxisType, FeatureInfo>;
    Rtree               mRtreeRoot;
    RtreeFeature_s      mRtreeFeature_s;
    CentroidListPtr     mCentroidListPtr;

    CentroidRtreePtr    mCentroidRtreePtr;

    bool    m_cluster_state_successfull = false;
public:
    Classifier():
        mRtreeRoot(),
        mRtreeFeature_s()
    { }

    void append_feature(Feature feature, FeatureInfo featureInfo) {
        auto tmp_RtreeFeature = RtreeFeature(feature);
        tmp_RtreeFeature.setInfo(featureInfo);
        mRtreeFeature_s.push_back(tmp_RtreeFeature);
        mRtreeRoot.insert({mRtreeFeature_s.back(), mRtreeRoot.size()});
    }

    auto query_data(ClusterModel<AxisType, FeatureInfo> & cm, AxisType box_size) -> DataQueryReturnType{

        RtreeFeature rtfeat(cm.mean());

        std::vector<RtreeValue> result_s;
        mRtreeRoot.query(
                bgi::intersects(rtfeat.createBox(box_size)), back_inserter(result_s)
        );

        DataQueryReturnType result_collections;
        for(auto &iter: result_s){
            result_collections.push_back(
                std::make_tuple(
                    iter.second,
                    cm.eval(mRtreeFeature_s[iter.second].getFeature()),
                    mRtreeFeature_s[iter.second].getInfo(),
                    cm.model_key(),
                    mRtreeFeature_s[iter.second].getFeature()
                )
            );
        }
        return result_collections;
    }

    auto run_cluster(
        AxisType centroid_distance,
        uint16_t minimal_count = 1,
        AxisType ratio_of_minimum_diff = 0.01,
        bool     using_centroid_distance = true
        ) -> std::shared_ptr<ScorerSetBase<AxisType, FeatureInfo>> {

        // --------------------------------------------------------
        // --- If the entire process go through the end------------
        // ----m_cluster_state_successfull will set to be true. ---
        // --------------------------------------------------------
        m_cluster_state_successfull = false;
        // --------------------------------------------------------
        // --------------------------------------------------------
        // --------------------------------------------------------

        mCentroidListPtr = std::make_shared<CentroidList>();

        const auto minimum_diff = centroid_distance * ratio_of_minimum_diff;

        // --------------------------------------------------------
        // First part: Generate centroid grid from a given distance
        // --------------------------------------------------------
        CentroidHashSet    centroidHashSet;
        std::hash<std::string>   stringToHash;
        auto createCentroidFromRtreeFeature = [&](RtreeFeature & rt_feature){
            auto centroidKeyStr = std::string();
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
                centroidKeyStr += fmt_string(nKey, true)+":";
            }
            for(size_t i=0; i<Dimension; i++){
                centroidKeyStr += +":"+fmt_string(centroid_distance);
            }

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
                if( needCleanTooFewCountCentroid && c_iter->count() < minimal_count){
                    auto tmp_iter = c_iter;
                    c_iter++;
                    centroidRtreeRoot.remove({tmp_iter->getPoint(), tmp_iter});
                    mCentroidListPtr->erase(tmp_iter);
                    continue;
                }

                std::vector<CentroidRtreeValue> result_s;
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
            centroid.count_in_range_feature_s(centroidRtreeRoot, using_centroid_distance);
        }
        
        for(auto & centroid: *mCentroidListPtr){
            centroid.updateCovariantMatrix(centroidRtreeRoot, mCentroidListPtr);
        }

        // Build the Centroid Rtree(Ptr) from the Centroid List(Ptr) 
        mCentroidRtreePtr = std::make_shared<CentroidRtree>();
        for(auto c_iter = mCentroidListPtr->begin(); c_iter != mCentroidListPtr->end(); c_iter++){
            mCentroidRtreePtr->insert({c_iter->getPoint(), c_iter});
        }

        // --------------------------------------------------------
        // --------------------------------------------------------
        // --------------------------------------------------------
        m_cluster_state_successfull = true;
        // --------------------------------------------------------
        // ---Ending the cluster analysis, setting up--------------
        // --- m_cluster_state_successfull true. ------------------
        // --------------------------------------------------------

        return std::make_shared<ScorerSet>(ScorerSet(mCentroidListPtr, mCentroidRtreePtr));
    }
};
////-----------------------------------
////--- Implementation of Classifier
////===================================
//
template<typename AxisType, typename FeatureInfo, size_t Dimension>
class ClassifierFactory{
    typedef std::shared_ptr<ClassifierBase<AxisType, FeatureInfo>> ClassifierBasePtr;

    template<size_t x, size_t to>
    struct t_static_for {
        void create(std::vector<ClassifierBasePtr> &classifierBaseSet) {
            classifierBaseSet.push_back(std::make_shared<Classifier<AxisType,x,FeatureInfo>>());
            t_static_for<x+1,to>().create(classifierBaseSet);
        }
    };

    template<size_t to>
    struct t_static_for<to, to> {
        void create(std::vector<ClassifierBasePtr> &classifierBaseSet) {}
    };

    typedef t_static_for<1, Dimension+1>    static_for;

public:
    auto create() -> std::vector<ClassifierBasePtr> {
        std::vector<ClassifierBasePtr> classifierBaseSet;

        static_for  m_static_for;
        m_static_for.create(classifierBaseSet);
        return classifierBaseSet;
    }

};



};

#endif //LOCAL_GAUSSIAN_MODEL_H
