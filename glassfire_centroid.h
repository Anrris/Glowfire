
#ifndef CENTROID_H
#define CENTROID_H

#include "glassfire_base.h"

namespace glassfire{

//===================================
//--- Implementation of Centroid ----
//-----------------------------------
template<typename AxisType, size_t Dimension, typename FeatureInfo>
class GlassfireType<AxisType, Dimension, FeatureInfo>::Centroid: public RtreeFeature
{
public:
    typedef list<Centroid>                                  CentroidList;
    typedef shared_ptr<CentroidList>                        CentroidListPtr;
    typedef typename CentroidList::iterator                 CentroidListIterator;
    typedef pair<RtreePoint, CentroidListIterator>          CentroidRtreeValue;
    typedef bgi::rtree<CentroidRtreeValue, bgi::linear<16>> CentroidRtree;

private:
    Rtree &             mRtree_ref;
    RtreeFeature_s &    mRtreeFeature_s_ref;
    AxisType            centroid_distance;
    string              mCentroidKeyStr;

    AxisType            N_soft_k;
    AxisType            N_k;

    Matrix              mCmat;
    Matrix              mInvCmat;
    AxisType            mCmatDet;
    AxisType            m_diff = numeric_limits<AxisType>::max();
    AxisType            m_minimum_difference;

    size_t              m_occupation_count = 0;

    Feature_s           m_in_range_feature_s;

    const AxisType  _pi_ = 3.141592653589793;
public:
    Centroid(
        Rtree &             _Rtree_ref,
        RtreeFeature_s&     _RtreeFeature_s_ref,
        Feature             _Mean,
        AxisType            _centroid_distance,
        string              _centroid_key_str
        ):      
            RtreeFeature        (_Mean),
            mRtree_ref          (_Rtree_ref),
            mRtreeFeature_s_ref (_RtreeFeature_s_ref),
            centroid_distance   (_centroid_distance),
            mCentroidKeyStr     (_centroid_key_str)
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
    auto getKeyStr() -> string {return mCentroidKeyStr;}
    auto needUpdateMean() -> bool { return m_diff > m_minimum_difference; }
    auto printMean() -> string {
        stringstream ss;
        for(size_t i=0 ; i<this->getFeature().size(); i++){
            ss << this->at(i) << " ";
        }
        return ss.str();
    }

    auto scoreOfFeature(const Feature &feature) -> AxisType {
        auto rowVec = feature_sub_mean_to_rowVector(feature);
        auto colVec = feature_sub_mean_to_colVector(feature);
        AxisType mahalanDistance = (rowVec * mInvCmat * colVec)(0, 0);
        AxisType result = (AxisType)exp( -0.5* mahalanDistance) / sqrt( pow(2 * _pi_, Dimension) * mCmatDet);
        return result;
    }

    auto count_in_range_feature_s(CentroidRtree & centroidRtreeRef) -> void {
        vector<CentroidRtreeValue> nearest_result;
        centroidRtreeRef.query(bgi::nearest((RtreePoint)*this, 2), std::back_inserter(nearest_result));

        AxisType nearest_neighbor_distance = this->distance_to(
                *nearest_result.front().second
        ) * 1.2;

        vector<RtreeValue> result_s;
        mRtree_ref.query(
                bgi::intersects(this->createBox(nearest_neighbor_distance)), back_inserter(result_s)
        );

        for (auto &it : result_s) {
            Feature feature = mRtreeFeature_s_ref[it.second].getFeature();
            if (this->distance_to(feature) < nearest_neighbor_distance)
                m_in_range_feature_s.push_back(feature);
        }
    }

    auto get_in_range_feature_s() -> const Feature_s &{ return m_in_range_feature_s; }

    auto updateCovariantMatrix(CentroidRtree & centroidRtreeRef, CentroidListPtr centroidListPtr) -> AxisType {

        Matrix tmpCmat, diffCmat;
        bool is_fresh_start = true;

        auto iterate_parameters = [&]() {

            // Zerolize tmpCmat
            for (int idx_r = 0; idx_r < tmpCmat.rows(); idx_r++)
            for (int idx_c = 0; idx_c < tmpCmat.cols(); idx_c++)
                tmpCmat(idx_r, idx_c) = 0;

            AxisType p_denumerator = 0;

            list<AxisType> weightList;
            for (auto & feature : m_in_range_feature_s) {
                if(is_fresh_start){
                    weightList.push_back(1.0);
                }
                else{
                    weightList.push_back(scoreOfFeature(feature));
                }
                p_denumerator += weightList.back();
            }
            for (auto & feature : m_in_range_feature_s) {

                auto colMeanVec = feature_sub_mean_to_colVector(feature);
                auto weight = weightList.front()/p_denumerator;
                weightList.pop_front();

                for (int idx_r = 0; idx_r < tmpCmat.rows(); idx_r++)
                for (int idx_c = 0; idx_c < tmpCmat.rows(); idx_c++) {
                    if(is_fresh_start){
                        tmpCmat(idx_r, idx_c) += colMeanVec(idx_r, 0) * colMeanVec(idx_c, 0)/
                                                    m_in_range_feature_s.size();
                    }
                    else{
                        tmpCmat(idx_r, idx_c) += (weight/(weight+0.000001)) * 
                                                 colMeanVec(idx_r, 0) * colMeanVec(idx_c, 0)/
                                                    m_in_range_feature_s.size();
                    }
                }
            }
            diffCmat = mCmat - tmpCmat;

            mCmat = tmpCmat;
            mInvCmat = mCmat.inverse();
            mCmatDet = mCmat.determinant();
        };

        AxisType max_diff = 0;

        auto self_consistnt = [&](){
            iterate_parameters();

            if( is_fresh_start ){
                is_fresh_start = false;
                return;
            }

            max_diff = 0;
            for (int idx_r = 0; idx_r < tmpCmat.rows(); idx_r++)
            for (int idx_c = 0; idx_c < tmpCmat.cols(); idx_c++) {
                max_diff = max(max_diff, abs(diffCmat(idx_r, idx_c)) );
            }
        };

        // Perform single shot calculation
        for(int i=0; i<2; i++){
            self_consistnt();
        }
        do{
            self_consistnt();
        } while(max_diff > 0.001);

        return max_diff;
    }



    auto get_model() -> ClusterModel {
        return ClusterModel(this->getFeature(), getCovariantMatrix(), mCentroidKeyStr);
    }

    auto getCovariantMatrix() -> const Matrix & { return mCmat; }

private:
    auto feature_sub_mean_to_rowVector(const Feature &feature) -> RowVector {
        RowVector rowVector;
        for(size_t i=0; i<feature.size(); i++){
            rowVector(0,i) = feature[i] - this->at(i);
        }
        return rowVector;
    }
    auto feature_sub_mean_to_colVector(const Feature &feature) -> ColVector {
        ColVector colVector;
        for(size_t i=0; i<feature.size(); i++){
            colVector(i,0) = feature[i] - this->at(i);
        }
        return colVector;
    }
};
//-----------------------------------
//--- Implementation of Centroid ----
//===================================

}



#endif