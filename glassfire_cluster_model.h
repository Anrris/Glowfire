#ifndef CENTROID_MODEL_H
#define CENTROID_MODEL_H

#include "glassfire_base.h"

namespace glassfire{

template<typename AxisType, size_t Dimension, typename FeatureInfo>
class GlassfireType<AxisType, Dimension, FeatureInfo>::ClusterModel{
    AxisType  _pi_ = 3.141592653589793;
public:
    Feature             mMean;
    Matrix              mCmat;
    Matrix              mInvCmat;
    AxisType            mCmatDet;

    string              mModelKey;

    ClusterModel(const Feature & mean, const Matrix & cmat, const string & modelKey):
        mMean(mean),
        mCmat(cmat),
        mModelKey(modelKey)
    {
        set_covariant_matrix(mCmat);
    }
    auto set_covariant_matrix(const Matrix & cmat) -> void {
        mCmat = cmat;
        mInvCmat = mCmat.inverse();
        mCmatDet = mCmat.determinant();
    }
    auto eval(Feature feature) -> AxisType {
        auto rowVec = feature_sub_mean_to_rowVector(feature);
        auto colVec = feature_sub_mean_to_colVector(feature);
        AxisType mahalanDistance = (rowVec * mInvCmat * colVec)(0, 0);
        AxisType result = (AxisType)exp( -0.5* mahalanDistance) / sqrt( 2 * _pi_ * mCmatDet );
        return result;
    }
    auto mean() -> const Feature & {return mMean;}
    auto covariant_matrix() -> const Matrix & {return mCmat;}

private:
    auto feature_sub_mean_to_rowVector(const Feature &feature) -> RowVector {
        RowVector rowVector;
        for(size_t i=0; i<feature.size(); i++){
            rowVector(0,i) = feature[i] - mMean[i];
        }
        return rowVector;
    }
    auto feature_sub_mean_to_colVector(const Feature &feature) -> ColVector {
        ColVector colVector;
        for(size_t i=0; i<feature.size(); i++){
            colVector(i,0) = feature[i] - mMean[i];
        }
        return colVector;
    }
};

}


#endif