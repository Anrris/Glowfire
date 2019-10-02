#ifndef CENTROID_MODEL_H
#define CENTROID_MODEL_H

#include "glassfire_base.h"

namespace glassfire{

template<typename AxisType, typename FeatureInfo>
class ClusterModel{
public:
    typedef Eigen::Matrix<AxisType, Eigen::Dynamic, Eigen::Dynamic>               Matrix;
    typedef Eigen::Matrix<AxisType, 1, Eigen::Dynamic>               RowVector;
    typedef Eigen::Matrix<AxisType, Eigen::Dynamic, 1>               ColVector;
    typedef std::vector<AxisType>    Feature;
private:
    AxisType  _pi_ = 3.141592653589793;
    RowVector   mRowMean;
    Feature     mMean;

    Matrix      mCmat;
    Matrix      mInvCmat;
    AxisType    mCmatDet;
    std::string mModelKey;
    size_t      Dimension;

    auto set_mean_cmat(const Feature & mean, const Matrix & cmat) -> void {
        Dimension = mean.size();
        mRowMean = RowVector(1, Dimension);
        for(size_t i=0; i<mean.size(); i++)
            mRowMean(0, i) = mean[i];
        mMean = mean;
        mCmat = cmat;
        mInvCmat = mCmat.inverse();
        mCmatDet = mCmat.determinant();
    }
public:
    ClusterModel(){}

    ClusterModel(const Feature & mean, const Matrix & cmat, const std::string & modelKey):
        mMean(mean),
        mCmat(cmat),
        mModelKey(modelKey)
    {
        set_mean_cmat(mean, mCmat);
    }

    auto eval(Feature feature) -> AxisType {
        auto rowVec = feature_sub_mean_to_rowVector(feature);
        auto colVec = feature_sub_mean_to_colVector(feature);
        AxisType mahalanDistance = (rowVec * mInvCmat * colVec)(0, 0);
        AxisType result = (AxisType)exp( -0.5* mahalanDistance) / sqrt( pow(2 * _pi_, feature.size()) * mCmatDet );
        return result;
    }
    auto cmean()-> const RowVector & {return mRowMean;}
    auto mean() -> const Feature & {return mMean;}
    auto cov_mat() -> const Matrix & {return mCmat;}
    auto model_key() -> const string & {return mModelKey;}

private:
    auto feature_sub_mean_to_rowVector(const Feature &feature) -> Matrix {
        Matrix rowVector(1, Dimension);
        for(size_t i=0; i<feature.size(); i++){
            rowVector(0,i) = feature[i] - mMean[i];
        }
        return rowVector;
    }
    auto feature_sub_mean_to_colVector(const Feature &feature) -> Matrix {
        Matrix colVector(Dimension, 1);
        for(size_t i=0; i<feature.size(); i++){
            colVector(i,0) = feature[i] - mMean[i];
        }
        return colVector;
    }
};

}


#endif