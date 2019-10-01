#ifndef CENTROID_MODEL_H
#define CENTROID_MODEL_H

#include "glassfire_base.h"

namespace glassfire{

template<typename AxisType, typename FeatureInfo>
class ClusterModel{
    AxisType  _pi_ = 3.141592653589793;
public:
    typedef Eigen::Matrix<AxisType, 1, Eigen::Dynamic>               RowVector;
    typedef Eigen::Matrix<AxisType, Eigen::Dynamic, 1>               ColVector;
    typedef std::vector<AxisType>    Feature;

    Feature             mMean;
    Eigen::MatrixXd     mCmat;
    Eigen::MatrixXd     mInvCmat;
    AxisType            mCmatDet;

    std::string              mModelKey;

    size_t              Dimension;

    ClusterModel(const Feature & mean, const Eigen::MatrixXd & cmat, const std::string & modelKey):
        mMean(mean),
        mCmat(cmat),
        mModelKey(modelKey)
    {
        Dimension = mean.size();
        set_mean_cmat(mean, mCmat);
    }
    auto set_mean_cmat(const Feature & mean, const Eigen::MatrixXd & cmat) -> void {
        mMean = mean;
        mCmat = cmat;
        mInvCmat = mCmat.inverse();
        mCmatDet = mCmat.determinant();
    }
    auto eval(Feature feature) -> AxisType {
        auto rowVec = feature_sub_mean_to_rowVector(feature);
        auto colVec = feature_sub_mean_to_colVector(feature);
        AxisType mahalanDistance = (rowVec * mInvCmat * colVec)(0, 0);
        AxisType result = (AxisType)exp( -0.5* mahalanDistance) / sqrt( pow(2 * _pi_, feature.size()) * mCmatDet );
        return result;
    }
    auto mean() -> const Feature & {return mMean;}
    auto covariant_matrix() -> const Eigen::MatrixXd & {return mCmat;}

private:
    auto feature_sub_mean_to_rowVector(const Feature &feature) -> Eigen::MatrixXd {
        Eigen::MatrixXd rowVector(1, Dimension);
        for(size_t i=0; i<feature.size(); i++){
            rowVector(0,i) = feature[i] - mMean[i];
        }
        return rowVector;
    }
    auto feature_sub_mean_to_colVector(const Feature &feature) -> Eigen::MatrixXd {
        Eigen::MatrixXd colVector(Dimension, 1);
        for(size_t i=0; i<feature.size(); i++){
            colVector(i,0) = feature[i] - mMean[i];
        }
        return colVector;
    }
};

}


#endif