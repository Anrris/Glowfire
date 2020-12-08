#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <iostream>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;

#include "src/glassfire.h"

class Cluster {
public:
    typedef double AxisType;
    typedef std::string FeatureInfo;

    typedef std::shared_ptr<glassfire::ClassifierBase<AxisType, FeatureInfo>> ClassifierBasePtr;
    typedef std::shared_ptr<glassfire::ScorerSetBase<AxisType, FeatureInfo>> ScorerSetBasePtr;
    typedef glassfire::ClusterModel<AxisType, FeatureInfo> ClusterModel;
    typedef glassfire::ClassifierBase<AxisType, FeatureInfo>::DataQueryReturnType DataQueryReturnType;
private:
    std::vector<ClassifierBasePtr> mClassifierBaseSet =
        glassfire::ClassifierFactory<double, FeatureInfo, 10>().create();

    ScorerSetBasePtr mScorerSetBase;
    ssize_t m_dimension;
    auto _cluster() -> ClassifierBasePtr {
        return mClassifierBaseSet[m_dimension-1];
    }

    bool m_been_executed = false;

public:

    Cluster(Eigen::MatrixXd input, std::vector<FeatureInfo> info, std::string loading_type){
        if(loading_type != "row" && loading_type != "col"){
            throw std::runtime_error(":loading_type: (='row' or ='col) must be given.");
        }
        if(input.size() == 0)
            throw std::runtime_error(":input: is an empty array.");
        
        if (info.size() !=0 ){
            if(loading_type == "row" && (size_t)input.rows() != info.size()){
                throw std::runtime_error(":info: dimension does not match :input:.");
            }
            if(loading_type == "col" && (size_t)input.cols() != info.size()){
                throw std::runtime_error(":info: dimension does not match :input:.");
            }
        }

        ssize_t OUTER_LOOP=0, INNER_LOOP=0;
        if(loading_type == "row"){
            OUTER_LOOP = input.rows();
            INNER_LOOP = input.cols();
            m_dimension = input.cols();
        }
        else if(loading_type == "col"){
            OUTER_LOOP = input.cols();
            INNER_LOOP = input.rows();
            m_dimension = input.rows();
        }

        if(info.size() == 0){
            for(int i=0; i<OUTER_LOOP; i++) info.push_back("NONE.");
        }

        if ((size_t)INNER_LOOP > mClassifierBaseSet.size())
            throw std::runtime_error("Cannot handle feature size over: "+std::to_string(mClassifierBaseSet.size()));

        for(ssize_t i=0; i<OUTER_LOOP; i++){
            vector<double> feature(INNER_LOOP);
            for(ssize_t j=0; j<INNER_LOOP; j++){
                if(loading_type == "row") feature[j] = input(i,j);
                if(loading_type == "col") feature[j] = input(j,i);
            }
            _cluster()->append_feature(feature, info[i]);
        }
    }

    auto run(AxisType box_len, bool using_box_len) -> void {
        m_been_executed = true;
        mScorerSetBase = _cluster()->run_cluster(box_len, 1, 0.0001, using_box_len);
    }
    auto get_models(AxisType regularize) -> std::vector<ClusterModel> {
        if(!m_been_executed)
            throw std::runtime_error("Cluster.run(...) never been executed before using 'get_models'.");
        return mScorerSetBase->get_model_set(regularize);
    }
    auto query_model(std::vector<AxisType> & input, AxisType regularize)
            -> std::tuple<bool, AxisType, glassfire::ClusterModel<AxisType, FeatureInfo>, std::string> {
        if(!m_been_executed)
            throw std::runtime_error("Cluster.run(...) never been executed before using 'query_model'.");
        return mScorerSetBase->query(input, regularize, -1);
    }
    auto query_data(ClusterModel & model, AxisType box_len) -> DataQueryReturnType {
        if(!m_been_executed)
            throw std::runtime_error("Cluster.run(...) never been executed before using 'query_data'.");
        return _cluster()->query_data(model, box_len);
    }
};

PYBIND11_MODULE(glassfire, m) {
    py::class_<Cluster::ClusterModel>(m, "ClusterModel")
        .def("eval", &Cluster::ClusterModel::eval)
        .def("count", &Cluster::ClusterModel::count)
        .def("cmean", &Cluster::ClusterModel::cmean, py::return_value_policy::reference_internal)
        .def("cov_mat", &Cluster::ClusterModel::cov_mat, py::return_value_policy::reference_internal)
        .def("model_key", &Cluster::ClusterModel::model_key, py::return_value_policy::reference_internal)
        .def("get_data_index", &Cluster::ClusterModel::get_data_index, py::return_value_policy::reference_internal)
        .def("set_regularize", &Cluster::ClusterModel::set_regularize)
        ;

    py::class_<Cluster>(m, "Cluster")
        .def(py::init<Eigen::MatrixXd, std::vector<std::string>, std::string>(),
          py::arg("input"),
          py::arg("info")=std::vector<Cluster::FeatureInfo>(),
          py::arg("loading_type")="row"
        )
        .def("run", &Cluster::run,
          py::arg("box_len"),
          py::arg("using_box_len")
        )
        .def("query_model", &Cluster::query_model,
          py::arg("input"),
          py::arg("regularize")
        )
        .def("query_data", &Cluster::query_data,
          py::arg("model"),
          py::arg("box_len")
        )
        .def("get_models", &Cluster::get_models,
          py::arg("regularize")
        )
        ;
}
