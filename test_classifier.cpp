//
// Created by Yuan Yen Tai on 11/23/17.
// All rights reserved.
//

#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
using namespace std;

#include "glassfire.h"

int main(){

    typedef glassfire::Classifier<double,2, std::string> Classifier;
    auto feature_s  = Classifier::Feature_s();
    auto classifier = Classifier();

    size_t Index;
    double dummy, axis_0, axis_1;

    cout << "-- Load feature from file ..." << endl;
    auto infile = ifstream("rand.csv");
    while( infile >> Index >> dummy >> axis_0 >> axis_1 )
        feature_s.push_back({axis_0, axis_1});
    infile.close();

    cout << "Append features ..." << endl;
    for(size_t i=0; i<feature_s.size(); i++){
        classifier.append_feature(feature_s[i], "test"+std::to_string(i));
    }

    cout << "Run cluster algorithm ..." << endl;
    auto scorerSet = classifier.run_cluster(8.0, 1, 0.0001);
    scorerSet->get_model_set(1);

    vector<string> cluster_id_s;
    Classifier::Feature_s result_s;

    //cout << "Predict scores from the scorer ..." << endl;
    //for(auto & feature: feature_s){
    //    auto score_dict = scorerSet->calc_scores(feature);

    //    cluster_id_s.push_back(score_dict.begin()->second);

    //    Classifier::Feature predict = {score_dict.begin()->first};
    //    for(auto & elem: feature){
    //        predict.push_back(elem);
    //    }

    //    result_s.push_back(predict);
    //}

    //cout << "Save result to a file ..." << endl;
    //auto outfile = ofstream("api-demo_predict.csv");
    //for(size_t i=0; i<cluster_id_s.size(); i++){
    //    outfile << cluster_id_s[i] << " ";
    //    auto & predict = result_s[i];
    //    for(auto & elem: predict) outfile << elem << " ";
    //    outfile << endl;
    //}
    //outfile.close();

    //cout << "Print Centroid id/positions ..." << endl;
    //outfile = ofstream("api-demo_centroid.csv");
    //auto centroids = scorerSet->get_model_set();
    //for(auto iter: centroids){
    //    auto feature = iter.mean();
    //    cout << iter.model_key() << " " << iter.cmean() << endl;
    //    cout << iter.cov_mat() << endl;
    //    outfile << iter.model_key() << " " << iter.cmean() << endl;
    //}
    //outfile.close();

    //auto result = scorerSet->query({23.8, 52.7}, 20);
    //cout << std::get<0>(result)<<endl;
    //cout << std::get<1>(result)<<endl;
    //cout << std::get<2>(result).model_key()<<endl;

    //auto tmp = classifier.query_data(std::get<2>(result), 0.5);

    //for(auto &iter: tmp){
    //    std::cout << std::get<0>(iter) << std::endl;
    //    std::cout << std::get<1>(iter) << std::endl;
    //    std::cout << std::get<2>(iter) << std::endl;
    //    std::cout << std::get<4>(iter)[0] << std::endl;
    //    std::cout << "-------" << std::endl;
    //}

    return 0;
}
