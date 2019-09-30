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

    typedef glassfire::Classifier<double,2, size_t> Classifier;
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
        classifier.append_feature(feature_s[i], i);
    }

    cout << "Run cluster algorithm ..." << endl;
    auto scorerSet = classifier.run_cluster(8.0, 1, 0.0001);

    vector<string> cluster_id_s;
    Classifier::Feature_s result_s;

    cout << "Predict scores from the scorer ..." << endl;
    for(auto & feature: feature_s){
        auto score_dict = scorerSet.calc_scores(feature);

        cluster_id_s.push_back(score_dict.begin()->second);

        Classifier::Feature predict = {score_dict.begin()->first};
        for(auto & elem: feature){
            predict.push_back(elem);
        }

        result_s.push_back(predict);
    }

    cout << "Save result to a file ..." << endl;
    auto outfile = ofstream("api-demo_predict.csv");
    for(size_t i=0; i<cluster_id_s.size(); i++){
        outfile << cluster_id_s[i] << " ";
        auto & predict = result_s[i];
        for(auto & elem: predict) outfile << elem << " ";
        outfile << endl;
    }
    outfile.close();

    cout << "Print Centroid id/positions ..." << endl;
    outfile = ofstream("api-demo_centroid.csv");
    auto centroids = scorerSet.get_centroids();
    for(auto iter: centroids){
        auto feature = iter.getFeature();
        cout << iter.getKeyStr() << " " << feature[0] << " " << feature[1] << endl;
        outfile << iter.getKeyStr() << " " << feature[0] << " " << feature[1] << endl;
    }
    outfile.close();

    auto result = scorerSet.query_centroids({23.8, 52.7}, 20);
    auto model = result.second;
    model = result.second;
    cout << model.mModelKey << endl;
    cout << model.eval({24.0265, 52.9805}) << endl;

    return 0;
}
