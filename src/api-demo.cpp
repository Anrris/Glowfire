//
// Created by Yuan Yen Tai on 11/23/17.
// All rights reserved.
//

#include <iostream>
#include <fstream>
#include <chrono>
using namespace std;

#include "glassfire.h"

int main(){

    typedef glassfire::Classifier<float,2, size_t> Classifier;
    auto feature_s  = Classifier::Feature_s();
    auto classifier = Classifier();

    size_t Index;
    float dummy, axis_0, axis_1;

    cout << "-- Load feature from file ..." << endl;
    auto infile = ifstream("rand.csv");
    while( infile >> Index >> dummy >> axis_0 >> axis_1 )
        feature_s.push_back({axis_0, axis_1});
    infile.close();

    cout << "Append features ..." << endl;
    for(size_t i=0; i<feature_s.size(); i++){
        classifier.append_feature(feature_s[i], i);
    }

    cout << "Run lgm cluster algorithm ..." << endl;
    auto scorer = classifier.run_cluster(8.0);

    vector<string> cluster_id_s;
    Classifier::Feature_s result_s;

    cout << "Predict scores from the scorer ..." << endl;
    for(auto & feature: feature_s){
        auto score_dict = scorer.calc_scores(feature);

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

    return 0;
}
