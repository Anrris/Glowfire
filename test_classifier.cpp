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

    cout << "Print Centroid id/positions ..." << endl;
    auto centroids = scorer.get_centroids();
    for(auto iter: centroids){
        auto feature = iter.getFeature();
        cout << iter.getKeyStr() << " " << feature[0] << " " << feature[1] << endl;
    }

    auto result = scorer.query_centroids({23.8, 52.7}, 20);
    cout << endl;
    cout << result.second.mModelKey << endl;
    cout << endl;

    auto model = result.second;
    model = result.second;
    cout << model.mModelKey << endl;
    cout << model.eval({23, 52.7}) << endl;
    cout << model.eval({23, 55.7}) << endl;
    cout << model.mCmatDet << endl;

    result = scorer.query_centroids({20, 50}, 20);
    cout << endl;
    cout << result.second.mModelKey << endl;
    cout << result.first;
    cout << endl;



    return 0;
}
