//
// Created by Tai, Yuan-yen on 11/11/17.
// All rights reserved.
//

#include <iostream>
#include <fstream>
using namespace std;

#include "lgm.h"

int main(){

    auto feature_s  = LGM::LgmClassifier<double,2>::Feature_s();

    size_t Index;
    double axis_0, axis_1;

    // Load features from a given file.
    auto infile = ifstream("rand.csv");
    while( infile >> Index >> axis_0 >> axis_1 )
        feature_s.push_back({axis_0, axis_1});
    infile.close();

    // Create instance of LgmClassifier
    auto classifier = LGM::LgmClassifier<double,2>();

    // Append features
    for(auto & feature: feature_s)
        classifier.append_feature(feature);

    // Run the cluster algorithm with separation = 7.0.
    classifier.run_cluster(7.0);

    // Save classified cluster (with cluster-id) to file.
    auto outfile = ofstream("predict.csv");
    for(auto & feature: feature_s){
        auto score_dict = classifier.calc_score(feature);
        outfile << score_dict.begin()->second<< "\t\t";
        outfile << score_dict.begin()->first<< " ";
        for(auto & elem: feature){
            outfile << elem << "  ";
        }
        outfile << endl;
    }
    outfile.close();

    return 0;
}