//
// Created by Tai, Yuan-yen on 11/11/17.
// All rights reserved.
//

#include <iostream>
#include <fstream>
#include <chrono>
using namespace std;

#include "lgm.h"

template<typename TimeT = std::chrono::milliseconds>
struct measure
{
    template<typename F, typename ...Args>
    static typename TimeT::rep execution(F&& func, Args&&... args)
    {
        auto start = std::chrono::steady_clock::now();
        std::forward<decltype(func)>(func)(std::forward<Args>(args)...);
        auto duration = std::chrono::duration_cast< TimeT>
                (std::chrono::steady_clock::now() - start);
        return duration.count();
    }
};

int main(){
    /*
     * The perf-test.cpp is testing performances of the LGM algorithm.
     * The execution time will be evaluated inside measure<>::execution
     * with a lambda input.
     * */



    typedef LGM::LgmClassifier<float,2> Classifier;
    auto feature_s  = Classifier::Feature_s();
    auto classifier = Classifier();

    Classifier::VarType cluster_distance;
    cout << "Please enter the cluster distance:";
    cin  >> cluster_distance;

    size_t Index;
    Classifier::VarType dummy, axis_0, axis_1;


    //-------------------------------------------------
    cout << "Feature loading : ";
    cout << measure<>::execution([&](){

        auto infile = ifstream("rand.csv");
        while( infile >> Index >> dummy >> axis_0 >> axis_1 )
            feature_s.push_back({axis_0, axis_1});
        infile.close();

    });
    cout << " msec"<< endl;

    cout << "Total: " << feature_s.size() << " of features." << endl;

    //-------------------------------------------------
    cout << "Append feature to rtree : ";
    cout << measure<>::execution([&](){
        // Append features
        for(size_t i=0; i<feature_s.size(); i++){
            classifier.append_feature(feature_s[i]);
        }

    });
    cout << " msec"<< endl;

    //-------------------------------------------------
    Classifier::Scorer scorer;
    cout << "Execute clustering algorithm : ";
    cout << measure<>::execution([&](){
        // run cluster algorithm
        cout << "Found: "<< classifier.run_cluster(cluster_distance) << " clusters."<<endl;
        scorer = classifier.create_scorer();
    });
    cout << " msec"<< endl;

    //-------------------------------------------------
    vector<size_t> cluster_id_s;
    Classifier::Feature_s result_s;
    cout << "Predict score and classify : ";
    cout << measure<>::execution([&](){

        for(auto & feature: feature_s){
            auto score_dict = scorer.calc_score(feature);

            cluster_id_s.push_back(score_dict.begin()->second);

            Classifier::Feature predict = {score_dict.begin()->first};
            for(auto & elem: feature){
                predict.push_back(elem);
            }
            result_s.push_back(predict);
        }
    });
    cout << " msec"<< endl;

    //-------------------------------------------------
    cout << "Save predictions to file : ";
    cout << measure<>::execution([&](){
        auto outfile = ofstream("predict.csv");

        for(size_t i=0; i<cluster_id_s.size(); i++){
            outfile << cluster_id_s[i] << " ";
            auto & predict = result_s[i];
            for(auto & elem: predict) outfile << elem << " ";
            outfile << endl;
        }

        for(auto & centroid: scorer.get_centroids()){
            outfile << "C ";
            for(auto & coord: centroid.getFeature()){
                outfile << coord << " ";
            }
            outfile << centroid.count() << endl;
        }

        outfile.close();
    });
    cout << " msec"<< endl;

    return 0;
}