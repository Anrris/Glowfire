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

    typedef LGM::LgmClassifier<double,2>::Scorer Scorer;
    auto feature_s  = LGM::LgmClassifier<double,2>::Feature_s();
    auto classifier = LGM::LgmClassifier<double,2>();

    size_t Index;
    double dummy, axis_0, axis_1;


    cout << "Feature loading : ";
    cout << measure<>::execution([&](){

        auto infile = ifstream("rand.csv");
        while( infile >> Index >> dummy >> axis_0 >> axis_1 )
            feature_s.push_back({axis_0, axis_1});
        infile.close();

    });
    cout << " msec"<< endl;

    cout << "Total: " << feature_s.size() << " of features." << endl;

    cout << "Append feature to rtree : ";
    cout << measure<>::execution([&](){
        // Append features
        //for(auto & feature: feature_s){
        for(size_t i=0; i<feature_s.size(); i++){
            if(i%10==0) classifier.append_feature(feature_s[i]);
        }

    });
    cout << " msec"<< endl;


    Scorer scorer1, scorer2;
    cout << "Execute clustering algorithm : ";
    cout << measure<>::execution([&](){
        // run cluster algorithm
        classifier.run_cluster(10.0);
        scorer1 = classifier.create_scorer();
        classifier.run_cluster(7.0);
        scorer2 = classifier.create_scorer();
    });
    cout << " msec"<< endl;

    auto saveToFile = [&](LGM::LgmClassifier<double,2>::Scorer & scorer, string filename){
        vector<size_t> cluster_id_s;
        vector<vector<double>> result_s;
        cout << "Predict score and classify : ";
        cout << measure<>::execution([&](){

            for(auto & feature: feature_s){
                auto score_dict = scorer.calc_score(feature);

                cluster_id_s.push_back(score_dict.begin()->second);

                vector<double> predict = {score_dict.begin()->first};
                for(auto & elem: feature){
                    predict.push_back(elem);
                }

                result_s.push_back(predict);
            }

        });
        cout << " msec"<< endl;

        cout << "Save to file : ";
        cout << measure<>::execution([&](){
            auto outfile = ofstream(filename);

            for(size_t i=0; i<cluster_id_s.size(); i++){
                outfile << cluster_id_s[i] << " ";
            auto & predict = result_s[i];
            for(auto & elem: predict) outfile << elem << " ";
            outfile << endl;
        }

        outfile.close();
        });
        cout << " msec"<< endl;
    };

    saveToFile(scorer1, "predict1.csv");
    saveToFile(scorer2, "predict2.csv");

    return 0;
}