//
// Created by Tai, Yuan-yen on 11/11/17.
//

#ifndef LOCALGAUSSIANMODEL_LGM_FACTORY_H
#define LOCALGAUSSIANMODEL_LGM_FACTORY_H

#include "lgm_classifier.h"

namespace LGM{

    // A factory class which can generate a set of N dimensional classifier from N=x... to N=y.
    //
    // Attension, Not fully tested, use at your own risk!!!!
    //

    template<int From, int To>
    class LgmClassifierFactory{
        template<int x, int to>
        struct static_for
        {
            void operator()(vector<LgmBase<double>::Ptr> & refProduct) {
                refProduct.push_back(
                        make_shared<LgmClassifier<double, x>>()
                );
                static_for<x+1,to>()(refProduct);
            }
        };
        template<int to>
        struct static_for<to,to>
        {
            void operator()(vector<LgmBase<double>::Ptr> &refProduct) {}
        };

        vector<LgmBase<double>::Ptr> factoryProduct;
    public:
        LgmClassifierFactory(){
            static_for<From, To>()(factoryProduct);
        }

        LgmBase<double>::Ptr getProduct(int product_index){
            return factoryProduct[product_index];
        }
    };

}
#endif //LOCALGAUSSIANMODEL_LGM_FACTORY_H
