# README #

This repository is the implementation of the following arXiv article:
https://arxiv.org/abs/1709.08470

Where I introduced a novel clustering algorithm, named Local Gaussian Model (LGM), that depends only on a single parameter to describe the approximate separation distance between clusters.

Academic research of using this algorithm including the source code must reference this arXiv article.


### How to use LGM ###

The followings are the major API component:

Create a classifier instance:
```cpp
auto classifier = LGM::LgmClassifier<double,2>();
```

Append feature to the cluster:
```cpp
// Usually, one need to append many feature as possible.
classifier.append_feature({num1, num2});
```

Append feature to the cluster:
```cpp
classifier.run_cluster(7.0);
```

Return the score of a feature from the classifier:
```cpp
auto score_dict = classifier.calc_score(feature);
```

A typicle example can be found in ```lgmbin/api-demo.cpp```

### How to generate a set of testing data ? ###
Execute the following python script which located inside the py folder:
```bash
$ ./cluster_generator.py
```
Will generate ```rand.csv``` and ```rand.mc``` files.
Where "rand.csv" is a set of several gaussian clusters and "rand.mc" is the covariant matrix of them.

One can also use the ```plot_scatter.py``` to plot the classified data set with a given cluster-id.
```bash
$ ./plot_scatter.py  api-demo_predict.csv
```
where the api-demo_predict.csv can be generated from the exectuable of ```api-demo.cpp```.

### How do I get set up? ###

* This library depends on the "boost" and "Eigen" C++ library. Please download these library yourself and setup correspond library path for them.

* The source code is entirely written in C++ header files, one just need to include the ```"lgm.h"``` file from the ```include``` directory of this repository.

