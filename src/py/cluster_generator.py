#!/usr/bin/env python
#
# Created by Tai, Yuan-yen on 11/17/17.
# All rights reserved.
#
import numpy as np


class NDRandClusterGenerator(object):
    def __init__(self, dim=1):
        if dim < 1:
            raise Exception("Dimension should be grater than 0.")

        self.dimension = dim
        self.point_collections = []
        self.mean_cov_collections = []
        self.cluster_count = 0

    def random_symm_matrix(self, diag, offdiag):
        asymm = np.random.rand(self.dimension, self.dimension)
        asymm = np.dot(asymm, asymm.T)
        symm = (asymm + asymm.T) / 2
        (col, row) = symm.shape

        (diag_base, diag_scal) = diag
        (offdiag_base, offdiag_scal) = offdiag
        for i in xrange(col):
            for j in xrange(i + 1, row):
                flag = np.random.random_integers(0, 1, 1)
                m = 1
                if flag == 0: m = -1
                symm[i, j] = m * (offdiag_base + symm[i, j] * offdiag_scal)
                symm[j, i] = m * (offdiag_base + symm[j, i] * offdiag_scal)
            symm[i, i] = diag_base + symm[i, i] * diag_scal
        return symm

    def seed(self, mean, cov, count):
        self.mean_cov_collections.append((mean, cov))
        for elem in np.random.multivariate_normal(mean, cov, count).tolist():
            self.point_collections.append((self.cluster_count, elem))
        self.cluster_count += 1

    def seed_in_range(self, mean, count, diag=0, offdiag=1):
        self.seed(mean, self.random_symm_matrix(diag, offdiag), count)

    def save_data(self, filename):
        # Save the distributed data point to file
        out = open(filename+".csv", 'w')
        for row in self.point_collections:
            (cluster, data) = row
            out.write(str(cluster)+" ")
            for elem in data:
                out.write(str(elem) + " ")
            out.write('\n')
        out.close()

        # Save the Gaussian distribution to file
        out = open(filename+".mc", 'w')
        out.write("Dimension = "+str(self.dimension))
        for item in self.mean_cov_collections:
            out.write(">>>")
            (mean_mat, cov_mat) = item
            for mean_elem in mean_mat:
                out.write(str(mean_elem) + " ")
            out.write('\n')

            (row, col) = cov_mat.shape
            for i in xrange(row):
                for j in xrange(col):
                    out.write(str(cov_mat[i,j]) + " ")
                out.write('\n')
            out.write('\n')
        out.close()



if __name__ == "__main__":
    test = NDRandClusterGenerator(2)
    test.seed_in_range(mean=(0, 0),  count=30000, diag=(3, 3), offdiag=(1, 1))
    test.seed_in_range(mean=(10, 0), count=30000, diag=(3, 3), offdiag=(1, 2))
    test.seed_in_range(mean=(0, 10), count=30000, diag=(3, 3), offdiag=(2, 1))
    test.seed_in_range(mean=(10, 10),count=30000, diag=(3, 3), offdiag=(1, 2))
    test.seed_in_range(mean=(20, 15),count=30000, diag=(3, 3), offdiag=(1, 2))
    test.seed_in_range(mean=(20, 30),count=30000, diag=(3, 3), offdiag=(1, 2))
    test.seed_in_range(mean=(20, 40),count=30000, diag=(6, 6), offdiag=(3, 2))
    test.seed_in_range(mean=(10, 40),count=30000, diag=(6, 6), offdiag=(3, 2))
    test.seed_in_range(mean=(44, 40),count=30000, diag=(4, 7), offdiag=(3, 2))
    test.seed_in_range(mean=(44, 32),count=30000, diag=(4, 7), offdiag=(3, 2))

    test.save_data("rand")
