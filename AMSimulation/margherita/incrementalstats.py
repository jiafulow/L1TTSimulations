#!/usr/bin/env python

import numpy as np

# ______________________________________________________________________________
# Reference
#   https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
#   http://www.johndcook.com/blog/standard_deviation/
#   http://www.boost.org/doc/libs/1_61_0/boost/accumulators/statistics/mean.hpp

# ______________________________________________________________________________
class IncrementalStats:
    def __init__(self, d=1, cov_d=None, p=0.5, cache_size=1, verbose=0):
        self.d = d
        if cov_d is None:
            cov_d = d
        self.cov_d = cov_d
        self.p = p
        self.cache_size = cache_size
        self.verbose = verbose
        self.cache = []
        self.cnt = 0

        # Sums
        self.sumw = 0.0
        self.sumw2 = 0.0

        # Parametric
        self.v_mean = np.zeros(d)
        self.v_variance = np.zeros(d)
        self.v_mean2 = np.zeros(cov_d)
        self.m_covariance = np.zeros((d,cov_d))

        # Non-parametric
        self.v_minimum = np.ones(d) * float('inf')
        self.v_maximum = np.ones(d) * float('-inf')

    def add(self, variables, covariables=None, weight=None):
        assert(variables.shape == (self.d,))
        if covariables is None:
            covariables = variables
        assert(covariables.shape == (self.cov_d,))

        self.cnt += 1

        # Sums
        if weight is None:
            w = 1.0
        else:
            w = weight
        self.sumw += w
        self.sumw2 += w

        # Parametric
        if weight is None:
            v_delta = (variables - self.v_mean)
            #self.v_mean = (self.v_mean * (self.cnt-1) + variables) / float(self.cnt)
            self.v_mean += v_delta / float(self.cnt)
            #if self.cnt > 1:
                #v_delta = (variables - self.v_mean)
                #self.v_variance += (v_delta * v_delta) / float(self.cnt-1) - self.v_variance / float(self.cnt)
                #self.m_covariance += np.outer(v_delta, v_delta) / float(self.cnt-1) - self.m_covariance / float(self.cnt)
            self.v_variance += v_delta * (variables - self.v_mean)

            v_delta2 = (covariables - self.v_mean2)
            self.v_mean2 += v_delta2 / float(self.cnt)
            self.m_covariance += np.outer(v_delta, v_delta2) * (self.cnt-1) / float(self.cnt)
        else:
            v_delta = (variables - self.v_mean)
            self.v_mean += v_delta * (w/self.sumw)
            self.v_variance += v_delta * (variables - self.v_mean) * w
            v_delta2 = (covariables - self.v_mean2)
            self.v_mean2 += v_delta2 * (w/self.sumw)
            self.m_covariance += np.outer(v_delta, v_delta2) * (w * (self.sumw-w) / self.sumw)

        # Non-parametric
        self.v_minimum = np.minimum(self.v_minimum, variables)
        self.v_maximum = np.maximum(self.v_maximum, variables)
        if len(self.cache) < self.cache_size:
            self.cache.append(variables)
        else:
            pass
        return

    def dim(self):
        return self.d

    def count(self):
        return self.cnt

    def mean(self):
        return self.v_mean

    def variance(self, ddof=0):
        if self.cnt < 2:
            return float('NaN')
        else:
            #return self.v_variance
            return self.v_variance / float(self.sumw - ddof)

    def stdev(self, ddof=0):
        if self.cnt < 2:
            return float('NaN')
        else:
            return np.sqrt(self.variance(ddof=ddof))

    def covariance(self, ddof=0):
        if self.cnt < 2:
            return float('NaN')
        else:
            #return self.m_covariance
            return self.m_covariance / float(self.sumw - ddof)

    def minimum(self):
        return self.v_minimum

    def maximum(self):
        return self.v_maximum

    def quantile(self, p=None):
        if not p:
            p = self.p
        return np.percentile(self.cache, int(p*100), axis=0)


# ______________________________________________________________________________
if __name__ == '__main__':

    d = 6
    nevents = 10000
    myarray = []
    myarray_w = []

    stat = IncrementalStats(d=6, p=0.5, cache_size=nevents/10)

    # Event loop
    for i in xrange(nevents):
        mu, sigma = 0., 1.
        variables = np.random.normal(mu, sigma, d)
        weights = np.random.uniform(0.5, 2.0, 1)
        stat.add(variables, weight=weights[0])

        # Cross check
        #print i, variables
        myarray.append(variables)
        myarray_w.append(weights[0])

    myarray = np.array(myarray)
    myarray_w = np.array(myarray_w)

    print stat.count()
    print stat.mean()
    print stat.variance()
    print stat.stdev()
    print stat.covariance()
    print stat.minimum()
    print stat.maximum()
    print stat.quantile()

    print "Cross check"
    print len(myarray)
    print np.mean(myarray, axis=0)
    print np.var(myarray, axis=0)
    print np.average(myarray, axis=0, weights=myarray_w)
    print np.average(np.square(myarray-np.average(myarray, axis=0, weights=myarray_w)), axis=0, weights=myarray_w)
