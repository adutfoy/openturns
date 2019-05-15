#! /usr/bin/env python

from __future__ import print_function
from openturns import *

TESTPREAMBLE()
RandomGenerator.SetSeed(0)

try:
    # Instanciate one distribution object
    R = CorrelationMatrix(3)
    R[0, 1] = 0.5
    R[0, 2] = 0.25
    d1 =  ot.ComposedDistribution([ot.Uniform(),  ot.Triangular(1.0, 2.0, 3.0)], ot.NormalCopula(R))
    d2 = ot.Normal(2)
    collection = [d1, d2]
    dist = ot.BlockIndependentDistribution(collection)

    print("Distribution ", repr(dist))
    print("Distribution ", Distribution)

    # Is this copula elliptical ?
    print("Elliptical distribution= ", dist.isElliptical())

    # Is this copula continuous ?
    print("Continuous = ", dist.isContinuous())

    # Is this copula elliptical ?
    print("Elliptical = ", dist.hasEllipticalCopula())

    # Is this copula independent ?
    print("Independent = ", dist.hasIndependentCopula())

    # Test for realization of copula
    oneRealization = dist.getRealization()
    print("oneRealization=", repr(oneRealization))

    # Test for sampling
    size = 10000
    oneSample = dist.getSample(size)
    print("oneSample first=", repr(
        oneSample[0]), " last=", repr(oneSample[size - 1]))
    print("mean=", repr(oneSample.computeMean()))
    print("covariance=", repr(oneSample.computeCovariance()))

    # Define a point
    point = Point(dist.getDimension(), 0.6)
    print("Point= ", repr(point))

    # Show PDF and CDF of point
    cdfPoint = dist.computeCDF(point)
    pdfPoint = dist.computePDF(point)
    print('CDF(point)=', repr(cdfPoint))
    print('PDF(point)=', repr(pdfPoint))
    cdfPointRef = d1.computeCDF(point.getMarginal([0,1,2])) * d2.computeCDF(point.getMarginal([3,4]))
    print('CDF ref (point)=', repr(cdfPointRef))
    pdfPointRef = d1.computePDF(point.getMarginal([0,1,2])) * d2.computePDF(point.getMarginal([3,4]))
    print('PDF ref (point)=', repr(pdfPointRef))

    # Quantile
    quantile = dist.computeQuantile(0.95)
    print("quantile=", repr(quantile))
    print("cdf(quantile)=%.6f" % dist.computeCDF(quantile))
    
    # Get 95% survival function
    inverseSurvival = Point(dist.computeInverseSurvivalFunction(0.95))
    print("InverseSurvival=", repr(inverseSurvival))
    print("Survival(inverseSurvival)=%.6f" %
          dist.computeSurvivalFunction(inverseSurvival))
    print("entropy=%.6f" % dist.computeEntropy())

    # Mean, parameters
    mean = dist.getMean()
    mean1 = d1.getMean()
    mean2 = d2.getMean()
    print("mean=", repr(mean))
    print("mean ref=", repr(mean1.stack(mean2)))
    
    parameters = dist.getParametersCollection()
    print("parameters=", repr(parameters))
    correlation = dist.getCorrelation()
    print("correlation=", correlation)
    spearman = dist.getSpearmanCorrelation()
    print("spearman=", spearman)
    kendall = dist.getKendallTau()
    print("kendall=", kendall)
    
    # Extract a 4-D marginal
    dim = 4
    point = Point(dim, 0.25)
    indices = Indices(dim, 0)
    indices[0] = 1
    indices[1] = 2
    indices[2] = 4
    indices[3] = 5
    print("indices=", repr(indices))
    margins = dist.getMarginal(indices)
    print("margins=", repr(margins))
    print("margins PDF=%.6f" % margins.computePDF(point))
    print("margins CDF=%.6f" % margins.computeCDF(point))
    quantile = margins.computeQuantile(0.95)
    print("margins quantile=", repr(quantile))
    print("margins CDF(quantile)=%.6f" % margins.computeCDF(quantile))
    
    print("margins realization=", repr(margins.getRealization()))

except:
    import sys
    print("t_ComposedCopula.py", sys.exc_info()[0], sys.exc_info()[1])
