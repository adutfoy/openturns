//                                               -*- C++ -*-
/**
 *  @brief The test file of class Triangular for standard methods
 *
 *  Copyright 2005-2025 Airbus-EDF-IMACS-ONERA-Phimeca
 *
 *  This library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include "openturns/OT.hxx"
#include "openturns/OTtestcode.hxx"

using namespace OT;
using namespace OT::Test;

class TestObject : public Triangular
{
public:
  TestObject() : Triangular(-0.5, 1.5, 2.5) {}
  virtual ~TestObject() {}
};


int main(int, char *[])
{
  TESTPREAMBLE;
  OStream fullprint(std::cout);

  try
  {
    // Test basic functionnalities
    checkClassWithClassName<TestObject>();

    Collection< Triangular > coll;
    coll.add(Triangular(-0.5,  1.5, 2.5));
    coll.add(Triangular(-0.5, -0.5, 2.5));
    coll.add(Triangular(-0.5,  2.5, 2.5));
    coll.add(Triangular(-2.5,  0.0, 2.5));
    Point u = {0.1, 0.01, 0.001, 0.0001, 0.00001};
    Collection< Collection<Complex> > refValues(4);
    refValues[0].add(Complex(9.9127099972903484510e-01, 1.1618827274767113648e-01));
    refValues[0].add(Complex(9.9991250210240596186e-01, 1.1666187507734305978e-02));
    refValues[0].add(Complex(9.9999912500021024303e-01, 1.1666661875000773437e-03));
    refValues[0].add(Complex(9.9999999125000002102e-01, 1.1666666618750000077e-04));
    refValues[0].add(Complex(9.9999999991250000000e-01, 1.1666666666187500000e-05));
    refValues[1].add(Complex(9.9625754367833246793e-01, 4.9821102073080346159e-02));
    refValues[1].add(Complex(9.9996250075519992255e-01, 4.9998208360230421163e-03));
    refValues[1].add(Complex(9.9999962500007552082e-01, 4.9999982083336023065e-04));
    refValues[1].add(Complex(9.9999999625000000755e-01, 4.9999999820833333602e-05));
    refValues[1].add(Complex(9.9999999996250000000e-01, 4.9999999998208333333e-06));
    refValues[2].add(Complex(9.8629664648967592135e-01, 1.4909782248906273663e-01));
    refValues[2].add(Complex(9.9986250467180267611e-01, 1.4999095853247538925e-02));
    refValues[2].add(Complex(9.9999862500046718743e-01, 1.4999990958335324777e-03));
    refValues[2].add(Complex(9.9999998625000004672e-01, 1.4999999909583333532e-04));
    refValues[2].add(Complex(9.9999999986250000000e-01, 1.4999999999095833333e-05));
    refValues[3].add(Complex(9.9480250525936690737e-01, 0.0000000000000000000e+00));
    refValues[3].add(Complex(9.9994791775172400105e-01, 0.0000000000000000000e+00));
    refValues[3].add(Complex(9.9999947916677517360e-01, 0.0000000000000000000e+00));
    refValues[3].add(Complex(9.9999999479166667752e-01, 0.0000000000000000000e+00));
    refValues[3].add(Complex(9.9999999994791666667e-01, 0.0000000000000000000e+00));

    for (UnsignedInteger nTriangular = 0; nTriangular < coll.getSize(); ++nTriangular)
    {
      // Instantiate one distribution object
      Triangular distribution(coll[nTriangular]);
      fullprint << "Distribution " << distribution << std::endl;
      std::cout << "Distribution " << distribution << std::endl;

      // Is this distribution elliptical ?
      fullprint << "Elliptical = " << (distribution.isElliptical() ? "true" : "false") << std::endl;

      // Is this distribution continuous ?
      fullprint << "Continuous = " << (distribution.isContinuous() ? "true" : "false") << std::endl;

      // Test for realization of distribution
      Point oneRealization = distribution.getRealization();
      fullprint << "oneRealization=" << oneRealization << std::endl;

      // Test for sampling
      UnsignedInteger size = 10000;
      Sample oneSample = distribution.getSample( size );
      fullprint << "oneSample first=" << oneSample[0] << " last=" << oneSample[size - 1] << std::endl;
      fullprint << "mean=" << oneSample.computeMean() << std::endl;
      fullprint << "covariance=" << oneSample.computeCovariance() << std::endl;
      size = 100;
      for (UnsignedInteger i = 0; i < 2; ++i)
      {
        fullprint << "Kolmogorov test for the generator, sample size=" << size << " is " << (FittingTest::Kolmogorov(distribution.getSample(size), distribution).getBinaryQualityMeasure() ? "accepted" : "rejected") << std::endl;
        size *= 10;
      }

      // Define a point, left part of the support
      Point point( distribution.getDimension(), 1.0 );
      fullprint << "Point= " << point << std::endl;

      // Show PDF and CDF of point
      Scalar eps = 1e-5;
      Point DDF = distribution.computeDDF( point );
      fullprint << "ddf     =" << DDF << std::endl;
      Scalar LPDF = distribution.computeLogPDF( point );
      fullprint << "log pdf=" << LPDF << std::endl;
      Scalar PDF = distribution.computePDF( point );
      fullprint << "pdf     =" << PDF << std::endl;
      fullprint << "pdf (FD)=" << (distribution.computeCDF( point + Point(1, eps) ) - distribution.computeCDF( point  + Point(1, -eps) )) / (2.0 * eps) << std::endl;
      Scalar CDF = distribution.computeCDF( point );
      fullprint << "cdf=" << CDF << std::endl;
      Scalar CCDF = distribution.computeComplementaryCDF( point );
      fullprint << "ccdf=" << CCDF << std::endl;
      Scalar Survival = distribution.computeSurvivalFunction( point );
      fullprint << "survival=" << Survival << std::endl;
      Point InverseSurvival = distribution.computeInverseSurvivalFunction(0.95);
      fullprint << "Inverse survival=" << InverseSurvival << std::endl;
      fullprint << "Survival(inverse survival)=" << distribution.computeSurvivalFunction(InverseSurvival) << std::endl;
      Complex CF = distribution.computeCharacteristicFunction( point[0] );
      fullprint << "characteristic function=" << CF << std::endl;
      Complex LCF = distribution.computeLogCharacteristicFunction( point[0] );
      fullprint << "log characteristic function=" << LCF << std::endl;
      for (UnsignedInteger j = 0; j < refValues[nTriangular].getSize(); ++j)
        assert_almost_equal(distribution.computeCharacteristicFunction(u[j]), refValues[nTriangular][j]);
      try
      {
        Point PDFgr = distribution.computePDFGradient( point );
        fullprint << "pdf gradient     =" << PDFgr << std::endl;
        Point PDFgrFD(3);
        PDFgrFD[0] = (Triangular(distribution.getA() + eps, distribution.getM(), distribution.getB()).computePDF(point) -
                      Triangular(distribution.getA() - eps, distribution.getM(), distribution.getB()).computePDF(point)) / (2.0 * eps);
        PDFgrFD[1] = (Triangular(distribution.getA(), distribution.getM() + eps, distribution.getB()).computePDF(point) -
                      Triangular(distribution.getA(), distribution.getM() - eps, distribution.getB()).computePDF(point)) / (2.0 * eps);
        PDFgrFD[2] = (Triangular(distribution.getA(), distribution.getM(), distribution.getB() + eps).computePDF(point) -
                      Triangular(distribution.getA(), distribution.getM(), distribution.getB() - eps).computePDF(point)) / (2.0 * eps);
        fullprint << "pdf gradient (FD)=" << PDFgrFD << std::endl;
        Point CDFgr = distribution.computeCDFGradient( point );
        fullprint << "cdf gradient     =" << CDFgr << std::endl;
        Point CDFgrFD(3);
        CDFgrFD[0] = (Triangular(distribution.getA() + eps, distribution.getM(), distribution.getB()).computeCDF(point) -
                      Triangular(distribution.getA() - eps, distribution.getM(), distribution.getB()).computeCDF(point)) / (2.0 * eps);
        CDFgrFD[1] = (Triangular(distribution.getA(), distribution.getM() + eps, distribution.getB()).computeCDF(point) -
                      Triangular(distribution.getA(), distribution.getM() - eps, distribution.getB()).computeCDF(point)) / (2.0 * eps);
        CDFgrFD[2] = (Triangular(distribution.getA(), distribution.getM(), distribution.getB() + eps).computeCDF(point) -
                      Triangular(distribution.getA(), distribution.getM(), distribution.getB() - eps).computeCDF(point)) / (2.0 * eps);
        fullprint << "cdf gradient (FD)=" << CDFgrFD << std::endl;
      }
      catch(const NotDefinedException &)
      {
      }
      Point quantile = distribution.computeQuantile( 0.25 );
      fullprint << "quantile=" << quantile << std::endl;
      fullprint << "cdf(quantile)=" << distribution.computeCDF(quantile) << std::endl;
      // Confidence regions
      Scalar threshold;
      fullprint << "Minimum volume interval=" << distribution.computeMinimumVolumeIntervalWithMarginalProbability(0.95, threshold) << std::endl;
      fullprint << "threshold=" << threshold << std::endl;
      Scalar beta;
      LevelSet levelSet(distribution.computeMinimumVolumeLevelSetWithThreshold(0.95, beta));
      fullprint << "Minimum volume level set=" << levelSet << std::endl;
      fullprint << "beta=" << beta << std::endl;
      fullprint << "Bilateral confidence interval=" << distribution.computeBilateralConfidenceIntervalWithMarginalProbability(0.95, beta) << std::endl;
      fullprint << "beta=" << beta << std::endl;
      fullprint << "Unilateral confidence interval (lower tail)=" << distribution.computeUnilateralConfidenceIntervalWithMarginalProbability(0.95, false, beta) << std::endl;
      fullprint << "beta=" << beta << std::endl;
      fullprint << "Unilateral confidence interval (upper tail)=" << distribution.computeUnilateralConfidenceIntervalWithMarginalProbability(0.95, true, beta) << std::endl;
      fullprint << "beta=" << beta << std::endl;
      // Define a point, right part of the support
      point = Point( distribution.getDimension(), 2.0 );
      fullprint << "Point= " << point << std::endl;

      // Show PDF and CDF of point
      DDF = distribution.computeDDF( point );
      fullprint << "ddf     =" << DDF << std::endl;
      fullprint << "ddf (FD)=" << Point(1, (distribution.computePDF( point + Point(1, eps) ) - distribution.computePDF( point  + Point(1, -eps) )) / (2.0 * eps)) << std::endl;
      PDF = distribution.computePDF( point );
      fullprint << "pdf     =" << PDF << std::endl;
      fullprint << "pdf (FD)=" << (distribution.computeCDF( point + Point(1, eps) ) - distribution.computeCDF( point  + Point(1, -eps) )) / (2.0 * eps) << std::endl;
      CDF = distribution.computeCDF( point );
      fullprint << "cdf=" << CDF << std::endl;
      try
      {
        Point PDFgr = distribution.computePDFGradient( point );
        fullprint << "pdf gradient     =" << PDFgr << std::endl;
        Point PDFgrFD(3);
        PDFgrFD[0] = (Triangular(distribution.getA() + eps, distribution.getM(), distribution.getB()).computePDF(point) -
                      Triangular(distribution.getA() - eps, distribution.getM(), distribution.getB()).computePDF(point)) / (2.0 * eps);
        PDFgrFD[1] = (Triangular(distribution.getA(), distribution.getM() + eps, distribution.getB()).computePDF(point) -
                      Triangular(distribution.getA(), distribution.getM() - eps, distribution.getB()).computePDF(point)) / (2.0 * eps);
        PDFgrFD[2] = (Triangular(distribution.getA(), distribution.getM(), distribution.getB() + eps).computePDF(point) -
                      Triangular(distribution.getA(), distribution.getM(), distribution.getB() - eps).computePDF(point)) / (2.0 * eps);
        fullprint << "pdf gradient (FD)=" << PDFgrFD << std::endl;
        Point CDFgr = distribution.computeCDFGradient( point );
        fullprint << "cdf gradient     =" << CDFgr << std::endl;
        Point CDFgrFD(3);
        CDFgrFD[0] = (Triangular(distribution.getA() + eps, distribution.getM(), distribution.getB()).computeCDF(point) -
                      Triangular(distribution.getA() - eps, distribution.getM(), distribution.getB()).computeCDF(point)) / (2.0 * eps);
        CDFgrFD[1] = (Triangular(distribution.getA(), distribution.getM() + eps, distribution.getB()).computeCDF(point) -
                      Triangular(distribution.getA(), distribution.getM() - eps, distribution.getB()).computeCDF(point)) / (2.0 * eps);
        CDFgrFD[2] = (Triangular(distribution.getA(), distribution.getM(), distribution.getB() + eps).computeCDF(point) -
                      Triangular(distribution.getA(), distribution.getM(), distribution.getB() - eps).computeCDF(point)) / (2.0 * eps);
        fullprint << "cdf gradient (FD)=" << CDFgrFD << std::endl;
      }
      catch(const NotDefinedException &)
      {
      }
      quantile = distribution.computeQuantile( 0.95 );
      fullprint << "quantile=" << quantile << std::endl;
      fullprint << "cdf(quantile)=" << distribution.computeCDF(quantile) << std::endl;
      fullprint << "entropy=" << distribution.computeEntropy() << std::endl;
      fullprint << "entropy (MC)=" << -distribution.computeLogPDF(distribution.getSample(1000000)).computeMean()[0] << std::endl;
      // Moments
      Point mean = distribution.getMean();
      fullprint << "mean=" << mean << std::endl;
      Point standardDeviation = distribution.getStandardDeviation();
      fullprint << "standard deviation=" << standardDeviation << std::endl;
      Point skewness = distribution.getSkewness();
      fullprint << "skewness=" << skewness << std::endl;
      Point kurtosis = distribution.getKurtosis();
      fullprint << "kurtosis=" << kurtosis << std::endl;
      CovarianceMatrix covariance = distribution.getCovariance();
      fullprint << "covariance=" << covariance << std::endl;
      CovarianceMatrix correlation = distribution.getCorrelation();
      fullprint << "correlation=" << correlation << std::endl;
      CovarianceMatrix spearman = distribution.getSpearmanCorrelation();
      fullprint << "spearman=" << spearman << std::endl;
      CovarianceMatrix kendall = distribution.getKendallTau();
      fullprint << "kendall=" << kendall << std::endl;
      Triangular::PointWithDescriptionCollection parameters = distribution.getParametersCollection();
      fullprint << "parameters=" << parameters << std::endl;
      fullprint << "Standard representative=" << distribution.getStandardRepresentative().__str__() << std::endl;
    }
  }
  catch (const TestFailed & ex)
  {
    std::cerr << ex << std::endl;
    return ExitCode::Error;
  }


  return ExitCode::Success;
}
