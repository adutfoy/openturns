//                                               -*- C++ -*-
/**
 *  @brief The MeixnerDistribution distribution
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
#include <cmath>
#include "openturns/MeixnerDistribution.hxx"
#include "openturns/RandomGenerator.hxx"
#include "openturns/SpecFunc.hxx"
#include "openturns/PersistentObjectFactory.hxx"
#include "openturns/TNC.hxx"
#include "openturns/Brent.hxx"
#include "openturns/GaussKronrod.hxx"
#include "openturns/FunctionImplementation.hxx"
#include "openturns/OptimizationAlgorithm.hxx"
#include "openturns/Distribution.hxx"

BEGIN_NAMESPACE_OPENTURNS

CLASSNAMEINIT(MeixnerDistribution)

static const Factory<MeixnerDistribution> Factory_MeixnerDistribution;

/* Default constructor */
MeixnerDistribution::MeixnerDistribution()
  : DistributionImplementation()
  , solver_(new TNC())
  , beta_(0.0)
  , alpha_(0.0)
  , delta_(0.0)
  , gamma_(0.0)
  , logNormalizationFactor_(0.0)
{
  setName("MeixnerDistribution");
  // Create the optimization solver parameters using the parameters in the ResourceMap
  initializeOptimizationAlgorithmParameter();
  setBetaAlphaDelta(1.0, 0.0, 1.0);
  setDimension(1);
}

/* Parameters constructor */
MeixnerDistribution::MeixnerDistribution(const Scalar beta,
    const Scalar alpha,
    const Scalar delta,
    const Scalar mu)
  : DistributionImplementation()
  , solver_(new TNC())
  , beta_(0.0)
  , alpha_(0.0)
  , delta_(0.0)
  , gamma_(mu)
  , logNormalizationFactor_(0.0)
{
  setName("MeixnerDistribution");
  // Create the optimization solver parameters using the parameters in the ResourceMap
  initializeOptimizationAlgorithmParameter();
  setBetaAlphaDelta(beta, alpha, delta);
  setDimension(1);
}

/* Initialize optimization solver parameter using the ResourceMap */
void MeixnerDistribution::initializeOptimizationAlgorithmParameter()
{
  solver_.setMaximumAbsoluteError(ResourceMap::GetAsScalar("MeixnerDistribution-MaximumAbsoluteError"));
  solver_.setMaximumRelativeError(ResourceMap::GetAsScalar("MeixnerDistribution-MaximumRelativeError"));
  solver_.setMaximumResidualError(ResourceMap::GetAsScalar("MeixnerDistribution-MaximumObjectiveError"));
  solver_.setMaximumConstraintError(ResourceMap::GetAsScalar("MeixnerDistribution-MaximumConstraintError"));
}

/* Comparison operator */
Bool MeixnerDistribution::operator ==(const MeixnerDistribution & other) const
{
  if (this == &other) return true;
  return (beta_ == other.beta_) && (alpha_ == other.alpha_) && (delta_ == other.delta_) && (gamma_ == other.gamma_);
}

Bool MeixnerDistribution::equals(const DistributionImplementation & other) const
{
  const MeixnerDistribution* p_other = dynamic_cast<const MeixnerDistribution*>(&other);
  return p_other && (*this == *p_other);
}

/* String converter */
String MeixnerDistribution::__repr__() const
{
  OSS oss(true);
  oss << "class=" << MeixnerDistribution::GetClassName()
      << " name=" << getName()
      << " dimension=" << getDimension()
      << " beta=" << beta_
      << " alpha=" << alpha_
      << " delta=" << delta_
      << " gamma=" << gamma_
      << " logNormalizationFactor=" << logNormalizationFactor_
      << " b=" << b_
      << " c=" << c_
      << " dc=" << dc_;
  return oss;
}

String MeixnerDistribution::__str__(const String & ) const
{
  OSS oss(false);
  oss << getClassName() << "(beta = " << beta_ << ", alpha = " << alpha_ << ", delta = " << delta_ << ", gamma = " << gamma_ << ")";
  return oss;
}

/* Beta accessor */
void MeixnerDistribution::setBeta(const Scalar beta)
{
  if (!(beta > 0.0)) throw InvalidArgumentException(HERE) << "Beta MUST be positive";
  if (beta != beta_)
  {
    beta_ = beta;
    update();
  }
}

Scalar MeixnerDistribution::getBeta() const
{
  return beta_;
}


/* Alpha accessor */
void MeixnerDistribution::setAlpha(const Scalar alpha)
{
  if (std::abs(alpha) >= M_PI) throw InvalidArgumentException(HERE) << "Alpha MUST be in (-pi, pi)";
  if (alpha != alpha_)
  {
    alpha_ = alpha;
    update();
  }
}

Scalar MeixnerDistribution::getAlpha() const
{
  return alpha_;
}

/* Delta accessor */
void MeixnerDistribution::setDelta(const Scalar delta)
{
  if (!(delta > 0.0)) throw InvalidArgumentException(HERE) << "Delta MUST be positive";
  if (delta != delta_)
  {
    delta_ = delta;
    update();
  }
}

Scalar MeixnerDistribution::getDelta() const
{
  return delta_;
}

void MeixnerDistribution::setBetaAlphaDelta(const Scalar beta,
    const Scalar alpha,
    const Scalar delta)
{
  if (!(beta > 0.0)) throw InvalidArgumentException(HERE) << "Beta MUST be positive";
  if (std::abs(alpha) >= M_PI) throw InvalidArgumentException(HERE) << "Alpha MUST be in (-pi, pi)";
  if (!(delta > 0.0)) throw InvalidArgumentException(HERE) << "Delta MUST be positive";
  if ((beta != beta_) || (alpha != alpha_) || (delta != delta_))
  {
    beta_ = beta;
    alpha_ = alpha;
    delta_ = delta;
    update();
  }
}

/* Mu accessor */
void MeixnerDistribution::setGamma(const Scalar gamma)
{
  if (gamma_ != gamma)
  {
    gamma_ = gamma;
    computeRange();
  }
}

Scalar MeixnerDistribution::getGamma() const
{
  return gamma_;
}

/* Virtual constructor */
MeixnerDistribution * MeixnerDistribution::clone() const
{
  return new MeixnerDistribution(*this);
}

/* Compute the numerical range of the distribution given the parameters values */
void MeixnerDistribution::computeRange()
{
  const Point mu(getMean());
  const Point sigma(getStandardDeviation());
  const Scalar logPDF = computeLogPDF(mu);
  const Scalar logPDFEpsilon = std::log(getPDFEpsilon());
  Point lowerBound(mu);
  // Find the numerical lower bound based on the PDF value
  Scalar logPDFLower = logPDF;
  while (logPDFLower > logPDFEpsilon)
  {
    lowerBound -= sigma;
    logPDFLower = computeLogPDF(lowerBound);
  }
  // Find the numerical upper bound based on the PDF value
  Point upperBound(mu);
  Scalar logPDFUpper = logPDF;
  while (logPDFUpper > logPDFEpsilon)
  {
    upperBound += sigma;
    logPDFUpper = computeLogPDF(upperBound);
  }
  const Interval::BoolCollection finiteLowerBound(1, false);
  const Interval::BoolCollection finiteUpperBound(1, false);
  setRange(Interval(lowerBound, upperBound, finiteLowerBound, finiteUpperBound));
}


namespace
{

class MeixnerBoundB: public FunctionImplementation
{
public:
  MeixnerBoundB(const DistributionImplementation::Implementation & p_distribution)
    : FunctionImplementation()
    , p_distribution_(p_distribution)
  {
    // Nothing to do
  }

  MeixnerBoundB * clone() const override
  {
    return new MeixnerBoundB(*this);
  }

  Point operator() (const Point & point) const override
  {
    return Point(1, p_distribution_->computePDF(point));
  }

  Matrix gradient(const Point & point) const override
  {
    const Point value = p_distribution_->computeDDF(point);
    return MatrixImplementation(getInputDimension(), getOutputDimension(), value);
  }

  UnsignedInteger getInputDimension() const override
  {
    return p_distribution_->getDimension();
  }

  UnsignedInteger getOutputDimension() const override
  {
    return 1;
  }

private:
  const DistributionImplementation::Implementation p_distribution_;
}; // class MeixnerBoundB

class MeixnerBoundCD: public FunctionImplementation
{
public:
  MeixnerBoundCD(const DistributionImplementation::Implementation & p_distribution)
    : FunctionImplementation()
    , p_distribution_(p_distribution)
  {
    // Nothing to do
  }

  MeixnerBoundCD * clone() const override
  {
    return new MeixnerBoundCD(*this);
  }

  Point operator() (const Point & point) const override
  {
    const Scalar pdf = p_distribution_->computePDF(point);
    return Point(1, point[0] * std::sqrt(pdf));
  }

  Matrix gradient(const Point & point) const override
  {
    const Scalar sqrtPDF = std::sqrt(p_distribution_->computePDF(point));
    if (sqrtPDF <= 0.0) return MatrixImplementation(1, 1);
    const Point ddf(p_distribution_->computeDDF(point));
    const Point value(1, sqrtPDF + 0.5 * point[0] * ddf[0] / sqrtPDF);
    return MatrixImplementation(1, 1, value);
  }

  UnsignedInteger getInputDimension() const override
  {
    return p_distribution_->getDimension();
  }

  UnsignedInteger getOutputDimension() const override
  {
    return 1;
  }

private:
  const DistributionImplementation::Implementation p_distribution_;
}; // class MeixnerBoundCD

// Here we inherit from EvaluationImplementation in order to benefit from the inherited evaluation over a sample
class MeixnerPDFGradientWrapper: public EvaluationImplementation
{
public:
  MeixnerPDFGradientWrapper(const DistributionImplementation::Implementation & p_distribution)
    : EvaluationImplementation()
    , p_distribution_(p_distribution)
  {
    // Nothing to do
  }

  MeixnerPDFGradientWrapper * clone() const override
  {
    return new MeixnerPDFGradientWrapper(*this);
  }

  Point operator() (const Point & point) const override
  {
    return p_distribution_->computePDFGradient(point);
  }

  UnsignedInteger getInputDimension() const override
  {
    return 1;
  }

  UnsignedInteger getOutputDimension() const override
  {
    return 4;
  }

private:
  const DistributionImplementation::Implementation p_distribution_;
}; // class MeixnerBoundCD


} // anonymous namespace

/* Update the derivative attributes */
void MeixnerDistribution::update()
{
  // First, the parameters of the distribution
  logNormalizationFactor_ = 2.0 * delta_ * std::log(2.0 * std::cos(0.5 * alpha_)) - std::log(2.0 * M_PI * beta_) - SpecFunc::LogGamma(2.0 * delta_);
  computeRange();
  // Second, the moments
  isAlreadyComputedMean_ = false;
  isAlreadyComputedCovariance_ = false;
  // Third, the CDF and complementary CDF
  Collection<PiecewiseHermiteEvaluation> CDFCCDF(interpolateCDF(ResourceMap::GetAsUnsignedInteger("MeixnerDistribution-CDFDiscretization")));
  cdfApproximation_ = CDFCCDF[0];
  ccdfApproximation_ = CDFCCDF[1];
  // Fourth, the random generator

  // define objectives functions
  MeixnerBoundB fB(clone());
  MeixnerBoundCD fCD(clone());

  // Define Optimization problem1 : maximization fB
  OptimizationProblem problem1(fB);
  problem1.setMinimization(false);
  problem1.setObjective(fB);
  problem1.setBounds(getRange());
  OptimizationAlgorithm solver(solver_);
  solver.setStartingPoint(getMean());
  solver.setProblem(problem1);
  solver.run();
  b_ = std::sqrt(solver.getResult().getOptimalValue()[0]);

  // Define Optimization problem2 : minimization fCD
  OptimizationProblem problem2(fCD);
  problem2.setMinimization(true);
  solver.setProblem(problem2);
  solver.run();
  c_ = solver.getResult().getOptimalValue()[0];

  // Define Optimization problem3 : maximization fCD
  OptimizationProblem problem3(fCD);
  problem3.setMinimization(false);
  solver.setProblem(problem3);
  solver.run();
  dc_ = solver.getResult().getOptimalValue()[0] - c_;
}

/* Build a C1 interpolation of the CDF function */
Collection<PiecewiseHermiteEvaluation> MeixnerDistribution::interpolateCDF(const UnsignedInteger n)
{
  const Function pdf(getPDF());
  const Scalar xMin = getRange().getLowerBound()[0];
  const Scalar xMax = getRange().getUpperBound()[0];
  const Scalar mu = getMean()[0];
  Point locationsCDF(n);
  Point locationsCCDF(n);
  Point valuesCDF(n);
  Point valuesCCDF(n);
  Point derivativesCDF(n);
  Point derivativesCCDF(n);
  Scalar xCDFOld = xMin;
  Scalar xCCDFOld = xMax;
  locationsCDF[0] = xMin;
  locationsCCDF[n - 1] = xMax;
  GaussKronrod algo;
  const Scalar stepCDF = (mu - xMin) / (n - 1.0);
  const Scalar stepCCDF = (xMax - mu) / (n - 1.0);
  for (UnsignedInteger i = 1; i < n; ++i)
  {
    const Scalar xCDF = xMin + i * stepCDF;
    const Scalar xCCDF = xMax - i * stepCCDF;
    locationsCDF[i] = xCDF;
    locationsCCDF[n - i - 1] = xCCDF;
    Point ai;
    Point bi;
    Sample fi;
    Point ei;
    Scalar error = -1.0;
    valuesCDF[i] = valuesCDF[i - 1] + algo.integrate(pdf, xCDFOld, xCDF, error, ai, bi, fi, ei)[0];
    valuesCCDF[n - i - 1] = valuesCCDF[n - i] + algo.integrate(pdf, xCCDF, xCCDFOld, error, ai, bi, fi, ei)[0];
    derivativesCDF[i] = computePDF(xCDF);
    derivativesCCDF[n - i - 1] = -computePDF(xCCDF);
    xCDFOld = xCDF;
    xCCDFOld = xCCDF;
  }
  Collection<PiecewiseHermiteEvaluation> coll(2);
  coll[0] = PiecewiseHermiteEvaluation(locationsCDF, valuesCDF, derivativesCDF);
  coll[1] = PiecewiseHermiteEvaluation(locationsCCDF, valuesCCDF, derivativesCCDF);
  return coll;
}

/* Get one realization of the distribution
   We use the ratio-of-uniform method:
   if b=\sup_x \sqrt{p(x)}, c=\inf_x x\sqrt{p(x)}, d=\sup_x x\sqrt{p(x)}
   then \{(u,v)\,|\,0\leq u\leq b,c\leq v\leq d\} is an enclosing region for
   \{(u,v)\,|\,0\leq u\leq\sqrt{p(u/v)}\}
*/
Point MeixnerDistribution::getRealization() const
{
  while (true)
  {
    const Scalar u = b_ * RandomGenerator::Generate();
    if (u == 0.0) continue;
    const Scalar v = c_ + dc_ * RandomGenerator::Generate();
    const Scalar rho = v / u;
    if (2.0 * std::log(u) <= computeLogPDF(rho)) return Point(1, rho);
  }
}

/* Get the PDF of the distribution */
Scalar MeixnerDistribution::computePDF(const Point & point) const
{
  if (point.getDimension() != 1) throw InvalidArgumentException(HERE) << "Error: the given point must have dimension=1, here dimension=" << point.getDimension();

  return std::exp(computeLogPDF(point));
}

Scalar MeixnerDistribution::computeLogPDF(const Point & point) const
{
  if (point.getDimension() != 1) throw InvalidArgumentException(HERE) << "Error: the given point must have dimension=1, here dimension=" << point.getDimension();

  const Scalar z = (point[0] - gamma_) / beta_;
  return logNormalizationFactor_ + alpha_ * z + 2.0 * SpecFunc::LogGamma(Complex(delta_, z)).real();
}

/* Get the CDF of the distribution */
Scalar MeixnerDistribution::computeCDF(const Point & point) const
{
  if (point.getDimension() != 1) throw InvalidArgumentException(HERE) << "Error: the given point must have dimension=1, here dimension=" << point.getDimension();
  if (point[0] <= getMean()[0]) return cdfApproximation_(point)[0];
  else return 1.0 - ccdfApproximation_(point)[0];
}

/* Get the complementary CDF of the distribution */
Scalar MeixnerDistribution::computeComplementaryCDF(const Point & point) const
{
  if (point.getDimension() != 1) throw InvalidArgumentException(HERE) << "Error: the given point must have dimension=1, here dimension=" << point.getDimension();
  if (point[0] <= getMean()[0]) return 1.0 - cdfApproximation_(point)[0];
  else return ccdfApproximation_(point)[0];
}

/* Compute the quantile function of the distribution */
Scalar MeixnerDistribution::computeScalarQuantile(const Scalar prob,
    const Bool tail) const
{
  if (!((prob >= 0.0) && (prob <= 1.0)))
    throw InvalidArgumentException(HERE) << "computeScalarQuantile expected prob to belong to [0,1], but is " << prob;
  const Scalar a = getRange().getLowerBound()[0];
  const Scalar b = getRange().getUpperBound()[0];
  if (prob <= 0.0) return (tail ? b : a);
  if (prob >= 1.0) return (tail ? a : b);
  const UnsignedInteger n = cdfApproximation_.getLocations().getSize();
  if (tail)
  {
    // Here we have to solve ComplementaryCDF(x) = prob which is mathematically
    // equivalent to CDF(x) = 1 - prob, but numerically different with an
    // accuracy that depends on prob.
    // The cut-off is around the mean value
    if (prob <= ccdfApproximation_.getValues()(0, 0)) return Brent(quantileEpsilon_, cdfEpsilon_, cdfEpsilon_, quantileIterations_).solve(ccdfApproximation_, prob, ccdfApproximation_.getLocations()[0], ccdfApproximation_.getLocations()[n - 1], ccdfApproximation_.getValues()(0, 0), ccdfApproximation_.getValues()(n - 1, 0));
    return Brent(quantileEpsilon_, cdfEpsilon_, cdfEpsilon_, quantileIterations_).solve(cdfApproximation_, 1.0 - prob, cdfApproximation_.getLocations()[0], cdfApproximation_.getLocations()[n - 1], cdfApproximation_.getValues()(0, 0), cdfApproximation_.getValues()(n - 1, 0));
  }
  // Here we have to solve CDF(x) = prob which is mathematically
  // equivalent to ComplementaryCDF(x) = 1 - prob, but numerically
  // different with an accuracy that depends on prob.
  // The cut-off is around the mean value
  if (prob <= cdfApproximation_.getValues()(n - 1, 0)) return Brent(quantileEpsilon_, cdfEpsilon_, cdfEpsilon_, quantileIterations_).solve(cdfApproximation_, prob, cdfApproximation_.getLocations()[0], cdfApproximation_.getLocations()[n - 1], cdfApproximation_.getValues()(0, 0), cdfApproximation_.getValues()(n - 1, 0));
  return Brent(quantileEpsilon_, cdfEpsilon_, cdfEpsilon_, quantileIterations_).solve(ccdfApproximation_, 1.0 - prob, ccdfApproximation_.getLocations()[0], ccdfApproximation_.getLocations()[n - 1], ccdfApproximation_.getValues()(0, 0), ccdfApproximation_.getValues()(n - 1, 0));
}

/* Get the characteristic function of the distribution, i.e. phi(u) = E(exp(I*u*X)) */
Complex MeixnerDistribution::computeCharacteristicFunction(const Scalar x) const
{
  return std::exp(computeLogCharacteristicFunction(x));
}

Complex MeixnerDistribution::computeLogCharacteristicFunction(const Scalar x) const
{
  return Complex(2.0 * delta_ * std::log(std::cos(0.5 * alpha_)), gamma_ * x) - 2.0 * delta_ * std::log(std::cosh(Complex(0.5 * beta_ * x, -0.5 * gamma_)));
}

/* Compute the mean of the distribution */
void MeixnerDistribution::computeMean() const
{
  mean_ = Point(1, beta_ * delta_ * std::tan(0.5 * alpha_) + gamma_);
  isAlreadyComputedMean_ = true;
}

/* Get the PDFGradient of the distribution */
Point MeixnerDistribution::computePDFGradient(const Point & point) const
{
  if (point.getDimension() != 1) throw InvalidArgumentException(HERE) << "Error: the given point must have dimension=1, here dimension=" << point.getDimension();

  const Scalar z = (point[0] - gamma_) / beta_;
  const Complex zPlus(delta_, z);
  const Complex zMinus(delta_, -z);
  const Complex gammaPlus(SpecFunc::Gamma(zPlus));
  const Complex gammaMinus(SpecFunc::Gamma(zMinus));
  const Complex gammaProduct(gammaPlus * gammaMinus);
  const Complex psiPlus(SpecFunc::DiGamma(zPlus));
  const Complex psiMinus(SpecFunc::DiGamma(zMinus));
  const Complex iZ(0.0, z);
  const Complex i(0.0, 1.0);
  const Scalar cosHalfAlpha = std::cos(0.5 * alpha_);
  const Scalar powCosHalfAlpha = std::pow(cosHalfAlpha, 2.0 * delta_);
  const Scalar sinHalfAlpha = std::sin(0.5 * alpha_);
  const Scalar expAlphaZ = std::exp(alpha_ * z);
  const Scalar powTwo = std::pow(2.0, 2.0 * delta_ - 1.0);
  const Scalar gammaTwoDelta = SpecFunc::Gamma(2.0 * delta_);
  const Scalar denominator = M_PI * beta_ * gammaTwoDelta;
  const Complex numerator(powTwo * expAlphaZ * powCosHalfAlpha * gammaProduct);
  Point pdfGradient(4, 0.0);
  pdfGradient[0] = -(numerator * (iZ * (psiPlus - psiMinus) + 1.0 + alpha_ * z)).real() / (beta_ * denominator);
  pdfGradient[1] = -numerator.real() * (sinHalfAlpha * delta_ - cosHalfAlpha * z) / (cosHalfAlpha * denominator);
  pdfGradient[2] = -(numerator * (2.0 * SpecFunc::Psi(2.0 * delta_) - 2.0 * std::log(cosHalfAlpha) - psiMinus - psiPlus - 2.0 * M_LN2)).real() / denominator;
  pdfGradient[3] =  (numerator * i * (psiMinus - psiPlus + i * alpha_)).real() / (beta_ * denominator);
  return pdfGradient;
}


/* Get the CDFGradient of the distribution */
Point MeixnerDistribution::computeCDFGradient(const Point & point) const
{
  if (point.getDimension() != 1) throw InvalidArgumentException(HERE) << "Error: the given point must have dimension=1, here dimension=" << point.getDimension();

  return GaussKronrod().integrate(MeixnerPDFGradientWrapper(clone()), Interval(range_.getLowerBound(), point));
}

/* Get the standard deviation of the distribution */
Point MeixnerDistribution::getStandardDeviation() const
{
  return Point(1, beta_ * std::sqrt(delta_ / (1.0 + std::cos(alpha_))));
}

/* Get the skewness of the distribution */
Point MeixnerDistribution::getSkewness() const
{
  return Point(1, std::sin(0.5 * alpha_) * std::sqrt(2.0 / delta_));
}

/* Get the kurtosis of the distribution */
Point MeixnerDistribution::getKurtosis() const
{
  return Point(1, 3.0 + (2.0 - std::cos(alpha_)) / delta_);
}

/* Get the standard representative in the parametric family, associated with the standard moments */
Distribution MeixnerDistribution::getStandardRepresentative() const
{
  MeixnerDistribution standard(1.0, alpha_, delta_, 0.0);
  standard.setDescription(getDescription());
  return standard;
}

/* Compute the covariance of the distribution */
void MeixnerDistribution::computeCovariance() const
{
  covariance_ = CovarianceMatrix(1);
  covariance_(0, 0) = beta_ * beta_ * delta_ / (1.0 + std::cos(alpha_));
  isAlreadyComputedCovariance_ = true;
}

/* Parameters value and description accessor */
Point MeixnerDistribution::getParameter() const
{
  return {beta_, alpha_, delta_, gamma_};
}

void MeixnerDistribution::setParameter(const Point & parameter)
{
  if (parameter.getSize() != 4) throw InvalidArgumentException(HERE) << "Error: expected 4 values, got " << parameter.getSize();
  const Scalar w = getWeight();
  *this = MeixnerDistribution(parameter[0], parameter[1], parameter[2], parameter[3]);
  setWeight(w);
}

/* Parameters description accessor */
Description MeixnerDistribution::getParameterDescription() const
{
  return {"beta", "alpha", "delta", "gamma"};
}

/* Check if the distribution is elliptical */
Bool MeixnerDistribution::isElliptical() const
{
  return alpha_ == 0.0;
}

/* Method save() stores the object through the StorageManager */
void MeixnerDistribution::save(Advocate & adv) const
{
  DistributionImplementation::save(adv);
  adv.saveAttribute( "beta_", beta_ );
  adv.saveAttribute( "alpha_", alpha_ );
  adv.saveAttribute( "delta_", delta_ );
  adv.saveAttribute( "gamma_", gamma_ );
  adv.saveAttribute( "logNormalizationFactor_", logNormalizationFactor_ );
}

/* Method load() reloads the object from the StorageManager */
void MeixnerDistribution::load(Advocate & adv)
{
  DistributionImplementation::load(adv);
  adv.loadAttribute( "beta_", beta_ );
  adv.loadAttribute( "alpha_", alpha_ );
  adv.loadAttribute( "delta_", delta_ );
  if (adv.hasAttribute("mu_")) // old parameter set
  {
    adv.loadAttribute( "mu_", gamma_ );

    // swap alpha, beta
    const Scalar tmp = beta_;
    beta_ = alpha_;
    alpha_ = tmp;
  }
  else
    adv.loadAttribute( "gamma_", gamma_ );
  adv.loadAttribute( "logNormalizationFactor_", logNormalizationFactor_ );
  update();
}

/* Optimization solver accessor */
OptimizationAlgorithm MeixnerDistribution::getOptimizationAlgorithm() const
{
  return solver_;
}

void MeixnerDistribution::setOptimizationAlgorithm(const OptimizationAlgorithm & solver)
{
  solver_ = solver;
}

END_NAMESPACE_OPENTURNS
