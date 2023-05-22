"""
Estimate a GEV on the Fremantle sea levels data
===============================================
"""
# %%
# In this example, we illustrate various techniques of extreme value modeling applied
# on the annual maximum sea-levels recorded at Fremantle, near Perth, western Austria. 
# Readers should refer to [coles2001]_ to get more details. 
#
# We illustrate techniques to:
#
# - model a stationary and non stationary GEV,
# - estimate a return level,
#
# using:
#
# - the log-likelihhod function,
# - the profile log-likelihhod function.
#
# %%
# Load the Fremantle dataset of the annual maximum sea-levels and look at them through time.
import openturns as ot
import openturns.viewer as otv
import openturns.experimental as otexp
from openturns.usecases import coles

data = coles.Coles().fremantle
print(data[:5])
cloud = ot.Cloud(data[:, :2])
cloud.setColor("red")
graph.add(cloud)
view = otv.View(graph)

# %%
# Select the sea levels column:
sample = data[:, 1]

# %%
# **Stationary GEV modeling via the log-likelihood function**
# We first assume that the dependence through time is negligible.
# Westimate the parameters of the GEV by maximizing the log-likehood of the data
# and compute the resulting GEV model:
factory = ot.GeneralizedExtremeValueFactory()
resultLL = factory.buildMethodOfLikelihoodMaximizationEstimator(sample)

# %%
# We get the estimated parameter of :math:`(\mu, \sigma, \xi)`:
estimate = resultLL.getDistribution()
desc = estimate.getParameterDescription()
p = estimate.getParameter()
print(", ".join([f"{param}: {value:.3f}" for param, value in zip(desc, p)]))

# %%
# We get the asymptotic distribution of the estimator :math:`(\hat{\mu}, \hat{\sigma}, \hat{\xi})`.
# In that case, the asymptotic distribution is normal.
parameterEstimate = resultLL.getParameterDistribution()
print(parameterEstimate)

# %%
# We get the covariance matrix  and the standard deviation of :math:`(\mu, \sigma, \xi)`:
print('Cov matrix = ', parameterEstimate.getCovariance())
print('standard dev = ', parameterEstimate.getStandardDeviation())

# %%
# We get the marginal confidence intervals of order 0.95:
for i in range(3):
    ci = parameterEstimate.getMarginal(i).computeBilateralConfidenceInterval(0.95)
    print(desc[i] + ":", ci)

# %%
# At last, we can validate the inference result thanks the 4 usual diagnostic plots:
validation = otexp.GeneralizedExtremeValueValidation(resultLL, sample)
graph = validation.drawDiagnosticPlot()
view = otv.View(graph)

# %%
# **Stationary GEV modeling via the profile log-likelihood function**
# Now, we use the profile log-likehood function rather than log-likehood function  to estimate the parameters of the GEV:
resultPLL = factory.buildMethodOfProfileLikelihoodMaximizationEstimator(sample)

# %%
# That plot allows to get the profile log-likelihood plot and the confidence interval for :math:`\xi` of order 0.95. 
# The graph also indicates the optimal value of :math:`\xi`, the maximum profile log-likelihood and
# the confidence interval.
view = otv.View(resultPLL.drawProfileLikelihoodFunction())

# %%
# We can get the numerical values of the confidence interval: it appears to be a bit smaller
# with the profile log-likelihood function than with the log-likelihood function.
resultPLL.setConfidenceLevel(0.95)
print(resultPLL.getParameterConfidenceInterval())

# %%
# **Return level estimate from the estimated stationary GEV**
# We estimate the :math:`m`-block return level :math:`z_m`: it is computed as a particular quantile of the
# GEV model estimated using the log-likelihood function. We just have to use the maximum log-likelihood
# estimator built in the previous section.
# As the data are annual sea-levels, each block corresponds to one year: the 10-year return level
# corresponds to :math:`m=10` and the 100-year return level corresponds to :math:`m=100`.
# The method also provides the asymptotic distribution of the estimator :math:`\hat{z}_m`.
zm10 = factory.buildReturnLevelEstimator(resultLL, 10.0)
return_level10 = zm10.getMean()
print(f"10 years return level={return_level10}")
return_level_ci10 = zm10.computeBilateralConfidenceInterval(0.95)
print(f"CI = {return_level_ci10}")

zm100 = factory.buildReturnLevelEstimator(resultLL, 100.0)
return_level100 = zm100.getMean()
print(f"100 years return level={return_level100}")
return_level_ci100 = zm100.computeBilateralConfidenceInterval(0.95)
print(f"CI = {return_level_ci100}")

# %%
# **Return level estimate via the profile log-likelihood function of a stationary GEV**
# We can estimate the :math:`m`-block return level :math:`z_m` directly from the data using the profile
# likelihood with respect to :math:`z_m`.
result_rl10_prof = factory.buildReturnLevelProfileLikelihoodEstimator(sample, 10.0)
zm = result_rl10_prof.getParameter()
print(f"10 years return level (profile) = {zm}")

# %%
# We can get the confidence interval of :math:`z_m`:  one cmore, it appears to be a bit smaller
# with the profile log-likelihood function than with the log-likelihood function
result_rl10_prof.setConfidenceLevel(0.95)
return_level_ci10 = result_rl10_prof.getParameterConfidenceInterval()
print(f"CI={return_level_ci10}")

# %%
# We can also plot the profile log-likelihood and read the confidence interval, the optimal value
# of :math:`z_m` and its confidnece interval:
view = otv.View(result_rl10_prof.drawProfileLikelihoodFunction())

# %%
# **Non stationary GEV modeling via the log-likelihood function**
# If we look at the data carefully, we see that the pattern of variation has not remained constant over
# the observation period. There is an increase in the data through time.
# We want to model this trend becaucse a slight increase in extreme sea-levels could have
# significant impact on the safety of coastal flood defenses.
# 
# First we need to get the grid of time values (in years here):
mesh = ot.Mesh(data[:, 0])

# %%
# Then, we define the functional basis for each parameter of the GEV model. Even if we have
# the possibility to affect a time-varying model to each of the 3 parameters :math:`(\mu, \sigma, \xi)`,
# it is strongly recommended not to vary the parameter :math:`\xi`.
# We suppose that :math:`\mu` is linear with time, and that the other parameters remain constant:
# 
# ..math::
#    nowrap
#
#    \begin{align*}
#    \mu(t) & = \beta_1 + \beta_2t \\
#    \sigma(t) & = \beta_3 \\
#    \xi(t) & = \beta_4
#    \end{align*}
constant = ot.SymbolicFunction(["t"], ["1.0"])
basis_mu = ot.Basis([constant, ot.SymbolicFunction(["t"], ["t"])]) 
basis_sigma = ot.Basis([constant])  
basis_xi = ot.Basis([constant])  
basis_coll = [basis_mu, basis_sigma, basis_xi]

# %%
# We can now estimate the list of coefficients :math:`(\beta_1, \beta_2, \beta_3, \beta_4)` using the log-likelihood of the data:
resultNonStatLL = factory.buildTimeVarying(sample, mesh, basis_coll)
beta = resultNonStatLL.getOptimalParameter()
beta_1, beta_2 = beta[:2]
print('beta1, beta2, beta3, beta4 = ', beta)
print(f"mu(t) = {beta_1:.4f} + {beta_2:.4f} * t")

# %%
# We get the asymptotic distribution of :math:`\vect{\beta}` to compute some confidence intervals of
# the estimates, for exemple of order :math:`p = 0.95`.
dist_beta = resultNonStatLL.getParameterDistribution()
condifence level = 0.95
for i in range(beta.getSize()):
    lower_bound = dist_beta.getMarginal(i).computeQuantile((1-condifence level)/2)[0]
    upper_bound = dist_beta.getMarginal(i).computeQuantile((1+condifence level)/2)[0]
    print('Conf interval for beta_' + str(i+1) + ' = [' + str(lower_bound) + '; ' + str(upper_bound) + ']')   

# %%
# In order to compare different modelings, we get the optimal log-likelihood of the data:
optim_LL = resultNonStatLL.getLogLikelihood()
print('Max log-likelihood = ', optim_LL)

# %%
# We can draw the estimated trend  :math:`t \mapsto \mu(t)` with the data: the graph confirms the increase of the annual maximum sea-levels through time.
graph = resultNonStatLL.drawParameterFunction(0)
cloud = ot.Cloud(data[:, :2])
cloud.setColor("red")
graph.add(cloud)
view = otv.View(graph)

# %%
# In order to test the validity of the model with time varying parameters (model :math:`\mathcal{M}_1`) relative to the stationary model (model :math:`\mathcal{M}_0 \subset \mathcal{M}_1`), we use the Likelihood Ratio test.
# The null hypothesis is the stationary model :math:`\mathcal{M}_0`. The Type 1 error :math:`\alpha`
# is taken equal to 0.05. 
# This test confirms that the dependence through time is not negligible: it means that the linear trend component explains a large variation in the data.
llh_LL = resultLL.getLogLikelihood()
llh_NonStatLL = resultNonStatLL.getLogLikelihood()
resultLikRatioTest = ot.HypothesisTest.LikelihoodRatioTest(llh_LL, llh_NonStatLL, 0.05)
accepted = resultLikRatioTest.getBinaryQualityMeasure()
print(f"hypothesis H0(stationary model) accepted={accepted}")

# %%
# We detail the statistics of the Likelihood Ratio test: the deviance statistics :math:`\mathcal{D}_p` follows a :math:`\cih^2_1` distribution.
# The model :math:`\mathcal{M}_0` is rejected of the deviance statistics estimated on the data is greater than the p-value :math:`c_{\alpha}` corresponding to :math:`\alpha`.
print(f"Dp={result4.getStatistic():.2f}")
print(f"cAlpha={result4.getThreshold():.2f}")

# %%
# Readers should do the same study with a quadratic trend in :math:`\mu` or a linear trend in
# :math:`\sigma`.
#
# We can also draw the function :math:`t \mapsto q_p(t)` where :math:`q_p(t)` is the quantile of order :math:`p` of the GEV distribution at time :math:`t`. 
# Here, :math:`mu(t)` is a linear function and the other parameters are constant, so the quantile function is also a linear function in time. 
graph = resultNonStatLL.drawQuantileFunction(0.5)
view = otv.View(graph)

# %%
otv.View.ShowAll()
