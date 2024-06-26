.. _sensitivity_ancova:

Sensivity analysis with correlated inputs
-----------------------------------------

The ANCOVA (ANalysis of COVAriance) method, is a variance-based method
generalizing the ANOVA (ANalysis Of VAriance) decomposition for models
with correlated input parameters.

Let us consider a model :math:`Y = \model(\inputRV)` without making any
hypothesis on the dependence structure of
:math:`\inputRV = \{X^1, \ldots, X^\inputDim\}`, a :math:`d`-dimensional
random vector. The covariance decomposition requires a functional
decomposition of the model. Thus the model response :math:`Y` is
expanded as a sum of functions of increasing dimension as follows:

.. math::
  :label: Model

       \model(\inputRV) = \model_0 + \sum_{u\subseteq\{1,\dots,d\}} \model_u(X_u)

:math:`\model_0` is the mean of :math:`Y`. Each function :math:`\model_u`
represents, for any non empty set :math:`u\subseteq\{1, \dots, d\}`,
the combined contribution of the variables :math:`X_u` to :math:`Y`.

Using the properties of the covariance, the variance of :math:`Y` can be
decomposed into a variance part and a covariance part as follows:

.. math::

   \begin{aligned}
       Var[Y] &=& Cov\left[\model_0 + \sum_{u\subseteq\{1,\dots,d\}} \model_u(X_u), \model_0 + \sum_{u\subseteq\{1,\dots,n\}} \model_u(X_u)\right] \\
              &=& \sum_{u\subseteq\{1,\dots,d\}} Cov\left[\model_u(X_u), \sum_{u\subseteq\{1,\dots,d\}} \model_u(X_u)\right] \\
              &=& \sum_{u\subseteq\{1,\dots,d\}} \left[Var[\model_u(X_u)] + Cov[\model_u(X_u), \sum_{v\subseteq\{1,\dots,d\}, v\cap u=\varnothing} \model_v(X_v)]\right]
     \end{aligned}

The total part of variance of :math:`Y` due to :math:`X_u` reads:

.. math:: S_u = \frac{Cov[Y, \model_u(X_u)]}{Var[Y]}

The variance formula described above enables to define each sensitivity
measure :math:`S_u` as the sum of a :math:`\mathit{physical}` (or
:math:`\mathit{uncorrelated}`) part and a :math:`\mathit{correlated}`
part such as:

.. math:: S_u = S_u^U + S_u^C

where :math:`S_u^U` is the uncorrelated part of variance of :math:`Y`
due to :math:`X_u`:

.. math:: S_u^U = \frac{Var[\model_u(X_u)]}{Var[Y]}

and :math:`S_u^C` is the contribution of the correlation of :math:`X_u`
with the other parameters:

.. math:: S_u^C = \frac{Cov[\model_u(X_u), \displaystyle \sum_{v\subseteq\{1,\dots,d\}, v\cap u=\varnothing} \model_v(X_v)]}{Var[Y]}

As the computational cost of the indices with the numerical model
:math:`h` can be very high, it is suggested to approximate the model
response with a polynomial chaos expansion. However, for the sake of
computational simplicity, the latter is constructed considering
:math:`\mathit{independent}` components :math:`\{X^1,\dots,X^\inputDim\}`.
Thus the chaos basis is not orthogonal with respect to the correlated
inputs under consideration, and it is only used as a metamodel to
generate approximated evaluations of the model response and its summands
in :eq:`Model`.

.. math:: Y \simeq \hat{h} = \sum_{j=0}^{P-1} \alpha_j \Psi_j(x)

Then one may identify the component functions. For instance, for
:math:`u = \{1\}`:

.. math:: \model_1(X_1) = \sum_{\alpha | \alpha_1 \neq 0, \alpha_{i \neq 1} = 0} y_{\alpha} \Psi_{\alpha}(\inputRV)

where :math:`\alpha` is a set of degrees associated to the :math:`d`
univariate polynomial :math:`\psi_i^{\alpha_i}(X_i)`.

Then the model response :math:`Y` is evaluated using a sample
:math:`X=\{x_k, k=1,\dots,\sampleSize\}` of the correlated joint distribution.
Finally, the several indices are computed using the model response and
its component functions that have been identified on the polynomial
chaos.


.. topic:: API:

    - See :class:`~openturns.ANCOVA`


.. topic:: Examples:

    - See :doc:`/auto_reliability_sensitivity/sensitivity_analysis/plot_sensitivity_ancova`


.. topic:: References:

    - [caniou2012]_

