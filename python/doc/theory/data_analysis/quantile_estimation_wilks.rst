.. _quantile_estimation_wilks:

Estimation of a quantile upper bound by Wilks' method
-----------------------------------------------------

We consider a random variable :math:`X` of dimension 1 and its unknown quantile :math:`x_{\alpha}` of order :math:`\alpha`.
We seek to evaluate an upper bound of :math:`x_{\alpha}` with a confidence greater than :math:`\beta`, using a given order statistics.

We denote by :math:`(X_1, \dots, X_\sampleSize)` some independent and identically distributed variables according to :math:`X`.

We denote by :math:`(X_1, \dots, X_\sampleSize)` some independent and identically distributed variables according to :math:`X`.

Let :math:`X_{(k)}` be the :math:`k` -th order statistics of :math:`(X_1, \dots, X_\sampleSize)` which means that
:math:`X_{(k)}` is the :math:`k` -th maximum of :math:`(X_1, \dots, X_\sampleSize)` for :math:`1 \leq k \leq \sampleSize`. For
example, :math:`X_{(1)} = \min (X_1, \dots, X_\sampleSize)` is the minimum
and :math:`X_{(\sampleSize)} = \max (X_1, \dots, X_\sampleSize)` is the maximum. We have:

.. math::

    X_{(1)} \leq X_{(2)} \leq \dots \leq X_{(\sampleSize)}


Estimation of the order statistics that will provide an upper bound to :math:`x_{\alpha}`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We assume that we have a sample of size :math:`\sampleSize` of :math:`X`, denoted by :math:`(x_1, \dots, x_\sampleSize)`.

Given :math:`\alpha`, :math:`\beta`, and :math:`\sampleSize`, we seek for the smallest :math:`k` such that:

.. math::
    :label: EqOrderStat

    \Prob{X_{(k)} \geq x_{\alpha}} \geq \beta

The order statistics :math:`X_{(k)}` follows the distribution which probability density function and cumulated distribution function
are defined by:

.. math::
    :label: DistOrderStat

    F_{X_{(k)}}(x) & = \sum_{i=k}^{\sampleSize} \binom{\sampleSize}{i}\left[F(x)\right]^i \left[1-F(x)
    \right]^{\sampleSize-i} \\
    p_{X_{(k)}}(x) & = (\sampleSize-k+1)\binom{\sampleSize}{k-1}\left[F(x)\right]^{k-1} \left[1-F(x)
    \right]^{\sampleSize-k} p(x)

We notice that :math:`F_{X_{(k)}}(x) = \overline{F}_{(\sampleSize,F(x))}(k-1)` where :math:`F_{(\sampleSize,F(x))}` is the cumulated
distribution function of the Binomial distribution :math:`\cB(\sampleSize,F(x))`.

Then, we have:

.. math::

    F_{X_{(k)}}(x_{\alpha}) = \sum_{i=k}^{\sampleSize} \binom{\sampleSize}{i} \alpha^i (1-\alpha)^{\sampleSize-i}
    = \overline{F}_{(\sampleSize,\alpha)}(k-1)

and relation :eq:`EqOrderStat` can be written as:

.. math::
    :label: EqOrderStat2

    1-F_{X_{(k)}}(x_{\alpha})\geq \beta \\
    F_{\sampleSize, \alpha}(k-1)\geq \beta

The smallest order :math:`k_{sol}` is defined by:

.. math::

    k_{sol} & = \min \{ k \in [1, n] \, | \, F_{\sampleSize, \alpha}(k-1)\geq \beta \}\\
            & = 1 +  \min \{ k \in [1, n] \, | \, F_{\sampleSize, \alpha}(k)\geq \beta \}

An upper bound of  :math:`x_{\alpha}` is estimated by the value of :math:`X_{(k_{sol})}`  on the sample
:math:`(x_1, \dots, x_\sampleSize)`.

Estimation of the smallest sample size that ensures that a given order statistics will provide an upper bound
to :math:`x_{\alpha}`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~

Given :math:`\alpha`, :math:`\beta`, and :math:`k`, we seek for the smallest sample size :math:`\sampleSize`
that verifies :eq:`EqOrderStat`. This time, :eq:`EqOrderStat2` is solved wrt :math:`\sampleSize`.

Once the smallest size :math:`\sampleSize`  has been estimated, a sample of size :math:`\sampleSize` can be
generated from
:math:`X` and an upper bound of :math:`x_{\alpha}` is estimated by the value of the value of
:math:`X_{(\sampleSize-i)}` on the sample
:math:`(x_1, \dots, x_\sampleSize)`.


.. topic:: API:

    - See :class:`~openturns.Wilks`

.. topic:: Examples:

    - See :doc:`/auto_data_analysis/manage_data_and_samples/plot_quantile_estimation_wilks`

.. topic:: References:

    - Wilks, S.S. (1962). "Mathematical Statistics", New York-London
    - Robert C.P., Casella G. (2004). Monte-Carlo Statistical Methods, Springer, ISBN 0-387-21239-6, 2nd ed.
    - Rubinstein R.Y. (1981). Simulation and The Monte-Carlo methods, John Wiley & Sons
