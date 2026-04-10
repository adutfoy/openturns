.. _enumeration_strategy:

Multivariate indices enumeration functions
------------------------------------------


Enumeration functions are bijections from :math:`\Nset` to :math:`\Nset^{\inputDim}`. We detail here some particular bijections:

- Linear enumeration function
- Hyperbolic enumeration function
- Anisotropic hyperbolic enumeration function
- Infinity norm enumeration function

All these enumeration functions are illustrated in :doc:`/auto_surrogate_modeling/polynomial_chaos/plot_enumeratefunction`.

A possible use is to build a multivariate basis as the tensorization of univariate basis: this is the
case for example in the functional chaos expansion setting (refer to  :ref:`functional_chaos` and
:ref:`enumeration_multivariate_basis`).

Let  :math:`\vect{\alpha}` be a multi-index is defined by:

.. math::

    \vect{\alpha} = (\alpha_1, \dots, \alpha_{\inputDim}) \in \Nset^{\inputDim}


An enumeration function :math:`\tau` is a bijection from :math:`\Nset` to :math:`\Nset^{\inputDim}`,
which creates a one-to-one mapping between an integer :math:`j` and a multi-index :math:`\vect{\alpha}`. The function :math:`\tau` is defined by:

.. math::

   \begin{array}{llcl}
         \tau \, : & \Nset & \longrightarrow & \Nset^{\inputDim} \\
         &  j & \longmapsto & \vect{\alpha}(j) = \{\alpha_1(j),\dots, \alpha_{\inputDim}(j)\}
    \end{array}


Let the *length* of any multi-index :math:`{\vect{\alpha}} \in {\Nset}^{\inputDim}` be defined by:

.. math::

    |{\vect{\alpha}}| = \sum_{i=1}^{\inputDim} \alpha_i


Linear enumeration function
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The linear enumeration function :math:`\tau : \Nset \longrightarrow \Nset^{\inputDim}` is a function:

.. math::
    :label: linearEnumFct

    \tau(j) = \vect{\alpha}(j) = (\alpha_1(j),\dots, \alpha_{\inputDim}(j))

for :math:`j \in \Nset` such that:

- first, :math:`\vect{\alpha}(0) = (0,\dots,0)`,
- for any :math:`k \in \Nset` and any :math:`j \in \{1, ..., k - 1\}`, we have:

.. math::

    \vect{\alpha}(j) < \vect{\alpha}(k)

which means that:

- either the length of :math:`\vect{\alpha}(j)` is strictly lower than :math:`\vect{\alpha}(k)`:

.. math::
    :label: cond_i

    \left|\vect{\alpha}(j)\right| < \left|\vect{\alpha}(k)\right|


- or the length of :math:`\vect{\alpha}(j)` equal to the length of :math:`\vect{\alpha}(k)`:.

.. math::
    :label: cond_ii

    \left|\vect{\alpha}(j)\right| = \left|\vect{\alpha}(k)\right|


and there exists :math:`m \in \{1,\dots,\inputDim\}` such that:

.. math::

    \begin{array}{ll}
    & \alpha_1(j) = \alpha_1(k) \\
    & \alpha_2(j) = \alpha_2(k) \\
    & \vdots \\
    & \alpha_{m - 1}(j) = \alpha_{m - 1}(k) \\
    & \alpha_m(j) < \alpha_m(k).
    \end{array}

Both conditions :eq:`cond_i` and :eq:`cond_ii` ensure that the mapping :math:`\tau`
implies a strict order on the set :math:`{\vect{\alpha}} \in {\Nset}^{\inputDim}`:
Condition :eq:`cond_i` states that the two multi-indices :math:`\vect{\alpha}_j` and :math:`\vect{\alpha}_k`
are not on the same strata; Condition :eq:`cond_ii` states that, if the two
multi-indices :math:`\vect{\alpha}_j` and :math:`\vect{\alpha}_k` are on the same strata,
then at least one of the component (denoted by :math:`m` in the definition) is different while the previous components are equal.

This numeration function is depicted in  :doc:`/auto_surrogate_modeling/polynomial_chaos/plot_enumeratefunction`.

Hyperbolic enumeration function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For any real number :math:`q` in :math:`(0,1]`, let
:math:`q`-*hyperbolic norm* (or :math:`q`-*norm* for short) of a
multi-index :math:`\vect{\alpha}` be defined by:

.. math::
   :label: hyperBolEnumFct

    \|\vect{\alpha}\|_{q} \, \, = \, \, \left(\sum_{i=1}^{\inputDim} \; \alpha_i^q \right)^{1/q}

The operator :math:`\|\cdot\|_q` is  a norm if anf only if :math:`q \geq 1` and is
a *pseudo-norm* if :math:`0 < q < 1` since it does not satisfy the triangular
inequality. However this abuse of language will be used in the
following. Note that the case :math:`q=1` corresponds to the
definition of the length of :math:`\vect{\alpha}`.

Let :math:`\lambda` be a real positive number. Let :math:`\cA_{\lambda}` be the set of
multi-indices with :math:`q`-norm not greater than :math:`\lambda` as
follows:

.. math::
    :label: eq_q_set

    \cA_{\lambda} \, \, = \, \, \{\vect{\alpha} \in \Nset^{\inputDim} \, : \, \|\vect{\alpha}\|_q \, \leq \lambda \}.

Moreover, let the *front* of :math:`\cA_{\lambda}` be defined by:

.. math::

   \partial \cA_{\lambda} \, \, = \, \, \left\{\vect{\alpha} \in \cA_{\lambda} \, : \, \exists \; i \; \in \; \{1,\dots,\inputDim\} \, , \, \, \vect{\alpha} \, + \, \vect{e_i} \, \notin \, \cA_{\lambda} \right\}

where :math:`\vect{e_i}` is a multi-index with a unit :math:`i`-entry
and zero :math:`k`-entries, :math:`k\neq i`.

We also define the set of *candidates* from the elements of :math:`\cA_\lambda`. The set of
the candidates is denoted by :math:`\cC_\lambda` and is defined by:

  .. math::

     \cC_\lambda\, \, = \, \, \left\{\vect{\alpha} \, + \, \vect{e_i} \, : \,
     \vect{\alpha} \in \partial \cA_{\lambda} \, , \,
     \vect{\alpha} + \, \vect{e_i} \notin  \cA_{\lambda} \, , \,
     1 \leq i \leq \inputDim, \right\}


We note that for all :math:`\lambda`, :math:`\cC_\lambda \neq \emptyset` because for any :math:`\lambda \in \Rset^+`,
there exists :math:`\vect{\alpha} \in \Nset^\inputDim` sich that  :math:`\|\vect{\alpha}\|_{q} > \lambda`.

The principle consists in exploring the space :math:`\Nset^{\inputDim}` through the
:math:`q`-norm of its elements. In this purpose, we define an appropriate
increasing sequence :math:`(\lambda_n)_{n \in \Nset}`  as follows:

.. math::

     \left\{
       \begin{array}{l}
         \lambda_0  =  0 \\
         \lambda_{n+1}  =  \min_{\vect{\alpha} \in \cC_{\lambda_n}}  \left\{ \|\vect{\alpha}\|_{q} \right\}
       \end{array}
     \right.

The sequence is well defined because by definition, all the elements of :math:`\cC_{\lambda_n}`
have a :math:`q`-*norm* strictly greater than :math:`\lambda_n`.

From the sequence :math:`(\lambda_n)_{n \in \Nset}`, we call :math:`(\cA_{\lambda_n})_{n \in \Nset}` the sequence of *cumulated strata*.
The sequence :math:`(\Delta_n)_{n \in \Nset}` of the strata is defined by:

.. math::

   \left\{
       \begin{array}{l}
         \Delta_0 =  \cA_{\lambda_{0}} \\
         \Delta_{n+1} =  \cA_{\lambda_{n+1}}  \setminus  \cA_{\lambda_n}
       \end{array}
    \right.
    
Note that we have :math:`\cA_{\lambda_{n+1}} \subset \cA_{\lambda_n} \cup \cC_n` and that
:math:`\Delta_n =  \left\{ \vect{\alpha} \in \Nset^\inputDim, \|\vect{\alpha}\|_{q} = \lambda_n \right\}`.

The sequence  of strata is depicted in  :doc:`/auto_surrogate_modeling/polynomial_chaos/plot_enumeratefunction`.

Anisotropic hyperbolic enumeration function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We consider enumeration functions based on an
*anisotropic* hyperbolic norm defined by:

.. math::
    :label: anisotropEnumFct

    \|\vect{\alpha}\|_{\vect{w},q} \, \, = \, \, \left(\sum_{i=1}^{\inputDim} \; w_i\; \alpha_i^q \right)^{1/q}

where the weights :math:`w_i` are real positive numbers. They enable to weight
some specific marginal indices.

In this setup, we consider both schemes outlined in the previous paragraph: it is only necessary to
replace the isotropic :math:`q`-norm in :eq:`eq_q_set` with the
:math:`(\vect{w},q)`-anisotropic one.

This enumerate function emphasizes multi-indices whose components are larger
when the associated weights are smaller.

The sequence  of strata is depicted in  :doc:`/auto_surrogate_modeling/polynomial_chaos/plot_enumeratefunction`.

Infinity norm enumeration function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We consider the enumeration function based on the infinite norm:

.. math::
     :label: infEnumFct

     \|\vect{\alpha}\|_{\infty} \, \, = \, \, \max_{1 \leq i \leq \inputDim} \; \alpha_i

This enumeration function is depicted in  :doc:`/auto_surrogate_modeling/polynomial_chaos/plot_enumeratefunction`.

Link between enumeration functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Note the following links.

If :math:`q=1`, the :math:`\|\vect{\alpha}\|_{1}` is the length of the multi_index. The hyperbolic enumeration
function becomes the linear enumeration function
for which the sequence becomes :math:`\lambda_n = n`.

If :math:`q \rightarrow 0`, then :math:`\|\vect{\alpha}\|_{q} \rightarrow \|\vect{\alpha}\|_{\infty}`
defined by :eq:`infEnumFct`. In that case, we have :math:`\lambda_n = n`.


.. topic:: API:

    - See :class:`~openturns.LinearEnumerateFunction`
    - See :class:`~openturns.HyperbolicAnisotropicEnumerateFunction`
    - See :class:`~openturns.NormInfEnumerateFunction`

.. topic:: Examples:

    - See :doc:`/auto_surrogate_modeling/polynomial_chaos/plot_functional_chaos`
    - See :doc:`/auto_surrogate_modeling/polynomial_chaos/plot_enumeratefunction`
    - See :doc:`/auto_surrogate_modeling/fields_surrogate_models/plot_fieldfunction_metamodel`

.. topic:: References:

    - [blatman2009]_
