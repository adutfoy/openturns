Create an optimisation problem with inequality constraints
==========================================================

Create the following optimization problem:


.. math::

    \min 2(x_1-2)^2+3(x_2-2)^2\\
    \vect{x}\in [-5, 5]^2 \\
    x_1^2+x_2^2 \leq 4

.. literalinclude:: t_OptimizationInequality.py
