Examples
========

Example : Resistance Stress model 
---------------------------------

1. Problem statement
````````````````````

The example is a resistance minus stress performance function used with normal
random variable :

.. math::

    g(\uX) = R - S

with

.. math::
    
    \begin{aligned}
    R \sim \mathcal N (7, 1) \\
    S \sim \mathcal N (2, 1)
    \end{aligned}

The failure situation is considered when the performance function :math:`g` is
being negative. The computing probability reads :

.. math::
    
    \Prob{ g(\uX) \leq q}

The goal is to find :math:`q` such as the probability is equal to :math:`2.10^{-4}`.


2. Resolution
`````````````

.. literalinclude:: t_SubsetInverseSampling_R-S.py

3. Result
`````````

.. literalinclude:: t_SubsetInverseSampling_R-S.expout

