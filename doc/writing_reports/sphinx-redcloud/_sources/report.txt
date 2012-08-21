.. Automatically generated reST file from Doconce source
   (http://code.google.com/p/doconce/)


Experiments with Schemes for Exponential Decay
==============================================

:Author: Hans Petter Langtangen (hpl at simula.no)

:Date: August 20, 2012

*Summary.* This report investigates the accuracy of three finite difference
schemes for the ordinary differential equation :math:`u'=-au` with the
aid of numerical experiments. Numerical artifacts are in particular
demonstrated.






.. Section with multi-line equation.


.. _math:problem:


Mathematical problem
====================


.. index:: model problem

.. index:: exponential decay


We address the initial-value problem


.. math::
        
        u'(t) &= -au(t), \quad t \in (0,T], \\
        u(0)  &= I,                         
        

where :math:`a`, :math:`I`, and :math:`T` are prescribed parameters, and :math:`u(t)` is
the unknown function to be estimated. This mathematical model
is relevant for physical phenomena featuring exponential decay
in time.

.. Section with single-line equation and a bullet list.


.. _numerical:problem:

Numerical solution method
=========================


.. index:: mesh in time

.. index:: :math:`\theta`-rule

.. index:: numerical scheme


.. index:: finite difference scheme


We introduce a mesh in time with points :math:`0= t_0< t_1 \cdots < t_N=T`.
For simplicity, we assume constant spacing :math:`\Delta t` between the
mesh points: :math:`\Delta t = t_{n}-t_{n-1}`, :math:`n=1,\ldots,N`. Let
:math:`u^n` be the numerical approximation to the exact solution at :math:`t_n`.

The :math:`\theta`-rule is used to solve (:ref:`ode`) numerically:


.. math::
        
        u^{n+1} = \frac{1 - (1-\theta) a\Delta t}{1 + \theta a\Delta t}u^n,
        

for :math:`n=0,1,\ldots,N-1`. This scheme corresponds to

  * The Forward Euler scheme when :math:`\theta=0`

  * The Backward Euler scheme when :math:`\theta=1`

  * The Crank-Nicolson scheme when :math:`\theta=1/2`
.. Section with computer code taken from a part of

.. a file. The fromto: f@t syntax copies from the

.. regular expression f up to the line, but not

.. including, the regular expression t.


Implementation
==============

The numerical method is implemented in a Python function:


.. code-block:: python

        def theta_rule(I, a, T, dt, theta):
            """Solve u'=-a*u, u(0)=I, for t in (0,T] with steps of dt."""
            N = int(round(T/float(dt)))  # no of intervals
            u = zeros(N+1)
            t = linspace(0, T, N+1)
        
            u[0] = I
            for n in range(0, N):
                u[n+1] = (1 - (1-theta)*a*dt)/(1 + theta*dt*a)*u[n]
            return u, t


.. Section with figures.


Numerical experiments
=====================


.. index:: numerical experiments


We define a set of numerical experiments where :math:`I`, :math:`a`, and :math:`T` are
fixed, while :math:`\Delta t` and :math:`\theta` are varied.
In particular, :math:`I=1`, :math:`a=2`, :math:`\Delta t = 1.25, 0.75, 0.5, 0.1`.



.. Subsection with inline figure (figure without caption).


The Backward Euler method
-------------------------


.. index:: BE



.. figure:: BE.png
   :width: 800





.. Subsection with inline figure (figure without caption).


The Crank-Nicolson method
-------------------------


.. index:: CN



.. figure:: CN.png
   :width: 800





.. Subsection with inline figure (figure without caption).


The Forward Euler method
------------------------


.. index:: FE



.. figure:: FE.png
   :width: 800





Error vs :math:`\Delta t`
-------------------------

.. Exemplify referring to a figure with label and caption.



.. index:: error vs time step


How :math:`E` varies with :math:`\Delta t` for :math:`\theta=0,0.5,1`
is shown in Figure :ref:`fig:E`.


.. _fig:E:

.. figure:: error.png
   :width: 400

   *Error versus time step*
