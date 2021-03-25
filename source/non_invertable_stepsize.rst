Non-invertable Step Size
==================


F-based
-------------

Assume the current optimziation variable is :math:`x`, the searching direction is :math:`\Delta x`.

Denote

.. math::

   F(\alpha) = F(x + \alpha \Delta x)

We would like to solve the following equation:

.. math::

    \det F(\alpha) = \epsilon \det F(0)

The equation can be converted to the following form:

.. math::

    \det (I + \alpha A) = \epsilon


FEM
^^^^^^^^^^

.. math::
    \begin{split}
        F(\alpha) &= (T + \alpha \Delta T) B^{-1} \\
        F(\alpha)F(0)^{-1} &= (T + \alpha \Delta T) T^{-1} = I + \alpha \Delta T T^{-1}
    \end{split}

MPM
^^^^^^^^^^

.. math::

    \begin{split}
    F(\alpha) &= (\sum_i (x_i + \alpha \Delta x_i ) (\nabla \omega_{ip}^n)^T)F^n\\
            &= F(0) + \alpha (\sum_i \Delta x_i (\nabla \omega_{ip}^n)^T)F^n\\
    F(\alpha)F(0)^{-1} &= I + \alpha (\sum_i \Delta x_i (\nabla \omega_{ip}^n)^T)F^nF(0)^{-1}
    \end{split}


Use the following code to solve the equation:

.. code-block:: C++

    T a, b, c, d;
    if constexpr (dim == 2) {
        a = 0;
        b = A.determinant();
    }
    else { // dim == 3
        a = A.determinant();
        b = A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0) + A(0, 0) * A(2, 2) - A(0, 2) * A(2, 0) + A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1);
    }
    c = A.diagonal().sum();
    d = 1 - eps;
    alpha = Math::get_smallest_positive_real_cubic_root(a, b, c, d);