# polynomial_solver_degree_smaller_than_4
This code enables to find the roots of low degree polynomials.

This is an implementation of exact root expressions from Cardano's formulas (degree 3 poly.) and Euler formulas (degree 4 poly.). Degree 2, 1 and 0 are also dealt with.

The main difference with numpy.roots is that it computes the roots with floats. It may be a very (extremely ?) small difference.

The problem I encountered (and where this all code comes from incidentally) is that numpy.roots seems to return 0.0 results for very small coefficients whereas, using floats, you do get usually non zero results. Indeed the result are approximate as always.

