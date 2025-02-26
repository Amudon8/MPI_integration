# MPI_integration

This is a simple program for doing numerical integration using trapezoidal rule. The issue faced in doing the integration is the edge of the node in MPI. This has been resolved in the program by simply implementing the following equation,
$$
I = \int_{a}^{b} f(x) dx \approx \frac{h}{2} \left[ f(a) + 2*\sum_{i = 1}^{n} f(x_{i}) + f(b) \right]
$$

You can defined your own profile which is assigned to the variable 'local_a'. Here I am using only sine.
