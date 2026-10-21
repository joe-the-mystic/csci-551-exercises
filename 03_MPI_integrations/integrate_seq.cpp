//////////////////////////////////////////////////////////////////////////////
// Project: Integrate
// Author: Josiah Lawrence
// Date: 28 September 2026
// Description: This project numerically integrates a function over a given
// interval. It was developed at CSU, Chico. This version uses a constant step
// size.
//////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include <iomanip>
#include <iostream>

#define STEP_SIZE 0.0001
#define OUTPUT_PRECISION 8

// Train Problem
// #define INTERVAL_FROM_ZERO (1800)
// #define TSCALE (1800 / (2 * M_PI))
// #define ASCALE 0.236589076381454

// Basic interval from 0 to pi
#define INTERVAL_FROM_ZERO (M_PI)
#define TSCALE 1
#define ASCALE 1

// Interval from 0 to 10
// #define INTERVAL_FROM_ZERO 10
// #define TSCALE 1
// #define ASCALE 1

using namespace std;

/* Defines the function to be integrated */
double
integrand(double x)
{
  // return 10.0;
  return (sin(x / TSCALE) * ASCALE); // acceleration
}

/* Defines the antiderivative for the purpose of testing the numerical 
 * integration. */
double
antiderivative(double x)
{
  // return 10 * x;
  return (1.0 - cos(x / TSCALE)) * ASCALE * TSCALE; // velocity
}

/* Computes the middle riemann sum over the given interval */
double
riemann_sum(double start, double end)
{
  int nstep = lround((end - start) / STEP_SIZE);
  double sum = 0.0;
  double dx = (end - start) / nstep;
  
  cout << "Step size = " << dx << " for " << nstep
  << " steps over interval " << start << " to " << end << "." << endl;
  
  for (int idx = 0; idx < nstep; idx++) {
    double x = 0.0, fx = 0.0;
    x = start + (idx * dx) + (dx / 2);
    fx = integrand(x);

    // Add the rectangle at the middle of each subinterval.
    sum += fx * dx;
  }

  return sum;
}

/* Computes the trapezoidal sum over the given interval */
double
trapezoidal_sum(double start, double end)
{
  int nstep = lround((end - start) / STEP_SIZE);
  double sum = 0.0;
  double dx = (end - start) / nstep;
  
  cout << "Step size = " << dx << " for " << nstep
  << " steps over interval " << start << " to " << end << "." << endl;
  
  for (int idx = 0; idx < nstep; idx++) {
    double x = 0.0;
    x = start + (idx * dx);

    // Add the area of the trapezoid on the current subinterval.
    sum += ( ((integrand(x) + integrand(x+dx)) / 2) * dx );
  }

  return sum;
}

int
main(int argc, char* argv[])
{
  double start = 0.0;
  double end = INTERVAL_FROM_ZERO;
  double result = 0.0;
  cout << std::fixed << std::setprecision(OUTPUT_PRECISION);

  /* Use a Riemann sum to estimate the definite integral using the area under
   * the curve.*/
  result = riemann_sum(start, end);

  cout << "Riemann sum: " << result << endl;

  /* Use the trapezoid rule to estimate the definite integral using the area under
   * the curve.*/
  result = trapezoidal_sum(start, end);

  cout << "Trapezoidal sum: " << result << endl;

  /* Check the numerical integration using the 2nd fundamental theorem of
   * Calculus: "The definite integral of a function over a closed interval where
   * the function is Riemann integrable is equal to the difference between the
   * values of an antiderivative evaluated at the endpoints of the interval.""
   * -- Wikipedia
   */

  cout << "Verify by antiderivative: " << antiderivative(end) - antiderivative(start) << endl;
}
