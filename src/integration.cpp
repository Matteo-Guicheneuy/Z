// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <string>       // Strings                      //
// -------- Classes ----------------------------------- //
#include "process.h"            // Process definition   //
#include <gsl/gsl_monte_vegas.h>// GSL                  //
// -------- Functions --------------------------------- //
void DisplayXsec(const double& res, const double& err, const std::string& label);
double Tot(double*,size_t,void*);


// ************************************************************************* //
//  Helper Function for VEGAS                                                //
// ************************************************************************* //
void PerformIntegration(double& res, double& err, double& chi, Process* proc, size_t dim)
{
  // Initialize random number generator
  gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(r, time(NULL));

  // Define the integrand
  gsl_monte_function I = { &Tot, dim, proc };

  // Integration bounds
  double xmin[dim], xmax[dim];
  for (size_t i = 0; i < dim; ++i) { xmin[i] = 0.0; xmax[i] = 1.0; }

  // Display initialization
  DisplayXsec(res, err, "-init");

  // Allocate and initialize VEGAS state
  gsl_monte_vegas_state* s = gsl_monte_vegas_alloc(dim);

  // Warm-up
  gsl_monte_vegas_integrate(&I, xmin, xmax, dim, 10000 / 50, r, s, &res, &err);
  chi = s->chisq;
  DisplayXsec(res, err, "Warm-up");

  // Precision target and refinement loop
  double precision = 1e9;
  int counter = 1;
  size_t calls = 10000;
  while(precision>1e-2 && counter<=10)
  {
    std::string label = "Refine-" + std::to_string(counter);
    gsl_monte_vegas_integrate(&I, xmin, xmax, dim, calls, r, s, &res, &err);
    precision = std::abs(err / res);
    DisplayXsec(res, err, label);
    counter++;
    calls *= 5;
  }

  // Final display
  chi = s->chisq;
  DisplayXsec(res, err, "final");

  // Cleanup
  gsl_monte_vegas_free(s);
  gsl_rng_free(r);
}
