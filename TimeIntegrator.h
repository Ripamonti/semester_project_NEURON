/*
 * © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
 * Laboratory ANMC, 2016
 * See the LICENSE.TXT file for more details.
 * 
 * Authors
 * The original Fortran codes were developed by Assyr Abdulle (for ROCK2) and 
 * Sommeijer, Shampine, Verwer (for RKC). Translation to C++ by Giacomo Rosilho
 * de Souza.
 * 
 * The ROCK2 code is described in 
 * 
 * Abdulle, Assyr, and Alexei A. Medovikov. 
 * "Second order Chebyshev methods based on orthogonal polynomials." 
 * Numerische Mathematik 90.1 (2001): 1-18.	
 * 
 * and the RKC code in
 * 
 * Sommeijer, B. P., L. F. Shampine, and J. G. Verwer. 
 * "RKC: An explicit solver for parabolic PDEs." 
 * Journal of Computational and Applied Mathematics 88.2 (1998): 315-326.
 * 
 * Please cite these articles in any publication describing research
 * performed using the software.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef TIMEINTEGRATOR_H
#define	TIMEINTEGRATOR_H

#include <limits>
#include <vector>
#include <iostream>
#include <cmath>
#include "Ode.h"
#include "operation.h" 

using namespace std;

template<class ODE> class TimeIntegrator {
  public:
    TimeIntegrator(bool onestep=true, bool verb=true, bool dtadap=true, double atol=1e-2, double rtol=1e-2, bool intrho=false, bool scalartol=true);
    void advance(ODE* ode, double& h, int& idid);	// Compute next step
    int check_correctness(double& h); 		// Check if the input parameters are ok
    void print_info();				// Starting infos
    void print_statistics();			// Statistics at the end

  protected:
    virtual void rtstep(ODE* ode, const double t, const double& h, vector<double>& y, vector<double>& yn)=0; // Do the stages
    void accepted_step(double& t, double& h, vector<double>& y);
    void rejected_step(double& t, double& h, vector<double>& y);
    void compute_hnew(double& h);		
    virtual void update_n_stages(double& h)=0; 			// Given rho computes the number of stages
    void update_rho(ODE* ode, int& idid);				// Updates the spectral radius
    void rho(double& eigmax, ODE* ode, int& idid);	// Power method approximating rho
  public:
  /*----------------------------------------- 
    Integration settings 
  -----------------------------------------*/
    bool one_step;		// If true just one step is executed, to do the next one just recall advance with the same arguments. If false we integrate until end.
    bool dt_adaptivity;	// If false the time step is fixed
    bool internal_rho;	// If true the ODE's spectral radius is computed by the internal power like method. Else a method rho is provided by the ODE class.
    bool scalar_tol;	// Tells if the tolerances are scalars or if every solution's component has its own
    double r_tol;		// Scalar relative tolerance
    double a_tol;		// Scalar absolute tolerance
    bool verbose;		// If true prints info about the timestep

  /*----------------------------------------- 
    Integration statics 
  -----------------------------------------*/
    int max_rho;		// Maximal spectral radius computed
    int min_rho;		// Minimal spectral radius computed
    int max_s;		// Maximal number of stages used
    double dt_max;		// Maximal time step used
    double dt_min;		// Minimal time step used
    int n_f_eval_rho;	// Number of right end side evaluations for the spectral radius computation in the power method
    int n_f_eval;		// Number of right end side evaluations for time integration
    int n_steps;		// Number of time steps
    int acc_steps;		// Number of accepted steps
    int rej_steps;		// Number of rejected steps
 
  protected:
  /*----------------------------------------- 
    Step variables
  -----------------------------------------*/ 
    double eigmax;		// Spectral radius of current time step
    double uround;		// Minimal time step allowed
    double err;		// Local error estimation
    double errp;		// Local error estimated at the previuos step
    double fac;		// Factor used in the choice of the new time step
    double facp;		// Factor used in the previous step
    double facmax;		// Maximal allowed factor
    double hp;		// Previous time step size
    double hnew;		// Next time step size
    double told;		// Starting time of previous time step. Used to detect multiple rejections
    int s;			// Number od stages for RKC. Number of stages minus 2 for ROCK2
    int sold;		// Previous number of stages
    int nrej;		// Number of consecutive rejections
    int nrho;		// Number of steps passed after last spectral radius estimation
    bool last;		// Is true if the current step is the last one
    bool reject;		// Is true if last step has been rejected
    bool fn_uptodate;	// Tells if the right hand side has already been computed for the current step

};

#endif	/* TIMEINTEGRATOR_H */
