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

template<class ODE> class TimeIntegrator {
  public:
    // Constructor
    TimeIntegrator(bool onestep=true, bool verb=true, bool dtadap=true, double atol=1e-2, double rtol=1e-2, bool intrho=false, bool scalartol=true);
    // Method to advance with the next step
    void advance(ODE* ode, double& h, int& idid);
    // Method to check the correctness of problem parameters
    int check_correctness(double& h);
    // Print info of the methos
    void print_info();
    // Print info of the simulation at the end
    void print_statistics();		
  protected:
    // Compute a step
    virtual void rtstep(ODE* ode, const double t, const double& h, std::vector<double>& y, std::vector<double>& yn)=0;
    // Method that manages a succesful step
    void accepted_step(double& t, double& h, std::vector<double>& y);
    // Method that reset the simulation at the previous step if the actual one is rejected
    void rejected_step(double& t, double& h, std::vector<double>& y);
    // Compute the next time step (need adaptive option)
    void compute_hnew(double& h);		
    // Update the coefficients
    virtual void update_n_stages(double& h)=0;
    // Method to update internally the spectral radius
    //  N.B. The possibility to compute the spectral radius inside the ODE class is not yet 
    //  implemented
    void update_rho(ODE* ode, int& idid);
    // Compute the spectral radius
    void rho(double& eigmax, ODE* ode, int& idid);
  public:
  /*----------------------------------------- 
    Integration settings 
  -----------------------------------------*/
    bool one_step;		// If true just one step is executed, to do the next one just recall advance with the same arguments. If false we integrate until end.
    bool dt_adaptivity;	        // If false the time step is fixed
    bool internal_rho;	        // If true the ODE's spectral radius is computed by the internal power like method. Else a method rho is provided by the ODE class (TODO: implement this second possibility)
    bool scalar_tol;	        // Tells if the tolerances are scalars or if every solution's component has its own
    double r_tol;		// Scalar relative tolerance
    double a_tol;		// Scalar absolute tolerance
    bool verbose;		// If true prints info about the timestep
  /*----------------------------------------- 
    Integration statics 
  -----------------------------------------*/
    double max_rho;		// Maximal spectral radius computed
    double min_rho;		// Minimal spectral radius computed
    int max_s;		        // Maximal number of stages used
    double dt_max;		// Maximal time step used
    double dt_min;		// Minimal time step used
    int n_f_eval_rho;	        // Number of right end side evaluations for the spectral radius computation in the power method
    int n_f_eval;		// Number of right end side evaluations for time integration
    int n_steps;		// Number of time steps
    int acc_steps;		// Number of accepted steps
    int rej_steps;		// Number of rejected steps
  protected:
  /*----------------------------------------- 
    Step variables
  -----------------------------------------*/
protected: 
    double eigmax;		// Spectral radius of current time step
    double uround;		// Minimal time step allowed
public:
    double err;		        // Local error estimation
protected:
    double errp;		// Local error estimated at the previuos step
    double fac;		        // Factor used in the choice of the new time step
    double facp;		// Factor used in the previous step
    double facmax;		// Maximal allowed factor
    double hp;		        // Previous time step size
    double hnew;		// Next time step size
    double told;		// Starting time of previous time step. Used to detect multiple rejections
    int s;			// Number od stages for RKC. Number of stages minus 2 for ROCK2
    int sold;		        // Previous number of stages
    int nrej;		        // Number of consecutive rejections
    int nrho;		        // Number of steps passed after last spectral radius estimation
    bool last;		        // Is true if the current step is the last one
    bool reject;		// Is true if last step has been rejected
    bool fn_uptodate;	        // Tells if the right hand side has already been computed for the current step

};

template <class ODE>
TimeIntegrator<ODE>::TimeIntegrator(bool onestep, bool verb, bool dtadap, double atol, double rtol, bool intrho, bool scalartol)
{
  uround=std::numeric_limits<double>::epsilon();
  errp=0.0;
  nrho=1;
  facmax=5.0;
  s=0;
  sold=0;
  nrej=0;
  
  one_step=onestep;
  verbose=verb;
  dt_adaptivity=dtadap;
  internal_rho=intrho;
  scalar_tol=scalartol;
  r_tol=rtol;
  a_tol=atol;

  max_rho=std::numeric_limits<int>::min();
  min_rho=std::numeric_limits<int>::max();
  max_s=0;
  n_f_eval_rho=0;
  n_f_eval=0;
  n_steps=0;
  acc_steps=0;
  rej_steps=0;
  dt_max=0.;
  dt_min=std::numeric_limits<double>::max();
}

template <class ODE>
void TimeIntegrator<ODE>::advance(ODE* ode, double& h, int& idid)
{
  // Define a pointer to the actual time in ODE class: if it is modified here it is modified as well 
  // as in the class
  double& t=ode->time;
  double tend=ode->tend;
  // Define pointers to the vectors in the ode class for the solutions at the previous and actual time steps
  std::vector<double>& y = ode->ynpu;	
  std::vector<double>& yn = ode->yn; 

  told=t-1.;	 // Different from t
  last=false;    // A priori is not the last step
  reject=false;	 // If we enter advance then the previous step has been accepted
  idid=1;

  // If it is the first step then we:
  // i)  recover the initial value from ode
  // ii) compute the spectral radius of the system
  if (n_steps==0)
  {
    yn=ode->get_un();	  // Takes the initial value from ODE
    fn_uptodate=false;	  // We have a new yn, so fn=f(yn) has to be updated
    update_rho(ode,idid); // Computes spectral radius    
  }
  // else the last solution is already contained in yn and we already have an approximation of the spectral radius

  /*----------------------------------------- 
    Integration loop 
  -----------------------------------------*/
  for(;;)
  {
    n_steps++;	// New time step
    // If we are very close to end or t+h>tend then adjust h
    if (1.1*std::abs(h)>=std::abs(tend-t))	// h can be negative (when t>tend)
    {
      h=tend-t;
      last=true; // It will be the last step (if accepted)
    }
    if (h<10.0*uround)	// Too small h
    {
      std::cout << "Tolerances are too small." << std::endl;
      idid=-2;	// idid=-2 means that there's a problem
      return; 	// Exit
    }

    // Do some statistics
    dt_max=std::max(dt_max,h);
    dt_min=std::min(dt_min,h);
     
    // Compute spectral radius every 25 steps if it isn't constant
    if (!ode->cte_rho && nrho==0) // nrho==0 every 25 steps
      update_rho(ode,idid);

    // The number of stages chosen depending on h and spectral radius
    update_n_stages(h);

    // Computation of an integration step.
    // This function is implemented in derived classes
    // y will contain the new value, ode provides the right hand side function
    rtstep(ode,t,h,y,yn);

    // Store the last used number of stages
    sold=s;

    if (dt_adaptivity)  // If time step adaptivity enabled
      compute_hnew(h);  // Computed new h
    else                // We are not using time step adaptivity
      hnew=h;           // keeps same h

    // Accepted step or without time step adaptivity
    if (err<1.0 || !dt_adaptivity)
    {
      accepted_step(t,h,y);
      
      yn.swap(y);        // Put the solution in the right space  
      fn_uptodate=false; // fn not up to date anymore

      if (last||one_step)  // We have finished
      {
        // remember that idid=1
        ode->set_un(yn); // Pass value to ODE, for outputs
        if (one_step)
          idid=2;
        break;
      }
    } else if(one_step) {
      ode->set_un(yn);	//pass value to ode, for outputs
      idid=2;           //says to outside world that we didn't reach the end
      break;
    } else
    rejected_step(t,h,y);
    }
    return;
}

template <class ODE>
void TimeIntegrator<ODE>::update_rho(ODE* ode, int& idid)
{
  // Computed externally by ODE::rho
  //  N.B. Actually it does nothing.... Chose always the second option
  if (!internal_rho)
    ode->rho(eigmax);
  // Computed internally by TimeIntegrator::rho
  else
    this->rho(eigmax,ode,idid);
  
  // Recover statistics
  if ((int)eigmax>max_rho) max_rho=(int)eigmax;
  if ((int)eigmax<min_rho) min_rho=(int)eigmax;

  // Print info
  // TODO: Since this method is used to solve only  part of the system the following output 
  // gives a not so good layout. We keep the next 3 lines of code in case this function will be used 
  // separately
  /*
  std::cout<<"-------- Spectral radius estimations --------"<<std::endl;
  std::cout<<"Spectral radius = "<<(int)eigmax<<std::endl;
  std::cout<<"---------------------------------------------\n"<<std::endl;   
  */
}


template <class ODE>
void TimeIntegrator<ODE>::compute_hnew(double& h)
{
  /*-------------------------------------------
    Estimation of the optimal time step
  -------------------------------------------*/
  // In the first time step we do not have a previous error estimation so errp=0

  // We use this estimation if previous error not zero and previous step has been accepted
  if (!reject && errp)
  {
    facp=1.0/err;
    fac=facp*(h/hp);
    fac=errp*fac*fac;
    fac=std::min(facp,fac);
    fac=sqrt(fac);
  } 
  else // In first step we go here, or if the previous step has been rejected
    fac=sqrt(1.0/err);
    fac=std::min(facmax,std::max(0.1,0.8*fac));  
    hnew=h*fac;
}

template <class ODE>
void TimeIntegrator<ODE>::accepted_step(double& t, double& h, std::vector<double>& y)
{
  // Called when the step is accepted. Does output, sets the new time step size and update some variables
  // Define delta character and red color
  std::string delta=u8"\u0394";
  std::string bcol="\033[31;1m";
  std::string ecol="\033[0m";
  if (dt_adaptivity && verbose)
    std::cout<<"Step t = "<<t<<", "<<delta<<"t = "<<h
    <<", s = "<<s<<". Acceptd with ||e_n||_L2 = "<<err<<std::endl;
    // TODO: infinity norm need to be implemented
    // <<" and ||y_n||_Linf = "<<y 
  else if (verbose)
    std::cout<<"Step t = "<<t<<", "<<delta<<"t = "<<h
    <<", s = "<<s<<","<<(err>1? bcol:"")<<" ||e_n||_L2 = "<<err<<(err>1? ecol:"")<<std::endl;
  
  acc_steps++;
  facmax=2.0;  // h can grow 2 times
  t=t+h; 
  errp=err;
  if (reject)  // The previous time step has been rejected so a smaller h is chosen
  {
    hnew=h>0.0 ? std::min(hnew,h):std::max(hnew,h);
    reject=false;  // This step has been accepted
    nrej=0;        // Set the consecutive rejections to 0
  }
  hp=h;			// Previous h
  h=hnew;		// Next h
  nrho=nrho+1;		// Consecutive steps without computing the spectral radius
  nrho=(nrho+1)%25;	// Set to 0 every 25 steps
}

template <class ODE>
void TimeIntegrator<ODE>::rejected_step(double& t, double& h, std::vector<double>& y)
{
  // Called when the step is rejected. Does output, choses the new time step and updates some variables
  // Define delta character and red color
  std::string delta=u8"\u0394";
  std::string bcol="\033[31;1m";
  std::string ecol="\033[0m";
  if (verbose)
   std::cout<<"Step t = "<<t<<", "<<delta<<"t = "<<h
   <<", s = "<<s<<"."<<bcol<<" Refected with ||e_n||_L2 = "<<err<<ecol<<std::endl;

  rej_steps++;
  reject=true;
  last=false;	// The step is rejected, it can't be the last one
  h=0.8*hnew;
  facmax=2.0;   // Next step h can at most double
  if (n_steps==0)    // Here we are in the first step, starting time was too large
    h=0.1*h;
  if(told==t)	// The previous step was also rejected
  {
    nrej=nrej+1;   // Consecutive rejections
    if (nrej==10)   // after 10 consecutive rejections ..
      h=h*1e-3;
  }
  else
    told=t;	// First rejection

  // The spectral radius is recomputed after a step failure
  if (nrho)
    nrho=0;	// Spectral radius will be computed
  else
    nrho=1;	// Was already computed
}

template <class ODE>
int TimeIntegrator<ODE>::check_correctness(double& h)
{
  // Some parameters checking before starting integration
  // Test the initial step size and tolerances
  if (h<10.0*uround)
  {
    std::cout<<"Initial step-size is too small."<<std::endl;
    return 0;
  }
  if (scalar_tol)	// Scalar tolerances
  {
    if (a_tol<=0.0 || r_tol<=10.0*uround)
    {
      std::cout<<"Tolerances are too small."<<std::endl;
      return 0;
    }
  }
  else
  {
    std::cerr<<"NON SCALAR TOLERANCES NOT IMPLEMENTED YET"<<std::endl;
    return 0;
  }
  return 1;
}

template <class ODE>
void TimeIntegrator<ODE>::rho(double& eigmax, ODE* ode, int& idid)
{
  double eigmaxo,sqrtu,znor,ynor,quot,dzyn,dfzfn;
  const int maxiter=50;
  const double safe=1.2;
  const double tol=0.01;
  double t=ode->time;

  sqrtu=sqrt(uround);

// ------ The initial vectors for the power method are yn --------
//       and yn+c*f(v_n), where vn=f(yn) a perturbation of yn 
//       (if iwork(6)=0) or a perturbation of the last computed
//       eigenvector (if iwork(6).neq.0).

    std::vector<double>& yn = ode->yn;
    std::vector<double>& fn = ode->fn;
    std::vector<double>& dz = ode->dz1;
    std::vector<double>& z  = ode->tmp1;

    if (n_steps==0)
    {
        ode->rhs(t,yn,fn);
        n_f_eval++;

        ode->set_Dirichlet(t,fn);
        ode->rhs(t,fn,z);
        ode->set_Dirichlet(t,z,true);
        n_f_eval_rho++;
        
        ode->set_Dirichlet(t,fn,true);
        fn_uptodate=true;
    }
    else
        z=dz;

// ------ Perturbation.--------
    ode->set_Dirichlet(t,yn,true);
    ynor= l2_norm(yn);
    ode->set_Dirichlet(t,yn);
    znor= l2_norm(z);
    if(ynor!=0.0 && znor!=0.0)
    {
        dzyn=ynor*sqrtu;
        quot=dzyn/znor;
        scale(z,quot);
        add(z,yn);
    }
    else if(ynor!=0.0)
    {
        dzyn=ynor*sqrtu;
        z=yn;
        scale(z,1.+sqrtu);
        ode->set_Dirichlet(t,z);
    }
    else if(znor!=0.0)
    {
        dzyn=sqrtu;
        quot=dzyn/znor;
        scale(z,quot);
        ode->set_Dirichlet(t,z);
    }
    else
    {
        dzyn=sqrtu*sqrt(z.size());
        add(z,sqrtu);
        ode->set_Dirichlet(t,z);
    }
    //here dzyn=||z-yn|| and z=yn+(small perturbation)
    //dzyn=||z-yn|| will be always true, even with new z in the loop

    //Start the power method for non linear operator rhs
    eigmax=0.0;
    for(int iter=1;iter<=maxiter;iter++)
    {
        ode->rhs(t,z,dz);
        n_f_eval_rho++;
        
        add(dz,-1.,fn); //dz is the new perturbation, not normalized yet
        ode->set_Dirichlet(t,dz,true);
        dfzfn= l2_norm(dz);
        
        eigmaxo=eigmax;
        eigmax=dfzfn/dzyn; //approximation of the Rayleigh quotient (not with dot product but just norms)
        eigmax=safe*eigmax;
        
        if (iter>=2 && std::abs(eigmax-eigmaxo)<= eigmax*tol)
        {
            //The last perturbation is stored. It will very likely be a
            // good starting point for the next rho call.
            dz=z;
            add(dz,-1.,yn);    
            return;
        }
        if (dfzfn!=0.0)
        {
            quot=dzyn/dfzfn;
            z=dz;
            scale(z,quot);
            add(z,yn); //z is built so that dzyn=||z-yn|| still true
        }
        else
            break;
    }
    std::cout<<"ERROR: Convergence failure in the spectral radius computation."<<std::endl;
    idid=-3;
}

template <class ODE>
void TimeIntegrator<ODE>::print_info()
{
  /*----------------------------------------- 
    Integration Parameters
  -----------------------------------------*/
  std::cout<<"\n-------   Integration parameters   ---------------"<<std::endl;
  std::cout<<"Time-step adaptivity "<<(dt_adaptivity ? "enabled.":"disabled.")<<std::endl;
  std::cout<<"Spectral radius computed "<<(internal_rho ? "internally.":"externally.")<<std::endl;
  std::cout<<"Absolute tolerance = "<<a_tol<<std::endl;
  std::cout<<"Relative tolerance = "<<r_tol<<std::endl;
  std::cout<<"--------------------------------------------------\n"<<std::endl;
}

template <class ODE> 
void TimeIntegrator<ODE>::print_statistics()
{
  /*----------------------------------------- 
    Integration Statistics
  -----------------------------------------*/
  std::cout<<"\n--------   Integration Statistics   ----------------"<<std::endl;
  std::cout<<"Max estimation of the spectral radius = "<<max_rho<<std::endl;
  std::cout<<"Min estimation of the spectral radius = "<<min_rho<<std::endl;
  std::cout<<"Max number of stages used = "<<max_s<<std::endl;
  std::cout<<"Number of f total evaluations = "<<n_f_eval<<std::endl;
  std::cout<<"Number of f eval. for the spectr. radius = "<<n_f_eval_rho<<std::endl;
  std::cout<<"Maximal time step used: "<<dt_max<<std::endl;
  std::cout<<"Steps: "<<n_steps<<std::endl;
  std::cout<<"Accepted steps: "<<acc_steps<<std::endl;
  std::cout<<"Rejected steps: "<<rej_steps<<std::endl; 
  std::cout<<"----------------------------------------------------\n"<<std::endl;
}

#endif	/* TIMEINTEGRATOR_H */
