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

#include "TimeIntegrator.h"

using namespace std;

template <class ODE>
TimeIntegrator<ODE>::TimeIntegrator(bool onestep, bool verb, bool dtadap, double atol, double rtol, bool intrho, bool scalartol)
{
  uround=numeric_limits<double>::epsilon();
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

  max_rho=numeric_limits<int>::min();
  min_rho=numeric_limits<int>::max();
  max_s=0;
  n_f_eval_rho=0;
  n_f_eval=0;
  n_steps=0;
  acc_steps=0;
  rej_steps=0;
  dt_max=0.;
  dt_min=numeric_limits<double>::max();
}

template <class ODE>
void TimeIntegrator<ODE>::advance(ODE* ode, double& h, int& idid)
{
  double told;
  vector<double> *swap_ptr;

  double& t=ode->time;	// If t is modified then ode->time as well
  double tend=ode->tend;

  vector<double>& y = ode->ynpu;	// Space for next computed solution
  vector<double>& yn = ode->yn;  // Space for initial value 

  told=t-1.;			// Different from t
  last=false;		// A priori is not the last step
  reject=false;		// If we enter advance then the previous step has been accepted
  idid=1;

  // If it is the first step then we:
  // i)  recover the initial value from ode
  // ii) compute the spectral radius of the system

  if (n_steps==0)
  {
    yn=ode->get_un();		// Takes the initial value from ODE (da sistemare)
    fn_uptodate=false;		// We have a new yn, so fn=f(yn) has to be updated
    update_rho(ode,idid);	// Computes spectral radius
    
  }
  // else the last solution is already contained in yn and we already have an approximation of the spectral radius

  /*----------------------------------------- 
    Integration loop 
  -----------------------------------------*/
  for(;;)
  {
    n_steps++;			// New time step
    // If we are very close to end or t+h>tend then adjust h
    if (1.1*abs(h)>=abs(tend-t))	// h can be negative (when t>tend)
    {
      h=tend-t;
      last=true;		// It will be the last step (if accepted)
    }
    if (h<10.0*uround)		// Too small h
    {
      cout << "Tolerances are too small." << endl;
      idid=-2;			// idid=-2 means that there's a problem
      return; 			// Exit
    }

    // Do some statistics
    dt_max=max(dt_max,h);
    dt_min=min(dt_min,h);
     
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

    if (dt_adaptivity) // If time step adaptivity enabled
      compute_hnew(h);  // Computed new h
    else // We are not using time step adaptivity
      hnew=h; // keeps same h

    // Accepted step or without time step adaptivity
    if (err<1.0 || !dt_adaptivity)
    {
      accepted_step(t,h,y);
      
      /*swap_ptr=yn;
      yn=y;			// yn takes the new value
      y=swap_ptr;*/		// and y will be free free memory
      yn.swap(y);
      fn_uptodate=false;	// fn not up to date anymore

      if (last)  // We have finished
      {
        // remember that idid=1
        ode->set_un(yn);	// Pass value to ODE, for outputs
        cout<<endl;
        break;
      }
    } else if(one_step) {
      ode->set_un(yn);		//pass value to ode, for outputs
      idid=2;                   //says to outside world that we didn't reach the end
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
  if (!internal_rho)
    ode->rho(eigmax);
  // Computed internally by TimeIntegrator::rho
  else
    this->rho(eigmax,ode,idid);
  
  // Recover statistics
  if ((int)eigmax>max_rho) max_rho=(int)eigmax;
  if ((int)eigmax<min_rho) min_rho=(int)eigmax;

  // Print info
/*
  cout<<"-------- Spectral radius estimations --------"<<endl;
  cout<<"Spectral radius = "<<(int)eigmax<<endl;
  cout<<"---------------------------------------------\n"<<endl;   
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
    fac=min(facp,fac);
    fac=sqrt(fac);
  } 
  else // In first step we go here, or if the previous step has been rejected
    fac=sqrt(1.0/err);
  fac=min(facmax,max(0.1,0.8*fac));  
  hnew=h*fac;
}

template <class ODE>
void TimeIntegrator<ODE>::accepted_step(double& t, double& h, vector<double>& y)
{
  // Called when the step is accepted. Does output, sets the new time step size and update some variables
  // Define delta character and red color
  std::string delta=u8"\u0394";
  string bcol="\033[31;1m";
  string ecol="\033[0m";
  if (dt_adaptivity && verbose)
    cout<<"Step t = "<<t<<", "<<delta<<"t = "<<h
    <<", s = "<<s<<". Acceptd with ||e_n||_L2 = "<<err<<endl;
    /* <<" and ||y_n||_Linf = "<<y */
  else if (verbose)
    cout<<"Step t = "<<t<<", "<<delta<<"t = "<<h
    <<", s = "<<s<<","<<(err>1? bcol:"")<<" ||e_n||_L2 = "<<err<<(err>1? ecol:"")<<endl;
    /* si dovrebbe mettere la norma infinito.. */
  
  acc_steps++;
  facmax=2.0;		// h can grow 2 times
  t=t+h;
  errp=err;
  if (reject)  // The previous time step has been rejected so a smaller h is chosen
  {
    hnew=h>0.0 ? min(hnew,h):max(hnew,h);
    reject=false;  // This step has been accepted
    nrej=0;		 // Set the consecutive rejections to 0
  }
  hp=h;			// Previous h
  h=hnew;		// Next h
  nrho=nrho+1;		// Consecutive steps without computing the spectral radius
  nrho=(nrho+1)%25;	// Set to 0 every 25 steps
}

template <class ODE>
void TimeIntegrator<ODE>::rejected_step(double& t, double& h, vector<double>& y)
{
  // Called when the step is rejected. Does output, choses the new time step and updates some variables
  // Define delta character and red color
  std::string delta=u8"\u0394";
  string bcol="\033[31;1m";
  string ecol="\033[0m";
  if (verbose)
   cout<<"Step t = "<<t<<", "<<delta<<"t = "<<h
   <<", s = "<<s<<"."<<bcol<<" Refected with ||e_n||_L2 = "<<err<<ecol<<endl;
   /* anche qui bisognerebbe mettere la norma infinito */ 

  rej_steps++;
  reject=true;
  last=false;	// The step is rejected, it can't be the last one
  h=0.8*hnew;
  facmax=2.0;		// Next step h can at most double
  if (n_steps==0)	// Here we are in the first step, starting time was too large
    h=0.1*h;
  if(told==t)		// The previous step was also rejected
  {
    nrej=nrej+1;	// Consecutive rejections
    if (nrej==10)	// after 10 consecutive rejections ..
      h=h*1e-3;
  }
  else
    told=t;		// First rejection

  // The spectral radius is recomputed after a step failure
  if (nrho)
    nrho=0;		// Spectral radius will be computed
  else
    nrho=1;		// Was already computed
}

template <class ODE>
int TimeIntegrator<ODE>::check_correctness(double& h)
{
  // Some parameters checking before starting integration
  
  // Test the initial step size and tolerances
  if (h<10.0*uround)
  {
    cout<<"Initial step-size is too small."<<endl;
    return 0;
  }
  if (scalar_tol)	// Scalar tolerances
  {
    if (a_tol<=0.0 || r_tol<=10.0*uround)
    {
      cout<<"Tolerances are too small."<<endl;
      return 0;
    }
  }
  else
  {
    cerr<<"NON SCALAR TOLERANCES NOT IMPLEMENTED YET"<<endl;
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

    vector<double>& yn = ode->yn;
    vector<double>& fn = ode->fn;
    vector<double>& dz = ode->dz1;
    vector<double>& z  = ode->tmp1;

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

    // D A  S I S T E M A R E
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
        
        if (iter>=2 && abs(eigmax-eigmaxo)<= eigmax*tol)
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
    cout<<"ERROR: Convergence failure in the spectral radius computation."<<endl;
    idid=-3;
}

template <class ODE>
void TimeIntegrator<ODE>::print_info()
{
  /*----------------------------------------- 
    Integration Parameters
  -----------------------------------------*/
  cout<<"\n-------   Integration parameters   ---------------"<<endl;
  cout<<"Time-step adaptivity "<<(dt_adaptivity ? "enabled.":"disabled.")<<endl;
  cout<<"Spectral radius computed "<<(internal_rho ? "internally.":"externally.")<<endl;
  cout<<"Absolute tolerance = "<<a_tol<<endl;
  cout<<"Relative tolerance = "<<r_tol<<endl;
  cout<<"--------------------------------------------------\n"<<endl;
}

template <class ODE> 
void TimeIntegrator<ODE>::print_statistics()
{
  /*----------------------------------------- 
    Integration Statistics
  -----------------------------------------*/
  cout<<"\n--------   Integration Statistics   ----------------"<<endl;
  cout<<"Max estimation of the spectral radius = "<<max_rho<<endl;
  cout<<"Min estimation of the spectral radius = "<<min_rho<<endl;
  cout<<"Max number of stages used = "<<max_s<<endl;
  cout<<"Number of f total evaluations = "<<n_f_eval<<endl;
  cout<<"Number of f eval. for the spectr. radius = "<<n_f_eval_rho<<endl;
  cout<<"Maximal time step used: "<<dt_max<<endl;
  cout<<"Steps: "<<n_steps<<endl;
  cout<<"Accepted steps: "<<acc_steps<<endl;
  cout<<"Rejected steps: "<<rej_steps<<endl; 
  cout<<"----------------------------------------------------\n"<<endl;
}
