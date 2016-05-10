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

#include<ImplicitIntegrator.h>
#include<Ode.h>
#include<problem.h>
#include<vector>
#include <cmath>

CN::CN(CABLE* ode, bool onestep, bool verb):  n_L1(ode->limit_branch_A), n_L2(ode->limit_branch_B), n_L3(ode->limit_branch_C), n(ode->n), m(ode->m), h(ode->h), y(&(ode->un)), verbose(verb), one_step(onestep), s_step(ode->s_step)
{
  //  N.B. The upper and lower diagonals have one less element than the main diaognal 
  //       and the known term
  upper_diagonal.resize(ode->un.size()-1);
  lower_diagonal.resize(ode->un.size()-1);
  mid_diagonal.resize(ode->un.size());
  known_term.resize(ode->un.size());
  
  // Initialize the upper and lower diagonal of the system to be solved
  // They are constant through all the simulation

  // Neumann condiiton branch B2
  upper_diagonal[0]=-2/s_step;
  // Point near the ending node in branch B2
  upper_diagonal[1]=-(10^4)*a/(3*R*s_step*s_step);
  lower_diagonal[0]=-2*(10^4)*a/(3*R*s_step*s_step);
  // Points in branch B2 far from boundaries
  for (std::size_t i=2; i<(n_L1-2); i++)
  {
    upper_diagonal[i]=-(10^4)*a/(4*R*s_step*s_step);
    lower_diagonal[i-1]=-(10^4)*a/(4*R*s_step*s_step);  
  }
  // Point near the branching node in branch B2
  upper_diagonal[n_L1-2]=-2*(10^4)*a/(3*R*s_step*s_step); // <====== Attention: we put this element, which should be far from upper diagonal, in this position in order to store everything in the same data structure
  lower_diagonal[n_L1-3]=-(10^4)*a/(3*R*s_step*s_step);
  // Neumann condition in branch B1
  upper_diagonal[n_L1-1]=-2/s_step;
  lower_diagonal[n_L1-2]=-2/s_step; // <====== Attention: we put this element, which should be far from lower diagonal, in this posiiton in order to store everythin in the same data structure
  // Point near the ending node in branch B1
  upper_diagonal[n_L1]=-(10^4)*a/(3*R*s_step*s_step);
  lower_diagonal[n_L1-1]=-2*(10^4)*a/(3*R*s_step*s_step);
  // Points in branch B1 far from boundaries
  for (std::size_t i=n_L1+1; i<(n_L1+n_L2-3); i++)
  {
    upper_diagonal[i]=-(10^4)*a/(4*R*s_step*s_step);
    lower_diagonal[i-1]=-(10^4)*a/(4*R*s_step*s_step);     
  }
  // Point near the branching node in branch B1
  upper_diagonal[n_L1+n_L2-3]=-2*(10^4)*a/(3*R*s_step*s_step);
  lower_diagonal[n_L1+n_L2-4]=-(10^4)*a/(3*R*s_step*s_step);
  // Branching node
  upper_diagonal[n_L1+n_L2-2]=2/s_step;
  lower_diagonal[n_L1+n_L2-3]=-2/s_step;
  // Point near the branching node in branch A
  upper_diagonal[n_L1+n_L2-1]=-(10^4)*a/(3*R*s_step*s_step);
  lower_diagonal[n_L1+n_L2-2]=-2*(10^4)*a/(3*R*s_step*s_step);
  // Points in branch A far from boundaries
  for (std::size_t i=n_L1+n_L2; i<(n_L1+n_L2+n_L3-4); i++)
  {
    upper_diagonal[i]=-(10^4)*a/(4*R*s_step*s_step);
    lower_diagonal[i-1]=-(10^4)*a/(4*R*s_step*s_step);     
  } 
  // Point near the ending node in branch A
  upper_diagonal[n_L1+n_L2+n_L3-4]=-2*(10^4)*a/(3*R*s_step*s_step);
  lower_diagonal[n_L1+n_L2+n_L3-5]=-(10^4)*a/(3*R*s_step*s_step);
  // Neumann condition in branch A
  lower_diagonal[n_L1+n_L2+n_L3-4]=-2/s_step;   
}

void CN::build_diag(CABLE* ode, double dt)
{

  // Neumann condition in branch B2
  mid_diagonal[0]=2/s_step;
  known_term[0]=0;
  // Point near the ending node in branch B2
  mid_diagonal[1]=(10^(-3))*c_m/dt
                 +(10^4)*a/(R*s_step*s_step)
                 +0.5*g_K*((*n)[1]*(*n)[1]*(*n)[1]*(*n)[1])
                 +0.5*g_Na*((*m)[1]*(*m)[1]*(*m)[1])*(*h)[1]
                 +0.5*g_l;
  known_term[1]=(10^(-3))*c_m/dt*(*y)[1]
               +(10^4)*a/(6*R*s_step*s_step)*(2*(*y)[2]-6*(*y)[1]+4*(*y)[0])
               -g_K*((*n)[1]*(*n)[1]*(*n)[1]*(*n)[1])*(0.5*(*y)[1]-e_K)
               -g_Na*((*m)[1]*(*m)[1]*(*m)[1])*(*h)[1]*(0.5*(*y)[1]-e_Na)
               -g_l*(0.5*(*y)[1]-e_l);
  // Points in branch B2 far from boundaries
  for (std::size_t i=2; i<(n_L1-2); i++)
  {
    mid_diagonal[i]=(10^(-3))*c_m/dt
                   +(10^4)*a/(2*R*s_step*s_step)   
                   +0.5*g_K*((*n)[i]*(*n)[i]*(*n)[i]*(*n)[i])
                   +0.5*g_Na*((*m)[i]*(*m)[i]*(*m)[i])*(*h)[i]
                   +0.5*g_l;             
    known_term[i]=(10^(-3))*c_m/dt*(*y)[i]
                 +(10^4)*a/(4*R*s_step*s_step)*((*y)[i+1]-2*(*y)[i]+(*y)[i-1])
                 -g_K*((*n)[i]*(*n)[i]*(*n)[i]*(*n)[i])*(0.5*(*y)[i]-e_K)
                 -g_Na*((*m)[i]*(*m)[i]*(*m)[i])*(*h)[i]*(0.5*(*y)[i]-e_Na)
                 -g_l*(0.5*(*y)[i]-e_l);
  }
  // Point near the branching node in branch B2
  mid_diagonal[n_L1-2]=(10^(-3))*c_m/dt
                      +(10^4)*a/(R*s_step*s_step)
                      +0.5*g_K*((*n)[n_L1-2]*(*n)[n_L1-2]*(*n)[n_L1-2]*(*n)[n_L1-2])
                      +0.5*g_Na*((*m)[n_L1-2]*(*m)[n_L1-2]*(*m)[n_L1-2])*(*h)[n_L1-2]
                      +0.5*g_l;
  known_term[n_L1-2]=(10^(-3))*c_m/dt*(*y)[n_L1-2]
                    +(10^4)*a/(6*R*s_step*s_step)*(4*(*y)[n_L1-1]-6*(*y)[n_L1-2]+2*(*y)[n_L1-3])
                    -g_K*((*n)[n_L1-2]*(*n)[n_L1-2]*(*n)[n_L1-2]*(*n)[n_L1-2])*(0.5*(*y)[n_L1-2]-e_K)
                    -g_Na*((*m)[n_L1-2]*(*m)[n_L1-2]*(*m)[n_L1-2])*(*h)[n_L1-2]*(0.5*(*y)[n_L1-2]-e_Na)
                    -g_l*(0.5*(*y)[n_L1-2]-e_l);
  // Neumann condition in branch B1
  mid_diagonal[n_L1-1]=2/s_step;
  known_term[n_L1-1]=0;
  // Point near the ending node in branch B1
  mid_diagonal[n_L1]=(10^(-3))*c_m/dt
                 +(10^4)*a/(R*s_step*s_step)
                 +0.5*g_K*((*n)[n_L1]*(*n)[n_L1]*(*n)[n_L1]*(*n)[n_L1])
                 +0.5*g_Na*((*m)[n_L1]*(*m)[n_L1]*(*m)[n_L1])*(*h)[n_L1]
                 +0.5*g_l;
  known_term[n_L1]=(10^(-3))*c_m/dt*(*y)[n_L1]
               +(10^4)*a/(6*R*s_step*s_step)*(2*(*y)[n_L1+1]-6*(*y)[n_L1]+4*(*y)[n_L1-1])
               -g_K*((*n)[n_L1]*(*n)[n_L1]*(*n)[n_L1]*(*n)[n_L1])*(0.5*(*y)[n_L1]-e_K)
               -g_Na*((*m)[n_L1]*(*m)[n_L1]*(*m)[n_L1])*(*h)[n_L1]*(0.5*(*y)[n_L1]-e_Na)
               -g_l*(0.5*(*y)[n_L1]-e_l);
  // Points in branch B1 far from boundaries
  for (std::size_t i=n_L1+1; i<(n_L1+n_L2-3); i++)
  {
    mid_diagonal[i]=(10^(-3))*c_m/dt
                   +(10^4)*a/(2*R*s_step*s_step)   
                   +0.5*g_K*((*n)[i]*(*n)[i]*(*n)[i]*(*n)[i])
                   +0.5*g_Na*((*m)[i]*(*m)[i]*(*m)[i])*(*h)[i]
                   +0.5*g_l
                   +((ode->time)>t_off)*(i==(n_L1-2+round(n_L2/2)))*(10*10)*g_syn/(2*M_PI*a*s_step)*((ode->time)-t_off)/tau*exp(-((ode->time)-t_off-tau)/tau); // <====== spike           
    known_term[i]=(10^(-3))*c_m/dt*(*y)[i]
                 +(10^4)*a/(4*R*s_step*s_step)*((*y)[i+1]-2*(*y)[i]+(*y)[i-1])
                 -g_K*((*n)[i]*(*n)[i]*(*n)[i]*(*n)[i])*(0.5*(*y)[i]-e_K)
                 -g_Na*((*m)[i]*(*m)[i]*(*m)[i])*(*h)[i]*(0.5*(*y)[i]-e_Na)
                 -g_l*(0.5*(*y)[i]-e_l)
                 -((ode->time-dt)>t_off)*(i==(n_L1-2+round(n_L2/2)))*(10*10)*g_syn/(2*M_PI*a*s_step)*((ode->time-dt)-t_off)/tau*exp(-((ode->time -dt)-t_off-tau)/tau)*(0.5*((*y)[i]-e_syn)); // <====== spike
  }
  // Point near the branching node in branch B1
  mid_diagonal[n_L1+n_L2-3]=(10^(-3))*c_m/dt
                           +(10^4)*a/(R*s_step*s_step)
                           +0.5*g_K*((*n)[n_L1+n_L2-3]*(*n)[n_L1+n_L2-3]*(*n)[n_L1+n_L2-3]*(*n)[n_L1+n_L2-3])
                           +0.5*g_Na*((*m)[n_L1+n_L2-3]*(*m)[n_L1+n_L2-3]*(*m)[n_L1+n_L2-3])*(*h)[n_L1+n_L2-3]
                           +0.5*g_l;
  known_term[n_L1+n_L2-3]=(10^(-3))*c_m/dt*(*y)[n_L1+n_L2-3]
                         +(10^4)*a/(6*R*s_step*s_step)*(4*(*y)[n_L1+n_L2-2]-6*(*y)[n_L1+n_L2-3]+2*(*y)[n_L1+n_L2-4])
                         -g_K*((*n)[n_L1+n_L2-3]*(*n)[n_L1+n_L2-3]*(*n)[n_L1+n_L2-3]*(*n)[n_L1+n_L2-3])*(0.5*(*y)[n_L1+n_L2-3]-e_K)
                         -g_Na*((*m)[n_L1+n_L2-3]*(*m)[n_L1+n_L2-3]*(*m)[n_L1+n_L2-3])*(*h)[n_L1+n_L2-3]*(0.5*(*y)[n_L1+n_L2-3]-e_Na)
                         -g_l*(0.5*(*y)[n_L1+n_L2-3]-e_l);
  // Branching node
  mid_diagonal[n_L1+n_L2-2]=2/s_step;
  known_term[n_L1+n_L2-2]=0; 
  // Point near the branching node in branch A
  mid_diagonal[n_L1+n_L2-1]=(10^(-3))*c_m/dt
                 +(10^4)*a/(R*s_step*s_step)
                 +0.5*g_K*((*n)[n_L1+n_L2-1]*(*n)[n_L1+n_L2-1]*(*n)[n_L1+n_L2-1]*(*n)[n_L1+n_L2-1])
                 +0.5*g_Na*((*m)[n_L1+n_L2-1]*(*m)[n_L1+n_L2-1]*(*m)[n_L1+n_L2-1])*(*h)[n_L1+n_L2-1]
                 +0.5*g_l;
  known_term[n_L1+n_L2-1]=(10^(-3))*c_m/dt*(*y)[n_L1+n_L2-1]
               +(10^4)*a/(6*R*s_step*s_step)*(2*(*y)[n_L1+n_L2]-6*(*y)[n_L1+n_L2-1]+4*(*y)[n_L1+n_L2-2])
               -g_K*((*n)[n_L1+n_L2-1]*(*n)[n_L1+n_L2-1]*(*n)[n_L1+n_L2-1]*(*n)[n_L1+n_L2-1])*(0.5*(*y)[n_L1+n_L2-1]-e_K)
               -g_Na*((*m)[n_L1+n_L2-1]*(*m)[n_L1+n_L2-1]*(*m)[n_L1+n_L2-1])*(*h)[n_L1+n_L2-1]*(0.5*(*y)[n_L1+n_L2-1]-e_Na)
               -g_l*(0.5*(*y)[n_L1+n_L2-1]-e_l); 
  // Points in branch A far from boundaries
  for (std::size_t i=n_L1+n_L2; i<(n_L1+n_L2+n_L3-4); i++)
  {
    mid_diagonal[i]=(10^(-3))*c_m/dt
                   +(10^4)*a/(2*R*s_step*s_step)   
                   +0.5*g_K*((*n)[i]*(*n)[i]*(*n)[i]*(*n)[i])
                   +0.5*g_Na*((*m)[i]*(*m)[i]*(*m)[i])*(*h)[i]
                   +0.5*g_l;             
    known_term[i]=(10^(-3))*c_m/dt*(*y)[i]
                 +(10^4)*a/(4*R*s_step*s_step)*((*y)[i+1]-2*(*y)[i]+(*y)[i-1])
                 -g_K*((*n)[i]*(*n)[i]*(*n)[i]*(*n)[i])*(0.5*(*y)[i]-e_K)
                 -g_Na*((*m)[i]*(*m)[i]*(*m)[i])*(*h)[i]*(0.5*(*y)[i]-e_Na)
                 -g_l*(0.5*(*y)[i]-e_l);
  }  
  // Point near the ending node in branch A
  mid_diagonal[n_L1+n_L2+n_L3-4]=(10^(-3))*c_m/dt
                           +(10^4)*a/(R*s_step*s_step)
                           +0.5*g_K*((*n)[n_L1+n_L2+n_L3-4]*(*n)[n_L1+n_L2+n_L3-4]*(*n)[n_L1+n_L2+n_L3-4]*(*n)[n_L1+n_L2+n_L3-4])
                           +0.5*g_Na*((*m)[n_L1+n_L2+n_L3-4]*(*m)[n_L1+n_L2+n_L3-4]*(*m)[n_L1+n_L2+n_L3-4])*(*h)[n_L1+n_L2+n_L3-4]
                           +0.5*g_l;
  known_term[n_L1+n_L2+n_L3-4]=(10^(-3))*c_m/dt*(*y)[n_L1+n_L2+n_L3-4]
                         +(10^4)*a/(6*R*s_step*s_step)*(4*(*y)[n_L1+n_L2+n_L3-3]-6*(*y)[n_L1+n_L2+n_L3-4]+2*(*y)[n_L1+n_L2+n_L3-5])
                         -g_K*((*n)[n_L1+n_L2+n_L3-4]*(*n)[n_L1+n_L2+n_L3-4]*(*n)[n_L1+n_L2+n_L3-4]*(*n)[n_L1+n_L2+n_L3-4])*(0.5*(*y)[n_L1+n_L2+n_L3-4]-e_K)
                         -g_Na*((*m)[n_L1+n_L2+n_L3-4]*(*m)[n_L1+n_L2+n_L3-4]*(*m)[n_L1+n_L2+n_L3-4])*(*h)[n_L1+n_L2+n_L3-4]*(0.5*(*y)[n_L1+n_L2+n_L3-4]-e_Na)
                         -g_l*(0.5*(*y)[n_L1+n_L2+n_L3-4]-e_l);
  // Neumann condition in branch A
  mid_diagonal[n_L1+n_L2+n_L3-3]=2/s_step;
  known_term[n_L1+n_L2+n_L3-3]=0;
}

void CN::advance(CABLE* ode, double& h, int& idid)
{
  ode->time = ode->time+h;
  build_diag(ode, h);
  modified_Thomas();
  
  // Part to be fixed: now it has the same structure of the explicit method, but while in that case mantaining all these copies of the solution is necessary, here is only a waste of time
  ode->un=mid_diagonal;
  ode->yn=mid_diagonal;
  ode->ynpu=mid_diagonal;
}

void CN::modified_Thomas()
{
  // Upper triangularization
  for (std::size_t i=1; i<n_L1-1; i++)
  {
    mid_diagonal[i]=mid_diagonal[i]-upper_diagonal[i-1]*(lower_diagonal[i-1]/mid_diagonal[i-1]);
    known_term[i]=known_term[i]-known_term[i-1]*(lower_diagonal[i-1]/mid_diagonal[i-1]);
  }
  mid_diagonal[n_L1+n_L2-2]=mid_diagonal[n_L1+n_L2-2]-upper_diagonal[n_L1-2]*(lower_diagonal[n_L1-2]/mid_diagonal[n_L1-2]);
  known_term[n_L1+n_L2-2]=known_term[n_L1+n_L2-2]-known_term[n_L1-2]*(lower_diagonal[n_L1-2]/mid_diagonal[n_L1-2]);
 for (std::size_t i=n_L1; i<n_L1+n_L2+n_L3-2; i++)
 {
    mid_diagonal[i]=mid_diagonal[i]-upper_diagonal[i-1]*(lower_diagonal[i-1]/mid_diagonal[i-1]);
    known_term[i]=known_term[i]-known_term[i-1]*(lower_diagonal[i-1]/mid_diagonal[i-1]); 
 }
 // Lower triangularization
  mid_diagonal[n_L1+n_L2+n_L3-3]=known_term[n_L1+n_L2+n_L3-3]/mid_diagonal[n_L1+n_L2+n_L3-3];
  for (std::size_t i=n_L1+n_L2+n_L3-4; i>n_L1-2; i--)
  {
    mid_diagonal[i]=(known_term[i]-mid_diagonal[i+1]*upper_diagonal[i])/mid_diagonal[i];
  }
  mid_diagonal[n_L1-2]=(known_term[n_L1-2]-mid_diagonal[n_L1+n_L2-2]*upper_diagonal[n_L1-2])/mid_diagonal[n_L1-2]; 
   for (std::size_t i=n_L1-3; i>=1; i--)
  {
    mid_diagonal[i]=(known_term[i]-mid_diagonal[i+1]*upper_diagonal[i])/mid_diagonal[i];
  }
  mid_diagonal[0]=(known_term[0]-mid_diagonal[1]*upper_diagonal[0])/mid_diagonal[0];
}
