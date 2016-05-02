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

#include "problem.h"
#include <cmath>
// Constructor of CABLE object: call constructor of ODE for one variable
CABLE::CABLE(Mesh& mesh,bool intrho,double t,std::vector<double>& init):ODE(mesh,intrho,t,init),n(NULL),m(NULL),h(NULL),n_L1(mesh.n_L1),n_L2(mesh.n_L2),n_L3(mesh.n_L3)
{
}

void CABLE::rhs(double t, std::vector<double>& y, std::vector<double>& f)
{
  // Point near the ending node in branch A
  f[1]=(10^3)/c_m*(((10^4)*a)/(3*s_step*s_step*R)*(2*y[2]-6*y[1]+4*y[0])	// membrane current density
       -g_K*((*n)[1]*(*n)[1]*(*n)[1]*(*n)[1])*(y[1]-e_K)			// potassium current contribution
       -g_Na*((*m)[1]*(*m)[1]*(*m)[1])*(*h)[1]*(y[1]-e_Na)			// sodium current contribution
       -g_l*(y[1]-e_l));							// leak current contribution
  for (std::size_t i=2; i<(n_L1-2); i++)
  {
  // Points in branch A far from boundaries
    f[i]=(10^3)/c_m*(((10^4)*a)/(2*s_step*s_step*R)*(y[i+1]-2*y[i]+y[i-1])	// membrane current density
       -g_K*((*n)[i]*(*n)[i]*(*n)[i]*(*n)[i])*(y[i]-e_K)			// potassium current contribution
       -g_Na*((*m)[i]*(*m)[i]*(*m)[i])*(*h)[i]*(y[i]-e_Na)			// sodium current contribution
       -g_l*(y[i]-e_l));							// leak current contribution         
  }
  // Point near the branching node in branch A
  f[n_L1-2]=(10^3)/c_m*(((10^4)*a)/(3*s_step*s_step*R)*(4*y[n_L1+n_L2-2]-6*y[n_L1-2]+2*y[n_L1-3])	// membrane current density
       -g_K*((*n)[n_L1-2]*(*n)[n_L1-2]*(*n)[n_L1-2]*(*n)[n_L1-2])*(y[n_L1-2]-e_K)			// potassium current contribution
       -g_Na*((*m)[n_L1-2]*(*m)[n_L1-2]*(*m)[n_L1-2])*(*h)[n_L1-2]*(y[n_L1-2]-e_Na)			// sodium current contribution
       -g_l*(y[n_L1-2]-e_l));										// leak current contribution
  // Point near the ending node in branch B
  f[n_L1]=(10^3)/c_m*(((10^4)*a)/(3*s_step*s_step*R)*(2*y[n_L1+1]-6*y[n_L1]+4*y[n_L1-1])		// membrane current density
       -g_K*((*n)[n_L1]*(*n)[n_L1]*(*n)[n_L1]*(*n)[n_L1])*(y[n_L1]-e_K)					// potassium current contribution
       -g_Na*((*m)[n_L1]*(*m)[n_L1]*(*m)[n_L1])*(*h)[n_L1]*(y[n_L1]-e_Na)				// sodium current contribution
       -g_l*(y[n_L1]-e_l));										// leak current contribution
  for (std::size_t i=n_L1+1; i<(n_L1+n_L2-3); i++)
  {
  // Points in branch B far from boundaries
    f[i]=(10^3)/c_m*(((10^4)*a)/(2*s_step*s_step*R)*(y[i+1]-2*y[i]+y[i-1])		// membrane current density
       -g_K*((*n)[i]*(*n)[i]*(*n)[i]*(*n)[i])*(y[i]-e_K)				// potassium current contribution
       -g_Na*((*m)[i]*(*m)[i]*(*m)[i])*(*h)[i]*(y[i]-e_Na)				// sodium current contribution
       -g_l*(y[i]-e_l)									// leak current contribution 
       -(t>t_off)*(i==(n_L1-2+round(n_L2/2)))*(10*10)*g_syn/(2*M_PI*a*s_step)*(t-t_off)/tau*exp(-(t-t_off-tau)/tau)*(y[i]-e_syn));	// spike contribution		<=========================
  }
  // Point near the branching node in branch B
  f[n_L1+n_L2-3]=(10^3)/c_m*(((10^4)*a)/(3*s_step*s_step*R)*(4*y[n_L1+n_L2-2]-6*y[n_L1+n_L2-3]+2*y[n_L1+n_L2-4])		// membrane current density
       -g_K*((*n)[n_L1+n_L2-3]*(*n)[n_L1+n_L2-3]*(*n)[n_L1+n_L2-3]*(*n)[n_L1+n_L2-3])*(y[n_L1+n_L2-3]-e_K)			// potassium current contribution
       -g_Na*((*m)[n_L1+n_L2-3]*(*m)[n_L1+n_L2-3]*(*m)[n_L1+n_L2-3])*(*h)[n_L1+n_L2-3]*(y[n_L1+n_L2-3]-e_Na)			// sodium current contribution
       -g_l*(y[n_L1+n_L2-3]-e_l));												// leak current contribution
  // Point near the branching node in branch C
  f[n_L1+n_L2-1]=(10^3)/c_m*(((10^4)*a)/(3*s_step*s_step*R)*(2*y[n_L1+n_L2]-6*y[n_L1+n_L2-1]+4*y[n_L1+n_L2-2])			// membrane current density
       -g_K*((*n)[n_L1+n_L2-1]*(*n)[n_L1+n_L2-1]*(*n)[n_L1+n_L2-1]*(*n)[n_L1+n_L2-1])*(y[n_L1+n_L2-1]-e_K)			// potassium current contribution
       -g_Na*((*m)[n_L1+n_L2-1]*(*m)[n_L1+n_L2-1]*(*m)[n_L1+n_L2-1])*(*h)[n_L1+n_L2-1]*(y[n_L1+n_L2-1]-e_Na)			// sodium current contribution
       -g_l*(y[n_L1+n_L2-1]-e_l));												// leak current contribution 
  for (std::size_t i=n_L1+n_L2; i<(n_L1+n_L2+n_L3-4); i++)
  {
  // Points in branch C far from boundaries
    f[i]=(10^3)/c_m*(((10^4)*a)/(2*s_step*s_step*R)*(y[i+1]-2*y[i]+y[i-1])	// membrane current density
       -g_K*((*n)[i]*(*n)[i]*(*n)[i]*(*n)[i])*(y[i]-e_K)			// potassium current contribution
       -g_Na*((*m)[i]*(*m)[i]*(*m)[i])*(*h)[i]*(y[i]-e_Na)			// sodium current contribution
       -g_l*(y[i]-e_l));							// leak current contribution 
  } 
  // Point near the ending node in branch C
  f[n_L1+n_L2+n_L3-4]=(10^3)/c_m*(((10^4)*a)/(3*s_step*s_step*R)*(4*y[n_L1+n_L2+n_L3-3]-6*y[n_L1+n_L2+n_L3-4]+2*y[n_L1+n_L2+n_L3-5])				// membrane current density
       -g_K*((*n)[n_L1+n_L2+n_L3-4]*(*n)[n_L1+n_L2+n_L3-4]*(*n)[n_L1+n_L2+n_L3-4]*(*n)[n_L1+n_L2+n_L3-4])*(y[n_L1+n_L2+n_L3-4]-e_K)				// potassium current contribution
       -g_Na*((*m)[n_L1+n_L2+n_L3-4]*(*m)[n_L1+n_L2+n_L3-4]*(*m)[n_L1+n_L2+n_L3-4])*(*h)[n_L1+n_L2+n_L3-41]*(y[n_L1+n_L2+n_L3-4]-e_Na)				// sodium current contribution
       -g_l*(y[n_L1+n_L2+n_L3-4]-e_l));																// leak current contribution 
}

void CABLE::get_gate_state(std::vector<double>& sodium,std::vector<double>& potassium,std::vector<double>& leak)
{
  n=&sodium;
  m=&potassium;
  h=&leak;
}

void CABLE::set_Dirichlet(double t, std::vector<double>& y, bool use)
{
    // Neumann boundary condition
    y[0]=y[1]; 
    y[n_L1-1]=y[n_L1];
    y[n_L1+n_L2+n_L3-3]=y[n_L1+n_L2+n_L3-4];
    // Branching boundary condition
    y[n_L1+n_L2-2]=y[n_L1-1]+y[n_L1+n_L2-3]-y[n_L1+n_L2-1];
}

GATE_N::GATE_N(Mesh& mesh,bool intrho,double t,std::vector<double>& init):ODE(mesh,intrho,t,init),alpha_n(0),beta_n(0)
{
}

void GATE_N::get_potential(std::vector<double>& y)
{
  v=&y;
}

void GATE_N::set_Dirichlet(double t, std::vector<double>& y, bool use)
{
}

void GATE_N::rhs(double t, std::vector<double>& y, std::vector<double>& f)
{
  for (unsigned int i=0; i<f.size(); i++)
    f[i]=0.01*(-((*v)[i]+55.0))/(exp(-((*v)[i]+55.0)/10.0)-1.0)	// alpha_n
         *(1.0-y[i])						// *(1-n)
         -0.125*exp(-((*v)[i]+65.0)/80.0)			// -beta_n
         *y[i];							// *n
}

GATE_M::GATE_M(Mesh& mesh,bool intrho,double t,std::vector<double>& init):ODE(mesh,intrho,t,init),alpha_m(0),beta_m(0)
{
}

void GATE_M::get_potential(std::vector<double>& y)
{
  v=&y;
}

void GATE_M::set_Dirichlet(double t, std::vector<double>& y, bool use)
{
}

void GATE_M::rhs(double t, std::vector<double>& y, std::vector<double>& f)
{
  for (unsigned int i=0; i<f.size(); i++)
    f[i]=0.1*(-((*v)[i]+40.0))/(exp(-((*v)[i]+40.0)/10.0)-1.0)	// alpha_m
         *(1.0-y[i])						// *(1-m)
         -4.0*exp(-((*v)[i]+65.0)/20.0)				// -beta_m
         *y[i];							// *m
}

GATE_H::GATE_H(Mesh& mesh,bool intrho,double t,std::vector<double>& init):ODE(mesh,intrho,t,init),alpha_h(0),beta_h(0)
{
}

void GATE_H::get_potential(std::vector<double>& y)
{
  v=&y;
}

void GATE_H::set_Dirichlet(double t, std::vector<double>& y, bool use)
{
}

void GATE_H::rhs(double t, std::vector<double>& y, std::vector<double>& f)
{
  for (unsigned int i=0; i<f.size(); i++)
    f[i]=0.07*exp(-((*v)[i]+65.0)/20.0)				// alpha_h
         *(1-y[i])						// *(1-h)
         -1.0/(exp(-((*v)[i]+35.0)/10.0)+1.0)			// -beta_h
         *y[i];							// *h
}
