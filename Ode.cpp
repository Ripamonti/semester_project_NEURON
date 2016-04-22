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

#include "Ode.h"
#include <cmath>

ODE::ODE(Mesh& mesh,bool intrho,double t) : time(0) , tend(t), cte_rho(true), s_step(mesh.h_space), limit_branch_A(mesh.n_L1), limit_branch_B(mesh.n_L2)
, limit_branch_C(mesh.n_L3){
  // Reserve space and define the initial condition
  un.reserve(mesh.n_elem);
  for (int i=0; i<mesh.n_elem; i++)
  {
    un.push_back(10*sin(mesh.grid[i]*2*M_PI/25));
  }
  // Reserve space for auxiliary variables used during the time step
  yn.resize(mesh.n_elem);
  ynpu.resize(mesh.n_elem);

  // Reserve space for auxiliary variables used during spectral radius computation
  dz1.resize(mesh.n_elem);
  z.resize(mesh.n_elem);
  fnpu.resize(mesh.n_elem);
  tmp1.resize(mesh.n_elem); // Probailmente da togliere..

  // Space for the evaluation
  fn.resize(mesh.n_elem);
}

void ODE::print_info()
{
  // TODO
}

std::vector<double>& ODE::get_un()
{
  return un;
}

void ODE::set_un(std::vector<double>& v)
{
  un=v;
}

void ODE::rho(double& eigmax)
{
  // TODO: For the future implement a method for approximate the spectral radius of the problem
  // For now, give a (rough) theoretical approximation
  eigmax= 4/(s_step*s_step);
}
 
void ODE::rhs(double t, std::vector<double>& y, std::vector<double>& f)
{
  // TODO: Implement a more efficient way to compute the rhs
  
  for (int i=1; i<(limit_branch_A-2); i++)
    f[i]=(y[i+1]-2*y[i]+y[i-1])/(s_step*s_step);
  f[limit_branch_A-2]=(y[limit_branch_A+limit_branch_B-2]-2*y[limit_branch_A-2]+y[limit_branch_A-3])/(s_step*s_step);
  for (int i=limit_branch_A; i<(limit_branch_A+limit_branch_B+limit_branch_C-2); i++)
    f[i]=(y[i+1]-2*y[i]+y[i-1])/(s_step*s_step);
}

void ODE::set_Dirichlet(double t, std::vector<double>& y, bool use)
{
  // TODO: Implement more general BD. Ask Giacomo the meaning of the variable use
  if (!use)
  {
    y[0]=0; 
    y[y.size()-1]=0;
    y[limit_branch_A-1]=0;
    y[limit_branch_A+limit_branch_B-2]=y[limit_branch_A+limit_branch_B-1]+y[limit_branch_A+limit_branch_B-3]-y[limit_branch_A-2];
  } else {
    y[0]=0; 
    y[y.size()-1]=0;
    y[limit_branch_A-2]=0;
  }
}

double ODE::normalized_L2_norm(std::vector<double>& u)
{
  // TODO: Ask how to compute this norm (discretize, continuous, meaning of normalized)
  // For now use L2 continuous norm
    double accum = 0.;
    for (double x : u) {
        accum += x * x;
    }
    return sqrt(s_step*accum);
}
/*
ODE::~ODE(void)
{
  
}
*/
