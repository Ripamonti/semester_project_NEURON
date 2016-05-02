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

#include <algorithm>
#include "operation.h"
#include <numeric>

using namespace std;

void add(vector<double>& v, double c, vector<double>& w)
{
  transform(v.begin(),v.end(),w.begin(),v.begin(),[c](double v, double w){return v+c*w;});
}

void add(vector<double>& v, double c)
{
transform(v.begin(), v.end(), v.begin(),bind2nd(std::plus<double>(), c));
}

void add(vector<double>& v,vector<double>& w)
{
  transform(v.begin(),v.end(),w.begin(),v.begin(),[](double v, double w){return v+w;});
}

void scale(vector<double>& v, double c)
{
  transform(v.begin(),v.end(),v.begin(),bind1st(multiplies<double>(),c));
}

double l2_norm(vector<double> const& u) {
return sqrt(std::inner_product(u.begin(), u.end(),u.begin(),0.0));
}

double l2_cont_norm(vector<double> const& u, double c) {
return sqrt(c*std::inner_product(u.begin(), u.end(),u.begin(),0.0));
}
