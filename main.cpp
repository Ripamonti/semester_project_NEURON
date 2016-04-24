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

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdio.h>

#include "TimeIntegrator.h"
#include "ChebyshevIntegrators.h"
#include "Ode.h"
#include "Mesh.h"
#include "problem.h"

using namespace std;

int main(int argc, char** argv)
{
//----------------    DEFAULT ACCURACY PARAMETERS   ----------------------------
    int n_ref = 2;    		//mesh size or mesh name depending if non local or local
    double dt =  0.025;		//initial step size
    double rtol= 1.0e-2;	//relative tolerance
    double atol= rtol; 	      	//absolute tolerance

//----------------    DEFAULT INTEGRATION OPTIONS   ----------------------------
    bool dt_adaptivity = false; 	//time step adaptivity enabled/disabled
    bool one_step=true;         	//advance one time step per ROCK2 call
    bool intrho=true;			//spectral radius computed internally in the time integrator
    double tend=10;			//last time

//-----    READ OPTIONAL INPUT ACCURACY/INTEGRATION PARAMETERS   ---------------
    if(argc>1){
        istringstream Nbuf(argv[1]);
        istringstream dtbuf(argv[2]);
        Nbuf>>n_ref;
        dtbuf>>dt;
        dt_adaptivity= (*argv[3]=='1');
        one_step= (*argv[4]=='1');
    }

//-------------------    MESH INITIALIZATION  ----------------------------------
    double dx = 0.5;			//grid space
    Mesh mesh(dx);			//initialize the mesh
    mesh.print_info();

//-----------------    INITIAL DATA DEFINITION  --------------------------------
    vector<double> pot_initial(mesh.n_elem,-64.974);
    double alpha_n(0.01*(-(-64.974+55.0))/(exp(-(-64.974+55)/10)-1));
    double beta_n(0.125*exp(-(-64.974+65)/80));
    double alpha_m(0.1*(-(64.974+40))/(exp(-(-64.974+40)/10)-1));
    double beta_m(0.4*exp(-(-64.974+65)/18));
    double alpha_h(0.07*exp(-(-64.974+65)/20));
    double beta_h(1/(exp(-(-64.974+35)/10)+1));
    vector<double> n_initial(mesh.n_elem,alpha_n/(alpha_n+beta_n));
    vector<double> m_initial(mesh.n_elem,alpha_m/(alpha_m+beta_m));
    vector<double> h_initial(mesh.n_elem,alpha_h/(alpha_h+beta_h));

//-------------------    PROBLEM INITIALIZATION  -------------------------------
    CABLE* cable = new CABLE(mesh,intrho,tend,pot_initial);
    GATE_N* gate_n = new GATE_N(mesh,intrho,tend,n_initial);
    GATE_M* gate_m = new GATE_M(mesh,intrho,tend,m_initial);
    GATE_H* gate_h = new GATE_H(mesh,intrho,tend,h_initial);
    gate_n->get_potential(cable->un);
    gate_m->get_potential(cable->un);
    gate_h->get_potential(cable->un);
    cable->get_gate_state(gate_n->un,gate_m->un,gate_h->un);

//---------------------    ROCK2/RKC INITIALIZATION  ---------------------------
    bool verbose=true;
    ROCK2 rock_cable(one_step, verbose, dt_adaptivity, atol, rtol, intrho); 
    rock_cable.print_info();
    verbose=true;
    ROCK2 rock_gate_n(one_step, verbose, dt_adaptivity, atol, rtol, intrho);  
    ROCK2 rock_gate_m(one_step, verbose, dt_adaptivity, atol, rtol, intrho);
    ROCK2 rock_gate_h(one_step, verbose, dt_adaptivity, atol, rtol, intrho);   
//-------------------------   TIME LOOP   --------------------------------------

    if(rock_cable.check_correctness(dt)==0 && rock_gate_n.check_correctness(dt)==0
       && rock_gate_m.check_correctness(dt)==0 && rock_gate_h.check_correctness(dt)==0) //checks if given parameters are ok
        return 0;//something is wrong
    int idid = 2;
    dt=dt/2;
    rock_gate_n.advance(gate_n,dt,idid);
    rock_gate_m.advance(gate_m,dt,idid);
    rock_gate_h.advance(gate_h,dt,idid);
    dt=dt*2;    
    idid=2;
fflush(stdout);

    for(int iter=0;idid==2||cable->time<tend;++iter)
    {
        rock_cable.advance(cable,dt,idid);
        rock_gate_n.advance(gate_n,dt,idid);
        rock_gate_m.advance(gate_m,dt,idid);
        rock_gate_h.advance(gate_h,dt,idid);
        //printf(" %f ",cable->time);fflush(stdout);
    }   
      rock_gate_n.print_statistics();
//---------------------------   WRITE ON SCREEN   ----------------------------------------
   bool print_GNU=true;
   if (print_GNU)
   {
   FILE * branchA;
   FILE * branchB;
   FILE * branchC;
   branchA = fopen ("branchA.txt","w");
   branchB = fopen ("branchB.txt","w");
   branchC = fopen ("branchC.txt","w");
    for (int i=0; i<(mesh.n_L1-1); i++){
      fprintf(branchA,"%.8f \t %.8f \n",mesh.grid[i],cable->un[i]);
    }
   fclose (branchA);
    for (int i=(mesh.n_L1-1); i<(mesh.n_L1+mesh.n_L2-2); i++){
      fprintf(branchB,"%.8f \t %.8f \n",mesh.grid[i],cable->un[i]);
    }
   fclose (branchB);
    for (int i=(mesh.n_L1+mesh.n_L2-2); i<mesh.n_elem; i++){
      fprintf(branchC,"%.8f \t %.8f \n",mesh.grid[i],cable->un[i]);
    }
   fclose (branchC);
    }
//------------------------   DEALLOCATE   -------------------------------------    
//    delete heat;
    return 0;
}

