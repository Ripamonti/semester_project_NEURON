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

using namespace std;

int main(int argc, char** argv)
{
//----------------    DEFAULT ACCURACY PARAMETERS   ----------------------------
    int n_ref = 2;    		//mesh size or mesh name depending if non local or local
    double dt =  0.025;		//initial step size
    double rtol= 1.0e-2;	//relative tolerance
    double atol= rtol; 	      	//absolute tolerance
    bool print_GNU=true;

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

//-------------------    PROBLEM INITIALIZATION  -------------------------------
    ODE* cable = new ODE(mesh,intrho,tend);
    ODE* gates = new ODE(mesh,intrho,tend);

/*
//---------------------    ROCK2/RKC INITIALIZATION  ---------------------------
    bool verbose=true;
    if (print_GNU)
      verbose=false;
    //RKC rock(one_step, verbose, dt_adaptivity, atol, rtol, intrho);
    ROCK2 rock(one_step, verbose, dt_adaptivity, atol, rtol, intrho);    
    if (!print_GNU)
      rock.print_info();


//-------------------------   TIME LOOP   --------------------------------------

    if(rock.check_correctness(dt)==0) //checks if given parameters are ok
        return 0;//something is wrong

    int idid = 2;
    for(int iter=0;idid==2;++iter)
    {
        rock.advance(heat,dt,idid);
    }   
    if (!print_GNU)
      rock.print_statistics();

//---------------------------   WRITE ON SCREEN   ----------------------------------------
    if (print_GNU)
    {
   FILE * branchA;
   FILE * branchB;
   FILE * branchC;
   branchA = fopen ("branchA.txt","w");
   branchB = fopen ("branchB.txt","w");
   branchC = fopen ("branchC.txt","w");
    for (int i=0; i<(mesh.n_L1-1); i++){
      fprintf(branchA,"%.8f \t %.8f \n",mesh.grid[i],heat->yn[i]);
    }
   fclose (branchA);
    for (int i=(mesh.n_L1-1); i<(mesh.n_L1+mesh.n_L2-2); i++){
      fprintf(branchB,"%.8f \t %.8f \n",mesh.grid[i],heat->yn[i]);
    }
   fclose (branchB);
    for (int i=(mesh.n_L1+mesh.n_L2-2); i<mesh.n_elem; i++){
      fprintf(branchC,"%.8f \t %.8f \n",mesh.grid[i],heat->yn[i]);
    }
   fclose (branchC);
    }
//------------------------   DEALLOCATE   -------------------------------------    
//    delete heat;
*/
    return 0;
}

