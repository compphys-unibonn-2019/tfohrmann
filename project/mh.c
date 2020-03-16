/*  
 *  Metropolis-Hastings-Algorithm for scalar field theory using the parameters chosen below.
 *  This program will store the lattice as .txt-files.
 *  Compile using: gcc -Wall -pedantic -std=c99 mh.c -o mh -lgsl -lm
 *  Execute with: ./mh
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <sys/time.h>

double delta_s(double**** lattice, int T, int N, double change, int pos[4], double mass, double field);
double sweep(double**** lattice, int T, int N, double mass, double field, double w);
void save_to_txt(double**** lattice, int T, int N, char* name);

int main(){
    // Simulation parameters
    const double mass = 1.0;
    const double field = 0.0;
    
    const int N = 16; //size of lattice in "time-direction"
    const int T = 32; //size of lattice in "space-direction"
    
    const double w = 1.0; //width of (uniform) proposal density
    const int n_sweeps = 30000; //number of sweeps to be calculated in total
    
    
    // Random generator initilisation
    const gsl_rng_type * type;
    gsl_rng * r;
    gsl_rng_env_setup();
    type = gsl_rng_ranlxd2;
    r = gsl_rng_alloc(type);
    
    struct timeval tv;
    gettimeofday(&tv, 0);
    unsigned long seed = tv.tv_sec + tv.tv_usec;
    gsl_rng_set(r, seed);
    
    
    
    // Set up simulation
    double AR = 0.0; //acceptance rate

    // Set up lattice:
    double**** lattice = (double****)malloc(T * sizeof(double***));
    for(int i=0;  i<T;  i++){
		lattice[i] = (double***)malloc(N * sizeof(double**));
        
		for(int j=0;  j<N;  j++){
			lattice[i][j] = (double**)malloc(N * sizeof(double*));
            
            for(int k=0; k<N; k++){
                lattice[i][j][k] = (double*)malloc(N * sizeof(double));
            }
		}
	}
    
    
    // Initilaise hot start:
    double total = 0;
    for(int i=0; i<T;  i++)
	{
		for(int j=0;  j<N;  j++)
		{
			for(int k=0;  k<N;  k++)
			{
                for(int l=0; l<N; l++)
                {
                    lattice[i][j][k][l] = 2.0*gsl_rng_uniform(r)-1.0;
                    total += lattice[i][j][k][l];
                }	
			}
		}
	}
    for(int i=0; i<T;  i++)
	{
		for(int j=0;  j<N;  j++)
		{
			for(int k=0;  k<N;  k++)
			{
                for(int l=0; l<N; l++)
                {
                    //lattice[i][j][k][l] /= total;
                    lattice[i][j][k][l] = 0;
                }	
			}
		}
	}
    gsl_rng_free(r);
    
    /*
    // Thermalise
    for(int i=0; i<100; i++){
        AR += sweep(lattice, T, N, mass, field, w);
    }
    printf("Thermalised!\n");
    */
    // Perform n_sweeps sweeps of the lattice
    char name[30];
    for(int i=0; i<n_sweeps; i++){
        AR += sweep(lattice, T, N, mass, field, w);
        if(i%100==0){
            sprintf(name, "./m1f0/sweep_%d.txt", i/100);
            save_to_txt(lattice, T, N, name);
        }
        
        // Fortschrittsanzeige
		printf("\b\b\b\b%3.0lf%%", ((double)i)/((double)n_sweeps)*100.0);
		fflush(stdout);
    }
    AR /= (double)(n_sweeps+100);
    printf("\n--- Done ---\n");
    printf("AR: %f\n", AR);
    
    
    
    return 0;
}


double delta_s(double**** lattice, int T, int N, double change, int pos[4], double mass, double field){
    // Calculates the difference in the action given a proposed change to the lattice at position pos.
    double delta = 0;
    double sp = lattice[pos[0]][pos[1]][pos[2]][pos[3]] + change;
    
    // After change:
    // Time
    delta += pow(sp - lattice[(pos[0]-1+T)%T][pos[1]][pos[2]][pos[3]], 2) /2.0; //before
    delta += pow(lattice[(pos[0]+1)%T][pos[1]][pos[2]][pos[3]] - sp, 2) /2.0; //after
    
    // x-Dim
    delta += pow(sp - lattice[pos[0]][(pos[1]-1+N)%N][pos[2]][pos[3]], 2) /2.0; //before
    delta += pow(lattice[pos[0]][(pos[1]+1)%N][pos[2]][pos[3]] - sp, 2) /2.0; //after
    
    // y-Dim
    delta += pow(sp - lattice[pos[0]][pos[1]][(pos[2]-1+N)%N][pos[3]], 2) /2.0; //before
    delta += pow(lattice[pos[0]][pos[1]][(pos[2]+1)%N][pos[3]] - sp, 2) /2.0; //after
    
    // z-Dim
    delta += pow(sp - lattice[pos[0]][pos[1]][pos[2]][(pos[3]-1+N)%N], 2) /2.0; //before
    delta += pow(lattice[pos[0]][pos[1]][pos[2]][(pos[3]+1)%N] - sp, 2) /2.0; //after
    
    // Pure Field
    delta += -pow(mass * sp, 2)/2.0 - field/24.0 * pow(sp, 4);
    
    
    // Before change:
    sp -= change;
    // Time
    delta -= pow(sp - lattice[(pos[0]-1+T)%T][pos[1]][pos[2]][pos[3]], 2) /2.0; //before
    delta -= pow(lattice[(pos[0]+1)%T][pos[1]][pos[2]][pos[3]] - sp, 2) /2.0; //after
    
    // x-Dim
    delta -= pow(sp - lattice[pos[0]][(pos[1]-1+N)%N][pos[2]][pos[3]], 2) /2.0; //before
    delta -= pow(lattice[pos[0]][(pos[1]+1)%N][pos[2]][pos[3]] - sp, 2) /2.0; //after
    
    // y-Dim
    delta -= pow(sp - lattice[pos[0]][pos[1]][(pos[2]-1+N)%N][pos[3]], 2) /2.0; //before
    delta -= pow(lattice[pos[0]][pos[1]][(pos[2]+1)%N][pos[3]] - sp, 2) /2.0; //after
    
    // z-Dim
    delta -= pow(sp - lattice[pos[0]][pos[1]][pos[2]][(pos[3]-1+N)%N], 2) /2.0; //before
    delta -= pow(lattice[pos[0]][pos[1]][pos[2]][(pos[3]+1)%N] - sp, 2) /2.0; //after
    
    // Pure Field
    delta -= -pow(mass * sp, 2)/2.0 - field/24.0 * pow(sp, 4);
    
    
    return delta;
}

double sweep(double**** lattice, int T, int N, double mass, double field, double w){
    // Sequentially makes update to the entire lattice.
    // Random generator init
    const gsl_rng_type * type;
    gsl_rng * s;
    gsl_rng_env_setup();
    type = gsl_rng_ranlxd2;
    s = gsl_rng_alloc(type);
    
    struct timeval tv;
    gettimeofday(&tv, 0);
    unsigned long seed = tv.tv_sec + tv.tv_usec;
    gsl_rng_set(s, seed);
    
    
    // Loop through whole lattice
    double AR = 0.0;
    double proposal;
    double del_S;
    int pos[4];
    double norm = 1.0/pow(mass, 4);
    for(int i=0; i<T;  i++)
	{
        pos[0] = i;
		for(int j=0;  j<N;  j++)
		{
            pos[1] = j;
			for(int k=0;  k<N;  k++)
			{
                pos[2] = k;
                for(int l=0; l<N; l++)
                {
                    pos[3] = l;
                    proposal = w * gsl_rng_uniform(s) - w/2.0;
                    del_S = delta_s(lattice, T, N, proposal, pos, mass, field);
                    
                    if(gsl_rng_uniform(s) < exp(-del_S*norm)){
                        lattice[i][j][k][l] += proposal;
                        AR += 1.0;
                    }
                }	
			}
		}
    }
    
    gsl_rng_free(s);
    return(AR/(T*pow(N,3)));
}

void save_to_txt(double**** lattice, int T, int N, char* name){
    // Saves a lattice configuration as a series of field points.
    // Open file
	FILE* output;
	output = fopen(name, "w");
	if(output == NULL){
		printf("Couldn't create/open '%s'\n", name);
		exit(129);
	}
    
    for(int i=0; i<T;  i++)
	{
		for(int j=0;  j<N;  j++)
		{
			for(int k=0;  k<N;  k++)
			{
                for(int l=0; l<N; l++)
                {
                    fprintf(output, "%f,", lattice[i][j][k][l]);
                }	
			}
		}
	}
    fprintf(output, "373.0"); //this is just to not have a "," at the end

    // Close file
	if(fclose(output)){
		printf("'%s' couldn't be closed!\n", name);
		printf("Exiting program.\n");
		exit(130);
	}
    /*
	else{
		printf("'%s' written\n", name);
	}
    */
}
