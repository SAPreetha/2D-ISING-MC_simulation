//
//  Preetha_Pr1_wolf.cpp
//  
//
//  Created by Preetha Saha on 3/15/17.
//
//


#include <stdio.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <list>
#include <tuple>

using namespace std;

const double J = +1; // ferromagnetic interaction in nearest neighbours

const int N=4;
const int Ns=N*N;// number of spins in 2D square lattice
int **s;   // the spins
double T;  // temperature

//int steps;
const int MCSteps=500000; // No of monte carlo steps
const int EQSteps=5000; // No of monte carlo steps

void initialize_lattice( ){
    s = new int* [N];
    for (int i = 0; i < N; i++)
        s[i] = new int [N];
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++){
            s[i][j] = (float)rand()/RAND_MAX < 0.5 ? +1 : -1;
            printf("%d\n",s[i][j]);
        }
    
}

bool **cluster;    // array of boolean variables cluster[i][j] = true if i,j belongs
double addProb;    // 1 - e^(-2J/T)

void initializeClusterVariables(double T) {
    
    
    cluster = new bool* [N];
    for (int i = 0; i < N; i++)
        cluster[i] = new bool [N];
    
    // probability to add a like spin to the cluster
    addProb = 1 - exp(-2*(J/T));
    
    
}


// declare functions to implement Wolff algorithm
void growCluster(int i, int j, int clusterSpin,double T);
void tryAdd(int i, int j, int clusterSpin,double T);



void oneMonteCarlo(double T) {
    
    // All cluster assigned to false to indicate starting lattice has no spin in cluster
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            cluster[i][j] = false;
    
    // choosing a random spin and grow a cluster
    int i = int(rand()%N);
    int j = int(rand()%N);
    growCluster(i, j, s[i][j],T);
    
    //++steps;
}

void growCluster(int i, int j, int clusterSpin,double T) {
    
    // starting at (i.j) spin as belonging to the cluster(variable set to true), and spin value is flipped
    cluster[i][j] = true;
    s[i][j] = -s[i][j];
    
    // 4 nearest neighbors of each (i,j)spin and implementing periodic boundary conditions
    
    int iPrev = i == 0    ? N-1 : i-1;
    int iNext = i == N-1 ? 0    : i+1;
    int jPrev = j == 0    ? N-1 : j-1;
    int jNext = j == N-1 ? 0    : j+1;
    
    // if the neighbor spin does not belong to the cluster we try to add it to the cluster
    
    if (!cluster[iPrev][j]){
        tryAdd(iPrev, j, clusterSpin,T);
        
    }
    if (!cluster[iNext][j]){
        tryAdd(iNext, j, clusterSpin,T);
        
    }
    if (!cluster[i][jPrev])
        tryAdd(i, jPrev, clusterSpin,T);
    if (!cluster[i][jNext])
        tryAdd(i, jNext, clusterSpin,T);
}

//tryAdd and growcluster are recursively applied on each spin which is successfully becomes a part of thr lattice

void tryAdd(int i, int j, int clusterSpin,double T) {
    if (s[i][j] == clusterSpin)
        if ((double)rand()/RAND_MAX < addProb)
            growCluster(i, j, clusterSpin,T);
}


// Calculating Magnetization of the lattice
float Magnetisation(){
    double Mag=0;
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            Mag =Mag + s[i][j];          }
    }
    return Mag;
    
}

// Calculating Energy of the lattice
float Energy(){
    float Ene=0;
    //float energy;
    int i;
    int j;
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            Ene =Ene + (-J*s[i][j]*(s[(i+1)%N][j]+s[i][(j+1)%N]));
        }
        
    }
    //printf("Ene is %f",Ene);
    return Ene;
    
    
}

// Getting the system in Equilibrium before taking measurements
void Thermalise(double T){
    for(int i=0;i<EQSteps;i++){
        oneMonteCarlo(T);
    }
}

//Looping over Teperature
void Measure_observables(){
    FILE *g = fopen("wolf.txt", "w+");
    
    for(double T=.01;T<6;T+=0.04){
        initializeClusterVariables(T);
        Thermalise(T);
        
        double Enesum=0;
        double Magsum=0;
        double Ene_sq_sum=0;
        double count=0;
        
        for(int i=0;i<MCSteps;i++){
            oneMonteCarlo(T);
            if(i%100==0){
                count +=1;
                
                double Ene =  Energy();
                double m = Magnetisation();
                
                
                Enesum = Enesum+Ene;
                Magsum = Magsum+m;
                Ene_sq_sum += Ene*Ene;
                
                
                
            }
        }
        double E_avg=Enesum/(count*Ns);
        double Cv_avg=((Ene_sq_sum/(count)-((Enesum/count)*(Enesum/count)))/(T*T))/Ns;
        //printf("Final Enesum is %lf\n",Enesum/(count*Ns));
        fprintf(g,"%f\t%f\t%f\n",T,E_avg,Cv_avg);

    }
}

int main() {
    
    double Enesum;
    double count;
    srand(time(NULL));
    
    
    
    initialize_lattice();
    Measure_observables();
    return 0;
}










