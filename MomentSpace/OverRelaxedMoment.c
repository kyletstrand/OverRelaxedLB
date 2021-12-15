//OverRelaxedLB.c
//Looking at OverRelaxation in diffusive LB to compare to IntegerLG version.

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <mygraph.h>

#define xdim 320
#define V 3

double f[V][xdim];
double feq[V][xdim];
double rho[xdim];
double w[V];
double tau[V] = {1., 1., 1.};
double n0 = 100;
double omega = 1.;
double theta = 1./3.;

int next = 0;
int Pause = 1;
int done = 0;
int repeat = 1;
int iterations;

void SetWeights() {
	w[0] = 1. - theta;
	w[1] = w[2] = theta/2.;
}	

void SetEqDist() {
	for (int a = 0; a < V; a++) {
		for (int i = 0; i < xdim; i++) {
			feq[a][i] = rho[i] * w[a];
		}
	}
}

void Initialize() {
	
	iterations = 0;
	SetWeights();
	for (int i = 0; i < xdim; i++) {
		rho[i] = n0 * (1. + sin(2.*M_PI * i/xdim)); //sine wave
		for (int a = 0; a < V; a++) {
			f[a][i] = rho[i] * w[a];
		}
	}
	SetEqDist();
}

void Stream() {
	double tmp1 = f[1][xdim-1];
	double tmp2 = f[2][0];
	memmove(&f[1][1], &f[1][0], (xdim-1)*sizeof(double));
	memmove(&f[2][0], &f[2][1], (xdim-1)*sizeof(double));
	f[1][0] = tmp1;
	f[2][xdim-1] = tmp2;
}

void Iteration() {

	double M[V][xdim];
	double Meq[V][xdim];
	//double Lambda[3] = {1., 2. - omega, 2. - omega}; // It works with this as standard. 
	double Lambda[3] = {1., omega - 2., 2. - omega}; // This fails
	double m[3][3] = {   //This might be wrong...double check
		{1., 1., 1.},
		{0., sqrt(1./theta), -sqrt(1./theta)},
		{-sqrt(theta/(1.-theta)), sqrt((1.-theta)/theta), sqrt((1.-theta)/theta)}
	};

	for (int i = 0; i < xdim; i++) {
		rho[i] = f[0][i] + f[1][i] + f[2][i];
	}
	SetEqDist();
	for (int i = 0; i < xdim; i++) {	
		for (int v = 0; v < V; v++) {
  	  M[v][i] = 0;
			Meq[v][i] = 0;
  	}

		//forward xform
		for (int a = 0; a < V; a++) {
			for (int b = 0; b < V; b++) {
				M[a][i] += m[a][b] * f[b][i];
				Meq[a][i] += m[a][b] * feq[b][i];
			}
		}

		for (int a = 0; a < V; a++) {
			//M[a][i] = (1. - Lambda[a]) * M[a][i];  //This could also be wrong. Verify! Clear over-relaxation at tau=0.6. what is that in omega?
			//M[a][i] = 1./tau[a] * (Meq[a][i] * (tau[a] - 1.) * M[a][i]);
			M[a][i] = Lambda[a] * (Meq[a][i] + (1./Lambda[a] - 1.) * M[a][i]);
		}

		f[0][i] = f[1][i] = f[2][i] = 0;

		//backward xform
		for (int a = 0; a < V; a++) {
			for (int b = 0; b < V; b++) {
				f[a][i] += m[b][a] * M[b][i] * w[a];
			}
		}	
	}

	Stream();
	iterations++;
}

void GUI() {
	static int Xdim = xdim;

	DefineGraphN_R("rho", &rho[0], &Xdim, NULL);
	StartMenu("D1Q3 Moment Space Diffusion", 1);
		DefineInt("Iterations", &iterations);
		DefineDouble("Omega", &omega);
		DefineFunction("Initialize", &Initialize);
		DefineGraph(curve2d_, "Graphs");
		DefineInt("Repeat", &repeat);
		DefineBool("Next", &next);
		DefineBool("Pause", &Pause);
		DefineBool("Quit", &done);
	EndMenu();
}

int main(int argc, char *argv[]) {
	int newdata = 1;
	int i;

	Initialize();
	GUI();

	while (done == 0) {
		Events(newdata);
		DrawGraphs();
		if (next || !Pause) {
			newdata = 1;
			next = 0;
			for (int i = 0; i < repeat; i++) {
				Iteration();
			}
		}
		else sleep(1);
	}

	return 0;
}
