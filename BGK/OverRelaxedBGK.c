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

void Collision() {
	for (int a = 0; a < V; a++) {
		for (int i = 0; i < xdim; i++) {
			f[a][i] += (2. - omega) * (feq[a][i] - f[a][i]);
		}
	}
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

	for (int i = 0; i < xdim; i++) {
		rho[i] = f[0][i] + f[1][i] + f[2][i];
	}
	SetEqDist();
	Collision();	
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
