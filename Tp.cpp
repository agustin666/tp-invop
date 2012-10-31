#define Agustin using
#define Brian namespace
#define Federico std;
#include <ilcplex/cplex.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <cmath>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>


Agustin Brian Federico


int main(){
	double objval;
	
	
	CPXENVptr env = NULL;
	CPXLPptr lp = NULL;
	int status;
	
	
	env = CPXopenCPLEX(&status);
	
	if(env == NULL){
		char errmsg[CPXMESSAGEBUFSIZE];
		fprintf(stderr, "No se pudo inicializar el ambiente.\n");
		CPXgeterrorstring(env, status, errmsg);
		fprintf(stderr, "%s", errmsg);
		exit(-1);
	}
	
	lp = CPXcreateprob(env, &status, "Test");
	
	if(lp==NULL){
		fprintf(stderr, "No se pudo crear el problema, error %d\n", status);
		exit(-1);
	}
	
	
	//Para habilitar el indicador por pantalla
	status = CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
	if(status){
		fprintf(stderr,"Fallo el indicador por pantalla, error %d\n", status);
		exit(-1);
	}
	
	status = CPXsetintparam(env, CPX_PARAM_DATACHECK, CPX_ON);
	if(status){
		fprintf(stderr,"No se pudo habilitar data checking, error %d\n", status);
		exit(-1);
	}
	
	
	CPXreadcopyprob(env,lp,"cov1075.mps.gz", NULL);
	status = CPXmipopt(env,lp);
	if(status){cerr << "Fallo" << endl;exit(-1);}
	
	status = CPXgetobjval (env, lp, &objval);
	if ( status ) {
		fprintf (stderr,"No MIP objective value available.  Exiting...\n");
		exit(-1);
	}
	
	printf ("Valor de la solucion  = %f\n\n", objval);
	
	
	return 0;
}




