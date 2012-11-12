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
#include <algorithm>


Agustin Brian Federico


static int losCutCallbacks(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p);
void theCovernCuts(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p);
void theGreatGomoryCuts(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p);
void ohMyCliqueCuts(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p);


double* xestrella;
double* aMochila;
double bMochila;
int cantVar;

int main(){
	double objval;
	
	
	
	////////////////////////////////////////////////////////////////////
	//Arrancamos CPLEX
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
	
	///////////////////////////////////////////////////////////////////////////////////
	//Parametros de cplex
	//
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
	
	//Con este parametro en off, siempre no estamos refiriendo al problema original
	status = CPXsetintparam (env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
   	if ( status )  exit(-1);
	
	
	/* Turn on traditional search for use with control callbacks */
	//~ status = CPXsetintparam (env, CPX_PARAM_MIPSEARCH, CPX_MIPSEARCH_TRADITIONAL);
   	//~ if ( status ) exit(-1);
	
	
	///////////////////////////////////////////////////////////////////////////////////
	
	
	
	
	
	
	//Lectura de instancia
	CPXreadcopyprob(env,lp,"cov1075.mps.gz", NULL);
	
	cantVar = CPXgetnumcols(env,lp);
	xestrella = (double*)malloc(sizeof(double)*cantVar);
	aMochila = (double*)malloc(sizeof(double)*cantVar);
	
	
	
	
	////////////////////////////////////////////////////////////////////
	//Optimiçao na pesquisa operacional
	
	//Seteamos nuestra funcion de callback
	status = CPXsetusercutcallbackfunc(env, losCutCallbacks, NULL);
	
	
	
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


///////////////////////////////////////
static int losCutCallbacks (CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p){
	
	theCovernCuts(env, cbdata, wherefrom, cbhandle, useraction_p);
	//~ theGreatGomoryCuts(env, cbdata, wherefrom, cbhandle, useraction_p);
	//~ ohMyCliqueCuts(env, cbdata, wherefrom, cbhandle, useraction_p);
	return 0;
} 

void theCovernCuts(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p){
	*useraction_p = CPX_CALLBACK_DEFAULT;
	
	int status;
	
	//Traemos la solucion de la relajacion
	status = CPXgetcallbacknodex(env,cbdata,wherefrom,xestrella,0,cantVar-1);
	if(status){cerr << "no pude traer x actual" << endl;exit(-1);}
	
	//Aca tenemos las aMochila
	//~ memcpy(aMochila, nosededonde, sizeof aMochila);
	//~ memset(aMochila, 63, sizeof aMochila);
	//Aca tenemos la bMochila
	//~ bMochila = vemos de donde te traemos;
	//~ bMochila = -100;
	
	
	
	//Variable para ver que heuristica o exacto corremos
	int tengoCover = 0;
	
	
	//Greedy 1 (A_j/(1-x_j)
	vector<pair<double, int> > vec;
	for(int i=0;i<cantVar;++i){
		vec.push_back(pair<double,int>(aMochila[i]/(1-xestrella[i]),i));
	}
	sort(vec.begin(),vec.end(), greater<pair<double,int> >());
	
	double acumA=0.0, acumPesos=0.0;
	
	//Ahora sumamos las a_j hasta pasarnos del bMochila+1
	int i;
	for(i=0;i<cantVar&&acumA<bMochila+0.5;++i){
		acumA += aMochila[vec[i].second];
		acumPesos += (1.-xestrella[vec[i].second]);
	}
	
	//Preguntamos si los pesos suman menos que 1 (si es asi, hay un corte violado)
	if(acumPesos < 1.){
		cerr << "Cover violado" << endl;
		double rhs = i-1.;
		int noZeroCount = i;
		
		int *cut_matind=(int*) calloc(noZeroCount,sizeof(int));
		double *cut_matval=(double*) calloc(noZeroCount,sizeof(double));
		for(int j=0;j<i;++j){
			cut_matind[j] = vec[j].second;
			cut_matval[j] = 1.0;
		}
		status = CPXcutcallbackadd(env, cbdata, wherefrom, noZeroCount, rhs, 'L', cut_matind, cut_matval, 0);
		
		tengoCover++;
		
	}
	
	
	
	//~ cerr << "I tried to dis covern!!!" << endl;
	
	if(tengoCover)*useraction_p = CPX_CALLBACK_SET;
	return;
}


void theGreatGomoryCuts(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p){
	*useraction_p = CPX_CALLBACK_DEFAULT;
	cerr << "a little something to make me gomory" << endl;


	return;
}


void ohMyCliqueCuts(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p){
	*useraction_p = CPX_CALLBACK_DEFAULT;
	cerr << "oh baby refrain, cuz you are breaking my clique" << endl;
	

	return;
}

