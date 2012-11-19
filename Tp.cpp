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

#define SINPARAMETROS
#define COVERS

Agustin Brian Federico


static int losCutCallbacks(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p);
void theCovernCuts(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p);
void theGreatGomoryCuts(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p);
void ohMyCliqueCuts(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p);
int greedyForCovern(int fila,CPXCENVptr env,void *cbdata,int wherefrom);
int dpForCovern(int fila, CPXCENV ptr env, void *cbdate, int wherefrom);



double* xestrella;
double* aMochila;
double* rhsMochila;
int* rmatbegMochila;
int cantFilasMochila;
double* rmatvalMochila;
char* senseMochila;
int* rmatindMochila;
vector<int> lasMochila;
double bMochila;
int cantVar;
int noZeroCountMochila;

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
	
	
	
	//De Juanjo
	status = CPXsetintparam(env, CPX_PARAM_PREIND, 0);
	status = CPXsetintparam(env, CPX_PARAM_PRELINEAR, 0);
	status = CPXsetintparam(env, CPX_PARAM_EACHCUTLIM, 0);
    status = CPXsetintparam(env, CPX_PARAM_CUTPASS, 0);
    status = CPXsetintparam(env, CPX_PARAM_FRACCUTS, -1);
    status = CPXsetintparam(env, CPX_PARAM_HEURFREQ, -1);
    status = CPXsetintparam(env, CPX_PARAM_RINSHEUR, -1);
    status = CPXsetintparam(env, CPX_PARAM_REDUCE, 0);
    status = CPXsetintparam(env, CPX_PARAM_IMPLBD, -1);
    status = CPXsetintparam(env, CPX_PARAM_MCFCUTS, -1);
    status = CPXsetintparam(env, CPX_PARAM_ZEROHALFCUTS, -1);
    status = CPXsetintparam(env, CPX_PARAM_MIRCUTS, -1);
    status = CPXsetintparam(env, CPX_PARAM_GUBCOVERS, -1);
    status = CPXsetintparam(env, CPX_PARAM_FLOWPATHS, -1);
    status = CPXsetintparam(env, CPX_PARAM_FLOWCOVERS, -1);
    status = CPXsetintparam(env, CPX_PARAM_DISJCUTS, -1);
    status = CPXsetintparam(env, CPX_PARAM_COVERS, -1);
    status = CPXsetintparam(env, CPX_PARAM_CLIQUES, -1);
    status = CPXsetintparam(env, CPX_PARAM_THREADS, 1);
    status = CPXsetintparam(env, CPX_PARAM_MIPSEARCH, 1);
    status = CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
	
	
	
	
	///////////////////////////////////////////////////////////////////////////////////
	
	
	
	
	
	
	//Lectura de instancia
	//~ CPXreadcopyprob(env,lp,"neos-820146.mps.gz", NULL);
	CPXreadcopyprob(env,lp,"nuestraInstancia.lp", NULL);
	 //~ macrophage.mps.gz
	//~ bley_xl1.mps.gz

	cantVar = CPXgetnumcols(env,lp);
	xestrella = (double*)malloc(sizeof(double)*cantVar);
	aMochila = (double*)malloc(sizeof(double)*cantVar);
	cantFilasMochila = CPXgetnumrows(env,lp);
	int surplus;
	//~ int rmatspaceMochila = cantVar*cantFilasMochila;
	int rmatspaceMochila = 1000000;
	cerr << "llega " <<cantVar<<' ' << cantFilasMochila << ' ' << rmatspaceMochila<< endl;
	rhsMochila = (double*)malloc(sizeof(double)*cantFilasMochila);
	rmatbegMochila = (int*)malloc(sizeof(int)*cantFilasMochila);
	rmatindMochila = (int*)malloc(sizeof(int)*rmatspaceMochila);
	rmatvalMochila = (double*)malloc(sizeof(double)*rmatspaceMochila);
	senseMochila = (char*)  malloc(sizeof(char)*cantFilasMochila);
	
	cerr << rhsMochila << ' ' << rmatbegMochila << ' ' << rmatindMochila << ' '<<rmatvalMochila << endl;

	status = CPXgetsense(env, lp, senseMochila, 0, cantFilasMochila-1);
	if(status){exit(-1);}
	status = CPXgetrhs(env, lp, rhsMochila, 0, cantFilasMochila-1);
	if(status){exit(-1);}
	status = CPXgetrows(env, lp, &noZeroCountMochila, rmatbegMochila, rmatindMochila, rmatvalMochila, rmatspaceMochila, &surplus, 0, cantFilasMochila-1);
	if(surplus<0){cerr << "No teniamos espacio para traer las filas" << endl; exit(-1);}
	if(status){exit(-1);}
	
	
	
	//Chequeamos todas las filas y indicamos cuales son mochila
	for(int i=0;i<cantFilasMochila;++i){
		if(senseMochila[i]!='L') continue;
		if(rhsMochila[i]<0)continue;
		if(fabs(rhsMochila[i]-floor(rhsMochila[i]))>10e-6)continue;
		bool laMeto = true;
		for(int j=rmatbegMochila[i];j<(i==cantFilasMochila-1?noZeroCountMochila:rmatbegMochila[i+1]);++j){
			if(rmatvalMochila[j]<0)laMeto = false;
			if(fabs(rmatvalMochila[j]-floor(rmatvalMochila[j]))>10e-6)laMeto = false;
		}
		if(laMeto)lasMochila.push_back(i);
	}
	
	cerr << endl << endl << endl << "Total de filas: " << cantFilasMochila << endl << "Total mochila: " << lasMochila.size() << endl << endl;
	
	
	
	
	
	//~ for(int i=0;i<cantFilas;++i){
		//~ cerr << rhs[i] << endl;
	//~ }
	
	
	
	
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
	
	status = CPXgetx(env, lp, xestrella, 0, cantVar-1);
	for(int i=0;i<cantVar;++i){
		cerr << xestrella[i] << endl;
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
	
	//Variable para ver que heuristica o exacto corremos
	int tengoCover = 0;
	
	
	for(int i=0;i<(int)lasMochila.size();++i){
		tengoCover += greedyForCovern(lasMochila[i], env, cbdata, wherefrom);
	}
	
	for(int i=0;!tengoCover&&i<(int)lasMochila.size();++i){
		tengoCover += dpForCovern(lasMochila[i], env, cbdata, wherefrom);
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


int dpForCovern(int fila, CPXCENV ptr env, void *cbdate, int wherefrom){
	int tengoCover = 0;
	bMochila = rhsMochila[fila]++;
	if(bmochila>1000){
		cerr << "Muy grande el b como para hacer dp" << endl;
		return 0;
	}
	
	double dpMochila[2][bMochila+5];
	int dpPrevMochila[(fila==cantFilasMochila-1?noZeroCountMochila:rmatbegMochila[fila+1])-rmatbegMochila[fila]][bMochila+5];
	int act = 0;
	dpMochila[0][0] = 0;
	dpPrevMochila[0][0] = 0;
	for(int j=1;j<bMochila;++j){
		if(j<=rmatvalMochila[rmatbegMochila[fila]]){
			dpMochila[0][j] = (1-xestrella[rmatindMochila[rmatbegMochila[fila]]]);
			dpPrevMochila[0][j] = 1;
		}else{
			dpMochila[0][j] = 10e100;
			dpPrevMochila[0][j] = 0;
		}
	}
	for(int i=rmatbegMochila[fila]+1;i<(fila==cantFilasMochila-1?noZeroCountMochila:rmatbegMochila[fila+1]);++i){
		for(int j=0;j<bMochila;++j){
			//~ dpMochila[act%2][j] = min((1-xestrella[rmatindMochila[i]])+dpMochila[(act-1)%2][max(0,j-(rmatvalMochila[i]))],dpMochila[(act-1)%2][j]);
			if((1-xestrella[rmatindMochila[i]])+dpMochila[(act-1)%2][max(0,j-(rmatvalMochila[i]))] < dpMochila[(act-1)%2][j]){
				dpMochila[act%2][j] = (1-xestrella[rmatindMochila[i]])+dpMochila[(act-1)%2][max(0,j-(rmatvalMochila[i]))];
				dpPrevMochila[i-rmatbegMochila[fila]][j] = 1;
			}else{
				dpMochila[act%2][j] = dpMochila[(act-1)%2][j];
				dpPrevMochila[i-rmatbegMochila[fila]][j] = 0;
			}
		}
		act++;
	}
	
	
	
	if(dpMochila[(fila==cantFilasMochila-1?noZeroCountMochila:rmatbegMochila[fila+1])-rmatbegMochila[fila]-1][bMochila] < 0.99){
		cerr << "Cover violado por DP" << endl;
		//~ exit(-1);
		
		//Hay que recuperar la solucion y ahi generar el corte
		
		
		double rhs = i-1.;
		//~ cerr << "Dios mio alto corte: " << endl << endl;
		//~ cerr << rhs << endl;
		int noZeroCount = i;
		
		int *cut_matind=(int*) calloc(noZeroCount,sizeof(int));
		double *cut_matval=(double*) calloc(noZeroCount,sizeof(double));
		for(int j=0;j<i;++j){
			cut_matind[j] = rmatindMochila[vec[j].second];
			cut_matval[j] = 1.0;
			cerr << rmatindMochila[vec[j].second] << endl;
		}
		CPXcutcallbackadd(env, cbdata, wherefrom, noZeroCount, rhs, 'L', cut_matind, cut_matval, 0);
		tengoCover++;
	}
	
	return tengoCover;
}

int greedyForCovern(int fila,CPXCENVptr env,void *cbdata,int wherefrom ){
	int tengoCover = 0;
	
	
	bMochila = rhsMochila[fila];
	
	
	//Greedy 1 (A_j/(1-x_j)
	vector<pair<double, int> > vec;
	for(int i=rmatbegMochila[fila];i<(fila==cantFilasMochila-1?noZeroCountMochila:rmatbegMochila[fila+1]);++i){
		vec.push_back(pair<double,int>(rmatvalMochila[i]/(1-xestrella[rmatindMochila[i]]),i));
		cerr << "Toda la info: " << "el i " << i << endl << "el valor: " << (rmatvalMochila[i]/(1-xestrella[rmatindMochila[i]])) << endl;
		cerr << "rmatvalMochila: " << rmatvalMochila[i] << " lo otro " << 1-xestrella[rmatindMochila[i]] << endl;
	}
	sort(vec.begin(),vec.end(), greater<pair<double,int> >());
	
	double acumA=0.0, acumPesos=0.0;
	
	//Ahora sumamos las a_j hasta pasarnos del bMochila+1
	int i;
	for(i=0;i<(int)vec.size()&&acumA<bMochila+0.5;++i){
		acumA += rmatvalMochila[vec[i].second];
		acumPesos += (1.-xestrella[rmatindMochila[vec[i].second]]);
	}
	
	//Preguntamos si los pesos suman menos que 1 (si es asi, hay un corte violado)
	if(acumPesos < .99 && acumA>bMochila+0.5){
		cerr << "Cover violado" << endl;
		//~ exit(-1);
		double rhs = i-1.;
		cerr << "Dios mio alto corte: " << endl << endl;
		cerr << rhs << endl;
		int noZeroCount = i;
		
		int *cut_matind=(int*) calloc(noZeroCount,sizeof(int));
		double *cut_matval=(double*) calloc(noZeroCount,sizeof(double));
		for(int j=0;j<i;++j){
			cut_matind[j] = rmatindMochila[vec[j].second];
			cut_matval[j] = 1.0;
			cerr << rmatindMochila[vec[j].second] << endl;
		}
		CPXcutcallbackadd(env, cbdata, wherefrom, noZeroCount, rhs, 'L', cut_matind, cut_matval, 0);
		tengoCover++;
	}
	return tengoCover;
}
