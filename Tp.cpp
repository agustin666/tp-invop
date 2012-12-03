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
#include "input.h"

Agustin Brian Federico 

#define COVERGREEDY
//#define COVERDP

//#define INSTANCIA "./miplib2010-benchmark/bab5.mps.gz"
//#define INSTANCIA "./miplib2010-benchmark/neos-1337307.mps.gz"
//#define INSTANCIA "./miplib2010-benchmark/neos-1396125.mps.gz"
//#define INSTANCIA "./miplib2010-benchmark/reblock67.mps.gz"
//#define INSTANCIA "./miplib2010-benchmark/ns1830653.mps.gz"
#define INSTANCIA "reblock67.mps.gz"
#define JOYA if(status)errorHandler(status,env)
#define ULTIMO(fila) (fila==cantFilasMochila-1?noZeroCountMochila:rmatbegMochila[fila+1])

static int losCutCallbacks(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p);
void theCovernCuts(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p);
void theGreatGomoryCuts(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p);
void ohMyCliqueCuts(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p);
int greedyForCovern(int fila,CPXCENVptr env,void *cbdata,int wherefrom);
int dpForCovern(int fila, CPXCENVptr env, void *cbdate, int wherefrom);



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
int totCortesGreedy=0;


char errmsg[CPXMESSAGEBUFSIZE];
void errorHandler(int status, CPXENVptr env){
    CPXgeterrorstring(env, status, errmsg);
    fprintf(stderr, "%s", errmsg);
    exit(-1);
}

int* last,*next,*adj;
int ejes; 
int* lastSepClique,*nextSepClique,*adjSepClique;
int ejesSepClique; 

void add_edge(int u,int v){
    next[ejes]=last[u]; adj[ejes]=v; last[u]=ejes++;
    swap(u,v);
    next[ejes]=last[u]; adj[ejes]=v; last[u]=ejes++;
}

int main(int argc, char *argv[]){
  
  string archivo;
  if(isParam("-f", argv, argc)){
    archivo = getParam("-f", argv, argc);
  }else{
    archivo = INSTANCIA;
  }
  
  
  ////////////////////////////////////////////////////////////////////
  //Arrancamos CPLEX
  CPXENVptr env = NULL;
  CPXLPptr lp = NULL;
  int status;
  
  
  env = CPXopenCPLEX(&status);
  
  if(env == NULL){
    fprintf(stderr, "No se pudo inicializar el ambiente.\n");
    errorHandler(status,env);
  }
  
  lp = CPXcreateprob(env, &status, "Test");
  
  if(lp==NULL){
    fprintf(stderr, "No se pudo crear el problema\n", status);
    errorHandler(status,env);
  }
  
  ///////////////////////////////////////////////////////////////////////////////////
  //Parametros de cplex
  //
  status = CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
  if(status){
    fprintf(stderr,"Fallo el indicador por pantalla\n", status);
    errorHandler(status,env);
  }
  
  status = CPXsetintparam(env, CPX_PARAM_DATACHECK, CPX_ON);
  if(status){
    fprintf(stderr,"No se pudo habilitar data checking\n", status);
    errorHandler(status,env);
  }
  
  //Con este parametro en off, siempre no estamos refiriendo al problema original
  status = CPXsetintparam (env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
  JOYA;
  
  
  /* Turn on traditional search for use with control callbacks */
  //~ status = CPXsetintparam (env, CPX_PARAM_MIPSEARCH, CPX_MIPSEARCH_TRADITIONAL);
  //~ if ( status ) exit(-1);
  
  
  if(isParam("-c", argv, argc)){
    cerr << " USANDO EL VERDADERO PODER DE CPLEX " << endl;
  }else{
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
  }
  
  //Lectura de instancia
  CPXreadcopyprob(env,lp,archivo.c_str(), NULL);

  cantVar = CPXgetnumcols(env,lp);
  xestrella = (double*)malloc(sizeof(double)*cantVar);
  aMochila = (double*)malloc(sizeof(double)*cantVar);
  cantFilasMochila = CPXgetnumrows(env,lp);
  int surplus;
  int rmatspaceMochila = min( cantVar*cantFilasMochila , 1000000);
  rhsMochila = (double*)malloc(sizeof(double)*cantFilasMochila);
  rmatbegMochila = (int*)malloc(sizeof(int)*cantFilasMochila);
  rmatindMochila = (int*)malloc(sizeof(int)*rmatspaceMochila);
  rmatvalMochila = (double*)malloc(sizeof(double)*rmatspaceMochila);
  senseMochila = (char*)  malloc(sizeof(char)*cantFilasMochila);
  status = CPXgetsense(env, lp, senseMochila, 0, cantFilasMochila-1);
  JOYA;
  status = CPXgetrhs(env, lp, rhsMochila, 0, cantFilasMochila-1);
  JOYA;
  status = CPXgetrows(env, lp, &noZeroCountMochila, rmatbegMochila, rmatindMochila, rmatvalMochila, rmatspaceMochila, &surplus, 0, cantFilasMochila-1);
  if(surplus<0){cerr << "No teniamos espacio para traer las filas" << endl; errorHandler(status,env);}
  JOYA;

  
  //Chequeamos todas las filas y indicamos cuales son mochila. Excluimos cover y cambiamos signo
  int covers=0;
  int cantG=0;
  // Inicializamos el grafo
  cerr << cantFilasMochila << endl;
  last=(int*)malloc(2*cantVar*sizeof(int));
  adj=(int*)malloc(min(1000000,4*cantVar*cantVar)*sizeof(int));
  next=(int*)malloc(min(1000000,4*cantVar*cantVar)*sizeof(int));
  lastSepClique=(int*)malloc(2*cantVar*sizeof(int));
  adjSepClique=(int*)malloc(min(1000000,4*cantVar*cantVar)*sizeof(int));
  nextSepClique=(int*)malloc(min(1000000,4*cantVar*cantVar)*sizeof(int));
  int cantFilasClique=cantFilasMochila;
  ejes=0;
  memset(last,-1,sizeof last);
  for(int i=0;i<cantVar;++i){
      for(int j=i+1;j<cantVar;++j){
		bool serompe;
		for(int a=0;a<4;++a){
		  int probi=a&1;
		  int probj=a&2;
		  serompe=false;
		  
		  for(int k=0;!serompe && k<cantFilasMochila;++k){
			  if(senseMochila[i]=='E') continue;
			  int cota=0;
			  
			  for(int l=rmatbegMochila[k];l<(k==cantFilasMochila-1?noZeroCountMochila:rmatbegMochila[k+1]);++l){
				if(rmatindMochila[l]!=i && rmatindMochila[l]!=j){
				  if(senseMochila[i]=='L' && rmatvalMochila[l]<0) cota+=rmatvalMochila[l];
				  if(senseMochila[i]=='G' && rmatvalMochila[l]>0) cota+=rmatvalMochila[l];
				}else{
					if(rmatindMochila[l]==i) cota+=probi*rmatvalMochila[l];
					if(rmatindMochila[l]==j) cota+=probj*rmatvalMochila[l];
				}
			  }
			  if(senseMochila[k]=='L' && cota>rhsMochila[k]) serompe=true;
			  if(senseMochila[k]=='G' && cota<rhsMochila[k]) serompe=true;
			  
		  }
		  
		  if(serompe){
			  add_edge(i*2+1-probi,j*2+1-probj);
		  }
		}   
      }
  }
  cerr <<cantVar <<' '<< ejes << endl; 
  exit(-1);
  //Chequeamos todas las filas y indicamos cuales son mochila
  for(int i=0;i<cantFilasMochila;++i){
    if(senseMochila[i]=='E') continue;
    int sgn=1;
    if(senseMochila[i]=='G') sgn=-1;
    
    if(rhsMochila[i]<0)continue;
    if(fabs(rhsMochila[i]-floor(rhsMochila[i]))>10e-6)continue;
    
    bool laMeto = true;
    bool esCover= true;
    for(int j=rmatbegMochila[i];j<(i==cantFilasMochila-1?noZeroCountMochila:rmatbegMochila[i+1]);++j){
      if(rmatvalMochila[j]*sgn<0)laMeto = false;
      if(rmatvalMochila[j]!=sgn) esCover=false;
      if(fabs(rmatvalMochila[j]-floor(rmatvalMochila[j]))>10e-6)laMeto = false;
    }
    if(!esCover && laMeto && sgn==-1){
      cantG++;
      rhsMochila[i]*=sgn;
      for(int j=rmatbegMochila[i];j<(i==cantFilasMochila-1?noZeroCountMochila:rmatbegMochila[i+1]);++j) rmatvalMochila[j]*=sgn;
    }
    if(!esCover && laMeto)lasMochila.push_back(i);
    if(esCover) covers++;
  }
  
  cerr <<  endl << "Total de filas: " << cantFilasMochila << endl ;
  cerr << "Mochilas por mayor " << cantG << endl;
  cerr << "Total mochila: " << lasMochila.size() << endl;
  cerr << "Total covers: " << covers << endl;
  
  
  //~ for(int i=0;i<cantFilas;++i){
    //~ cerr << rhs[i] << endl;
  //~ }
  
  

  
  ////////////////////////////////////////////////////////////////////
  //Optimiçao na pesquisa operacional
  
  //Seteamos nuestra funcion de callback
  status = CPXsetusercutcallbackfunc(env, losCutCallbacks, NULL);
  
  
  
  status = CPXmipopt(env,lp);
  if(status){cerr << "Fallo" << endl;exit(-1);}
  
  double objval;
  status = CPXgetobjval (env, lp, &objval);
  if ( status ) {
    fprintf (stderr,"No MIP objective value available.  Exiting...\n");
    exit(-1);
  }
  
  status = CPXgetx(env, lp, xestrella, 0, cantVar-1);
  //for(int i=0;i<cantVar;++i){
    //cerr << xestrella[i] << endl;
  //}
  printf ("Valor de la solucion  = %f\n\n", objval);
  
  
  return 0;
}


///////////////////////////////////////
static int losCutCallbacks (CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p){
	int status;
	//Traemos la solucion de la relajacion
	status = CPXgetcallbacknodex(env,cbdata,wherefrom,xestrella,0,cantVar-1);
	if(status){cerr << "no pude traer x actual" << endl;exit(-1);}
	theCovernCuts(env, cbdata, wherefrom, cbhandle, useraction_p);
	//~ theGreatGomoryCuts(env, cbdata, wherefrom, cbhandle, useraction_p);
	ohMyCliqueCuts(env, cbdata, wherefrom, cbhandle, useraction_p);
	return 0;
} 

void theCovernCuts(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p){
  *useraction_p = CPX_CALLBACK_DEFAULT;
  int status;
  
  
  
  //Variable para ver que heuristica o exacto corremos
  int tengoCover = 0;
  
  #ifdef COVERGREEDY
    for(int i=0;i<(int)lasMochila.size();++i){
      tengoCover += greedyForCovern(lasMochila[i], env, cbdata, wherefrom);
    }
  #endif
  
  #ifdef COVERDP
    for(int i=0;!tengoCover&&i<(int)lasMochila.size();++i){
      tengoCover += dpForCovern(lasMochila[i], env, cbdata, wherefrom);
    }
  #endif
  
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
  
  	vector< pair<double, int> > xestrellaOrdenado;
	//Armamos subgrafo inducido por xestrella
	memset(lastSepClique, -1 , sizeof lastSepClique);
	int ejesActual = 0;
	for(int i=0;i<cantVar;++i){
		if(xestrella[i]<=0.01)continue;
		int ejeActual = last[i];
		while(ejeActual != -1){
			if(xestrella[adj[ejeActual]]>0.01){
				nextSepClique[ejesActual]=lastSepClique[i]; adjSepClique[ejes]=adj[ejeActual]; lastSepClique[i]=ejesActual++;
			}
			ejeActual = next[ejeActual];
		}
		xestrellaOrdenado.push_back(pair<double,int>(xestrella[i],i));
	}
	sort(xestrellaOrdenado.begin(),xestrellaOrdenado.end());
	
	//heuristica de busqueda clique
	
	return;
}


int dpForCovern(int fila, CPXCENVptr env, void *cbdata, int wherefrom){
  int tengoCover = 0;
  int bMochila = rhsMochila[fila] + 1;
  if(bMochila>2005){
    cerr << "Muy grande el b como para hacer dp " << bMochila << endl;
    return 0;
  }
  
  
  double dpMochila[2][bMochila+5];
  int dpPrevMochila[ULTIMO(fila)-rmatbegMochila[fila]][bMochila+5];
  memset(dpPrevMochila,-1,sizeof dpPrevMochila);
  int ant=1,act=0;
  
  dpMochila[0][0] = 0;
  dpMochila[0][1] = 0;
  
  for(int j=1;j<=bMochila;++j){
    dpMochila[ant][j] = 1e100;
    
    if(j<=rmatvalMochila[rmatbegMochila[fila]]){
      dpMochila[act][j] = (1-xestrella[rmatindMochila[rmatbegMochila[fila]]]);
      dpPrevMochila[0][j] = 0;
    }else{
      dpMochila[act][j] = 1e100;
    }
  }
  
  for(int i=rmatbegMochila[fila]+1;i<ULTIMO(fila);++i){
    swap(act,ant);
    for(int j=0;j<=bMochila;++j) dpMochila[act][j]=dpMochila[ant][j];
    
    for(int j=0;j<=bMochila;++j){  // Recorro todos los valores que tenia
      if( (1-xestrella[rmatindMochila[i]]) + dpMochila[ant][j] <= 
           dpMochila[act][min(bMochila,j+((int)rmatvalMochila[i]))])
      {
        dpMochila[act][min(bMochila,j+((int)rmatvalMochila[i]))] = (1-xestrella[rmatindMochila[i]]) + dpMochila[ant][j];
        dpPrevMochila[i-rmatbegMochila[fila]][min(bMochila,j+((int)rmatvalMochila[i]))] = j;
      }
    }
  
  }
  
  
  
  if(dpMochila[act][bMochila] < 0.99){
    cerr << "Cover violado por DP" << endl;
    //for(int i=0;i<ULTIMO(fila)-rmatbegMochila[fila];cerr<<endl,i++){
      //for(int j=0;j<=bMochila;j++) cerr << dpPrevMochila[i][j] << ' ' ;
    //}
    
    //Hay que recuperar la solucion y ahi generar el corte
    vector<int> usados;
    int columna=bMochila;
    int pfila=ULTIMO(fila)-rmatbegMochila[fila]-1;
    while(pfila>=0){
      if(dpPrevMochila[pfila][columna]==-1) pfila--;
      else{
        usados.push_back(pfila);
        columna=dpPrevMochila[pfila][columna];
        pfila--;
      }
    }
       
    
    double rhs = usados.size()-1.;
    int noZeroCount = usados.size();
    
    int *cut_matind=(int*) calloc(noZeroCount,sizeof(int));
    double *cut_matval=(double*) calloc(noZeroCount,sizeof(double));
    for(int j=0;j<(int)usados.size();++j){
      cut_matind[j] = rmatindMochila[usados[j]+rmatbegMochila[fila]];
      cut_matval[j] = 1.0;
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
  for(int i=rmatbegMochila[fila];i<ULTIMO(fila);++i){
    vec.push_back(pair<double,int>(rmatvalMochila[i]/(1-xestrella[rmatindMochila[i]]),i));
    //cerr << "Toda la info: " << "el i " << i << endl << "el valor: " << (rmatvalMochila[i]/(1-xestrella[rmatindMochila[i]])) << endl;
    //cerr << "rmatvalMochila: " << rmatvalMochila[i] << " lo otro " << 1-xestrella[rmatindMochila[i]] << endl;
  }
  sort(vec.begin(),vec.end(), greater<pair<double,int> >());
  
  double acumA=0.0, acumPesos=0.0;
  
  //Ahora sumamos las a_j hasta pasarnos del bMochila+1
  int i;
  
  //Para extender la cover!
  int maximoValor=0; 
  
  for(i=0;i<(int)vec.size()&&acumA<bMochila+0.5;++i){
    maximoValor=max(maximoValor,(int)rmatvalMochila[vec[i].second]);
    acumA += rmatvalMochila[vec[i].second];
    acumPesos += (1.-xestrella[rmatindMochila[vec[i].second]]);
  }
  
  //Preguntamos si los pesos suman menos que 1 (si es asi, hay un corte violado)
  if(acumPesos < .98 && acumA>bMochila+0.5){
    double rhs = i-1.;
    int noZeroCount = i;
    
    
    //Vamos a extender la mochila con los que tienen A mayor o igual al maximo
    int extendemos=0;
    for(int j=rmatbegMochila[fila];j<ULTIMO(fila);j++){
      if(rmatvalMochila[j]>=maximoValor) extendemos++;
    }
    
    
    int *cut_matind=(int*) calloc(noZeroCount+extendemos,sizeof(int));
    double *cut_matval=(double*) calloc(noZeroCount+extendemos,sizeof(double));
    
    int ptr=0;
    for(int j=0;j<i;++j){
      cut_matind[ptr] = rmatindMochila[vec[j].second];
      cut_matval[ptr++] = 1.0;
    }
    
    // Las extendidas
    for(int j=rmatbegMochila[fila];j<ULTIMO(fila);j++){
      if(rmatvalMochila[j]>=maximoValor){
        cut_matind[ptr]  = rmatindMochila[j];
        cut_matval[ptr++]= 1.0; 
      }
    }
    sort(cut_matind,cut_matind+ptr);
    int distintas=unique(cut_matind,cut_matind+ptr)-cut_matind;
    
    cerr << "Cover violado. Lo extendimos por " << distintas-noZeroCount<<".Van "<<(++totCortesGreedy)<<" !"  << endl;
    
    CPXcutcallbackadd(env, cbdata, wherefrom, distintas, rhs, 'L', cut_matind, cut_matval, 0);
    tengoCover++;
  }
  return tengoCover;
}

/*

  28335 18180  -116599.4537   486  -103157.5365  -119789.1446  1001733   16.12%
Elapsed real time = 5025.96 sec. (tree size = 886.63 MB, solutions = 7)
Nodefile size = 759.10 MB (435.38 MB after compression)
  28541 18245  -117307.9430   541  -103157.5365  -119763.6160  1011636   16.10%







    
Elapsed real time = 945.76 sec. (tree size = 177.72 MB, solutions = 0)
Nodefile size = 50.65 MB (30.58 MB after compression)
   4343  3405   -98593.7680    86                -122137.0054   185937         
   4572  3544  -116444.7307   474                -122065.4128   194610         
Cover violado
   4797  3678  -121166.6174   448                -121992.2699   204321         
   4996  3814  -111350.2882   419                -121925.6843   213555         
   5305  4117  -104746.6039    96                -121925.6843   220075         
*  5477  4255      integral     0   -95024.2072  -121925.6843   222425   28.31%
   5674  4360  -121856.2725   597   -95024.2072  -121847.2772   231677   28.23%
   5872  4533  -107251.7208   286   -95024.2072  -121837.6007   241987   28.22%
Cover violado









*/
