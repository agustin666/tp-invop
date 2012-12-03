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
 
//#define COVERGREEDY
//#define COVERDP

#define INSTANCIA "miplib3/egout.mps"
#define JOYA if(status)errorHandler(status,env)
#define ULTIMO(fila) (fila==cantFilasMochila-1?noZeroCountMochila:rmatbegOriginal[fila+1])

static int losCutCallbacks(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p);
void theCovernCuts(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p);
void theGreatGomoryCuts(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p);
void ohMyCliqueCuts(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p);
int greedyForCovern(int fila,CPXCENVptr env,void *cbdata,int wherefrom);
int dpForCovern(int fila, CPXCENVptr env, void *cbdate, int wherefrom);



double* xestrella;
double* aOriginal;
double* rhsOriginal;
int* rmatbegOriginal;
int cantFilasMochila;
double* rmatvalOriginal;
char* senseOriginal;
int* rmatindOriginal;
vector<int> lasMochila;
double bMochila;
int cantVar;
int noZeroCountMochila;


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
  aOriginal = (double*)malloc(sizeof(double)*cantVar);
  cantFilasMochila = CPXgetnumrows(env,lp);
  int surplus;
  int rmatspaceOriginal = min( cantVar*cantFilasMochila , 1000000);
  rhsOriginal = (double*)malloc(sizeof(double)*cantFilasMochila);
  rmatbegOriginal = (int*)malloc(sizeof(int)*cantFilasMochila);
  rmatindOriginal = (int*)malloc(sizeof(int)*rmatspaceOriginal);
  rmatvalOriginal = (double*)malloc(sizeof(double)*rmatspaceOriginal);
  senseOriginal = (char*)  malloc(sizeof(char)*cantFilasMochila);
  status = CPXgetsense(env, lp, senseOriginal, 0, cantFilasMochila-1);
  JOYA;
  status = CPXgetrhs(env, lp, rhsOriginal, 0, cantFilasMochila-1);
  JOYA;
  status = CPXgetrows(env, lp, &noZeroCountMochila, rmatbegOriginal, rmatindOriginal, rmatvalOriginal, rmatspaceOriginal, &surplus, 0, cantFilasMochila-1);
  if(surplus<0){cerr << "No teniamos espacio para traer las filas" << endl; errorHandler(status,env);}
  JOYA;
  
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
  memset(last,-1,2*cantVar*sizeof(int));
  for(int i=0;i<cantVar;++i){
      for(int j=i+1;j<cantVar;++j){
		bool serompe;
		for(int a=0;a<4;++a){
		  int probi=a&1;
		  int probj=(a&2)?1:0;
		  serompe=false;
		  
		  
		  //cantFilasMochila
		  for(int k=0;!serompe && k<cantFilasMochila;++k){
			  if(senseOriginal[i]=='E') continue;
			  int cota=0;
			  
			  for(int l=rmatbegOriginal[k];l<(k==cantFilasMochila-1?noZeroCountMochila:rmatbegOriginal[k+1]);++l){
				if(rmatindOriginal[l]!=i && rmatindOriginal[l]!=j){
				  if(senseOriginal[i]=='L' && rmatvalOriginal[l]<0) cota+=rmatvalOriginal[l];
				  if(senseOriginal[i]=='G' && rmatvalOriginal[l]>0) cota+=rmatvalOriginal[l];
				}else{
					if(rmatindOriginal[l]==i) cota+=probi*rmatvalOriginal[l];
					if(rmatindOriginal[l]==j) cota+=probj*rmatvalOriginal[l];
				}
			  }
			  if(senseOriginal[k]=='L' && cota>rhsOriginal[k]) serompe=true;
			  if(senseOriginal[k]=='G' && cota<rhsOriginal[k]) serompe=true;
			  
		  }
		  
		  if(serompe){
			  add_edge(i*2+1-probi,j*2+1-probj);
		  }
		}   
      }
  }
  
  cerr << "EL NEXT DE CERO ES " << next[0] << endl;
  cerr <<"Ejes en el grafo de conflicto "<< ejes << endl; 
  //~ exit(-1);
  
  //Chequeamos todas las filas y indicamos cuales son Original
  for(int i=0;i<cantFilasMochila;++i){
    if(senseOriginal[i]!='L') continue;
    if(rhsOriginal[i]<0)continue;
    if(fabs(rhsOriginal[i]-floor(rhsOriginal[i]))>10e-6)continue;
    bool laMeto = true;
    for(int j=rmatbegOriginal[i];j<(i==cantFilasMochila-1?noZeroCountMochila:rmatbegOriginal[i+1]);++j){
      if(rmatvalOriginal[j]<0)laMeto = false;
      if(fabs(rmatvalOriginal[j]-floor(rmatvalOriginal[j]))>10e-6)laMeto = false;
    }
    if(laMeto)lasMochila.push_back(i);
  }
  
  cerr <<  endl << "Total de filas: " << cantFilasMochila << endl << "Total Original: " << lasMochila.size() << endl << endl;
  
  
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

#define COTACLIQUE 20
int szclique;
int miclique[1010];
int cantcliques;

#define DBG

void genclique(int nodo, int pos,vector< pair<double, int> > &xestrella,double acum,void *cbdata,int wherefrom,CPXCENVptr env){
	int entro=false; szclique++;
	
	#ifdef DBG
		cerr << nodo << ' ' << acum << endl;
	#endif 
	
	miclique[szclique-1]=nodo;
	
	for(int i=pos+1;i<xestrella.size() && cantcliques<COTACLIQUE;++i){
		bool puedo=true;
		
		
		for(int v=xestrella[i].second*2;v<=xestrella[i].second*2+1;v++){

			
			for(int j=0;j<szclique;++j){
				bool hayeje = false;
				for(int k=last[v];~k && !hayeje;k=next[k]){
					if(adj[k]==miclique[j]) hayeje=true;
				}
				if(!hayeje) puedo=false;
			}
			if(puedo){
				entro=true;
				genclique(v,i,xestrella, acum+(v&1?0:xestrella[i].first),cbdata,wherefrom,env);
			}
		}
		
	}
	
	
	
	
	if(!entro && cantcliques!=COTACLIQUE){ // es maximal
		cantcliques++;
		
		if(acum>1.05){
			cerr<<"Metimos clique cut "<< acum << endl;
			
			double rhs = 1.;
			
			int noZeroCount = szclique;
			
			int *cut_matind=(int*) calloc(noZeroCount,sizeof(int));
			double *cut_matval=(double*) calloc(noZeroCount,sizeof(double));
			
			for(int j=0;j<szclique;++j){
			  cut_matind[j] = miclique[j]/2;
			  cut_matval[j] = 1.0;
			  if(miclique[j]&1) {  // Es el complemento
				  rhs-=1.;
				  cut_matval[j] = -1.0;
			  }
			}
			
			CPXcutcallbackadd(env, cbdata, wherefrom, noZeroCount, rhs, 'L', cut_matind, cut_matval, 0);
			
		}
		if(cantcliques==COTACLIQUE) return;
	}
	
	szclique--;
}

void ohMyCliqueCuts(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p){
	*useraction_p = CPX_CALLBACK_DEFAULT;
	//~ cerr << "oh baby refrain, cuz you are breaking my clique" << endl;
  
  	vector< pair<double, int> > xestrellaOrdenado;
	//Armamos subgrafo inducido por xestrella
	//~ memset(lastSepClique, -1 , sizeof lastSepClique);
	//~ int ejesActual = 0;
	for(int i=0;i<cantVar;++i){
		if(xestrella[i]<=0.01)continue;
		//~ int ejeActual = last[i];
		//~ while(ejeActual != -1){
			//~ if(xestrella[adj[ejeActual]]>0.01){
				//~ nextSepClique[ejesActual]=lastSepClique[i]; adjSepClique[ejesActual]=adj[ejeActual]; lastSepClique[i]=ejesActual++;
			//~ }
			//~ ejeActual = next[ejeActual];
		//~ }
		xestrellaOrdenado.push_back(pair<double,int>(xestrella[i],i));
	}
	
	
	
	sort(xestrellaOrdenado.begin(),xestrellaOrdenado.end());
	cantcliques=0;
	szclique=0;
	cerr << "Armo subgrafo " << endl;
	
	for(int i=0;i<cantVar;++i){
		genclique(xestrellaOrdenado[i].second*2,i,xestrellaOrdenado,xestrellaOrdenado[i].first,cbdata,wherefrom,env);
	}
	
	return;
}


int dpForCovern(int fila, CPXCENVptr env, void *cbdata, int wherefrom){
  int tengoCover = 0;
  int bMochila = rhsOriginal[fila] + 1;
  if(bMochila>2005){
    cerr << "Muy grande el b como para hacer dp " << bMochila << endl;
    return 0;
  }
  
  
  double dpOriginal[2][bMochila+5];
  int dpPrevOriginal[ULTIMO(fila)-rmatbegOriginal[fila]][bMochila+5];
  memset(dpPrevOriginal,-1,sizeof dpPrevOriginal);
  int ant=1,act=0;
  
  dpOriginal[0][0] = 0;
  dpOriginal[0][1] = 0;
  
  for(int j=1;j<=bMochila;++j){
    dpOriginal[ant][j] = 1e100;
    
    if(j<=rmatvalOriginal[rmatbegOriginal[fila]]){
      dpOriginal[act][j] = (1-xestrella[rmatindOriginal[rmatbegOriginal[fila]]]);
      dpPrevOriginal[0][j] = 0;
    }else{
      dpOriginal[act][j] = 1e100;
    }
  }
  
  for(int i=rmatbegOriginal[fila]+1;i<ULTIMO(fila);++i){
    swap(act,ant);
    for(int j=0;j<=bMochila;++j) dpOriginal[act][j]=dpOriginal[ant][j];
    
    for(int j=0;j<=bMochila;++j){  // Recorro todos los valores que tenia
      if( (1-xestrella[rmatindOriginal[i]]) + dpOriginal[ant][j] <= 
           dpOriginal[act][min(bMochila,j+((int)rmatvalOriginal[i]))])
      {
        dpOriginal[act][min(bMochila,j+((int)rmatvalOriginal[i]))] = (1-xestrella[rmatindOriginal[i]]) + dpOriginal[ant][j];
        dpPrevOriginal[i-rmatbegOriginal[fila]][min(bMochila,j+((int)rmatvalOriginal[i]))] = j;
      }
    }
  
  }
  
  
  
  if(dpOriginal[act][bMochila] < 0.99){
    cerr << "Cover violado por DP" << endl;
    //for(int i=0;i<ULTIMO(fila)-rmatbegOriginal[fila];cerr<<endl,i++){
      //for(int j=0;j<=bMochila;j++) cerr << dpPrevOriginal[i][j] << ' ' ;
    //}
    
    //Hay que recuperar la solucion y ahi generar el corte
    vector<int> usados;
    int columna=bMochila;
    int pfila=ULTIMO(fila)-rmatbegOriginal[fila]-1;
    while(pfila>=0){
      if(dpPrevOriginal[pfila][columna]==-1) pfila--;
      else{
        usados.push_back(pfila);
        columna=dpPrevOriginal[pfila][columna];
        pfila--;
      }
    }
       
    
    double rhs = usados.size()-1.;
    int noZeroCount = usados.size();
    
    int *cut_matind=(int*) calloc(noZeroCount,sizeof(int));
    double *cut_matval=(double*) calloc(noZeroCount,sizeof(double));
    for(int j=0;j<(int)usados.size();++j){
      cut_matind[j] = rmatindOriginal[usados[j]+rmatbegOriginal[fila]];
      cut_matval[j] = 1.0;
    }
    CPXcutcallbackadd(env, cbdata, wherefrom, noZeroCount, rhs, 'L', cut_matind, cut_matval, 0);
    tengoCover++;
  }
  
  return tengoCover;
}

int greedyForCovern(int fila,CPXCENVptr env,void *cbdata,int wherefrom ){
  int tengoCover = 0;  
  bMochila = rhsOriginal[fila];
  
  
  //Greedy 1 (A_j/(1-x_j)
  vector<pair<double, int> > vec;
  for(int i=rmatbegOriginal[fila];i<ULTIMO(fila);++i){
    vec.push_back(pair<double,int>(rmatvalOriginal[i]/(1-xestrella[rmatindOriginal[i]]),i));
    //cerr << "Toda la info: " << "el i " << i << endl << "el valor: " << (rmatvalOriginal[i]/(1-xestrella[rmatindOriginal[i]])) << endl;
    //cerr << "rmatvalOriginal: " << rmatvalOriginal[i] << " lo otro " << 1-xestrella[rmatindOriginal[i]] << endl;
  }
  sort(vec.begin(),vec.end(), greater<pair<double,int> >());
  
  double acumA=0.0, acumPesos=0.0;
  
  //Ahora sumamos las a_j hasta pasarnos del bMochila+1
  int i;
  for(i=0;i<(int)vec.size()&&acumA<bMochila+0.5;++i){
    acumA += rmatvalOriginal[vec[i].second];
    acumPesos += (1.-xestrella[rmatindOriginal[vec[i].second]]);
  }
  
  //Preguntamos si los pesos suman menos que 1 (si es asi, hay un corte violado)
  if(acumPesos < .99 && acumA>bMochila+0.5){
    cerr << "Cover violado" << endl;
    double rhs = i-1.;
    //cerr << "Dios mio alto corte: " << endl << endl;
    //cerr << rhs << endl;
    int noZeroCount = i;
    
    int *cut_matind=(int*) calloc(noZeroCount,sizeof(int));
    double *cut_matval=(double*) calloc(noZeroCount,sizeof(double));
    for(int j=0;j<i;++j){
      cut_matind[j] = rmatindOriginal[vec[j].second];
      cut_matval[j] = 1.0;
      //cerr << rmatindOriginal[vec[j].second] << endl;
    }
    CPXcutcallbackadd(env, cbdata, wherefrom, noZeroCount, rhs, 'L', cut_matind, cut_matval, 0);
    tengoCover++;
  }
  return tengoCover;
}
