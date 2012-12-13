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

// Opciones para correr el programa
int XTIEMPO=0;
int VERBOSE=0;
int COVERGREEDY=0;
int COVERDP=0;
int CORTECLIQUE=0;
int CUTANDBRANCH=0;


//Instancia Default si no se indica el archivo
#define INSTANCIA "miplib3/harp2.mps"

//Definimos unas macros para mejorar un poco la legibilidad
#define ULTIMO(fila) (fila==cantFilas-1?noZeroCountMochila:rmatbegOriginal[fila+1])
#define JOYA if(status){CPXgeterrorstring(env, status, errmsg);fprintf(stderr, "%s", errmsg);exit(-1);}
char errmsg[CPXMESSAGEBUFSIZE];


static int losCutCallbacks(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p);
void theCovernCuts(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p);
void ohMyCliqueCuts(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p);
int greedyForCovern(int fila,CPXCENVptr env,void *cbdata,int wherefrom);
int dpForCovern(int fila, CPXCENVptr env, void *cbdate, int wherefrom);


//Vamos a guardarnos todo el LP original en memoria
int cantVar;
int cantFilas;
double* xestrella;
double* aOriginal;
double* rhsOriginal;
int* rmatbegOriginal;
double* rmatvalOriginal;
char* senseOriginal;
int* rmatindOriginal;
int noZeroCountMochila;

//En un vector nos guardamos el indice de las restricciones que son mochila, y no son cover
vector<int> lasMochila;
double bMochila;
int totCortesGreedy=0;
int totCortesClique=0;

//Estructura de grafo
int* last,*next,*adj;
int ejes; 

//Arreglos para meter cortes
int *cut_matind;
double *cut_matval;
int *ord;

void add_edge(int u,int v){
    next[ejes]=last[u]; adj[ejes]=v; last[u]=ejes++;
    swap(u,v);
    next[ejes]=last[u]; adj[ejes]=v; last[u]=ejes++;
}

int main(int argc, char *argv[]){
  string archivo;
  if(argc==1 || isParam("-h",argv,argc)){
    puts("usage: ./Tp [opciones]");
    puts("");
    puts("-f [archivo] usamos al archivo de entrada");
    puts("-c  Usa todos los cortes de CPLEX");
    puts("-G  Usa cortes cover greedy");
    puts("-D  Usa cortes cover con DP");
    puts("-K  Usa cortes clique");
    puts("-X  Corta en el primer nodo (Cut & Branch)");
    puts("-v  Indica cuando se realizan los cortes");
    puts("-t  Limite de tiempo 5 minutos");
    exit(0);
  }
  if(isParam("-G",argv,argc)) COVERGREEDY=1;
  if(isParam("-D",argv,argc)) COVERDP=1;
  if(isParam("-X",argv,argc)) CUTANDBRANCH=1;
  if(isParam("-K",argv,argc)) CORTECLIQUE=1;
  if(isParam("-v",argv,argc)) VERBOSE=1;
  if(isParam("-t",argv,argc)) XTIEMPO=1;
  
  if(isParam("-f", argv, argc)){
    archivo = getParam("-f", argv, argc);
  }else{
    archivo = INSTANCIA;
  }
  
  
  ////////////////////////////////////////////////////////////////////
  //Arrancamos CPLEX
  int status;
  CPXENVptr env = NULL;
  CPXLPptr lp = NULL;
  env = CPXopenCPLEX(&status);
  
  if(env == NULL){
    fprintf(stderr, "No se pudo inicializar el ambiente.\n");JOYA;
  }
  
  lp = CPXcreateprob(env, &status, "Test");
  
  if(lp==NULL){
    fprintf(stderr, "No se pudo crear el problema\n", status);JOYA;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////
  //Parametros de cplex
  //
  status = CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
  JOYA;
  
  status = CPXsetintparam(env, CPX_PARAM_DATACHECK, CPX_ON);
  JOYA;
  
  //Con este parametro en off, siempre no estamos refiriendo al problema original
  status = CPXsetintparam (env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
  JOYA;
  
  //Hay una opcion para correr cplex con todo  
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
  }
  
  //Lectura de instancia
  CPXreadcopyprob(env,lp,archivo.c_str(), NULL);
  
  /* Con esto nos quedamos con la informacion del nodo raiz
  status = CPXsetintparam(env, CPX_PARAM_NODELIM, 0);
	if(status)exit(-1);
  */

  // Copiamos el LP a memoria global
  cantVar = CPXgetnumcols(env,lp);
  xestrella = (double*)malloc(sizeof(double)*cantVar);
  aOriginal = (double*)malloc(sizeof(double)*cantVar);

  cantFilas = CPXgetnumrows(env,lp);
  int surplus;
  int rmatspaceOriginal = min( cantVar*cantFilas , 1000000);
  rhsOriginal = (double*)malloc(sizeof(double)*cantFilas);
  rmatbegOriginal = (int*)malloc(sizeof(int)*cantFilas);
  rmatindOriginal = (int*)malloc(sizeof(int)*rmatspaceOriginal);
  rmatvalOriginal = (double*)malloc(sizeof(double)*rmatspaceOriginal);
  senseOriginal = (char*)  malloc(sizeof(char)*cantFilas);
  status = CPXgetsense(env, lp, senseOriginal, 0, cantFilas-1);
  JOYA;
  status = CPXgetrhs(env, lp, rhsOriginal, 0, cantFilas-1);
  JOYA;
  status = CPXgetrows(env, lp, &noZeroCountMochila, rmatbegOriginal, rmatindOriginal, rmatvalOriginal, rmatspaceOriginal, &surplus, 0, cantFilas-1);
  if(surplus<0){cerr << "No teniamos espacio para traer las filas" << endl;}
  JOYA;

  
  // Inicializamos el grafo, usamos la misma estructura que utiliza el paper de listas de ejes
  last=(int*)malloc(2*cantVar*sizeof(int));
  adj=(int*)malloc(min(1000000,4*cantVar*cantVar)*sizeof(int));
  next=(int*)malloc(min(1000000,4*cantVar*cantVar)*sizeof(int));
  int cantFilasClique=cantFilas;
  ejes=0;
  memset(last,-1,2*cantVar*sizeof(int));
  
  if(CORTECLIQUE){
    for(int i=0;i<cantVar;++i){
      add_edge(i*2,i*2+1);          // i y ~i son conflicto
      for(int j=i+1;j<cantVar;++j){ 
        bool serompe;
        for(int a=0;a<4;++a){ // Probamos las 4 posiblidades 00,01,10,11
          int probi=a&1;      
          int probj=(a&2)?1:0;
          serompe=false;
          
          for(int k=0;!serompe && k<cantFilas;++k){  // Tenemos que buscar alguna restriccion infactible
            if(senseOriginal[i]=='E') continue;
            double cota=0;
            
            for(int l=rmatbegOriginal[k];l<(k==cantFilas-1?noZeroCountMochila:rmatbegOriginal[k+1]);++l){
              if(rmatindOriginal[l]!=i && rmatindOriginal[l]!=j){
                if(senseOriginal[k]=='L' && rmatvalOriginal[l]<0) cota+=rmatvalOriginal[l];
                if(senseOriginal[k]=='G' && rmatvalOriginal[l]>0) cota+=rmatvalOriginal[l];
              }else{
                if(rmatindOriginal[l]==i) cota+=probi*rmatvalOriginal[l]; // Si probi es 1, uso la variable, sino, no la uso
                if(rmatindOriginal[l]==j) cota+=probj*rmatvalOriginal[l];
              }
            }
            if(senseOriginal[k]=='L' && cota>rhsOriginal[k]+0.001) serompe=true;
            if(senseOriginal[k]=='G' && cota<rhsOriginal[k]-0.001) serompe=true;
          }
          if(serompe){          
            add_edge(i*2+1-probi,j*2+1-probj);
          }
        }   
      }
    }
    cerr <<"Ejes en el grafo de conflicto "<< ejes << endl; 
  }
  
  //Chequeamos todas las filas y indicamos cuales son mochila. Excluimos cover y cambiamos signo
  int covers=0;
  int cantG=0;
  for(int i=0;i<cantFilas;++i){
    if(senseOriginal[i]=='E') continue;
    int sgn=1;
    if(senseOriginal[i]=='G') sgn=-1;
    
    if(rhsOriginal[i]<0)continue;
    if(fabs(rhsOriginal[i]-floor(rhsOriginal[i]))>10e-6)continue;
    
    bool laMeto = true;
    bool esCover= true;
    for(int j=rmatbegOriginal[i];j<(i==cantFilas-1?noZeroCountMochila:rmatbegOriginal[i+1]);++j){
      if(rmatvalOriginal[j]*sgn<0)laMeto = false;
      if(rmatvalOriginal[j]!=sgn) esCover=false;
      if(fabs(rmatvalOriginal[j]-floor(rmatvalOriginal[j]))>10e-6)laMeto = false;
    }
    if(!esCover && laMeto && sgn==-1){
      cantG++;
      rhsOriginal[i]*=sgn;
      for(int j=rmatbegOriginal[i];j<(i==cantFilas-1?noZeroCountMochila:rmatbegOriginal[i+1]);++j) rmatvalOriginal[j]*=sgn;
    }
    if(!esCover && laMeto)lasMochila.push_back(i);
    if(esCover) covers++;
   }
  
  if(VERBOSE){
    cerr <<  endl << "Total de variables: " << cantVar << endl ;
    cerr << "Total de filas: " << cantFilas << endl ;
    cerr << "Mochilas por mayor " << cantG << endl;
    cerr << "Total mochila: " << lasMochila.size() << endl;
    cerr << "Total covers: " << covers << endl;
  }
  //Para resultados del informe usamos limite de tiempo
  if(XTIEMPO){
    status = CPXsetdblparam(env, CPX_PARAM_TILIM, 6400);
    JOYA;
  }
  
  //Memoria para los cuts
  cut_matind=(int*) calloc(cantVar,sizeof(int));
  cut_matval=(double*) calloc(cantVar,sizeof(double));
  
  //En los cortes cliques voy a ordenar las variables
  ord=(int*) calloc(cantVar,sizeof(int));
  for(int i=0;i<cantVar;i++) ord[i]=i;
  
  ////////////////////////////////////////////////////////////////////
  //Optimiçao na pesquisa operacional
  //Seteamos nuestra funcion de callback
  status = CPXsetusercutcallbackfunc(env, losCutCallbacks, NULL);
  
  double start,end;
  CPXgettime(env,&start);
  status = CPXmipopt(env,lp);
  CPXgettime(env,&end);
  double tardo=end-start;
  
  double objval,gap;
  int nodecount=CPXgetnodecnt(env,lp);
  status = CPXgetobjval (env, lp, &objval);
  status = CPXgetmiprelgap (env, lp, &gap);
  cerr << "Total covers encontradas " << totCortesGreedy << endl;
  cerr << "Total clique encontradas " << totCortesClique << endl;
  fprintf(stderr,"Tiempo  = %lf\n", tardo);
  fprintf(stderr,"#nodos  = %d\n", nodecount);
  fprintf(stderr,"Gap  = %.2lf%%\n", gap*100);
  fprintf(stderr,"Valor de la solucion  = %lf\n", objval);
  
  return 0;
}


///////////////////////////////////////
static int losCutCallbacks (CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p){
  int status;
  
  //Cut-and-Branch
  if(CUTANDBRANCH){
    CPXINT result;
    status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_DEPTH, &result);
    if(result){return 0;}
  }
  ///////////////////
  
  //Traemos la solucion de la relajacion
  status = CPXgetcallbacknodex(env,cbdata,wherefrom,xestrella,0,cantVar-1);
  JOYA;
  
  if(COVERDP || COVERGREEDY) theCovernCuts(env, cbdata, wherefrom, cbhandle, useraction_p);
  if(CORTECLIQUE) ohMyCliqueCuts(env, cbdata, wherefrom, cbhandle, useraction_p);
  
  return 0;
} 

void theCovernCuts(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p){
  *useraction_p = CPX_CALLBACK_DEFAULT;
  int status;
  
  //Variable para ver que heuristica o exacto corremos
  int tengoCover = 0;
  
  if(COVERGREEDY){
    for(int i=0;i<(int)lasMochila.size()&&tengoCover==0;i++){
      tengoCover += greedyForCovern(lasMochila[i], env, cbdata, wherefrom);
    }
  }
  
  if(COVERDP){
    for(int i=0;tengoCover&&i<(int)lasMochila.size();++i){
      tengoCover += dpForCovern(lasMochila[i], env, cbdata, wherefrom);
    }
  }
  
  if(tengoCover) *useraction_p = CPX_CALLBACK_SET;
  return;
}


// Maxima cantidad de cliques que vamos a explorar
#define COTACLIQUE 10
int miclique[1000];
int szclique;
int cantcliques;

void genclique(int nodo, int pos,double acum,void *cbdata,int wherefrom,CPXCENVptr env){
  int entro=false; 
  miclique[szclique]=nodo;
  szclique++;
  
  for(int i=pos+1;i<cantVar && cantcliques<COTACLIQUE;++i){
    for(int v=ord[i]*2;v<=ord[i]*2+1  && cantcliques<COTACLIQUE;v++){
      int j;
      
      for(j=0;j<szclique;++j){
        bool hayeje = false;
        for(int k=last[v];~k && !hayeje;k=next[k]){
          if(adj[k]==miclique[j]) hayeje=true;
        }
        if(!hayeje) break;
      }
      if(j==szclique){    // hay eje con todas, puedo extender la clique
        entro=true;
        genclique(v,i,acum+(v&1?0:xestrella[ord[i]]),cbdata,wherefrom,env);
      }
    }
  }
  
  //Si no pudo extender clique era maximal
  if(!entro && cantcliques!=COTACLIQUE){ 
    cantcliques++;
    if(acum>=1.001){
      if(VERBOSE){
        cerr<<"Metimos clique cut de tamanyo "<< szclique<< endl;
      }
      totCortesClique++;
      
      double rhs = 1.;
      int noZeroCount = szclique;
      
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

// Funcion de comparacion de variables por su valor en la relajacion
bool byXestrella(int a,int b){
  return xestrella[a]>xestrella[b];
}

void ohMyCliqueCuts(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p){
  *useraction_p = CPX_CALLBACK_DEFAULT;
  
  // Vamos a tratar de armar cliques priorizando los ejes 'mas presentes'
  sort(ord,ord+cantFilas,byXestrella);
  
  cantcliques=0;
  szclique=0;
  for(int i=0;i<cantVar && cantcliques<COTACLIQUE;++i){
    genclique(ord[i]*2,i,xestrella[ord[i]],cbdata,wherefrom,env);
  }
  
  return;
}


// Como inicialmente el dp no parecia muy prometedor, ni tan confiable, porque el algo pseudo polinomial puede requerir mucha memoria
// quedo medio incompleto, por ejemplo, asignar arreglos dinamicos fuera de la funcion y extender los covers
int dpForCovern(int fila, CPXCENVptr env, void *cbdata, int wherefrom){
  int tengoCover = 0;
  int bMochila = rhsOriginal[fila] + 1;
  if(bMochila>2005){
    return 0;
  }
  
  double dpMochila[2][bMochila+5];
  int dpPrevMochila[ULTIMO(fila)-rmatbegOriginal[fila]][bMochila+5];
  memset(dpPrevMochila,-1,sizeof dpPrevMochila);
  int ant=1,act=0;
  
  dpMochila[0][0] = 0;
  dpMochila[0][1] = 0;
  
  for(int j=1;j<=bMochila;++j){
    dpMochila[ant][j] = 1e100;
    
    if(j<=rmatvalOriginal[rmatbegOriginal[fila]]){
      dpMochila[act][j] = (1-xestrella[rmatindOriginal[rmatbegOriginal[fila]]]);
      dpPrevMochila[0][j] = 0;
    }else{
      dpMochila[act][j] = 1e100;
    }
  }
  
  for(int i=rmatbegOriginal[fila]+1;i<ULTIMO(fila);++i){
    swap(act,ant);
    for(int j=0;j<=bMochila;++j) dpMochila[act][j]=dpMochila[ant][j];
    
    for(int j=0;j<=bMochila;++j){  // Recorro todos los valores que tenia
      if( (1-xestrella[rmatindOriginal[i]]) + dpMochila[ant][j] <= 
           dpMochila[act][min(bMochila,j+((int)rmatvalOriginal[i]))])
      {
        dpMochila[act][min(bMochila,j+((int)rmatvalOriginal[i]))] = (1-xestrella[rmatindOriginal[i]]) + dpMochila[ant][j];
        dpPrevMochila[i-rmatbegOriginal[fila]][min(bMochila,j+((int)rmatvalOriginal[i]))] = j;
      }
    }
  }
  
  if(dpMochila[act][bMochila] < 0.99){
    if(VERBOSE){
      cerr << "Cover violado por DP" << endl;
    }
    
    //Hay que recuperar la solucion y ahi generar el corte
    vector<int> usados;
    int columna=bMochila;
    int pfila=ULTIMO(fila)-rmatbegOriginal[fila]-1;
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
    for(int j=0;j<(int)usados.size();++j){
      cut_matind[j] = rmatindOriginal[usados[j]+rmatbegOriginal[fila]];
      cut_matval[j] = 1.0;
    }
    CPXcutcallbackadd(env, cbdata, wherefrom, noZeroCount, rhs, 'L', cut_matind, cut_matval, 0);
    tengoCover++;
  }
  
  return tengoCover;
}


// Maxima cantidad de variables que vamos a manejar
#define COTATAM 15000
pair<double,int> vec[COTATAM];
int greedyForCovern(int fila,CPXCENVptr env,void *cbdata,int wherefrom ){
  int tengoCover = 0;  
  bMochila = rhsOriginal[fila];
  
  //Greedy 1 (A_j/(1-x_j)
  //vector<pair<double, int> > vec;
  int len=0;
  
  for(int i=rmatbegOriginal[fila];i<ULTIMO(fila);++i){
    vec[len++]=pair<double,int>(rmatvalOriginal[i]/(1-xestrella[rmatindOriginal[i]]),i);
  }
  sort(vec,vec+len, greater<pair<double,int> >());
  double acumA=0.0, acumPesos=0.0;
  
  //Ahora sumamos las a_j hasta pasarnos del bMochila+1
  int i;
  
  //Para extender la cover!
  int maximoValor=0; 
  
  for(i=0;i<len&&acumA<bMochila+0.5;++i){
    maximoValor=max(maximoValor,(int)rmatvalOriginal[vec[i].second]);
    acumA += rmatvalOriginal[vec[i].second];
    acumPesos += (1.-xestrella[rmatindOriginal[vec[i].second]]);
  }
  
  //Preguntamos si los pesos suman menos que 1 (si es asi, hay un corte violado)
  if(acumPesos < .98 && acumA>bMochila+0.5){
    totCortesGreedy++;
    double rhs = i-1.;
    int noZeroCount = i;
    
    //Vamos a extender la mochila con los que tienen A mayor o igual al maximo
    int extendemos=0;
    for(int j=rmatbegOriginal[fila];j<ULTIMO(fila);j++){
      if(rmatvalOriginal[j]>=maximoValor) extendemos++;
    }    
    
    int ptr=0;
    for(int j=0;j<i;++j){
      cut_matind[ptr] = rmatindOriginal[vec[j].second];
      cut_matval[ptr++] = 1.0;
    }
    
    // Las extendidas
    for(int j=rmatbegOriginal[fila];j<ULTIMO(fila);j++){
      if(rmatvalOriginal[j]>=maximoValor){
        cut_matind[ptr]  = rmatindOriginal[j];
        cut_matval[ptr++]= 1.0; 
      }
    }
    sort(cut_matind,cut_matind+ptr);
    int distintas=unique(cut_matind,cut_matind+ptr)-cut_matind;
    
    CPXcutcallbackadd(env, cbdata, wherefrom, distintas, rhs, 'L', cut_matind, cut_matval, 0);
    return 1;
  }
  return tengoCover;
}

