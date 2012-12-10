#include <iostream>
#include <cstdio>
#include <vector>
#include <cstring>
#include <algorithm>
using namespace std;

int N,R;
int A[20][20];
int rhs[20];

int generarclique(){
  N=random()%6+4;
  
  R=random()%N+2;
  
  memset(A,0,sizeof A);
  
  int i,j,k;
  
  int mx=-1;
  //Todas los nodos estan en alguna constint
  for(i=0;i<N;i++){
    int id=random()%R;
    A[id][i]=random()%50+1;
    mx=max(mx,A[id][i]);      // No quiero que ninguna variable me vuelva infeasible si 
  }                           // esta prendida entonces fuerzo el rhs a ser mayor al maximo
  
  
  // Gernero RHS random. Deberia poder meter cualquier eje al grafo de conflicto
  for(i=0;i<R;i++) rhs[i]=random()%1000+mx+2;
  
  //cerr << "llega " << endl;
  // Trato de forzar el random, para cada par de nodos
  // Alguna restriccion impide que esten ambas
  for(i=1;i<N;i++){
    for(j=0;j<i;j++){
      int mask1=0;  // mask1 me indica las restricciones que tienen a i, y mask2 a j
      int mask2=0;
      for(k=0;k<R;k++){
        if(A[k][i]>0) mask1|=(1<<k);
        if(A[k][j]>0) mask2|=(1<<k);
      }
      
      if((mask1^mask2)==0){ // En todas las resticciones estan las dos, o no esta ninguna
        // No se que hacer, lo que haya en alguna restriccion lo incrementamos para que se rompa
        // Voy subiendo de a poco hasta que se rompa
        int P=random()%R;
        for(k=0;A[P][i]+A[P][j]<=rhs[P];k++){
          if(k&1){
            if(A[P][i]+1<=rhs[P]) A[P][i]+=1;
          }else{
            if(A[P][j]+1<=rhs[P]) A[P][j]+=1;
          }
        }
        //se rompio
        //si inicialmente las dos valian 0, ahora valen cerca de rhs[0]/2. por ahi cuentita andaba....
        
      }else{
        if(random()&1){ // meto x_i en una restriccion con x_j o al reves?
          f1:;
          int mask=(mask1^mask2)&mask1; // ahora la mascara indica donde i tiene algo que j no
          if(mask==0) goto f2;          // quizas, no habia ninguna restriccion asi, pero seguro en alguna hay algo que la otra no tiene
          
          //elijo un bit al azar
          int P=random()%(__builtin_popcount(mask))+1;
          
          for(k=0;P;k++) if(mask&(1<<k)) P--;
          //en k tengo el Pesimo bit (el peor de todos)
          //o sea que voy a poner a j en la kesima restriccion
          //Quiero A[k][j] tal que A[k][i] + A[k][j] > rhs[k]
          // A[k][j]>=1+rhs[k]-A[k][i]
          // A[k][j]<=rhs[k]  para que no se haga infeasible
          
          if(A[k][i]==1) return -1;
          A[k][j] = rhs[k]-random()%(A[k][i]-1);
          
          
        }else{
          f2:;
          int mask=(mask1^mask2)&mask2;
          if(mask==0) goto f1;
          int P=random()%(__builtin_popcount(mask));
          for(k=0;P;k++) if(mask&(1<<k)) P--;
          
          if(A[k][j]==1) return -1;
          A[k][i] = rhs[k]-random()%(A[k][j]-1);
          
        }
        if(A[k][i]+A[k][j]<=rhs[k]) return -1;
      }
    }
  }
  
  //Termina
  //cerr << N << ' ' << R << endl;
  //
  //for(i=0;i<R;i++){
    //for(j=0;j<N;j++) cerr << A[i][j] << ' ' ;
    //cerr << "   " << rhs[i] << endl;
  //}
  return 0;
}

int main(){
  srand(time(0));
  //No tengo idea como hacer esto
  //Trato de hacer algo que arme cliques violables, lo hago muchas veces para tener algo
  //que el arbol completo de B&B sea grande.
  //Primero armo el LP, despues me encargo de imprimirlo
  int cont=0;
  puts("Maximize");
  puts("Subject To");
  int i,j,k;
  for(i=0;i<80;i++){
    while(generarclique()==-1);
    
    for(j=0;j<R;j++){
      for(k=0;k<N;k++)
        if(k==0)
          printf("%d x%d ",A[j][k],cont+k);
        else
          printf("+ %d x%d ",A[j][k],cont+k);
          
      printf(" <= %d\n",rhs[j]);
    }
    
    cont+=N;
  }
  
  printf("obj: ");
  for(i=0;i<cont;i++) printf("%ld x%d + ",random()%10,i);
  puts("");
  
  puts("Binaries");
  for(i=0;i<cont;i++) printf("x%d\n",i);
  puts("End");
}
