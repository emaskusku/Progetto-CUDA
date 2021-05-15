#include <iostream>
#include <cmath>
#include "mathlink.h"

void PivTriangolazione(double *m, double v[], int N){
  double max=m[0];
  double appo;
  double vappo[N];
  double coef;

  for(int i=0; i<N-1; i++){
    //cout<<"Passaggio "<<i+1<<" / "<<N-1<<endl;
    //seleziono e scambio le righe col massimo
    for(int j=i; j<N; j++){
      if(max<abs(m[j*N+i])){
        max=abs(m[j*N+i]);
        for(int k=i; k<N; k++){
          vappo[k]=m[j*N+k];
          m[j*N+k]=m[i*N+k];
          m[i*N+k]=vappo[k];
        }
        appo=v[j];
        v[j]=v[i];
        v[i]=appo;
      }
    }
    max=m[(i+1)*N+i+1];
    //annullo i coefficienti della i-esima colonna
    if(m[i*N+i]!=0){
      for(int l=i; l<N-1; l++){
        coef=-m[(l+1)*N+i]/m[i*N+i];
        v[l+1]+=coef*v[i];
        for(int q=i; q<N; q++){
          m[(l+1)*N+q]+=coef*m[i*N+q];
        }
      }
    }
  }
}

void solve(double *m, double termini[], double sol[], int N){
    sol[N-1]=termini[N-1]/m[(N-1)*N+N-1];
    for(int i=N-2; i>=0; i--){
      for(int j=N-1; j>i; j--){
        termini[i]-=m[i*N+j]*sol[j];
      }
      if(m[i*N+i]!=0){
        sol[i]=termini[i]/m[i*N+i];
      }
      else{
        sol[i]=0;
      }
    }
}

void solveC (double term[], long dim){
  long   *dimensions;
  char   **heads;
  long   depth;
  double *data;
  double sol[dim];

  MLGetDoubleArray(stdlink, &data, &dimensions, &heads, &depth);

  PivTriangolazione(data,term,dimensions[0]);
  solve (data,term,sol,dim);

  MLPutRealList(stdlink, sol, dim);

  MLDisownDoubleArray(stdlink, data, dimensions, heads, depth);
}
