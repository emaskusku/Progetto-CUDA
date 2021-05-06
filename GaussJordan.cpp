#include <iostream>
#include <cmath>

#define N 4


using namespace std;

void printmatrice(double m[N][N], double v[N]){
  for(int i=0; i<N; i++){
    cout << endl;
    for(int j=0; j<N; j++){
      cout <<"\t"<< m[i][j];
    }
    cout<<"\t\t"<<v[i];
  }

  cout << endl;
}

void Triangolazione(double m[N][N], double v[N]){
  double coef;

  for(int i=0; i<N-1; i++){
    for(int k=0; k<i+1; k++){
      if(m[k][k]!=0){
        coef=m[i+1][k]/m[k][k];
          for(int j=0; j<N; j++){
            m[i+1][j]+=-m[k][j]*coef;
          }
          v[i+1]+=-coef*v[k];
        }
      }
    }
  }


void PivTriangolazione(double m[N][N], double v[N]){
  double max=m[0][0];
  double appo;
  double vappo[N];
  double coef;

  for(int i=0; i<N-1; i++){
    //seleziono e scambio le righe col massimo
    for(int j=i; j<N; j++){
      if(max<abs(m[j][i])){
        max=abs(m[j][i]);
        for(int k=i; k<N; k++){
          vappo[k]=m[j][k];
          m[j][k]=m[i][k];
          m[i][k]=vappo[k];
        }
        appo=v[j];
        v[j]=v[i];
        v[i]=appo;
      }
    }
    max=m[i+1][i+1];
    //annullo i coefficienti della i-esima colonna
    if(m[i][i]!=0){
      for(int l=i; l<N-1; l++){
        coef=-m[l+1][i]/m[i][i];
        v[l+1]+=coef*v[i];
        for(int q=i; q<N; q++){
          m[l+1][q]+=coef*m[i][q];
        }
      }
    }
  }
}




int main(){
  double matrice [N][N]={8,1,0,1,4,5,5,7,9,11,7,5,7,4,6,8};
  double termini[N]={2,3,6,80}; //termini noti

  printmatrice(matrice,termini);
  cout <<endl<<"Triangolazione:"<<endl;

  //Triangolazione(matrice,termini);
  PivTriangolazione(matrice,termini);
  printmatrice(matrice,termini);


  return 0;
}
