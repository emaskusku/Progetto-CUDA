#include <iostream>
#include <cmath>
#include <fstream>
#include <ctime>

#define N 7000


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

void printSol(double v[N]){
  cout << "\nSoluzioni: "<<endl;
  for(int i=0; i<N; i++){
    cout<< v[i]<<endl;
  }
}

double determinant( double** matrix, int n) {
   int det = 0;
   double **submatrix;
   submatrix=new double *[N];
   for(int i=0; i<N; i++){
     submatrix[i]=new double [N];
   }
   if (n == 2)
   return ((matrix[0][0] * matrix[1][1]) - (matrix[1][0] * matrix[0][1]));
   else {
      for (int x = 0; x < n; x++) {
         int subi = 0;
         for (int i = 1; i < n; i++) {
            int subj = 0;
            for (int j = 0; j < n; j++) {
               if (j == x)
               continue;
               submatrix[subi][subj] = matrix[i][j];
               subj++;
            }
            subi++;
         }
         det = det + (pow(-1, x) * matrix[0][x] * determinant( submatrix, n - 1 ));
      }
   }
   return det;
}


void Triangolazione(double **m, double v[N]){
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


void PivTriangolazione(double **m, double v[N]){
  double max=m[0][0];
  double appo;
  double vappo[N];
  double coef;

  for(int i=0; i<N-1; i++){
    //cout<<"Passaggio "<<i+1<<" / "<<N-1<<endl;
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

void solve(double **m, double termini[N], double sol[N]){
    sol[N-1]=termini[N-1]/m[N-1][N-1];
    for(int i=N-2; i>=0; i--){
      for(int j=N-1; j>i; j--){
        termini[i]-=m[i][j]*sol[j];
      }
      if(m[i][i]!=0){
        sol[i]=termini[i]/m[i][i];
      }
      else{
        sol[i]=0;
      }
    }
}

void confrontaSol(double sol[N],double solmath[N]){
  cout<<"Soluzioni con differenza >= 0.00001 rispetto al risultato esatto: "<<endl;
  for(int i=0; i<N; i++){
	//cout<<i<<" "<< abs(sol[i]/solmath[i]-1)<<endl;
	if(abs(sol[i]/solmath[i]-1)>=0.00001)
	cout<<i<<" "<< abs(sol[i]/solmath[i]-1)<<endl;
  }
}

int main(){

  clock_t t;
  int dim;
  //double matrice [N][N];
  double termini[N]; //termini noti
  double soluzioni [N],soluzionimath[N];
  ifstream GetMatrix;
  ifstream GetTerm;
  ifstream GetSol;
  double **matrice;

  

  matrice=new double *[N];
  for(int i=0; i<N; i++){
    matrice[i]=new double [N];
  }

  GetMatrix.open("matrix.txt");
  if(GetMatrix.fail()){
     cout<< endl << "Problema apertura file di ingresso dati! Esco!";
     return 0;
  }

  GetTerm.open("term.txt");
  if(GetTerm.fail()){
     cout<< endl << "Problema apertura file di ingresso dati! Esco!";
     return 0;
  }

  GetMatrix>>dim;

  for(int i=0; i<N; i++){
    GetTerm>>termini[i];
    for(int j=0; j<N; j++){
      GetMatrix>>matrice[i][j];
    }
  }

  GetMatrix.close();
  GetTerm.close();


  t=clock();

  //printmatrice(matrice,termini);
  //cout <<endl<<"Triangolazione:"<<endl;

  Triangolazione(matrice,termini);
  //PivTriangolazione(matrice,termini);
  //printmatrice(matrice,termini);
  //cout << endl;
  //cout<<N<<endl;
  //cout<< "Faccio il determinante!"<<endl;
  //if(determinant(matrice,N)!=0){
    //cout << "Determinante > 0, ora risolvo!"<<endl;
  solve(matrice,termini,soluzioni);
  //printSol(soluzioni);
  //}
  //else{
    //cout << "Il sistema non ha un unica soluzione, non lo risolvo!"<<endl;
  //}

  t=clock()-t;

  cout<<"---------------------------------"<<endl;

  GetSol.open("solutions.txt");
  if(GetSol.fail()){
     cout<< endl << "Problema apertura file di ingresso dati! Esco!";
     return 0;
  }
  
  for(int i=0; i<N; i++){
     GetSol >> soluzionimath[i];
  }  

  GetSol.close();


  ofstream PrintSol;

  PrintSol.open("solcpp.txt");
  for( int i=0; i<N; i++){
	PrintSol<<soluzioni[i]<<endl;
  }

  
  //cout<< "Sol math "<<endl;
  //printSol(soluzionimath);
  //cout << "Sol cpp "<<endl;
  //printSol(soluzioni);

  //cout<<"Confronto la soluzione trovata con quella esatta: "<<endl;

  confrontaSol(soluzioni,soluzionimath);
  cout<<soluzioni[0]<<" "<<soluzionimath[0]<<endl;
  cout<<soluzioni[1]<<" "<<soluzionimath[1]<<endl;
  cout<<soluzioni[2]<<" "<<soluzionimath[2]<<endl;
  cout<<soluzioni[3]<<" "<<soluzionimath[3]<<endl; 
  cout<<soluzioni[4]<<" "<<soluzionimath[4]<<endl;
  
  cout << "Tempo impiegato: "<<((float)t)/CLOCKS_PER_SEC<<" sec"<<endl;


  for(int i=0; i<N; i++){
    delete [] matrice[i];
  }

  delete [] matrice;

  return 0;
}
