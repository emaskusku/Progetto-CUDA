#include <iostream> 
#include <fstream>
#include <cstdlib>
#include <cmath>
#define EPS 1e-6
#include <ctime>

#define N 4
#define THREADS 32

using namespace std;


typedef float dato;

void printmatrice(float *m){
  for(int i=0; i<N; i++){
    cout << endl;
    for(int j=0; j<N+1; j++){
      cout <<"\t"<< m[i*(N+1)+j];
    }
   }

  cout << endl;
}

__global__ void triangola (dato *a_dev){

	__shared__ float coef [THREADS];
	int ncicli= (N+1-blockIdx.x)/THREADS + 1;

	if(blockIdx.x==1){
		//a_dev[3]=100000;
	}
    
       for(int i=blockIdx.x+1; i<N ; i++){
	   if(a_dev[blockIdx.x*(N+1)+blockIdx.x]!=0){
	   for(int j=0; j<ncicli; j++){
	     int index = threadIdx.x + j*THREADS + (N+1)*blockIdx.x+blockIdx.x; //dovrebbe andare bene
	      if(index < N+1){
		coef[threadIdx.x]=a_dev[index]/a_dev[blockIdx.x*(N+1)+blockIdx.x];
		a_dev[i*(N+1)+threadIdx.x+j*THREADS+blockIdx.x]-= coef[threadIdx.x]*a_dev[i*(N+1)+threadIdx.x+j*THREADS+blockIdx.x];
	  }
	}
       }
     }
     __syncthreads();
    }

int main (){

  cudaSetDevice(0);	
	
  clock_t t;
  int dim;
  float *matrice;
  float termini[N]; //termini noti
  float soluzioni [N],soluzionimath[N];
  ifstream GetMatrix;
  ifstream GetTerm;
  ifstream GetSol;
  dato *a_dev;	
  int nblock = (N+THREADS-1)/THREADS;	
  int width = nblock*THREADS;
  size_t size = N*(N+1)*sizeof(dato);
  cudaMalloc((void**)&a_dev,size);

  matrice=new float [N*(N+1)];

  t=clock();

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
    GetTerm>>matrice[i*(N+1)+N];
    for(int j=0; j<N; j++){
      GetMatrix>>matrice[i*(N+1)+j];
    }
  }

  GetMatrix.close();
  GetTerm.close();


  printmatrice(matrice);

  cudaMemcpy(a_dev, matrice, size, cudaMemcpyHostToDevice);
  triangola<<<(dim-1),THREADS>>>(a_dev);
  cudaMemcpy(matrice, a_dev, size, cudaMemcpyDeviceToHost);

  printmatrice(matrice);
	
  //cudaEvent_t start,stop;
  //cudaEventCreate(&start);
  //cudaEventCreate(&stop);
  //cudaEventRecord(start,0);

  


  cudaFree(a_dev);
  return 0;
}


