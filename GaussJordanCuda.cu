#include <iostream> 
#include <fstream>
#include <cstdlib>
#include <cmath>
#define EPS 1e-6
#include <ctime>

#define N 500
#define THREADS 32

using namespace std;


typedef float dato;

void printmatrice(dato *m){
  for(int i=0; i<N; i++){
    cout << endl;
    for(int j=0; j<N+1; j++){
      cout <<"\t"<< m[i*(N+1)+j];
    }
   }

  cout << endl;
}

void confrontaSol(dato sol[N],dato solmath[N]){
  cout<<"Soluzioni con differenza > 0.0001 rispetto al risultato esatto: "<<endl;
  for(int i=0; i<N; i++){
	if(abs(sol[i]/solmath[i]-1)>0.0001)
	cout<<i<<" "<< abs(sol[i]/solmath[i]-1)<<endl;
  }
}

__global__ void completa_matrice(dato *a) {

	int i = blockIdx.y*blockDim.y + threadIdx.y;
	int j = blockIdx.x*blockDim.x + threadIdx.x;
	int offset = gridDim.x*blockDim.x*i + j;
	
	if(i >= N || j >= N+1) 
		a[offset] = 0.;


};



__global__ void triangolo(dato* matric_dev, int index){
	
	int i = blockIdx.y*blockDim.y + threadIdx.y;
	int j = blockIdx.x*blockDim.x + threadIdx.x;
	int offset = gridDim.x*blockDim.x*i + j;

	int row = index*gridDim.x*blockDim.x;
	int ind = row+j;
	//matric_dev[threadIdx.x]=2;

	
	if(i>index && j>index){
		if(matric_dev[row+index]!=0){
			dato coef = matric_dev[i*gridDim.x*blockDim.x+index]/matric_dev[row+index];

			matric_dev[offset] -= matric_dev[ind]*coef;
			//matric_dev[offset]=4;
		}
	}
}



__global__ void findmax (dato* matrice_dev, dato* pivot_dev, int col, int *Max){
	
	__shared__ dato cachecontrol [THREADS];
	__shared__ int indices [THREADS];
	
	int row = blockIdx.x*blockDim.x + threadIdx.x;
	indices[threadIdx.x]=row;

	if(row>=col){
		cachecontrol[threadIdx.x] = abs(matrice_dev[(row)*(gridDim.x*blockDim.x)+col]);
	} else {
		cachecontrol[threadIdx.x] = 0;
	}
	
	__syncthreads();

	int i=THREADS/2;
	if(threadIdx.x<i){
		while (i!=0){
			if(threadIdx.x<i){
				if(cachecontrol[threadIdx.x]<cachecontrol[threadIdx.x+i]){
					cachecontrol[threadIdx.x] = cachecontrol[threadIdx.x+i];
					indices[threadIdx.x]=indices[threadIdx.x+i];
				}
				i/=2;
				__syncthreads();
			} else {
				i=0;
				}
		}
	}

	if(threadIdx.x==0){
		pivot_dev[blockIdx.x]=cachecontrol[0];
		Max[blockIdx.x]=indices[0];
	}

	__syncthreads();

	if(threadIdx.x==0 && blockIdx.x==0){
		for( int i=0; i<gridDim.x; i++){
			if(pivot_dev[0]<pivot_dev[i]){
				pivot_dev[0]=pivot_dev[i];
				Max[0]=Max[i];
			}
		}
	}
}

__global__ void pivoting (dato* matrice_dev, int col, int* rowmax ){ // metti int rowmax
	
	
	int lenrow = gridDim.x*blockDim.x;
	int row = blockIdx.x*blockDim.x + threadIdx.x;

	dato appo = matrice_dev[lenrow*rowmax[0]+row];
	matrice_dev[lenrow*rowmax[0]+row] = matrice_dev[lenrow*col+row];
	matrice_dev[lenrow*col+row]=appo;
}


void solve(dato *m, dato sol[]){
    sol[N-1]=m[N*(N+1)-1]/m[N*(N+1)-2];
    for(int i=N-2; i>=0; i--){
      for(int j=N-1; j>i; j--){
	m[i*(N+1)+N] -= m[i*(N+1)+j]*sol[j];
      }
      if(m[i*(N+1)+i]!=0){
        sol[i]=m[i*(N+1)+N]/m[i*(N+1)+i];
      }
      else{
        sol[i]=0;
      }
   }
}


void triangolizza( dato *matrice ){
	
	dato* matrice_dev; 
	dato *pivot_dev;
	int *Max;
	//int *Maxcpu;
	int nblock1 = N/THREADS + 1; // y, colonne (j)
	int nblock2 = (N+1)/THREADS +1; // x, righe (i)
	int width1= nblock1*THREADS;
	int width2= nblock2*THREADS;
	
	cudaMalloc((void**)&matrice_dev,width1*width2*sizeof(dato));
	cudaMalloc((void**)&pivot_dev,nblock1*sizeof(dato));
	cudaMalloc((void**)&Max,nblock1*sizeof(int));
	//Maxcpu = new int [nblock1];
	
	cout<< "numero blocchi: "<<nblock2<<" x "<< nblock1 <<" width: "<<width2<<" x "<<width1<<endl;

	cudaMemcpy2D(matrice_dev,width2*sizeof(dato),matrice,(N+1)*sizeof(dato),(N+1)*sizeof(dato),N,cudaMemcpyHostToDevice);

	dim3 c_threads(THREADS,THREADS);
	dim3 c_blocks(nblock2,nblock1);

	completa_matrice<<<c_blocks,c_threads>>>(matrice_dev);

	for( int i=0; i<N-1; i++){
		findmax<<<nblock1,THREADS>>>(matrice_dev, pivot_dev, i, Max);
		//cudaMemcpy(Maxcpu,Max,nblock1*sizeof(int),cudaMemcpyDeviceToHost);
		//cout<<Maxcpu[0]<<endl;
		pivoting<<<nblock1,THREADS>>>(matrice_dev, i, Max);
		triangolo<<<c_blocks,c_threads>>>(matrice_dev,i);
	}
	
	cudaMemcpy2D(matrice,(N+1)*sizeof(dato),matrice_dev,width2*sizeof(dato),(N+1)*sizeof(dato),N,cudaMemcpyDeviceToHost);

	cudaFree(pivot_dev);
	cudaFree(matrice_dev);
  }




int main (){

  cudaSetDevice(0);	
	
  //clock_t t;
  int dim;
  dato *matrice;
  dato soluzioni [N], soluzionimath[N];

  ifstream GetMatrix;
  ifstream GetTerm;
  ifstream GetSol;

  matrice=new dato [N*(N+1)];

  //t=clock();

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

  GetSol.open("solutions.txt");
  if(GetSol.fail()){
     cout<< endl << "Problema apertura file di ingresso dati! Esco!";
     return 0;
  }


  for(int i=0; i<N; i++){
	GetSol>> soluzionimath[i];
  }
  

  GetMatrix.close();
  GetTerm.close();


  //printmatrice(matrice);
  triangolizza(matrice);	
  //printmatrice(matrice);
  solve(matrice, soluzioni);

  cout<< "Soluzioni CUDA: "<<endl;
  for(int i=0; i<N; i++){
	cout<<soluzioni[i]<<endl;
  }
  
  cout<< "Soluzioni mathematica: "<<endl;
  for(int i=0; i<N; i++){
	cout<<soluzionimath[i]<<endl;
  }
  
  confrontaSol(soluzioni, soluzionimath);
	
  //cudaEvent_t start,stop;
  //cudaEventCreate(&start);
  //cudaEventCreate(&stop);
  //cudaEventRecord(start,0);

  
  return 0;
}


