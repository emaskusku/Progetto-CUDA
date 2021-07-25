#include <iostream> 
#include <fstream>
#include <cstdlib>
#include <cmath>
#define EPS 1e-6
#include <ctime>

#define N 7000
#define THREADS 32
#define THREADS1 1024

using namespace std;


typedef double dato;

void printmatrice(dato *m){
  for(int i=0; i<N; i++){
    cout << endl;
    for(int j=0; j<N+1; j++){
      cout <<"\t"<< m[i*(N+1)+j];
    }
   }

  cout << endl;
}

void confrontaSol(dato *sol,dato *solmath){
  cout<<"Soluzioni con differenza > 0.00001 rispetto al risultato esatto: "<<endl;
  for(int i=0; i<N; i++){
	if(abs(sol[i]/solmath[i]-1)>0.00001)
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
	__syncthreads();
}



__global__ void findmax (dato* matrice_dev, dato* pivot_dev, int col, int *Max){
	
	int row = blockIdx.x*blockDim.x + threadIdx.x;
	int len = int((N+1)/THREADS+1)*THREADS;

	__shared__ dato cachecontrol [THREADS1];
	__shared__ int indices [THREADS1];
	
	indices[threadIdx.x]=row;

	if(row>=col && row < len){
	cachecontrol[threadIdx.x] = abs(matrice_dev[(row)*(len)+col]);
	} else {
		cachecontrol[threadIdx.x] = 0;
	}
	
	__syncthreads();

	int i=blockDim.x/2;

	while (i!=0){
		if(threadIdx.x<i){
		if(cachecontrol[threadIdx.x]<cachecontrol[threadIdx.x+i]){
			cachecontrol[threadIdx.x] = cachecontrol[threadIdx.x+i];
			indices[threadIdx.x]=indices[threadIdx.x+i];
		}
		}

		i/=2;
		__syncthreads();

	}

	if(threadIdx.x==0){
		pivot_dev[blockIdx.x]=cachecontrol[0];
		Max[blockIdx.x]=indices[0];
	}

	__syncthreads();
}


__global__ void getmax(dato * pivot_dev, int * Max, int nblock3){

	if(threadIdx.x==0 && blockIdx.x==0){
		for( int i=0; i<nblock3; i++){
			if(pivot_dev[0]<pivot_dev[i]){
				pivot_dev[0]=pivot_dev[i];
				Max[0]=Max[i];
			}
		}
	}
	__syncthreads();
}


__global__ void findmax1(dato *matrice_dev, dato* pivrow, int col, int *Maxind){ // da mandare con 1 blocco e THREADS*nblock threads

	Maxind[threadIdx.x]=threadIdx.x;

	if(threadIdx.x>=col){
		pivrow[threadIdx.x] = matrice_dev[threadIdx.x*blockDim.x+col];
	} else {
		pivrow[threadIdx.x] = 0;
	}
	
	__syncthreads();

	int i=blockDim.x/2;

	while (i!=0){
		if(threadIdx.x<i){
		if(pivrow[threadIdx.x] < pivrow[threadIdx.x+i]){
			pivrow[threadIdx.x] = pivrow[threadIdx.x+i];
			Maxind[threadIdx.x] = Maxind[threadIdx.x+i];
		}
		}

		i/=2;
		__syncthreads();

	}
}


__global__ void pivoting (dato* matrice_dev, dato* pivot_dev, int col, int* Max ){

	if(pivot_dev[0]==0){
		__syncthreads();
		return;
	}

	int lenrow = int((N+1)/THREADS+1)*THREADS;
	int row = blockIdx.x*blockDim.x + threadIdx.x;

	if(row < int((N+1)/THREADS+1)*THREADS){
		dato appo = matrice_dev[lenrow*Max[0]+row];
		matrice_dev[lenrow*Max[0]+row] = matrice_dev[lenrow*col+row];
		matrice_dev[lenrow*col+row]=appo;
		//Max[blockIdx.x]=0;
		//pivot_dev[blockIdx.x]=0;
	}
	__syncthreads();
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


__global__ void solveCuda(dato*m, dato *sol, int block1){

    int len = block1*THREADS;
    sol[N-1]=m[(N-1)*(len)+N]/m[(N-1)*(len)+N-1];

    for(int i=N-2; i>=0; i--){
      for(int j=N-1; j>i; j--){
	m[i*(len)+N] -= m[i*(len)+j]*sol[j];
      }
      if(m[i*(len)+i]!=0){
        sol[i]=m[i*(len)+N]/m[i*(len)+i];
      }
      else{
        sol[i]=0;
      }
   }

   
}


/*
__global__ void triangolocuda(int nblock1, dato* matrice_dev, dato* pivot_dev, int* Max){

	dim3 c_threads(THREADS,THREADS);
	dim3 c_blocks(nblock2,nblock1);
	
	for( int i=0; i<N-1; i++){
		findmax<<<nblock1,THREADS>>>(matrice_dev, pivot_dev, i, Max);
		getmax<<<1,THREADS>>>(pivot_dev, Max, nblock1);
		pivoting<<<nblock1,THREADS>>>(matrice_dev, pivot_dev, i, Max);
		triangolo<<<c_blocks,c_threads>>>(matrice_dev,i);
		cudaDeviceSynchronize();
	}
}
*/

dato* triangolizza( dato *matrice ){
	
	dato* matrice_dev; 
	dato *pivot_dev;
	dato* solcuda;
	dato *solcpu;
	int *Max;
	int nblock1 = N/THREADS + 1; // y, colonne (j)
	int nblock2 = (N+1)/THREADS + 1; // x, righe (i)
	int nblock3 = (N+1)/(THREADS1) + 1;
	nblock1=nblock2;
	int width1= nblock1*THREADS;
	int width2= nblock2*THREADS;
	
	solcpu = new dato [N];
	cudaMalloc((void**)&matrice_dev,width2*width1*sizeof(dato));
	cudaMalloc((void**)&pivot_dev,nblock3*sizeof(dato));
	cudaMalloc((void**)&Max,nblock3*sizeof(int));
	cudaMalloc((void**)&solcuda,N*sizeof(dato));

	cudaEvent_t start,stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start,0);
	
	cout<< "numero blocchi kernel 2D: "<<nblock2<<" x "<< nblock1 <<", width: "<<width2<<" x "<<width1<<endl;

	cudaMemcpy2D(matrice_dev,width2*sizeof(dato),matrice,(N+1)*sizeof(dato),(N+1)*sizeof(dato),N,cudaMemcpyHostToDevice);


	dim3 c_threads(THREADS,THREADS);
	dim3 c_blocks(nblock2,nblock1);

	completa_matrice<<<c_blocks,c_threads>>>(matrice_dev);

	//triangolocuda<<<1,1>>>(nblock1, matrice_dev, pivot_dev, Max);
	cout<< "numero blocchi kernel 1D: "<<nblock3<<", width: "<<THREADS1*nblock3<<endl;

	for( int i=0; i<N-1; i++){

		findmax<<<nblock3,THREADS1>>>(matrice_dev, pivot_dev, i, Max);
		//cudaDeviceSynchronize();
		getmax<<<1,1>>>(pivot_dev, Max, nblock3);
		//cudaDeviceSynchronize();
		pivoting<<<nblock3,THREADS1>>>(matrice_dev, pivot_dev, i, Max);
		//cudaDeviceSynchronize();
		triangolo<<<c_blocks,c_threads>>>(matrice_dev,i);
		//cudaDeviceSynchronize();
	}

	solveCuda<<<1,1>>>(matrice_dev, solcuda, nblock1);

	cudaMemcpy(solcpu, solcuda, N*sizeof(dato), cudaMemcpyDeviceToHost);
	//cudaMemcpy2D(matrice,(N+1)*sizeof(dato),matrice_dev,width2*sizeof(dato),(N+1)*sizeof(dato),N,cudaMemcpyDeviceToHost);


	//printmatrice(matrice);

	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	float elapsed;
	cudaEventElapsedTime(&elapsed,start,stop);
	cout<<"tempo: " << elapsed/1000. << endl;
	
	//dealloco 
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	cudaFree(solcuda);
	cudaFree(pivot_dev);
	cudaFree(matrice_dev);
	return solcpu;
  }




int main (){

  cudaSetDevice(2);	
	
  //clock_t t;
  int dim;
  dato *matrice;
  //dato soluzioni [N];
  dato soluzionimath[N];
  dato* solfromcuda;
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

  GetSol.open("solutions.txt"); //occhio da quale file prendi
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
  solfromcuda = triangolizza(matrice);	
  //printmatrice(matrice);
  //solve(matrice, soluzioni);

/*
  cout<< "Soluzioni CUDA: "<<endl;
  for(int i=0; i<N; i++){
	cout<<soluzioni[i]<<endl;
  }
  
  cout<< "Soluzioni mathematica: "<<endl;
  for(int i=0; i<N; i++){
	cout<<soluzionimath[i]<<endl;
  }
 */

  cout<<solfromcuda[0]<<" "<<soluzionimath[0]<<endl;
  cout<<solfromcuda[1]<<" "<<soluzionimath[1]<<endl;
  cout<<solfromcuda[2]<<" "<<soluzionimath[2]<<endl;
  cout<<solfromcuda[3]<<" "<<soluzionimath[3]<<endl; 
  cout<<solfromcuda[4]<<" "<<soluzionimath[4]<<endl;

/*
  cout<<soluzioni[0]<<" "<<soluzionimath[0]<<endl;
  cout<<soluzioni[1]<<" "<<soluzionimath[1]<<endl;
  cout<<soluzioni[2]<<" "<<soluzionimath[2]<<endl;
  cout<<soluzioni[3]<<" "<<soluzionimath[3]<<endl; 
  cout<<soluzioni[4]<<" "<<soluzionimath[4]<<endl;
*/
  //confrontaSol(soluzioni, soluzionimath);
  confrontaSol(solfromcuda, soluzionimath);
	
  //cudaEvent_t start,stop;
  //cudaEventCreate(&start);
  //cudaEventCreate(&stop);
  //cudaEventRecord(start,0);

  
  return 0;
}

