#include "mathlink.h"

double prodscalar(double v1[], long dim1, double v2[], long dim2){

	double prod=0;
	for(int i=0; i<dim1; i++){
		prod+=v1[i]*v2[i];
	}
	//MLPutReal(stdlink, prod);
	return prod;
}

int main(int argc, char* argv[]) {
        return MLMain(argc, argv);
}
