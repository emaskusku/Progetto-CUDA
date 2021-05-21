#include "mathlink.h"

void addtwoM(void) {
        long i, j, sum;
        MLGetLongInteger(stdlink, &i);
        MLGetLongInteger(stdlink, &j);
        sum = i + j;
        MLPutLongInteger(stdlink, sum);
}

void main (int argc, char *argv[]){

	addtwoM();
	return;
}
