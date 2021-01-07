#include <stdio.h>
#include <mpi.h>
#include <vector>

int main(){
    
    const double * a;
    const double * b;
    const double * c;
    const double * d;

    double response = predicates::incircle(a,b,c,d);
    return 0;
}