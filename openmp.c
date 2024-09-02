#include<stdio.h>
#include<omp.h>

int main(void){

omp_set_num_threads(50);

    #pragma omp parallel
{
    int id = omp_get_thread_num();
    int mt = omp_get_num_threads();

    printf("\n Hilo: %d de un total de: %d\n", id, mt);

}

return 0;
}