#include <cstdio>
#include <omp.h>

#define MAX 10

int count = 0;

int main(void){
    int hilo = 0;
    #pragma omp parallel num_threads(MAX)
    {
        hilo = omp_get_thread_num();
        #pragma omp atomic  //investigara este pragma
        //#pragma omp barrier  //investigara este pragma
        count++;
        printf("Secuencia: %d hilo: %d\n",count, hilo);
    }
    return 0;
}
