#include <cstdio>
#include <omp.h>

int main(void){
    #pragma omp parallel
    {
        #pragma omp single //Investigar este pragma
        {
            printf("Solo entramos una vez a el hilo: %d\n", omp_get_thread_num());
        }
        printf("Thread %d\n",omp_get_thread_num());
    }



}