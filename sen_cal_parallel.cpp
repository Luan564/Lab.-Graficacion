#include <cstdio>
#include <cmath>
#include <ctime>
#include <omp.h>

long getSegundos(){
   return time(NULL);
}

double f(const double &x){
    return sin(x);
}


int main(){
    double xMin = 0.0;
    double xMax = M_PI;

    double x = 0;
    double resultado = 0.0;

    unsigned long muestras = 1000000000;
    unsigned long regionSize = 1000000;
    unsigned long regiones = muestras / regionSize;

    long tiempoInicial = 0;
    long tiempoFinal = 0;
    double delta = (xMax -xMin) / muestras;
    tiempoInicial = getSegundos();

    #pragma omp parallel for

    for(unsigned long i = 0; i < regiones; i++){

        if(i == 0){
            printf("Hilos: %d\n",omp_get_num_threads());
        }

        unsigned long desde = i * regionSize;
        unsigned long hasta = (i + 1) * regionSize;

        double resultadoParcial = 0.0;

        for(unsigned long j  = desde; j < hasta; j++){
            double x = xMin + delta * j;
            resultadoParcial += f(x) * delta;
        }
        #pragma omp atomic
        resultado += resultadoParcial;
    }

    tiempoFinal = getSegundos();
    tiempoFinal -= tiempoInicial;
    printf("El resultado de integrar f(x) en [%f,%f] es: %f en: %ld\n",xMin,xMax,resultado,tiempoFinal);

}