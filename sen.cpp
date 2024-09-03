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

int main(void){
    double xMin = 0.0;
    double xMax = M_PI;

    double x = 0;
    double resultado = 0;

    unsigned long muestras = 1000000000;
    long tiempoInicial = 0;
    long tiempoFinal = 0;
    double delta = (xMax -xMin) / muestras;
    
    #pragma parallel for schedule(dinamic,1)
    {
        for(x = xMin; x <= xMax; x += delta){
            resultado += f(x)*delta;
        }
    }

    
    tiempoFinal = getSegundos();
    tiempoFinal -= tiempoInicial;
    printf("La f(x) = %f Se realizo en %ld segundos\n", resultado, tiempoFinal);


}
