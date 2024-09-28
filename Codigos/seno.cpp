/*
Debemos encontrar el area bajo la curva del seno

I = Integral[ f(x)dx = (1/N)*(SUMATORIA de i = 0 hasta N-1)*(f(xi)/p(xi))   ]


xi = a + rand*(b-a)

p(xi) = 1/(b-a)

(Función de distribución de probabilidad) PDF

(Area del cuadrado) A = base*altura = 1

xi es tomada aleatoriamente del dominio D (de 0 a pi), y tiene una probabilidad p(xi)

*/


#include <iostream>
#include <cstdlib>
#include <cmath>

class MonteCarlo {
private:
    double a, b;
    double I;
    int n;

public:
    MonteCarlo(double a, double b, int n) : a(a), b(b), n(n){}

    double f(double xi) {
        return sin(xi);  // La función f(x)   
    }

    double getResultado() {
        return I;
    }

    double Xi() {
        return a + static_cast<double>(rand()) / RAND_MAX * (b - a);  // Generar un xi entre a y b
    }

    double pdf(){
        return (1/ (b-a));
    }

    double integrar() {
        
        for(unsigned int i = 0; i<n; i++){
            double xi = Xi();
            double px = pdf();
            I += (f(xi)/px);
        }


        I = I/n;
        return I;
    }
};

int main() {
    MonteCarlo mc(0, M_PI, 100000000);  
    double res = mc.integrar();
    printf("El resultado de la integral es: %f\n", res);

    return 0;
}
