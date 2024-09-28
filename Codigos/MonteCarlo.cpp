#include <iostream>
#include <cstdlib>
#include <cmath>

class MonteCarlo {
private:
    double a, b;
    unsigned int N;
    unsigned int aciertos;
    double resultado;

public:
    MonteCarlo(double a, double b, unsigned int N) : a(a), b(b), N(N), aciertos(0), resultado(0) {}

    double f(double x) {
        return (x*x*x*x) + 1;  // La función f(x)   
    }

    double getResultado() {
        return resultado;
    }

    double Xi() {
        return a + static_cast<double>(rand()) / RAND_MAX * (b - a);  // Generar un xi entre a y b
    }

    double Yi(double M) {
        return static_cast<double>(rand()) / RAND_MAX * M;
    }

    double calcularM() {
        return f(b);
    }

    void integrar() {
        double M = calcularM();  // Se obtiene el valor máximo de la función
        for (unsigned int i = 0; i < N; i++) {
            double xi = Xi();
            double yi = Yi(M);
            double fx = f(xi);

            if (yi <= fx) {
                aciertos++;
            }
        }
        resultado = (b - a) * calcularM() * static_cast<double>(aciertos) / N;  // Estimación del área
    }
};

int main() {
    MonteCarlo mc(1, 3, 1000);  
    mc.integrar();
    printf("El resultado de la integral es: %f\n", mc.getResultado());
    return 0;
}
