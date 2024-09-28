#include <iostream>
#include <cstdlib>
#include <cmath>

class circulo {
private:
    double a, b;
    unsigned int N;
    unsigned int adentro;
    unsigned int afuera;
    double resultado;

public:
    circulo(double a, double b, unsigned int N) : a(a), b(b), N(N), adentro(0), resultado(0) {}

    double f(double x, double y) {
        return (x*x + y*y);  // La función f(x)   
    }

    double getResultado() {
        return resultado;
    }

    double Xi() {
        return 1 - static_cast<double>(rand()) / RAND_MAX * (b - a);  // Generar un xi entre a y b
    }

    double Yi(double M) {
        return 1 - static_cast<double>(rand()) / RAND_MAX * M;
    }

    double calcularM() {
        // El valor máximo de la función en el intervalo [a, b]. En este caso, f(x) = x^5 alcanza su máximo en x = 3
        return f(b,a);
    }

    void integrar() {
        double M = calcularM();  // Se obtiene el valor máximo de la función
        for (unsigned int i = 0; i < N; i++) {
            double xi = Xi();
            double yi = Yi(M);
            double fx = f(xi,yi);

            if (fx <= 1) {
                adentro++;
            }
            else{
                afuera++;
            }
        }
        resultado = (b - a) * calcularM() * static_cast<double>(adentro) / N;  // Estimación del área
    }

    int getAdentro() {
        return adentro;
    }
    int getAfuera() {
        return afuera;
    }

};

int main() {
    circulo mc(-1, 1, 1000000);  
    mc.integrar();
    std::cout << "Puntos adentro son: " << mc.getAdentro() << std::endl;
    std::cout << "Puntos afuera son: " << mc.getAfuera() << std::endl;
    std::cout << "El resultado de la integral es: " << mc.getResultado() << std::endl;
    return 0;
}
