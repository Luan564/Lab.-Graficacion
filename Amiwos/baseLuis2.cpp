#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <limits>
#include <iostream>
#include <iomanip>
using namespace std;

class Vector
{
public:
	double x, y, z; // coordenadas x,y,z

	// Constructor del vector, parametros por default en cero
	Vector(double x_ = 0, double y_ = 0, double z_ = 0)
	{
		x = x_;
		y = y_;
		z = z_;
	}

	// operador para suma y resta de vectores
	Vector operator+(const Vector &b) const { return Vector(x + b.x, y + b.y, z + b.z); }
	Vector operator-(const Vector &b) const { return Vector(x - b.x, y - b.y, z - b.z); }
	// operator multiplicacion vector y escalar
	Vector operator*(double b) const { return Vector(x * b, y * b, z * b); }

	// operator % para producto cruz
	Vector operator%(Vector &b) { return Vector(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }
	Vector operator%(const Vector &b) const { return Vector(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }

	// producto punto con vector b
	double dot(const Vector &b) const { return x * b.x + y * b.y + z * b.z; }

	// producto elemento a elemento (Hadamard product)
	Vector mult(const Vector &b) const { return Vector(x * b.x, y * b.y, z * b.z); }

	// normalizar vector
	Vector &normalize() { return *this = *this * (1.0 / sqrt(x * x + y * y + z * z)); }
};
typedef Vector Point;
typedef Vector Color;

class Ray
{
public:
	Point o;
	Vector d;								   // vec y direcccion del rayo
	Ray(Point o_, Vector d_) : o(o_), d(d_) {} // constructor
};

class Sphere
{
public:
	double r; // radius
	Point p;  // position
	Color c;  // color
	Color e;  // radiancia

	Sphere(double r_, Point p_, Color c_, Color e_) : r(r_), p(p_), c(c_), e(e_) {}

	// interseccion del rayo con la esfera
	double intersect(const Ray &ray) const
	{
		Vector op = p - ray.o; //vector de origen desde el centro de la esfera
		double t, epsilon = 1e-4;
		double b = op.dot(ray.d);
		double det = b * b - op.dot(op) + r * r;
		if (det < 0)
			return 0;
		else
			det = sqrt(det);
		return (t = b - det) > epsilon ? t : ((t = b + det) > epsilon ? t : 0);
	}
};

Sphere spheres[] = {
	// Escena: radio, posicion, color, radiancia
	Sphere(1e5, Point(-1e5 - 49, 0, 0), Color(.75, .25, .25), Color()),	  // pared izq
	Sphere(1e5, Point(1e5 + 49, 0, 0), Color(.25, .25, .75), Color()),	  // pared der
	Sphere(1e5, Point(0, 0, -1e5 - 81.6), Color(.25, .75, .25), Color()), // pared detras
	Sphere(1e5, Point(0, -1e5 - 40.8, 0), Color(.25, .75, .75), Color()), // suelo
	Sphere(1e5, Point(0, 1e5 + 40.8, 0), Color(.75, .75, .25), Color()),  // techo
	Sphere(16.5, Point(-23, -24.3, -34.6), Color(.2, .3, .4), Color()),	  // esfera abajo-izq
	Sphere(16.5, Point(23, -24.3, -3.6), Color(.4, .3, .2), Color()),	  // esfera abajo-der
	Sphere(10.5, Point(0, 24.3, 0), Color(1, 1, 1), Color(10, 10, 10))	  // esfera arriba
};

// limita el valor de x a [0,1]
inline double clamp(const double x)
{
	if (x < 0.0)
		return 0.0;
	else if (x > 1.0)
		return 1.0;
	return x;
}

// convierte un valor de color en [0,1] a un entero en [0,255]
inline int toDisplayValue(const double x)
{
	return int(pow(clamp(x), 1.0 / 2.2) * 255 + .5);
}

// PROYECTO 1
/* calcular la intersección del rayo r con todas las esferas
regresar true si hubo una intersección, falso de otro modo
almacenar en t la distancia sobre el rayo en que sucede la interseccion
almacenar en id el indice de spheres[] de la esfera cuya interseccion es mas cercana */

inline bool intersect(const Ray &r, double &t, int &id)
{
	bool interseccion = false;
	double cercana = std::numeric_limits<double>::max();

	for (int i = 0; i < sizeof(spheres) / sizeof(spheres[0]); ++i)
	{
		double actual = spheres[i].intersect(r);

		// verifica si la interseccion es valida y mas cercana que la ultima interseccion
		if (actual > 0 && actual < cercana)
		{
			cercana = actual;
			id = i;
			interseccion = true;
		}
	}

	t = cercana;
	return interseccion;
}

//Obtencion de coordenadas locales
void coordinateSystem(const Vector &n, Vector &s, Vector &t)
{
	if (fabs(n.x) > fabs(n.y))
	{
		float invLen = 1.0f / sqrt(n.x * n.x + n.z * n.z);
		t = Vector(n.z * invLen, 0.0f, -n.x * invLen);
	}
	else
	{
		float invLen = 1.0f / sqrt(n.y * n.y + n.z * n.z);
		t = Vector(0.0f, n.z * invLen, -n.y * invLen);
	}
	s = n % t;
}

Color MonteCarloEsferico(const Ray &r, int muestreos)
{
    double t;
    int id = 0;

    if (!intersect(r, t, id))
        return Color();

    const Sphere &obj = spheres[id];

    if (id == 7)
    {
        return obj.e;
    }

    Color radiancia_acum = Color();

    for (int i = 0; i < muestreos; i++)
    {
		//interseccion
        Point x = r.o + r.d * t;
        Vector n = (x - obj.p).normalize();

        // Genera un sistema de coordenadas locales
        Vector s, tVec;
        coordinateSystem(n, s, tVec);

        double r1 = 2 * M_PI * rand() / (double)RAND_MAX; 
        double r2 = rand() / (double)RAND_MAX; 
        double r2s = sqrt(r2); 

        // Direcciones ortogonales en el sistema de coordenadas local
        Vector w = n;
        Vector u = ((fabs(w.x) > .1 ? Vector(0, 1) : Vector(1)) % w).normalize();
        Vector v = w % u;

        // Dirección del muestreo en coordenadas esféricas
        Vector direc = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).normalize();

        double t2;
        int id2 = 0;

        if (!intersect(Ray(x, direc), t2, id2))
            continue;

        const Sphere &obj2 = spheres[id2];

        // Radiancia de la esfera intersectada
        Color radiancia = obj2.e;
        // BRDF difusa
        Color brdf = obj.c * (1 / M_PI);
        // Coseno del ángulo entre la normal y la dirección del muestreo
        double coseno = n.dot(direc);

        // Calcula la radiancia para este muestreo
        Color radiancia_muestra = radiancia.mult(brdf) * coseno;

        // Acumula la radiancia del muestreo
        radiancia_acum = radiancia_acum + radiancia_muestra;
    }

    return radiancia_acum * (1.0 / muestreos);
}


Color MonteCarloHemisferico(const Ray &r, int muestreos) {
    double t;
    int id = 0;

    if (!intersect(r, t, id)) return Color(); 
    
    const Sphere &obj = spheres[id];

    if (id == 7) {
        return obj.e; 
    }

    Point x = r.o + r.d * t;
    Vector n = (x - obj.p).normalize(); 

    // Genera un sistema de coordenadas locales
    Vector s, tVec;
    coordinateSystem(n, s, tVec);

    Color color_acumulado = Color();
	
    for (int i = 0; i < muestreos; i++) {

        double r1 = rand() / (double)RAND_MAX;
        double r2 = rand() / (double)RAND_MAX;

        // Calcula los ángulos theta y phi para el muestreo hemisférico
        double theta = acos(sqrt(r1));
        double phi = 2 * M_PI * r2; 

        // Calcula la dirección del muestreo en el sistema de coordenadas local
        Vector direc = s * (sin(theta) * cos(phi)) + tVec * (sin(theta) * sin(phi)) + n * cos(theta);

        double t2;
        int id2 = 0;

        if (!intersect(Ray(x, direc), t2, id2)) continue; 

        const Sphere &obj2 = spheres[id2];

        // Radiancia de la esfera intersectada
        Color radiancia = obj2.e;

        // BRDF difusa
        Color brdf = obj.c * (1 / M_PI);

        // Coseno del ángulo entre la normal y la dirección del muestreo
        double coseno = n.dot(direc);

        color_acumulado = color_acumulado + radiancia.mult(brdf) * coseno;
    }

    return color_acumulado * (1.0 / muestreos);
}


Color shadeE(const Ray &r, int muestreos)
{
	return MonteCarloEsferico(r, muestreos);
}

Color shadeH(const Ray &r, int muestreos)
{
	return MonteCarloHemisferico(r, muestreos);
}

int main(int argc, char *argv[])
{
int opcion;
#pragma omp parallel
	opcion = 0;
	do
	{
		int opcion;
		std::cout << "\nMENU:" << std::endl;
		std::cout << "1. Muestreo uniforme esferico" << std::endl;
		std::cout << "2. Muestreo uniforme hemisferico" << std::endl;
		std::cout << "3. Salir" << std::endl;
		std::cout << "Elija una opcion: ";
		std::cin >> opcion;

		switch (opcion)
		{
		case 1:
		{
			int w = 1024, h = 768; // image resolution

			// fija la posicion de la camara y la dirección en que mira
			Ray camera(Point(0, 11.2, 214), Vector(0, -0.042612, -1).normalize());

			// parametros de la camara
			Vector cx = Vector(w * 0.5095 / h, 0., 0.);
			Vector cy = (cx % camera.d).normalize() * 0.5095;

			// auxiliar para valor de pixel y matriz para almacenar la imagen
			Color *pixelColors = new Color[w * h];

			int muestreos;
			std::cout << "\nIngresa el numero de muestras: ";
			std::cin >> muestreos;

			if (muestreos <= 0)
			{
				std::cout << "Elige un numero de muestreos valido.\n" << std::endl;
				break;
			}

			// PROYECTO 1
			// usar openmp para paralelizar el ciclo: cada hilo computara un renglon (ciclo interior),
			
			#pragma omp parallel for
			for (int y = 0; y < h; y++)
			{
				// recorre todos los pixeles de la imagen
				fprintf(stderr, "\r%5.2f%%", 100. * y / (h - 1));
				for (int x = 0; x < w; x++)
				{
					int idx = (h - y - 1) * w + x; // index en 1D para una imagen 2D x,y son invertidos
					Color pixelValue = Color();	   // pixelValue en negro por ahora
					// para el pixel actual, computar la dirección que un rayo debe tener
					Vector cameraRayDir = cx * (double(x) / w - .5) + cy * (double(y) / h - .5) + camera.d;

					// computar el color del pixel para el punto que intersectó el rayo desde la camara
					pixelValue = shadeE(Ray(camera.o, cameraRayDir.normalize()), muestreos);

					// limitar los tres valores de color del pixel a [0,1]
					pixelColors[idx] = Color(clamp(pixelValue.x), clamp(pixelValue.y), clamp(pixelValue.z));
					
				}
			}

			fprintf(stderr, "\n");

			// PROYECTO 1
			// Investigar formato ppm
			FILE *f = fopen("image.ppm", "w");
			// escribe cabecera del archivo ppm, ancho, alto y valor maximo de color
			fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
			for (int p = 0; p < w * h; p++)
			{ // escribe todos los valores de los pixeles
				fprintf(f, "%d %d %d ", toDisplayValue(pixelColors[p].x), toDisplayValue(pixelColors[p].y),
						toDisplayValue(pixelColors[p].z));
			}
			fclose(f);

			delete[] pixelColors;
			break;
		}
		case 2:
		{
			int w = 1024, h = 768; // image resolution

			// fija la posicion de la camara y la dirección en que mira
			Ray camera(Point(0, 11.2, 214), Vector(0, -0.042612, -1).normalize());

			// parametros de la camara
			Vector cx = Vector(w * 0.5095 / h, 0., 0.);
			Vector cy = (cx % camera.d).normalize() * 0.5095;

			// auxiliar para valor de pixel y matriz para almacenar la imagen
			Color *pixelColors = new Color[w * h];

			int muestreos2;
			std::cout << "\nIngresa el numero de muestras: ";
			std::cin >> muestreos2;

			if (muestreos2 <= 0)
			{
				std::cout << "Elige un numero de muestreos valido.\n" << std::endl;
				break;
			}

			// PROYECTO 1
			// usar openmp para paralelizar el ciclo: cada hilo computara un renglon (ciclo interior),
			#pragma omp parallel for
			for (int y = 0; y < h; y++)
			{
				// recorre todos los pixeles de la imagen
				fprintf(stderr, "\r%5.2f%%", 100. * y / (h - 1));
				for (int x = 0; x < w; x++)
				{
					int idx = (h - y - 1) * w + x; // index en 1D para una imagen 2D x,y son invertidos
					Color pixelValue = Color();	   // pixelValue en negro por ahora
					// para el pixel actual, computar la dirección que un rayo debe tener
					Vector cameraRayDir = cx * (double(x) / w - .5) + cy * (double(y) / h - .5) + camera.d;

					// computar el color del pixel para el punto que intersectó el rayo desde la camara
					pixelValue = shadeH(Ray(camera.o, cameraRayDir.normalize()), muestreos2);

					// limitar los tres valores de color del pixel a [0,1]
					pixelColors[idx] = Color(clamp(pixelValue.x), clamp(pixelValue.y), clamp(pixelValue.z));
				}
			}

			fprintf(stderr, "\n");

			// PROYECTO 1
			// Investigar formato ppm
			FILE *f = fopen("image.ppm", "w");
			// escribe cabecera del archivo ppm, ancho, alto y valor maximo de color
			fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
			for (int p = 0; p < w * h; p++)
			{ // escribe todos los valores de los pixeles
				fprintf(f, "%d %d %d ", toDisplayValue(pixelColors[p].x), toDisplayValue(pixelColors[p].y),
						toDisplayValue(pixelColors[p].z));
			}
			fclose(f);

			delete[] pixelColors;
			break;
		}
		case 3:
			std::cout << "Saliendo del programa...\n" << std::endl;
			opcion = 3;
			return 0;
		default:
			std::cout << "Elige una opcion correcta...\n" << std::endl;

			std::cin.clear();
			std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			break;
		}
	} while (opcion != 3);

	return 0;
}
