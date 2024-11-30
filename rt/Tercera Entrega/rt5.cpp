// rt: un lanzador de rayos minimalista
 // g++ -O3 -fopenmp rt.cpp -o rt
#include <math.h>
#include <stdlib.h>
#include <stdio.h>  
#include <omp.h>
#include <time.h>

class Vector {
public:        
	double x, y, z; // coordenadas x,y,z 
  
	// Constructor del vector, parametros por default en cero
	Vector(double x_= 0, double y_= 0, double z_= 0){ x=x_; y=y_; z=z_; }
  
	// operador para suma y resta de vectores
	Vector operator+(const Vector &b) const{ 
		return Vector(x + b.x, y + b.y, z + b.z); 
	}
	Vector operator-(const Vector &b) const{ 
		return Vector(x - b.x, y - b.y, z - b.z); 
	}

	// operator multiplicacion vector y escalar 
	Vector operator*(double b) const{ 
		return Vector(x * b, y * b, z * b); 
	}
  
	// operator % para producto cruz
	Vector operator%(const Vector &b){
		return Vector(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
	}
	
	// producto punto con vector b
	double dot(const Vector &b) const{ 
		return x * b.x + y * b.y + z * b.z; 
	}

	// producto elemento a elemento (Hadamard product)
	Vector mult(const Vector &b) const{ 
		return Vector(x * b.x, y * b.y, z * b.z); 
	}
	
	// normalizar vector 
	Vector& normalize(){ 
		return *this = *this * (1.0 / sqrt(x * x + y * y + z * z)); 
	}

	Vector cross(const Vector& other) const {
    return Vector(
        y * other.z - z * other.y,  // Componente x
        z * other.x - x * other.z,  // Componente y
        x * other.y - y * other.x   // Componente z
    );
}
};
typedef Vector Point;
typedef Vector Color;

class Ray { 
public:
	Point o;
	Vector d; // origen y direcccion del rayo
	Ray(Point o_, Vector d_) : o(o_), d(d_) {} // constructor
};

class Sphere {
public:
	double r;	// radio de la esfera
	Point p;	// posicion
	Color c;	// color  
	Color emLuz;	//Si emite o no Luz

	Sphere(double r_, Point p_, Color c_, Color L_): r(r_), p(p_), c(c_), emLuz(L_) {}
  
	// PROYECTO 1
	// [1]	determina si el rayo intersecta a esta esfera
	// [2]	regresar distancia si hay intersección
	// [3]	regresar 0.0 si no hay interseccion
	double intersect(const Ray &ray) const {
		double discriminante = 0;

		Vector aux = ray.o.operator-(p);					//Vector    (o-p)
		double a = ray.d.dot(ray.d);						//double	(d.d)
		double b = 2.0*(aux.dot(ray.d));					//double	2*[(o-p).d]
		double c = aux.dot(aux) - (r*r); 					//double	(o-p).(o-p)-r²
		discriminante = (b*b)-4.0*(a*c);

		if(discriminante < 0){
			return 0.0;
		}
		else{
			double t = (-b - sqrt(discriminante)) / (2.0 * a);
        	if(t > 0) {
            	return t;
        	}
			else{
            	return 0; 
        	}
		}
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
inline double clamp(const double x) { 
	if(x < 0.0)
		return 0.0;
	else if(x > 1.0)
		return 1.0;
	return x;
}

// convierte un valor de color en [0,1] a un entero en [0,255]
inline int toDisplayValue(const double x) {
	return int( pow( clamp(x), 1.0/2.2 ) * 255 + .5); 
}

// PROYECTO 1
// [1]	calcular la intersección del rayo r con todas las esferas
// [2]	regresar true si hubo una intersección, falso de otro modo
// [3]	almacenar en t la distancia sobre el rayo en que sucede la interseccion
// [4]	almacenar en id el indice de spheres[] de la esfera cuya interseccion es mas cercana
inline bool intersect(const Ray &r, double &t, int &id) {
	bool interseccion = false;
    double cercana = std::numeric_limits<double>::max(); 	//Numero maximo para un double 

    for (int i = 0; i <= 8; ++i) {
        double actual = spheres[i].intersect(r); 

        // verifica si la interseccion es valida y mas cercana que la ultima interseccion
        if(actual > 0 && actual < cercana){
            cercana = actual; 
            id = i; 
            interseccion = true; 
        }
    }
    t = cercana; 
    return interseccion; 
}

Vector localesToGlobales(const Vector &n, Vector &s, Vector &t, Vector &local){
	Vector global{ 
		(s.x * local.x + t.x * local.y + n.x * local.z), 
		(s.y * local.x + t.y * local.y + n.y * local.z),
		(s.z * local.x + t.z * local.y + n.z * local.z)
	};

	return global;
}

Vector globalesToLocales(const Vector &n, Vector &s, Vector &t, Vector &global){
	Vector local{
		(s.x * global.x + s.y * global.y + s.z * global.z), 
		(t.x * global.x + t.y * global.y + t.z * global.z),
		(n.z * global.x + n.z * global.y + n.z * global.z)
	};

	return local;
}

// Calcula el valor de color para el rayo dado
Color shade(const Ray &r) {
	double t;
	int id = 0;
	// determinar que esfera (id) y a que distancia (t) el rayo intersecta
	if (!intersect(r, t, id)){
		return Color();	// el rayo no intersecto objeto, return Vector() == negro
	}
	const Sphere &obj = spheres[id];
	
	// PROYECTO 1
	// determinar coordenadas del punto de interseccion
	Point x;
	x = r.o + r.d.mult(t);		// Origen del rayo + (t * direccion)

	// determinar la dirección normal en el punto de interseccion

	Vector n;		//Definiremos el vector desde el centro de la esfera hasta el punto de intersección calculado
	n = x.operator-(spheres[id].p);		//n = Punto De Intersección - Origen de la Esfera
	n.normalize();						//Normalizamos "n" con su función ya definida

	// determinar el color que se regresara
	Color colorValue;
	colorValue = spheres[id].c;

	return colorValue; 
	
}


class FuenteArea {
public:
    Point centro;    // Centro de la fuente
    Vector u, v;     // Vectores para definir el área (dimensiones y orientación, solo si es necesario)
    Color emLuz;     // Radiancia emitida por la fuente
    double radio;    // Radio del disco

    FuenteArea(Point c_, double radio_, Color L_) 
        : centro(c_), radio(radio_), emLuz(L_) {}

    // Genera un punto aleatorio dentro del disco
    Point randomPoint() const {
        // Genera un ángulo aleatorio en [0, 2*pi]
        double angulo = 2 * M_PI * (rand() / (double)RAND_MAX);
        
        // Genera una distancia aleatoria dentro del radio (rango [0, radio])
        double distancia = sqrt(rand() / (double)RAND_MAX) * radio; 

        // Coordenadas polares a cartesianas para el punto en el disco
        double x = centro.x + distancia * cos(angulo);
        double z = centro.z + distancia * sin(angulo);

        // El disco está alineado con el plano XZ, por lo que la coordenada Y permanece constante
        return Point(x, centro.y, z);
    }
};

// Instancia de fuente de área en lugar de la esfera luminosa (ahora un disco)
FuenteArea fuenteLuminosa(
    Point(0, 23.4, -10),    // Centro del disco
    15.0,                 // Radio del disco
    Color(10, 10, 10)     // Radiancia
);

Color monteCarloArea(const Ray &r, int muestreos) {
    double t;
    int id = 0;
    Color radiancia;

    // Intersección del rayo con la escena
    if (!intersect(r, t, id)) return Color(); // No hay intersección
    if (id == 7) return spheres[id].emLuz; // Fuente de luz visible

    const Sphere &obj = spheres[id];
    Point x = r.o + (r.d * t); // Punto de intersección
    Vector n = (x - obj.p).normalize(); // Normal en el punto de intersección

    // BRDF simplificada
    Vector brdf = obj.c * (1.0 / M_PI);

    for (int i = 0; i < muestreos; i++) {
        Point lightPoint = fuenteLuminosa.randomPoint(); // Punto aleatorio en el disco
        Vector dir_to_light = (lightPoint - x).normalize();

        double cos_theta = n.dot(dir_to_light);
        if (cos_theta <= 0) continue; // Ignorar si está en dirección opuesta

        Ray shadowRay(x, dir_to_light);
        double t2;
        int id2;

        // Verificar si el rayo alcanza la fuente sin obstrucciones
        if (!intersect(shadowRay, t2, id2) || id2 == 7) { // Asegurarse de que los rayos lleguen correctamente
            radiancia = radiancia + (fuenteLuminosa.emLuz.mult(brdf) * cos_theta);
        }
    }

    // Normalización de la radiancia acumulada con un factor de atenuación
    return radiancia * (1.0 / muestreos) * 0.1; 
}
 

// Define una fuente de ángulo sólido
class FuenteAnguloSolido {
public:
    Point posicion;       // Posición de la fuente
    Vector direccion;     // Dirección principal de la fuente
    double thetaMax;      // Ángulo máximo en radianes
    Color emision;        // Color (radiancia) emitida

    FuenteAnguloSolido(Point p, Vector d, double maxTheta, Color L)
        : posicion(p), direccion(d.normalize()), thetaMax(maxTheta), emision(L) {}

    // Genera una dirección dentro del ángulo sólido
    Vector muestrearDireccion() const {
        double r1 = rand() / (double)RAND_MAX; // Aleatorio entre 0.0 y 1.0
        double r2 = rand() / (double)RAND_MAX;

        // Muestreo dentro del cono definido por thetaMax
        double cosTheta = 1.0 - r1 * (1.0 - cos(thetaMax));
        double sinTheta = sqrt(1.0 - cosTheta * cosTheta);
        double phi = 2.0 * M_PI * r2;

        // Coordenadas locales
        Vector localDir(sinTheta * cos(phi), sinTheta * sin(phi), cosTheta);

        // Marco de referencia local a global
        Vector u = (fabs(direccion.x) > 0.1 ? Vector(0, 1, 0) : Vector(1, 0, 0)) % direccion;
        u = u.normalize();
        Vector v = direccion.cross(u);

        // Transformar a coordenadas globales
        return Vector(
            u.x * localDir.x + v.x * localDir.y + direccion.x * localDir.z,
            u.y * localDir.x + v.y * localDir.y + direccion.y * localDir.z,
            u.z * localDir.x + v.z * localDir.y + direccion.z * localDir.z
        ).normalize();
    }
};

// Ejemplo de inicialización de una fuente de ángulo sólido
FuenteAnguloSolido fuente(
	Point(0, 24.3, 0), 
	Vector(0, 0, 0), M_PI, 
	Color(10, 10, 10));

// Función de cálculo con fuente de ángulo sólido
Color monteCarloAnguloSolido(const Ray &r, int muestreos, double theta_max) {
    double t;
    int id = 0;
    Color colorValue;
    Color radiancia;

    if (!intersect(r, t, id))
        return Color(); // El rayo no intersectó ningún objeto, return negro.

    const Sphere &obj = spheres[id];
    if (id == 7)
        return obj.emLuz; // La fuente emite luz directamente.

    Point x = r.o + (r.d * t);         // Primer punto de intersección.
    Vector n = (x - obj.p).normalize(); // Vector normal a la superficie.

    // BRDF
    Vector brdf = obj.c * (1.0 / M_PI);

    // Muestreo en ángulo sólido
    for (int z = 0; z < muestreos; z++) {
        double r1 = rand() / (double)RAND_MAX;
        double r2 = rand() / (double)RAND_MAX;

        double theta = acos(1 - r1 * (1 - cos(theta_max))); // Restricción en theta.
        double phi = 2 * M_PI * r2;

        // Conversión esféricas a cartesianas para la dirección muestreada.
        Vector dir_local(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
        dir_local = dir_local.normalize();

        // Marco de referencia local a partir del vector normal `n`.
        float invLen = 1.0f / sqrt(n.x * n.x + n.z * n.z);
        Vector te = Vector(n.z * invLen, 0.0f, -n.x * invLen);
        Vector s = te % n;

        Vector dir_global = localesToGlobales(n, s, te, dir_local);
        dir_global = dir_global.normalize();

        // Coseno theta
        double cos_theta = n.dot(dir_global);

        // Probabilidad
        double prob = 1.0 / (2 * M_PI * (1 - cos(theta_max))); // PDF para ángulo sólido.

        // Rayo con origen fijo y dirección muestreada.
        Ray ray2(x, dir_global);

        // Segunda intersección.
        int id2 = 0;
        double t2;
        intersect(ray2, t2, id2);
        const Sphere &obj2 = spheres[id2];
        if (id2 == 7) { // Intersecta la fuente de luz.
            Vector emision = obj2.emLuz;
            radiancia = radiancia + ((emision.mult(brdf)) * cos_theta * (1.0 / prob));
        }
        if (z == muestreos - 1) {
            colorValue = radiancia * (1.0 / muestreos);
        }
    }
    return colorValue;
}

int main(int argc, char *argv[]) {
	clock_t inicio; 
	int op = 0;
	int muestreos = 0;
	int w = 1024, h = 768; // image resolution
  
	// fija la posicion de la camara y la dirección en que mira
	Ray camera( Point(0, 11.2, 214), Vector(0, -0.042612, -1).normalize() );

	// parametros de la camara
	Vector cx = Vector( w * 0.5095 / h, 0., 0.); 
	Vector cy = (cx % camera.d).normalize() * 0.5095;
  
	// auxiliar para valor de pixel y matriz para almacenar la imagen
	Color *pixelColors = new Color[w * h];

		printf("------ Menú -----\n");
		printf("1.- Muestreo con fuente de área\n");
        printf("2.- Muestreo Angulo Solido\n");
		printf("3.- Salir\n");
		printf("Seleccione una opción: ");
		scanf("%d", &op);

		if (op == 3) {
			return 0; 
		}





		printf("Numero de muestras: ");
		scanf("%d",&muestreos);


		switch (op){
		case 1:{
			#pragma omp parallel for
            for (int y = 0; y < h; y++) { 
                fprintf(stderr, "\r%5.2f%%", 100. * y / (h - 1));
                for (int x = 0; x < w; x++) {
                    int idx = (h - y - 1) * w + x;
                    Vector cameraRayDir = cx * (double(x) / w - .5) + cy * (double(y) / h - .5) + camera.d;
                    pixelColors[idx] = monteCarloArea(Ray(camera.o, cameraRayDir.normalize()), muestreos);
                }
            }
		}break;

		case 2: { // Fuente de ángulo sólido
			#pragma omp parallel for
			for (int y = 0; y < h; y++) { 
				fprintf(stderr, "\r%5.2f%%", 100. * y / (h - 1));
				for (int x = 0; x < w; x++) {
					int idx = (h - y - 1) * w + x;
					Vector cameraRayDir = cx * (double(x) / w - .5) + cy * (double(y) / h - .5) + camera.d;
					double theta_max = M_PI; // Definir un valor para el ángulo máximo.
					Ray cameraRay(camera.o, cameraRayDir.normalize());
					pixelColors[idx] = monteCarloAnguloSolido(cameraRay, muestreos, theta_max);
				}
			}
		} break;


		default:
			break;
		}





	// // PROYECTO 1
	// // usar openmp para paralelizar el ciclo: cada hilo computara un renglon (ciclo interior),
	// #pragma omp parallel for
	// for(int y = 0; y < h; y++) { 
	// 	// recorre todos los pixeles de la imagen
	// 	fprintf(stderr,"\r%5.2f%%",100.*y/(h-1));
	// 	for(int x = 0; x < w; x++ ) {
	// 		int idx = (h - y - 1) * w + x; // index en 1D para una imagen 2D x,y son invertidos
	// 		Color pixelValue = Color(); // pixelValue en negro por ahora
	// 		// para el pixel actual, computar la dirección que un rayo debe tener
	// 		Vector cameraRayDir = cx * ( double(x)/w - .5) + cy * ( double(y)/h - .5) + camera.d;
			
	// 		// computar el color del pixel para el punto que intersectó el rayo desde la camara
	// 		pixelValue = shadeU( Ray(camera.o, cameraRayDir.normalize()), muestreos );

	// 		// limitar los tres valores de color del pixel a [0,1]
	// 		pixelColors[idx] = Color(clamp(pixelValue.x), clamp(pixelValue.y), clamp(pixelValue.z));
	// 	}
	// }

	fprintf(stderr,"\n");

	// PROYECTO 1
	// Investigar formato ppm
	FILE *f = fopen("image3.ppm", "w");
	// escribe cabecera del archivo ppm, ancho, alto y valor maximo de color
	fprintf(f, "P3\n%d %d\n%d\n", w, h, 255); 
	for (int p = 0; p < w * h; p++) 
	{ // escribe todos los valores de los pixeles
    		fprintf(f,"%d %d %d ", toDisplayValue(pixelColors[p].x), toDisplayValue(pixelColors[p].y), 
				toDisplayValue(pixelColors[p].z));
  	}
  	fclose(f);

  	delete[] pixelColors;

	// Obtener el tiempo de fin
    clock_t fin = clock();

    // Calcular la duración en segundos
    double duracion = (double)(fin - inicio) / CLOCKS_PER_SEC;

    printf("Tiempo de ejecución: %.6f segundos\n", duracion);

	return 0;
}
