// rt: un lanzador de rayos minimalista
 // g++ fuente.cpp -o compilado -fopenmp
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>  
#include <omp.h>
#include <cstdlib> 
#include <ctime>

using namespace::std;

const double PI =  3.14159265358979323846;
#define N 32    //Número de muestras para el montecarlo 


class Vector 
{
public:        
	double x, y, z; // coordenadas x,y,z 
  
	// Constructor del vector, parametros por default en cero
	Vector(double x_= 0, double y_= 0, double z_= 0){ x=x_; y=y_; z=z_; }
  
	// operador para suma y resta de vectores
	Vector operator+(const Vector &b) const { return Vector(x + b.x, y + b.y, z + b.z); }
	Vector operator-(const Vector &b) const { return Vector(x - b.x, y - b.y, z - b.z); }
	// operator multiplicacion vector y escalar 
	Vector operator*(double b) const { return Vector(x * b, y * b, z * b); }
  
	// operator % para producto cruz
	Vector operator%(Vector&b){return Vector(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);}
	
	// producto punto con vector b
	double dot(const Vector &b) const { return x * b.x + y * b.y + z * b.z; }

	// producto elemento a elemento (Hadamard product)
	Vector mult(const Vector &b) const { return Vector(x * b.x, y * b.y, z * b.z); }
	
	// normalizar vector 
	Vector& normalize(){ return *this = *this * (1.0 / sqrt(x * x + y * y + z * z)); }
};
typedef Vector Point;   //Un punto tiene 3 coordenadas al igual que un vector
typedef Vector Color;   //El color es un vector rgb de 3 valores al igual que el vector

class Ray 
{ 
public:
	Point o;
	Vector d; // origen y direcccion del rayo
	Ray(Point o_, Vector d_) : o(o_), d(d_) {} // constructor
};

class Sphere 
{
public:
	double r;	// radio de la esfera
	Point p;	// posicion
	Color c;	// color  
    Color l;

	Sphere(double r_, Point p_, Color c_, Color l_): r(r_), p(p_), c(c_), l(l_) {}
  
	// PROYECTO 1
	// determina si el rayo intersecta a esta esfera
	double intersect(const Ray &ray) const {
		// regresar distancia si hay intersección
		// regresar 0.0 si no hay interseccion       
        double dis;
        double terminos[3]={0,0,0};     //Almacena los terminos dentro de la raiz        
        double t1,t2;
        //dot no es un operador, no se sobrecarga, se define
        terminos[0] = pow( (ray.o - this -> p).dot(ray.d) , 2);
        terminos[1] = ( (ray.o - this -> p).dot(ray.o - this -> p) * -1);
        terminos[2] = pow(this -> r,2);
        dis = terminos[0] + terminos[1] + terminos[2];
        if(dis > 0)
        {
            //Se calcula la menor distancia de intersección t1 y t2
            t1 = ( (ray.o - this -> p).dot(ray.d)  * -1 ) + sqrt(dis);
            t2 = ( (ray.o - this -> p).dot(ray.d)  * -1 ) - sqrt(dis);
            if(t1 > t2) 
                return t2;
            else
                return t1;
        }
        else
        {
		    return 0.0;     //En caso de que el rayo no intercepte la esfera
	    }
    }
};

Sphere spheres[] = {
	//Escena: radio, posicion, color 
	Sphere(1e5,  Point(-1e5 - 49, 0, 0),   Color(.75, .25, .25), Color()), // pared izq
	Sphere(1e5,  Point(1e5 + 49, 0, 0),    Color(.25, .25, .75), Color()), // pared der
	Sphere(1e5,  Point(0, 0, -1e5 - 81.6), Color(.25, .75, .25), Color()), // pared detras
	Sphere(1e5,  Point(0, -1e5 - 40.8, 0), Color(.25, .75, .75), Color()), // suelo
	Sphere(1e5,  Point(0, 1e5 + 40.8, 0),  Color(.75, .75, .25), Color()), // techo
	Sphere(16.5, Point(-23, -24.3, -34.6), Color(.2, .3, .4), Color()), // esfera abajo-izq
	Sphere(16.5, Point(23, -24.3, -3.6),   Color(.4, .3, .2), Color()), // esfera abajo-der
	Sphere(10.5, Point(0, 24.3, 0),        Color(1, 1, 1), Color(10,10,10)) // esfera arriba
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
// calcular la intersección del rayo r con todas las esferas
// regresar true si hubo una intersección, falso de otro modo
// almacenar en t la distancia sobre el rayo en que sucede la interseccion
// almacenar en id el indice de spheres[] de la esfera cuya interseccion es mas cercana

inline bool intersect(const Ray &r, double &t, int &id) {

    double array_t[8] = {0,0,0,0,0,0,0,0};  //Para almacenar las distancias a las esferas para despues ordenarlas
    double t_min;
    int t_min_index;
    int interseccion = 0;
    //Se calcula la intersección de r con las esferas segun el discriminante
    for(int i = 0; i < 8; i++)
    {
        array_t[i] = spheres[i].intersect(r);
        if(array_t[i] != 0)
        {   
            interseccion++; //El rayo interceptó alguna esfera
        }
    }
    //Se busca la distancia más pequeña de las intercepciones
    t_min = INFINITY;   //Para que la primera distancia > 0 se tome la menor sin importar su valor
    t_min_index = 0;
    for(int i = 0; i < 8; i++)
    {
        //Se salta las esferas que no son interseptadas
        if(array_t[i] > 0 && array_t[i] < t_min)
        {
            t_min = array_t[i];
            t_min_index = i;
        }
    }

    if(interseccion != 0)
    {
        //cout<<"La esfera mas cercana es: "<<t_min_index<<" que esta a: "<<t_min<<endl;
        id = t_min_index;   //Esfera más cercana
        t = t_min;
        return true;
    }
    else    //La esfera no tocó o rozó la esfera
    {
        return false;
    }
}

// Calcula el valor de color para el rayo dado
Color shade(const Ray &r) {
    double t;
    int id = 0;
    Color colorValue;   //Color del pixel actual, falta ve como la fuente modifica ese color
    Color radiancia;
    // determinar que esfera (id) y a que distancia (t) el rayo intersecta
    if (!intersect(r, t, id))
        return Color();	// el rayo no intersecto objeto, return Vector() == negro
    else
    {
        const Sphere &obj = spheres[id];
        if(id == 7)
        {
            colorValue = obj.l;
        }
        else
        {
            //MONTE CARLOS
            for(int z = 0; z < N; z++)
            {
                Point x = r.o + (r.d * t);  //primer punto intersección 
                Vector n =  (x - obj.p).normalize();    //Vector normal a la superficie

                //BRDF
                Vector brdf = obj.c * (1.0/ M_PI);

                //Se lanza un segundo rayo desde el punto de intersección ¿Con qué dirección?

                //muestreo uniforme esférico 

                double aleatorio1 = drand48();   //Aleatorio entre 0.0 y 1.0
                double theta = acos(1 - (2 * aleatorio1));
                double aleatorio2 = drand48();
                double phi = 2 * M_PI * aleatorio2;

                //Conversión esféricas-cartesianas para la dirección muestreada
                Vector dir_local(sin(theta)*cos(phi) , sin(theta)*sin(phi), cos(theta) );
                dir_local = dir_local.normalize();

                //Marco de referencia local a partir del vector n
                float invLen = 1.0f / sqrt(n.x * n.x + n.z * n.z);
                Vector te = Vector(n.z * invLen, 0.0f, -n.x * invLen);
                Vector s = te % n;


                //componentes vector global (local-global)
                double c1_x = s.x*dir_local.x + te.x*dir_local.y + n.x*dir_local.z;
                double c2_y =s.y*dir_local.x + te.y*dir_local.y + n.y*dir_local.z; 
                double c3_z =s.z*dir_local.x + te.z*dir_local.y + n.z*dir_local.z; 
                Vector dir_global(c1_x, c2_y ,c3_z);
                dir_global = dir_global.normalize();

                //Coseno theta
                double cos_theta = n.dot(dir_global);
                
                //Probabilidad
                double prob = 1.0 / (4 * M_PI); 

                //Rayo con origen fijo y dirección muestreada para calcular la segundo intersección
                Ray ray2(x, dir_global);

                //Se calcula la segunda intersección
                int id2 = 0;
                double t2;
                intersect(ray2, t2, id2);
                const Sphere obj2 = spheres[id2];   //Segunda esfera con la que intercepta
                if(id2 == 7)    //Conecta con la fuente
                {
                    Vector emision = obj2.l;
                    radiancia = radiancia + ( ((emision.mult(brdf)) * cos_theta) * (1.0 / prob) );
                }
                if(z == N-1)
                {
                    colorValue = radiancia * (1.0/N);
                }
            }   //Llave del for
        }
    }
    return colorValue;
}


int main(int argc, char *argv[]) {
	int w = 1024, h = 768; // image resolution
  
    srand48(time(nullptr));     //Planta la semilla
	// fija la posicion de la camara y la dirección en que mira
	Ray camera( Point(0, 11.2, 214), Vector(0, -0.042612, -1).normalize() );

	// parametros de la camara
	Vector cx = Vector( w * 0.5095 / h, 0., 0.); 
	Vector cy = (cx % camera.d).normalize() * 0.5095;
  
	// auxiliar para valor de pixel y matriz para almacenar la imagen
	Color *pixelColors = new Color[w * h];

	// PROYECTO 1
	// usar openmp para paralelizar el ciclo: cada hilo computara un renglon (ciclo interior),
	for(int y = 0; y < h; y++) 
	{ 
		// recorre todos los pixeles de la imagen
		fprintf(stderr,"\r%5.2f%%",100.*y/(h-1));
        #pragma parallel for 
        for(int x = 0; x < w; x++ ) {
			int idx = (h - y - 1) * w + x; // index en 1D para una imagen 2D x,y son invertidos
			Color pixelValue = Color(); // pixelValue en negro por ahora
			// para el pixel actual, computar la dirección que un rayo debe tener
			Vector cameraRayDir = cx * ( double(x)/w - .5) + cy * ( double(y)/h - .5) + camera.d;
		
			// computar el color del pixel para el punto que intersectó el rayo desde la camara
            pixelValue = shade( Ray(camera.o, cameraRayDir.normalize()) );
			// limitar los tres valores de color del pixel a [0,1]
			pixelColors[idx] = Color(clamp(pixelValue.x), clamp(pixelValue.y), clamp(pixelValue.z));
		}
	}

	fprintf(stderr,"\n");

	// PROYECTO 1
	// Investigar formato ppm
	FILE *f = fopen("image.ppm", "w");
	// escribe cabecera del archivo ppm, ancho, alto y valor maximo de color
	fprintf(f, "P3\n%d %d\n%d\n", w, h, 255); 
	for (int p = 0; p < w * h; p++) 
	{ // escribe todos los valores de los pixeles
    		fprintf(f,"%d %d %d ", toDisplayValue(pixelColors[p].x), toDisplayValue(pixelColors[p].y), 
				toDisplayValue(pixelColors[p].z));
  	}
  	fclose(f);

  	delete[] pixelColors;

	return 0;
}
