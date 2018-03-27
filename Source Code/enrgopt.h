/* Funciones para un optimizador de energia */
#ifndef __ENRGOPT_H__
#define __ENRGOPT_H__

/* Funcion para calcular la potencia instantanea */
float inst_power(float wbase, int alfa, int beta, float ** voltajes, float ** frecuencias, 
				int * modos, float p_idle, int cores, float coreC);

/* Funcion para calcular la potencia dinamica */
float p_din(float c, float v, float f);

/* Funcion que lee los datos de un fichero de configuracion y devuelve distintos punteros 
	1) numero de cores
	2) el vector voltajes
	3) el vector de frecuencias
	4) numero de p-estados
	5) tiempo de ejecucion secuencial
	6) grado de paralelizacion
	7) division de trabajo en cada hebra
	8) numero de p_estados
*/
void read_data(char* file_conf, int* cores, float* wbase, float* wcoreinactivo, float* coreC, 
				float** voltajes, float** frequencias, int* timepo_secuencial, int* grado_paralelizacion, 
				float** division_trabajo, int* p_estados);

//float energy_cons(float p, float time);

/* Funcion que calcula el aumento de tiempo en funcion del tiempo base */
float * slow_down_time(float * frecs, float t_base);

/* Funcion para facilitar la lectura de datos */
int getCharCount(char tok, char* string);

/* Funcion de comparacion para usar el quicksort() */
int comp(const void *a,const void *b);

/* Funcion que cuenta cuantos elementos son menores o iguales en que un elem en array */
int count_leq(float elem, float *array, int grado_paralelizacion);

/* Funcion que calcula el consumo total de energia dados los paramentros */
float energy_cons_total(int * modos, int grado_paralelizacion, float wbase, float coreC, float ** frequencias, 
				float ** voltajes, float wcoreinactivo, float * tiempos, int cores);

/* Funcion que busca el valor maximo de un array */
float armax(float * array);

/* Funcion que ordena los p-estado en funcion de los tiempos optimizados */
void sort_states(int * modes, float * time_old, float * time_cur, int grado_paralelizacion);

/* Funcion que optimiza los modos de funcionamiento de los CPUs */
void modes_optimizer(int * modos, float * tiempos, float ** frecuencias, int p_estados, 
				int grado_paralelizacion);

#endif