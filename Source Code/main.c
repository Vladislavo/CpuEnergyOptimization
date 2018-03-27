#include "enrgopt.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv){
	int i;

	int cores;
	float wbase;
	float wcoreinactivo;
	float coreC;
	float ** frecuencias = (float **) malloc(sizeof(float*));
	float ** voltajes = (float **) malloc(sizeof(float*));
	int tiempo_secuencial;
	int grado_paralelizacion;
	float ** division_trabajo = (float **) malloc(sizeof(float*));
	int p_estados;
	int * modos;

	read_data(argv[1], &cores, &wbase, &wcoreinactivo, &coreC, 
				voltajes, frecuencias, &tiempo_secuencial, &grado_paralelizacion, 
				division_trabajo, &p_estados);

	printf("\t\t Datos leidos desde %s :\n", argv[1]);
	printf("Numero de cores = %u\n", cores);
	printf("Potencia en watios base procesador = %.3f W\n", wbase);
	printf("Potencia en watios core inactivo = %.3f W\n", wcoreinactivo);
	printf("Constante de C por core = %.3f nF\n\n", coreC);

	printf("\tFrecuencias(GHz)\tVoltajes(V)\n");
	for(i = 0; i < p_estados; i++)
		printf("P%d \t %.3f \t\t\t %.3f\n", i, (*frecuencias)[i], (*voltajes)[i]);

	printf("\nDuracion en segundos aplicacion secuencial = %u s\n", tiempo_secuencial);
	printf("Indice de paralelizacion = %d hebras\n", grado_paralelizacion);
	printf("\n\tDivision de trabajo entre hebras: \n");
	for(i = 0; i < grado_paralelizacion; i++)
		printf("Hebra %d - %.2f%%\n", i, (*division_trabajo)[i]);

	modos = (int *) malloc(sizeof(int)*cores);
	for(i = 0; i < cores; i++)
		modos[i] = p_estados-1;

	float * tiempos = (float*) malloc(sizeof(float)*grado_paralelizacion);

	printf("\nTiempos dedicados a las tareas:\n");
	for(i = 0; i < grado_paralelizacion; i++){
		tiempos[i] = tiempo_secuencial*(*division_trabajo)[i]/100;
		printf("Hebra %d %.2f s\n", i, tiempos[i]);
	}
	qsort(tiempos, grado_paralelizacion, sizeof(*tiempos), comp);

	float total_energy = energy_cons_total(modos, grado_paralelizacion, wbase, coreC, frecuencias, 
				voltajes, wcoreinactivo, tiempos, cores);
	
	printf("\nCon ajuste: \n");
	for(i = 0; i < grado_paralelizacion; i++)
		printf("Core %d P%d\n", i, modos[i]);
	printf("el consumo de energia sin optimizacion es de %.3f Wh\n", total_energy);

	modes_optimizer(modos, tiempos, frecuencias, p_estados, grado_paralelizacion);

	printf("\nEl ajuste optimizado:\n");
	for(i = 0; i < grado_paralelizacion; i++)
		printf("Core %d P%d\n", i, modos[i]);

	float * times_old = (float *) malloc(sizeof(float)*grado_paralelizacion);
	memcpy(times_old, tiempos, sizeof(float)*grado_paralelizacion);

	qsort(tiempos, grado_paralelizacion, sizeof(*tiempos), comp);
	
	sort_states(modos, times_old, tiempos, grado_paralelizacion);

	float total_energy_opt = energy_cons_total(modos, grado_paralelizacion, wbase, coreC, frecuencias, 
				voltajes, wcoreinactivo, tiempos, cores);
	printf("con el consumo de energia es de %.3f Wh\n", total_energy_opt);

	printf("\nLa mejora obtenida = %.3f Wh\n", total_energy - total_energy_opt);

	printf("\nTiempos optimizados:\n");
	for(i = 0; i < grado_paralelizacion; i++){
		printf("Hebra %d - %.2f s - P%d\n", i, tiempos[i], modos[i]);
	}

	free(*frecuencias);
	free(*voltajes);
	free(frecuencias);
	free(voltajes);
	free(*division_trabajo);
	free(division_trabajo);
	free(modos);
	free(tiempos);
	free(times_old);

	return 0;
}