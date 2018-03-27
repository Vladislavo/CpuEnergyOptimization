#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef __ENRGOPT_H__
#include "enrgopt.h"
#endif

/* Funcion para calcular la potencia instantanea */
float inst_power(float wbase, int alfa, int beta, float ** voltajes, float ** frecuencias, 
				int * modos, float p_idle, int grado_paralelizacion, float coreC){
	float ipwr = wbase;
	int i = grado_paralelizacion-1, j;

	for(j = 0; j < alfa; j++, i--){
		ipwr += p_din(coreC, (*voltajes)[modos[i]], (*frecuencias)[modos[i]]);
		//printf("modo[%d] %d v = %.2f f = %.2f\n", i, modos[i], (*voltajes)[modos[i]], (*frecuencias)[modos[i]]);
		//printf("ipwr = %.2f\n", p_din(coreC, (*voltajes)[modos[i]], (*frecuencias)[modos[i]]));
	}
	ipwr += beta*p_idle;
	//printf("STOP inst_power\n");

	return ipwr;
}

/* Funcion para calcular la potencia dinamica */
float p_din(float c, float v, float f){
	return c*pow(v,2)*f;
}

/* Funcion que calcula el consiumo de energia dados la potencia y el tiempo (Wh) */
float energy_cons(float p, float time){
	return p*time/3600;
}

int getCharCount(char tok, char* string){
	int i = 0, ret = 0;
	for(; i < strlen(string); i++)
		if(string[i] == tok)
			ret++;
	return ret;
}

/* Funcion que lee los datos de un fichero y devuelve distintos punteros */
void read_data(char* file_conf, int* cores, float* wbase, float* wcoreinactivo, float* coreC, 
				float** voltajes, float** frequencias, int* timepo_secuencial, int* grado_paralelizacion, 
				float** division_trabajo, int *p_estados){
	int read, i;
	FILE *file = fopen(file_conf, "r");

	char *buffer = NULL;
	char *ops, value[128], conf_data[9][128];
	size_t len = 0;

	i = 0;
	while((read = getline(&buffer, &len, file)) != -1){
        /* saltar comentarios */
        if(buffer[0] == '#')
        	continue;

        ops = strtok(buffer,"=");
		while(ops != NULL){
			strcpy(value, ops);
			ops = strtok(NULL, "=");
		}

		strcpy(conf_data[i++], value);
	}

	close(file);

	/* numero de cores */
	*cores = atoi(conf_data[0]);

	*wbase = atof(conf_data[1]);

	*wcoreinactivo = atof(conf_data[2]);

	*coreC = atof(conf_data[3]);

	i = 0;
	*p_estados = getCharCount(';',conf_data[4])+1;
	*frequencias = (float*) malloc(sizeof(float)*(*p_estados));
	ops = strtok(conf_data[4], ";");
	while(ops != NULL){
		(*frequencias)[i] = atof(ops);
		ops = strtok(NULL, ";");
		i++;
	}
	
	i = 0;
	*voltajes = (float*) malloc(sizeof(float)*(*p_estados));
	ops = strtok(conf_data[5], ";");
	while(ops != NULL){
		(*voltajes)[i] = atof(ops);
		ops = strtok(NULL, ";");
		i++;
	}

	*timepo_secuencial = atoi(conf_data[6]);

	*grado_paralelizacion = atoi(conf_data[7]);

	i = 0;
	*division_trabajo = (float*) malloc(sizeof(float)*(*grado_paralelizacion));
	ops = strtok(conf_data[8], ";");
	while(ops != NULL){
		(*division_trabajo)[i] = atof(ops);
		ops = strtok(NULL, ";");
		i++;
	}
}

int comp(const void *a,const void *b) {
	float *x = (float *) a;
	float *y = (float *) b;

	return (*x < *y)? -1 : ((*x > *y)? 1 : 0);

	//if (*x < *y) return -1;
	//else if (*x > *y) return 1; return 0;
}

int count_leq(float elem, float *array, int grado_paralelizacion){
	int res = 0, i;
	for(i = 0; i < grado_paralelizacion; i++){
		if(array[i] <= elem)
			res++;
	}
	return res;
}

float energy_cons_total(int * modos, int grado_paralelizacion, float wbase, float coreC, float ** frecuencias, 
				float ** voltajes, float wcoreinactivo, float * tiempos, int cores){
	float energy, pot_inst, t = 0;
	int c = grado_paralelizacion, i = 0;
	for(i = 0; i < grado_paralelizacion; i++){
		pot_inst = inst_power(wbase, c, (cores-c), voltajes, frecuencias, modos, wcoreinactivo, grado_paralelizacion, coreC);
		energy += energy_cons(pot_inst, (tiempos[i] - t));
		//printf("pot_inst = %.3f, energy = %.3f, alfa = %d, time = %.2f, t = %.2f i = %d\n", pot_inst, energy, c, tiempos[i], t, i);
		c = grado_paralelizacion - count_leq(tiempos[i], tiempos, grado_paralelizacion);
		//printf("count_leq = %d\n", count_leq(tiempos[i], tiempos, grado_paralelizacion));
		i = count_leq(tiempos[i], tiempos, grado_paralelizacion)-1;
		t += tiempos[i] - t;
	}
	return energy;
}

float armax(float * array){
	int i, max = array[0];
	for(i = 1; i < sizeof(array); i++){
		if(array[i] > max)
			max = array[i];
	}
	return max;
}

void sort_states(int * modes, float * time_old, float * time_cur, int grado_paralelizacion){
	int * modes_cpy = (int *) malloc(sizeof(int)*grado_paralelizacion);
	int i, j;
	memcpy(modes_cpy, modes, sizeof(int)*grado_paralelizacion);

	for(i = 0; i < grado_paralelizacion; i++){
		for(j = 0; j < grado_paralelizacion; j++){
			if(time_cur[i] == time_old[j]){
				modes[i] = modes_cpy[j];
			}
		}
	}
	//for(i = 0; i < grado_paralelizacion; i++)
	//	printf("P%d -> P%d\n", modes_cpy[i], modes[i]);

	free(modes_cpy);
}

void modes_optimizer(int * modos, float * tiempos, float ** frecuencias, int p_estados, 
				int grado_paralelizacion){
	float max;
	int found, j, i;
	max = armax(tiempos);
	for(j = 0; j < grado_paralelizacion; j++){
		i = p_estados-1;
		found = 0;
		//printf("tiempo origial = %.2f\n", tiempos[j]);
		if(tiempos[j] < max){
			while(!found && i >= 0){
				//printf("%.3f para P%d t[%d], f[m[%d]] = %.3f, f = %.3f\n", tiempos[j]*(*frecuencias)[modos[j]]/(*frecuencias)[i], j, j, i, (*frecuencias)[modos[j]], (*frecuencias)[i]);
				
				if(tiempos[j]*(*frecuencias)[modos[j]]/(*frecuencias)[i] > max) {// + max*admitancia){
					found = 1;
					tiempos[j] = tiempos[j]*(*frecuencias)[modos[j]]/(*frecuencias)[i+1];
					modos[j] = i + 1;
				//	printf("FOUND %.2f\n", tiempos[j]);
				}
				i--;
			}
			if(!found){
				tiempos[j] = tiempos[j]*(*frecuencias)[modos[j]]/(*frecuencias)[0];
				modos[j] = 0;
			}
		}

		//printf("Modo elejido para core %d es P%d\n", j, modos[j]);
	}
}