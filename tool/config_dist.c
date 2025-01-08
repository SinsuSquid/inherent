#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define LINESIZE 256
#define MAXTIMESTEP 100000
#define NUMK 6
#define NUMBINS 10000
#define DELTA_T 1

int numTraj;
int timestep[MAXTIMESTEP];
int numAtoms[MAXTIMESTEP];
double box[MAXTIMESTEP][3][2];
int ***atom;
double ***coord;
double binsize;
int numParticles[3];
double box_x, box_y, box_z;

int readTraj(void);
void config_distance(void);

FILE *fp_in;
FILE *fp_out;

int main(int argc, char *argv[]){
	if (argc != 3){
		printf("USAGE : ./scattering.x ***.lammpstrj configurational_distance.out\n");
		exit(1);
	}
	fp_in = fopen(argv[1], "r");
	fp_out = fopen(argv[2], "w");

	atom = (int***)malloc(sizeof(int**) * MAXTIMESTEP);
	coord = (double***)malloc(sizeof(double**) * MAXTIMESTEP);

	readTraj();

	printf("\n\t\t< FILE SUMMARY >\n\n");
	printf("\tNumber of Timesteps : %d\n", numTraj);
	printf("\tNumber of Atoms : %d\n\n", numAtoms[0]);

	box_x = (box[0][0][1] - box[0][0][0]) / 2.0;
	box_y = (box[0][1][1] - box[0][1][0]) / 2.0;
	box_z = (box[0][2][1] - box[0][2][0]) / 2.0;

	binsize = sqrt(box_x * box_x + box_y * box_y + box_z * box_z) * 5.0 / NUMBINS;

	numParticles[0] = 0; numParticles[1] = 0; numParticles[2] = 0;
	for (int i = 0; i < numAtoms[0]; i++){
		switch(atom[0][i][1]){
			case 1: numParticles[0] += 1; break;
			case 2: numParticles[1] += 1; break;
			case 3: numParticles[2] += 1; break;
		}
	}
	printf("box_x = %lf\n", box_x * 2.0);
	printf("box_y = %lf\n", box_y * 2.0);
	printf("box_z = %lf\n", box_z * 2.0);
	printf("numLi = %d\n", numParticles[0]);
	printf("numCl = %d\n", numParticles[1]);
	printf("numAl = %d\n", numParticles[2]);
	printf("\n");

	// Add a function here.
	config_distance();

	free(atom);
	free(coord);

	printf("\tAll Tasks are Done ! >:D\n");

	return 0;
}

void config_distance(void){
	printf("\tNow calculating configurational space distances ...\n");
	fprintf(fp_out, "#\tt\t\\Delta R(t)^I\n");

	for (int i = 0; i < numTraj - 1; i++){
		double distance = 0;
		for (int ii = 0; ii < numAtoms[i]; ii++){
			double dx = coord[i+1][ii][0] - coord[i][ii][0];
			double dy = coord[i+1][ii][1] - coord[i][ii][1];
			double dz = coord[i+1][ii][2] - coord[i][ii][2];

			distance += sqrt(dx * dx + dy * dy + dz * dz);
		}
		fprintf(fp_out, "%f\t%f\n", (double)(i+0.5), distance);
	}
	return;
}

int readTraj(void){
	char *iostat;
	char line[LINESIZE];

	while(1){
		iostat = fgets(line, LINESIZE, fp_in);
		// printf("s", line);

		if (!iostat) break;
		else if (strcmp(line, "ITEM: TIMESTEP\n") == 0){
			int temp;
			fscanf(fp_in, "%d", &temp);
			timestep[numTraj] = temp;
			// printf("timestep : %d\n", timestep[numTraj]);
		}
		else if (strcmp(line, "ITEM: NUMBER OF ATOMS\n") == 0){
			int temp;
			fscanf(fp_in, "%d", &temp);
			numAtoms[numTraj] = temp;
			// printf("numAtoms : %d\n", numAtoms[numTraj]);
		}
		else if (strcmp(line, "ITEM: BOX BOUNDS pp pp pp\n") == 0){
			double temp1, temp2;
			for (int i = 0; i < 3; i++){
				fscanf(fp_in, "%lf %lf", &temp1, &temp2);
				box[numTraj][i][0] = temp1;
				box[numTraj][i][1] = temp2;
			}
			/*
			printf("box :\n%f %f\n%f %f\n%f %f\n",
			        box[numTraj][0][0], box[numTraj][0][1],
			        box[numTraj][1][0], box[numTraj][1][1],
			        box[numTraj][2][0], box[numTraj][2][1]);
			*/
		}
		else if (strcmp(line, "ITEM: ATOMS id type xu yu zu\n") == 0){
			double x, y, z;
			int type, id, ix, iy, iz;

			double boxlength[3] = {box[numTraj][0][1] - box[numTraj][0][0],
					       box[numTraj][1][1] - box[numTraj][1][0],
					       box[numTraj][2][1] - box[numTraj][2][0]};

			int **atomPerTraj;
			atomPerTraj = (int**)malloc(sizeof(int*) * numAtoms[numTraj]);
			double **coordPerTraj;
			coordPerTraj = (double**)malloc(sizeof(double*) * numAtoms[numTraj]);

			for (int i = 0; i < numAtoms[numTraj]; i++){
				int *atomTemp;
				atomTemp = (int*)malloc(sizeof(int) * 2);
				double *coordTemp;
				coordTemp = (double*)malloc(sizeof(double) * 3);

				fscanf(fp_in, "%d %d %lf %lf %lf",
				       &id, &type, &x, &y, &z);

				atomTemp[0] = id;
				atomTemp[1] = type;
				atomPerTraj[i] = atomTemp;

				coordTemp[0] = x;
				coordTemp[1] = y;
				coordTemp[2] = z;
				coordPerTraj[i] = coordTemp;

				/*
				printf("timestep : %d id : %d, type : %d\n",
				       timestep[numTraj], atomPerTraj[i][0], atomPerTraj[i][1]);
				printf("\tx : %f y : %f z : %f\n",
				       coordPerTraj[i][0],
				       coordPerTraj[i][1],
				       coordPerTraj[i][2]);
				*/
			}
			atom[numTraj] = atomPerTraj;
			coord[numTraj] = coordPerTraj;
			numTraj += 1;
		}
	}
	fclose(fp_in);
	return 0;
}
