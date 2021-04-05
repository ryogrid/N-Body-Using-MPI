#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include <mpi.h>
#include "hacked_mpi.h"
#include "definitions.hpp"
#include "iofunc.hpp"
#include "ogl.hpp"

MPI_Datatype mpi_star;

typedef struct Star {
  int m;
  float x;
  float y;
  float sx;
  float sy;
} Star;

void moveGalaxy(Star *galaxy, int nbStars, int id, int split) {

  int i,j;
  long double dx, dy, dst, Fx = 0, Fy = 0, ax, ay;

  for(i = id * split; i < (id + 1) * split && i < nbStars; i++) {
    for(j = 0; j < nbStars; j++) {
      if(i == j) continue;
      //We calculate the direction vector;
      dx = galaxy[j].x - galaxy[i].x;
      dy = galaxy[j].y - galaxy[i].y;
      dst = dx * dx + dy * dy;
      //We calculate the acceleration
      Fx += (/*CST_G * */((double) /*galaxy[i].m * */galaxy[j].m) / dst) * dx;
      Fy += (/*CST_G * */((double) /*galaxy[i].m * */galaxy[j].m) / dst) * dy;
    }
    ax = CST_G * Fx;// / ((double) galaxy[i].m);
    ay = CST_G * Fy;// / ((double) galaxy[i].m);
    //We update the speed
    galaxy[i].sx += ax * DELTA_T;
    galaxy[i].sy += ay * DELTA_T;
    Fx = Fy = 0;
  }

  for(i = id * split; i < (id + 1) * split && i < nbStars; i++) {
    //We update the position
    galaxy[i].x += galaxy[i].sx * DELTA_T;
    galaxy[i].y += galaxy[i].sy * DELTA_T;
  }
}


int main(int c,char **v) {

  if(c != 4 && c != 2) { printf("Usage : program [InputFile] [OutputFile] [Iterations] OR program [OpenGLInputFile].\n"); exit(WRONG_USAGE); }

  int nbIterations = atoi(v[3]);

  if(nbIterations < 1) { printf("Invalid number of iterations.\n"); exit(UNSUPPORTED); }

  MPI_Init(&c,&v);

  int id, size, nbStars;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&id);

  //We define in MPI the new type
  MPI_Datatype types[2] = {MPI_INT, MPI_FLOAT};
  int len[2] = {1,4};
  MPI_Aint offset[2] = {(MPI_Aint)offsetof(Star, m),(MPI_Aint)offsetof(Star, x)};
  MPI_Type_create_struct(2, len, offset, types, &mpi_star);
  MPI_Type_commit(&mpi_star);

  Star *galaxy = NULL;
  FILE *f = NULL;

  if(id == 0) {
    galaxy = loadGalaxy(v[1], &nbStars,NULL);

    //We send the necessary data
    MPI_Bcast(&nbStars,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(galaxy,nbStars,mpi_star,0,MPI_COMM_WORLD);

    f = initStorage(v[2],nbStars, nbIterations);
    storeGalaxy(f,galaxy,nbStars);
  } else {
    MPI_Bcast(&nbStars,1,MPI_INT,0,MPI_COMM_WORLD);
    galaxy = (Star*)malloc(nbStars * sizeof(Star));
    MPI_Bcast(galaxy,nbStars,mpi_star,0,MPI_COMM_WORLD);
  }

  int i, split = nbStars / size;
  for(i = 0; i < nbIterations; i++) {
    moveGalaxy(galaxy,nbStars,id,split);
    if(id == 0)
      moveGalaxy(galaxy,nbStars,size,split);
    MPI_Allgather(&(galaxy[split * id]),split,mpi_star,galaxy,split,mpi_star,MPI_COMM_WORLD);
    if(id == 0)
      storeGalaxy(f,galaxy,nbStars);
  }

  free(galaxy);
  fclose(f);
  MPI_Finalize();
  return EXIT_SUCCESS;

}
