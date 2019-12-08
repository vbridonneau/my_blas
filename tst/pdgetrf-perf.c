#include <mpi.h>
#include <stdlib.h>
#include "dgetrf.h"
#include "util.h"

void check_pdgetrf_split() {
  int M = 9, N = 9;
  for (int Mb = 0; Mb < 10; ++Mb) {
  double *A = alloc_matrix(M, N);
  rnd_matrix_buff(A, 1, 10, M * N, 1);

  /* Datatypes */
  MPI_Datatype band_type, last_band_type;

  /* Datatype creation */
  MPI_Type_vector( block_size, M, M, MPI_DOUBLE, &band_type );
  MPI_Type_vector( N % block_size, M, M, MPI_DOUBLE, &last_band_type );

  /* Commit types */
  MPI_Type_commit( &band_type );
  MPI_Type_commit( &last_band_type );

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  /* Info for bandes per process */
  int Ndb              = (N + block_size - 1) / block_size;
  int NT               = (Ndb + size - 1) / size;
  int rlast            = ((Ndb % size) + size - 1) % size;
  int is_last          = (rank == rlast);
  int band_size        = block_size * M;
  int last_band_size   = (N % block_size) * M;
  int is_last_complete = ((N % block_size) == 0);

  /* Matrices */
  int NT_treated = (rank <= rlast) ? NT : NT - 1;
  int matrix_size;
  if (is_last && !is_last_complete) {
    matrix_size = (NT_treated - 1) * band_size + last_band_size;
  } else {
    matrix_size = NT_treated * band_size;
  }
  double *Asub = malloc(sizeof(double) * matrix_size);

  /* Scatter */
  my_pdgetrf_scatter(M, N, A, M, Asub, &band_type, &last_band_type);

  /* Check everyone has the correct bands */
  double *my_A = A + rank * band_size;
  for (int band = 0; band < NT_treated; ++band) {
    double *shiftAsub = Asub + band * band_size;
    double *shiftA    = my_A + band * band_size * size;
    for (int j = 0; j < block_size; ++j) {
      int is_ok = 1;
      for (int i = 0; i < M; ++i) {
        if (!eq_double(shiftA[j*M + i], shiftAsub[j*M + i], 1e-16)) {
          fprintf(stderr, "[%d] : %d,%d : A : %lf -- Asub : %lf\n",
            rank, i, j, shiftA[j*M + i], shiftAsub[j*M + i]);
          is_ok = 0;
          break;
        }
      }
      if(!is_ok) {
        fprintf(stderr, "[%d] : Bad scatter for %dx%d matrix with %d block size\n",
        rank, M, N, block_size);
      }
    }
  }

  /* Gather */
  my_pdgetrf_gather(M, N, A, M, Asub, &band_type, &last_band_type);

  /* Check that 0 and 1 has the same matrix */
  double *A1;
  if(rank == 0) {
    A1 = alloc_matrix(M, N);
    MPI_Recv(A1, M*N, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int elt = 0; elt < M*N; ++elt) {
      if(!eq_double(A1[elt], A[elt], 1e-16)) {
        fprintf(stderr, "[0] : Bad gather for %dx%d matrix with %d block size\n",
          M, N, block_size);
        break;
      }  
    }
    
    free(A1);
  } else if (rank == 1) {
    MPI_Send(A, M * N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }

  free(Asub);

  /* Free types*/
  MPI_Type_free( &band_type );
  MPI_Type_free( &last_band_type );
  }
}

void check_pdgetrf() {
  int M = 9, N = 9;
  double *A = alloc_matrix(M, N);
  rnd_matrix_buff(A, 1, 10, M * N, 1);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  my_pdgetrf(CblasColMajor, M, N, A, M);

  free(A);
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  check_pdgetrf_split();
  MPI_Finalize();
  return EXIT_SUCCESS;
}