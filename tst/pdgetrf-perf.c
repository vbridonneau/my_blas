#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include "dgetrf.h"
#include "util.h"


void check_pdgetrf_split() {
  for (int Mb = 1; Mb < 10; ++Mb) {
    for (int M = 9; M < 100; M += 5) {
      for(int N = 9; N < 100; N += 5) {
	double *A = alloc_matrix(M, N);
	rnd_matrix_buff(A, 1, 10, M * N, 1);

	/* Datatypes */
	MPI_Datatype band_type, last_band_type;

	/* Datatype creation */
	MPI_Type_vector( Mb, M, M, MPI_DOUBLE, &band_type );
	MPI_Type_vector( N % Mb, M, M, MPI_DOUBLE, &last_band_type );

	/* Commit types */
	MPI_Type_commit( &band_type );
	MPI_Type_commit( &last_band_type );

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	/* Info for bandes per process */
	int Ndb              = (N + Mb - 1) / Mb;
	int NT               = (Ndb + size - 1) / size;
	int rlast            = ((Ndb % size) + size - 1) % size;
	int is_last          = (rank == rlast);
	int band_size        = Mb * M;
	int last_band_size   = (N % Mb) * M;
	int is_last_complete = ((N % Mb) == 0);

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
	my_pdgetrf_scatter(M, N, Mb, A, M, Asub, &band_type, &last_band_type);

	/* Check everyone has the correct bands */
	double *my_A = A + rank * band_size;
	for (int band = 0; band < NT_treated; ++band) {
	  int is_ok = 1;
	  double *shiftAsub = Asub + band * band_size;
	  double *shiftA    = my_A + band * band_size * size;
	  int maxj = (band == NT_treated - 1 && is_last && !is_last_complete) ? (N % Mb) : Mb;
	  for (int j = 0; j < maxj; ++j) {
	    for (int i = 0; i < M; ++i) {
	      if (!eq_double(shiftA[j*M + i], shiftAsub[j*M + i], 1e-16)) {
		fprintf(stderr, "[%d] : %d,%d,%d : A : %lf -- Asub : %lf\n",
			rank, band, i, j, shiftA[j*M + i], shiftAsub[j*M + i]);
		is_ok = 0;
		break;
	      }
	    }
	    if(!is_ok) {
	      fprintf(stderr, "[%d] : Bad scatter for %dx%d matrix with %d block size\n",
		      rank, M, N, Mb);
	      break;
	    }
	  }
	  if(!is_ok)break;
	}

	/* Gather */
	my_pdgetrf_gather(M, N, Mb, A, M, Asub, &band_type, &last_band_type);

	/* Check that 0 and 1 has the same matrix */
	double *A1;
	if(rank == 0) {
	  A1 = alloc_matrix(M, N);
	  MPI_Recv(A1, M*N, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	  for (int elt = 0; elt < M*N; ++elt) {
	    if(!eq_double(A1[elt], A[elt], 1e-16)) {
	      fprintf(stderr, "[0] : Bad gather for %dx%d matrix with %d block size\n",
		      M, N, Mb);
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
  }
}

void check_pdgetrf() {
	for (int Mb = 1; Mb < 10; ++Mb)
	{
		for (int M = 10; M < 100; M += 10)
		{
			for (int N = 10; N < 100; N += 10)
			{
				double *A = alloc_matrix(M, N);
				double *my_LU = alloc_matrix(M, N);

				int rank;
				MPI_Comm_rank(MPI_COMM_WORLD, &rank);

				if (rank == 0) {
					rnd_matrix_buff(A, 1, 10, M*N, 1);

					int maxiter = min(M, N);
					for (int i = 0; i < maxiter; ++i)
					{
						A[i + i *maxiter] += 50;
					}

					memcpy(my_LU , A, sizeof(double) * M * N);
				}


				MPI_Barrier(MPI_COMM_WORLD);
				my_pdgetrf(CblasColMajor, M, N, Mb, my_LU, M);
				MPI_Barrier(MPI_COMM_WORLD);

				if(rank == 0) {
					my_dgetrf(CblasColMajor, M, N, A, M);
					int is_ok = 1;
					for (int j = 0; j < N; ++j) {
						for (int i = 0; i < M; ++i) {
							if (!eq_double(my_LU[j*M + i], A[j*M + i], 1e-16)) {
								fprintf(stderr, "%d,%d : A : %lf -- LU : %lf\n",
									i, j, A[j*M + i], my_LU[j*M + i]);
								is_ok = 0;
								break;
							}
						}
						if(!is_ok) {
							fprintf(stderr, "Bad result for %dx%d matrix with %d block size\n",
								M, N, Mb);
							break;
						}
					}
				// printf("myLU mpi\n");
				// affiche(M, N, my_LU, M, stdout);
				// printf("myLU seq\n");
				// affiche(M, N, A, M, stdout);
				}
				free(my_LU);
				free(A);
				MPI_Barrier(MPI_COMM_WORLD);
			}
		}	
	}
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  check_pdgetrf();
  MPI_Finalize();
  return EXIT_SUCCESS;
}
