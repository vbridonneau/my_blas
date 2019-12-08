#include <mpi.h>
#include "dgetrf.h"
#include "algonum.h"
#include "util.h"
#include <string.h>

const int block_size = 1;

void my_pdgetrf(const CBLAS_LAYOUT Order, int M, int N, double* A, int lda ) {
  MPI_Datatype band_type, last_band_type;

  /* Treat too small matrix size */
  if (block_size >= min( M, N )) {
    my_dgetf2( Order, M, N, A, lda );
  }

  /* Get ranks */
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

  /***********************************************/
  /* FIXME : Create a data type for band to send */
  /* by A and one to send by Asub.               */
  /* Or think about array of array.              */
  /***********************************************/

  /* Matrices */
  int NT_treated = (rank <= rlast) ? NT : NT - 1;
  int matrix_size;
  if (is_last && !is_last_complete) {
    matrix_size = (NT_treated - 1) * band_size + last_band_size;
  } else {
    matrix_size = NT_treated * band_size;
  }
  double *Asub = malloc(sizeof(double) * matrix_size);

  /* Datatype creation */
  MPI_Type_vector( block_size, M, lda, MPI_DOUBLE, &band_type );
  MPI_Type_vector( N % block_size, M, lda, MPI_DOUBLE, &last_band_type );

  /* Commit types */
  MPI_Type_commit( &band_type );
  MPI_Type_commit( &last_band_type );

  /* Split */

  /* Begin by scattering complete band to all processes */
  for (int i = 0; i < NT - 1; ++i) {
    MPI_Scatter(A + i * size * block_size * lda, size, band_type, Asub + i * block_size * M, 1, band_type, 0, MPI_COMM_WORLD);
  }

  /* Then scatter last incomplete bands to last processes */
  if (NT == NT_treated) {
    if (rank == 0) {
      int shiftA = (NT - 1) * size * block_size * lda;

      /* 0 is the last, so it receives the last band in 'Asub' */
      if (rlast == 0) {
        double *recvAsub  = Asub + (NT - 1) * band_size;
        double *ptrAshift = A + shiftA;
        int     last_size = (is_last_complete) ? band_size : last_band_size;
        memcpy(recvAsub, ptrAshift, size);
      } 
      /* Ohterwise 0 sends to all untill the last one a complete band
       * and it manages the last band case */
      else {
        /* rlast - 1 complete band */
        for(int rnk = 1; rnk < rlast; ++rnk) {
          MPI_Send(A + shiftA + rnk * block_size * lda, 1, band_type, rnk, 0, MPI_COMM_WORLD);
        }
        /* Last block info to send it last band properly */
        MPI_Datatype ltype = (is_last_complete) ? band_type : last_band_type;
        double *sendAptr   = A + shiftA + rlast * block_size * lda;
        int last_size      = (is_last_complete) ? band_size : last_band_size;
        MPI_Send(sendAptr, 1, ltype, rlast, 0, MPI_COMM_WORLD);
      }
    }
    /* Other processes receive their band */
    else if (rank <= rlast) {
      MPI_Datatype ltype = (!is_last || is_last_complete) ? band_type : last_band_type;
      double *recvAsub   = Asub + (NT - 1) * block_size * M;
      int last_size      = (!is_last || is_last_complete) ? band_size : last_band_size;
      MPI_Recv(recvAsub, last_size, ltype, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  }

  /* Computation */

  /* Gather data */
  for (int i = 0; i < NT - 1; ++i) {
    MPI_Gather(Asub + i * band_size, 1, band_type, A + i * size * band_size, size, band_type, 0, MPI_COMM_WORLD);
  }

  free(Asub);

  /* Free types*/
  MPI_Type_free( &band_type );
  MPI_Type_free( &last_band_type );
}
