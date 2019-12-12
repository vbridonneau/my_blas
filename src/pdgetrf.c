#include <mpi.h>
#include "dgetrf.h"
#include "algonum.h"
#include "util.h"
#include "dtrsm.h"
#include "dgemm.h"
#include <string.h>
#include <stdio.h>
#include <unistd.h>

void my_pdgetrf_scatter(int M, int N, int block_size, double* A, int lda, double *Asub, MPI_Datatype *pband_type, MPI_Datatype *plast_band_type) {
  /* MPI info */
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* Info for bandes per process */
  int Ndb              = (N + block_size - 1) / block_size;
  int NT               = (Ndb + size - 1) / size;
  int rlast            = ((Ndb % size) + size - 1) % size;
  int is_last          = (rank == rlast);
  int band_size        = block_size * lda;
  int NT_treated       = (rank <= rlast) ? NT : NT - 1;
  int last_band_size   = (N % block_size) * lda;
  int is_last_complete = ((N % block_size) == 0);

  /* Types */
  MPI_Datatype band_type      = *pband_type;
  MPI_Datatype last_band_type = *plast_band_type;

  /* Begin by scattering complete band to all processes */
  for (int i = 0; i < NT - 1; ++i) {
    MPI_Scatter(A + i * size * band_size, 1, band_type, Asub + i * band_size, 1, band_type, 0, MPI_COMM_WORLD);
  }

  /* Then scatter last incomplete bands to last processes */
  if (NT == NT_treated) {
    if (rank == 0) {
      int shiftA = (NT - 1) * size * band_size;

      /* 0 receives its the last band in 'Asub' */
      double *recvAsub  = Asub + (NT - 1) * band_size;
      double *ptrAshift = A + shiftA;
      int     last_size = (!is_last || is_last_complete) ? band_size : last_band_size;
      memcpy(recvAsub, ptrAshift, last_size * sizeof(double));
    
      /* 0 sends to all untill the last complete band
       * and it manages the last band case */
      /* rlast - 1 complete band */
      for(int rnk = 1; rnk < rlast; ++rnk) {
        MPI_Send(A + shiftA + rnk * band_size, 1, band_type, rnk, 0, MPI_COMM_WORLD);
      }
      if(rlast != 0) {
        /* Last band info to send it properly */
        MPI_Datatype ltype = (is_last_complete) ? band_type : last_band_type;
        double *sendAptr   = A + shiftA + rlast * block_size * lda;
        int last_size      = (is_last_complete) ? band_size : last_band_size;
        MPI_Send(sendAptr, 1, ltype, rlast, 0, MPI_COMM_WORLD);
      }
    }
    /* Other processes receive their band */
    else if (rank <= rlast) {
      MPI_Datatype ltype = (!is_last || is_last_complete) ? band_type : last_band_type;
      double *recvAsub   = Asub + (NT - 1) * band_size;
      int last_size      = (!is_last || is_last_complete) ? band_size : last_band_size;
      MPI_Recv(recvAsub, last_size, ltype, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  }
}

void my_pdgetrf_gather(int M, int N, int block_size, double* A, int lda, double *Asub, MPI_Datatype *pband_type, MPI_Datatype *plast_band_type) {
  /* MPI info */
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* Info for bandes per process */
  int Ndb              = (N + block_size - 1) / block_size;
  int NT               = (Ndb + size - 1) / size;
  int rlast            = ((Ndb % size) + size - 1) % size;
  int is_last          = (rank == rlast);
  int band_size        = block_size * lda;
  int NT_treated       = (rank <= rlast) ? NT : NT - 1;
  int last_band_size   = (N % block_size) * lda;
  int is_last_complete = ((N % block_size) == 0);

  /* Types */
  MPI_Datatype band_type      = *pband_type;
  MPI_Datatype last_band_type = *plast_band_type;

  for (int i = 0; i < NT - 1; ++i) {
    MPI_Gather(Asub + i * band_size, 1, band_type, A + i * size * band_size, 1, band_type, 0, MPI_COMM_WORLD);
  }

  /* Then gather by hand from processes that have one more band */
  if (NT == NT_treated) {
    /* 0 receives data from all other processes */
    if (rank == 0) {
      int shiftA = (NT - 1) * size * band_size;

      /* First 0 takes its data from Asub */
      double *ptrAshift = A + shiftA;
      double *sendAsub  = Asub + (NT - 1) * band_size;
      int     last_size = (!is_last || is_last_complete) ? band_size : last_band_size;
      memcpy(ptrAshift, sendAsub, last_size * sizeof(double));

      /* Receive from all others until rlast excluded */
      for (int rnk = 1; rnk < rlast; ++rnk) {
        MPI_Recv(A + shiftA + rnk * band_size, 1, band_type, rnk, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
      /* Treat last band case */
      if (rlast != 0) {
        MPI_Datatype ltype = (is_last_complete) ? band_type : last_band_type;
        double *recvAptr   = A + shiftA + rlast * band_size;
        int     last_size  = (is_last_complete) ? band_size : last_band_size;
        MPI_Recv(recvAptr, 1, ltype, rlast, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }      
    }
    /* Others send their last band to 0 */
    else if (rank <= rlast) {
      MPI_Datatype rtype     = (is_last && !is_last_complete) ? last_band_type : band_type;
      double      *sendAsub  = Asub + (NT - 1) * band_size;
      int          last_size = (is_last && !is_last_complete) ? last_band_size : band_size;
      MPI_Send(sendAsub, 1, rtype, 0, 0, MPI_COMM_WORLD);
    }
  }
}

void my_pdgetrf(const CBLAS_LAYOUT Order, int M, int N, int block_size, double* A, int lda ) {
  MPI_Datatype band_type, last_band_type;

  /* Get ranks */
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  /* Treat too small matrix size */
  if (block_size >= min( M, N )) {
    if(rank == 0)my_dgetf2( Order, M, N, A, lda );
    return;
  }

  /* Info for bandes per process */
  int Ndb              = (N + block_size - 1) / block_size;
  int NT               = (Ndb + size - 1) / size;
  int rlast            = ((Ndb % size) + size - 1) % size;
  int is_last          = (rank == rlast);
  int band_size        = block_size * lda;
  int last_band_size   = (N % block_size) * lda;
  int is_last_complete = ((N % block_size) == 0);

  /***********************************************/
  /* FIXME : Create a data type for band to send */
  /* by A and one to send by Asub.               */
  /* Or think about array of array.              */
  /* Warning : Not relevant in first             */
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

  /**************************** Split *********************/
  my_pdgetrf_scatter(M, N, block_size, A, lda, Asub, &band_type, &last_band_type);

  /********************** Computation ************************/
  double *placeholder, *tmp_band = malloc(sizeof(double) * band_size);
  int Nmax    = (is_last && !is_last_complete) ? 
    (NT_treated - 1) * block_size + (N % block_size)
    : NT_treated * block_size;
  int maxiter = min( M, N );
  int jb;
  int jbb;
  for(int j = 0, bandid = 0; j < maxiter; j += block_size, ++ bandid) {
    jb  = min( min(M, N) - j, block_size );
    jbb = (bandid == Ndb-1 && !is_last_complete) ? (N % block_size) : block_size;

    MPI_Datatype btype = (bandid == Ndb-1 && !is_last_complete) ?
      last_band_type : band_type;
    /* Root perform dgetrf */
    if( rank == (bandid % size) ) {
      int local_col_shift = (bandid / size) * band_size;
      placeholder = Asub + local_col_shift;
      my_dgetf2( Order, M - j, jbb, placeholder + j, lda );
    } else {
      placeholder = tmp_band;
    }
    if (bandid != Ndb-1) {
      MPI_Bcast( placeholder, 1, btype, (bandid % size), MPI_COMM_WORLD );
      /* Perform a lot of dtrsm and dgemm instead of 1 */
      for(int b = (bandid / size) + (bandid % size >= rank); b < NT_treated; ++b) {
        /* Pointer for B in dtrsm arg's name */
        double *Aptr = Asub + b * band_size;
        int jj  = b * block_size;
        int sub_blocksize = (b == NT_treated - 1 && is_last && !is_last_complete) ? (N % block_size) : block_size;//min ( min(M, N) - jj, block_size );
        if (jj + sub_blocksize <= Nmax) {
          my_dtrsm( Order, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, jb, sub_blocksize, 1., placeholder + j, lda, Aptr + j, lda );
          if (j + jb < M) {
            my_dgemm( Order, CblasNoTrans, CblasNoTrans, M - j - jb, sub_blocksize, jbb, -1., placeholder + j + jb, lda, Aptr + j, lda, 1., Aptr + j + jb, lda );
          }
        }
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  /********************** Gather data ************************/
  my_pdgetrf_gather(M, N, block_size, A, lda, Asub, &band_type, &last_band_type);

  /*********************** Free Section ***********************/
  free(tmp_band);
  free(Asub);

  /* Free types*/
  MPI_Type_free( &band_type );
  MPI_Type_free( &last_band_type );
}
