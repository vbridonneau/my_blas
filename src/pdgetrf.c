#include <mpi.h>
#include "dgetrf.h"
#include "algonum.h"

const int block_size = 1;

void my_pdgetrf(const CBLAS_LAYOUT Order, int M, int N, double* A, int lda ) {
  MPI_Datatype band_type, last_band_type;

  /* Treat too small matrix size */
  if (block_size >= min( M, N )) {
    my_dgetf2( Order, M, N, A, lda );
  }

  /* Datatype creation */
  MPI_Type_vector( block_size, M, lda, MPI_DOUBLE, &band_type );
  MPI_Type_vector( N % block_size, M, lda, MPI_DOUBLE, &last_band_type );

  /* Commit types */
  MPI_Type_commit( &band_type );
  MPI_Type_commit( &last_band_type );

  /* Split */

  /* Computation */

  /* Gather data */

  /* Free types*/
  MPI_Type_free( &band_type );
  MPI_Type_free( &last_band_type );
}
