#include "algonum.h"
#include "codelets.h"

void
my_dgetrf_tiled_starpu( CBLAS_LAYOUT layout,
                        int M, int N, int b, double **A )
{
    starpu_data_handle_t *handlesA;
    starpu_data_handle_t hAkk, hAkn, hAmk, hAmn;

    /* Let's compute the total number of tiles with a *ceil* */
    int MT = my_iceil( M, b );
    int NT = my_iceil( N, b );
    int KT = my_imin( MT, NT );
    int m, n, k;

    handlesA = calloc( MT * NT, sizeof(starpu_data_handle_t) );

    /* TODO */

    unregister_starpu_handle( MT * NT, handlesA );

    /* Let's wait for the end of all the tasks */
    starpu_task_wait_for_all();
#if defined(ENABLE_MPI)
    starpu_mpi_barrier(MPI_COMM_WORLD);
#endif

    free( handlesA );
}
