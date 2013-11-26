#include <bfam_subdomain_dgx.h>
#include <bfam_jacobi.h>
#include <bfam_kron.h>
#include <bfam_log.h>
#include <bfam_util.h>
#include <bfam_vtk.h>

#ifndef BFAM_DGX_DIMENSION
#define BFAM_DGX_DIMENSION
#define USE_GENERIC_DGX_DIMENSION
#else
#define DIM (BFAM_DGX_DIMENSION)
#endif

static void
bfam_subdomain_dgx_geo(int N, bfam_locidx_t K, int Np, int ***gmask,
                        const int *restrict Ng, const int *restrict Ngp,
                        int num_Vi,
                    /*const*/ bfam_long_real_t **restrict xi,
                    /*const*/ bfam_long_real_t  *restrict Dr,
                              bfam_long_real_t **restrict Jrx,
                              bfam_long_real_t  *restrict J,
                              bfam_long_real_t **restrict ni,
                              bfam_long_real_t *restrict sJ,
                              const int inDIM)
{
#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_geo");
  const int DIM = inDIM;
#endif
  BFAM_ABORT_IF(num_Vi != DIM && !(J == NULL && ni == NULL && sJ == NULL),
      "[J,ni,sJ] != NULL not implemented for this case");
  BFAM_ABORT_IF_NOT((J == NULL) == (ni == NULL) && (J == NULL) == (sJ == NULL),
      "[J,ni,sJ] must all be NULL or not NULL");

  const int Nfaces = Ng[0];
  const int Nfp    = Ngp[0];

  /*
  BFAM_ASSUME_ALIGNED( x, 32);
  BFAM_ASSUME_ALIGNED( y, 32);
  BFAM_ASSUME_ALIGNED(Dr, 32);
  BFAM_ASSUME_ALIGNED(Jrx, 32);
  BFAM_ASSUME_ALIGNED(Jsx, 32);
  BFAM_ASSUME_ALIGNED(Jry, 32);
  BFAM_ASSUME_ALIGNED(Jsy, 32);
  BFAM_ASSUME_ALIGNED( J, 32);
  */

  BFAM_ASSERT(DIM > 0);
  for(bfam_locidx_t k = 0, vsk = 0, fsk = 0; k < K; ++k)
  {
    if(DIM == 1)
      for(int v = 0; v < num_Vi; v++)
      {
        for(int i = 0; i < N+1; i++) Jrx[v][vsk+i] = 0;
        bfam_util_mvmult(N+1,N+1,Dr,N+1,xi[v]+vsk,Jrx[v]+vsk);
      }
    else if(DIM == 2)
      for(int v = 0; v < num_Vi; v++)
      {
        BFAM_KRON_IXA(N+1, Dr, xi[v] + vsk, Jrx[0+2*v] + vsk); /* xr */
        BFAM_KRON_AXI(N+1, Dr, xi[v] + vsk, Jrx[1+2*v] + vsk); /* xs */
      }
    else if(DIM == 3)
      for(int v = 0; v < num_Vi; v++)
      {
        BFAM_KRON_IXIXA(N+1, Dr, xi[v] + vsk, Jrx[0+3*v] + vsk); /* xi_r */
        BFAM_KRON_IXAXI(N+1, Dr, xi[v] + vsk, Jrx[1+3*v] + vsk); /* xi_s */
        BFAM_KRON_AXIXI(N+1, Dr, xi[v] + vsk, Jrx[2+3*v] + vsk); /* xi_t */
      }
    else BFAM_ABORT("Cannot handle dim = %d",DIM);

    if(J)
    {
      if(DIM == 1)
      {
        for(int n = 0; n < Np; ++n)
          J[n+vsk] = BFAM_LONG_REAL_ABS(Jrx[0][n+vsk]);
        for(int n = 0; n < Nfp; ++n)
        {
          const bfam_locidx_t fidx0 = fsk + 0 * Nfp + n;
          const bfam_locidx_t fidx1 = fsk + 1 * Nfp + n;

          const bfam_locidx_t vidx0 = vsk + gmask[0][0][n];
          const bfam_locidx_t vidx1 = vsk + gmask[0][1][n];

          /* face 0 */
          ni[0][fidx0] = -Jrx[0][vidx0]; /* -sy */

          /* face 1 */
          ni[0][fidx1] =  Jrx[0][vidx1]; /*  sy */


          sJ[fidx0] = BFAM_LONG_REAL_ABS(ni[0][fidx0]);
          sJ[fidx1] = BFAM_LONG_REAL_ABS(ni[0][fidx1]);

          ni[0][fidx0] /= sJ[fidx0];
          ni[0][fidx1] /= sJ[fidx1];
        }
      }
      else if(DIM == 2)
      {
        for(int n = 0; n < Np; ++n)
        {
          bfam_locidx_t idx = n + vsk;

          const bfam_long_real_t xr = Jrx[0][idx];
          const bfam_long_real_t xs = Jrx[1][idx];
          const bfam_long_real_t yr = Jrx[2][idx];
          const bfam_long_real_t ys = Jrx[3][idx];

          /* xr*ys - xs*yr */
          J[idx] = xr*ys - xs*yr;

          /* J*rx = ys */
          Jrx[0][idx] =  ys;

          /* J*ry = -xs */
          Jrx[1][idx] = -xs;

          /* J*sx = -yr */
          Jrx[2][idx] = -yr;

          /* J*sy = xr */
          Jrx[3][idx] =  xr;
        }
        for(int n = 0; n < Nfp; ++n)
        {
          const bfam_locidx_t fidx0 = fsk + 0 * Nfp + n;
          const bfam_locidx_t fidx1 = fsk + 1 * Nfp + n;
          const bfam_locidx_t fidx2 = fsk + 2 * Nfp + n;
          const bfam_locidx_t fidx3 = fsk + 3 * Nfp + n;

          const bfam_locidx_t vidx0 = vsk + gmask[0][0][n];
          const bfam_locidx_t vidx1 = vsk + gmask[0][1][n];
          const bfam_locidx_t vidx2 = vsk + gmask[0][2][n];
          const bfam_locidx_t vidx3 = vsk + gmask[0][3][n];

          /* rx = 0; ry = 1; sx = 2; sy = 3 */

          /* face 0 */
          ni[0][fidx0] = -Jrx[0][vidx0]; /* -Jrx/sJ */
          ni[1][fidx0] = -Jrx[1][vidx0]; /* -Jry/sJ */

          /* face 1 */
          ni[0][fidx1] =  Jrx[0][vidx1]; /*  Jrx/sJ */
          ni[1][fidx1] =  Jrx[1][vidx1]; /*  Jry/sJ */

          /* face 2 */
          ni[0][fidx2] = -Jrx[2][vidx2]; /* -Jsx/sJ */
          ni[1][fidx2] = -Jrx[3][vidx2]; /* -Jsy/sJ */

          /* face 3 */
          ni[0][fidx3] =  Jrx[2][vidx3]; /*  Jsx/sJ */
          ni[1][fidx3] =  Jrx[3][vidx3]; /*  Jsy/sJ */

          sJ[fidx0] = BFAM_LONG_REAL_HYPOT(ni[0][fidx0],ni[1][fidx0]);
          sJ[fidx1] = BFAM_LONG_REAL_HYPOT(ni[0][fidx1],ni[1][fidx1]);
          sJ[fidx2] = BFAM_LONG_REAL_HYPOT(ni[0][fidx2],ni[1][fidx2]);
          sJ[fidx3] = BFAM_LONG_REAL_HYPOT(ni[0][fidx3],ni[1][fidx3]);

          ni[0][fidx0] /= sJ[fidx0];
          ni[1][fidx0] /= sJ[fidx0];
          ni[0][fidx1] /= sJ[fidx1];
          ni[1][fidx1] /= sJ[fidx1];
          ni[0][fidx2] /= sJ[fidx2];
          ni[1][fidx2] /= sJ[fidx2];
          ni[0][fidx3] /= sJ[fidx3];
          ni[1][fidx3] /= sJ[fidx3];
        }
      }
      else if(DIM == 3)
      {
        for(int n = 0; n < Np; ++n)
        {
          bfam_locidx_t idx = n + vsk;

          const bfam_long_real_t xr = Jrx[0][idx];
          const bfam_long_real_t xs = Jrx[1][idx];
          const bfam_long_real_t xt = Jrx[2][idx];
          const bfam_long_real_t yr = Jrx[3][idx];
          const bfam_long_real_t ys = Jrx[4][idx];
          const bfam_long_real_t yt = Jrx[5][idx];
          const bfam_long_real_t zr = Jrx[6][idx];
          const bfam_long_real_t zs = Jrx[7][idx];
          const bfam_long_real_t zt = Jrx[8][idx];

          /*       xr*(ys*zt-yt*zs) - xs*(yr*zt-yt*zr) + xt*(yr*zs-ys*zr) */
          J[idx] = xr*(ys*zt-yt*zs) - xs*(yr*zt-yt*zr) + xt*(yr*zs-ys*zr);

          /* J*rx     =  ( ys * zt - yt * zs ) */
          Jrx[0][idx] =  ( ys * zt - yt * zs );

          /* J*ry     = -( xs * zt - xt * zs ) */
          Jrx[1][idx] = -( xs * zt - xt * zs );

          /* J*rz     =  ( xs * yt - xt * ys ) */
          Jrx[2][idx] =  ( xs * yt - xt * ys );

          /* J*sx     = -( yr * zt - yt * zr ) */
          Jrx[3][idx] = -( yr * zt - yt * zr );

          /* J*sy     =  ( xr * zt - xt * zr ) */
          Jrx[4][idx] =  ( xr * zt - xt * zr );

          /* J*sz     = -( xr * yt - xt * yr ) */
          Jrx[5][idx] = -( xr * yt - xt * yr );

          /* J*tx     =  ( yr * zs - ys * zr ) */
          Jrx[6][idx] =  ( yr * zs - ys * zr );

          /* J*ty     = -( xr * zs - xs * zr ) */
          Jrx[7][idx] = -( xr * zs - xs * zr );

          /* J*tz     =  ( xr * ys - xs * yr ) */
          Jrx[8][idx] =  ( xr * ys - xs * yr );
        }
        for(int n = 0; n < Nfp; ++n)
        {
          const bfam_locidx_t fidx0 = fsk + 0 * Nfp + n;
          const bfam_locidx_t fidx1 = fsk + 1 * Nfp + n;
          const bfam_locidx_t fidx2 = fsk + 2 * Nfp + n;
          const bfam_locidx_t fidx3 = fsk + 3 * Nfp + n;
          const bfam_locidx_t fidx4 = fsk + 4 * Nfp + n;
          const bfam_locidx_t fidx5 = fsk + 5 * Nfp + n;

          const bfam_locidx_t vidx0 = vsk + gmask[0][0][n];
          const bfam_locidx_t vidx1 = vsk + gmask[0][1][n];
          const bfam_locidx_t vidx2 = vsk + gmask[0][2][n];
          const bfam_locidx_t vidx3 = vsk + gmask[0][3][n];
          const bfam_locidx_t vidx4 = vsk + gmask[0][4][n];
          const bfam_locidx_t vidx5 = vsk + gmask[0][5][n];

          /* face 0 */
          ni[0][fidx0] = -Jrx[0][vidx0]; /* -Jrx/sJ */
          ni[1][fidx0] = -Jrx[1][vidx0]; /* -Jry/sJ */
          ni[2][fidx0] = -Jrx[2][vidx0]; /* -Jrz/sJ */

          /* face 1 */
          ni[0][fidx1] =  Jrx[0][vidx1]; /*  Jrx/sJ */
          ni[1][fidx1] =  Jrx[1][vidx1]; /*  Jry/sJ */
          ni[2][fidx1] =  Jrx[2][vidx1]; /*  Jrz/sJ */

          /* face 2 */
          ni[0][fidx2] = -Jrx[3][vidx2]; /* -Jsx/sJ */
          ni[1][fidx2] = -Jrx[4][vidx2]; /* -Jsy/sJ */
          ni[2][fidx2] = -Jrx[5][vidx2]; /* -Jsz/sJ */

          /* face 3 */
          ni[0][fidx3] =  Jrx[3][vidx3]; /*  Jsx/sJ */
          ni[1][fidx3] =  Jrx[4][vidx3]; /*  Jsy/sJ */
          ni[2][fidx3] =  Jrx[5][vidx3]; /*  Jsz/sJ */

          /* face 4 */
          ni[0][fidx4] = -Jrx[6][vidx4]; /* -Jtx/sJ */
          ni[1][fidx4] = -Jrx[7][vidx4]; /* -Jty/sJ */
          ni[2][fidx4] = -Jrx[8][vidx4]; /* -Jtz/sJ */

          /* face 5 */
          ni[0][fidx5] =  Jrx[6][vidx5]; /*  Jtx/sJ */
          ni[1][fidx5] =  Jrx[7][vidx5]; /*  Jty/sJ */
          ni[2][fidx5] =  Jrx[8][vidx5]; /*  Jtz/sJ */

          sJ[fidx0] = BFAM_LONG_REAL_HYPOT(ni[0][fidx0],
                        BFAM_LONG_REAL_HYPOT(ni[1][fidx0],ni[2][fidx0]));
          sJ[fidx1] = BFAM_LONG_REAL_HYPOT(ni[0][fidx1],
                        BFAM_LONG_REAL_HYPOT(ni[1][fidx1],ni[2][fidx1]));
          sJ[fidx2] = BFAM_LONG_REAL_HYPOT(ni[0][fidx2],
                        BFAM_LONG_REAL_HYPOT(ni[1][fidx2],ni[2][fidx2]));
          sJ[fidx3] = BFAM_LONG_REAL_HYPOT(ni[0][fidx3],
                        BFAM_LONG_REAL_HYPOT(ni[1][fidx3],ni[2][fidx3]));
          sJ[fidx4] = BFAM_LONG_REAL_HYPOT(ni[0][fidx4],
                        BFAM_LONG_REAL_HYPOT(ni[1][fidx4],ni[2][fidx4]));
          sJ[fidx5] = BFAM_LONG_REAL_HYPOT(ni[0][fidx5],
                        BFAM_LONG_REAL_HYPOT(ni[1][fidx5],ni[2][fidx5]));

          ni[0][fidx0] /= sJ[fidx0];
          ni[1][fidx0] /= sJ[fidx0];
          ni[2][fidx0] /= sJ[fidx0];
          ni[0][fidx1] /= sJ[fidx1];
          ni[1][fidx1] /= sJ[fidx1];
          ni[2][fidx1] /= sJ[fidx1];
          ni[0][fidx2] /= sJ[fidx2];
          ni[1][fidx2] /= sJ[fidx2];
          ni[2][fidx2] /= sJ[fidx2];
          ni[0][fidx3] /= sJ[fidx3];
          ni[1][fidx3] /= sJ[fidx3];
          ni[2][fidx3] /= sJ[fidx3];
          ni[0][fidx4] /= sJ[fidx4];
          ni[1][fidx4] /= sJ[fidx4];
          ni[2][fidx4] /= sJ[fidx4];
          ni[0][fidx5] /= sJ[fidx5];
          ni[1][fidx5] /= sJ[fidx5];
          ni[2][fidx5] /= sJ[fidx5];
        }
      }
    }

    vsk += Np;
    fsk += Nfaces*Nfp;
  }
}

static void
bfam_subdomain_dgx_field_init(bfam_subdomain_t *subdomain,
    const char *name, bfam_real_t time, bfam_subdomain_init_field_t init_field,
    void *arg)
{
  bfam_subdomain_dgx_t *s = (bfam_subdomain_dgx_t*) subdomain;

  bfam_real_t *field = bfam_dictionary_get_value_ptr(&s->base.fields,name);

  BFAM_ABORT_IF(field==NULL, "Init: Field %s not found in subdomain %s",
      name, subdomain->name);

  size_t fieldLength = s->Np*s->K;

  bfam_real_t *restrict x0 =
    bfam_dictionary_get_value_ptr(&subdomain->fields, "_grid_x0");
  bfam_real_t *restrict x1 =
    bfam_dictionary_get_value_ptr(&subdomain->fields, "_grid_x1");
  bfam_real_t *restrict x2 =
    bfam_dictionary_get_value_ptr(&subdomain->fields, "_grid_x2");
  /* this is just to check is we need to update the init function */
  BFAM_ASSERT(NULL == bfam_dictionary_get_value_ptr(&subdomain->fields,
                                                      "_grid_x3"));

  init_field(fieldLength, name, time, x0,x1,x2, subdomain, arg, field);
}

static void
bfam_subdomain_dgx_vtk_interp(bfam_locidx_t K,
    int N_d,       bfam_real_t * restrict d,
    int N_s, const bfam_real_t * restrict s,
    const bfam_real_t *restrict interp)
{
  BFAM_ABORT("interp not implemented");
}

static int
bfam_subdomain_dgx_vtk_write_vtu_piece(bfam_subdomain_t *subdomain,
    FILE *file, bfam_real_t time, const char **scalars, const char **vectors,
    const char **components, int writeBinary, int writeCompressed,
    int rank, bfam_locidx_t id, int Np_write)
{
  bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t*) subdomain;

  const char *format;

  if(writeBinary)
    format = "binary";
  else
    format = "ascii";

  const bfam_locidx_t  K  = sub->K;
  int              N_vtk  = sub->N;
  int              Np_vtk = sub->Np;
  bfam_real_t *interp = NULL;

  bfam_real_t *restrict stor1 = NULL;
  bfam_real_t *restrict stor2 = NULL;
  bfam_real_t *restrict stor3 = NULL;

  if(Np_write > 0)
  {
    BFAM_ABORT_IF_NOT(Np_write > 1, "Np_write = %d is not valid",
        Np_write);

    N_vtk  = Np_write - 1;
    Np_vtk = Np_write;
    if(sub->dim > 1) Np_vtk *= Np_write;
    if(sub->dim > 2) Np_vtk *= Np_write;

    interp = bfam_malloc_aligned(sizeof(bfam_real_t)*(sub->N+1)*(N_vtk+1));

    bfam_long_real_t *cal_interp =
      bfam_malloc_aligned(sizeof(bfam_long_real_t)*(sub->N+1)*(N_vtk+1));
    bfam_long_real_t *lr =
      bfam_malloc_aligned(sizeof(bfam_long_real_t)*Np_write);

    for(int r = 0; r < Np_write; r++)
      lr[r] = -1 + 2*(bfam_long_real_t)r/(Np_write-1);

    bfam_jacobi_p_interpolation(0, 0, sub->N, Np_write, lr, sub->V,cal_interp);

    for(int n = 0; n < (sub->N+1)*(N_vtk+1); n++)
      interp[n] = (bfam_real_t)cal_interp[n];

    stor1 = bfam_malloc_aligned(sizeof(bfam_real_t)*Np_vtk*K);
    stor2 = bfam_malloc_aligned(sizeof(bfam_real_t)*Np_vtk*K);
    stor3 = bfam_malloc_aligned(sizeof(bfam_real_t)*Np_vtk*K);

    bfam_free_aligned(lr);
    bfam_free_aligned(cal_interp);
  }

  const int Ncorners = sub->Ng[sub->numg-1];

  bfam_locidx_t Ncells = K * N_vtk;
  if(sub->dim > 1) Ncells *= N_vtk;
  if(sub->dim > 2) Ncells *= N_vtk;
  const bfam_locidx_t Ntotal = K * Np_vtk;

  bfam_real_t *restrict x =
    bfam_dictionary_get_value_ptr(&subdomain->fields, "_grid_x0");
  bfam_real_t *restrict y =
    bfam_dictionary_get_value_ptr(&subdomain->fields, "_grid_x1");
  bfam_real_t *restrict z =
    bfam_dictionary_get_value_ptr(&subdomain->fields, "_grid_x2");

  if(interp == NULL)
  {
    stor1 = x;
    stor2 = y;
    stor3 = z;
  }
  else
  {
    bfam_subdomain_dgx_vtk_interp(K,N_vtk,stor1,sub->N,x,interp);
    bfam_subdomain_dgx_vtk_interp(K,N_vtk,stor2,sub->N,y,interp);
    bfam_subdomain_dgx_vtk_interp(K,N_vtk,stor3,sub->N,z,interp);
  }

  fprintf(file,
           "    <Piece NumberOfPoints=\"%jd\" NumberOfCells=\"%jd\">\n",
           (intmax_t) Ntotal, (intmax_t) Ncells);

  /*
   * Points
   */
  fprintf (file, "      <Points>\n");

  bfam_vtk_write_real_vector_data_array(file, "Position", writeBinary,
      writeCompressed, Ntotal, stor1, stor2, stor3);

  fprintf(file, "      </Points>\n");

  /*
   * Cells
   */
  fprintf(file, "      <Cells>\n");

  /*
   * Connectivity
   */
  fprintf(file, "        <DataArray type=\"%s\" Name=\"connectivity\""
          " format=\"%s\">\n", BFAM_LOCIDX_VTK, format);
  if(writeBinary)
  {
    size_t cellsSize = Ncells*Ncorners*sizeof(bfam_locidx_t);
    bfam_locidx_t *cells = bfam_malloc_aligned(cellsSize);

    if(sub->dim == 1)
      for(bfam_locidx_t k = 0, i = 0; k < K; ++k)
      {
        for(int n = 0; n < N_vtk; ++n)
        {
          cells[i++] = Np_vtk * k + (n + 0);
          cells[i++] = Np_vtk * k + (n + 1);
        }
      }
    else if(sub->dim == 2)
      for(bfam_locidx_t k = 0, i = 0; k < K; ++k)
      {
        for(int m = 0; m < N_vtk; ++m)
        {
          for(int n = 0; n < N_vtk; ++n)
          {
            cells[i++] = Np_vtk * k + (N_vtk + 1) * (m + 0) + (n + 0);
            cells[i++] = Np_vtk * k + (N_vtk + 1) * (m + 0) + (n + 1);
            cells[i++] = Np_vtk * k + (N_vtk + 1) * (m + 1) + (n + 0);
            cells[i++] = Np_vtk * k + (N_vtk + 1) * (m + 1) + (n + 1);
          }
        }
      }
    else if(sub->dim == 3)
      for(bfam_locidx_t k = 0, i = 0; k < K; ++k)
      {
        for(int l = 0; l < N_vtk; ++l)
        {
          for(int m = 0; m < N_vtk; ++m)
          {
            for(int n = 0; n < N_vtk; ++n)
            {
              cells[i++] = Np_vtk*k +
                (N_vtk+1)*(N_vtk+1)*(l+0) + (N_vtk+1)*(m+0) + (n+0);
              cells[i++] = Np_vtk*k +
                (N_vtk+1)*(N_vtk+1)*(l+0) + (N_vtk+1)*(m+0) + (n+1);
              cells[i++] = Np_vtk*k +
                (N_vtk+1)*(N_vtk+1)*(l+0) + (N_vtk+1)*(m+1) + (n+0);
              cells[i++] = Np_vtk*k +
                (N_vtk+1)*(N_vtk+1)*(l+0) + (N_vtk+1)*(m+1) + (n+1);
              cells[i++] = Np_vtk*k +
                (N_vtk+1)*(N_vtk+1)*(l+1) + (N_vtk+1)*(m+0) + (n+0);
              cells[i++] = Np_vtk*k +
                (N_vtk+1)*(N_vtk+1)*(l+1) + (N_vtk+1)*(m+0) + (n+1);
              cells[i++] = Np_vtk*k +
                (N_vtk+1)*(N_vtk+1)*(l+1) + (N_vtk+1)*(m+1) + (n+0);
              cells[i++] = Np_vtk*k +
                (N_vtk+1)*(N_vtk+1)*(l+1) + (N_vtk+1)*(m+1) + (n+1);
            }
          }
        }
      }
    else BFAM_ABORT("not implemented for dim = %d", sub->dim);

    fprintf(file, "          ");
    int rval =
      bfam_vtk_write_binary_data(writeCompressed, file, (char*)cells,
          cellsSize);
    fprintf(file, "\n");
    if(rval)
      BFAM_WARNING("Error encoding cells");

    bfam_free_aligned(cells);
  }
  else
  {
    if(sub->dim == 1)
      for(bfam_locidx_t k = 0; k < K; ++k)
        for(int n = 0; n < N_vtk; ++n)
          fprintf(file,
              "          %8jd %8jd\n",
              (intmax_t) Np_vtk * k + (n + 0),
              (intmax_t) Np_vtk * k + (n + 1));
    else if(sub->dim == 2)
      for(bfam_locidx_t k = 0; k < K; ++k)
        for(int m = 0; m < N_vtk; ++m)
          for(int n = 0; n < N_vtk; ++n)
            fprintf(file,
                "          %8jd %8jd %8jd %8jd\n",
                (intmax_t) Np_vtk * k + (N_vtk + 1) * (m + 0) + (n + 0),
                (intmax_t) Np_vtk * k + (N_vtk + 1) * (m + 0) + (n + 1),
                (intmax_t) Np_vtk * k + (N_vtk + 1) * (m + 1) + (n + 0),
                (intmax_t) Np_vtk * k + (N_vtk + 1) * (m + 1) + (n + 1));
    else if(sub->dim == 3)
      for(bfam_locidx_t k = 0; k < K; ++k)
        for(int l = 0; l < N_vtk; ++l)
          for(int m = 0; m < N_vtk; ++m)
            for(int n = 0; n < N_vtk; ++n)
              fprintf(file,
                  "          %8jd %8jd %8jd %8jd %8jd %8jd %8jd %8jd\n",
                  (intmax_t) Np_vtk*k
                      +(N_vtk+1)*(N_vtk+1)*(l+0) + (N_vtk+1)*(m+0) + (n+0),
                  (intmax_t) Np_vtk*k
                      +(N_vtk+1)*(N_vtk+1)*(l+0) + (N_vtk+1)*(m+0) + (n+1),
                  (intmax_t) Np_vtk*k
                      +(N_vtk+1)*(N_vtk+1)*(l+0) + (N_vtk+1)*(m+1) + (n+0),
                  (intmax_t) Np_vtk*k
                      +(N_vtk+1)*(N_vtk+1)*(l+0) + (N_vtk+1)*(m+1) + (n+1),
                  (intmax_t) Np_vtk*k
                      +(N_vtk+1)*(N_vtk+1)*(l+1) + (N_vtk+1)*(m+0) + (n+0),
                  (intmax_t) Np_vtk*k
                      +(N_vtk+1)*(N_vtk+1)*(l+1) + (N_vtk+1)*(m+0) + (n+1),
                  (intmax_t) Np_vtk*k
                      +(N_vtk+1)*(N_vtk+1)*(l+1) + (N_vtk+1)*(m+1) + (n+0),
                  (intmax_t) Np_vtk*k
                      +(N_vtk+1)*(N_vtk+1)*(l+1) + (N_vtk+1)*(m+1) + (n+1));
    else BFAM_ABORT("not implemented for dim = %d %d", sub->dim);
  }
  fprintf(file, "        </DataArray>\n");

  /*
   * Offsets
   */
  fprintf (file, "        <DataArray type=\"%s\" Name=\"offsets\""
           " format=\"%s\">\n", BFAM_LOCIDX_VTK, format);
  fprintf(file, "          ");
  if(writeBinary)
  {
    size_t offsetsSize = Ncells*sizeof(bfam_locidx_t);
    bfam_locidx_t *offsets = bfam_malloc_aligned(offsetsSize);

    for(bfam_locidx_t i = 1; i <= Ncells; ++i)
      offsets[i - 1] = Ncorners * i;

    int rval =
      bfam_vtk_write_binary_data(writeCompressed, file, (char*)offsets,
          offsetsSize);
    if(rval)
      BFAM_WARNING("Error encoding offsets");

    bfam_free_aligned(offsets);
  }
  else
  {
    for(bfam_locidx_t i = 1, sk = 1; i <= Ncells; ++i, ++sk)
    {
      fprintf(file, " %8jd", (intmax_t) (Ncorners * i));
      if(!(sk % 20) && i != Ncells)
        fprintf(file, "\n          ");
    }
  }
  fprintf(file, "\n");
  fprintf(file, "        </DataArray>\n");

  /*
   * Types
   */
  fprintf(file, "        <DataArray type=\"UInt8\" Name=\"types\""
           " format=\"%s\">\n", format);
  fprintf(file, "          ");
  if(writeBinary)
  {
    size_t typesSize = Ncells*sizeof(uint8_t);
    uint8_t *types = bfam_malloc_aligned(typesSize);

    if(sub->dim == 1)
      for(bfam_locidx_t i = 0; i < Ncells; ++i)
        types[i] = 3; /* VTK_LINE */
    else if(sub->dim == 2)
      for(bfam_locidx_t i = 0; i < Ncells; ++i)
        types[i] = 8; /* VTK_PIXEL */
    else if(sub->dim == 3)
      for(bfam_locidx_t i = 0; i < Ncells; ++i)
        types[i] = 11; /* VTK_VOXEL */
    else BFAM_ABORT("cannot handle dim = %d",sub->dim);

    int rval =
      bfam_vtk_write_binary_data(writeCompressed, file, (char*)types,
          typesSize);
    if(rval)
      BFAM_WARNING("Error encoding types");

    bfam_free_aligned(types);
  }
  else
  {
    for(bfam_locidx_t i = 0, sk = 1; i < Ncells; ++i, ++sk)
    {
      if(sub->dim==1) fprintf(file, " 3"); /* VTK_LINE */
      else if(sub->dim==2) fprintf(file, " 8"); /* VTK_PIXEL */
      else if(sub->dim==3) fprintf(file, " 11"); /* VTK_VOXEL */
      else BFAM_ABORT("cannot handle dim = %d",sub->dim);
      if (!(sk % 20) && i != (Ncells - 1))
        fprintf(file, "\n         ");
    }
  }
  fprintf(file, "\n");
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "      </Cells>\n");

  /*
   * Cell Data
   */
  fprintf(file, "      <CellData Scalars=\"time,mpirank,subdomain_id\">\n");
  fprintf(file, "        <DataArray type=\"%s\" Name=\"time\""
           " format=\"%s\">\n", BFAM_REAL_VTK, format);
  fprintf(file, "          ");
  if(writeBinary)
  {
    size_t timesize = Ncells*sizeof(bfam_real_t);
    bfam_real_t *times = bfam_malloc_aligned(timesize);

    for(bfam_locidx_t i = 0; i < Ncells; ++i)
      times[i] = time;

    int rval =
      bfam_vtk_write_binary_data(writeCompressed, file, (char*)times,
          timesize);
    if(rval)
      BFAM_WARNING("Error encoding times");

    bfam_free_aligned(times);
  }
  else
  {
    for(bfam_locidx_t i = 0, sk = 1; i < Ncells; ++i, ++sk)
    {
      fprintf(file, " %"BFAM_REAL_FMTe, time);
      if (!(sk % 8) && i != (Ncells - 1))
        fprintf(file, "\n         ");
    }
  }
  fprintf(file, "\n");
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"%s\" Name=\"mpirank\""
           " format=\"%s\">\n", BFAM_LOCIDX_VTK, format);
  fprintf(file, "          ");
  if(writeBinary)
  {
    size_t ranksSize = Ncells*sizeof(bfam_locidx_t);
    bfam_locidx_t *ranks = bfam_malloc_aligned(ranksSize);

    for(bfam_locidx_t i = 0; i < Ncells; ++i)
      ranks[i] = rank;

    int rval =
      bfam_vtk_write_binary_data(writeCompressed, file, (char*)ranks,
          ranksSize);
    if(rval)
      BFAM_WARNING("Error encoding ranks");

    bfam_free_aligned(ranks);
  }
  else
  {
    for(bfam_locidx_t i = 0, sk = 1; i < Ncells; ++i, ++sk)
    {
      fprintf(file, " %6jd", (intmax_t)rank);
      if (!(sk % 8) && i != (Ncells - 1))
        fprintf(file, "\n         ");
    }
  }
  fprintf(file, "\n");
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"%s\" Name=\"subdomain_id\""
           " format=\"%s\">\n", BFAM_LOCIDX_VTK, format);
  fprintf(file, "          ");
  if(writeBinary)
  {
    size_t idsSize = Ncells*sizeof(bfam_locidx_t);
    bfam_locidx_t *ids = bfam_malloc_aligned(idsSize);

    for(bfam_locidx_t i = 0; i < Ncells; ++i)
      ids[i] = id;

    int rval =
      bfam_vtk_write_binary_data(writeCompressed, file, (char*)ids,
          idsSize);
    if(rval)
      BFAM_WARNING("Error encoding ids");

    bfam_free_aligned(ids);
  }
  else
  {
    for(bfam_locidx_t i = 0, sk = 1; i < Ncells; ++i, ++sk)
    {
      fprintf(file, " %6jd", (intmax_t)id);
      if (!(sk % 8) && i != (Ncells - 1))
        fprintf(file, "\n         ");
    }
  }
  fprintf(file, "\n");
  fprintf(file, "        </DataArray>\n");

  fprintf(file, "      </CellData>\n");

  char pointscalars[BFAM_BUFSIZ];
  bfam_util_strcsl(pointscalars, scalars);

  char pointvectors[BFAM_BUFSIZ];
  bfam_util_strcsl(pointvectors, vectors);

  fprintf(file, "      <PointData Scalars=\"%s\" Vectors=\"%s\">\n",
      pointscalars, pointvectors);

  if(scalars)
  {
    for(size_t s = 0; scalars[s]; ++s)
    {
      bfam_real_t *sdata = bfam_dictionary_get_value_ptr(&subdomain->fields,
          scalars[s]);
      BFAM_ABORT_IF(sdata == NULL, "VTK: Field %s not in subdomain %s",
          scalars[s], subdomain->name);
      if(interp == NULL)
      {
        stor1 = sdata;
      }
      else
      {
        bfam_subdomain_dgx_vtk_interp(K,N_vtk,stor1,sub->N,sdata,interp);
      }

      bfam_vtk_write_real_scalar_data_array(file, scalars[s],
          writeBinary, writeCompressed, Ntotal, stor1);
    }
  }

  if(vectors)
  {
    for(size_t v = 0; vectors[v]; ++v)
    {

      bfam_real_t *v1 =
        bfam_dictionary_get_value_ptr(&subdomain->fields, components[3*v+0]);
      bfam_real_t *v2 =
        bfam_dictionary_get_value_ptr(&subdomain->fields, components[3*v+1]);
      bfam_real_t *v3 =
        bfam_dictionary_get_value_ptr(&subdomain->fields, components[3*v+2]);

      BFAM_ABORT_IF(v1 == NULL, "VTK: Field %s not in subdomain %s",
          components[3*v+0], subdomain->name);
      BFAM_ABORT_IF(v2 == NULL, "VTK: Field %s not in subdomain %s",
          components[3*v+1], subdomain->name);
      BFAM_ABORT_IF(v3 == NULL, "VTK: Field %s not in subdomain %s",
          components[3*v+2], subdomain->name);
      if(interp == NULL)
      {
        stor1 = v1;
        stor2 = v2;
        stor3 = v3;
      }
      else
      {
        bfam_subdomain_dgx_vtk_interp(K,N_vtk,stor1,sub->N,v1,interp);
        bfam_subdomain_dgx_vtk_interp(K,N_vtk,stor2,sub->N,v2,interp);
        bfam_subdomain_dgx_vtk_interp(K,N_vtk,stor3,sub->N,v3,interp);
      }

      bfam_vtk_write_real_vector_data_array(file, vectors[v],
          writeBinary, writeCompressed, Ntotal, stor1, stor2, stor3);
    }
  }

  fprintf(file, "      </PointData>\n");
  fprintf(file, "    </Piece>\n");

  if(interp != NULL)
  {
    bfam_free_aligned(interp);
    bfam_free_aligned(stor1);
    bfam_free_aligned(stor2);
    bfam_free_aligned(stor3);
  }
  return 1;
}

static int
bfam_subdomain_dgx_field_face_add(bfam_subdomain_t *subdomain,
    const char *name)
{
  bfam_subdomain_dgx_t *s = (bfam_subdomain_dgx_t*) subdomain;

  if(bfam_dictionary_get_value_ptr(&s->base.fields_face,name))
    return 1;

  size_t fieldSize = s->Ng[0]*s->Ngp[0]*s->K*sizeof(bfam_real_t);
  bfam_real_t *field = bfam_malloc_aligned(fieldSize);

  int rval = bfam_dictionary_insert_ptr(&s->base.fields_face, name, field);

  BFAM_ASSERT(rval != 1);

  if(rval == 0)
    bfam_free_aligned(field);

  return rval;
}

static int
bfam_subdomain_dgx_field_add(bfam_subdomain_t *subdomain, const char *name)
{
  bfam_subdomain_dgx_t *s = (bfam_subdomain_dgx_t*) subdomain;

  if(bfam_dictionary_get_value_ptr(&s->base.fields,name))
    return 1;

  size_t fieldSize = s->Np*s->K*sizeof(bfam_real_t);
  bfam_real_t *field = bfam_malloc_aligned(fieldSize);
#ifdef BFAM_DEBUG
  for(int i = 0; i < s->Np*s->K;i++) field[i] = bfam_real_nan("");
#endif

  int rval = bfam_dictionary_insert_ptr(&s->base.fields, name, field);

  BFAM_ASSERT(rval != 1);

  if(rval == 0)
    bfam_free_aligned(field);

  return rval;
}

static inline int***
BFAM_APPEND_EXPAND(bfam_subdomain_dgx_gmask_set_,BFAM_DGX_DIMENSION)
         (const int numg, const int N, int *Np, int *Ng, int *Ngp, int inDIM)
{
#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_gmask_set");
  const int DIM = inDIM;
#endif
  BFAM_ABORT_IF(DIM > 3 || DIM < 0,
      "bfam_subdomain_dgx_gmask_set cannot handle dim = %d",DIM);

  if(DIM == 0)
  {
    *Np = 1;
    return NULL;
  }


  /* this could probably be made generic for arbitrary dimensions, but until
   * that's needed... */
  switch(DIM)
  {
    case 1:
      *Np = N+1;

      /* just corners */
      Ng [0]  = 2;
      Ngp[0]  = 1;
      break;

    case 2:
      *Np = (N+1)*(N+1);

      /* edges */
      Ng [0]  = 4;
      Ngp[0]  = N+1;

      /* corners */
      Ng [1]  = 4;
      Ngp[1]  = 1;
      break;

    case 3:
      *Np = (N+1)*(N+1)*(N+1);

      /* faces */
      Ng [0]  = 6;
      Ngp[0]  = (N+1)*(N+1);

      /* edges */
      Ng [1]  = 12;
      Ngp[1]  = N+1;

      /* corners */
      Ng [2]  = 8;
      Ngp[2]  = 1;
      break;

    default:
      BFAM_ABORT("cannot handle dim = %d",DIM);
  }

  int ***gmask = bfam_malloc_aligned(numg * sizeof(int**));
  for(int g = 0; g < numg; g++)
  {
    gmask[g] = bfam_malloc_aligned(Ng[g] * sizeof(int*));
    for(int i = 0; i < Ng[g]; i++)
      gmask[g][i] = bfam_malloc_aligned(Ngp[g] * sizeof(int));
  }

  switch(DIM)
  {
    case 1:
      gmask[0][0][0] = 0;
      gmask[0][1][0] = N;
      break;

    case 2:
      /* edges */
      for(int i = 0; i < N+1; ++i) gmask[0][0][i] = i*(N+1);
      for(int i = 0; i < N+1; ++i) gmask[0][1][i] = (i+1)*(N+1)-1;
      for(int i = 0; i < N+1; ++i) gmask[0][2][i] = i;
      for(int i = 0; i < N+1; ++i) gmask[0][3][i] = (N+1)*N + i;

      /* corners */
      for(int j = 0; j < 2; ++j)
        for(int i = 0; i < 2; ++i)
          gmask[1][i+j*2][0] = i*N + j*(N+1);
      break;

    case 3:
      /* This could all probably be cleaned up... */

      /* faces */
      {
        int n,i,j,k,f=-1;

        n = 0; i = 0; f++;
        for(k = 0; k < N+1; k++) for(j = 0; j < N+1; j++)
          gmask[0][f][n++] = i+j*(N+1)+k*(N+1)*(N+1);

        n = 0; i = N; f++;
        for(k = 0; k < N+1; k++) for(j = 0; j < N+1; j++)
          gmask[0][f][n++] = i+j*(N+1)+k*(N+1)*(N+1);

        n = 0; j = 0; f++;
        for(k = 0; k < N+1; k++) for(i = 0; i < N+1; i++)
          gmask[0][f][n++] = i+j*(N+1)+k*(N+1)*(N+1);

        n = 0; j = N; f++;
        for(k = 0; k < N+1; k++) for(i = 0; i < N+1; i++)
          gmask[0][f][n++] = i+j*(N+1)+k*(N+1)*(N+1);

        n = 0; k = 0; f++;
        for(j = 0; j < N+1; j++) for(i = 0; i < N+1; i++)
          gmask[0][f][n++] = i+j*(N+1)+k*(N+1)*(N+1);

        n = 0; k = N; f++;
        for(j = 0; j < N+1; j++) for(i = 0; i < N+1; i++)
          gmask[0][f][n++] = i+j*(N+1)+k*(N+1)*(N+1);
      }

      /* edges */
      {
        int n,i,j,k,e = 0;

        for(k = 0; k < N+1;k+=N)
          for(j = 0; j < N+1;j+=N)
          {
            n = 0;
            for(i = 0; i < N+1; i++)
              gmask[1][e][n++] = i+j*(N+1)+k*(N+1)*(N+1);
            e++;
          }
        for(k = 0; k < N+1;k+=N)
          for(i = 0; i < N+1;i+=N)
          {
            n = 0;
            for(j = 0; j < N+1; j++)
              gmask[1][e][n++] = i+j*(N+1)+k*(N+1)*(N+1);
            e++;
          }
        for(j = 0; j < N+1;j+=N)
          for(i = 0; i < N+1;i+=N)
          {
            n = 0;
            for(k = 0; k < N+1; k++)
              gmask[1][e][n++] = i+j*(N+1)+k*(N+1)*(N+1);
            e++;
          }
      }

      /* corners */
      for(int k = 0, c = 0; k < N+1; k+=N)
        for(int j = 0; j < N+1;j+=N)
          for(int i = 0; i < N+1;i+=N)
            gmask[2][c++][0] = i+j*(N+1)+k*(N+1)*(N+1);

      break;

    default:
      BFAM_ABORT("cannot handle dim = %d",DIM);
  }

  return gmask;
}

static void
bfam_subdomain_dgx_buildmaps(bfam_locidx_t K, int Np, int Nfp, int Nfaces,
   const bfam_locidx_t *EToE, const int8_t *EToF, int ***gmask,
   bfam_locidx_t *restrict vmapP, bfam_locidx_t *restrict vmapM)
{
  for(bfam_locidx_t k1 = 0, sk = 0; k1 < K; ++k1)
  {
    for(int8_t f1 = 0; f1 < Nfaces; ++f1)
    {
      bfam_locidx_t k2 = EToE[Nfaces * k1 + f1];
      int8_t        f2 = EToF[Nfaces * k1 + f1] % Nfaces;
      int8_t        o  = EToF[Nfaces * k1 + f1] / Nfaces;

      for(int n = 0; n < Nfp; ++n)
      {
        vmapM[sk + n] = Np * k1 + gmask[0][f1][n];

        if(o)
          vmapP[sk + n] = Np * k2 + gmask[0][f2][Nfp-1-n];
        else
          vmapP[sk + n] = Np * k2 + gmask[0][f2][n];
      }

      sk += Nfp;
    }
  }
}

/* set all subdomain values to something logical */
static void
bfam_subdomain_dgx_null_all_values(bfam_subdomain_dgx_t *sub)
{
  sub->K       =  0;
  sub->N       = -1;
  sub->Np      =  0;
  sub->Ngp     = NULL;
  sub->numg    =  0;
  sub->Ng      = NULL;
  sub->r       = NULL;
  sub->w       = NULL;
  sub->wi      = NULL;
  sub->Dr      = NULL;
  sub->V       = NULL;
  sub->K       =  0;
  sub->vmapM   = NULL;
  sub->vmapP   = NULL;
  sub->gmask   = NULL;
}

void
BFAM_APPEND_EXPAND(bfam_subdomain_dgx_init_,BFAM_DGX_DIMENSION)(
                              bfam_subdomain_dgx_t *subdomain,
                        const bfam_locidx_t         id,
                        const char                 *name,
                        const int                   N,
                        const bfam_locidx_t         Nv,
                        const int                   num_Vi,
                        const bfam_long_real_t    **Vi,
                        const bfam_locidx_t         K,
                        const bfam_locidx_t        *EToV,
                        const bfam_locidx_t        *EToE,
                        const int8_t               *EToF,
                        const int                   inDIM)
{
#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_init");
  const int DIM = inDIM;
#endif

  BFAM_ASSERT(DIM == inDIM);
  BFAM_ABORT_IF(DIM < 0, "dimension %d is not possible in bfam",DIM);
  BFAM_ABORT_IF(DIM == 0 && N != 0,
                "if DIM < 1 then N must be zero (i.e., constant");

  bfam_subdomain_init(&subdomain->base, id, name);
  bfam_subdomain_add_tag(&subdomain->base, "_subdomain_dgx");
  char dim_str[BFAM_BUFSIZ];
  snprintf(dim_str,BFAM_BUFSIZ,"_dimension_%d",DIM);
  bfam_subdomain_add_tag(&subdomain->base, dim_str);

  bfam_subdomain_dgx_null_all_values(subdomain);

  subdomain->dim = DIM;

  subdomain->base.free =
              BFAM_APPEND_EXPAND(bfam_subdomain_dgx_free_,BFAM_DGX_DIMENSION);
  subdomain->base.vtk_write_vtu_piece =
    bfam_subdomain_dgx_vtk_write_vtu_piece;
  subdomain->base.field_add = bfam_subdomain_dgx_field_add;
  subdomain->base.field_face_add = bfam_subdomain_dgx_field_face_add;
  subdomain->base.field_init = bfam_subdomain_dgx_field_init;

  subdomain->numg = DIM;
  const int numg = subdomain->numg;

  int *Ng  = NULL;
  if(numg > 0) Ng  = bfam_malloc_aligned(sizeof(int)*numg);
  subdomain->Ng = Ng;

  int *Ngp  = NULL;
  if(numg > 0) Ngp = bfam_malloc_aligned(sizeof(int)*numg);
  subdomain->Ngp = Ngp;

  subdomain->gmask =
    BFAM_APPEND_EXPAND(bfam_subdomain_dgx_gmask_set_,BFAM_DGX_DIMENSION)(
                       numg, N, &subdomain->Np, Ng, Ngp, DIM);

  const int Np = subdomain->Np;

  if(DIM > 0)
  {
    const int Nrp = N+1;
    bfam_long_real_t *lr, *lw;
    lr = bfam_malloc_aligned(Nrp*sizeof(bfam_long_real_t));
    lw = bfam_malloc_aligned(Nrp*sizeof(bfam_long_real_t));

    bfam_jacobi_gauss_lobatto_quadrature(0, 0, N, lr, lw);

    bfam_long_real_t **lxi =
      bfam_malloc_aligned(num_Vi*sizeof(bfam_long_real_t*));

    for(int i = 0;i < num_Vi;i++)
      lxi[i] = bfam_malloc_aligned(K*Np*sizeof(bfam_long_real_t));

    /* Loop over all the elements and set up the grid*/
    int Ncorners = Ng[numg-1];
    for(bfam_locidx_t k = 0; k < K; ++k)
    {
      const bfam_locidx_t *v = EToV+Ncorners*k;
      bfam_long_real_t w[Ncorners];

      if(DIM == 1)
        for(int n = 0; n < Nrp; ++n)
        {
          int offset = n;
          w[0] = 1-lr[n];
          w[1] = 1+lr[n];
          for(int i = 0; i < num_Vi;i++)
          {
            lxi[i][Np*k + offset] = 0;
            for(int c = 0; c < Ncorners; c++)
              lxi[i][Np*k + offset] += w[c]*Vi[i][v[c]];
            lxi[i][Np*k + offset] /= Ncorners;
          }
        }
      else if(DIM == 2)
        for(int n = 0; n < Nrp; ++n)
          for(int m = 0; m < Nrp; ++m)
          {
            int offset = n*Nrp+m;
            w[0] = (1-lr[m])*(1-lr[n]);
            w[1] = (1+lr[m])*(1-lr[n]);
            w[2] = (1-lr[m])*(1+lr[n]);
            w[3] = (1+lr[m])*(1+lr[n]);
            for(int i = 0; i < num_Vi;i++)
            {
              lxi[i][Np*k + offset] = 0;
              for(int c = 0; c < Ncorners; c++)
                lxi[i][Np*k + offset] += w[c]*Vi[i][v[c]];
              lxi[i][Np*k + offset] /= Ncorners;
            }
        }
      else if(DIM == 3)
        for(int n = 0; n < Nrp; ++n)
          for(int m = 0; m < Nrp; ++m)
            for(int l = 0; l < Nrp; ++l)
            {
              int offset = n*Nrp*Nrp+m*Nrp+l;
              w[0] = (1-lr[l])*(1-lr[m])*(1-lr[n]);
              w[1] = (1+lr[l])*(1-lr[m])*(1-lr[n]);
              w[2] = (1-lr[l])*(1+lr[m])*(1-lr[n]);
              w[3] = (1+lr[l])*(1+lr[m])*(1-lr[n]);
              w[4] = (1-lr[l])*(1-lr[m])*(1+lr[n]);
              w[5] = (1+lr[l])*(1-lr[m])*(1+lr[n]);
              w[6] = (1-lr[l])*(1+lr[m])*(1+lr[n]);
              w[7] = (1+lr[l])*(1+lr[m])*(1+lr[n]);
              for(int i = 0; i < num_Vi;i++)
              {
                lxi[i][Np*k + offset] = 0;
                for(int c = 0; c < Ncorners; c++)
                  lxi[i][Np*k + offset] += w[c]*Vi[i][v[c]];
                lxi[i][Np*k + offset] /= Ncorners;
              }
            }
      else BFAM_ABORT("not setup of dim = %d",DIM);

    }

    subdomain->V = bfam_malloc_aligned(Nrp*Nrp*sizeof(bfam_long_real_t));
    bfam_long_real_t *restrict V = subdomain->V;

    bfam_jacobi_p_vandermonde(0, 0, N, Nrp, lr, V);

    bfam_long_real_t *restrict D =
      bfam_malloc_aligned(Nrp*Nrp*sizeof(bfam_long_real_t));

    bfam_jacobi_p_differentiation(0, 0, N, Nrp, lr, V, D);

    bfam_long_real_t *restrict M =
      bfam_malloc_aligned(Nrp*Nrp*sizeof(bfam_long_real_t));

    bfam_jacobi_p_mass(0, 0, N, V, M);

    bfam_long_real_t **lJrx =
      bfam_malloc_aligned(num_Vi*DIM*sizeof(bfam_long_real_t*));
    for(int n = 0; n < num_Vi*DIM; n++)
      lJrx[n] = bfam_malloc_aligned(K*Np*sizeof(bfam_long_real_t));

    bfam_long_real_t  *lJ  = NULL;
    bfam_long_real_t **lni  = NULL;
    bfam_long_real_t  *lsJ = NULL;
    if(DIM == num_Vi)
    {
      lJ = bfam_malloc_aligned(K*Np*sizeof(bfam_long_real_t));

      /* Ng[0] = number of faces, Ngp[0] = number of face points */
      lni = bfam_malloc_aligned(num_Vi*sizeof(bfam_long_real_t*));
      for(int n = 0; n < num_Vi; n++)
        lni[n] = bfam_malloc_aligned(K*Ng[0]*Ngp[0]*sizeof(bfam_long_real_t));

      lsJ = bfam_malloc_aligned(K*Ng[0]*Ngp[0]*sizeof(bfam_long_real_t));
    }

    bfam_subdomain_dgx_geo(N, K, Np, subdomain->gmask, Ng, Ngp, num_Vi,
        lxi, D, lJrx, lJ, lni, lsJ, DIM);

    subdomain->r  = bfam_malloc_aligned(Nrp*sizeof(bfam_real_t));
    subdomain->w  = bfam_malloc_aligned(Nrp*sizeof(bfam_real_t));
    subdomain->wi = bfam_malloc_aligned(Nrp*sizeof(bfam_real_t));

    for(int n = 0; n<Nrp; ++n)
    {
      subdomain->r[n]  = (bfam_real_t) lr[n];
      subdomain->w[n]  = (bfam_real_t) lw[n];
      subdomain->wi[n] = (bfam_real_t) (1.0l/lw[n]);
    }

    subdomain->K = K;

    /* store the volume stuff */
    subdomain->N = N;
    /* store the grid */
    for(int i = 0; i < num_Vi; i++)
    {
      char name[BFAM_BUFSIZ];
      snprintf(name,BFAM_BUFSIZ,"_grid_x%d",i);
      int rval = bfam_subdomain_dgx_field_add(&subdomain->base, name);
      BFAM_ABORT_IF_NOT(rval == 2, "Error adding %s",name);
      bfam_real_t *restrict xi =
        bfam_dictionary_get_value_ptr(&subdomain->base.fields, name);
      for(int n = 0; n < K*Np; ++n)
      {
        xi[n] = (bfam_real_t) lxi[i][n];
      }
    }

    /* store the metric stuff */
    if(lJ)
    {
      for(int d = 0; d < DIM; d++)
      {
        for(int v = 0; v < num_Vi; v++)
        {
          char name[BFAM_BUFSIZ];
          snprintf(name,BFAM_BUFSIZ,"_grid_Jr%dx%d",d,v);
          int rval = bfam_subdomain_dgx_field_add(&subdomain->base, name);
          BFAM_ABORT_IF_NOT(rval == 2, "Error adding %s",name);
          bfam_real_t *restrict Jrx =
            bfam_dictionary_get_value_ptr(&subdomain->base.fields, name);
          for(int n = 0; n < K*Np; ++n)
            Jrx[n] = (bfam_real_t) lJrx[v+d*num_Vi][n];
        }
      }
      {
        char name[] = "_grid_J";
        int rval = bfam_subdomain_dgx_field_add(&subdomain->base, name);
        BFAM_ABORT_IF_NOT(rval == 2, "Error adding %s",name);
        bfam_real_t *restrict J =
          bfam_dictionary_get_value_ptr(&subdomain->base.fields, name);
        for(int n = 0; n < K*Np; ++n)
          J[n] = (bfam_real_t) lJ[n];
      }
      BFAM_ASSERT(lni != NULL && lsJ != NULL);
      for(int v = 0; v < num_Vi;v++)
      {
        char name[BFAM_BUFSIZ];
        snprintf(name,BFAM_BUFSIZ,"_grid_nx%d",v);
        int rval = bfam_subdomain_dgx_field_face_add(&subdomain->base, name);
        BFAM_ABORT_IF_NOT(rval == 2, "Error adding %s",name);
        bfam_real_t *restrict ni =
          bfam_dictionary_get_value_ptr(&subdomain->base.fields_face, name);
        for(int n = 0; n < K*Ng[0]*Ngp[0]; ++n)
          ni[n] = (bfam_real_t) lni[v][n];
      }
      {
        char name[] = "_grid_sJ";
        int rval = bfam_subdomain_dgx_field_face_add(&subdomain->base, name);
        BFAM_ABORT_IF_NOT(rval == 2, "Error adding %s",name);
        bfam_real_t *restrict sJ =
          bfam_dictionary_get_value_ptr(&subdomain->base.fields_face, name);
        for(int n = 0; n < K*Ng[0]*Ngp[0]; ++n)
          sJ[n] = (bfam_real_t) lsJ[n];
      }
    }
    else
    {
      /* In this Jrx really has dx/dr */
      for(int v = 0; v < num_Vi; v++)
      {
        for(int d = 0; d < DIM; d++)
        {
          char name[BFAM_BUFSIZ];
          snprintf(name,BFAM_BUFSIZ,"_grid_x%dr%d",v,d);
          int rval = bfam_subdomain_dgx_field_add(&subdomain->base, name);
          BFAM_ABORT_IF_NOT(rval == 2, "Error adding %s",name);
          bfam_real_t *restrict xr =
            bfam_dictionary_get_value_ptr(&subdomain->base.fields, name);
          for(int n = 0; n < K*Np; ++n)
            xr[n] = (bfam_real_t) lJrx[d+v*DIM][n];
        }
      }
      BFAM_ASSERT(lni == NULL && lsJ == NULL);
    }

    /* store the face stuff */

    subdomain->Dr = bfam_malloc_aligned(Nrp * Nrp * sizeof(bfam_real_t));
    for(int n = 0; n < Nrp*Nrp; ++n) subdomain->Dr[n] = (bfam_real_t) D[n];

    subdomain->vmapP = bfam_malloc_aligned(K*Ngp[0]*Ng[0]*sizeof(bfam_locidx_t));
    subdomain->vmapM = bfam_malloc_aligned(K*Ngp[0]*Ng[0]*sizeof(bfam_locidx_t));

    bfam_subdomain_dgx_buildmaps(K, Np, Ngp[0], Ng[0], EToE, EToF,
        subdomain->gmask, subdomain->vmapP, subdomain->vmapM);

    /* free stuff */
    if(lsJ) bfam_free_aligned(lsJ);
    if(lni)
    {
      for(int n = 0; n < num_Vi; n++)
        bfam_free_aligned(lni[n]);
      bfam_free_aligned(lni);
    }

    if(lJ) bfam_free_aligned(lJ);

    for(int n = 0; n < num_Vi*DIM; n++)
      bfam_free_aligned(lJrx[n]);
    bfam_free_aligned(lJrx);

    bfam_free_aligned(D);
    bfam_free_aligned(M);

    bfam_free_aligned(lr);
    bfam_free_aligned(lw);

    for(int i = 0;i < num_Vi;i++) bfam_free_aligned(lxi[i]);
    bfam_free_aligned(lxi);
  }
}

bfam_subdomain_dgx_t*
BFAM_APPEND_EXPAND(bfam_subdomain_dgx_new_,BFAM_DGX_DIMENSION)(
                       const bfam_locidx_t      id,
                       const char              *name,
                       const int                N,
                       const bfam_locidx_t      Nv,
                       const int                num_Vi,
                       const bfam_long_real_t **Vi,
                       const bfam_locidx_t      K,
                       const bfam_locidx_t     *EToV,
                       const bfam_locidx_t     *EToE,
                       const int8_t            *EToF,
                       const int                inDIM)
{
#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_new");
  const int DIM = inDIM;
#endif
  BFAM_ASSERT(DIM == inDIM);

  bfam_subdomain_dgx_t* newSubdomain =
    bfam_malloc(sizeof(bfam_subdomain_dgx_t));

  BFAM_APPEND_EXPAND(bfam_subdomain_dgx_init_,BFAM_DGX_DIMENSION)(
                         newSubdomain, id, name, N, Nv, num_Vi, Vi, K, EToV,
                         EToE, EToF,DIM);
  return newSubdomain;
}

static int
bfam_subdomain_dgx_free_fields(const char * key, void *val,
    void *arg)
{
  bfam_free_aligned(val);

  return 1;
}

void
BFAM_APPEND_EXPAND(bfam_subdomain_dgx_free_,BFAM_DGX_DIMENSION)(
    bfam_subdomain_t *thisSubdomain)
{
  bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t*) thisSubdomain;

  bfam_dictionary_allprefixed_ptr(&sub->base.fields,"",
      &bfam_subdomain_dgx_free_fields,NULL);
  bfam_dictionary_allprefixed_ptr(&sub->base.fields_p,"",
      &bfam_subdomain_dgx_free_fields,NULL);
  bfam_dictionary_allprefixed_ptr(&sub->base.fields_m,"",
      &bfam_subdomain_dgx_free_fields,NULL);
  bfam_dictionary_allprefixed_ptr(&sub->base.fields_face,"",
      &bfam_subdomain_dgx_free_fields,NULL);

  bfam_subdomain_free(thisSubdomain);

  if(sub->Dr) bfam_free_aligned(sub->Dr); sub->Dr = NULL;
  if(sub->V ) bfam_free_aligned(sub->V ); sub->V  = NULL;

  if(sub->r ) bfam_free_aligned(sub->r ); sub->r  = NULL;
  if(sub->w ) bfam_free_aligned(sub->w ); sub->w  = NULL;
  if(sub->wi) bfam_free_aligned(sub->wi); sub->wi = NULL;

  if(sub->gmask)
  {
    for(int g = 0; g < sub->numg; g++)
    {
      for(int i = 0; i < sub->Ng[g]; i++)
        bfam_free_aligned(sub->gmask[g][i]);
      bfam_free_aligned(sub->gmask[g]);
    }
    bfam_free_aligned(sub->gmask);
    sub->gmask = NULL;
  }
  if(sub->Ng ) bfam_free_aligned(sub->Ng ); sub->Ng  = NULL;
  if(sub->Ngp) bfam_free_aligned(sub->Ngp); sub->Ngp = NULL;

  if(sub->vmapP) bfam_free_aligned(sub->vmapP); sub->vmapP = NULL;
  if(sub->vmapM) bfam_free_aligned(sub->vmapM); sub->vmapM = NULL;

  bfam_subdomain_dgx_null_all_values(sub);
}

bfam_subdomain_dgx_t*
BFAM_APPEND_EXPAND(bfam_subdomain_dgx_glue_new_,BFAM_DGX_DIMENSION)(
                             const bfam_locidx_t              id,
                             const char                      *name,
                             const int                        N_m,
                             const bfam_locidx_t              rank_m,
                             bfam_subdomain_dgx_t            *sub_m,
                             bfam_locidx_t                   *ktok_m,
                             const bfam_locidx_t              K,
                             const int                        inDIM)
{
#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_glue_new");
  const int DIM = inDIM;
#endif
  BFAM_ASSERT(DIM == inDIM);

  bfam_subdomain_dgx_t* newSubdomain =
    bfam_malloc(sizeof(bfam_subdomain_dgx_t));

  BFAM_APPEND_EXPAND(bfam_subdomain_dgx_glue_init_,BFAM_DGX_DIMENSION)(
      newSubdomain,id,name,N_m,rank_m,sub_m,ktok_m,K,DIM);
  return newSubdomain;
}

void
BFAM_APPEND_EXPAND(bfam_subdomain_dgx_glue_init_,BFAM_DGX_DIMENSION)(
                             bfam_subdomain_dgx_t            *subdomain,
                             const bfam_locidx_t              id,
                             const char                      *name,
                             const int                        N_m,
                             const bfam_locidx_t              rank_m,
                             bfam_subdomain_dgx_t            *sub_m,
                             bfam_locidx_t                   *ktok_m,
                             const bfam_locidx_t              K,
                             const int                        inDIM)
{
#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_glue_init");
  const int DIM = inDIM;
#endif
  BFAM_ASSERT(DIM == inDIM);

  BFAM_ABORT_IF(DIM < 0, "dimension %d is not possible in bfam",DIM);

  bfam_subdomain_init(&subdomain->base, id, name);
  bfam_subdomain_add_tag(&subdomain->base, "_subdomain_dgx");
  bfam_subdomain_add_tag(&subdomain->base, "_subdomain_dgx_glue");
  char dim_str[BFAM_BUFSIZ];
  snprintf(dim_str,BFAM_BUFSIZ,"_dimension_%d",DIM);
  bfam_subdomain_add_tag(&subdomain->base, dim_str);

  bfam_subdomain_dgx_null_all_values(subdomain);
}
