#include <bfam.h>
#include <bfam_subdomain_sbp.h>

void
zero_buffer(bfam_real_t* field,bfam_subdomain_sbp_t *sub)
{
  bfam_locidx_t ix[6] = {0,0,0,0,0,0};
  bfam_locidx_t bx[6] = {0,0,0,0,0,0};
  for(int d = 0;d < sub->dim;d++)
  {
    ix[2*d  ] = sub->buf_sz[2*d];
    ix[2*d+1] = sub->buf_sz[2*d]+sub->Nl[d];
    bx[2*d+1] = sub->buf_sz[2*d]+sub->Nl[d]+sub->buf_sz[2*d+1];
  }

  int nx = bx[1]-bx[0]+1;
  int ny = bx[3]-bx[2]+1;
  for(int k = 0; k <= bx[5]; k++)
    for(int j = 0; j <= bx[3]; j++)
      for(int i = 0; i <= bx[1]; i++)
        if(i<ix[0] && ix[2]<j && j<ix[3] && ix[4]<k && k<ix[5])
          field[i+nx*(j+ny*k)] = 0;
        else if(ix[1]<i && ix[2]<j && j<ix[3] && ix[4]<k && k<ix[5])
          field[i+nx*(j+ny*k)] = 1;
        else if(ix[0]<i && i<ix[1] && j<ix[2] && ix[4]<k && k<ix[5])
          field[i+nx*(j+ny*k)] = 2;
        else if(ix[0]<i && i<ix[1] && ix[3]<j && ix[4]<k && k<ix[5])
          field[i+nx*(j+ny*k)] = 3;
        else if(ix[0]<i && i<ix[1] && ix[2]<j && j<ix[3] && k<ix[4])
          field[i+nx*(j+ny*k)] = 4;
        else if(ix[0]<i && i<ix[1] && ix[2]<j && j<ix[3] && ix[5]<k)
          field[i+nx*(j+ny*k)] = 5;
}


static void
poly1_field(bfam_locidx_t npoints, bfam_real_t time, bfam_real_t *restrict x,
    bfam_real_t *restrict y, bfam_real_t *restrict z,
    struct bfam_subdomain *s, void *arg, bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(field, 32);

  if(z == NULL)
    for(bfam_locidx_t n=0; n < npoints; ++n)
      field[n] = -x[n] - y[n];
  else
    for(bfam_locidx_t n=0; n < npoints; ++n)
      field[n] = -x[n] - y[n] - z[n];

  zero_buffer(field,(bfam_subdomain_sbp_t*)s);
}

static void
poly1_field_check(bfam_locidx_t npoints, bfam_real_t time,
    bfam_real_t *restrict x, bfam_real_t *restrict y, bfam_real_t *restrict z,
    struct bfam_subdomain *s, void *arg, bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(field, 32);

  if(z == NULL)
    for(bfam_locidx_t n=0; n < npoints; ++n)
      BFAM_ABORT_IF_NOT(field[n] == -x[n] - y[n],
          "poly1_field_check failed on %d got %f expected %f",
          n, field[n],-x[n]-y[n]);
  else
    for(bfam_locidx_t n=0; n < npoints; ++n)
      BFAM_ABORT_IF_NOT(field[n] == -x[n] - y[n] - z[n],
          "poly1_field_check failed on %d got %f expected %f",
          n, field[n],-x[n]-y[n]-z[n]);

}

static void
poly2_field(bfam_locidx_t npoints, bfam_real_t time, bfam_real_t *restrict x,
    bfam_real_t *restrict y, bfam_real_t *restrict z,
    struct bfam_subdomain *s, void *arg, bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(field, 32);

  if(z == NULL)
    for(bfam_locidx_t n=0; n < npoints; ++n)
      field[n] = x[n] + y[n];
  else
    for(bfam_locidx_t n=0; n < npoints; ++n)
      field[n] = x[n] + y[n] + z[n];

  zero_buffer(field,(bfam_subdomain_sbp_t*)s);
}

static void
poly2_field_check(bfam_locidx_t npoints, bfam_real_t time,
    bfam_real_t *restrict x, bfam_real_t *restrict y, bfam_real_t *restrict z,
    struct bfam_subdomain *s, void *arg, bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(field, 32);

  if(z == NULL)
    for(bfam_locidx_t n=0; n < npoints; ++n)
      BFAM_ABORT_IF_NOT(field[n] == x[n] + y[n],
          "poly2_field_check failed on %d got %f expected %f",
          n,field[n],x[n]+y[n]);
  else
    for(bfam_locidx_t n=0; n < npoints; ++n)
      BFAM_ABORT_IF_NOT(field[n] == x[n] + y[n] + z[n],
          "poly2_field_check failed on %d got %f expected %f",
          n,field[n],x[n]+y[n]+z[n]);

  zero_buffer(field,(bfam_subdomain_sbp_t*)s);
}

static void
poly3_field(bfam_locidx_t npoints, bfam_real_t time, bfam_real_t *restrict x,
    bfam_real_t *restrict y, bfam_real_t *restrict z,
    struct bfam_subdomain *s, void *arg, bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(field, 32);

  BFAM_ASSERT(z != NULL);
  for(bfam_locidx_t n=0; n < npoints; ++n)
    field[n] = x[n] + y[n] - z[n];

  zero_buffer(field,(bfam_subdomain_sbp_t*)s);
}


static void
poly3_field_check(bfam_locidx_t npoints, bfam_real_t time,
    bfam_real_t *restrict x, bfam_real_t *restrict y, bfam_real_t *restrict z,
    struct bfam_subdomain *s, void *arg, bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(field, 32);

  BFAM_ASSERT(z != NULL);
  for(bfam_locidx_t n=0; n < npoints; ++n)
    BFAM_ABORT_IF_NOT(field[n] == x[n] + y[n] - z[n],
        "poly3_field_check failed on %d got %f expected %f",
        n,field[n],x[n]+y[n]-z[n]);

  zero_buffer(field,(bfam_subdomain_sbp_t*)s);
}


void
simple_partition_1d(bfam_locidx_t *Nl, bfam_gloidx_t *gx, bfam_locidx_t *Nb,
    bfam_gloidx_t N, bfam_locidx_t size, bfam_locidx_t rank,
    bfam_locidx_t bufsz)
{
  Nl[0] = (N+1)/size;
  gx[0] = Nl[0]*rank;

  bfam_locidx_t rem = (N+1) % size;
  if(rank >= rem)
  {
    Nl[0] -= 1;
    gx[0] += rem;
  }
  else
  {
    gx[0] += rank;
  }
  if(rank == 0) Nb[0] = 0;
  else Nb[0] = bufsz;
  if(rank == size-1) Nb[1] = 0;
  else Nb[1] = bufsz;

  BFAM_ABORT_IF(Nl[0] < bufsz, "fewer points than buffer size: %d %d",Nl[0],bufsz);
}

/** Simple domain partitioner for an sbp multi-block
 *
 * \param [out]     Nl          local ending index (size dim)
 * \param [out]     gx          my global starting index (size dim)
 * \param [out]     Nb          number of buffer cells in each dimension
 *                              (size 2*dim)
 * \param [out]     face_neigh  ranks of my face neighbors (size 2*dim)
 *                              \note not used is NULL
 * \param [out]     bx          my location in the subdomain partitioning
 * \param [in]      pd          processor dimension / layout (size dim)
 * \param [in]      N           global ending index
 * \param [in]      size_in     number of pieces to partion into
 * \param [in]      rank        my number in the partiion
 * \param [in]      bufsz       size of the buffer to create is necessary
 * \param[in]       dim         dimension of the of the parition
 *
 * \note if face_
 */
void
simple_partition(bfam_locidx_t *Nl, bfam_gloidx_t *gx, bfam_locidx_t *Nb,
    bfam_locidx_t *face_neigh, bfam_locidx_t *bx, const int *pd,
    const bfam_gloidx_t *N, bfam_locidx_t rank,
    bfam_locidx_t bufsz,int dim)
{
  bfam_locidx_t tmp_rank = rank;
  if(dim > 2)
  {
    bx[2] = tmp_rank/(pd[0]*pd[1]);
    simple_partition_1d(&Nl[2],&gx[2],&Nb[2*2],N[2],pd[2],
        bx[2],bufsz);
    tmp_rank = tmp_rank % (pd[0]*pd[1]);
  }

  if(dim > 1)
  {
    bx[1] = tmp_rank/pd[0];
    simple_partition_1d(&Nl[1],&gx[1],&Nb[1*2],N[1],pd[1],bx[1],bufsz);
  }

  bx[0] = tmp_rank%pd[0];
  simple_partition_1d(&Nl[0],&gx[0],&Nb[0*2],N[0],pd[0],bx[0],bufsz);

  /* now set up the neighboring processor info */
  if(face_neigh != NULL)
  {
    bfam_locidx_t offset = 1;
    for(int d = 0;d<dim;d++)
    {
      face_neigh[2*d] = rank;
      if(bx[d] != 0) face_neigh[2*d  ] -= offset;

      face_neigh[2*d+1] = rank;
      if(bx[d] != pd[d]-1) face_neigh[2*d+1] += offset;

      offset *= pd[d];
    }
  }
}

void
simple_load_balance(bfam_locidx_t *procs,const bfam_gloidx_t *N,
    int mpi_size,int num_blocks,int dim)
{
  if(mpi_size <= num_blocks)
  {
    for(int b = 0; b < num_blocks; b++)
    {
      procs[2*b  ] = 0;      // starting rank
      procs[2*b+1] = mpi_size-1; // ending rank
    }
  }
  else
  {
    /* get cost per block */
    bfam_real_t cost[num_blocks];
    bfam_locidx_t num_procs[num_blocks];
    for(int b = 0; b < num_blocks; b++)
    {
      num_procs[b] = 1;
      cost[b] = 1;
      for(int d = 0;d < dim;d++) cost[b] *= N[dim*b+1];
    }

    /* figured out number of procs per block */
    for(int p = num_blocks; p < mpi_size; p++)
    {
      bfam_locidx_t max_blck = 0;
      bfam_real_t   max_cost = cost[max_blck]/num_procs[max_blck];
      for(int b = 1; b < num_blocks; b++)
      {
        bfam_real_t test = cost[b]/num_procs[b];
        if(test > max_cost)
        {
          max_cost = test;
          max_blck = b;
        }
      }
      num_procs[max_blck]++;
    }

    bfam_locidx_t lp = -1;
    for(int b = 0; b < num_blocks; b++)
    {
      procs[2*b  ] = lp+1;
      procs[2*b+1] = lp+num_procs[b];
      lp           = procs[2*b+1];
    }
  }
}


void
setup_subdomains(bfam_domain_t *domain,
    int dim, int num_blocks, bfam_locidx_t bufsz,int rank,const char **names,
    const bfam_locidx_t *procs,
    const bfam_locidx_t *EToV,const bfam_locidx_t *EToE,const int8_t *EToF,
    const bfam_gloidx_t *N, const bfam_long_real_t* Vx,
    const bfam_long_real_t *Vy, const bfam_long_real_t *Vz)
{
  const int num_corners = (int)pow(2,dim);
  const int num_face = 2*dim;
  bfam_long_real_t *x = NULL;
  bfam_long_real_t *y = NULL;
  bfam_long_real_t *z = NULL;
  if(dim > 0) x = bfam_malloc_aligned(num_corners*sizeof(bfam_long_real_t));
  if(dim > 1) y = bfam_malloc_aligned(num_corners*sizeof(bfam_long_real_t));
  if(dim > 2) z = bfam_malloc_aligned(num_corners*sizeof(bfam_long_real_t));
  for(int b = 0; b < num_blocks;b++)
  {
    if(procs[2*b] <= rank && rank <= procs[2*b+1])
    {
      bfam_locidx_t Nl[  dim];
      bfam_gloidx_t gx[  dim];
      bfam_locidx_t Nb[2*dim];
      for(int d = 0;d<dim;d++)
      {
        Nl[d] = -1;
        gx[d] = -1;
        Nb[2*d] = 0;
        Nb[2*d+1] = 0;
      }

      bfam_locidx_t l_rank = rank-procs[2*b];
      bfam_locidx_t l_size = procs[2*b+1]-procs[2*b]+1;
      bfam_locidx_t face_neigh[2*dim];
      bfam_locidx_t bx[2*dim];

      int pd[dim];
      for(int d = 0; d < dim; d++) pd[d] = 0;

      BFAM_MPI_CHECK(MPI_Dims_create(l_size, dim, pd));

      simple_partition(Nl,gx,Nb,face_neigh,bx,pd,&N[dim*b],
          l_rank,bufsz,dim);

      for(int ix = 0; ix < num_corners; ix++)
      {
        int V = EToV[num_corners*b+ix];
        if(x != NULL) x[ix] = Vx[V];
        if(y != NULL) y[ix] = Vy[V];
        if(z != NULL) z[ix] = Vz[V];
      }

      /* make pd really the local top index */
      for(int d = 0; d < dim; d++) pd[d]--;

      /* add subdomain */
      bfam_subdomain_sbp_t *sub =
        bfam_subdomain_sbp_new(b,bx,pd,names[b],dim,
            &N[dim*b],Nl,Nb,gx,x,y,z);
      bfam_subdomain_add_tag((bfam_subdomain_t*)sub,"_volume");
      bfam_domain_add_subdomain(domain,(bfam_subdomain_t*)sub);

      /* add subdomain glue grid */
      for(int d = 0; d < dim; d++)
      {
        for(int face = 2*d;face<2*d+2;face++)
        {
          char glue_name[BFAM_BUFSIZ];
          snprintf(glue_name,BFAM_BUFSIZ,"%s_face_%d",names[b],face);
          bfam_locidx_t glue_id = b+num_blocks*face;
          if(face_neigh[face] != l_rank)
          {
            bfam_subdomain_sbp_intra_glue_t* glue =
              bfam_subdomain_sbp_intra_glue_new(glue_id,glue_name,
                  rank,face_neigh[face]+procs[2*b],sub,face);
            bfam_subdomain_add_tag((bfam_subdomain_t*)glue,"_intra_glue");
            bfam_domain_add_subdomain(domain,(bfam_subdomain_t*)glue);
          }
          else
          {
            bfam_locidx_t neigh  = EToE[b*num_face+face];
            int8_t n_face = EToF[b*num_face+face]%num_face;
            int8_t orient = EToF[b*num_face+face]/num_face;
            if(face == n_face && neigh == b)
            {
              //boundary
              BFAM_VERBOSE("BOUNDARY: %s",glue_name);
            }
            else
            {
              if(dim == 2)
              {
                int m_mask = 0;
                switch(face)
                {
                  case 0:
                  case 1:
                    m_mask = 1;
                    break;
                  case 2:
                  case 3:
                    m_mask = 0;
                    break;
                }
                int n_mask = 0;
                switch(n_face)
                {
                  case 0:
                  case 1:
                    n_mask = 1;
                    break;
                  case 2:
                  case 3:
                    n_mask = 0;
                    break;
                }
                int n_pd[2] = {0,0};
                bfam_locidx_t n_size = procs[2*neigh+1]-procs[2*neigh]+1;
                BFAM_MPI_CHECK(MPI_Dims_create(n_size, 2, n_pd));

                bfam_gloidx_t n_N = N[dim*neigh + n_mask];
                BFAM_ABORT_IF_NOT(N[dim*b+m_mask] == n_N,
                    "Cannot connect %3d face %3d of size %d with "
                    "%3d face %3d of size %d. "
                    "Non-conforming not implemented.",b,face,N[dim*b+m_mask],
                    neigh,n_face,n_N);

                /* m_ix are the smallest and largest indices I need from my
                 * neighbors */
                bfam_gloidx_t m_ix[2] = {gx[m_mask],gx[m_mask]+Nl[m_mask]};

                int n_pd_ix[2] = {0,0};
                switch(n_face)
                {
                  case 1:
                    n_pd_ix[0] = n_pd[0]-1;
                    break;
                  case 3:
                    n_pd_ix[1] = n_pd[1]-1;
                    break;
                }

                for(int nr = 0; nr < n_pd[n_mask];nr++)
                {
                  bfam_locidx_t n_Nl = -1;
                  bfam_gloidx_t n_gx = -1;
                  bfam_locidx_t n_Nb[2] = {-1,-1};

                  /* determine neighbors orientation and flip if needed */
                  simple_partition_1d(&n_Nl,&n_gx,n_Nb,n_N,n_pd[n_mask],nr,0);
                  if(orient==1)
                  {
                    n_gx = n_N - (n_gx+n_Nl);
                    bfam_locidx_t tmp = n_Nb[1];
                    n_Nb[1] = n_Nb[0];
                    n_Nb[0] = tmp;
                  }

                  /* determine overlap */
                  bfam_gloidx_t n_ix[2] = {n_gx,n_gx+n_Nl};
                  bfam_gloidx_t glue_ix[2] =
                    {BFAM_MAX(n_ix[0],m_ix[0]),BFAM_MIN(n_ix[1],m_ix[1])};
                  bfam_gloidx_t glue_size = glue_ix[1]+1-glue_ix[0];

                  if(glue_size > 0)
                  {
                    n_pd_ix[n_mask] = nr;
                    bfam_locidx_t n_rank = procs[2*neigh]
                                         + (n_pd_ix[0]+n_pd_ix[1]*n_pd[0]);
                    /* determine overlap */
                    char tmp_name[BFAM_BUFSIZ];
                    snprintf(tmp_name,BFAM_BUFSIZ,"%s_piece_%d",glue_name,nr);
                    bfam_subdomain_sbp_inter_glue_t* glue =
                      bfam_subdomain_sbp_inter_glue_new(glue_id,tmp_name,
                          rank,n_rank,sub,glue_ix,face,orient,
                          neigh,nr,n_face);
                    bfam_subdomain_add_tag((bfam_subdomain_t*)glue,
                                           "_inter_glue");
                    bfam_domain_add_subdomain(domain,(bfam_subdomain_t*)glue);
                  }
                }
              }
              else
              {
                int m_mask[2] = {0,0};
                switch(face)
                {
                  case 0:
                  case 1:
                    m_mask[0] = 1;
                    m_mask[1] = 2;
                    break;
                  case 2:
                  case 3:
                    m_mask[0] = 0;
                    m_mask[1] = 2;
                    break;
                  case 4:
                  case 5:
                    m_mask[0] = 0;
                    m_mask[1] = 1;
                    break;
                }
                int n_mask[2] = {0,0};
                switch(face)
                {
                  case 0:
                  case 1:
                    n_mask[0] = 1;
                    n_mask[1] = 2;
                    break;
                  case 2:
                  case 3:
                    n_mask[0] = 0;
                    n_mask[1] = 2;
                    break;
                  case 4:
                  case 5:
                    n_mask[0] = 0;
                    n_mask[1] = 1;
                    break;
                }
                bfam_gloidx_t n_N[2] =
                   {N[dim*neigh + n_mask[0]],N[dim*neigh + n_mask[1]]};
                bfam_gloidx_t m_N[2] =
                   {N[dim*b + m_mask[0]],N[dim*b + m_mask[1]]};

                int n_pd[3] = {0,0,0};
                bfam_locidx_t n_size = procs[2*neigh+1]-procs[2*neigh]+1;
                BFAM_MPI_CHECK(MPI_Dims_create(n_size, 3, n_pd));
                orient = BFAM_P8EST_ORIENTATION(face,n_face,orient);

                switch(orient)
                {
                  case 0:
                    BFAM_ABORT_IF_NOT(m_N[0] == n_N[0] && m_N[1] == n_N[1],
                        "Cannot connect %3d face %3d of size (%d,%d) with "
                        "%3d face %3d of size (%d,%3d) with orientation %d. "
                        "Non-conforming not implemented.",
                        b,face,m_N[0],m_N[1],
                        neigh,n_face,n_N[0],n_N[1],
                        orient);
                    break;
                  case 1:
                    BFAM_ABORT_IF_NOT(m_N[1] == n_N[0] && m_N[0] == n_N[1],
                        "Cannot connect %3d face %3d of size (%d,%d) with "
                        "%3d face %3d of size (%d,%3d) with orientation %d. "
                        "Non-conforming not implemented.",
                        b,face,m_N[0],m_N[1],
                        neigh,n_face,n_N[0],n_N[1],orient);
                    break;
                  case 3:
                    BFAM_ABORT_IF_NOT(m_N[1] == n_N[0] && m_N[0] == n_N[1],
                        "Cannot connect %3d face %3d of size (%d,%d) with "
                        "%3d face %3d of size (%d,%3d) with orientation %d. "
                        "Non-conforming not implemented.",
                        b,face,m_N[0],m_N[1],
                        neigh,n_face,n_N[0],n_N[1],orient);
                    break;
                  case 4:
                    BFAM_ABORT_IF_NOT(m_N[1] == n_N[0] && m_N[0] == n_N[1],
                        "Cannot connect %3d face %3d of size (%d,%d) with "
                        "%3d face %3d of size (%d,%3d) with orientation %d. "
                        "Non-conforming not implemented.",
                        b,face,m_N[0],m_N[1],
                        neigh,n_face,n_N[0],n_N[1],orient);
                    break;
                  default:
                    BFAM_ABORT("orientation %d: not implemented yet",orient);
                    break;
                }
              }
            }
          }
        }
      }
    }
  }

  if(x != NULL) bfam_free_aligned(x);
  if(y != NULL) bfam_free_aligned(y);
  if(z != NULL) bfam_free_aligned(z);
}

void
test_2d(int rank, int mpi_size,MPI_Comm mpicomm)
{
  /* setup the domain*/
  bfam_domain_t domain;
  bfam_domain_init(&domain,MPI_COMM_WORLD);

  /* connectivity */
  int num_blocks = 3;
  const char *names[] = {"sub0","sub1","sub2"};
  bfam_locidx_t bufsz = 3;

  int dim   =  2;
  bfam_locidx_t   EToV[] = {0,1,3,4,
                            1,2,4,5,
                            4,5,3,6};
  bfam_locidx_t   EToE[] = {0,1,0,2,
                            0,1,1,2,
                            0,2,1,2};
  int8_t          EToF[] = {0,0,2,0+4,
                            1,1,2,2,
                            3+4,1,3,3};
  bfam_gloidx_t      N[] = {100,110,
                            120,110,
                            120,100};
  bfam_long_real_t  Vx[] = {0,0.25,1,  0,0.25,0.5,0};
  bfam_long_real_t  Vy[] = {0,   0,0,0.5,0.25,0.5,1};
  bfam_long_real_t *Vz   = NULL;
  bfam_locidx_t  procs[2*num_blocks];

  /* load balance */
  simple_load_balance(procs,N,mpi_size,num_blocks,dim);

  /* setup subdomains */
  setup_subdomains(&domain,dim,num_blocks,bufsz,rank,
      names,procs,EToV,EToE,EToF,N,Vx,Vy,Vz);

  /* put in some fields */
  const char *volume[] = {"_volume", NULL};
  bfam_domain_add_field(&domain, BFAM_DOMAIN_OR, volume, "p1");
  bfam_domain_add_field(&domain, BFAM_DOMAIN_OR, volume, "p2");

  bfam_domain_init_field(&domain, BFAM_DOMAIN_OR, volume, "p1",
      0, poly1_field, NULL);
  bfam_domain_init_field(&domain, BFAM_DOMAIN_OR, volume, "p2",
      0, poly2_field, NULL);

  const char *intra_glue[] = {"_intra_glue", NULL};
  bfam_domain_add_minus_field(&domain, BFAM_DOMAIN_OR, intra_glue, "p1");
  bfam_domain_add_minus_field(&domain, BFAM_DOMAIN_OR, intra_glue, "p2");

  bfam_domain_add_plus_field(&domain, BFAM_DOMAIN_OR, intra_glue, "p1");
  bfam_domain_add_plus_field(&domain, BFAM_DOMAIN_OR, intra_glue, "p2");

  const char *inter_glue[] = {"_inter_glue", NULL};
  const char *fields[] = {"p1","p2","_grid_x","_grid_y",NULL};
  for(int f = 0;fields[f];f++)
  {
    bfam_domain_add_minus_field(&domain, BFAM_DOMAIN_OR, inter_glue, fields[f]);
    bfam_domain_add_plus_field(&domain, BFAM_DOMAIN_OR, inter_glue, fields[f]);
  }

  const char *glue[] = {"_inter_glue","_intra_glue", NULL};

  bfam_communicator_t* communicator =
    bfam_communicator_new(&domain, BFAM_DOMAIN_OR, glue,
        mpicomm, 10);

  /* start recv_send */
  bfam_communicator_start(communicator);

  /* finish recv */
  bfam_communicator_finish(communicator);

  /* check values */
  bfam_domain_init_field(&domain, BFAM_DOMAIN_OR, volume, "p1",
      0, poly1_field_check, NULL);
  bfam_domain_init_field(&domain, BFAM_DOMAIN_OR, volume, "p2",
      0, poly2_field_check, NULL);


  bfam_subdomain_sbp_inter_glue_t *subs[domain.numSubdomains];
  bfam_locidx_t numSubs;
  bfam_domain_get_subdomains(&domain,BFAM_DOMAIN_OR,inter_glue,
      domain.numSubdomains,(bfam_subdomain_t**) subs,&numSubs);
  for(int s = 0; s < numSubs;s++)
  {
    for(int f = 0;fields[f];f++)
    {
      bfam_real_t * f_m =
        bfam_dictionary_get_value_ptr(&subs[s]->base.fields_m, fields[f]);
      bfam_real_t * f_p =
        bfam_dictionary_get_value_ptr(&subs[s]->base.fields_p, fields[f]);
      for(int i = 0; i < (int)subs[s]->field_size_m;i++)
        /* for some reason valgrind does something to values that they are not
         * exactly the same, would prefer f_p[i] == f_m[i] */
        /* BFAM_ABORT_IF_NOT(f_p[i] == f_m[i],*/
        BFAM_ABORT_IF_NOT(BFAM_REAL_APPROX_EQ(f_p[i],f_m[i],10),
            "inter check: sub: %d face: %d field %s"
            " failed on %d got %f expected %f",
            subs[s]->sub_m->base.id,subs[s]->face,
            fields[f],i,f_m[i],f_p[i]);
    }
  }

  /* dump the entire mesh */
  const char *tags[] = {"_volume",NULL};
  const char *scalars[] = {"p1","p2",NULL};
  const char *vectors[] = {"p",NULL};
  const char *components[] = {"p1","p2","p3",NULL};
  bfam_vtk_write_struc_file(&domain,BFAM_DOMAIN_AND,
      tags,"sbp_fields_2d",scalars,vectors,components,0,1);

  /* clean up */
  bfam_communicator_free(communicator);
  bfam_free(communicator);
  bfam_domain_free(&domain);
}

void
test_3d(int rank, int mpi_size, MPI_Comm mpicomm)
{
  /* setup the domain*/
  bfam_domain_t domain;
  bfam_domain_init(&domain,MPI_COMM_WORLD);

  int dim   = 3;

  const char *names[] = {"sub0","sub1","sub2","sub3","sub4","sub5"};
  bfam_locidx_t bufsz = 3;

  /* connectivity from p8est_connectivity_new_rotcubes */
  bfam_locidx_t num_blocks = 6;
  const bfam_locidx_t EToV[6 * 8] = {
    0, 17, 3, 4, 15, 11, 13, 14,
    7, 2, 6, 17, 9, 12, 8, 11,
    2, 12, 5, 10, 17, 11, 4, 14,
    19, 13, 18, 14, 16, 15, 1, 11,
    14, 11, 21, 25, 18, 1, 22, 23,
    21, 20, 25, 24, 14, 10, 11, 12,
  };
  const bfam_gloidx_t N[6 * 3] = {
    11,11,11,
    11,11,11,
    11,11,11,
    11,11,11,
    11,11,11,
    11,11,11,
  };
  const bfam_locidx_t EToE[6 * 6] = {
    0, 2, 0, 0, 0, 3,
    1, 2, 1, 1, 1, 1,
    2, 5, 1, 2, 2, 0,
    3, 0, 3, 4, 3, 3,
    4, 4, 3, 4, 5, 4,
    4, 5, 5, 5, 5, 2,
  };
  const int8_t EToF[6 * 6] = {
    0, 5, 2, 3, 4, 13,
    0, 2, 2, 3, 4, 5,
    0, 23, 1, 3, 4, 1,
    0, 17, 2, 8, 4, 5,
    0, 1, 9, 3, 12, 5,
    16, 1, 2, 3, 4, 19,
  };

  const bfam_long_real_t Vx[26] =
    {0,1,2,0,1,2,1,2,1,2,2,1,2,0,1,0,0,1,1,0,2.5,2,2,2,2.5,2};
  const bfam_long_real_t Vy[26] =
    {0,0,0,1,1,1,-1,-1,-1,-1,1,0,0,1,1,0,0,0,1,1,1.5,1.5,1.5,.5,.5,.5};
  const bfam_long_real_t Vz[26] =
    {0,2,0,0,0,0,0,0,1,1,1,1,1,1,1,1,2,0,2,2,2,2,2.5,2.5,2,2};

  bfam_locidx_t  procs[2*num_blocks];

  /* load balance */
  simple_load_balance(procs,N,mpi_size,num_blocks,dim);

  /* setup subdomains */
  setup_subdomains(&domain,dim,num_blocks,bufsz,rank,
      names,procs,EToV,EToE,EToF,N,Vx,Vy,Vz);

  /* put in some fields */
  const char *volume[] = {"_volume", NULL};
  bfam_domain_add_field(&domain, BFAM_DOMAIN_OR, volume, "p1");
  bfam_domain_add_field(&domain, BFAM_DOMAIN_OR, volume, "p2");
  bfam_domain_add_field(&domain, BFAM_DOMAIN_OR, volume, "p3");

  bfam_domain_init_field(&domain, BFAM_DOMAIN_OR, volume, "p1",
      0, poly1_field, NULL);
  bfam_domain_init_field(&domain, BFAM_DOMAIN_OR, volume, "p2",
      0, poly2_field, NULL);
  bfam_domain_init_field(&domain, BFAM_DOMAIN_OR, volume, "p3",
      0, poly3_field, NULL);


  const char *glue[] = {"_intra_glue", NULL};
  bfam_domain_add_minus_field(&domain, BFAM_DOMAIN_OR, glue, "p1");
  bfam_domain_add_minus_field(&domain, BFAM_DOMAIN_OR, glue, "p2");
  bfam_domain_add_minus_field(&domain, BFAM_DOMAIN_OR, glue, "p3");

  bfam_domain_add_plus_field(&domain, BFAM_DOMAIN_OR, glue, "p1");
  bfam_domain_add_plus_field(&domain, BFAM_DOMAIN_OR, glue, "p2");
  bfam_domain_add_plus_field(&domain, BFAM_DOMAIN_OR, glue, "p3");

  bfam_communicator_t* communicator =
    bfam_communicator_new(&domain, BFAM_DOMAIN_OR, glue,
        mpicomm, 10);

  /* start recv_send */
  bfam_communicator_start(communicator);

  /* finish recv */
  bfam_communicator_finish(communicator);


  /* check the comm results */
  bfam_domain_init_field(&domain, BFAM_DOMAIN_OR, volume, "p1",
      0, poly1_field_check, NULL);
  bfam_domain_init_field(&domain, BFAM_DOMAIN_OR, volume, "p2",
      0, poly2_field_check, NULL);
  bfam_domain_init_field(&domain, BFAM_DOMAIN_OR, volume, "p3",
      0, poly3_field_check, NULL);

  /* dump the entire mesh */
  // const char *tags[] = {"_volume",NULL};
  // const char *scalars[] = {"p1","p2","p3",NULL};
  // const char *vectors[] = {"p",NULL};
  // const char *components[] = {"p1","p2","p3",NULL};
  // bfam_vtk_write_struc_file(&domain,BFAM_DOMAIN_AND,
  //     tags,"sbp_fields_3d",scalars,vectors,components,1,1);


  /* clean up */
  bfam_communicator_free(communicator);
  bfam_free(communicator);
  bfam_domain_free(&domain);
}

void
test_3d_2block(int rank, int mpi_size, MPI_Comm mpicomm,
    int face0,int face1,int rot)
{
  /* setup the domain*/
  bfam_domain_t domain;
  bfam_domain_init(&domain,MPI_COMM_WORLD);

  int dim   = 3;

  char *names[2];
  char sub0[BFAM_BUFSIZ];
  char sub1[BFAM_BUFSIZ];
  names[0] = sub0;
  names[1] = sub1;

  snprintf(names[0],BFAM_BUFSIZ,"sub0_face%d_rot%d",face0,rot);
  snprintf(names[1],BFAM_BUFSIZ,"sub1_face%d_rot%d",face1,rot);

  bfam_locidx_t bufsz = 3;

  bfam_locidx_t num_blocks = 2;
  const bfam_long_real_t Vx[12] = {
    0,1,2,
    0,1,2,
    0,1,2,
    0,1,2
  };
  const bfam_long_real_t Vy[12] = {
    0,0,0,
    1,1,1,
    0,0,0,
    1,1,1
  };
  const bfam_long_real_t Vz[12] = {
    0,0,0,
    0,0,0,
    1,1,1,
    1,1,1
  };

  bfam_locidx_t EToV[2 * 8];

  switch(face0)
  {
    case 0:
      EToV[0] =  4;
      EToV[1] =  3;
      EToV[2] =  1;
      EToV[3] =  0;
      EToV[4] = 10;
      EToV[5] =  9;
      EToV[6] =  7;
      EToV[7] =  6;
      break;
    case 1:
      EToV[0] =  0;
      EToV[1] =  1;
      EToV[2] =  3;
      EToV[3] =  4;
      EToV[4] =  6;
      EToV[5] =  7;
      EToV[6] =  9;
      EToV[7] = 10;
      break;
    case 2:
      EToV[0] =  1;
      EToV[1] =  4;
      EToV[2] =  0;
      EToV[3] =  3;
      EToV[4] =  7;
      EToV[5] = 10;
      EToV[6] =  6;
      EToV[7] =  9;
      break;
    case 3:
      EToV[0] =  3;
      EToV[1] =  0;
      EToV[2] =  4;
      EToV[3] =  1;
      EToV[4] =  9;
      EToV[5] =  6;
      EToV[6] = 10;
      EToV[7] =  7;
      break;
    case 4:
      EToV[0] =  1;
      EToV[1] =  7;
      EToV[2] =  4;
      EToV[3] = 10;
      EToV[4] =  0;
      EToV[5] =  6;
      EToV[6] =  3;
      EToV[7] =  9;
      break;
    case 5:
      EToV[0] =  3;
      EToV[1] =  9;
      EToV[2] =  0;
      EToV[3] =  6;
      EToV[4] =  4;
      EToV[5] = 10;
      EToV[6] =  1;
      EToV[7] =  7;
      break;
  }

  switch(face1)
  {
    case 1:
      EToV[8+0] =  4+1;
      EToV[8+1] =  3+1;
      EToV[8+2] =  1+1;
      EToV[8+3] =  0+1;
      EToV[8+4] = 10+1;
      EToV[8+5] =  9+1;
      EToV[8+6] =  7+1;
      EToV[8+7] =  6+1;
      break;
    case 0:
      EToV[8+0] =  0+1;
      EToV[8+1] =  1+1;
      EToV[8+2] =  3+1;
      EToV[8+3] =  4+1;
      EToV[8+4] =  6+1;
      EToV[8+5] =  7+1;
      EToV[8+6] =  9+1;
      EToV[8+7] = 10+1;
      break;
    case 3:
      EToV[8+0] =  1+1;
      EToV[8+1] =  4+1;
      EToV[8+2] =  0+1;
      EToV[8+3] =  3+1;
      EToV[8+4] =  7+1;
      EToV[8+5] = 10+1;
      EToV[8+6] =  6+1;
      EToV[8+7] =  9+1;
      break;
    case 2:
      EToV[8+0] =  3+1;
      EToV[8+1] =  0+1;
      EToV[8+2] =  4+1;
      EToV[8+3] =  1+1;
      EToV[8+4] =  9+1;
      EToV[8+5] =  6+1;
      EToV[8+6] = 10+1;
      EToV[8+7] =  7+1;
      break;
    case 5:
      EToV[8+0] =  1+1;
      EToV[8+1] =  7+1;
      EToV[8+2] =  4+1;
      EToV[8+3] = 10+1;
      EToV[8+4] =  0+1;
      EToV[8+5] =  6+1;
      EToV[8+6] =  3+1;
      EToV[8+7] =  9+1;
      break;
    case 4:
      EToV[8+0] =  3+1;
      EToV[8+1] =  9+1;
      EToV[8+2] =  0+1;
      EToV[8+3] =  6+1;
      EToV[8+4] =  4+1;
      EToV[8+5] = 10+1;
      EToV[8+6] =  1+1;
      EToV[8+7] =  7+1;
      break;
  }

  bfam_locidx_t EToE[2 * 6] = {
    0,0,0,0,0,0,
    1,1,1,1,1,1,
  };

  int8_t EToF[2 * 6] = {
    0,1,2,3,4,5,
    0,1,2,3,4,5,
  };

  bfam_gloidx_t N[2 * 3] = {
    50,75,100,
    125,150,175,
  };


  bfam_locidx_t  procs[2*num_blocks];

  /* load balance */
  simple_load_balance(procs,N,mpi_size,num_blocks,dim);

  /* setup subdomains */
  setup_subdomains(&domain,dim,num_blocks,bufsz,rank,
      (const char**)names,procs,EToV,EToE,EToF,N,Vx,Vy,Vz);

  //JK /* put in some fields */
  //JK const char *volume[] = {"_volume", NULL};
  //JK bfam_domain_add_field(&domain, BFAM_DOMAIN_OR, volume, "p1");
  //JK bfam_domain_add_field(&domain, BFAM_DOMAIN_OR, volume, "p2");
  //JK bfam_domain_add_field(&domain, BFAM_DOMAIN_OR, volume, "p3");

  //JK bfam_domain_init_field(&domain, BFAM_DOMAIN_OR, volume, "p1",
  //JK     0, poly1_field, NULL);
  //JK bfam_domain_init_field(&domain, BFAM_DOMAIN_OR, volume, "p2",
  //JK     0, poly2_field, NULL);
  //JK bfam_domain_init_field(&domain, BFAM_DOMAIN_OR, volume, "p3",
  //JK     0, poly3_field, NULL);


  //JK const char *glue[] = {"_intra_glue", NULL};
  //JK bfam_domain_add_minus_field(&domain, BFAM_DOMAIN_OR, glue, "p1");
  //JK bfam_domain_add_minus_field(&domain, BFAM_DOMAIN_OR, glue, "p2");
  //JK bfam_domain_add_minus_field(&domain, BFAM_DOMAIN_OR, glue, "p3");

  //JK bfam_domain_add_plus_field(&domain, BFAM_DOMAIN_OR, glue, "p1");
  //JK bfam_domain_add_plus_field(&domain, BFAM_DOMAIN_OR, glue, "p2");
  //JK bfam_domain_add_plus_field(&domain, BFAM_DOMAIN_OR, glue, "p3");

  //JK bfam_communicator_t* communicator =
  //JK   bfam_communicator_new(&domain, BFAM_DOMAIN_OR, glue,
  //JK       mpicomm, 10);

  //JK /* start recv_send */
  //JK bfam_communicator_start(communicator);

  //JK /* finish recv */
  //JK bfam_communicator_finish(communicator);


  //JK /* check the comm results */
  //JK bfam_domain_init_field(&domain, BFAM_DOMAIN_OR, volume, "p1",
  //JK     0, poly1_field_check, NULL);
  //JK bfam_domain_init_field(&domain, BFAM_DOMAIN_OR, volume, "p2",
  //JK     0, poly2_field_check, NULL);
  //JK bfam_domain_init_field(&domain, BFAM_DOMAIN_OR, volume, "p3",
  //JK     0, poly3_field_check, NULL);


  //JK /* clean up */
  //JK bfam_communicator_free(communicator);
  //JK bfam_free(communicator);
  //JK bfam_domain_free(&domain);
}

int
main (int argc, char *argv[])
{
  int rank, mpi_size;

  /* startup MPI */
  BFAM_MPI_CHECK(MPI_Init(&argc,&argv));
  BFAM_MPI_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));
  BFAM_MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
  bfam_log_init(rank,stdout,BFAM_LL_VERBOSE);

  /* test 2d */
  test_2d(rank,mpi_size,MPI_COMM_WORLD);

  /* test 3d */
  test_3d(rank,mpi_size,MPI_COMM_WORLD);

  for(int face0 = 0; face0 < 8; face0++)
    for(int face1 = face0; face1 < 8; face1++)
      for(int rot = 0; rot < 4; rot++)
        test_3d_2block(rank,mpi_size,MPI_COMM_WORLD,face0,face1,rot);

  /* stop MPI */
  BFAM_MPI_CHECK(MPI_Finalize());
}
