#ifdef BLADE_DGX_DIMENSION

/* START DIM 2 */
#if    BLADE_DGX_DIMENSION==2
#include "blade_dgx_rhs_2.h"

#define DIM 2
#define BLADE_D3_AP(A1,A2) (A1)
#define BLADE_D3_OP(A) BFAM_NOOP()


#define BLADE_DR1(A,B)     BFAM_KRON_IXA    (N+1, Dr , A, B)
#define BLADE_DR1_PE(A,B)  BFAM_KRON_IXA_PE (N+1, Dr , A, B)
#define BLADE_DR1T(A,B)    BFAM_KRON_IXAT   (N+1, Dr , A, B)
#define BLADE_DR1T_PE(A,B) BFAM_KRON_IXAT_PE(N+1, Dr , A, B)

#define BLADE_DR2(A,B)     BFAM_KRON_AXI    (N+1, Dr , A, B)
#define BLADE_DR2_PE(A,B)  BFAM_KRON_AXI_PE (N+1, Dr , A, B)
#define BLADE_DR2T(A,B)    BFAM_KRON_ATXI   (N+1, Dr , A, B)
#define BLADE_DR2T_PE(A,B) BFAM_KRON_ATXI_PE(N+1, Dr , A, B)

#define BLADE_DR3(A,B)     BFAM_NOOP()
#define BLADE_DR3_PE(A,B)  BFAM_NOOP()
#define BLADE_DR3T(A,B)    BFAM_NOOP()
#define BLADE_DR3T_PE(A,B) BFAM_NOOP()

/* START DIM 2 Generic */
#ifndef NORDER
#define NORDER
#define USE_GENERIC
#define GENERIC_INIT(INN,NFM)\
  const int N   = INN; \
  const int Np  = (INN+1)*(INN+1); \
  const int Nfp = (INN+1); \
  const int Nfaces  = 4; \
  BFAM_WARNING(\
      "using generic NFM for 2d of order %d with "\
      "Np = %d, Nfp = %d, and Nfaces = %d",\
      N, Np, Nfp, Nfaces);
#else
#define GENERIC_INIT(INN,NFM) BFAM_NOOP()
#define N    NORDER
#define Np  (NORDER+1)*(NORDER+1)
#define Np_BACK  (NORDER+1)*(NORDER+1)
#define Nfp (NORDER+1)
#define Nfaces  4
#endif
/* END DIM 2 Generic */

/* END DIM 2 */
/* START DIM 3 */
#elif  BLADE_DGX_DIMENSION==3
#include "blade_dgx_rhs_3.h"

#define DIM 3
#define BLADE_D3_AP(A1,A2) (A1 A2)
#define BLADE_D3_OP(A) A

#define BLADE_DR1(A,B)     BFAM_KRON_IXIXA    (N+1, Dr , A, B)
#define BLADE_DR1_PE(A,B)  BFAM_KRON_IXIXA_PE (N+1, Dr , A, B)
#define BLADE_DR1T(A,B)    BFAM_KRON_IXIXAT   (N+1, Dr , A, B)
#define BLADE_DR1T_PE(A,B) BFAM_KRON_IXIXAT_PE(N+1, Dr , A, B)

#define BLADE_DR2(A,B)     BFAM_KRON_IXAXI    (N+1, Dr , A, B)
#define BLADE_DR2_PE(A,B)  BFAM_KRON_IXAXI_PE (N+1, Dr , A, B)
#define BLADE_DR2T(A,B)    BFAM_KRON_IXATXI   (N+1, Dr , A, B)
#define BLADE_DR2T_PE(A,B) BFAM_KRON_IXATXI_PE(N+1, Dr , A, B)

#define BLADE_DR3(A,B)     BFAM_KRON_AXIXI   (N+1, Dr , A, B)
#define BLADE_DR3_PE(A,B)  BFAM_KRON_AXIXI_PE(N+1, Dr , A, B)
#define BLADE_DR3T(A,B)    BFAM_KRON_ATXIXI   (N+1, Dr , A, B)
#define BLADE_DR3T_PE(A,B) BFAM_KRON_ATXIXI_PE(N+1, Dr , A, B)

#ifndef NORDER
#define NORDER
#define USE_GENERIC
#define GENERIC_INIT(INN,NFM)\
  const int N   = INN; \
  const int Np  = (INN+1)*(INN+1)*(INN+1); \
  const int Nfp = (INN+1)*(INN+1); \
  const int Nfaces  = 6; \
  BFAM_WARNING(\
      "using generic NFM for 3d of order %d with "\
      "Np = %d, Nfp = %d, and Nfaces = %d",\
      N, Np, Nfp, Nfaces);
#else
#define GENERIC_INIT(INN,NFM) BFAM_NOOP()
#define N    NORDER
#define Np  (NORDER+1)*(NORDER+1)*(NORDER+1)
#define Np_BACK  (NORDER+1)*(NORDER+1)*(NORDER+1)
#define Nfp (NORDER+1)*(NORDER+1)
#define Nfaces  6
#endif

/* END DIM 3 */
/* START BAD DIMENSION */
#else
#error "bad dimension"
#endif
/* END BAD DIMENSION */

/* Some useful field macros */
#define BFAM_LOAD_FIELD_RESTRICT_ALIGNED(field,prefix,base,dictionary)         \
bfam_real_t *restrict field;                                                   \
{                                                                              \
  char bfam_load_field_name[BFAM_BUFSIZ];                                      \
  snprintf(bfam_load_field_name,BFAM_BUFSIZ,"%s%s",(prefix),(base));           \
  field = bfam_dictionary_get_value_ptr(dictionary, bfam_load_field_name);     \
  BFAM_ASSERT(field != NULL);                                                  \
}                                                                              \
BFAM_ASSUME_ALIGNED(field,32);

#define BFAM_LOAD_FIELD_ALIGNED(field,prefix,base,dictionary)                  \
bfam_real_t *field;                                                            \
{                                                                              \
  char bfam_load_field_name[BFAM_BUFSIZ];                                      \
  snprintf(bfam_load_field_name,BFAM_BUFSIZ,"%s%s",(prefix),(base));           \
  field = bfam_dictionary_get_value_ptr(dictionary, bfam_load_field_name);     \
  BFAM_ASSERT(field != NULL);                                                  \
}                                                                              \
BFAM_ASSUME_ALIGNED(field,32);


/* Function macros */
#define BLADE_APPEND_4(a,b,c,d) a ## b ## c ##d
#define BLADE_APPEND_EXPAND_4(a,b,c,d) BLADE_APPEND_4(a,b,c,d)

#define blade_dgx_scale_rates_advection \
  BLADE_APPEND_EXPAND_4(blade_dgx_scale_rates_advection_,DIM,_,NORDER)

#define blade_dgx_intra_rhs_advection \
  BLADE_APPEND_EXPAND_4(blade_dgx_intra_rhs_advection_,DIM,_,NORDER)

#define blade_dgx_inter_rhs_interface \
  BLADE_APPEND_EXPAND_4(blade_dgx_inter_rhs_interface_,DIM,_,NORDER)

#define blade_dgx_energy \
  BLADE_APPEND_EXPAND_4(blade_dgx_energy_,DIM,_,NORDER)

#define blade_dgx_add_rates_advection \
  BLADE_APPEND_EXPAND_4(blade_dgx_add_rates_advection_,DIM,_,NORDER)

#define HALF BFAM_REAL(0.5)

#define MASSPROJECTION   (1)
#define PROJECTION       (0)

void blade_dgx_scale_rates_advection(
    int inN, bfam_subdomain_dgx_t *sub, const char *rate_prefix,
    const bfam_long_real_t a)
{
  GENERIC_INIT(inN,blade_dgx_scale_rates_advection);

  const bfam_locidx_t num_pts = sub->K * Np;
  bfam_dictionary_t *fields = &sub->base.fields;

  const char *f_names[] = {"q",NULL};

  for(bfam_locidx_t f = 0; f_names[f]!=NULL; ++f)
  {
    BFAM_LOAD_FIELD_RESTRICT_ALIGNED(rate,rate_prefix,f_names[f],fields);
    for(bfam_locidx_t n = 0; n < num_pts; n++) rate[n] *= a;
  }
}

void blade_dgx_energy(
    int inN, bfam_real_t *energy_sq,
    bfam_subdomain_dgx_t *sub, const char *field_prefix)
{
  GENERIC_INIT(inN,blade_dgx_energy);

  /* get the fields we will need */
  bfam_dictionary_t *fields = &sub->base.fields;
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(q ,field_prefix,"q" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(J  ,"","_grid_J",fields);

  bfam_locidx_t K  = sub->K;

  bfam_real_t *w  = sub->w;
  BFAM_ASSUME_ALIGNED(w ,32);

  /* loop through all the elements */
  for(bfam_locidx_t e = 0; e < K;e++)
  {
    bfam_locidx_t n = e*Np;
#if DIM==3
    for(bfam_locidx_t k = 0; k < N+1; k++)
    {
      const bfam_real_t wk = w[k];
#endif
      for(bfam_locidx_t j = 0; j < N+1; j++)
      {
        const bfam_real_t wj = w[j];
        for(bfam_locidx_t i = 0; i < N+1; i++,n++)
        {
          const bfam_real_t wi = w[i];

          energy_sq[0] += BLADE_D3_AP(wi*wj,*wk)*J[n]*q[n]*q[n];
        }
      }
#if DIM==3
    }
#endif
  }
}

#define GAM (BFAM_REAL(0))
static inline void
blade_dgx_compute_flux(
    const bfam_real_t   u, const bfam_real_t   qM, const bfam_real_t qP,
          bfam_real_t  *uqS)
{
  /* compute the state u*q */
  *uqS  = HALF*u*qM; /* half due to the skew-symmetric splitting */
  *uqS -= HALF*u*(qM+qP);
}

static inline void
blade_dgx_add_flux(const bfam_real_t scale,
          bfam_real_t *dq, const bfam_locidx_t iM, const bfam_real_t uqS,
    const bfam_real_t  sJ, const bfam_real_t   JI, const bfam_real_t wi)
{
  dq[iM] += scale*wi*JI*sJ*uqS;
}

void blade_dgx_intra_rhs_advection(
    int inN, bfam_subdomain_dgx_t *sub, const char *rate_prefix,
    const char *field_prefix, const bfam_long_real_t t)
{
  GENERIC_INIT(inN,blade_dgx_intra_rhs_advection);

  /* get the fields we will need */
  bfam_dictionary_t *fields = &sub->base.fields;
  bfam_dictionary_t *fields_face = &sub->base.fields_face;

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(q  ,field_prefix,"q" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dq ,rate_prefix ,"q" ,fields);

  /* get the material properties and metric terms */
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(ux ,field_prefix,"ux" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(uy ,field_prefix,"uy" ,fields);
  BLADE_D3_OP(BFAM_LOAD_FIELD_RESTRICT_ALIGNED(uz,field_prefix,"uz",fields));

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(J    ,"","_grid_J"    ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(JI   ,"","_grid_JI"   ,fields);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Jr1x1,"","_grid_Jr0x0",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Jr1x2,"","_grid_Jr0x1",fields);
  BLADE_D3_OP(BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Jr1x3,"","_grid_Jr0x2",fields));

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Jr2x1,"","_grid_Jr1x0",fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Jr2x2,"","_grid_Jr1x1",fields);
  BLADE_D3_OP(BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Jr2x3,"","_grid_Jr1x2",fields));

  BLADE_D3_OP(BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Jr3x1,"","_grid_Jr2x0",fields));
  BLADE_D3_OP(BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Jr3x2,"","_grid_Jr2x1",fields));
  BLADE_D3_OP(BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Jr3x3,"","_grid_Jr2x2",fields));

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n1,"","_grid_nx0",fields_face);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n2,"","_grid_nx1",fields_face);
  BLADE_D3_OP(BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n3,"","_grid_nx2",fields_face));
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(sJ  ,"","_grid_sJ",fields_face);

  bfam_locidx_t K  = sub->K;

  BFAM_ALIGN(32) bfam_real_t aux0[Np];
  BFAM_ALIGN(32) bfam_real_t aux1[Np];
  BFAM_ALIGN(32) bfam_real_t aux2[Np];
  BLADE_D3_OP(BFAM_ALIGN(32) bfam_real_t aux3[Np]);

  bfam_real_t *Dr = sub->Dr;
  bfam_real_t *w  = sub->w;
  bfam_real_t *wi  = sub->wi;
  BFAM_ASSUME_ALIGNED(Dr,32);
  BFAM_ASSUME_ALIGNED(w ,32);
  BFAM_ASSUME_ALIGNED(wi ,32);

  bfam_locidx_t *vmapM = sub->vmapM;
  bfam_locidx_t *vmapP = sub->vmapP;


  /* loop through all the elements */
  for(bfam_locidx_t e = 0; e < K;e++)
  {
    bfam_locidx_t off = e*Np;

    /* do differential equation */
    /* Note: Matrices on the LHS will be handled in add rates routine */

    /****** X-DERIVATIVE ******/
    /* q -= JI*(Jrx*Dr*(ux*q) + Jsx*Ds*(ux*q) + Jtx*Dt*(ux*q)) */
    for(bfam_locidx_t i = 0; i < Np; i++) aux0[i] = ux[off+i] * q[off+i];
    BLADE_DR1(aux0, aux1); /* Dr */
    BLADE_DR2(aux0, aux2); /* Ds */
    BLADE_DR3(aux0, aux3); /* Dt */
    for(bfam_locidx_t i = 0; i < Np; i++)
          dq[off+i] -= HALF*JI[off+i]
                        *BLADE_D3_AP( Jr1x1[off+i]*aux1[i]
                                    + Jr2x1[off+i]*aux2[i],
                                    + Jr3x1[off+i]*aux3[i]);

    /* q += MI*JI*ux*(Dr'*M*Jrx*q + Ds'*M*Jsx*q + Dt'*M*Jtx*q) */
    bfam_locidx_t n = 0;
#if DIM==3
    for(bfam_locidx_t k = 0; k < N+1; k++)
#endif
      for(bfam_locidx_t j = 0; j < N+1; j++)
        for(bfam_locidx_t i = 0; i < N+1; i++,n++)
        {
          aux1[n] = BLADE_D3_AP(w[i]*w[j],*w[k])*Jr1x1[off+n]*q[off+n];
          aux2[n] = BLADE_D3_AP(w[i]*w[j],*w[k])*Jr2x1[off+n]*q[off+n];
          BLADE_D3_OP(aux3[n] = w[i]*w[j] *w[k] *Jr3x1[off+n]*q[off+n]);
        }
    BLADE_DR1T   (aux1, aux0);
    BLADE_DR2T_PE(aux2, aux0);
    BLADE_DR3T_PE(aux3, aux0);

    n = 0;
#if DIM==3
    for(bfam_locidx_t k = 0; k < N+1; k++)
#endif
      for(bfam_locidx_t j = 0; j < N+1; j++)
        for(bfam_locidx_t i = 0; i < N+1; i++,n++)
        {
          const bfam_real_t wi_ijk = BLADE_D3_AP(wi[i]*wi[j],*wi[k]);
          dq[off+n] += HALF*wi_ijk*JI[off+n]*ux[off+n]*aux0[n];
        }


    /****** Y-DERIVATIVE ******/
    /* q -= JI*(Jry*Dr*(uy*q) + Jsy*Ds*(uy*q) + Jty*Dt*(uy*q)) */
    for(bfam_locidx_t i = 0; i < Np; i++) aux0[i] = uy[off+i] * q[off+i];
    BLADE_DR1(aux0, aux1); /* Dr */
    BLADE_DR2(aux0, aux2); /* Ds */
    BLADE_DR3(aux0, aux3); /* Dt */
    for(bfam_locidx_t i = 0; i < Np; i++)
          dq[off+i] -= HALF*JI[off+i]
                        *BLADE_D3_AP( Jr1x2[off+i]*aux1[i]
                                    + Jr2x2[off+i]*aux2[i],
                                    + Jr3x2[off+i]*aux3[i]);

    /* q += MI*JI*uy*(Dr'*M*Jry*q + Ds'*M*Jsy*q + Dt'*M*Jty*q) */
    n = 0;
#if DIM==3
    for(bfam_locidx_t k = 0; k < N+1; k++)
#endif
      for(bfam_locidx_t j = 0; j < N+1; j++)
        for(bfam_locidx_t i = 0; i < N+1; i++,n++)
        {
          aux1[n] = BLADE_D3_AP(w[i]*w[j],*w[k])*Jr1x2[off+n]*q[off+n];
          aux2[n] = BLADE_D3_AP(w[i]*w[j],*w[k])*Jr2x2[off+n]*q[off+n];
          BLADE_D3_OP(aux3[n] = w[i]*w[j] *w[k] *Jr3x2[off+n]*q[off+n]);
        }
    BLADE_DR1T   (aux1, aux0);
    BLADE_DR2T_PE(aux2, aux0);
    BLADE_DR3T_PE(aux3, aux0);

    n = 0;
#if DIM==3
    for(bfam_locidx_t k = 0; k < N+1; k++)
#endif
      for(bfam_locidx_t j = 0; j < N+1; j++)
        for(bfam_locidx_t i = 0; i < N+1; i++,n++)
        {
          const bfam_real_t wi_ijk = BLADE_D3_AP(wi[i]*wi[j],*wi[k]);
          dq[off+n] += HALF*wi_ijk*JI[off+n]*uy[off+n]*aux0[n];
        }

#if DIM==3
    /****** Z-DERIVATIVE ******/
    /* q -= JI*(Jrz*Dr*(uz*q) + Jsz*Ds*(uz*q) + Jtz*Dt*(uz*q)) */
    for(bfam_locidx_t i = 0; i < Np; i++) aux0[i] = uz[off+i] * q[off+i];
    BLADE_DR1(aux0, aux1); /* Dr */
    BLADE_DR2(aux0, aux2); /* Ds */
    BLADE_DR3(aux0, aux3); /* Dt */
    for(bfam_locidx_t i = 0; i < Np; i++)
          dq[off+i] -= HALF*JI[off+i]
                        *(  Jr1x3[off+i]*aux1[i]
                          + Jr2x3[off+i]*aux2[i]
                          + Jr3x3[off+i]*aux3[i]);

    /* q += MI*JI*uz*(Dr'*M*Jrz*q + Ds'*M*Jsz*q + Dt'*M*Jtz*q) */
    n = 0;
    for(bfam_locidx_t k = 0; k < N+1; k++)
      for(bfam_locidx_t j = 0; j < N+1; j++)
        for(bfam_locidx_t i = 0; i < N+1; i++,n++)
        {
          aux1[n] = w[i]*w[j]*w[k]*Jr1x3[off+n]*q[off+n];
          aux2[n] = w[i]*w[j]*w[k]*Jr2x3[off+n]*q[off+n];
          aux3[n] = w[i]*w[j] *w[k]*Jr3x3[off+n]*q[off+n];
        }
    BLADE_DR1T   (aux1, aux0);
    BLADE_DR2T_PE(aux2, aux0);
    BLADE_DR3T_PE(aux3, aux0);

    n = 0;
    for(bfam_locidx_t k = 0; k < N+1; k++)
      for(bfam_locidx_t j = 0; j < N+1; j++)
        for(bfam_locidx_t i = 0; i < N+1; i++,n++)
        {
          const bfam_real_t wi_ijk = wi[i]*wi[j]*wi[k];
          dq[off+n] += HALF*wi_ijk*JI[off+n]*uz[off+n]*aux0[n];
        }
#endif

    /* loop over faces */
    for(bfam_locidx_t face = 0; face < Nfaces; face++)
    {
      for(bfam_locidx_t pnt = 0; pnt < Nfp; pnt++)
      {
        const bfam_locidx_t f = pnt + Nfp*(face + Nfaces*e);
        const bfam_locidx_t iM = vmapM[f];
        const bfam_locidx_t iP = vmapP[f];

        /* we use the average here in case there is a slight miss-match in the
         * velocities (they should be the same value though)
         */
        const bfam_real_t u = HALF*BLADE_D3_AP( (ux[iM] + ux[iP]) * n1[f]
                                              + (uy[iM] + uy[iP]) * n2[f],
                                              + (uz[iM] + uz[iP]) * n3[f]);

        /* intra */
        bfam_real_t uqS;
        blade_dgx_compute_flux(u,q[iM],q[iP],&uqS);
        blade_dgx_add_flux(1, dq, iM, uqS, sJ[f], JI[iM],wi[0]);
      }
    }
  }
}

static inline void
blade_dgx_remove_flux( const int inN,
    bfam_locidx_t face, bfam_locidx_t e, const bfam_locidx_t *vmapM,
    const bfam_real_t *n1, const bfam_real_t *n2,
#if DIM==3
    const bfam_real_t *n3,
#endif
    const bfam_real_t *ux, const bfam_real_t *uy,
#if DIM==3
    const bfam_real_t *uz,
#endif
    const bfam_real_t *sJ, const bfam_real_t *JI,   const bfam_real_t *wi,
    const bfam_real_t * q,
          bfam_real_t *dq)
{
  GENERIC_INIT(inN,blade_dgx_remove_flux);

  for(bfam_locidx_t pnt = 0; pnt < Nfp; pnt++)
  {
    bfam_locidx_t f = pnt + Nfp*(face + Nfaces*e);
    bfam_locidx_t iM = vmapM[f];
    const bfam_real_t u = BLADE_D3_AP( ux[iM] * n1[f]
                                     + uy[iM] * n2[f],
                                     + uz[iM] * n3[f]);
    bfam_real_t uqS;
    blade_dgx_compute_flux(u,q[iM],q[iM],&uqS);
    blade_dgx_add_flux(-1, dq, iM, uqS, sJ[f], JI[iM],wi[0]);
  }
}

static inline void
blade_project_flux(bfam_real_t *qs, const bfam_real_t *qg,
               bfam_locidx_t  inN,     bfam_locidx_t  Nrpg,
           const bfam_real_t * P1
#if DIM==3
          ,const bfam_real_t * P2
#endif
           )
{
  GENERIC_INIT(inN,blade_project_flux);

  /* zero out the components */
  for(int i = 0; i < Nfp; i++)
  {
    qs[i] = 0;
  }

  for(int m = 0, p = 0; m < Nrpg; m++)
#if DIM==3
    for(int l = 0; l < Nrpg; l++)
#endif
    {
      for(int j = 0, n = 0; j < N+1; j++)
#if DIM==3
        for(int i = 0; i < N+1; i++)
#endif
        {
          const bfam_real_t PR = BLADE_D3_AP(P1[j+m*(N+1)],
                                            *P2[i+l*(N+1)]);
          qs[n] += PR*qg[p];
          n++;
        }
      p++;
    }
}

static inline void
blade_massproject_flux(bfam_real_t *qs, const bfam_real_t *qg,
               bfam_locidx_t  inN,     bfam_locidx_t  Nrpg,
           const bfam_real_t * MP1,
#if DIM==3
           const bfam_real_t * MP2,
#endif
           const bfam_real_t  *wi)
{
  GENERIC_INIT(inN,blade_massproject_flux);

  /* zero out the components */
  for(int i = 0; i < Nfp; i++)
  {
    qs[i] = 0;
  }

  for(int m = 0, p = 0; m < Nrpg; m++)
#if DIM==3
    for(int l = 0; l < Nrpg; l++)
#endif
    {
      for(int j = 0, n = 0; j < N+1; j++)
#if DIM==3
        for(int i = 0; i < N+1; i++)
#endif
        {
          const bfam_real_t wi_MP = BLADE_D3_AP(wi[j]*MP1[j+m*(N+1)],
                                               *wi[i]*MP2[i+l*(N+1)]);
          qs[n] += wi_MP*qg[p];
          n++;
        }
      p++;
    }
}

void blade_dgx_inter_rhs_interface(
    int inN, bfam_subdomain_dgx_t *sub_g, const char *rate_prefix,
    const char *field_prefix, const bfam_long_real_t t)
{
  GENERIC_INIT(inN,blade_dgx_inter_rhs_interface);

  bfam_subdomain_dgx_t* sub_m =
    (bfam_subdomain_dgx_t*) sub_g->base.glue_m->sub_m;

  /* get the fields we will need */
  bfam_subdomain_dgx_glue_data_t* glue_m =
    (bfam_subdomain_dgx_glue_data_t*) sub_g->base.glue_m;
  bfam_subdomain_dgx_glue_data_t* glue_p =
    (bfam_subdomain_dgx_glue_data_t*) sub_g->base.glue_p;
  BFAM_ASSERT(glue_m != NULL);
  BFAM_ASSERT(glue_p != NULL);
  bfam_dictionary_t *fields_m    = &glue_m->base.fields;
  bfam_dictionary_t *fields_p    = &glue_p->base.fields;
  bfam_dictionary_t *fields_g    = &sub_g->base.fields;
  bfam_dictionary_t *fields      = &sub_m->base.fields;
  bfam_dictionary_t *fields_face = &sub_m->base.fields_face;

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(q  ,field_prefix,"q" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(q_M,field_prefix,"q" ,fields_m);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(q_P,field_prefix,"q" ,fields_p);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(u_g,field_prefix,"u" ,fields_g);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(dq ,rate_prefix ,"q" ,fields);

  /* get the material properties and metric terms */
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(ux ,field_prefix,"ux" ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(uy ,field_prefix,"uy" ,fields);
  BLADE_D3_OP(BFAM_LOAD_FIELD_RESTRICT_ALIGNED(uz,field_prefix,"uz",fields));
//JK  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(J   ,"","_grid_J"  ,fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(JI  ,"","_grid_JI" ,fields);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n1,"","_grid_nx0",fields_face);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n2,"","_grid_nx1",fields_face);
  BLADE_D3_OP(BFAM_LOAD_FIELD_RESTRICT_ALIGNED(n3,"","_grid_nx2",fields_face));
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(sJ  ,"","_grid_sJ",fields_face);

  bfam_real_t *wi  = sub_m->wi;
  BFAM_ASSUME_ALIGNED(wi ,32);

  BFAM_ASSERT(glue_p->EToEm);
  BFAM_ASSERT(glue_p->EToFm);
  BFAM_ASSERT(glue_p->EToHm);

  for(bfam_locidx_t le = 0; le < sub_g->K; le++)
  {
    bfam_locidx_t e = glue_p->EToEm[le];
    int8_t face = glue_p->EToFm[le];

    if(glue_p->EToHm[le] < 2)
      blade_dgx_remove_flux(N,face,e,sub_m->vmapM,
          n1,n2,
#if DIM==3
          n3,
#endif
          ux,uy,
#if DIM==3
          uz,
#endif
          sJ,JI,wi,
          q,dq);

#ifndef USE_GENERIC
#undef Np
#undef N
#endif
    bfam_locidx_t Np_g  = sub_g->Np;
    bfam_locidx_t Nrp_g = sub_g->N+1;
#ifndef USE_GENERIC
#define Np Np_BACK
#define N  NORDER
#endif

    bfam_real_t uqS_g[Np_g];

    /* Compute the flux on the glue grid */
    for(bfam_locidx_t pnt = 0; pnt < Np_g; pnt++)
    {
      bfam_locidx_t iG = pnt+le*Np_g;
      blade_dgx_compute_flux(u_g[iG], q_M[iG], q_P[iG], &uqS_g[pnt]);
    }


    /* these will be used to the store the projected values if we need them */
    BFAM_ALIGN(32) bfam_real_t uqS_M_STORAGE[Nfp];

    /* For now set the pointer to the storage */
    bfam_real_t *restrict uqS_M = uqS_M_STORAGE;

    /* check to see if projection */
    if(MASSPROJECTION && (!glue_p->same_order || glue_p->EToHm[le]
                                              || glue_p->EToHp[le]))
    {
      BFAM_ASSERT(glue_m->massprojection);
#if   DIM == 2
      const bfam_real_t *MP1 = glue_m->massprojection[glue_p->EToHm[le]];
#elif DIM == 3
      const int I1 = (glue_p->EToHm[le] == 0) ? 0 : (glue_p->EToHm[le]-1)/2+1;
      const int I2 = (glue_p->EToHm[le] == 0) ? 0 : (glue_p->EToHm[le]-1)%2+1;
      const bfam_real_t *MP1 = glue_m->massprojection[I1];
      const bfam_real_t *MP2 = glue_m->massprojection[I2];
#else
#error "Bad Dimension"
#endif
              blade_massproject_flux(uqS_M, uqS_g, N, Nrp_g,
                  MP1,
#if DIM==3
                  MP2,
#endif
                                     wi);
    }
    else if(PROJECTION  && (!glue_p->same_order || glue_p->EToHm[le] || glue_p->EToHp[le]))
    {
      BFAM_ASSERT(glue_m->projection);
#if   DIM == 2
      const bfam_real_t *P1 = glue_m->projection[glue_p->EToHm[le]];
#elif DIM == 3
      const int I1 = (glue_p->EToHm[le] == 0) ? 0 : (glue_p->EToHm[le]-1)/2+1;
      const int I2 = (glue_p->EToHm[le] == 0) ? 0 : (glue_p->EToHm[le]-1)%2+1;
      const bfam_real_t *P1 = glue_m->projection[I1];
      const bfam_real_t *P2 = glue_m->projection[I2];
#else
#error "Bad Dimension"
#endif
      blade_project_flux(uqS_M, uqS_g, N, Nrp_g,
          P1
#if DIM==3
          ,P2
#endif
          );
    }
    else
      uqS_M = uqS_g;

    for(bfam_locidx_t pnt = 0; pnt < Nfp; pnt++)
    {
      bfam_locidx_t f = pnt + Nfp*(face + Nfaces*e);
      bfam_locidx_t iM = sub_m->vmapM[f];
      blade_dgx_add_flux(1, dq, iM, uqS_M[pnt], sJ[f], JI[iM],wi[0]);
    }
  }
}

void blade_dgx_add_rates_advection(
    int inN, bfam_subdomain_dgx_t *sub, const char *field_prefix_lhs,
    const char *field_prefix_rhs, const char *rate_prefix,
    const bfam_long_real_t a)
{
  GENERIC_INIT(inN,blade_dgx_add_rates_advection);

  const char *f_names[] = {"q", NULL};

  const bfam_locidx_t num_pts = sub->K * Np;

  /* get the fields we will need */
  bfam_dictionary_t *fields = &sub->base.fields;

  for(bfam_locidx_t f = 0; f_names[f] != NULL; f++)
  {
    BFAM_LOAD_FIELD_ALIGNED(         lhs ,field_prefix_lhs,f_names[f],fields);
    BFAM_LOAD_FIELD_ALIGNED(         rhs ,field_prefix_rhs,f_names[f],fields);
    BFAM_LOAD_FIELD_RESTRICT_ALIGNED(rate,rate_prefix     ,f_names[f],fields);
    for(bfam_locidx_t n = 0; n < num_pts; n++) lhs[n] = rhs[n] + a*rate[n];
  }
}

#endif
