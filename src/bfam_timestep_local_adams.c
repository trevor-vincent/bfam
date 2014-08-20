#include <bfam_timestep_local_adams.h>
#include <bfam_log.h>

#define BFAM_LOCAL_ADAMS_PREFIX ("_local_adams_rate_")
#define BFAM_LOCAL_ADAMS_LVL_PREFIX "_local_adams_lvl_"
#define BFAM_LOCAL_ADAMS_COMM_LVL_PREFIX ("_local_adams_comm_lvl_")

/*
#define BFAM_LOCAL_ADAMS_ALWAYS_INTERP
*/

typedef struct bfam_ts_local_allprefix
{
  bfam_ts_local_adams_t *ts;
  bfam_locidx_t    step; /* what number step are we on */
  bfam_locidx_t    lvl;  /* this marks the level to the updated */
  bfam_long_real_t dt;   /* this is the fastest time step */
} bfam_ts_local_adams_allprefix_t;

void
bfam_ts_local_adams_fill_level_tag(char* tag, size_t buf_sz, int level)
{
  snprintf(tag,buf_sz,"%s%d",BFAM_LOCAL_ADAMS_LVL_PREFIX,level);
}

static int
get_tag_level_number(const char * key, void *val)
{
  BFAM_ASSERT(*(bfam_locidx_t*)val < 0);
  return sscanf(key, BFAM_LOCAL_ADAMS_LVL_PREFIX"%"BFAM_LOCIDX_PRId,
      (bfam_locidx_t*) val);
}

void
bfam_ts_local_adams_fill_comm_level_tag(char* tag, size_t buf_sz, int level)
{
  snprintf(tag,buf_sz,"%s%d",BFAM_LOCAL_ADAMS_COMM_LVL_PREFIX,level);
}

static void
comm_send_prefix(bfam_subdomain_t *sub,
    char *prefix, size_t buf_siz, void* user_data)
{
  BFAM_ASSERT(sub->glue_m);
  BFAM_ASSERT(sub->glue_p);

  /* Get the plus and minus side levels */
  bfam_locidx_t m_lvl = -1;
  bfam_critbit0_allprefixed(&sub->glue_m->tags,
      BFAM_LOCAL_ADAMS_LVL_PREFIX, get_tag_level_number,&m_lvl);
  BFAM_ASSERT(m_lvl>=0);

  bfam_locidx_t p_lvl = -1;
  bfam_critbit0_allprefixed(&sub->glue_p->tags,
      BFAM_LOCAL_ADAMS_LVL_PREFIX, get_tag_level_number,&p_lvl);
  BFAM_ASSERT(p_lvl>=0);

  bfam_ts_local_adams_allprefix_t* data =
    (bfam_ts_local_adams_allprefix_t*) user_data;
  BFAM_ASSERT(data);

#ifdef BFAM_LOCAL_ADAMS_ALWAYS_INTERP
  /*
   * If my level number of greater than (or equal to) my neighbors, then I need
   * to send rates to my neighbor because I step slower.
   *
   * (The one step is b/c the first time we want to exchange the initial
   * condition)
   */
  if(data->ts->numStepsArray[m_lvl] > 0 && m_lvl >= p_lvl)
#else
  /*
   * If my level number of greater than my neighbors and my level number of
   * greater than the level being updated (i.e., I am not being updated), then I
   * need to send rates to my neighbor because I step slower.
   */
  if(m_lvl > p_lvl && data->lvl < m_lvl)
#endif
  {
    /*
     *
     * We have to shift back on b/c the stage number has already been moved
     * forward by 1
     */
    snprintf(prefix,buf_siz,"%s%d_",BFAM_LOCAL_ADAMS_PREFIX,
        (data->ts->currentStageArray[m_lvl]+data->ts->nStages-1)
        %data->ts->nStages);
    BFAM_LDEBUG("send: lvl = %d :: m_lvl = %d :: p_lvl = %d prefix = %s",
        data->lvl, m_lvl, p_lvl, prefix);
  }
}

static void
comm_recv_prefix(bfam_subdomain_t *sub,
    char *prefix, size_t buf_siz, void* user_data)
{
  BFAM_ASSERT(sub->glue_m);
  BFAM_ASSERT(sub->glue_p);

  /* Get the plus and minus side levels */
  bfam_locidx_t m_lvl = -1;
  bfam_critbit0_allprefixed(&sub->glue_m->tags,
      BFAM_LOCAL_ADAMS_LVL_PREFIX, get_tag_level_number,&m_lvl);
  BFAM_ASSERT(m_lvl>=0);

  bfam_locidx_t p_lvl = -1;
  bfam_critbit0_allprefixed(&sub->glue_p->tags,
      BFAM_LOCAL_ADAMS_LVL_PREFIX, get_tag_level_number,&p_lvl);
  BFAM_ASSERT(p_lvl>=0);

  bfam_ts_local_adams_allprefix_t* data =
    (bfam_ts_local_adams_allprefix_t*) user_data;
  BFAM_ASSERT(data);

#ifdef BFAM_LOCAL_ADAMS_ALWAYS_INTERP
  /*
   * If my level number of smaller than (or equal to) my neighbors, then I need
   * to send rates to my neighbor because I step faster.
   *
   * (The one step is b/c the first time we want to exchange the initial
   * condition)
   */
  if(data->ts->numStepsArray[m_lvl] > 0 && m_lvl <= p_lvl)
#else
  /*
   * If my level number of smaller than my neighbors and my neighbors level
   * number of greater than the level number (i.e., my neighbor is not
   * updating), then I need to send rates to my neighbor because I step faster.
   */
  if(m_lvl < p_lvl && data->lvl < p_lvl)
#endif
  {
   /*
   * We have to shift back on b/c the stage number has already been moved
   * forward by 1
   */
    snprintf(prefix,buf_siz,"%s%d_",BFAM_LOCAL_ADAMS_PREFIX,
        (data->ts->currentStageArray[p_lvl]+data->ts->nStages-1)
        %data->ts->nStages);
    BFAM_LDEBUG("recv: lvl = %d :: m_lvl = %d :: p_lvl = %d prefix = %s",
        data->lvl, m_lvl, p_lvl, prefix);
  }
}

bfam_ts_local_adams_t*
bfam_ts_local_adams_new(bfam_domain_t* dom, bfam_ts_local_adams_method_t method,
    bfam_locidx_t num_lvl, bfam_domain_match_t subdom_match,
    const char** subdom_tags, bfam_domain_match_t comm_match, const char**
    comm_tags, MPI_Comm mpicomm, int mpitag,
    bfam_subdomain_comm_args_t *comm_data,
    void (*aux_rates) (bfam_subdomain_t *thisSubdomain, const char *prefix),
    void (*glue_rates) (bfam_subdomain_t *thisSubdomain, const char *prefix),
    void (*scale_rates) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const bfam_long_real_t a),
    void (*intra_rhs) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const char *field_prefix,
      const bfam_long_real_t t),
    void (*inter_rhs) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const char *minus_rate_prefix,
      const char *field_prefix, const bfam_long_real_t t),
    void (*add_rates) (bfam_subdomain_t *thisSubdomain,
      const char *field_prefix_lhs, const char *field_prefix_rhs,
      const char *rate_prefix, const bfam_long_real_t a),
    void (*add_rates_glue_p) (bfam_subdomain_t *thisSubdomain,
      const char *field_prefix_lhs, const char *field_prefix_rhs,
      const char *rate_prefix, const bfam_long_real_t a),
    const int RK_init)
{
  bfam_ts_local_adams_t* newTS = bfam_malloc(sizeof(bfam_ts_local_adams_t));
  bfam_ts_local_adams_init(newTS, dom, method, num_lvl, subdom_match,
      subdom_tags, comm_match, comm_tags, mpicomm, mpitag, comm_data, aux_rates,
      glue_rates, scale_rates,intra_rhs,inter_rhs,add_rates,add_rates_glue_p,
      RK_init);
  return newTS;
}

static int
bfam_ts_local_adams_intra_rhs(const char * key, void *val, void *arg)
{
  bfam_ts_local_adams_allprefix_t *data =
    (bfam_ts_local_adams_allprefix_t *) arg;
  bfam_subdomain_t* sub = (bfam_subdomain_t*) val;

  /* We have to determine if this domain is updated or not */
  bfam_locidx_t lvl = -1;
  bfam_critbit0_allprefixed(&sub->tags, BFAM_LOCAL_ADAMS_LVL_PREFIX,
      get_tag_level_number,&lvl);
  char *rate_prefix = NULL;
  char rate_prefix_storage[BFAM_BUFSIZ];
  if(lvl > -1 && lvl <= data->lvl)
  {
    rate_prefix = rate_prefix_storage;
    snprintf(rate_prefix_storage,BFAM_BUFSIZ,"%s%d_",BFAM_LOCAL_ADAMS_PREFIX,
        data->ts->currentStageArray[lvl]%data->ts->nStages);
    BFAM_LDEBUG("Local Adams intra: level %"BFAM_LOCIDX_PRId
        " using rate prefix %s",lvl,rate_prefix_storage);
    data->ts->scale_rates(sub, rate_prefix, 0);
    data->ts->intra_rhs(sub, rate_prefix, "", data->ts->t);
  }

  return 1;
}

static void
interp_fields(bfam_ts_local_adams_allprefix_t* data,
    bfam_subdomain_t* sub, bfam_locidx_t m_lvl, bfam_locidx_t p_lvl)
{
  /*
   * First we we interp, then we don't have our current fields so first we get
   * those
   */
  bfam_subdomain_comm_args_t* sub_comm_data =
    (bfam_subdomain_comm_args_t*) data->ts->comm_array[m_lvl]->user_args;
  BFAM_ASSERT(sub_comm_data);
  BFAM_ASSERT(sub_comm_data->user_data == NULL);
  BFAM_ASSERT(sub_comm_data->user_prefix_function == NULL);
  sub->glue_put_send_buffer(sub, NULL, 0, sub_comm_data);

  /* Determine my step size and my neighbors */
  const bfam_locidx_t m_sz = (1<<m_lvl);
  const bfam_locidx_t p_sz = (1<<p_lvl);

  /* Assert the minus step size is faster than the plus size */
#ifdef BFAM_LOCAL_ADAMS_ALWAYS_INTERP
  BFAM_ASSERT(m_lvl <= p_lvl);
#else
  BFAM_ASSERT(m_lvl < p_lvl && p_lvl > data->lvl);
#endif

  /* figure out when the minus and plus sides last updated */
  const bfam_locidx_t m_last_step = data->step - m_sz;
  const bfam_locidx_t p_last_step = p_sz*(((data->step+p_sz-1)/p_sz)-1);

  /* Now get the integration limits */
  const bfam_locidx_t a = m_last_step-p_last_step;
  const bfam_locidx_t b = a + m_sz;
  const bfam_locidx_t d = p_sz;

  BFAM_LDEBUG(
      "step %2"BFAM_LOCIDX_PRId
      " :: lvl %2"BFAM_LOCIDX_PRId
      " :: ngh lvl %2"BFAM_LOCIDX_PRId
      " :: last  %2"BFAM_LOCIDX_PRId
      " :: ngh last  %2"BFAM_LOCIDX_PRId
      " :: d = %"BFAM_LOCIDX_PRId
      " :: a = %"BFAM_LOCIDX_PRId
      " :: b = %"BFAM_LOCIDX_PRId,
      data->step, m_lvl, p_lvl, m_last_step, p_last_step, d,a,b);

  /*
   * we have to look at the p_lvl because these are the guys we are
   * interpolating
   */
  bfam_long_real_t A[4];
  bfam_locidx_t    num_stages = BFAM_MIN(data->ts->numStepsArray[p_lvl],
                                         data->ts->nStages);
  switch(num_stages)
  {
    case 1:
      {
        /* \int_{a}^{b} 1 dt */
        A[0] = (bfam_long_real_t)(b-a);
      }
      break;
    case 2:
      {
        /* Integrate[(t + d)/(0 + d), {t, a, b}] // Together */
        A[0] = (bfam_long_real_t)(-a*a + b*b - 2*a*d + 2*b*d) /
               (bfam_long_real_t)(2*d);

        /* Integrate[(t)/(-d), {t, a, b}] // Together */
        A[1] = (bfam_long_real_t)(a*a-b*b) / (bfam_long_real_t)(2*d);
      }
      break;
    case 3:
      {
        /* Integrate[(t + d) (t + 2 d)/((0 + d) (0 + 2 d)), {t, a, b}] //  Together */
        A[0] = (bfam_long_real_t)(-2*a*a*a + 2*b*b*b - 9*a*a*d
                                 + 9*b*b*d - 12*a*d*d + 12*b*d*d) /
               (bfam_long_real_t)(12*d*d);

        /* Integrate[(t) (t + 2 d)/((-d) (-d + 2 d)), {t, a, b}] // Together */
        A[1] = (bfam_long_real_t)(a*a*a - b*b*b + 3*a*a*d - 3*b*b*d) /
               (bfam_long_real_t)(3*d*d);

        /* Integrate[(t) (t + d)/((-2 d) (-2 d + d)), {t, a, b}] // Together */
        A[2] = (bfam_long_real_t)(-2*a*a*a + 2*b*b*b - 3*a*a*d + 3*b*b*d) /
               (bfam_long_real_t)(12*d*d);

        break;
      }
    case 4:
      {
        /* Integrate[(t + d) (t + 2 d) (t + 3 d)/((0 + d) (0 + 2 d) (0 + 3 d)), {t, a, b}] // Together */
        A[0] = (bfam_long_real_t)(-a*a*a*a + b*b*b*b - 8*a*a*a*d + 8*b*b*b*d
                                  - 22*a*a*d*d + 22*b*b*d*d
                                  - 24*a*d*d*d + 24*b*d*d*d) /
               (bfam_long_real_t)(24*d*d*d);

        /* Integrate[(t) (t + 2 d) (t + 3 d)/((-d) (-d + 2 d) (-d + 3 d)), {t, a, b}] // Together */
        A[1] = (bfam_long_real_t)(3*a*a*a*a - 3*b*b*b*b + 20*a*a*a*d
                                 - 20*b*b*b*d + 36*a*a*d*d - 36*b*b*d*d) /
               (bfam_long_real_t)(24*d*d*d);

        /* Integrate[(t) (t + d) (t + 3 d)/((-2 d) (-2 d + d) (-2 d + 3 d)), {t, a, b}] // Together */
        A[2] = (bfam_long_real_t)(-3*a*a*a*a + 3*b*b*b*b - 16*a*a*a*d
                                  + 16*b*b*b*d - 18*a*a*d*d + 18*b*b*d*d) /
               (bfam_long_real_t)(24*d*d*d);

        /* Integrate[(t) (t + d) (t + 2 d)/((-3 d) (-3 d + d) (-3 d + 2 d)), {t, a, b}] // Together */
        A[3] = (bfam_long_real_t)(a*a*a*a - b*b*b*b + 4*a*a*a*d - 4*b*b*b*d
                                  + 4*a*a*d*d - 4*b*b*d*d) /
               (bfam_long_real_t)(24*d*d*d);
      }
      break;
    default:
      BFAM_ABORT("Adams-Bashforth order %d not implemented",num_stages);
  }

  for(int k = 0; k < num_stages;k++)
  {
    char prefix[BFAM_BUFSIZ];
    /* minus 1 here b/c we need the last guys (not the current guys */
    snprintf(prefix,BFAM_BUFSIZ,"%s%d_",BFAM_LOCAL_ADAMS_PREFIX,
        (data->ts->currentStageArray[p_lvl]+data->ts->nStages-k-1)
        %data->ts->nStages);
    BFAM_LDEBUG("Adams interp step: stage %d of %d using prefix %s",k,
        num_stages, prefix);
    data->ts->add_rates_glue_p(sub, "", "", prefix, data->dt*A[k]);
  }
}

static int
bfam_ts_local_adams_inter_rhs(const char * key, void *val, void *arg)
{
  bfam_ts_local_adams_allprefix_t *data =
    (bfam_ts_local_adams_allprefix_t *) arg;
  bfam_subdomain_t* sub = (bfam_subdomain_t*) val;

  /* We have to determine if this domain is updated or not */
  bfam_locidx_t lvl = -1;
  bfam_critbit0_allprefixed(&sub->tags, BFAM_LOCAL_ADAMS_LVL_PREFIX,
      get_tag_level_number,&lvl);
  char *rate_prefix = NULL;
  char rate_prefix_storage[BFAM_BUFSIZ];
  if(lvl > -1 && lvl <= data->lvl)
  {
    rate_prefix = rate_prefix_storage;
    snprintf(rate_prefix_storage,BFAM_BUFSIZ,"%s%d_",BFAM_LOCAL_ADAMS_PREFIX,
        data->ts->currentStageArray[lvl]%data->ts->nStages);
    BFAM_LDEBUG("Local Adams inter: level %"BFAM_LOCIDX_PRId
        " using rate prefix %s",lvl,rate_prefix_storage);
  }

  /*
   * Determine if the minus side exists and whether interpolation needs to
   * occur on the plus side
   */
  bfam_locidx_t m_lvl = -1;
  if(sub->glue_m && sub->glue_m->sub_m)
    bfam_critbit0_allprefixed(&sub->glue_m->sub_m->tags,
        BFAM_LOCAL_ADAMS_LVL_PREFIX, get_tag_level_number,&m_lvl);
  char *minus_rate_prefix = NULL;
  char minus_rate_prefix_storage[BFAM_BUFSIZ];
  if(m_lvl > -1 && m_lvl <= data->lvl)
  {
    minus_rate_prefix = minus_rate_prefix_storage;
    snprintf(minus_rate_prefix_storage,BFAM_BUFSIZ,"%s%d_",
        BFAM_LOCAL_ADAMS_PREFIX,
        data->ts->currentStageArray[m_lvl]%data->ts->nStages);
    BFAM_LDEBUG("Local Adams inter: level %"BFAM_LOCIDX_PRId
        " using minus rate prefix %s",m_lvl,minus_rate_prefix_storage);

    bfam_locidx_t p_lvl = -1;
    if(sub->glue_p)
      bfam_critbit0_allprefixed(&sub->glue_p->tags,
          BFAM_LOCAL_ADAMS_LVL_PREFIX, get_tag_level_number,&p_lvl);

#ifdef BFAM_LOCAL_ADAMS_ALWAYS_INTERP
    /*
     * If my update rate is faster than or equal to my neighbors then
     * interpolate (assume we have at least exchanged the initial condition)
     */
    if(data->ts->numStepsArray[m_lvl] > 0 && m_lvl <= p_lvl)
#else
    /*
     * If my level is less than my neighbor and my neighbor is not being updated
     * then do the interpolation
     */
    if(m_lvl < p_lvl && p_lvl > data->lvl)
#endif
    {
      interp_fields(data,sub,m_lvl,p_lvl);
    }
  }

  if(rate_prefix || minus_rate_prefix)
    data->ts->inter_rhs(sub, rate_prefix, minus_rate_prefix, "", data->ts->t);

  return 1;
}

static inline int
bfam_ts_local_adams_do_update(bfam_subdomain_t* sub, const bfam_long_real_t* A,
    const bfam_ts_local_adams_t* ts, const bfam_long_real_t dt,
    const int nStages, const bfam_locidx_t lvl)
{
  BFAM_LDEBUG("BFAM_TS_LOCAL_ADAMS_DO_UPDATE");

  /* Loop through the stages to scale rates and add in */
  /*
   * nStages is the computing number of stages whereas ts->nStages is the
   * storage number of stages
   */
  for(int k = 0; k < nStages;k++)
  {
    char prefix[BFAM_BUFSIZ];
    snprintf(prefix,BFAM_BUFSIZ,"%s%d_",BFAM_LOCAL_ADAMS_PREFIX,
        (ts->currentStageArray[lvl]+ts->nStages-k)%ts->nStages);
    BFAM_LDEBUG("Adams step: stage %d of %d using prefix %s",k,nStages,prefix);
    ts->add_rates(sub, "", "", prefix, dt*A[k]);
  }
  return 1;
}

static int
bfam_ts_local_adams_update(const char * key, void *val, void *arg)
{
  bfam_ts_local_adams_allprefix_t *data =
    (bfam_ts_local_adams_allprefix_t *) arg;
  bfam_subdomain_t* sub = (bfam_subdomain_t*) val;
  bfam_ts_local_adams_t* ts = data->ts;

  bfam_locidx_t lvl = -1;
  bfam_critbit0_allprefixed(&sub->tags, BFAM_LOCAL_ADAMS_LVL_PREFIX,
      get_tag_level_number,&lvl);

  if(lvl >= 0 && !(data->step%(1<<lvl)))
    switch(BFAM_MIN(ts->numStepsArray[lvl]+1,ts->nStages))
    {
      case 1:
        {
          bfam_long_real_t A[1] = {BFAM_LONG_REAL(1.0)};
          bfam_ts_local_adams_do_update(sub, A, ts, (1<<lvl)*data->dt, 1, lvl);
        }
        break;
      case 2:
        {
          bfam_long_real_t A[2] = {
            BFAM_LONG_REAL( 3.0) / BFAM_LONG_REAL( 2.0),
            BFAM_LONG_REAL(-1.0) / BFAM_LONG_REAL( 2.0),
          };
          bfam_ts_local_adams_do_update(sub, A, ts, (1<<lvl)*data->dt, 2, lvl);
        }
        break;
      case 3:
        {
          bfam_long_real_t A[3] = {
            BFAM_LONG_REAL(23.0)/ BFAM_LONG_REAL(12.0),
            BFAM_LONG_REAL(-4.0)/ BFAM_LONG_REAL( 3.0),
            BFAM_LONG_REAL( 5.0)/ BFAM_LONG_REAL(12.0),
          };
          bfam_ts_local_adams_do_update(sub, A, ts, (1<<lvl)*data->dt, 3, lvl);
        }
        break;
      case 4:
        {
          bfam_long_real_t A[4] = {
            BFAM_LONG_REAL( 55.0)/ BFAM_LONG_REAL( 24.0),
            BFAM_LONG_REAL(-59.0)/ BFAM_LONG_REAL( 24.0),
            BFAM_LONG_REAL( 37.0)/ BFAM_LONG_REAL( 24.0),
            BFAM_LONG_REAL( -3.0)/ BFAM_LONG_REAL(  8.0),
          };
          bfam_ts_local_adams_do_update(sub, A, ts, (1<<lvl)*data->dt, 4, lvl);
        }
        break;
      default:
        BFAM_ABORT("Adams-Bashforth order %d not implemented",
            BFAM_MIN(ts->numStepsArray[lvl]+1,ts->nStages));
    }
  return 1;
}

static void
bfam_ts_local_adams_step(bfam_ts_t *a_ts, bfam_long_real_t dt)
{
  bfam_ts_local_adams_t* ts = (bfam_ts_local_adams_t*) a_ts;
  bfam_locidx_t num_steps = 1<<(ts->numLevels-1);
  BFAM_LDEBUG("Number of steps for the local time stepper %"BFAM_LOCIDX_PRId,
      num_steps);

  dt /= num_steps;
  bfam_ts_local_adams_allprefix_t data;
  data.ts = ts;
  data.dt = dt;
  data.lvl = -1;


  /* determine the level of comm to do: max of this update and last update */
  bfam_locidx_t last_lvl = 0;
  for(bfam_locidx_t step = 0; step < num_steps; step++)
  {
    data.step = step;

    BFAM_LDEBUG("local time step number %"BFAM_LOCIDX_PRId, step);
    BFAM_LDEBUG("Number of effective levels is %"BFAM_LOCIDX_PRId,
        ts->effNumLvls);
    data.lvl = 0;

    /* check to see what level we are stepping */
    for(bfam_locidx_t lvl = 0; lvl < ts->numLevels; lvl++)
    {
      bfam_locidx_t chk = 1 << lvl;
      if(!(step%chk))
      {
        BFAM_LDEBUG("level %"BFAM_LOCIDX_PRId" to be updated",lvl);
        data.lvl = BFAM_MAX(data.lvl, lvl);
      }
    }

    /* If effNumLvls is smaller than (or equal) to the level we are stepping by
     * above calculation, then we really want to step all the levels (since
     * slower guys are currently stepping faster). Thus, we set the data lvl to
     * the real number of levels minus 1
     */
    if(data.lvl >= ts->effNumLvls-1)
    {
      data.lvl = ts->numLevels-1;
    }

    /* comm level is the max of this level and the last level updated */
    bfam_locidx_t comm_lvl = BFAM_MAX(data.lvl, last_lvl);
    BFAM_LDEBUG("step %02"BFAM_LOCIDX_PRId
        ": update level %02"BFAM_LOCIDX_PRId
        ": and communication level %02"BFAM_LOCIDX_PRId,
        step, data.lvl, comm_lvl);

    /* start the communication for the communication level */
    bfam_subdomain_comm_args_t* sub_comm_data =
      (bfam_subdomain_comm_args_t*) ts->comm_array[comm_lvl]->user_args;

    /* This just makes sure we won't kill anything later, if we do then we need
     * to rethink what we are doing here
     */
    BFAM_ASSERT(sub_comm_data);
    BFAM_ASSERT(sub_comm_data->user_data == NULL);
    BFAM_ASSERT(sub_comm_data->user_prefix_function == NULL);

    sub_comm_data->user_prefix_function = comm_send_prefix;
    sub_comm_data->user_data = &data;
    bfam_communicator_start(ts->comm_array[comm_lvl]);

    /* Do the intra work for the levels to be updated */
    bfam_dictionary_allprefixed_ptr(&ts->elems,
        "",&bfam_ts_local_adams_intra_rhs,&data);

    /* finish the communication */
    sub_comm_data->user_prefix_function = comm_recv_prefix;
    bfam_communicator_finish(ts->comm_array[comm_lvl]);

    /* reset the pointer to NULL */
    sub_comm_data->user_data = NULL;
    sub_comm_data->user_prefix_function = NULL;

    /* Do the inter work for the levels to be updated */
    bfam_dictionary_allprefixed_ptr(&ts->elems,
        "",&bfam_ts_local_adams_inter_rhs,&data);

    if(ts->lsrk)
    {
      bfam_locidx_t stage = (ts->currentStageArray[0]+1)%ts->nStages;

      /*
       * If we are using RK, we want to tell the RK scheme to use the next rate
       * for storage (not the current rate which is valid)
       */
      ts->lsrk->t = ts->t;
      char rate_prefix[BFAM_BUFSIZ];
      /* THIS ASSUMES THAT ALL THE RATES ARE AT THE SAME LEVEL INITIALLY */
      snprintf(rate_prefix,BFAM_BUFSIZ,"%s%d_",BFAM_LOCAL_ADAMS_PREFIX,
          stage);
      BFAM_LDEBUG("Local Adams step: RK rate rate_prefix %s",rate_prefix);
      BFAM_ASSERT(ts->lsrk->step_extended);
      ts->lsrk->step_extended((bfam_ts_t*)ts->lsrk,dt,rate_prefix,"","");
      ts->numLSRKsteps++;

      /* Set the time step that each stage is at (or at least will be when the
       * next rate is calculated)
       */
      for(bfam_locidx_t k = 0; k < ts->numLevels;k++)
        ts->lvlStepArray[k][stage] = ts->numLSRKsteps;

      if(ts->numLSRKsteps+1 >= ts->nStages)
      {
        bfam_ts_lsrk_free(ts->lsrk);
        bfam_free(ts->lsrk);
        ts->lsrk = NULL;
      }

    }
    else
      /* Do the local update */
      bfam_dictionary_allprefixed_ptr(&ts->elems,
          "",&bfam_ts_local_adams_update,&data);

    /* set last update to next update */
    last_lvl = data.lvl;

    /* update the stage counters */
    for(bfam_locidx_t k=0; k <= data.lvl;k++)
    {
      ts->numStepsArray[k]++;
      ts->currentStageArray[k] = (ts->currentStageArray[k]+1)%ts->nStages;
    }
    ts->t += dt;
  }

}

void
bfam_ts_local_adams_init(
    bfam_ts_local_adams_t*       ts,
    bfam_domain_t*               dom,
    bfam_ts_local_adams_method_t method,
    bfam_locidx_t                num_lvl,
    bfam_domain_match_t          subdom_match,
    const char**                 subdom_tags,
    bfam_domain_match_t          comm_match,
    const char**                 comm_tags,
    MPI_Comm                     mpicomm,
    int                          mpitag,
    bfam_subdomain_comm_args_t*  comm_data,
    void (*aux_rates) (bfam_subdomain_t *thisSubdomain, const char *prefix),
    void (*glue_rates) (bfam_subdomain_t *thisSubdomain, const char *prefix),
    void (*scale_rates) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const bfam_long_real_t a),
    void (*intra_rhs) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix,
      const char *field_prefix, const bfam_long_real_t t),
    void (*inter_rhs) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const char *minus_rate_prefix,
      const char *field_prefix, const bfam_long_real_t t),
    void (*add_rates) (bfam_subdomain_t *thisSubdomain,
      const char *field_prefix_lhs, const char *field_prefix_rhs,
      const char *rate_prefix, const bfam_long_real_t a),
    void (*add_rates_glue_p) (bfam_subdomain_t *thisSubdomain,
      const char *field_prefix_lhs, const char *field_prefix_rhs,
      const char *rate_prefix, const bfam_long_real_t a),
    const int RK_init)
{
  BFAM_LDEBUG("LOCAL ADAMS INIT");
#ifdef BFAM_LOCAL_ADAMS_ALWAYS_INTERP
  BFAM_WARNING("BFAM_LOCAL_ADAMS_ALWAYS_INTERP is set which may lead to "
               "energy instabilities due to round-off error drift");
#endif

  /*
   * set up some preliminaries
   */
  bfam_ts_init(&ts->base, dom);
  bfam_dictionary_init(&ts->elems);
  ts->t  = BFAM_LONG_REAL(0.0);
  ts->base.step = &bfam_ts_local_adams_step;

  /*
   * store the function calls
   */
  ts->scale_rates      = scale_rates;
  ts->intra_rhs        = intra_rhs;
  ts->inter_rhs        = inter_rhs;
  ts->add_rates        = add_rates;
  ts->add_rates_glue_p = add_rates_glue_p;


  ts->lsrk         = NULL;
  ts->comm_array   = NULL;

  switch(method)
  {
    default:
      BFAM_WARNING("Invalid Adams scheme, using ADAMS_3");
    case BFAM_TS_LOCAL_ADAMS_3:
      ts->nStages = 3;
      /* if necessary initialize the RK scheme */
      if(RK_init)
      {
        ts->lsrk = bfam_ts_lsrk_new_extended(dom, BFAM_TS_LSRK_KC54,
            subdom_match, subdom_tags, comm_match, comm_tags, mpicomm, mpitag,
            comm_data, aux_rates, scale_rates, intra_rhs, inter_rhs, add_rates,
            0);
      }
      break;
    case BFAM_TS_LOCAL_ADAMS_1:
      ts->nStages = 1;
      /* if necessary initialize the RK scheme */
      if(RK_init)
      {
        ts->lsrk = bfam_ts_lsrk_new_extended(dom, BFAM_TS_LSRK_FE,
            subdom_match, subdom_tags, comm_match, comm_tags, mpicomm, mpitag,
            comm_data, aux_rates, scale_rates, intra_rhs, inter_rhs, add_rates,
            0);
      }
      break;
    case BFAM_TS_LOCAL_ADAMS_2:
      ts->nStages = 2;
      /* if necessary initialize the RK scheme */
      if(RK_init)
      {
        ts->lsrk = bfam_ts_lsrk_new_extended(dom, BFAM_TS_LSRK_HEUN,
            subdom_match, subdom_tags, comm_match, comm_tags, mpicomm, mpitag,
            comm_data, aux_rates, scale_rates, intra_rhs, inter_rhs, add_rates,
            0);
      }
      break;
    case BFAM_TS_LOCAL_ADAMS_4:
      ts->nStages = 4;
      /* if necessary initialize the RK scheme */
      if(RK_init)
      {
        ts->lsrk = bfam_ts_lsrk_new_extended(dom, BFAM_TS_LSRK_KC54,
            subdom_match, subdom_tags, comm_match, comm_tags, mpicomm, mpitag,
            comm_data, aux_rates, scale_rates, intra_rhs, inter_rhs, add_rates,
            0);
      }
      break;
  }

  ts->numLevels = num_lvl;
  ts->currentStageArray = bfam_malloc(ts->numLevels*sizeof(bfam_locidx_t));
  ts->numStepsArray     = bfam_malloc(ts->numLevels*sizeof(bfam_locidx_t));
  if(RK_init)
  {
    ts->effNumLvls = 1;
    ts->lvlStepArray = bfam_malloc(ts->numLevels*sizeof(bfam_locidx_t));
  }
  else
  {
    ts->lvlStepArray = NULL;
    ts->effNumLvls = num_lvl;
  }
  for(bfam_locidx_t k = 0; k < ts->numLevels; k++)
  {
    ts->numStepsArray[k] = 0;
    ts->currentStageArray[k] = 0;
    if(ts->lvlStepArray)
    {
      ts->lvlStepArray[k] = bfam_malloc(ts->numLevels*sizeof(bfam_locidx_t));
      for(bfam_locidx_t j = 0; j < ts->nStages; j++)
        ts->lvlStepArray[k][j] = 0;
    }
  }
  ts->numLSRKsteps = 0;


  /*
   * get the subdomains and create rates we will need
   */
  bfam_subdomain_t *subs[dom->numSubdomains+1];
  bfam_locidx_t numSubs = 0;
  bfam_domain_get_subdomains(dom,subdom_match,subdom_tags,
      dom->numSubdomains,subs,&numSubs);
  for(int s = 0; s < numSubs;s++)
  {
    int rval = bfam_dictionary_insert_ptr(&ts->elems,subs[s]->name,subs[s]);
    BFAM_ABORT_IF_NOT(rval != 1, "Issue adding subdomain %s", subs[s]->name);

    for(int n = 0; n < ts->nStages; n++)
    {
      char aux_rates_name[BFAM_BUFSIZ];
      snprintf(aux_rates_name,BFAM_BUFSIZ,"%s%d_",BFAM_LOCAL_ADAMS_PREFIX,n);
      aux_rates(subs[s],aux_rates_name);
      glue_rates(subs[s],aux_rates_name);
    }
  }

  numSubs = 0;
  bfam_domain_get_subdomains(dom,comm_match,comm_tags,
      dom->numSubdomains,subs,&numSubs);

  /* this tracks the max level that I know about */
#ifdef BFAM_DEBUG
  int local_max_levels = 0;
#endif

  /* find the plus and minus levels for the glue grids */
  for(int s = 0; s < numSubs;s++)
  {
    bfam_locidx_t m_lvl = -1;
    bfam_locidx_t p_lvl = -1;

    BFAM_ASSERT(subs[s]->glue_m);
    BFAM_ASSERT(subs[s]->glue_p);

    bfam_critbit0_allprefixed(&subs[s]->glue_p->tags,
        BFAM_LOCAL_ADAMS_LVL_PREFIX, get_tag_level_number,&p_lvl);
    bfam_critbit0_allprefixed(&subs[s]->glue_m->tags,
        BFAM_LOCAL_ADAMS_LVL_PREFIX, get_tag_level_number,&m_lvl);

    BFAM_ASSERT(m_lvl >= 0);
    BFAM_ASSERT(p_lvl >= 0);

#ifdef BFAM_DEBUG
    /* just a sanity check since these should match */
    BFAM_ASSERT(subs[s]->glue_m->sub_m);
    bfam_locidx_t s_lvl = -1;
    bfam_critbit0_allprefixed(&subs[s]->glue_m->sub_m->tags,
        BFAM_LOCAL_ADAMS_LVL_PREFIX, get_tag_level_number,&s_lvl);
    BFAM_ASSERT(s_lvl == m_lvl);
#endif

    /* communicate as infrequently as needed, so use the higher level */
    char comm_lvl_tag[BFAM_BUFSIZ];
    bfam_ts_local_adams_fill_comm_level_tag(comm_lvl_tag,BFAM_BUFSIZ,
        BFAM_MAX(m_lvl, p_lvl));
    bfam_subdomain_add_tag(subs[s],comm_lvl_tag);

#ifdef BFAM_DEBUG
    /* update the max level */
    local_max_levels = BFAM_MAX(local_max_levels,BFAM_MAX(m_lvl,p_lvl));
#endif
  }


#ifdef BFAM_DEBUG
  BFAM_ASSERT(local_max_levels <= num_lvl);
#endif

  /* loop through all possible commmunication tags */
  char *local_comm_tags[ts->numLevels+1];
  char tag_stor[BFAM_BUFSIZ*ts->numLevels];
  ts->comm_array = bfam_malloc(ts->numLevels*sizeof(bfam_ts_local_adams_t*));

  /* since we will use these make sure they are NULL */
  BFAM_ASSERT(comm_data->user_data == NULL);
  BFAM_ASSERT(comm_data->user_prefix_function == NULL);
  BFAM_ASSERT(comm_data->user_comm_info == NULL);
  BFAM_ASSERT(comm_data->user_put_send_buffer == NULL);
  BFAM_ASSERT(comm_data->user_get_recv_buffer == NULL);

  for(int k=0; k < ts->numLevels; k++)
  {
    local_comm_tags[k] = &tag_stor[BFAM_BUFSIZ*k];
    bfam_ts_local_adams_fill_comm_level_tag(&tag_stor[BFAM_BUFSIZ*k],
        BFAM_BUFSIZ, k);
    local_comm_tags[k+1] = NULL;

    /*
     * Set up the communicator we will use
     */
    ts->comm_array[k] = bfam_communicator_new(dom, BFAM_DOMAIN_OR,
        (const char**)local_comm_tags, mpicomm, mpitag, comm_data);
  }
}

void
bfam_ts_local_adams_free(bfam_ts_local_adams_t* ts)
{
  BFAM_LDEBUG("LOCAL ADAMS FREE");
  if(ts->lsrk != NULL)
  {
    bfam_ts_lsrk_free(ts->lsrk);
    bfam_free(ts->lsrk);
  }
  ts->lsrk = NULL;
  for(int k = 0; k < ts->numLevels; k++)
  {
    bfam_communicator_free(ts->comm_array[k]);
    bfam_free(ts->comm_array[k]);
    ts->comm_array[k] = NULL;
    if(ts->lvlStepArray)
      bfam_free(ts->lvlStepArray[k]);
  }
  bfam_free(ts->comm_array);
  bfam_dictionary_clear(&ts->elems);
  if(ts->lvlStepArray)
    bfam_free(ts->lvlStepArray);
  /*
  bfam_free_aligned(ts->A);
  ts->A = NULL;
  */
  ts->nStages = 0;
  bfam_free(ts->numStepsArray);
  bfam_free(ts->currentStageArray);
  ts->t  = NAN;
  bfam_ts_free(&ts->base);
}
