#ifndef BFAM_TIMESTEP_ADAMS_H
#define BFAM_TIMESTEP_ADAMS_H

#include <bfam_base.h>
#include <bfam_domain.h>
#include <bfam_timestep.h>
#include <bfam_timestep_lsrk.h>
#include <bfam_communicator.h>

/**
 * structure comtaining the necessary features of an explicit Adams method
 *
 * Each step is of the form
 * t_{n+1}  := t_n + dt
 * q_{n+1}  := q_{n} + dt \sum_{k=0}^{m} a_{k} dq_{n-k}
 * dq_{n+1} := RHS(q_{n+1},t_{n+1})
 *
 */
typedef struct bfam_ts_adams
{
  bfam_ts_t base;            /**< parent timestepper */
  // bfam_long_real_t* A;       /**< coefficients */
  int nStages;               /**< number of steps */
  int currentStage;          /**< current stage counter */
  int numSteps;              /**< number of steps completed */
  bfam_long_real_t  t;       /**< domain time */
  bfam_communicator_t *comm; /**< communicator I handle */
  bfam_dictionary_t elems;   /**< dictionary of subdomains I step */

  /* LSRK method for initialization */
  bfam_ts_lsrk_t* lsrk;

  /* scale rates function */
  void (*scale_rates) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const bfam_long_real_t a);

  /* compute rhs that does not require communication */
  void (*intra_rhs) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const char *minus_rate_prefix,
      const char *field_prefix, const bfam_long_real_t t);

  /* compute rhs that does require communication */
  void (*inter_rhs) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const char *minus_rate_prefix,
      const char *field_prefix, const bfam_long_real_t t);

  /* add the rates to the fields: q_lhs := q_rhs + a*dq */
  /* NOTE: should handle case of in place addition */
  void (*add_rates) (bfam_subdomain_t *thisSubdomain,
      const char *field_prefix_lhs, const char *field_prefix_rhs,
      const char *rate_prefix, const bfam_long_real_t a);
} bfam_ts_adams_t;

typedef enum bfam_ts_adams_method
{
  BFAM_TS_ADAMS_1,
  BFAM_TS_ADAMS_2,
  BFAM_TS_ADAMS_3,
  BFAM_TS_ADAMS_4,
  BFAM_TS_ADAMS_NOOP,
} bfam_ts_adams_method_t;

/** create an explicit Adams scheme
 *
 * \warning It is the callers responsibility to ensure that
 *          \a dom is freed after this Adams is
 *
 * \param [in]  dom              pointer to the domain
 * \param [in]  method           Adams tyoe we are using
 * \param [in]  subdom_match     match type for subdomains
 * \param [in]  subdom_tags      tags for the subdomains to time step
 * \param [in]  comm_match       match type for communication
 * \param [in]  comm_tags        tags for the communication required for RHS
 * \param [in]  mpicomm          MPI_Comm to use for communication
 * \param [in]  mpitag           tag to use for MPI communcation
 * \param [in]  comm_data        user data passed to the communicator new
 * \param [in]  aux_rates        create rate field with given prefix
 * \param [in]  scale_rates      scale rates function
 * \param [in]  intra_rhs        function handle to intra RHS routine
 * \param [in]  inter_rhs        function handle to inter RHS routine
 * \param [in]  add_rates        function handle to add rates routine
 * \param [in]  RK_init          boolean which if true signifies using LSRK to
 *                               init Adams method
 *
 * \return the newly created low storage RK time stepper
 */
bfam_ts_adams_t*
bfam_ts_adams_new(bfam_domain_t* dom, bfam_ts_adams_method_t method,
    bfam_domain_match_t subdom_match, const char** subdom_tags,
    bfam_domain_match_t comm_match, const char** comm_tags,
    MPI_Comm mpicomm, int mpitag, void * comm_data,
    void (*aux_rates) (bfam_subdomain_t *thisSubdomain, const char *prefix),
    void (*scale_rates) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const bfam_long_real_t a),
    void (*intra_rhs) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const char*minus_rate_prefix,
      const char *field_prefix, const bfam_long_real_t t),
    void (*inter_rhs) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const char *minus_rate_prefix,
      const char *field_prefix, const bfam_long_real_t t),
    void (*add_rates) (bfam_subdomain_t *thisSubdomain,
      const char *field_prefix_lhs, const char *field_prefix_rhs,
      const char *rate_prefix, const bfam_long_real_t a),
    const int RK_init);

/** initialize an Adams scheme
 *
 * \warning It is the callers responsibility to ensure that
 *          \a dom is freed after this Adams is
 *
 * \param [in,out]  ts           pointer to time stepper to initialize
 * \param [in]  dom              pointer to the domain
 * \param [in]  method           Low storage RK tyoe we are using
 * \param [in]  subdom_match     match type for subdomains
 * \param [in]  subdom_tags      tags for the subdomains to time step
 * \param [in]  comm_match       match type for communication
 * \param [in]  comm_tags        tags for the communication required for RHS
 * \param [in]  mpicomm          MPI_Comm to use for communication
 * \param [in]  mpitag           tag to use for MPI communcation
 * \param [in]  comm_data        user data passed to the communicator new
 * \param [in]  aux_rates        create rate field with given prefix
 * \param [in]  scale_rates      scale rates function
 * \param [in]  intra_rhs        function handle to intra RHS routine
 * \param [in]  inter_rhs        function handle to inter RHS routine
 * \param [in]  add_rates        function handle to add rates routine
 * \param [in]  RK_init          boolean which if true signifies using LSRK to
 *                               init Adams method
 */
void
bfam_ts_adams_init(bfam_ts_adams_t* ts,
    bfam_domain_t* dom, bfam_ts_adams_method_t method,
    bfam_domain_match_t subdom_match, const char** subdom_tags,
    bfam_domain_match_t comm_match, const char** comm_tags,
    MPI_Comm mpicomm, int mpitag, void * comm_data,
    void (*aux_rates) (bfam_subdomain_t *thisSubdomain, const char *prefix),
    void (*scale_rates) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const bfam_long_real_t a),
    void (*intra_rhs) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const char *minus_rate_prefix,
      const char *field_prefix, const bfam_long_real_t t),
    void (*inter_rhs) (bfam_subdomain_t *thisSubdomain,
      const char *rate_prefix, const char *minus_rate_prefix,
      const char *field_prefix, const bfam_long_real_t t),
    void (*add_rates) (bfam_subdomain_t *thisSubdomain,
      const char *field_prefix_lhs, const char *field_prefix_rhs,
      const char *rate_prefix, const bfam_long_real_t a),
    const int RK_init);

/** free an Adams scheme
 *
 * \param [in,out]  ts       pointer to time stepper to free
 */
void
bfam_ts_adams_free(bfam_ts_adams_t* ts);

/** set the time of the scheme
 *
 * \param [in,out]  ts       pointer to time stepper to set
 * \param [in]      time     time to set
 */
void
bfam_ts_adams_set_time(bfam_ts_adams_t* ts, bfam_long_real_t time);


/** get the time of the scheme
 *
 * \param [in]  ts       pointer to adams to get time
 */
bfam_long_real_t
bfam_ts_adams_get_time(bfam_ts_adams_t* ts);

#endif
