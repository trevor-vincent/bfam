#ifndef BFAM_TIMESTEP_LSRK_H
#define BFAM_TIMESTEP_LSRK_H

#include <bfam_base.h>
#include <bfam_domain.h>
#include <bfam_timestep.h>

/**
 * structure comtaining the necessary features of a low memory RK scheme
 *
 * Each stages the form
 * t  := t0       + c*dt
 * dq := RHS(q,t) + a*dq
 * q  := q        + b*dt*dq
 */
typedef struct bfam_ts_lsrk
{
  bfam_ts_t base;       /**< parent timestepper */
  bfam_long_real_t* A;  /**< low memory RK A: rate scale */
  bfam_long_real_t* B;  /**< low memory RK B: update scale */
  bfam_long_real_t* C;  /**< low memory RK C: time scale*/
  int nStages;          /**< number of stages */
  bfam_long_real_t  t;  /**< domain time */
  bfam_long_real_t  dt; /**< domain dt   */
  bfam_subdomain_t ** subdomains; /**< subdomains to handle */
  void             ** fields;     /**< pointers to the fields I handle */
  void             ** rates;      /**< pointers to the my rates */
  bfam_locidx_t    numSubdomains; /**< number of subdomains I handle */
  bfam_locidx_t    numFields;     /**< number of fields I handle */
  /**< Function pointer for scale rates. s is a null terminated array of
   * subdomains. dq := a*dq */
  void (*scale_rates) (bfam_subdomain_t** s,void** dq, const bfam_real_t a);
  /**< Function pointer for update rates. s is a null terminated array of
   * subdomains. dq := dq + RHS(q,t) */
  void (*update_rates) (bfam_subdomain_t** s,void** dq,
      const void** q, const bfam_real_t t);
  /**< Function pointer for scale and add rates. s is a null terminated array of
   * subdomains. q := q + b*dq */
  void (*scale_add_rates) (bfam_subdomain_t** s,void** q,
      const void** dq, const bfam_real_t b);
} bfam_ts_lsrk_t;

typedef enum bfam_ts_lsrk_method
{
  BFAM_TS_LSRK_KC54,
  BFAM_TS_LSRK_FE,
  BFAM_TS_LSRK_HEUN,
  BFAM_TS_LSRK_W33
} bfam_ts_lsrk_method_t;

/** create a low storage RK scheme
 *
 * \warning It is the callers responsibility to ensure that
 *          \a dom is freed after this LSRK is
 *
 * \param [in]  dom              pointer to the domain
 * \param [in]  method           Low storage RK type we are using
 * \param [in]  tags             Tags for the subdomains to be updated
 * \param [in]  match            type of match, \c BFAM_DOMAIN_OR will
 *                               match subdomains with any of the tags
 *                               and \c BFAM_DOMAIN_AND will match subdomains
 *                               with all of the tags.
 * \param [in]  fields           Fields I am in charge of updating
 * \param [in]  scale_rates      function pointer to scale rates function
 *                               (dq := a * dq)
 * \param [in]  update_rates     function pointer to update rates function
 *                               (dq := dq + RHS(q,t))
 * \param [in]  scale_add_rates  function pointer to update rates function
 *                               (q  := q + b*dq)
 *
 * \return the newly created low storage RK time stepper
 */
bfam_ts_lsrk_t*
bfam_ts_lsrk_new(bfam_domain_t* dom, bfam_ts_lsrk_method_t method,
    const char* tags[], bfam_domain_match_t match,
    const char* fields[],
    void (*scale_rates) (bfam_subdomain_t**,void**,const bfam_real_t),
    void (*update_rates) (bfam_subdomain_t**,void**,const void**,
      const bfam_real_t),
    void (*scale_add_rates) (bfam_subdomain_t**,void**,
      const void**, const bfam_real_t));

/** initialize a low storage RK scheme
 *
 * \warning It is the callers responsibility to ensure that
 *          \a dom is freed after this LSRK is
 *
 * \param [in,out]  ts               pointer to time stepper to initialize
 * \param [in]      dom              pointer to the domain
 * \param [in]      method           Low storage RK type we are using
 * \param [in]      tags             Tags for the subdomains to be updated
 * \param [in]      match            type of match, \c BFAM_DOMAIN_OR will
 *                                   match subdomains with any of the tags
 *                                   and \c BFAM_DOMAIN_AND will match
 *                                   subdomains with all of the tags.
 * \param [in]      fields           Fields I am in charge of updating
 * \param [in]      scale_rates      function pointer to scale rates function
 *                                   (dq := a * dq)
 * \param [in]      update_rates     function pointer to update rates function
 *                                   (dq := dq + RHS(q,t))
 * \param [in]      scale_add_rates  function pointer to update rates function
 *                                   (q  := q + b*dq)
 */
void
bfam_ts_lsrk_init(bfam_ts_lsrk_t* ts, bfam_domain_t* dom,
    bfam_ts_lsrk_method_t method,
    const char* tags[], bfam_domain_match_t match,
    const char* fields[],
    void (*scale_rates) (bfam_subdomain_t**,void**,const bfam_real_t),
    void (*update_rates) (bfam_subdomain_t**,void**,const void**,
      const bfam_real_t),
    void (*scale_add_rates) (bfam_subdomain_t**,void**,
      const void**, const bfam_real_t));

/** free a low storage RK scheme
 *
 * \param [in,out]  ts       pointer to time stepper to free
 */
void
bfam_ts_lsrk_free(bfam_ts_lsrk_t* ts);

/** set the time of the scheme
 *
 * \param [in,out]  ts       pointer to time stepper to set
 * \param [in]      time     time to set
 */
void
bfam_ts_lsrk_set_time(bfam_ts_lsrk_t* ts,bfam_long_real_t time);

/** set dt of the scheme
 *
 * \param [in,out]  ts       pointer to time stepper to set
 * \param [in]      dt       dt to set
 */
void
bfam_ts_lsrk_set_dt(bfam_ts_lsrk_t* ts,bfam_long_real_t dt);

/** get the time of the scheme
 *
 * \param [in]  ts       pointer to lsrk to get time
 */
bfam_long_real_t
bfam_ts_lsrk_get_time(bfam_ts_lsrk_t* ts);

/** get dt of the scheme
 *
 * \param [in]  ts       pointer to lsrk to get dt
 */
bfam_long_real_t
bfam_ts_lsrk_get_dt(bfam_ts_lsrk_t* ts);

#endif
