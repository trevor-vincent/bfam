#include <bfam.h>
#include <bfam_subdomain_dummy.h>

void scale_rates(bfam_subdomain_t **subdomains_,void **dq_, const
    bfam_real_t a)
{
  bfam_subdomain_dummy_t **subdomains = (bfam_subdomain_dummy_t **) subdomains_;
  for(bfam_locidx_t s=0;subdomains[s];s++)
  {
    bfam_real_t *dq = (bfam_real_t *) dq_[s];
    dq[0] *= a;
  }
}

void update_rates(bfam_subdomain_t **subdomains_,void **dq_, const void **q_,
    const bfam_real_t t)
{
  bfam_subdomain_dummy_t **subdomains = (bfam_subdomain_dummy_t **) subdomains_;
  for(bfam_locidx_t s=0;subdomains[s];s++)
  {
    bfam_subdomain_dummy_t *sub = subdomains[s];
    bfam_real_t *dq = (bfam_real_t *) dq_[s];
    dq[0] += sub->N*pow(t,sub->N-1.0);
  }
}

void scale_add_rates(bfam_subdomain_t **subdomains_,void **q_, const void **dq_,
    const bfam_real_t b)
{
  bfam_subdomain_dummy_t **subdomains = (bfam_subdomain_dummy_t **) subdomains_;
  for(bfam_locidx_t s=0;subdomains[s];s++)
  {
    const bfam_real_t *dq = (const bfam_real_t *) dq_[s];
    bfam_real_t *q  = (bfam_real_t *)  q_[s];
    q[0] += b*dq[0];
  }
}


void
test_integration(int N, bfam_ts_lsrk_method_t method)
{
  bfam_ts_lsrk_t *ts;
  bfam_domain_t *dom = bfam_domain_new(NULL);
  bfam_subdomain_dummy_t *subDom = bfam_subdomain_dummy_new("a",N);
  bfam_domain_add_subdomain(dom, (bfam_subdomain_t*) subDom);

  const char *tags[]   = {"_subdomain_dummy", NULL};
  const char *fields[] = {"q",NULL};

  bfam_domain_add_fields(dom,BFAM_DOMAIN_AND,tags,fields);

  ts = bfam_ts_lsrk_new(dom,method, tags,BFAM_DOMAIN_AND,fields, &scale_rates,
      &update_rates, &scale_add_rates);

  bfam_ts_lsrk_step(ts,1.0);

  bfam_ts_lsrk_free(ts);
  bfam_free(ts);
  bfam_domain_free(dom);
  bfam_free(dom);
}

int
main (int argc, char *argv[])
{

  test_integration(4,BFAM_TS_LSRK_KC54);
  test_integration(3,BFAM_TS_LSRK_W33 );
  test_integration(2,BFAM_TS_LSRK_HEUN);
  test_integration(1,BFAM_TS_LSRK_FE  );

  return EXIT_SUCCESS;
}
