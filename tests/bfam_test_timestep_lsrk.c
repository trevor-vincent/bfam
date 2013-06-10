#include <bfam.h>
#include <bfam_subdomain_dummy.h>

int
main (int argc, char *argv[])
{
  bfam_ts_lsrk_t* ts;
  bfam_domain_t* dom = bfam_domain_new(NULL);
  bfam_subdomain_dummy_t* subDom = bfam_subdomain_dummy_new("a",4);
  bfam_domain_add_subdomain(dom, (bfam_subdomain_t*) subDom);

  const char *tags[]   = {"_subdomain_dummy", NULL};
  const char *fields[] = {"q",NULL};

  // printf("%e\n",(double) bfam_subdomain_dummy_exact(subDom, 1.0));
  // printf("%e\n",(double) bfam_subdomain_dummy_exact(subDom, 0.2));

  /* test Kennedy-Carpenter(5,4) */
  ts = bfam_ts_lsrk_new(dom,BFAM_TS_LSRK_KC54,
      tags,BFAM_DOMAIN_AND,fields, NULL, NULL, NULL);
  bfam_ts_lsrk_free(ts);
  bfam_free(ts);

  /* test Williamson(3,3) */
  ts = bfam_ts_lsrk_new(dom,BFAM_TS_LSRK_W33,
      tags,BFAM_DOMAIN_AND,fields, NULL, NULL, NULL);
  bfam_ts_lsrk_free(ts);
  bfam_free(ts);

  /* test Heun */
  ts = bfam_ts_lsrk_new(dom,BFAM_TS_LSRK_HEUN,
      tags,BFAM_DOMAIN_AND,fields, NULL, NULL, NULL);
  bfam_ts_lsrk_free(ts);
  bfam_free(ts);

  /* test Forward Euler */
  ts = bfam_ts_lsrk_new(dom,BFAM_TS_LSRK_FE,
      tags,BFAM_DOMAIN_AND,fields, NULL, NULL, NULL);
  bfam_ts_lsrk_free(ts);
  bfam_free(ts);

  bfam_domain_free(dom);
  bfam_free(dom);

  return EXIT_SUCCESS;
}
