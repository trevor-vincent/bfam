#include <bfam.h>
#include <bfam_subdomain_dummy.h>

void
test_integration(int N, bfam_ts_lsrk_method_t method)
{
  bfam_ts_lsrk_t* ts;
  bfam_domain_t* dom = bfam_domain_new(NULL);
  bfam_subdomain_dummy_t* subDom = bfam_subdomain_dummy_new("a",N);
  bfam_domain_add_subdomain(dom, (bfam_subdomain_t*) subDom);

  const char *tags[]   = {"_subdomain_dummy", NULL};
  const char *fields[] = {"q",NULL};

  bfam_domain_add_fields(dom,BFAM_DOMAIN_AND,tags,fields);

  ts = bfam_ts_lsrk_new(dom,method,
      tags,BFAM_DOMAIN_AND,fields, NULL, NULL, NULL);

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
