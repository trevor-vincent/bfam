#include <bfam_subdomain_dummy.h>

bfam_subdomain_dummy_t*
bfam_subdomain_dummy_new(const char             *name,
                         const int               N)
{
  bfam_subdomain_dummy_t* newSubdomain =
    bfam_malloc(sizeof(bfam_subdomain_dummy_t));

  bfam_subdomain_dummy_init(newSubdomain, name, N);

  return newSubdomain;
}

static int
bfam_subdomain_dummy_field_add(bfam_subdomain_t *subdomain, const char *name)
{
  bfam_subdomain_dummy_t *s = (bfam_subdomain_dummy_t*) subdomain;

  if(bfam_dictionary_get_value_ptr(&s->base.fields,name))
    return 1;

  size_t fieldSize = sizeof(bfam_real_t);
  bfam_real_t *field = bfam_malloc_aligned(fieldSize);

  *field = 0;

  int rval = bfam_dictionary_insert_ptr(&s->base.fields, name, field);

  BFAM_ASSERT(rval != 1);

  if(rval == 0)
    bfam_free_aligned(field);

  return rval;
}

void
bfam_subdomain_dummy_init(bfam_subdomain_dummy_t       *subdomain,
                             const char                *name,
                             const int                  N)
{
  bfam_subdomain_init(&subdomain->base, name);
  bfam_subdomain_add_tag(&subdomain->base, "_subdomain_dummy");
  subdomain->base.free = bfam_subdomain_dummy_free;
  subdomain->N = N;
  subdomain->base.field_add = bfam_subdomain_dummy_field_add;
}

static int
bfam_subdomain_dummy_free_fields(const char * key, void *val, void *arg)
{
  bfam_free_aligned(val);

  return 1;
}


void
bfam_subdomain_dummy_free(bfam_subdomain_t *thisSubdomain)
{
  bfam_subdomain_dummy_t *sub = (bfam_subdomain_dummy_t*) thisSubdomain;
  bfam_dictionary_allprefixed_ptr(&sub->base.fields,"",
      &bfam_subdomain_dummy_free_fields,NULL);
  bfam_subdomain_free(thisSubdomain);
}

bfam_long_real_t
bfam_subdomain_dummy_exact(bfam_subdomain_dummy_t *subdom, bfam_long_real_t t)
{
  return pow(t,(bfam_long_real_t) subdom->N);
}
