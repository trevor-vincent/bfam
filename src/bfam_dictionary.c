#include <bfam.h>
#include <bfam_dictionary.h>

void bfam_dictionary_init(bfam_dictionary_t *d)
{
  d->num_entries = 0;
  d->t.root = NULL;
}

void bfam_dictionary_clear(bfam_dictionary_t *d)
{
  bfam_critbit0_clear(&(d->t));
  d->num_entries = 0;
}

int bfam_dictionary_contains_check(const char * value, void * arg)
{
  (*(int*)arg)++;
  return 1;
}

int bfam_dictionary_contains(bfam_dictionary_t *d, const char *key)
{
  const size_t keylen = strlen(key);

  char *u = bfam_malloc(sizeof(char)*(keylen+2));

  memcpy(u, key, keylen);
  u[keylen] = BFAM_KEYVALUE_SPLIT;
  u[keylen+1] = '\0';
  int found = 0;
  bfam_critbit0_allprefixed(&(d->t),u,
      &bfam_dictionary_contains_check,&found);
  bfam_free(u);
  return found;
}

int bfam_dictionary_insert(bfam_dictionary_t *d, const char *key,
                                  const char *val)
{
  const size_t keylen = strlen(key);
  const size_t vallen = strlen(val);

  char *keyval = bfam_malloc(sizeof(char)*(keylen+vallen+2));

  memcpy(keyval, key, keylen);

  if(bfam_dictionary_contains(d,key))
  {
    free(keyval);
    return 1;
  }

  keyval[keylen]   = BFAM_KEYVALUE_SPLIT;
  memcpy(&keyval[keylen+1], val, vallen);
  keyval[keylen+vallen+1] = '\0';

  int rval = bfam_critbit0_insert(&(d->t), keyval);

  if(rval == 2)
   ++d->num_entries;

  bfam_free(keyval);
  return rval;
}

int bfam_dictionary_insert_ptr(bfam_dictionary_t *d, const char *key,
                                  const void *val_ptr)
{
  char val_str[BFAM_PTR_STR_LEN+1];
  snprintf(val_str,BFAM_PTR_STR_LEN+1,"%p",val_ptr);
  return bfam_dictionary_insert(d,key,val_str);
}

int bfam_dictionary_insert_locidx(bfam_dictionary_t *d, const char *key,
                                  const bfam_locidx_t val)
{
  char val_str[BFAM_BUFSIZ];
  snprintf(val_str,BFAM_BUFSIZ,"%"BFAM_LOCIDX_PRId,val);
  return bfam_dictionary_insert(d,key,val_str);
}

int bfam_dictionary_get_value_handle(const char * keyval, void * arg)
{
  char * key = (char  *)((void **)arg)[0];
  char** val = (char **)((void **)arg)[1];
  int keylen = strlen(key);

  *val = (char*)&keyval[keylen];

  return 1;
}


char* bfam_dictionary_get_value(bfam_dictionary_t *d, const char *key)
{
  if(!bfam_dictionary_contains(d,key))
    return NULL;

  const size_t keylen = strlen(key);

  char *u = bfam_malloc(sizeof(char)*(keylen+2));

  memcpy(u, key, keylen);
  u[keylen] = BFAM_KEYVALUE_SPLIT;
  u[keylen+1] = '\0';

  char *value = NULL;
  void *arg[2];
  arg[0] = u;
  arg[1] = &value;

  bfam_critbit0_allprefixed(&(d->t),u,
      &bfam_dictionary_get_value_handle,arg);

  bfam_free(u);

  return value;
}

void* bfam_dictionary_get_value_ptr(bfam_dictionary_t *d, const char *key)
{
  char* val_str = bfam_dictionary_get_value(d,key);
  if(val_str == NULL) return NULL;
  void *val_ptr = NULL;
  sscanf(val_str,"%p",&val_ptr);
  return val_ptr;
}

int bfam_dictionary_get_value_locidx(bfam_dictionary_t *d, const char *key,
    bfam_locidx_t *val)
{
  char* val_str = bfam_dictionary_get_value(d,key);
  if(val_str == NULL)
    return 0;
  int n = sscanf(val_str,"%"BFAM_LOCIDX_SCNd,val);
  return n;
}


typedef struct {
  int (*handle) (const char *, const char *, void *);
  void *arg;
} bfam_dict_allprex;

int
bfam_dictionary_allprefixed_usercall(const char * keyval, void * arg)
{
  char * key = (char *) keyval;
  char* split = strchr(key,BFAM_KEYVALUE_SPLIT);
  *split = '\0';
  bfam_dict_allprex *s_arg = (bfam_dict_allprex *)arg;
  s_arg->handle(key,split+1,s_arg->arg);
  *split = BFAM_KEYVALUE_SPLIT;
  return 1;
}

int bfam_dictionary_allprefixed(bfam_dictionary_t *d, const char *prefix,
                              int (*handle) (const char *, const char*, void *),
                              void *arg)
{
  bfam_dict_allprex args = {0,0};
  args.handle = handle;
  args.arg = arg;
  return bfam_critbit0_allprefixed(&(d->t),prefix,
      &bfam_dictionary_allprefixed_usercall,&args);
}

typedef struct {
  int (*handle) (const char *, void *, void *);
  void *arg;
} bfam_dict_allprex_ptr;

int
bfam_dictionary_allprefixed_usercall_ptr(const char * keyval, void * arg)
{
  char * key = (char *) keyval;
  char* split = strchr(key,BFAM_KEYVALUE_SPLIT);
  *split = '\0';

  void *val_ptr = NULL;
  sscanf(split+1,"%p",&val_ptr);

  bfam_dict_allprex_ptr *s_arg = (bfam_dict_allprex_ptr *)arg;
  s_arg->handle(key,val_ptr,s_arg->arg);

  *split = BFAM_KEYVALUE_SPLIT;
  return 1;
}


int bfam_dictionary_allprefixed_ptr(bfam_dictionary_t *d, const char *prefix,
                              int (*handle) (const char *, void*, void *),
                              void *arg)
{
  bfam_dict_allprex_ptr args = {0,0};
  args.handle = handle;
  args.arg = arg;
  return bfam_critbit0_allprefixed(&(d->t),prefix,
      &bfam_dictionary_allprefixed_usercall_ptr,&args);
}

int bfam_dictionary_swap(bfam_dictionary_t *d, const char *key1,
                                               const char *key2)
{
  char* t_val1 = bfam_dictionary_get_value(d,key1);
  if(!t_val1) return 1;

  /* return done if the keys are the same */
  if(0==strcmp(key1,key2)) return 3;

  char* t_val2 = bfam_dictionary_get_value(d,key2);
  if(!t_val2) return 2;


  /* Delete the first key-value pair */
  const size_t k1len = strlen(key1);
  const size_t v1len = strlen(t_val1);
  char *kv1 = bfam_malloc(sizeof(char)*(k1len+v1len+2));
  memcpy(kv1, key1, k1len);
  kv1[k1len]   = BFAM_KEYVALUE_SPLIT;
  memcpy(&kv1[k1len+1], t_val1, v1len);
  kv1[k1len+v1len+1] = '\0';
  bfam_critbit0_delete(&(d->t), kv1);

  /* Delete the second key-value pair */
  const size_t k2len = strlen(key2);
  const size_t v2len = strlen(t_val2);
  char *kv2 = bfam_malloc(sizeof(char)*(k2len+v2len+2));
  memcpy(kv2, key2, k2len);
  kv2[k2len]   = BFAM_KEYVALUE_SPLIT;
  memcpy(&kv2[k2len+1], t_val2, v2len);
  kv2[k2len+v2len+1] = '\0';
  bfam_critbit0_delete(&(d->t), kv2);

  /* re-insert the swapped values */
  int rval = bfam_dictionary_insert(d,key1,&kv2[k2len+1]);
  BFAM_ASSERT(rval != 1);
  bfam_free(kv2);
  kv2 = NULL;
  if(rval == 0)
  {
    bfam_free(kv1);
    kv1 = NULL;
    return rval;
  }

  rval = bfam_dictionary_insert(d,key2,&kv1[k1len+1]);
  BFAM_ASSERT(rval != 1);
  bfam_free(kv1);
  kv1 = NULL;
  if(rval == 0) return rval;

  return 3;
}
