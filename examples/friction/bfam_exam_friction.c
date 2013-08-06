#include <bfam.h>

#define REAL_APPROX_EQ(x, y, K)                                              \
  BFAM_APPROX_EQ((x), (y), (K), BFAM_REAL_ABS, BFAM_REAL_EPS, BFAM_REAL_EPS)

static int
bfam_bracketed_newton(bfam_real_t *sol, bfam_real_t x, bfam_real_t xmin,
    bfam_real_t xmax,
    void (*friction_law)(bfam_real_t *fval, bfam_real_t *dval,
      bfam_real_t val, void* args),
    void* args)
{
  const int MAX_ITERATIONS = 200;

  const bfam_real_t ftol = 10*BFAM_REAL_EPS;

  bfam_real_t f    = 0;
  bfam_real_t df = 0;

  bfam_real_t fmin = 0;
  friction_law(&fmin,&df,x,args);

  int rval = 1;

  bfam_real_t dx = 0;
  for(int i = 0; i < MAX_ITERATIONS; ++i)
  {
    /* evaluate function and derivative */
    friction_law(&f,&df,x,args);

    /* check to see if we are converged */
    if(BFAM_REAL_ABS(x) < ftol)
    {
      rval = 0;
      break;
    }

    /* update the bracket */
    if(f*fmin > 0) xmin =x;
    else xmax = x;

    /* update the solution */
    dx = -f/df;
    x += dx;

    /* check that solution is in the bracket */
    if(x<xmin || x > xmax)
    {
      x = (xmin+xmax)/2.0;
      dx = xmax-xmin;
    }

  }

  *sol = x;
  return rval;
}

typedef struct prefs
{
  lua_State *L;

  bfam_locidx_t num_subdomains;
} prefs_t;

typedef struct exam
{
  MPI_Comm mpicomm;
  int      mpirank;
  int      mpisize;

  bfam_domain_t  *domain;
} exam_t;

static void
init_mpi(exam_t *exam, MPI_Comm mpicomm)
{
  exam->mpicomm = mpicomm;
  BFAM_MPI_CHECK(MPI_Comm_rank(mpicomm, &exam->mpirank));
  BFAM_MPI_CHECK(MPI_Comm_size(mpicomm, &exam->mpisize));
}

static void
init_domain(exam_t *exam, prefs_t *prefs)
{
  exam->domain = bfam_domain_new(exam->mpicomm);
}

static void
free_exam(exam_t *exam)
{
  bfam_domain_free(exam->domain);
  bfam_free(exam->domain);
}

static void
run(MPI_Comm mpicomm, prefs_t *prefs)
{
  exam_t exam;

  init_mpi(&exam, mpicomm);

  init_domain(&exam, prefs);

  free_exam(&exam);
}

static int
get_global_int(lua_State *L, const char *name, int def)
{
  lua_getglobal(L, name);
  int result = def;
  if(!lua_isnumber(L, -1))
    BFAM_ROOT_WARNING("`%s' not found, using default", name);
  else
    result = (int)lua_tonumber(L, -1);
  lua_pop(L, 1);
  return result;
}

static prefs_t *
new_prefs(const char *prefs_filename)
{
  prefs_t *prefs = bfam_malloc(sizeof(prefs_t));

  lua_State *L = luaL_newstate();
  luaL_openlibs(L);

  if(luaL_loadfile(L, prefs_filename) || lua_pcall(L, 0, 0, 0))
    BFAM_LERROR("cannot run configuration file: `%s'", lua_tostring(L, -1));

  prefs->L = L;

  BFAM_ASSERT(lua_gettop(L)==0);

  prefs->num_subdomains = get_global_int(L, "num_subdomains", 1);
  BFAM_ASSERT(lua_gettop(L)==0);

  lua_pop(L, 1);

  BFAM_ASSERT(lua_gettop(L)==0);

  return prefs;
}

static void
print_prefs(prefs_t *prefs)
{
  BFAM_ROOT_INFO("----------Preferences----------");
  BFAM_ROOT_INFO("num_subdomains=%"BFAM_LOCIDX_PRId, prefs->num_subdomains);
  BFAM_ROOT_INFO("-------------------------------");
}

static void
free_prefs(prefs_t *prefs)
{
  lua_close(prefs->L);
}

int
main(int argc, char *argv[])
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;

  void *options= bfam_gopt_sort(&argc, (const char**)argv,
      bfam_gopt_start(
        bfam_gopt_option('h', 0,
                         bfam_gopt_shorts('h', '?'),
                         bfam_gopt_longs("help", "HELP")),
        bfam_gopt_option('V', 0,
                         bfam_gopt_shorts('V'),
                         bfam_gopt_longs("version")),
        bfam_gopt_option('v', BFAM_GOPT_REPEAT,
                         bfam_gopt_shorts('v'),
                         bfam_gopt_longs("verbose"))
        )
      );

  const char *help_text =
  "  %s [options] prefs_file\n"
  "\n"
  "  there are four possible options to this program, some of which have\n"
  "  multiple names:\n"
  "\n"
  "    -h -? --help --HELP\n"
  "    -V --version\n"
  "    -v --verbose  (which may be repeated for more verbosity)\n"
  "\n";

  if(bfam_gopt(options, 'h'))
  {
    /*
     * if any of the help options was specified
     */
    BFAM_ROOT_INFO(help_text, argv[0]);
    exit(EXIT_SUCCESS);
  }

  if(bfam_gopt(options, 'V'))
  {
    BFAM_ROOT_INFO("BFAM Version: %s", bfam_version_get());
    BFAM_ROOT_INFO("BFAM Compile Info:\n" BFAM_COMPILE_INFO);
    exit( EXIT_SUCCESS );
  }

  int verbosity = bfam_gopt(options, 'v');

  BFAM_MPI_CHECK(MPI_Init(&argc,&argv));
  BFAM_MPI_CHECK(MPI_Comm_rank(comm, &rank));

  if(argc != 2)
  {
    BFAM_LERROR("Unexpected number of arguments.");
    BFAM_ROOT_INFO(help_text, argv[0]);
    exit(EXIT_FAILURE);
  }

  int logLevel = BFAM_MAX(BFAM_LL_INFO - verbosity, BFAM_LL_ALWAYS);

  bfam_log_init(rank, stdout, logLevel);
  bfam_signal_handler_set();

  prefs_t *prefs = new_prefs(argv[1]);
  print_prefs(prefs);

  run(comm, prefs);
  free_prefs(prefs);
  bfam_free(prefs);

  sc_finalize();
  BFAM_MPI_CHECK(MPI_Finalize());

  bfam_gopt_free(options);

  return EXIT_SUCCESS;
}
