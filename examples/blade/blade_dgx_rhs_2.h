#ifndef BLADE_DGX_RHS_2_H
#define BLADE_DGX_RHS_2_H

#include <bfam.h>

#define X(order)                                                     \
  void blade_dgx_scale_rates_advection_2_##order(int N,              \
      bfam_subdomain_dgx_t *sub, const char *rate_prefix,            \
      const bfam_long_real_t a);                                     \
  void blade_dgx_intra_rhs_advection_2_##order(int N,                \
      bfam_subdomain_dgx_t *sub, const char *rate_prefix,            \
      const char *field_prefix, const bfam_long_real_t t);           \
  void blade_dgx_add_rates_advection_2_##order(int N,                \
      bfam_subdomain_dgx_t *sub, const char *field_prefix_lhs,       \
      const char *field_prefix_rhs, const char *rate_prefix,         \
      const bfam_long_real_t a);                                     \
  void blade_dgx_energy_2_##order(int N,                             \
      bfam_real_t* energy_sq, bfam_subdomain_dgx_t *sub,             \
      const char *field_prefix);
BFAM_LIST_OF_DGX_NORDERS
#undef X

void blade_dgx_scale_rates_advection_2_(int N,
    bfam_subdomain_dgx_t *sub, const char *rate_prefix,
    const bfam_long_real_t a);
void blade_dgx_intra_rhs_advection_2_(int N,
    bfam_subdomain_dgx_t *sub, const char *rate_prefix,
    const char *field_prefix, const bfam_long_real_t t);
void blade_dgx_add_rates_advection_2_(int N,
    bfam_subdomain_dgx_t *sub, const char *field_prefix_lhs,
    const char *field_prefix_rhs, const char *rate_prefix,
    const bfam_long_real_t a);
void blade_dgx_energy_2_(int N, bfam_real_t *energy_sq,
    bfam_subdomain_dgx_t *sub, const char *field_prefix);

#endif
