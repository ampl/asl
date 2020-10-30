/*******************************************************************
Copyright (C) 2016 AMPL Optimization, Inc.; written by David M. Gay.

Permission to use, copy, modify, and distribute this software and its
documentation for any purpose and without fee is hereby granted,
provided that the above copyright notice appear in all copies and that
both that the copyright notice and this permission notice and warranty
disclaimer appear in supporting documentation.

The author and AMPL Optimization, Inc. disclaim all warranties with
regard to this software, including all implied warranties of
merchantability and fitness.  In no event shall the author be liable
for any special, indirect or consequential damages or any damages
whatsoever resulting from loss of use, data or profits, whether in an
action of contract, negligence or other tortious action, arising out
of or in connection with the use or performance of this software.
*******************************************************************/

#ifndef JACPDIM_H_included
#define JACPDIM_H_included
#ifndef ASL_PFGH_included
#include "asl_pfgh.h"
#endif

#define conpival conpival_ew_ASL
#define conpgrd  conpgrd_ew_ASL
#define conpval  conpval_ew_ASL
#define jacpval  jacpval_ew_ASL
#define lconpval lconpval_ew_ASL
#define objpgrd  objpgrd_ew_ASL
#define objpval  objpval_ew_ASL
#define hvpcomp  hvpcomp_ew_ASL
#define hvpcompd hvpcompd_ew_ASL
#define hvpcompde hvpcompde_ew_ASL
#define hvpcompe hvpcompe_ew_ASL
#define hvpcomps hvpcomps_ew_ASL
#define hvpcompse hvpcompse_ew_ASL
#define xp2known xp2known_ew_ASL

#ifdef __cplusplus
extern "C" {
#endif
 extern real conpival(EvalWorkspace*, int nc, real *X, fint *ne);
 extern void conpgrd(EvalWorkspace*, int nc, real *X, real *G, fint *nerror);
 extern void conpval(EvalWorkspace*, real *X, real *F, fint *nerror);
 extern real eval2_ASL(int *o, EvalWorkspace*);
 extern void jacpval(EvalWorkspace*, real *X, real *JAC, fint *nerror);
 extern int  lconpval(EvalWorkspace*, int nc, real *X, fint *ne);
 extern void objpgrd(EvalWorkspace*, int nobj, real *X, real *G, fint *nerror);
 extern real objpval(EvalWorkspace*, int nobj, real *X, fint *nerror);
 extern void hvpcomp(EvalWorkspace*, real *hv, real *p, int nobj, real *ow, real *y);
 extern void hvpcompd(EvalWorkspace*,real *hv, real *p, int co);
 extern void hvpcompde(EvalWorkspace*,real *hv, real *p, int co, fint*);
 extern void hvpcompe(EvalWorkspace*, real *hv, real *p, int nobj, real *ow, real *y, fint*);
 extern varno_t hvpcomps(EvalWorkspace*, real *hv, real *p, int co, varno_t nz, varno_t *z);
 extern varno_t hvpcompse(EvalWorkspace*, real *hv, real *p, int co, varno_t nz, varno_t *z, fint*);

 extern int xp2known(EvalWorkspace*, real*, fint*);
#ifdef __cplusplus
	}
#endif
#endif /* JACPDIM_H_included */
