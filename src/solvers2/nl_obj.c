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

/* nl_obj(n) = 1 if objective n is nonlinear, 0 otherwise. */
#define SKIP_NL2_DEFINES
#include "nlp.h"
#include "nlp2.h"
#include "asl_pfg.h"
#include "asl_pfgh.h"
#include "obj_adj.h"

enum {nOPRET = 0};

 int
#ifdef KR_headers
nl_obj_ASL(asl, n) ASL *asl; int n;
#else
nl_obj_ASL(ASL *asl, int n)
#endif
{
	Objrep *od, **pod;
	int *o;
	static char who[] = "nl_obj";

	if (!asl)
		badasl_ASL(asl,0,who);
	else if (asl->i.ASLtype < ASL_read_f || asl->i.ASLtype >ASL_read_pfgh)
		badasl_ASL(asl,ASL_read_f,who);

	if (n >= 0 && n < n_obj) {
		if ((pod = asl->i.Or) && (od = pod[n])) {
			n = od->ico;
			switch(asl->i.ASLtype) {
			  case ASL_read_fgh:
				o = (((ASL_fgh*)asl)->I.con2_de_ + n)->o.e;
				break;
			  case ASL_read_pfg:
				o = (((ASL_pfg*)asl)->I.con_de_ + n)->o.e;
				break;
			  case ASL_read_pfgh:
				o = (((ASL_pfgh*)asl)->I.con2_de_ + n)->o.e;
				break;
			  default:
				o = (((ASL_fg*)asl)->I.con_de_ + n)->o.e;
			  }
			}
		else {
			switch(asl->i.ASLtype) {
			  case ASL_read_fgh:
				o = (((ASL_fgh*)asl)->I.obj2_de_ + n)->o.e;
				break;
			  case ASL_read_pfg:
				o = (((ASL_pfg*)asl)->I.obj_de_ + n)->o.e;
				break;
			  case ASL_read_pfgh:
				o = (((ASL_pfgh*)asl)->I.obj2_de_ + n)->o.e;
				break;
			  default:
				o = (((ASL_fg*)asl)->I.obj_de_ + n)->o.e;
			  }
			}
		if (o && (*o != nOPRET || o[1] >= 0))
			return 1;
		}
	return 0;
	}
