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

/* For LPs (and IPs and MIPs), objconst(n) is the constant term
 * for objective n (first objective has n = 0).
 */
#include "asl.h"

 real
objconst_ASL(ASL *asl, int n)
{
	real *oc;
	static char who[] = "objconst";

	if (!asl)
		badasl_ASL(asl,0,who);
	else if (asl->i.ASLtype < ASL_read_f || asl->i.ASLtype > ASL_read_pfgh)
		badasl_ASL(asl,ASL_read_f,who);

	if (n >= 0 && n < n_obj && (oc = asl->i.objconst))
		return oc[n];
	return 0.;
	}
