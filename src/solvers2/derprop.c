/****************************************************************
Copyright (C) 2014 AMPL Optimization, Inc.

Permission to use, copy, modify, and distribute this software and its
documentation for any purpose and without fee is hereby granted,
provided that the above copyright notice appear in all copies and that
both that the copyright notice and this permission notice and warranty
disclaimer appear in supporting documentation.

AMPL Optimization disclaims all warranties with regard to this
software, including all implied warranties of merchantability and
fitness.  In no event shall AMPL Optimization or any of its entities
be liable for any special, indirect or consequential damages or any
damages whatsoever resulting from loss of use, data or profits,
whether in an action of contract, negligence or other tortious action,
arising out of or in connection with the use or performance of this
software.
****************************************************************/

#include "asl.h"

 void
derprop(derpblock *db, real *s, real *w, real f)
{
	derp *d, *de;	
	size_t n;

	d = db->d0;
	s[d->b] = f;
	for(;;) {
		for(de = db->de; d < de; ++d)
			s[d->a] += s[d->b] * w[d->c];
		if ((n = db->nxt))
			db = *(derpblock**)&w[n];
		else if (!(db = db->next))
			break;
		d = db->d0;
		}
	}
