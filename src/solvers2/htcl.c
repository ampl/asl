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

#include "asl.h"

size_t mblk_mem_ASL;

 int
htcl_ASL(size_t x)
{
	int k = 0;
	size_t L = sizeof(void *);

	while(L < x) {
		++k;
		if (!(L <<= 1))
			memfailure_ASL("htcl", "arg too big", x);
		}
	return k;
	}

 void *
new_mblk_ASL(ASL *asl, uint k)
{
	MBavail *a;
	MBFree *mb;
	size_t L;
	uint na;

	if (k >= MBLK_KMAX_ASL)
		memfailure_ASL("new_mblk", "arg too big", k);
	mb = &asl->i.mblk_free[k];
	ACQUIRE_MBLK_LOCK(&asl->i, mb->Lock);
	na = ++mb->nalloc;
	if ((a = mb->a))
		mb->a = a->mbnext;
	FREE_MBLK_LOCK(&asl->i, mb->Lock);
	if (!a) {
		mblk_mem_ASL += L = sizeof(MBhead) + (sizeof(void*)<<k);
		a = (MBavail*)mem_ASL(asl, L);
		a->h.klass = k;
		}
	a->h.nalloc = na;
	return (void*)&a->mbnext;
	}

 void
Del_mblk_ASL(ASL *asl, void *x)
{
	MBavail *a;
	MBFree *mb;
	uint k;

	a = (MBavail*)((char*)x - sizeof(MBhead));
	if ((k = a->h.klass) >= MBLK_KMAX_ASL)
		memfailure_ASL("Del_mblk", "invalid klass", (size_t)k);
	mb = &asl->i.mblk_free[k];
	ACQUIRE_MBLK_LOCK(&asl->i, mb->Lock);
	a->mbnext = mb->a;
	mb->a = a;
	FREE_MBLK_LOCK(&asl->i, mb->Lock);
	}
