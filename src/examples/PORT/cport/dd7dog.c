/* dd7dog.f -- translated by f2c (version 20160102).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Subroutine */ int dd7dog_(doublereal *dig, integer *lv, integer *n, 
	doublereal *nwtstp, doublereal *step, doublereal *v)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal t, t1, t2, cfact, relax, cnorm, gnorm, rlambd, ghinvg, 
	    femnsq, ctrnwt, nwtnrm;


/*  ***  COMPUTE DOUBLE DOGLEG STEP  *** */

/*  ***  PARAMETER DECLARATIONS  *** */


/*  ***  PURPOSE  *** */

/*        THIS SUBROUTINE COMPUTES A CANDIDATE STEP (FOR USE IN AN UNCON- */
/*     STRAINED MINIMIZATION CODE) BY THE DOUBLE DOGLEG ALGORITHM OF */
/*     DENNIS AND MEI (REF. 1), WHICH IS A VARIATION ON POWELL*S DOGLEG */
/*     SCHEME (REF. 2, P. 95). */

/* --------------------------  PARAMETER USAGE  -------------------------- */

/*    DIG (INPUT) DIAG(D)**-2 * G -- SEE ALGORITHM NOTES. */
/*      G (INPUT) THE CURRENT GRADIENT VECTOR. */
/*     LV (INPUT) LENGTH OF V. */
/*      N (INPUT) NUMBER OF COMPONENTS IN  DIG, G, NWTSTP,  AND  STEP. */
/* NWTSTP (INPUT) NEGATIVE NEWTON STEP -- SEE ALGORITHM NOTES. */
/*   STEP (OUTPUT) THE COMPUTED STEP. */
/*      V (I/O) VALUES ARRAY, THE FOLLOWING COMPONENTS OF WHICH ARE */
/*             USED HERE... */
/* V(BIAS)   (INPUT) BIAS FOR RELAXED NEWTON STEP, WHICH IS V(BIAS) OF */
/*             THE WAY FROM THE FULL NEWTON TO THE FULLY RELAXED NEWTON */
/*             STEP.  RECOMMENDED VALUE = 0.8 . */
/* V(DGNORM) (INPUT) 2-NORM OF DIAG(D)**-1 * G -- SEE ALGORITHM NOTES. */
/* V(DSTNRM) (OUTPUT) 2-NORM OF DIAG(D) * STEP, WHICH IS V(RADIUS) */
/*             UNLESS V(STPPAR) = 0 -- SEE ALGORITHM NOTES. */
/* V(DST0) (INPUT) 2-NORM OF DIAG(D) * NWTSTP -- SEE ALGORITHM NOTES. */
/* V(GRDFAC) (OUTPUT) THE COEFFICIENT OF  DIG  IN THE STEP RETURNED -- */
/*             STEP(I) = V(GRDFAC)*DIG(I) + V(NWTFAC)*NWTSTP(I). */
/* V(GTHG)   (INPUT) SQUARE-ROOT OF (DIG**T) * (HESSIAN) * DIG -- SEE */
/*             ALGORITHM NOTES. */
/* V(GTSTEP) (OUTPUT) INNER PRODUCT BETWEEN G AND STEP. */
/* V(NREDUC) (OUTPUT) FUNCTION REDUCTION PREDICTED FOR THE FULL NEWTON */
/*             STEP. */
/* V(NWTFAC) (OUTPUT) THE COEFFICIENT OF  NWTSTP  IN THE STEP RETURNED -- */
/*             SEE V(GRDFAC) ABOVE. */
/* V(PREDUC) (OUTPUT) FUNCTION REDUCTION PREDICTED FOR THE STEP RETURNED. */
/* V(RADIUS) (INPUT) THE TRUST REGION RADIUS.  D TIMES THE STEP RETURNED */
/*             HAS 2-NORM V(RADIUS) UNLESS V(STPPAR) = 0. */
/* V(STPPAR) (OUTPUT) CODE TELLING HOW STEP WAS COMPUTED... 0 MEANS A */
/*             FULL NEWTON STEP.  BETWEEN 0 AND 1 MEANS V(STPPAR) OF THE */
/*             WAY FROM THE NEWTON TO THE RELAXED NEWTON STEP.  BETWEEN */
/*             1 AND 2 MEANS A TRUE DOUBLE DOGLEG STEP, V(STPPAR) - 1 OF */
/*             THE WAY FROM THE RELAXED NEWTON TO THE CAUCHY STEP. */
/*             GREATER THAN 2 MEANS 1 / (V(STPPAR) - 1) TIMES THE CAUCHY */
/*             STEP. */

/* -------------------------------  NOTES  ------------------------------- */

/*  ***  ALGORITHM NOTES  *** */

/*        LET  G  AND  H  BE THE CURRENT GRADIENT AND HESSIAN APPROXIMA- */
/*     TION RESPECTIVELY AND LET D BE THE CURRENT SCALE VECTOR.  THIS */
/*     ROUTINE ASSUMES DIG = DIAG(D)**-2 * G  AND  NWTSTP = H**-1 * G. */
/*     THE STEP COMPUTED IS THE SAME ONE WOULD GET BY REPLACING G AND H */
/*     BY  DIAG(D)**-1 * G  AND  DIAG(D)**-1 * H * DIAG(D)**-1, */
/*     COMPUTING STEP, AND TRANSLATING STEP BACK TO THE ORIGINAL */
/*     VARIABLES, I.E., PREMULTIPLYING IT BY DIAG(D)**-1. */

/*  ***  REFERENCES  *** */

/* 1.  DENNIS, J.E., AND MEI, H.H.W. (1979), TWO NEW UNCONSTRAINED OPTI- */
/*             MIZATION ALGORITHMS WHICH USE FUNCTION AND GRADIENT */
/*             VALUES, J. OPTIM. THEORY APPLIC. 28, PP. 453-482. */
/* 2. POWELL, M.J.D. (1970), A HYBRID METHOD FOR NON-LINEAR EQUATIONS, */
/*             IN NUMERICAL METHODS FOR NON-LINEAR EQUATIONS, EDITED BY */
/*             P. RABINOWITZ, GORDON AND BREACH, LONDON. */

/*  ***  GENERAL  *** */

/*     CODED BY DAVID M. GAY. */
/*     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH SUPPORTED */
/*     BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS MCS-7600324 AND */
/*     MCS-7906671. */

/* ------------------------  EXTERNAL QUANTITIES  ------------------------ */

/*  ***  INTRINSIC FUNCTIONS  *** */
/* /+ */
/* / */
/* --------------------------  LOCAL VARIABLES  -------------------------- */


/*  ***  V SUBSCRIPTS  *** */


/*  ***  DATA INITIALIZATIONS  *** */

/* /6 */
/*     DATA HALF/0.5D+0/, ONE/1.D+0/, TWO/2.D+0/, ZERO/0.D+0/ */
/* /7 */
/* / */

/* /6 */
/*     DATA BIAS/43/, DGNORM/1/, DSTNRM/2/, DST0/3/, GRDFAC/45/, */
/*    1     GTHG/44/, GTSTEP/4/, NREDUC/6/, NWTFAC/46/, PREDUC/7/, */
/*    2     RADIUS/8/, STPPAR/5/ */
/* /7 */
/* / */

/* +++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++ */

    /* Parameter adjustments */
    --v;
    --step;
    --nwtstp;
    --dig;

    /* Function Body */
    nwtnrm = v[3];
    rlambd = 1.;
    if (nwtnrm > 0.) {
	rlambd = v[8] / nwtnrm;
    }
    gnorm = v[1];
    ghinvg = v[6] * 2.;
    v[45] = 0.;
    v[46] = 0.;
    if (rlambd < 1.) {
	goto L30;
    }

/*        ***  THE NEWTON STEP IS INSIDE THE TRUST REGION  *** */

    v[5] = 0.;
    v[2] = nwtnrm;
    v[4] = -ghinvg;
    v[7] = v[6];
    v[46] = -1.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L20: */
	step[i__] = -nwtstp[i__];
    }
    goto L999;

L30:
    v[2] = v[8];
/* Computing 2nd power */
    d__1 = gnorm / v[44];
    cfact = d__1 * d__1;
/*     ***  CAUCHY STEP = -CFACT * G. */
    cnorm = gnorm * cfact;
    relax = 1. - v[43] * (1. - gnorm * cnorm / ghinvg);
    if (rlambd < relax) {
	goto L50;
    }

/*        ***  STEP IS BETWEEN RELAXED NEWTON AND FULL NEWTON STEPS  *** */

    v[5] = 1. - (rlambd - relax) / (1. - relax);
    t = -rlambd;
    v[4] = t * ghinvg;
    v[7] = rlambd * (1. - rlambd * .5) * ghinvg;
    v[46] = t;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L40: */
	step[i__] = t * nwtstp[i__];
    }
    goto L999;

L50:
    if (cnorm < v[8]) {
	goto L70;
    }

/*        ***  THE CAUCHY STEP LIES OUTSIDE THE TRUST REGION -- */
/*        ***  STEP = SCALED CAUCHY STEP  *** */

    t = -v[8] / gnorm;
    v[45] = t;
    v[5] = cnorm / v[8] + 1.;
    v[4] = -v[8] * gnorm;
/* Computing 2nd power */
    d__1 = v[44] / gnorm;
    v[7] = v[8] * (gnorm - v[8] * .5 * (d__1 * d__1));
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L60: */
	step[i__] = t * dig[i__];
    }
    goto L999;

/*     ***  COMPUTE DOGLEG STEP BETWEEN CAUCHY AND RELAXED NEWTON  *** */
/*     ***  FEMUR = RELAXED NEWTON STEP MINUS CAUCHY STEP  *** */

L70:
    ctrnwt = cfact * relax * ghinvg / gnorm;
/*     *** CTRNWT = INNER PROD. OF CAUCHY AND RELAXED NEWTON STEPS, */
/*     *** SCALED BY GNORM**-1. */
/* Computing 2nd power */
    d__1 = cfact;
    t1 = ctrnwt - gnorm * (d__1 * d__1);
/*     ***  T1 = INNER PROD. OF FEMUR AND CAUCHY STEP, SCALED BY */
/*     ***  GNORM**-1. */
/* Computing 2nd power */
    d__1 = cfact;
    t2 = v[8] * (v[8] / gnorm) - gnorm * (d__1 * d__1);
    t = relax * nwtnrm;
    femnsq = t / gnorm * t - ctrnwt - t1;
/*     ***  FEMNSQ = SQUARE OF 2-NORM OF FEMUR, SCALED BY GNORM**-1. */
/* Computing 2nd power */
    d__1 = t1;
    t = t2 / (t1 + sqrt(d__1 * d__1 + femnsq * t2));
/*     ***  DOGLEG STEP  =  CAUCHY STEP  +  T * FEMUR. */
    t1 = (t - 1.) * cfact;
    v[45] = t1;
    t2 = -t * relax;
    v[46] = t2;
    v[5] = 2. - t;
/* Computing 2nd power */
    d__1 = gnorm;
    v[4] = t1 * (d__1 * d__1) + t2 * ghinvg;
/* Computing 2nd power */
    d__1 = v[44] * t1;
    v[7] = -t1 * gnorm * ((t2 + 1.) * gnorm) - t2 * (t2 * .5 + 1.) * ghinvg - 
	    d__1 * d__1 * .5;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L80: */
	step[i__] = t1 * dig[i__] + t2 * nwtstp[i__];
    }

L999:
    return 0;
/*  ***  LAST LINE OF DD7DOG FOLLOWS  *** */
} /* dd7dog_ */

