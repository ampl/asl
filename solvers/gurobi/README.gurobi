The solver "gurobi" uses Gurobi (a trademark of Gurobi Optimization,
Inc.; see http://www.gurobi.com/) to solve integer, mixed-integer, and
linear programming problems.  Normally gurobi is invoked by AMPL's
solve command, which gives the invocation

     gurobi stub -AMPL

in which stub.nl is an AMPL generic output file (possibly written
by "ampl -obstub" or "ampl -ogstub").  After solving the problem,
gurobi writes a stub.sol file for use by ampl's solve and solution
commands.  When you run ampl, this all happens automatically if you
give the AMPL commands

     option solver gurobi;
     solve;

You can control gurobi by setting the environment variable gurobi_options
appropriately (either by using ampl's option command, or by using the
shell's set and export commands before you invoke ampl).  You can put
one or more (white-space separated) phrases in $gurobi_options.  To see
the possibilities, invoke

        gurobi '-='

----------
INSTALLING
==========

On Linux systems, libgurobi*.so (where the value of "*" depends
on the current version of Gurobi) and the libgurobi.so.* to which
it points need to appear in the current directory when gurobi
itself appears there, or in one of the standard places (specified by
/etc/ld.so.conf on some systems), or in a directory named in
$LD_LIBRARY_PATH.  An alternative is to add a short shell script,
such as

        #!/bin/sh
        LD_LIBRARY_PATH=/usr/local/lib
        export LD_LIBRARY_PATH
        exec /usr/local/bin/gurobix "$@"

to a directory in your usual $PATH (and mark the script executable
with, e.g., "chmod +x gurobi").  The above script assumes that the
true "gurobi" binary has been moved to /usr/local/bin/gurobix and that
the libgurobi* files have been moved to /usr/local/lib.

MacOSX systems are similar to Linux systems, but with DYLD_LIBRARY_PATH
in place of LD_LIBRARY_PATH.  Starting 20150225, gurobi binaries
for MacOSX should find the appropriate libgurobi.so.* if it appears
in the same directory as "gurobi".

On MS Windows systems, gurobi.exe and the relevant gurobi*.dll must
appear somewhere in your usual search $PATH (or in the current
directory).

AIX systems are similar to Linux systems, but with LIBPATH in place of
LD_LIBRARY_PATH and libgurobbi*.a in place of libgurobi*.so*.


-----------------------
solve_result_num values
=======================

Here is a table of solve_result_num values that "gurobi" can return
to an AMPL session, along with an indication of the text that appears
in the associated solve_message.

        Value   Message

          0     optimal solution
          1     optimal solution with integer variables rounded to integers
          2     optimal solution with nonintegral "integer" variables
        100     suboptimal: could not satisfy optimaliter tolerances
        101     bestbndstop reached
        102     bestobjstop reached
        103     bestobjstop or bestbndstop reached
        104     bestobjstop or bestbndstop reached with no solution available
        200     infeasible [IIS computation not attempted]
        201     infeasible [IIS returned]
        202     infeasible [IIS finder failed]
        203     infeasible; .dunbdd returned [IIS computation not attempted]
        204     infeasible; .dunbdd returned [IIS also returned]
        205     infeasible; .dunbdd returned [IIS finder failed]
        300     unbounded
        301     unbounded [unbounded or infeasible; IIS finder failed]
        302     unbounded; .unbdd returned
        303     unbounded; .unbdd returned [IIS finder failed]
        400     objective cutoff
        401     iteration limit with a feasible solution
        402     node limit with a feasible solution
        403     time limit with a feasible solution
        404     solution limit
        405     interrupted with a feasible solution
        411     iteration limit without a feasible solution
        412     node limit without a feasible solution
        413     time limit without a feasible solution
        415     interrupted without a feasible solution
        500     Could not create the gurobi environment
        501     Gurobi call failed [message gives routine name]
        502     misc. failure [message gives details]
        503     Bad $gurobi_options
        504     Surprise VBasis[...] = ...
        505     Surprise CBasis[...] = ...
        506     Gurobi set/get parameter failed [message gives more details]
        510     cannot open logfile (specified in $gurobi_options or command line)
        511     cannot open paramfile (specified in $gurobi_options or command line)
        512     missing value in paramfile
        513     extra text in paramfile
        514     invalid parameter name in paramfile
        520     numeric error
        521     nonlinear objective
        522     nonlinear constraint
        523     quadratic objective involving division by 0
        524     indefinite quadratic objective or constraint
        525     quadratic constraint involving division by 0
        530     could not open serverlic file
        531     error in serverlic file
        532     error while tuning
        540     Gurobi Compute Server not reached or bad pool_... settings
        541     Rejected by Gurobi Compute Server
        542     Feature not supported by Gurobi Compute Server
        543     Feature not supported
        563     logical constraint is not an indicator constraint
        564     bad suffixes for multiple objectives
        565     bug? Error return from named routine
        567     complementarity constraint
        570     solution found but not available (Gurobi bug?)
        571     expected just one solution when problem has no integer variables
        601     could not talk to Gurobi compute server
        602     job rejected by Gurobi compute server
        603     no license for specified gurobi compute server
        604     surprise return while trying to use Gurobi compute server
        605     bad value for cloudid or cloudkey, or Gurobi Cloud out of reach

Values 521-524 only arise in Gurobi versions >= 4.0.
Values 203-205 and 302-303 only arise in Gurobi versions >= 4.5.
Values 530-531 and 601-603 only arise in Gurobi versions >= 5.5.
Values 540-542 only arise in Gurobi versions >= 6.5.
Value  605 only arises in Gurobi versions >= 7.0.

*************************

If you have questions about or find bugs with this stuff,
please contact:

     David M. Gay
     dmg@ampl.com
