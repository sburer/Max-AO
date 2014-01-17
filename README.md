Max-AO
======

A heuristic for finding large cliques and stable sets

 ==============================================
                 a heuristic for finding large
                   cliques and stable sets
  ** Max-AO **   
                  Sam Burer, Renato Monteiro
                   Georgia Tech, v0.120301
 ==============================================

===> Web Page and Contacts <===

http://dollar.biz.uiowa.edu/~sburer

Samuel Burer (samuel-burer@uiowa.edu)
Renato D.C. Monteiro (monteiro@isye.gatech.edu)


===> Compiling <===

cc -O2 -o max-ao Max-AO_v0.120301.c -lm     (or something similar)


===> Running <===

max-ao <graph file> <parameter file>        (try "max-ao brock200_1 p.c11")


===> Graph File <===

Max-AO accepts graphs encoded in a standard text file. The encoding
is as follows:

  * The first line of the text file contains "n m" where n is the
    number of vertices and m is the number of edges.

  * The next m lines contain "i j w" where (i,j) represents an edge
    and w is a weight associated with the edge. (The weight is
    ignored by Max-AO.)

For example, the cycle on five vertices could be encoded as

  5 5
  1 2 1
  2 3 1
  3 4 1
  4 5 1
  5 1 1

Note that Max-AO does not worry whether i > j or j > i but does
assume that the vertices are labelled 1 through n (as opposed to 
0 through n-1, for example).


===> Parameter File <===

The parameter file consists of 7 lines of text with an integer on
each line. The parameters controlled by each line are as follows:

 1. Whether to find stable sets or cliques; 0 for stable set and
    1 for clique.

 2. Which rank to use; either 1 or 2.

 3. How many escapes to try; any positive integer. It is suggested
    that 1 is used if rank is 1.

 4. The maximum number of nonlinear iterations allowed; any
    positive integer. It is suggested to choose a large number,
    say 100000.

 5. The maximum time (in seconds) to allow; any positive integer.

 6. The number of limited-memory BFGS updates to store; any positive
    integer. It is suggested to choose 2 or 3.

 7. Whether or not to save the best stable set or clique found;
    0 for no and 1 for yes.

Six example parameter files have been provided. They differ only in
the first three parameters and are categorized as follows:

  p.<mode><rank><escapes>

     mode = c or s (clique or stable set)
     rank = 1 or 2
  escapes = 1 or 5


===> Output <===

Max-AO will first print some information to the screen regarding
the input graph as well as the options that have been selected.
During the running of the heuristic, it will print out a
line of information each time a new, better stable set or clique
has been found; the line is in the form "<size> <time found>".

If the parameter for saving the best stable set or clique to a
file has been set to 1, then Max-AO will create a file with an
appropriate suffix in which the vertices of the best stable set
or clique are listed on consecutive lines.


===> Miscellaneous <===

Max-AO currently does very little error checking, and so it is
important that the user verify his or her input.

The code contains very few comments but is not particularly
complicated or long. Thus, it should be possible for a user
to alter the code for their own purposes.

