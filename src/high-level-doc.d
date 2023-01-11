
/*!

\mainpage High-Level LARC Documentation

\section intro_sec Introduction

LARC is a software package developed at IDA/CCS that
stores matrices in a compressed format and performs operations on
matrices and vectors while they are compressed.  LARC assigns each
scalar, vector, and matrix a unique \a MatrixID.  These MatrixIDs
are used with two hash tables, the \a MatrixStore and the
\a OperationStore, to ensure that we never store the same matrix
twice or carry out the same operation twice.
We have developed several new techniques using a locality-sensitive
scalar hash which allow LARC to reduce numerical precision issues and
to preserve selected math identities.
LARC is optimized for
use on matrices with power-of-2 dimensions and is most effective when
there is submatrix reuse, which is common in recursive matrices,
self-similar matrices, matrices with tensor structure, and matrices
with a limited number of distinct scalars (e.g., sparse matrices).


\section data_sec Data Structure

In LARC, each scalar, vector, and matrix is assigned a unique MatrixID
and given a \a MatrixRecord that is stored in the MatrixStore hash table.
Inside the MatrixRecord, a matrix is represented recursively 
by a list of the MatrixIDs from its four quadrant submatrices;
we call this list its \a MatrixValue.  If the matrix is 1 by 1,
then its MatrixValue is simply the scalar itself.

For non-scalar matrices, each MatrixRecord in the
MatrixStore is indexed by a hash of its
four quadrant submatrix MatrixIDs from its recursive definition.
In order to not store the same matrix twice, we use the MatrixStore
hash table to confirm whether a matrix is new or has already been stored.


\section scalar_type Scalar Types

At compile time, LARC users select the underlying data type for scalars
inside of LARC matrices by compiling with a command such as:

    make TYPE=INTEGER

The scalarType also determines the appropriate scalar arithmetic. The
scalarTypes currently supported by LARC are:

scalarType | Implementation
---------- | --------------
Real (default) | C long double
Integer | C int64 t
Complex | C long double complex
MPInteger | GMP multiprecision mpz_t
MPRational | GMP multiprecision mpq_t
MPReal | GMP multiprecision mpfr_t with 256 bits of precision
MPComplex | GMP multiprecision mpc_t with 256 bits of precision
MPRatComplex | a specialized LARC structure with two fields, real and imag, each of type mpq_t


\section hash_info Locality-Sensitive Hashing

We use a locality-sensitive hash to handle finite precision issues
and to mimic symbolic/exact computation.

Finite precision can lead to math identities going wrong, such as
1−(1−a) &ne; a.
Without some way to compensate for this in LARC, this would lead to creation
of two or more scalars that differ only by some small amount due to numerical
precision issues. This could become a large problem: not only would we need to
store multiple copies of what should be the same scalar, we would also need to
store many matrices which would be equal but for these small errors.

We have developed a technique involving locality-sensitive hashing.
Locality-sensitive hashing is an algorithmic technique that hashes similar input items
into the same hash buckets with high probability. In our case ’similar’ means
that two scalars are close to equal. The locality-sensitive hash ensures that
scalars that would be in the same region hash to the same value, making it
fast to find a previously stored representative of that region or to determine
that no such representative exists.

When attempting to store a scalar value in LARC’s MatrixStore, we first
want to determine whether a sufficiently nearby scalar has already been
stored, that we should use instead. We detect this with high probability by
dividing the space of scalars into small regions, defined by a user-controlled
parameter, and hashing all scalars that would be in the same region to the
same hash bucket. The first scalar found in a particular region is stored in
the MatrixStore and we call this scalar the representative of that region.
Attempting to store any new value that would lie in an already occupied region
will result in LARC returning the MatrixID of the original representative
and not storing the new value.

This technique could clearly cause additional issues with failing identities
when substitute scalars are used in mathematical operations. For example, if
we have already stored the number 22 in LARC, and we try to store the scalar
s = 7π ≈ 21.9911485751, and we have set the size of regions to be large e.g.
1/32, then LARC will return the MatrixID for 22 instead of storing s. This
can lead to further errors (see examples of how choices of initial parameters
can affect this in MyPyLARC/Tutorial/preloading mults pi.py).
Two parameters, `regionbitparam` and `zeroregionbitparam`, are passed to LARC
during the initialization routine, and determine the size of the regions
used in the hash function.

We have an additional technique called preloading which we use to guarantee
that preselected identities are satisfied. For example, since zero is
particularly important to mathematical identities,
it is the first thing that we store
when initializing LARC. This ensures that zero is the representative of its
region, and any scalar within the region containing zero will be treated as
zero. Depending on the application, we may also choose to preload a
select set of scalars such as the n-th roots of unity, to ensure they become the
representatives of their regions and that mathematical identities with these
scalars hold (see examples of this preloading to ensure identities in
MyPyLARC/Tutorial/examples roots unity.py and fft paramFile init.py).

Because of the importance of identities involving zero, LARC gives the user
the choice to make the region around zero larger than the other regions used
by the locality sensitive hash.


\section matrixop_sec Matrix Operations

LARC can carry out linear algebra operations on compressed matrices
(leaving them in the compressed format), as long as the algebraic
operation can be described recursively in terms of quadrant
submatrices.
LARC currently contains several basic unitary and binary operations
including: matrix addition, matrix-matrix multiplication,
scalar-matrix multiplication, 
Kronecker product, adjoint (complex conjugate
transpose), and some matrix norms.

A list of common matrix operations can be found
<A HREF="matmath_8h.html">here</A>.

Operations are memoized and saved as an \a OperationRecord
in the OperationStore.
Just as with the MatrixStore, LARC implements this as a hash table.
The OperationStore hash function uses the type of operation and the
MatrixIDs of the input matrices.

LARC takes advantage of mathematical identities in several ways.
MatrixRecords of identity matrices and all zero matrices contain
flags marking them as such.  Then when a matrix \a M is added to a
zero matrix, or multiplied by an identity matrix, the 
MatrixID of \a M is returned without any computational work.
These zero and identity flags are used for initial checks in
many LARC operations, e.g. adjoint of these matrices is trivial.
Also, since matrix addition is commutative, we use the trick
of sorting the MatrixIDs of the inputs before including them
in the hash computation for the OperationStore.  Thus having
memoized the operation \a A + \a B serves to memoize \a B + \a A
as well.

When LARC is asked to perform an operation, it uses any identity
short cuts, then checks 
the
OperationStore to see if the operation has already been performed, and
if so returns the MatrixID of the result found in the OperationRecord.
If the operation has not been memoized, then LARC carries out the
operation.  Addition, multiplication and all other matrix operations
in LARC are carried out recursively in a way that produces the four
MatrixIDs of the quadrant submatrices of the output matrix.  As with
any matrix LARC attempts to store, it hashes the list of quadrant
MatrixIDs and then looks in that hash chain of the MatrixStore
to see if a matrix with
those quadrant submatrices has already been stored. If the matrix
is found its MatrixID is returned, otherwise a new MatrixRecord
is created and its MatrixID is returned.
Finally, the operation is now memoized in the OperationStore.

Even if the top level of an operation
has not been memoized, LARC may have previously memoized
operations further down the quadtree recursion, allowing LARC
to skip the work for all branches below the point of that memoization.
This cut down of work is similar to the cut down of memory
that occurred in the MatrixStore when a submatrix of a larger
matrix was previously stored.
The combination of compressed storage, carrying out operations in
compressed format, and memoizing operations allows LARC, and
application packages which use LARC,
to carry out computations that would otherwise seem to be intractable.


\section io_sec Input and Output

LARC has its own compressed format to describe a single matrix; this
can be used for input from a file to LARC and output from LARC to a
file. These files use a JSON container and have a separate line for
each unique MatrixRecord containing the MatrixID, the level, and
the MatrixValue (which is either the list of MatrixIDs of the four
quadrants or for scalar matrices, the value of the scalar).
These LARC compressed matrix files can also contain metadata
which LARC allows the user to record in an InfoStore indexed by the
MatrixID. The metadata in the InfoStore allows a user to track
parameters of interest such as scalar data type used, locality-sensitive
hash parameter used, and comments.  We designed the InfoStore to
make logging of experimental parameters easy for input/output. Since
there is no compression in the InfoStore, it is not intended for
tracking information on very large sets of matrices.

The ability
to store compressed matrices in files allows us to checkpoint
important intermediate results in long runs. We can also write staged
algorithms that output results and then start a new run loading the
results from the previous program. This has the advantage of removing
unneeded records and reducing memory requirements.

In addition to the LARC compressed matrix format, 
LARC has the capability to read files containing matrices expressed in
row-major format (header with dimensions, followed by items listed
one at a time reading along rows).
Using row-major format only makes sense for small matrices.
LARC can also read sparse matrices in the 
Matrix Market Exchange Format (header with various info including
number of nonzero items, followed by a line for each nonzero item
with its coordinates and content).
When LARC reads in a matrix, as always, it does not store copies of
any matrix or submatrix that is already in the MatrixStore.

LARC can output small matrices in row-major format and larger matrices
in LARC compressed format.  LARC currently does not output matrices
in Matrix Market Exchange Format.


\section mem_challenge Addressing Memory Challenges

The LARC package speeds computations by storing matrices and memoizing
operations when the matrices involved have repeated submatrices. However,
the memory requirements for some computations can become very large, so
it is sometimes necessary to give up some previously stored matrices and
memoized operations to permit continued computation.

If a computation can be divided into several programs which each culminate
in one or more matrices containing the results from that stage, then LARC
can take advantage of this staged computation. Each stage starts by
initializing LARC, and reading one or more files with the results from the last
stage; each stage ends by writing its results into one or more files. In this
way, unneeded MatrixRecords and OperationRecords from previous stages
do not occupy memory.

Another way to reduce memory requirements is to use the routines that
LARC has available for cleaning. By cleaning, we mean the removal of certain
matrices from the MatrixStore and removal of those OperationRecords from
the OperationStore that contain removed matrices. The cleaning routines
track recursive dependencies in the MatrixStore so that no MatrixRecord is
removed if another MatrixRecord refers to it. This means that a request to
remove a matrix will only be carried out if it does not violate this
recursive closure condition. When a matrix is removed, LARC also removes any
submatrices that are only used by this matrix.

When a matrix has been removed from the MatrixStore, it may still be
referred to in an OperationStore record. These obsolete OperationRecords are
removed whenever they are touched during the traversal of a OperationStore
hash chain while searching for some memoized operation. There is also a
cleaning option to empty the entire OperationStore if desired; this is most
useful if the computation is entering a new stage in which few previously
memoized operations are likely to be repeated.


\section clean_store Cleaning the Stores

The LARC package speeds computations by storing matrices and memoizing
operations when the matrices involved have repeated submatrices. However,
the memory requirements for some computations can become very large, so
it is sometimes necessary to free the memory of previously stored matrices
and memoized operations to permit continued computation.

One way to do this is to remove matrices from the MatrixStore hash table, as
well as removing those operations from the OperationStore hash table which
recorded operations on those matrices which have been removed. LARC
provides various routines to clean these stores. (In general, the InfoStore does
not take up a significant amount of memory, but users are provided routines
which allow them to clean it in the same way that the OperationStore is
cleaned.)

Removal of a matrix from the MatrixStore entails: freeing the memory
allocated for its MatrixRecord; removing the hash node in the MatrixStore hash
chain which referred to that record; updating relevant statistics that LARC
keeps about the contents of the MatrixStore; and indicating in a table
relating MatrixIDs to matrix pointers that the pointer for that MatrixID is
now invalid. (MatrixIDs are guaranteed to be unique. If a matrix is removed
from the MatrixStore and later recreated, it will not receive the old MatrixID
value but will instead be given a new one, since we deleted the information
that would confirm its pre-existence.)

We must be careful not to remove any matrix which is a submatrix of some
larger matrix still in the MatrixStore, as such loss of information would be
catastrophic to any operation performed on the larger matrix. To ensure this
does not happen, each MatrixRecord contains a counter which indicates the
number of matrices for which it is a quadrant submatrix. If this value is
greater than zero, the cleaning routine will not remove the matrix from the
MatrixStore. However, if it is zero (and other factors discussed in the next
paragraph do not apply), the matrix will be removed, and all of its quadrant
submatrices will have their own submatrix counters decremented; each of
these is then recursively checked to see if they in turn may be removed from
the MatrixStore.

To complement this submatrix counting mechanism, we also provide lock and
hold functions. Matrices that are either locked or held are exempted from
cleaning, even if their submatrix count is zero. (All submatrices of locked
and held matrices are guaranteed to have submatrix count of at least 1, and
are therefore also exempt from cleaning.) The lock flag for a matrix is set
to 1 when that matrix should never be cleaned, and once locked a matrix
may not be unlocked. LARC’s initialization locks certain important scalars
such as zero and one, as well as all-zero and identity matrices of various
sizes. Users may lock scalars or matrices that are particularly important to
their applications. In contrast, the hold field in a MatrixRecord is a counter
which is incremented and decremented as holds are placed and released; this
mechanism is used to guarantee a matrix is not cleaned until the user knows
it is no longer needed.

Cleaning of the MatrixStore is integrated into LARC’s matrix multiply
routine. While the standard matrix mult() routine does no cleaning, it is built
on top of a routine called matrix mult clean(). This latter routine allows
the user to clean matrices which are not part of the final result once they
have been used. The amount of memory saved, and the consequent increase
in time as operations are repeated, depends on the application and on how
aggressively cleaning is pursued.

The OperationStore is cleaned with a simpler algorithm, which only needs to
free the memory allocated for the OperationRecord and remove the node in
the OperationStore hash chain which referred to that record. As the absence
of a stored operation merely means that the operation must be performed,
removing an operation cannot cause an error in any calculation. Consequently,
much of the machinery devoted to preventing the incorrect deletion of a
matrix is unnecessary when cleaning operations.

LARC contains routines to clean a specified OperationStore hash chain, to
clean all hash chains in the OperationStore, and even to completely remove all
nodes from the entire OperationStore hash table. (While emptying the entire
OperationStore is safe, it may significantly increase the time needed for an
application, as any operations that had been in the store before the cleaning
are no longer memoized until they are recalculated.) When the user specifies
a hash chain to be cleaned, the hash chain is traversed; any node containing
an operation which references a matrix with invalid ID is removed, and its
OperationRecord is freed. Some cleaning is also performed automatically
during LARC operations when checking to see if a particular operation has
already been stored in an OperationRecord. In the process of trying to find
such an OperationRecord, we must look at each node of a particular hash
chain (corresponding to the hashed value of the operation). This traversal
ends by either finding the desired OperationRecord or by reaching the last
node in the chain without finding it. As we look at each OperationRecord in
turn, LARC checks the input matrices and the output matrix to see whether
they are still in the MatrixStore. If any of their matrix pointers are now equal
to NULL, then this OperationRecord is no longer valid and is removed from
the hash chain (as is typical, the nodes on either side of the removed node are
linked). Only that part of the OperationsStore hash chain that is traversed
during the search for an OperationRecord is cleaned.


\section python_sec Python Interface

In addition to a base in C code, there is a Python interface.
This interface is automatically generated by SWIG.
There are additional routines in pylarc.py that facilitate using LARC
from Python code.

There is a tutorial called MyPyLARC that is available from github.com/LARCmath/MyPyLARC.git.
It has sample applications that can be used as a starting point for your applications.

Additionally, if you are making any changes to your copy of LARC, you should be aware that
there are unit tests (not fully complete) you can run by issuing the command:

    make unittests

This command runs unit tests written in python from the LARC:tests/python directory.


*/


