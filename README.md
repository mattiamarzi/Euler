The code is written in Fortran and can be launched using the Makefile in the following way:

make

To solve the differential equations using Lax-Wendroff scheme:
./Euler L

To solve the differential equations using Roe scheme:
./Euler R

To solve the differential equations using MUSCL scheme with MINMOD flux-limiter:
./Euler M MIN

To solve the differential equations using MUSCL scheme with SUPERBEE flux-limiter:
./Euler M SUP

The plot file can be used to produce .gif representations of the solutions obtained at different times using any of the above schemes.
