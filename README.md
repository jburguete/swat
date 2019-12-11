SWAT
====

An updated SWAT 2012 revision 670 code

Objectives
----------

* Standard indentation and translation to Fortran 90 by using
[findent](https://sourceforge.net/projects/findent). See the
translate-fortran90.pl perl script file
* Exhaustive use of the "implicit none" directive to detect bad variable usage
* Generate a GNU [make](http://www.gnu.org/software/make) makefile and compile
with GNU [gfortran](https://gcc.gnu.org/fortran). See the gernerate-makefile.pl perl script file
* Remove non-used variables and format labels
* Detect and solve all uninitialized variables
* Remove unneeded variable initializations as:
> j=0 ! this line is not necessary
>
> j=ihru

Required tools
--------------

* [gfortran](https://gcc.gnu.org/fortran) (to compile the source code)
* [make](https://www.gnu.org/software/make) (to build the executable file)
* [perl](https://www.perl.org) (optional: to execute the perl scripts to 
update the makefile or to translate original files to Fortran 90)
* [findent](https://sourceforge.net/projects/findent) (optional: to translate
original files to Fortran 90 with a standard indentation)
* On Microsoft Windows systems you have to install
[MSYS2](http://sourceforge.net/projects/msys2) and the required
libraries and utilities. You can follow detailed instructions in
[install-unix](https://github.com/jburguete/install-unix/blob/master/tutorial.pdf)

Generation of an executable to test
-----------------------------------

* In UNIX type operative systems:
> $ make

* Cross-compiling a 32 bits Microsoft Window executable in a UNIX type operative
system:
> $ make clean
>
> $ prefix="i686-w64-mingw32-" EXE=".exe" LDFLAGS="-static" make

* Cross-compiling a 64 bits Microsoft Window executable in a UNIX type operative
system:
> $ make clean
>
> $ prefix="x86\_64-w64-mingw32-" EXE=".exe" LDFLAGS="-static" make


Generation of an optimized executable file
------------------------------------------

* In UNIX type operative systems:
> $ make clean
>
> $ CFLAGS="-march=native" make strip

* Cross-compiling a 32 bits Microsoft Window executable in a UNIX type operative
system:
> $ make clean
>
> $ prefix="i686-w64-mingw32-" EXE=".exe" LDFLAGS="-static" make strip

* Cross-compiling a 64 bits Microsoft Window executable in a UNIX type operative
system:
> $ make clean
>
> $ prefix="x86\_64-w64-mingw32-" EXE=".exe" LDFLAGS="-static" make strip

