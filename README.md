SWAT
====

An updated SWAT 2012 revision 670 code

Objectives
----------

* Standard indentation and translation to Fortran 90 by using
[findent](https://sourceforge.net/projects/findent). See the
translate-fortran90.pl perl script file (:heavy_check_mark:)
* Exhaustive use of the "implicit none" directive to detect bad variable usage
(:heavy_check_mark:)
* Generate a GNU [make](http://www.gnu.org/software/make) makefile and compile
with GNU [gfortran](https://gcc.gnu.org/fortran). See the gernerate-makefile.pl
perl script file (:heavy_check_mark:)
* Remove non-used variables and format labels (:heavy_check_mark:)
* Detect and solve all uninitialized variables (:heavy_check_mark:, but some
solutions can be incorrect) 
* Remove unneeded variable initializations (:heavy_check_mark:) as:

`j=0 ! this line is not necessary`

`j=ihru`

* Remove redundant code (:heavy_check_mark:)
* Exhaustive use of the "parameter" directive on constants (:heavy_check_mark:)
* Remove obsolete commented code (:x:)
* Update variable descriptions in comments (:x:)
* Standardize comments by using Doxygen style on order to generate documentation
(:x:)

Required tools
--------------

* [gfortran](https://gcc.gnu.org/fortran) (to compile the source code)
* [make](https://www.gnu.org/software/make) (to build the executable file)
* [perl](https://www.perl.org) (optional: to execute the perl scripts to 
update the makefile or to translate original files to Fortran 90)
* [findent](https://sourceforge.net/projects/findent) (optional: to translate
original files to Fortran 90 with a standard indentation)
* On Microsoft Windows systems you have to install
[MSYS2](http://sourceforge.net/projects/msys2) and the required utilities
(gfortran and make). You can follow detailed instructions in
[install-unix](https://github.com/jburguete/install-unix/blob/master/tutorial.pdf)

Instructions to generate an executable to test
----------------------------------------------

* In UNIX type operative systems:
> $ make

* Cross-compiling a 32 bits Microsoft Window executable in a UNIX type operative
system:
> $ prefix="i686-w64-mingw32-" EXE=".exe" LDFLAGS="-static" make

* Cross-compiling a 64 bits Microsoft Window executable in a UNIX type operative
system:
> $ prefix="x86\_64-w64-mingw32-" EXE=".exe" LDFLAGS="-static" make


Instructions to generate of an optimized executable file
--------------------------------------------------------

* In UNIX type operative systems:
> $ CFLAGS="-march=native -flto" LDFLAGS="-flto" make strip

* Cross-compiling a 32 bits Microsoft Window executable in a UNIX type operative
system:
> $ prefix="i686-w64-mingw32-" EXE=".exe" CFLAGS="-flto" LDFLAGS="-flto -static"
> make strip

* Cross-compiling a 64 bits Microsoft Window executable in a UNIX type operative
system:
> $ prefix="x86\_64-w64-mingw32-" EXE=".exe" CFLAGS="-flto"
> LDFLAGS="-flto -static" make strip

