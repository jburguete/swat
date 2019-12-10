SWAT
====

An updated SWAT 2012 revision 670 code

Objectives
----------

* Standard indentation and translation to Fortran 90 by using
[findent](https://sourceforge.net/projects/findent). See the
translate-fortran90.pl perl script file
* Exhaustive use of the "implicit none" directive to detect bad variable usage
* Generation of a GNU [make](http://www.gnu.org/software/make) makefile and
compilation with GNU [gfortran](https://gcc.gnu.org/fortran)
* Remove non-used variables

Required tools
--------------

* [gfortran](https://gcc.gnu.org/fortran) (to compile the source code)
* [make](http://www.gnu.org/software/make) (to build the executable file)
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
> $ prefix="i686-w64-mingw32-" EXE=".exe" LDFLAGS="-static" make

* Cross-compiling a 64 bits Microsoft Window executable in a UNIX type operative
system:
> $ make clean
> $ prefix="x86\_64-w64-mingw32-" EXE=".exe" LDFLAGS="-static" make


Generation of an optimized executable file
------------------------------------------

* In UNIX type operative systems:
> $ make clean
> $ CFLAGS="-march=native" make strip

* Cross-compiling a 32 bits Microsoft Window executable in a UNIX type operative
system:
> $ make clean
> $ prefix="i686-w64-mingw32-" EXE=".exe" LDFLAGS="-static" make strip

* Cross-compiling a 64 bits Microsoft Window executable in a UNIX type operative
system:
> $ make clean
> $ prefix="x86\_64-w64-mingw32-" EXE=".exe" LDFLAGS="-static" make strip

