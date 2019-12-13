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
* Generate a GNU [Make](http://www.gnu.org/software/make) makefile and compile
with GNU [GFortran](https://gcc.gnu.org/fortran). See the gernerate-makefile.pl
perl script file (:heavy_check_mark:)
* Remove non-used variables and format labels (:heavy_check_mark:)
* Detect and solve all uninitialized variables (:heavy_check_mark:
:construction:, some proposed solutions could be incorrect) 
* Remove unneeded variable initializations (:heavy_check_mark:) as:

  `j=0 ! this line is not necessary`

  `j=ihru`

* Remove redundant code (:heavy_check_mark:)
* Exhaustive use of the "parameter" directive on constants (:heavy_check_mark:)
* Generate a detailed list of issues detected in the original code
(:heavy_check_mark:, see at the end of this README)
* Remove obsolete commented code (:x:)
* Update variable descriptions in comments (:construction:, a lot of work)
* Standardize comments by using Doxygen style in order to generate
documentation. See at latex/refman.pdf (:construction:, a lot of work)

Required tools
--------------

* [GFortran](https://gcc.gnu.org/fortran) (to compile the source code)
* [Make](https://www.gnu.org/software/make) (to build the executable file)
* [Perl](https://www.perl.org) (optional: to execute the perl scripts to 
update the makefile or to translate original files to Fortran 90)
* [Findent](https://sourceforge.net/projects/findent) (optional: to translate
original files to Fortran 90 with a standard indentation)
* [Doxygen](http://www.doxygen.nl) (optional: to generate a reference
programming manual from source code)
* [TeX Live](https://www.tug.org/texlive) or [MiKTeX](https://miktex.org/)
(optional: to generate a reference programming manual from source code)
* On Microsoft Windows systems you have to install
[MSYS2](http://sourceforge.net/projects/msys2) and the required utilities
([GFortran](https://gcc.gnu.org/fortran) and
[Make](https://www.gnu.org/software/make)). You can follow detailed instructions
in
[install-unix](https://github.com/jburguete/install-unix/blob/master/tutorial.pdf)

Instructions to generate Fortran 90 style code from original code
-----------------------------------------------------------------
In order to generate Fortran 90 style code with standard indentation from
original code you have to type on a UNIX type terminal (you need
[Perl](https://www.perl.org) and
[Findent](https://sourceforge.net/projects/findent)):
> $ perl translate-fortran90.pl

Instructions to generate an initial GNU make Makefile
-----------------------------------------------------

Type on the UNIX type terminal, when translated the original code to Fortran 90
style (you need [Perl](https://www.perl.org)):
> $ perl generate-makefile.pl

Instructions to generate an executable to test
----------------------------------------------

Type on the UNIX type terminal (you need [GFortran](https://gcc.gnu.org/fortran)
and [Make](https://www.gnu.org/software/make))

* In UNIX type operative systems:
> $ make

* In a [MSYS2](http://sourceforge.net/projects/msys2) terminal in Microsoft
Windows:
> $ EXE=".exe" LDFLAGS="-static" make

* Cross-compiling a 32 bits Microsoft Windows executable in a UNIX type
operative system:
> $ prefix="i686-w64-mingw32-" EXE=".exe" LDFLAGS="-static" make

* Cross-compiling a 64 bits Microsoft Windows executable in a UNIX type
operative system:
> $ prefix="x86\_64-w64-mingw32-" EXE=".exe" LDFLAGS="-static" make


Instructions to generate an optimized executable file
-----------------------------------------------------

Type on the UNIX type terminal (you need [GFortran](https://gcc.gnu.org/fortran)
and [Make](https://www.gnu.org/software/make))

* In UNIX type operative systems:
> $ CFLAGS="-march=native -flto" LDFLAGS="-flto" make strip

* In a [MSYS2](http://sourceforge.net/projects/msys2) terminal in Microsoft
Windows:
> $ EXE=".exe" CFLAGS="-flto" LDFLAGS="-flto -static" make strip

* Cross-compiling a 32 bits Microsoft Windows executable in a UNIX type
operative system:
> $ prefix="i686-w64-mingw32-" EXE=".exe" CFLAGS="-flto" LDFLAGS="-flto -static"
> make strip

* Cross-compiling a 64 bits Microsoft Windows executable in a UNIX type
operative system:
> $ prefix="x86\_64-w64-mingw32-" EXE=".exe" CFLAGS="-flto"
> LDFLAGS="-flto -static" make strip

Instructions to generate a reference programming manual from source code
------------------------------------------------------------------------

Type on the UNIX type terminal (you need [Doxygen](http://www.doxygen.nl) and
[TeX Live](https://www.tug.org/texlive) or [MiKTeX](https://miktex.org/)):
> $ doxygen
>
> $ cd latex
>
> $ make

The reference programming manual file latex/refman.pdf is generated from source
code in PDF format

Issues in the original source code
----------------------------------

This is a list of possible issues detected in the original source code. These
issues have been mostly detected by the [GFortran](https://gcc.gnu.org/fortran)
(to compile the source code) compiler warnings. Some of them could not arise
because the logic of the variables is not possible.

* In biofilm.f:
  - "dcoef" is not defined.
    `dcoef=3` as in watqual.f?
    Then, I propose at beginning:
    `real*8, parameter :: dcoef = 3.`

* In bmp\_ri\_pond.f:
  - "qseep" and "qet" could be used not initialized at lines 133 and 134. 
    However the problem only arises for `nstep<1`

* In bmp\_sand\_filter.f:
  - "sed\_removed" at line 342 could be used not initialized if `sfsedstdev<=0`

* In bpm\_sed\_pond.f:
  - `bmp_sed _pond` seems to be `bmp_sed_pond` at line 186

* In bmp\_wet\_pond.f:
  - "hvol" could be used not initialized in "ext\_dpth" subroutine at line 267
    in first bucle iteration

* In clicon.f:
  - "tmxbsb", "tmnbsb", "rbsb", "rstpbsb", "rhdbsb", "rabsb", "rmxbsb",
    "daylbsb", "fradbsb" and "u10bsb" could be used not initialized at 186-207
    lines

* In conapply.f:
  - "k" and "kk" could be used not initialized at 121-122 lines if
    `iday_pest(j)/=ipst_freq(j)` and `curyr>nyskip`

* In confert.f:
  - "ifrt" seems to be "it" at line 214

* In curno.f:
  - "smxold" could be used not initialized if `cn1(h)<=1.e-6` and `curyr/=0` at
    line 96

* In drains.f:
  - "nlayer" could be used not initialized at line 23. However, the problem only
    arises if it is not set in the previous bucle (`mlyr<=1` or
    `sol_z(j1,j)<=0`)

* In etact.f:
  - "sev" could be used not initialized at line 286 if `dep>=esd` and `ly==2`

* In filter.f:
  - "remove21" seems to be "remove2" at line 316

* In grass\_wway.f:
  - "sf\_depth" and "sf\_sed" could be used not initialized at lines 133 and 137
    if `sf_area>0` and `sf_area<=1.e-6`

* In hhnoqual.f:
  - "algon" seems to be "algcon" at line 190

* In hhwatqual.f
  - "orgnpin" seems to be "orgpin" at line 278
  - `thour=1.0` at line 377 overwrites previous "thour" calculation. It is wrong

* In hmeas.f:
  - "rhdbsb" could be used not initialized at line 84

* In killop.f:
  - "ff1" and "ff2" are not defined at lines 167 and 267.
    They are set in harvkillop.f file (lines 257-258).
    They have to be included in modparm.f to share harvkillop.f values?
    or they have to be redefined as in harvkillop.f?

* In NCsed\_leach.f90:
  - "perc\_clyr" could be used not initialized at line 221 if `sol_nly(j)<2`

* In nrain.f:
  - "no2pcp" seems to be "no3pcp" at line 72

* In pmeas.f:
  - "rbsb" could be used not initialized at line 143
  - "flag" could be used not initialized if `a==' '` at line 210
  - "rainsb" could be used not initialized, however only if `nstep<=0`

* In pminrl2.f:
  - at line 95 a comma is necessary between "base" and "vara"
  - "ssp" could be used not initialized at line 196 if `xx<=1.e-6`

* In pothole.f:
  - "solp\_tileo" could be used not initialized at line 593 if
    `pot_vol(j)<=1.e-6` or `potvol_tile<=1.e-6`

* In potholehr.f:
  - "potflow" seems to be "potflwo" at line 447

* In readatmodep.f:
  - `momax=12*nbyr` is defined at line 65 but not used.
    It has to be "mo\_max"? but then, it overwrites the file read

* In readops.f:
  - `year = 0.` seems to be `iyear = 0` at line 98
  - "mg13" seems to be "mgt13" at line 206

* In readpnd.f:
  - "vselsetlpnd" seems to be "velsetlpnd" at line 279

* In readru.f:
  - "tck" is used but not initialized at line 79

* In readsepticbz.f:
  - at line 135 `4. e-8` seems to be `4.e-8`

* In rewind\_init.f:
  - "orig\_tnylda" is used but not initialized at line 174

* In routels.f:
  - "dstor" is used but not initialized at line 134.
    It has to be calculated as in watbal.f?
    or as in the commented line 109?
  - "latqout" and "gwqout" could be used not initialized at lines 142-143

* In rtbact.f:
  - "netwtr" could be used not initialized at line 124, however only if
    `nstep<1`

* In rthpest.f:
  - `thour=1.0` at line 183 overwrites previous "thour" calculation. It is wrong
  - "frsol" and "frsrb" could be used not initialized at lines 289-290 if
    `hrtwtr(ii)>0.001` and `hrtwtr(ii)/(idt*60)<=0.01`

* In rtpest.f:
  - `tday=1.0` at line 180 overwrites previous "tday" calculation. It is wrong

* In sched\_mgt.f:
  - `< =` seems to be `<=` at 202 line
  - "husc" and "igrow" at lines 264-265 are used but not initialized.
    "husc" has to be `phu_op(iop,ihru)` has in readmgt.f?
    "igrow" has to be `igro(ihru)` has in readmgt.f?

* In smeas.f:
  - "rabsb" could be used not initialized at line 86

* In sweep.f:
  - "fr\_curb" is used but not initialized at line 56.
    It has to be added to modparm.f to share result with sched\_mgt.f?
    or it has to be `mgt5op(nop(ihru),ihru)` as in sched\_mgt.f?

* In tmeas.f:
  - "tmxbsb" and "tmnbsb" could be used not initialized at lines 109-110

* In transfer.f:
  - "ratio", "xx" and "ratio1" could be used not initialized at lines 236, 239
    and 241 if `ihout==2`

* In wmeas.f:
  - "u10bsb" could be used not initialized at line 85

* In zero0.f:
  - "sol\_sumn03" seems to be "sol\_sumno3" at line 508

* In zero\_urbn.f:
  - "stp\_stagdis" seems to be "dtp\_stagdis" at line 84
  - "subdr\_kg" seems to be "subdr\_km" at line 149
  - "spl\_eros" is not defined at line 21, it could be "eros\_spl"?
