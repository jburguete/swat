@sources = split( /[\r\n+]/, `ls *.f90`);
open( DATA, ">Makefile" );
print DATA ".PHONY: strip clean tar\n".
	"cflags = -c -O3 -Wall -Wextra \$(CFLAGS)\n".
	"cc = \$(prefix)gfortran\n".
	"strip = \$(prefix)strip\n".
	"sources =";
foreach $source (@sources)
{
	print DATA " ".$source;
}
print DATA "\nobjs =";
foreach $source (@sources)
{
	if ( $source ne "modparm.f90" )
	{	
		$obj = $source;
		$obj =~ s/\.f90/\.o/g;
		print DATA " ".$obj;
	}
}
print DATA "\nmods = parm.mod\n".
	"\nswat\$(EXE): \$(mods) \$(objs)\n".
	"\t\$(cc) -O3 \$(LDFLAGS) \$(objs) -o swat\$(EXE)\n".
	"\nstrip:\n\tmake\n\t\$(strip) swat\$(EXE)\n".
	"\nclean:\n\trm -rf *.mod *.o swat* latex html\n".
	"\ntar:\n\ttar cJf swat.tar.xz *.f90 *.pl *.sh Makefile\n".
	"\nparm.mod: modparm.f90 main.f90 Makefile\n".
	"\t\$(cc) \$(cflags) main.f90 -o main.o\n\ttouch parm.mod\n";
foreach $source (@sources)
{
	if ($source eq "main.f90")
	{
		print DATA "\nmain.o:\n".
			"\t\$(cc) \$(cflags) ".$source." -o main.o\n";
	}
	elsif (($source ne "modparm.f90"))
	{
		$obj = $source;
		$obj =~ s/\.f90/\.o/g;
		print DATA "\n".$obj.": ".$source." parm.mod\n".
			"\t\$(cc) \$(cflags) ".$source." -o ".$obj."\n";
	}
}
print DATA
	"\nlatex/refman.pdf: \$(sources) bib.bib README.md Makefile Doxyfile\n".
	"\tdoxygen\n\tcd latex; make";
close DATA;
@sources = (
	"ascrv\.f90",
	"aunif\.f90",
	"aveval\.f90",
	"caps\.f90",
	"clgen\.f90",
	"dstn1\.f90",
	"ee\.f90",
	"estimate_ksat\.f90",
	"expo\.f90",
	"jdt\.f90",
	"log_normal\.f90",
	"nuts\.f90",
	"qman\.f90",
	"ran1\.f90",
	"theta\.f90",
	"vbl\.f90"
);
foreach $source (@sources)
{
	system "sed -i \"s/".$source." parm\.mod/".$source."/g\" Makefile";
}
system "sed -i \"s/) carbon_z/) -ffree-line-length-0 carbon_z/g\" Makefile"; 
system "sed -i \"s/) lay/) -ffree-line-length-0 lay/g\" Makefile"; 
