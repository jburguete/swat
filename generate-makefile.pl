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
	"\t\$(cc) \$(LDFLAGS) \$(objs) -o swat\$(EXE)\n".
	"\nstrip:\n\tmake\n\t\$(strip) swat\$(EXE)\n".
	"\nclean:\n\trm *.o swat*\n".
	"\ntar:\n\ttar cJf swat.tar.xz *.f90 *.pl *.sh Makefile\n".
	"\nparm.mod: modparm.f90 main.f90 Makefile\n".
	"\t\$(cc) \$(cflags) main.f90 -o main.o\n\ttouch parm.mod\n";
foreach $source (@sources)
{
	if (($source ne "modparm.f90") and ($source ne "main.f90"))
	{
		$obj = $source;
		$obj =~ s/\.f90/\.o/g;
		print DATA "\n".$obj.": ".$source." parm.mod Makefile\n".
			"\t\$(cc) \$(cflags) ".$source." -o ".$obj."\n"
	}
}
