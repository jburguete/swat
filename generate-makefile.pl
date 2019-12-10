@sources = split( /[\r\n+]/, `ls *.f90`);
open( DATA, ">Makefile" );
print DATA ".PHONY: clean tar\n".
	"cflags = -c -O3 -Wall -Wextra \$(CFLAGS)\n".
	"cc = \$(prefix)gfortran \$(LDFLAGS)\n".
	"strip = \$(prefix)strip\n".
	"sources =";
foreach $source (@sources)
{
	print DATA " ".$source;
}
print DATA "\nobjs =";
foreach $source (@sources)
{
	$obj = $source;
	$obj =~ s/\.f90/\.o/g;
	print DATA " ".$obj;
}
print DATA
	"\n\nswat\$(EXE): \$(objs)\n\t\$(cc) -static \$(objs) -o swat\$(EXE)\n".
	"\t\$(strip) swat\$(EXE)\n".
	"\nclean:\n\trm *.o swat*\n".
	"\ntar:\n\ttar cJf swat.tar.xz *.f90 *.pl *.sh Makefile\n".
	"\nmodparm.o: modparm.f90 Makefile\n".
	"\t\$(cc) \$(cflags) modparm.f90 -o modparm.o\n";
foreach $source (@sources)
{
	if ($source ne "modparm.f90")
	{
		$obj = $source;
		$obj =~ s/\.f90/\.o/g;
		print DATA "\n".$obj.": ".$source." modparm.o\n".
			"\t\$(cc) \$(cflags) ".$source." -o ".$obj."\n"
	}
}
