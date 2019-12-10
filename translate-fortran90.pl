@ls = split(/[\r\n]+/,`ls *.f90`);
foreach $l (@ls)
{
	system "findent -ofree < ".$l." > a";
	system "sed -i \"s/\t/ /g\" a";
	system "mv a ".$l;
}
@ls = split(/[\r\n]+/,`ls *.f`);
foreach $l (@ls)
{
	system "findent -ofree < ".$l." > ".$l."90";
	system "sed -i \"s/\t/ /g\" ".$l."90";
	system "rm ".$l;
}
