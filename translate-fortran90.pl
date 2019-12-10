@ls = split(/[\r\n]+/,`ls *.f90`);
foreach $l (@ls)
{
	system "findent -ofree < ".$l." > a";
	system "vm a ".$l;
}
@ls = split(/[\r\n]+/,`ls *.f`);
foreach $l (@ls)
{
	system "findent -ofree < ".$l." > ".$l."90";
	system "rm ".$l;
}
