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
system "sed -i \"s/include 'modparm.f'/include 'modparm.f90'/g\" main.f90";
system "sed -i \"s/use parm/use parm, only: Aunif/g\" atri.f90";
system "sed -i \"s/bmp_sed _pond/bmp_sed_pond/g\" bmp_sed_pond.f90";
system "sed -i \"s/USE PARM/USE PARM, ONLY: ch_w, ch_n, qcap, chxp, rchx, rcss, ch_s, chxa/g\" HQDAV.f90";
system "sed -i \"s/use parm/use parm, only: ihru, iseptic, sol_awc, sol_cbn, sol_bd, sol_cal, sol_clay, sol_ec, sol_k, sol_mc, sol_mn, sol_mp, sol_n, sol_nly, sol_no3, sol_orgn, sol_orgp, sol_ph, sol_rock, sol_sand, sol_silt, sol_solp, sol_z/g\" layersplit.f90";
system "sed -i \"s/use parm/use parm, only: cdn, sol_cbn, sol_no3/g\" ndenit.f90";
system "sed -i \"s/base vara/base, vara/g\" pminrl2.f90";
system "sed -i \"s/ALOG/LOG/g\" pminrl2.f90";
system "sed -i \"s/ e-/e-/g\" readsepticbz.f90";
system "sed -i \"s/use parm/use parm, only: fimp, hru_km, hru_sub, ihru, ireg, precipday, urblu/g\" regres.f90";
system "sed -i \"s/use parm/use parm, only: rch_dakm, rchaao, subgis, subtot/g\" rsedaa.f90";
system "sed -i \"s/< =/<=/g\" sched_mgt.f90";
system "sed -i \"s/use parm/use parm, only: tmn, tmp_hi, tmp_lo, tmx/g\" tair.f90";
system "sed -i \"s/use parm//g\" vbl.f90";
