.PHONY: strip clean tar
cflags = -c -O3 -Wall -Wextra $(CFLAGS)
cc = $(prefix)gfortran
strip = $(prefix)strip
sources = addh.f90 albedo.f90 allocate_parms.f90 alph.f90 anfert.f90 apex_day.f90 apply.f90 ascrv.f90 atri.f90 aunif.f90 autoirr.f90 aveval.f90 bacteria.f90 biozone.f90 bmp_det_pond.f90 bmpfixed.f90 bmpinit.f90 bmp_ri_pond.f90 bmp_sand_filter.f90 bmp_sed_pond.f90 bmp_wet_pond.f90 buffer.f90 burnop.f90 canopyint.f90 caps.f90 carbon_new.f90 carbon_zhang2.f90 cfactor.f90 chkcst.f90 clgen.f90 clicon.f90 command.f90 conapply.f90 confert.f90 crackflow.f90 crackvol.f90 curno.f90 dailycn.f90 decay.f90 depstor.f90 distrib_bmps.f90 dormant.f90 drains.f90 dstn1.f90 ee.f90 eiusle.f90 enrsb.f90 estimate_ksat.f90 etact.f90 etpot.f90 expo.f90 fert.f90 filter.f90 filtw.f90 finalbal.f90 gcycl.f90 getallo.f90 grass_wway.f90 graze.f90 grow.f90 gwmod_deep.f90 gwmod.f90 gw_no3.f90 gwnutr.f90 h2omgt_init.f90 harvestop.f90 harvgrainop.f90 harvkillop.f90 headout.f90 hhnoqual.f90 hhwatqual.f90 hmeas.f90 HQDAV.f90 hruaa.f90 hruallo.f90 hruday.f90 hrumon.f90 hrupond.f90 hrupondhr.f90 hruyr.f90 hydroinit.f90 icl.f90 impndaa.f90 impndday.f90 impnd_init.f90 impndmon.f90 impndyr.f90 irrigate.f90 irr_rch.f90 irr_res.f90 irrsub.f90 jdt.f90 killop.f90 lakeq.f90 latsed.f90 layersplit.f90 lid_cistern.f90 lid_greenroof.f90 lidinit.f90 lid_porpavement.f90 lid_raingarden.f90 lids.f90 log_normal.f90 lwqdef.f90 main.f90 modparm.f90 NCsed_leach.f90 ndenit.f90 newtillmix.f90 nfix.f90 nitvol.f90 nlch.f90 nminrl.f90 noqual.f90 npup.f90 nrain.f90 nup.f90 nuts.f90 openwth.f90 operatn.f90 orgncswat.f90 orgn.f90 origtile.f90 ovr_sed.f90 percmacro.f90 percmain.f90 percmicro.f90 pestlch.f90 pestw.f90 pesty.f90 pgen.f90 pgenhr.f90 pkq.f90 plantmod.f90 plantop.f90 pmeas.f90 pminrl2.f90 pminrl.f90 pond.f90 pondhr.f90 pothole.f90 potholehr.f90 print_hyd.f90 psed.f90 qman.f90 ran1.f90 rchaa.f90 rchday.f90 rchinit.f90 rchmon.f90 rchuse.f90 rchyr.f90 reachout.f90 readatmodep.f90 readbsn.f90 readchm.f90 readcnst.f90 readfcst.f90 readfert.f90 readfig.f90 readfile.f90 readgw.f90 readhru.f90 readinpt.f90 readlup.f90 readlwq.f90 readmgt.f90 readmon.f90 readops.f90 readpest.f90 readplant.f90 readpnd.f90 readres.f90 readrte.f90 readru.f90 readsdr.f90 readsepticbz.f90 readseptwq.f90 readsno.f90 readsol.f90 readsub.f90 readswq.f90 readtill.f90 readurban.f90 readwgn.f90 readwus.f90 readwwq.f90 readyr.f90 reccnst.f90 recday.f90 rechour.f90 recmon.f90 recyear.f90 regres.f90 resbact.f90 resetlu.f90 res.f90 reshr.f90 resinit.f90 resnut.f90 rewind_init.f90 rhgen.f90 rootfr.f90 route.f90 routels.f90 routeunit.f90 routres.f90 rsedaa.f90 rseday.f90 rsedmon.f90 rsedyr.f90 rtbact.f90 rtday.f90 rteinit.f90 rthmusk.f90 rthpest.f90 rthsed.f90 rthvsc.f90 rtmusk.f90 rtout.f90 rtpest.f90 rtsed_bagnold.f90 rtsed.f90 rtsed_kodatie.f90 rtsed_Molinas_Wu.f90 rtsed_yangsand.f90 sat_excess.f90 saveconc.f90 save.f90 sched_mgt.f90 schedule_ops.f90 sim_initday.f90 sim_inityr.f90 simulate.f90 slrgen.f90 smeas.f90 snom.f90 soil_chem.f90 soil_phys.f90 soil_write.f90 solp.f90 solt.f90 std1.f90 std2.f90 std3.f90 stdaa.f90 storeinitial.f90 structure.f90 subaa.f90 subbasin.f90 subday.f90 submon.f90 substor.f90 sub_subbasin.f90 subwq.f90 subyr.f90 sumhyd.f90 sumv.f90 surface.f90 surfst_h2o.f90 surfstor.f90 surq_daycn.f90 surq_greenampt.f90 swbl.f90 sweep.f90 swu.f90 tair.f90 tgen.f90 theta.f90 tillfactor.f90 tmeas.f90 tran.f90 transfer.f90 tstr.f90 ttcoef.f90 ttcoef_wway.f90 urban.f90 urbanhr.f90 urb_bmp.f90 varinit.f90 vbl.f90 virtual.f90 volq.f90 washp.f90 watbal.f90 water_hru.f90 watqual2.f90 watqual.f90 wattable.f90 watuse.f90 weatgn.f90 wetlan.f90 wmeas.f90 wndgen.f90 writeaa.f90 writea.f90 writed.f90 writem.f90 xmon.f90 ysed.f90 zero0.f90 zero1.f90 zero2.f90 zeroini.f90 zero_urbn.f90
objs = addh.o albedo.o allocate_parms.o alph.o anfert.o apex_day.o apply.o ascrv.o atri.o aunif.o autoirr.o aveval.o bacteria.o biozone.o bmp_det_pond.o bmpfixed.o bmpinit.o bmp_ri_pond.o bmp_sand_filter.o bmp_sed_pond.o bmp_wet_pond.o buffer.o burnop.o canopyint.o caps.o carbon_new.o carbon_zhang2.o cfactor.o chkcst.o clgen.o clicon.o command.o conapply.o confert.o crackflow.o crackvol.o curno.o dailycn.o decay.o depstor.o distrib_bmps.o dormant.o drains.o dstn1.o ee.o eiusle.o enrsb.o estimate_ksat.o etact.o etpot.o expo.o fert.o filter.o filtw.o finalbal.o gcycl.o getallo.o grass_wway.o graze.o grow.o gwmod_deep.o gwmod.o gw_no3.o gwnutr.o h2omgt_init.o harvestop.o harvgrainop.o harvkillop.o headout.o hhnoqual.o hhwatqual.o hmeas.o HQDAV.o hruaa.o hruallo.o hruday.o hrumon.o hrupond.o hrupondhr.o hruyr.o hydroinit.o icl.o impndaa.o impndday.o impnd_init.o impndmon.o impndyr.o irrigate.o irr_rch.o irr_res.o irrsub.o jdt.o killop.o lakeq.o latsed.o layersplit.o lid_cistern.o lid_greenroof.o lidinit.o lid_porpavement.o lid_raingarden.o lids.o log_normal.o lwqdef.o main.o NCsed_leach.o ndenit.o newtillmix.o nfix.o nitvol.o nlch.o nminrl.o noqual.o npup.o nrain.o nup.o nuts.o openwth.o operatn.o orgncswat.o orgn.o origtile.o ovr_sed.o percmacro.o percmain.o percmicro.o pestlch.o pestw.o pesty.o pgen.o pgenhr.o pkq.o plantmod.o plantop.o pmeas.o pminrl2.o pminrl.o pond.o pondhr.o pothole.o potholehr.o print_hyd.o psed.o qman.o ran1.o rchaa.o rchday.o rchinit.o rchmon.o rchuse.o rchyr.o reachout.o readatmodep.o readbsn.o readchm.o readcnst.o readfcst.o readfert.o readfig.o readfile.o readgw.o readhru.o readinpt.o readlup.o readlwq.o readmgt.o readmon.o readops.o readpest.o readplant.o readpnd.o readres.o readrte.o readru.o readsdr.o readsepticbz.o readseptwq.o readsno.o readsol.o readsub.o readswq.o readtill.o readurban.o readwgn.o readwus.o readwwq.o readyr.o reccnst.o recday.o rechour.o recmon.o recyear.o regres.o resbact.o resetlu.o res.o reshr.o resinit.o resnut.o rewind_init.o rhgen.o rootfr.o route.o routels.o routeunit.o routres.o rsedaa.o rseday.o rsedmon.o rsedyr.o rtbact.o rtday.o rteinit.o rthmusk.o rthpest.o rthsed.o rthvsc.o rtmusk.o rtout.o rtpest.o rtsed_bagnold.o rtsed.o rtsed_kodatie.o rtsed_Molinas_Wu.o rtsed_yangsand.o sat_excess.o saveconc.o save.o sched_mgt.o schedule_ops.o sim_initday.o sim_inityr.o simulate.o slrgen.o smeas.o snom.o soil_chem.o soil_phys.o soil_write.o solp.o solt.o std1.o std2.o std3.o stdaa.o storeinitial.o structure.o subaa.o subbasin.o subday.o submon.o substor.o sub_subbasin.o subwq.o subyr.o sumhyd.o sumv.o surface.o surfst_h2o.o surfstor.o surq_daycn.o surq_greenampt.o swbl.o sweep.o swu.o tair.o tgen.o theta.o tillfactor.o tmeas.o tran.o transfer.o tstr.o ttcoef.o ttcoef_wway.o urban.o urbanhr.o urb_bmp.o varinit.o vbl.o virtual.o volq.o washp.o watbal.o water_hru.o watqual2.o watqual.o wattable.o watuse.o weatgn.o wetlan.o wmeas.o wndgen.o writeaa.o writea.o writed.o writem.o xmon.o ysed.o zero0.o zero1.o zero2.o zeroini.o zero_urbn.o
mods = parm.mod

swat$(EXE): $(mods) $(objs)
	$(cc) $(LDFLAGS) $(objs) -o swat$(EXE)

strip:
	make
	$(strip) swat$(EXE)

clean:
	rm -rf *.mod *.o swat* latex html

tar:
	tar cJf swat.tar.xz *.f90 *.pl *.sh Makefile

parm.mod: modparm.f90 main.f90 Makefile
	$(cc) $(cflags) main.f90 -o main.o
	touch parm.mod

addh.o: addh.f90 parm.mod Makefile
	$(cc) $(cflags) addh.f90 -o addh.o

albedo.o: albedo.f90 parm.mod Makefile
	$(cc) $(cflags) albedo.f90 -o albedo.o

allocate_parms.o: allocate_parms.f90 parm.mod Makefile
	$(cc) $(cflags) allocate_parms.f90 -o allocate_parms.o

alph.o: alph.f90 parm.mod Makefile
	$(cc) $(cflags) alph.f90 -o alph.o

anfert.o: anfert.f90 parm.mod Makefile
	$(cc) $(cflags) anfert.f90 -o anfert.o

apex_day.o: apex_day.f90 parm.mod Makefile
	$(cc) $(cflags) apex_day.f90 -o apex_day.o

apply.o: apply.f90 parm.mod Makefile
	$(cc) $(cflags) apply.f90 -o apply.o

ascrv.o: ascrv.f90 Makefile
	$(cc) $(cflags) ascrv.f90 -o ascrv.o

atri.o: atri.f90 parm.mod Makefile
	$(cc) $(cflags) atri.f90 -o atri.o

aunif.o: aunif.f90 Makefile
	$(cc) $(cflags) aunif.f90 -o aunif.o

autoirr.o: autoirr.f90 parm.mod Makefile
	$(cc) $(cflags) autoirr.f90 -o autoirr.o

aveval.o: aveval.f90 Makefile
	$(cc) $(cflags) aveval.f90 -o aveval.o

bacteria.o: bacteria.f90 parm.mod Makefile
	$(cc) $(cflags) bacteria.f90 -o bacteria.o

biozone.o: biozone.f90 parm.mod Makefile
	$(cc) $(cflags) biozone.f90 -o biozone.o

bmp_det_pond.o: bmp_det_pond.f90 parm.mod Makefile
	$(cc) $(cflags) bmp_det_pond.f90 -o bmp_det_pond.o

bmpfixed.o: bmpfixed.f90 parm.mod Makefile
	$(cc) $(cflags) bmpfixed.f90 -o bmpfixed.o

bmpinit.o: bmpinit.f90 parm.mod Makefile
	$(cc) $(cflags) bmpinit.f90 -o bmpinit.o

bmp_ri_pond.o: bmp_ri_pond.f90 parm.mod Makefile
	$(cc) $(cflags) bmp_ri_pond.f90 -o bmp_ri_pond.o

bmp_sand_filter.o: bmp_sand_filter.f90 parm.mod Makefile
	$(cc) $(cflags) bmp_sand_filter.f90 -o bmp_sand_filter.o

bmp_sed_pond.o: bmp_sed_pond.f90 parm.mod Makefile
	$(cc) $(cflags) bmp_sed_pond.f90 -o bmp_sed_pond.o

bmp_wet_pond.o: bmp_wet_pond.f90 parm.mod Makefile
	$(cc) $(cflags) bmp_wet_pond.f90 -o bmp_wet_pond.o

buffer.o: buffer.f90 parm.mod Makefile
	$(cc) $(cflags) buffer.f90 -o buffer.o

burnop.o: burnop.f90 parm.mod Makefile
	$(cc) $(cflags) burnop.f90 -o burnop.o

canopyint.o: canopyint.f90 parm.mod Makefile
	$(cc) $(cflags) canopyint.f90 -o canopyint.o

caps.o: caps.f90 Makefile
	$(cc) $(cflags) caps.f90 -o caps.o

carbon_new.o: carbon_new.f90 parm.mod Makefile
	$(cc) $(cflags) carbon_new.f90 -o carbon_new.o

carbon_zhang2.o: carbon_zhang2.f90 parm.mod Makefile
	$(cc) $(cflags) -ffree-line-length-0 carbon_zhang2.f90 -o carbon_zhang2.o

cfactor.o: cfactor.f90 parm.mod Makefile
	$(cc) $(cflags) cfactor.f90 -o cfactor.o

chkcst.o: chkcst.f90 parm.mod Makefile
	$(cc) $(cflags) chkcst.f90 -o chkcst.o

clgen.o: clgen.f90 Makefile
	$(cc) $(cflags) clgen.f90 -o clgen.o

clicon.o: clicon.f90 parm.mod Makefile
	$(cc) $(cflags) clicon.f90 -o clicon.o

command.o: command.f90 parm.mod Makefile
	$(cc) $(cflags) command.f90 -o command.o

conapply.o: conapply.f90 parm.mod Makefile
	$(cc) $(cflags) conapply.f90 -o conapply.o

confert.o: confert.f90 parm.mod Makefile
	$(cc) $(cflags) confert.f90 -o confert.o

crackflow.o: crackflow.f90 parm.mod Makefile
	$(cc) $(cflags) crackflow.f90 -o crackflow.o

crackvol.o: crackvol.f90 parm.mod Makefile
	$(cc) $(cflags) crackvol.f90 -o crackvol.o

curno.o: curno.f90 parm.mod Makefile
	$(cc) $(cflags) curno.f90 -o curno.o

dailycn.o: dailycn.f90 parm.mod Makefile
	$(cc) $(cflags) dailycn.f90 -o dailycn.o

decay.o: decay.f90 parm.mod Makefile
	$(cc) $(cflags) decay.f90 -o decay.o

depstor.o: depstor.f90 parm.mod Makefile
	$(cc) $(cflags) depstor.f90 -o depstor.o

distrib_bmps.o: distrib_bmps.f90 parm.mod Makefile
	$(cc) $(cflags) distrib_bmps.f90 -o distrib_bmps.o

dormant.o: dormant.f90 parm.mod Makefile
	$(cc) $(cflags) dormant.f90 -o dormant.o

drains.o: drains.f90 parm.mod Makefile
	$(cc) $(cflags) drains.f90 -o drains.o

dstn1.o: dstn1.f90 Makefile
	$(cc) $(cflags) dstn1.f90 -o dstn1.o

ee.o: ee.f90 Makefile
	$(cc) $(cflags) ee.f90 -o ee.o

eiusle.o: eiusle.f90 parm.mod Makefile
	$(cc) $(cflags) eiusle.f90 -o eiusle.o

enrsb.o: enrsb.f90 parm.mod Makefile
	$(cc) $(cflags) enrsb.f90 -o enrsb.o

estimate_ksat.o: estimate_ksat.f90 Makefile
	$(cc) $(cflags) estimate_ksat.f90 -o estimate_ksat.o

etact.o: etact.f90 parm.mod Makefile
	$(cc) $(cflags) etact.f90 -o etact.o

etpot.o: etpot.f90 parm.mod Makefile
	$(cc) $(cflags) etpot.f90 -o etpot.o

expo.o: expo.f90 Makefile
	$(cc) $(cflags) expo.f90 -o expo.o

fert.o: fert.f90 parm.mod Makefile
	$(cc) $(cflags) fert.f90 -o fert.o

filter.o: filter.f90 parm.mod Makefile
	$(cc) $(cflags) filter.f90 -o filter.o

filtw.o: filtw.f90 parm.mod Makefile
	$(cc) $(cflags) filtw.f90 -o filtw.o

finalbal.o: finalbal.f90 parm.mod Makefile
	$(cc) $(cflags) finalbal.f90 -o finalbal.o

gcycl.o: gcycl.f90 parm.mod Makefile
	$(cc) $(cflags) gcycl.f90 -o gcycl.o

getallo.o: getallo.f90 parm.mod Makefile
	$(cc) $(cflags) getallo.f90 -o getallo.o

grass_wway.o: grass_wway.f90 parm.mod Makefile
	$(cc) $(cflags) grass_wway.f90 -o grass_wway.o

graze.o: graze.f90 parm.mod Makefile
	$(cc) $(cflags) graze.f90 -o graze.o

grow.o: grow.f90 parm.mod Makefile
	$(cc) $(cflags) grow.f90 -o grow.o

gwmod_deep.o: gwmod_deep.f90 parm.mod Makefile
	$(cc) $(cflags) gwmod_deep.f90 -o gwmod_deep.o

gwmod.o: gwmod.f90 parm.mod Makefile
	$(cc) $(cflags) gwmod.f90 -o gwmod.o

gw_no3.o: gw_no3.f90 parm.mod Makefile
	$(cc) $(cflags) gw_no3.f90 -o gw_no3.o

gwnutr.o: gwnutr.f90 parm.mod Makefile
	$(cc) $(cflags) gwnutr.f90 -o gwnutr.o

h2omgt_init.o: h2omgt_init.f90 parm.mod Makefile
	$(cc) $(cflags) h2omgt_init.f90 -o h2omgt_init.o

harvestop.o: harvestop.f90 parm.mod Makefile
	$(cc) $(cflags) harvestop.f90 -o harvestop.o

harvgrainop.o: harvgrainop.f90 parm.mod Makefile
	$(cc) $(cflags) harvgrainop.f90 -o harvgrainop.o

harvkillop.o: harvkillop.f90 parm.mod Makefile
	$(cc) $(cflags) harvkillop.f90 -o harvkillop.o

headout.o: headout.f90 parm.mod Makefile
	$(cc) $(cflags) headout.f90 -o headout.o

hhnoqual.o: hhnoqual.f90 parm.mod Makefile
	$(cc) $(cflags) hhnoqual.f90 -o hhnoqual.o

hhwatqual.o: hhwatqual.f90 parm.mod Makefile
	$(cc) $(cflags) hhwatqual.f90 -o hhwatqual.o

hmeas.o: hmeas.f90 parm.mod Makefile
	$(cc) $(cflags) hmeas.f90 -o hmeas.o

HQDAV.o: HQDAV.f90 parm.mod Makefile
	$(cc) $(cflags) HQDAV.f90 -o HQDAV.o

hruaa.o: hruaa.f90 parm.mod Makefile
	$(cc) $(cflags) hruaa.f90 -o hruaa.o

hruallo.o: hruallo.f90 parm.mod Makefile
	$(cc) $(cflags) hruallo.f90 -o hruallo.o

hruday.o: hruday.f90 parm.mod Makefile
	$(cc) $(cflags) hruday.f90 -o hruday.o

hrumon.o: hrumon.f90 parm.mod Makefile
	$(cc) $(cflags) hrumon.f90 -o hrumon.o

hrupond.o: hrupond.f90 parm.mod Makefile
	$(cc) $(cflags) hrupond.f90 -o hrupond.o

hrupondhr.o: hrupondhr.f90 parm.mod Makefile
	$(cc) $(cflags) hrupondhr.f90 -o hrupondhr.o

hruyr.o: hruyr.f90 parm.mod Makefile
	$(cc) $(cflags) hruyr.f90 -o hruyr.o

hydroinit.o: hydroinit.f90 parm.mod Makefile
	$(cc) $(cflags) hydroinit.f90 -o hydroinit.o

icl.o: icl.f90 parm.mod Makefile
	$(cc) $(cflags) icl.f90 -o icl.o

impndaa.o: impndaa.f90 parm.mod Makefile
	$(cc) $(cflags) impndaa.f90 -o impndaa.o

impndday.o: impndday.f90 parm.mod Makefile
	$(cc) $(cflags) impndday.f90 -o impndday.o

impnd_init.o: impnd_init.f90 parm.mod Makefile
	$(cc) $(cflags) impnd_init.f90 -o impnd_init.o

impndmon.o: impndmon.f90 parm.mod Makefile
	$(cc) $(cflags) impndmon.f90 -o impndmon.o

impndyr.o: impndyr.f90 parm.mod Makefile
	$(cc) $(cflags) impndyr.f90 -o impndyr.o

irrigate.o: irrigate.f90 parm.mod Makefile
	$(cc) $(cflags) irrigate.f90 -o irrigate.o

irr_rch.o: irr_rch.f90 parm.mod Makefile
	$(cc) $(cflags) irr_rch.f90 -o irr_rch.o

irr_res.o: irr_res.f90 parm.mod Makefile
	$(cc) $(cflags) irr_res.f90 -o irr_res.o

irrsub.o: irrsub.f90 parm.mod Makefile
	$(cc) $(cflags) irrsub.f90 -o irrsub.o

jdt.o: jdt.f90 Makefile
	$(cc) $(cflags) jdt.f90 -o jdt.o

killop.o: killop.f90 parm.mod Makefile
	$(cc) $(cflags) killop.f90 -o killop.o

lakeq.o: lakeq.f90 parm.mod Makefile
	$(cc) $(cflags) lakeq.f90 -o lakeq.o

latsed.o: latsed.f90 parm.mod Makefile
	$(cc) $(cflags) latsed.f90 -o latsed.o

layersplit.o: layersplit.f90 parm.mod Makefile
	$(cc) $(cflags) -ffree-line-length-0 layersplit.f90 -o layersplit.o

lid_cistern.o: lid_cistern.f90 parm.mod Makefile
	$(cc) $(cflags) lid_cistern.f90 -o lid_cistern.o

lid_greenroof.o: lid_greenroof.f90 parm.mod Makefile
	$(cc) $(cflags) lid_greenroof.f90 -o lid_greenroof.o

lidinit.o: lidinit.f90 parm.mod Makefile
	$(cc) $(cflags) lidinit.f90 -o lidinit.o

lid_porpavement.o: lid_porpavement.f90 parm.mod Makefile
	$(cc) $(cflags) lid_porpavement.f90 -o lid_porpavement.o

lid_raingarden.o: lid_raingarden.f90 parm.mod Makefile
	$(cc) $(cflags) lid_raingarden.f90 -o lid_raingarden.o

lids.o: lids.f90 parm.mod Makefile
	$(cc) $(cflags) lids.f90 -o lids.o

log_normal.o: log_normal.f90 Makefile
	$(cc) $(cflags) log_normal.f90 -o log_normal.o

lwqdef.o: lwqdef.f90 parm.mod Makefile
	$(cc) $(cflags) lwqdef.f90 -o lwqdef.o

NCsed_leach.o: NCsed_leach.f90 parm.mod Makefile
	$(cc) $(cflags) NCsed_leach.f90 -o NCsed_leach.o

ndenit.o: ndenit.f90 parm.mod Makefile
	$(cc) $(cflags) ndenit.f90 -o ndenit.o

newtillmix.o: newtillmix.f90 parm.mod Makefile
	$(cc) $(cflags) newtillmix.f90 -o newtillmix.o

nfix.o: nfix.f90 parm.mod Makefile
	$(cc) $(cflags) nfix.f90 -o nfix.o

nitvol.o: nitvol.f90 parm.mod Makefile
	$(cc) $(cflags) nitvol.f90 -o nitvol.o

nlch.o: nlch.f90 parm.mod Makefile
	$(cc) $(cflags) nlch.f90 -o nlch.o

nminrl.o: nminrl.f90 parm.mod Makefile
	$(cc) $(cflags) nminrl.f90 -o nminrl.o

noqual.o: noqual.f90 parm.mod Makefile
	$(cc) $(cflags) noqual.f90 -o noqual.o

npup.o: npup.f90 parm.mod Makefile
	$(cc) $(cflags) npup.f90 -o npup.o

nrain.o: nrain.f90 parm.mod Makefile
	$(cc) $(cflags) nrain.f90 -o nrain.o

nup.o: nup.f90 parm.mod Makefile
	$(cc) $(cflags) nup.f90 -o nup.o

nuts.o: nuts.f90 Makefile
	$(cc) $(cflags) nuts.f90 -o nuts.o

openwth.o: openwth.f90 parm.mod Makefile
	$(cc) $(cflags) openwth.f90 -o openwth.o

operatn.o: operatn.f90 parm.mod Makefile
	$(cc) $(cflags) operatn.f90 -o operatn.o

orgncswat.o: orgncswat.f90 parm.mod Makefile
	$(cc) $(cflags) orgncswat.f90 -o orgncswat.o

orgn.o: orgn.f90 parm.mod Makefile
	$(cc) $(cflags) orgn.f90 -o orgn.o

origtile.o: origtile.f90 parm.mod Makefile
	$(cc) $(cflags) origtile.f90 -o origtile.o

ovr_sed.o: ovr_sed.f90 parm.mod Makefile
	$(cc) $(cflags) ovr_sed.f90 -o ovr_sed.o

percmacro.o: percmacro.f90 parm.mod Makefile
	$(cc) $(cflags) percmacro.f90 -o percmacro.o

percmain.o: percmain.f90 parm.mod Makefile
	$(cc) $(cflags) percmain.f90 -o percmain.o

percmicro.o: percmicro.f90 parm.mod Makefile
	$(cc) $(cflags) percmicro.f90 -o percmicro.o

pestlch.o: pestlch.f90 parm.mod Makefile
	$(cc) $(cflags) pestlch.f90 -o pestlch.o

pestw.o: pestw.f90 parm.mod Makefile
	$(cc) $(cflags) pestw.f90 -o pestw.o

pesty.o: pesty.f90 parm.mod Makefile
	$(cc) $(cflags) pesty.f90 -o pesty.o

pgen.o: pgen.f90 parm.mod Makefile
	$(cc) $(cflags) pgen.f90 -o pgen.o

pgenhr.o: pgenhr.f90 parm.mod Makefile
	$(cc) $(cflags) pgenhr.f90 -o pgenhr.o

pkq.o: pkq.f90 parm.mod Makefile
	$(cc) $(cflags) pkq.f90 -o pkq.o

plantmod.o: plantmod.f90 parm.mod Makefile
	$(cc) $(cflags) plantmod.f90 -o plantmod.o

plantop.o: plantop.f90 parm.mod Makefile
	$(cc) $(cflags) plantop.f90 -o plantop.o

pmeas.o: pmeas.f90 parm.mod Makefile
	$(cc) $(cflags) pmeas.f90 -o pmeas.o

pminrl2.o: pminrl2.f90 parm.mod Makefile
	$(cc) $(cflags) pminrl2.f90 -o pminrl2.o

pminrl.o: pminrl.f90 parm.mod Makefile
	$(cc) $(cflags) pminrl.f90 -o pminrl.o

pond.o: pond.f90 parm.mod Makefile
	$(cc) $(cflags) pond.f90 -o pond.o

pondhr.o: pondhr.f90 parm.mod Makefile
	$(cc) $(cflags) pondhr.f90 -o pondhr.o

pothole.o: pothole.f90 parm.mod Makefile
	$(cc) $(cflags) pothole.f90 -o pothole.o

potholehr.o: potholehr.f90 parm.mod Makefile
	$(cc) $(cflags) potholehr.f90 -o potholehr.o

print_hyd.o: print_hyd.f90 parm.mod Makefile
	$(cc) $(cflags) print_hyd.f90 -o print_hyd.o

psed.o: psed.f90 parm.mod Makefile
	$(cc) $(cflags) psed.f90 -o psed.o

qman.o: qman.f90 Makefile
	$(cc) $(cflags) qman.f90 -o qman.o

ran1.o: ran1.f90 Makefile
	$(cc) $(cflags) ran1.f90 -o ran1.o

rchaa.o: rchaa.f90 parm.mod Makefile
	$(cc) $(cflags) rchaa.f90 -o rchaa.o

rchday.o: rchday.f90 parm.mod Makefile
	$(cc) $(cflags) rchday.f90 -o rchday.o

rchinit.o: rchinit.f90 parm.mod Makefile
	$(cc) $(cflags) rchinit.f90 -o rchinit.o

rchmon.o: rchmon.f90 parm.mod Makefile
	$(cc) $(cflags) rchmon.f90 -o rchmon.o

rchuse.o: rchuse.f90 parm.mod Makefile
	$(cc) $(cflags) rchuse.f90 -o rchuse.o

rchyr.o: rchyr.f90 parm.mod Makefile
	$(cc) $(cflags) rchyr.f90 -o rchyr.o

reachout.o: reachout.f90 parm.mod Makefile
	$(cc) $(cflags) reachout.f90 -o reachout.o

readatmodep.o: readatmodep.f90 parm.mod Makefile
	$(cc) $(cflags) readatmodep.f90 -o readatmodep.o

readbsn.o: readbsn.f90 parm.mod Makefile
	$(cc) $(cflags) readbsn.f90 -o readbsn.o

readchm.o: readchm.f90 parm.mod Makefile
	$(cc) $(cflags) readchm.f90 -o readchm.o

readcnst.o: readcnst.f90 parm.mod Makefile
	$(cc) $(cflags) readcnst.f90 -o readcnst.o

readfcst.o: readfcst.f90 parm.mod Makefile
	$(cc) $(cflags) readfcst.f90 -o readfcst.o

readfert.o: readfert.f90 parm.mod Makefile
	$(cc) $(cflags) readfert.f90 -o readfert.o

readfig.o: readfig.f90 parm.mod Makefile
	$(cc) $(cflags) readfig.f90 -o readfig.o

readfile.o: readfile.f90 parm.mod Makefile
	$(cc) $(cflags) readfile.f90 -o readfile.o

readgw.o: readgw.f90 parm.mod Makefile
	$(cc) $(cflags) readgw.f90 -o readgw.o

readhru.o: readhru.f90 parm.mod Makefile
	$(cc) $(cflags) readhru.f90 -o readhru.o

readinpt.o: readinpt.f90 parm.mod Makefile
	$(cc) $(cflags) readinpt.f90 -o readinpt.o

readlup.o: readlup.f90 parm.mod Makefile
	$(cc) $(cflags) readlup.f90 -o readlup.o

readlwq.o: readlwq.f90 parm.mod Makefile
	$(cc) $(cflags) readlwq.f90 -o readlwq.o

readmgt.o: readmgt.f90 parm.mod Makefile
	$(cc) $(cflags) readmgt.f90 -o readmgt.o

readmon.o: readmon.f90 parm.mod Makefile
	$(cc) $(cflags) readmon.f90 -o readmon.o

readops.o: readops.f90 parm.mod Makefile
	$(cc) $(cflags) readops.f90 -o readops.o

readpest.o: readpest.f90 parm.mod Makefile
	$(cc) $(cflags) readpest.f90 -o readpest.o

readplant.o: readplant.f90 parm.mod Makefile
	$(cc) $(cflags) readplant.f90 -o readplant.o

readpnd.o: readpnd.f90 parm.mod Makefile
	$(cc) $(cflags) readpnd.f90 -o readpnd.o

readres.o: readres.f90 parm.mod Makefile
	$(cc) $(cflags) readres.f90 -o readres.o

readrte.o: readrte.f90 parm.mod Makefile
	$(cc) $(cflags) readrte.f90 -o readrte.o

readru.o: readru.f90 parm.mod Makefile
	$(cc) $(cflags) readru.f90 -o readru.o

readsdr.o: readsdr.f90 parm.mod Makefile
	$(cc) $(cflags) readsdr.f90 -o readsdr.o

readsepticbz.o: readsepticbz.f90 parm.mod Makefile
	$(cc) $(cflags) readsepticbz.f90 -o readsepticbz.o

readseptwq.o: readseptwq.f90 parm.mod Makefile
	$(cc) $(cflags) readseptwq.f90 -o readseptwq.o

readsno.o: readsno.f90 parm.mod Makefile
	$(cc) $(cflags) readsno.f90 -o readsno.o

readsol.o: readsol.f90 parm.mod Makefile
	$(cc) $(cflags) readsol.f90 -o readsol.o

readsub.o: readsub.f90 parm.mod Makefile
	$(cc) $(cflags) readsub.f90 -o readsub.o

readswq.o: readswq.f90 parm.mod Makefile
	$(cc) $(cflags) readswq.f90 -o readswq.o

readtill.o: readtill.f90 parm.mod Makefile
	$(cc) $(cflags) readtill.f90 -o readtill.o

readurban.o: readurban.f90 parm.mod Makefile
	$(cc) $(cflags) readurban.f90 -o readurban.o

readwgn.o: readwgn.f90 parm.mod Makefile
	$(cc) $(cflags) readwgn.f90 -o readwgn.o

readwus.o: readwus.f90 parm.mod Makefile
	$(cc) $(cflags) readwus.f90 -o readwus.o

readwwq.o: readwwq.f90 parm.mod Makefile
	$(cc) $(cflags) readwwq.f90 -o readwwq.o

readyr.o: readyr.f90 parm.mod Makefile
	$(cc) $(cflags) readyr.f90 -o readyr.o

reccnst.o: reccnst.f90 parm.mod Makefile
	$(cc) $(cflags) reccnst.f90 -o reccnst.o

recday.o: recday.f90 parm.mod Makefile
	$(cc) $(cflags) recday.f90 -o recday.o

rechour.o: rechour.f90 parm.mod Makefile
	$(cc) $(cflags) rechour.f90 -o rechour.o

recmon.o: recmon.f90 parm.mod Makefile
	$(cc) $(cflags) recmon.f90 -o recmon.o

recyear.o: recyear.f90 parm.mod Makefile
	$(cc) $(cflags) recyear.f90 -o recyear.o

regres.o: regres.f90 parm.mod Makefile
	$(cc) $(cflags) regres.f90 -o regres.o

resbact.o: resbact.f90 parm.mod Makefile
	$(cc) $(cflags) resbact.f90 -o resbact.o

resetlu.o: resetlu.f90 parm.mod Makefile
	$(cc) $(cflags) resetlu.f90 -o resetlu.o

res.o: res.f90 parm.mod Makefile
	$(cc) $(cflags) res.f90 -o res.o

reshr.o: reshr.f90 parm.mod Makefile
	$(cc) $(cflags) reshr.f90 -o reshr.o

resinit.o: resinit.f90 parm.mod Makefile
	$(cc) $(cflags) resinit.f90 -o resinit.o

resnut.o: resnut.f90 parm.mod Makefile
	$(cc) $(cflags) resnut.f90 -o resnut.o

rewind_init.o: rewind_init.f90 parm.mod Makefile
	$(cc) $(cflags) rewind_init.f90 -o rewind_init.o

rhgen.o: rhgen.f90 parm.mod Makefile
	$(cc) $(cflags) rhgen.f90 -o rhgen.o

rootfr.o: rootfr.f90 parm.mod Makefile
	$(cc) $(cflags) rootfr.f90 -o rootfr.o

route.o: route.f90 parm.mod Makefile
	$(cc) $(cflags) route.f90 -o route.o

routels.o: routels.f90 parm.mod Makefile
	$(cc) $(cflags) routels.f90 -o routels.o

routeunit.o: routeunit.f90 parm.mod Makefile
	$(cc) $(cflags) routeunit.f90 -o routeunit.o

routres.o: routres.f90 parm.mod Makefile
	$(cc) $(cflags) routres.f90 -o routres.o

rsedaa.o: rsedaa.f90 parm.mod Makefile
	$(cc) $(cflags) rsedaa.f90 -o rsedaa.o

rseday.o: rseday.f90 parm.mod Makefile
	$(cc) $(cflags) rseday.f90 -o rseday.o

rsedmon.o: rsedmon.f90 parm.mod Makefile
	$(cc) $(cflags) rsedmon.f90 -o rsedmon.o

rsedyr.o: rsedyr.f90 parm.mod Makefile
	$(cc) $(cflags) rsedyr.f90 -o rsedyr.o

rtbact.o: rtbact.f90 parm.mod Makefile
	$(cc) $(cflags) rtbact.f90 -o rtbact.o

rtday.o: rtday.f90 parm.mod Makefile
	$(cc) $(cflags) rtday.f90 -o rtday.o

rteinit.o: rteinit.f90 parm.mod Makefile
	$(cc) $(cflags) rteinit.f90 -o rteinit.o

rthmusk.o: rthmusk.f90 parm.mod Makefile
	$(cc) $(cflags) rthmusk.f90 -o rthmusk.o

rthpest.o: rthpest.f90 parm.mod Makefile
	$(cc) $(cflags) rthpest.f90 -o rthpest.o

rthsed.o: rthsed.f90 parm.mod Makefile
	$(cc) $(cflags) rthsed.f90 -o rthsed.o

rthvsc.o: rthvsc.f90 parm.mod Makefile
	$(cc) $(cflags) rthvsc.f90 -o rthvsc.o

rtmusk.o: rtmusk.f90 parm.mod Makefile
	$(cc) $(cflags) rtmusk.f90 -o rtmusk.o

rtout.o: rtout.f90 parm.mod Makefile
	$(cc) $(cflags) rtout.f90 -o rtout.o

rtpest.o: rtpest.f90 parm.mod Makefile
	$(cc) $(cflags) rtpest.f90 -o rtpest.o

rtsed_bagnold.o: rtsed_bagnold.f90 parm.mod Makefile
	$(cc) $(cflags) rtsed_bagnold.f90 -o rtsed_bagnold.o

rtsed.o: rtsed.f90 parm.mod Makefile
	$(cc) $(cflags) rtsed.f90 -o rtsed.o

rtsed_kodatie.o: rtsed_kodatie.f90 parm.mod Makefile
	$(cc) $(cflags) rtsed_kodatie.f90 -o rtsed_kodatie.o

rtsed_Molinas_Wu.o: rtsed_Molinas_Wu.f90 parm.mod Makefile
	$(cc) $(cflags) rtsed_Molinas_Wu.f90 -o rtsed_Molinas_Wu.o

rtsed_yangsand.o: rtsed_yangsand.f90 parm.mod Makefile
	$(cc) $(cflags) rtsed_yangsand.f90 -o rtsed_yangsand.o

sat_excess.o: sat_excess.f90 parm.mod Makefile
	$(cc) $(cflags) sat_excess.f90 -o sat_excess.o

saveconc.o: saveconc.f90 parm.mod Makefile
	$(cc) $(cflags) saveconc.f90 -o saveconc.o

save.o: save.f90 parm.mod Makefile
	$(cc) $(cflags) save.f90 -o save.o

sched_mgt.o: sched_mgt.f90 parm.mod Makefile
	$(cc) $(cflags) sched_mgt.f90 -o sched_mgt.o

schedule_ops.o: schedule_ops.f90 parm.mod Makefile
	$(cc) $(cflags) schedule_ops.f90 -o schedule_ops.o

sim_initday.o: sim_initday.f90 parm.mod Makefile
	$(cc) $(cflags) sim_initday.f90 -o sim_initday.o

sim_inityr.o: sim_inityr.f90 parm.mod Makefile
	$(cc) $(cflags) sim_inityr.f90 -o sim_inityr.o

simulate.o: simulate.f90 parm.mod Makefile
	$(cc) $(cflags) simulate.f90 -o simulate.o

slrgen.o: slrgen.f90 parm.mod Makefile
	$(cc) $(cflags) slrgen.f90 -o slrgen.o

smeas.o: smeas.f90 parm.mod Makefile
	$(cc) $(cflags) smeas.f90 -o smeas.o

snom.o: snom.f90 parm.mod Makefile
	$(cc) $(cflags) snom.f90 -o snom.o

soil_chem.o: soil_chem.f90 parm.mod Makefile
	$(cc) $(cflags) soil_chem.f90 -o soil_chem.o

soil_phys.o: soil_phys.f90 parm.mod Makefile
	$(cc) $(cflags) soil_phys.f90 -o soil_phys.o

soil_write.o: soil_write.f90 parm.mod Makefile
	$(cc) $(cflags) soil_write.f90 -o soil_write.o

solp.o: solp.f90 parm.mod Makefile
	$(cc) $(cflags) solp.f90 -o solp.o

solt.o: solt.f90 parm.mod Makefile
	$(cc) $(cflags) solt.f90 -o solt.o

std1.o: std1.f90 parm.mod Makefile
	$(cc) $(cflags) std1.f90 -o std1.o

std2.o: std2.f90 parm.mod Makefile
	$(cc) $(cflags) std2.f90 -o std2.o

std3.o: std3.f90 parm.mod Makefile
	$(cc) $(cflags) std3.f90 -o std3.o

stdaa.o: stdaa.f90 parm.mod Makefile
	$(cc) $(cflags) stdaa.f90 -o stdaa.o

storeinitial.o: storeinitial.f90 parm.mod Makefile
	$(cc) $(cflags) storeinitial.f90 -o storeinitial.o

structure.o: structure.f90 parm.mod Makefile
	$(cc) $(cflags) structure.f90 -o structure.o

subaa.o: subaa.f90 parm.mod Makefile
	$(cc) $(cflags) subaa.f90 -o subaa.o

subbasin.o: subbasin.f90 parm.mod Makefile
	$(cc) $(cflags) subbasin.f90 -o subbasin.o

subday.o: subday.f90 parm.mod Makefile
	$(cc) $(cflags) subday.f90 -o subday.o

submon.o: submon.f90 parm.mod Makefile
	$(cc) $(cflags) submon.f90 -o submon.o

substor.o: substor.f90 parm.mod Makefile
	$(cc) $(cflags) substor.f90 -o substor.o

sub_subbasin.o: sub_subbasin.f90 parm.mod Makefile
	$(cc) $(cflags) sub_subbasin.f90 -o sub_subbasin.o

subwq.o: subwq.f90 parm.mod Makefile
	$(cc) $(cflags) subwq.f90 -o subwq.o

subyr.o: subyr.f90 parm.mod Makefile
	$(cc) $(cflags) subyr.f90 -o subyr.o

sumhyd.o: sumhyd.f90 parm.mod Makefile
	$(cc) $(cflags) sumhyd.f90 -o sumhyd.o

sumv.o: sumv.f90 parm.mod Makefile
	$(cc) $(cflags) sumv.f90 -o sumv.o

surface.o: surface.f90 parm.mod Makefile
	$(cc) $(cflags) surface.f90 -o surface.o

surfst_h2o.o: surfst_h2o.f90 parm.mod Makefile
	$(cc) $(cflags) surfst_h2o.f90 -o surfst_h2o.o

surfstor.o: surfstor.f90 parm.mod Makefile
	$(cc) $(cflags) surfstor.f90 -o surfstor.o

surq_daycn.o: surq_daycn.f90 parm.mod Makefile
	$(cc) $(cflags) surq_daycn.f90 -o surq_daycn.o

surq_greenampt.o: surq_greenampt.f90 parm.mod Makefile
	$(cc) $(cflags) surq_greenampt.f90 -o surq_greenampt.o

swbl.o: swbl.f90 parm.mod Makefile
	$(cc) $(cflags) swbl.f90 -o swbl.o

sweep.o: sweep.f90 parm.mod Makefile
	$(cc) $(cflags) sweep.f90 -o sweep.o

swu.o: swu.f90 parm.mod Makefile
	$(cc) $(cflags) swu.f90 -o swu.o

tair.o: tair.f90 parm.mod Makefile
	$(cc) $(cflags) tair.f90 -o tair.o

tgen.o: tgen.f90 parm.mod Makefile
	$(cc) $(cflags) tgen.f90 -o tgen.o

theta.o: theta.f90 Makefile
	$(cc) $(cflags) theta.f90 -o theta.o

tillfactor.o: tillfactor.f90 parm.mod Makefile
	$(cc) $(cflags) tillfactor.f90 -o tillfactor.o

tmeas.o: tmeas.f90 parm.mod Makefile
	$(cc) $(cflags) tmeas.f90 -o tmeas.o

tran.o: tran.f90 parm.mod Makefile
	$(cc) $(cflags) tran.f90 -o tran.o

transfer.o: transfer.f90 parm.mod Makefile
	$(cc) $(cflags) transfer.f90 -o transfer.o

tstr.o: tstr.f90 parm.mod Makefile
	$(cc) $(cflags) tstr.f90 -o tstr.o

ttcoef.o: ttcoef.f90 parm.mod Makefile
	$(cc) $(cflags) ttcoef.f90 -o ttcoef.o

ttcoef_wway.o: ttcoef_wway.f90 parm.mod Makefile
	$(cc) $(cflags) ttcoef_wway.f90 -o ttcoef_wway.o

urban.o: urban.f90 parm.mod Makefile
	$(cc) $(cflags) urban.f90 -o urban.o

urbanhr.o: urbanhr.f90 parm.mod Makefile
	$(cc) $(cflags) urbanhr.f90 -o urbanhr.o

urb_bmp.o: urb_bmp.f90 parm.mod Makefile
	$(cc) $(cflags) urb_bmp.f90 -o urb_bmp.o

varinit.o: varinit.f90 parm.mod Makefile
	$(cc) $(cflags) varinit.f90 -o varinit.o

vbl.o: vbl.f90 Makefile
	$(cc) $(cflags) vbl.f90 -o vbl.o

virtual.o: virtual.f90 parm.mod Makefile
	$(cc) $(cflags) virtual.f90 -o virtual.o

volq.o: volq.f90 parm.mod Makefile
	$(cc) $(cflags) volq.f90 -o volq.o

washp.o: washp.f90 parm.mod Makefile
	$(cc) $(cflags) washp.f90 -o washp.o

watbal.o: watbal.f90 parm.mod Makefile
	$(cc) $(cflags) watbal.f90 -o watbal.o

water_hru.o: water_hru.f90 parm.mod Makefile
	$(cc) $(cflags) water_hru.f90 -o water_hru.o

watqual2.o: watqual2.f90 parm.mod Makefile
	$(cc) $(cflags) watqual2.f90 -o watqual2.o

watqual.o: watqual.f90 parm.mod Makefile
	$(cc) $(cflags) watqual.f90 -o watqual.o

wattable.o: wattable.f90 parm.mod Makefile
	$(cc) $(cflags) wattable.f90 -o wattable.o

watuse.o: watuse.f90 parm.mod Makefile
	$(cc) $(cflags) watuse.f90 -o watuse.o

weatgn.o: weatgn.f90 parm.mod Makefile
	$(cc) $(cflags) weatgn.f90 -o weatgn.o

wetlan.o: wetlan.f90 parm.mod Makefile
	$(cc) $(cflags) wetlan.f90 -o wetlan.o

wmeas.o: wmeas.f90 parm.mod Makefile
	$(cc) $(cflags) wmeas.f90 -o wmeas.o

wndgen.o: wndgen.f90 parm.mod Makefile
	$(cc) $(cflags) wndgen.f90 -o wndgen.o

writeaa.o: writeaa.f90 parm.mod Makefile
	$(cc) $(cflags) writeaa.f90 -o writeaa.o

writea.o: writea.f90 parm.mod Makefile
	$(cc) $(cflags) writea.f90 -o writea.o

writed.o: writed.f90 parm.mod Makefile
	$(cc) $(cflags) writed.f90 -o writed.o

writem.o: writem.f90 parm.mod Makefile
	$(cc) $(cflags) writem.f90 -o writem.o

xmon.o: xmon.f90 parm.mod Makefile
	$(cc) $(cflags) xmon.f90 -o xmon.o

ysed.o: ysed.f90 parm.mod Makefile
	$(cc) $(cflags) ysed.f90 -o ysed.o

zero0.o: zero0.f90 parm.mod Makefile
	$(cc) $(cflags) zero0.f90 -o zero0.o

zero1.o: zero1.f90 parm.mod Makefile
	$(cc) $(cflags) zero1.f90 -o zero1.o

zero2.o: zero2.f90 parm.mod Makefile
	$(cc) $(cflags) zero2.f90 -o zero2.o

zeroini.o: zeroini.f90 parm.mod Makefile
	$(cc) $(cflags) zeroini.f90 -o zeroini.o

zero_urbn.o: zero_urbn.f90 parm.mod Makefile
	$(cc) $(cflags) zero_urbn.f90 -o zero_urbn.o

latex/refman.pdf: $(sources) bib.bib README.md Makefile Doxyfile
	doxygen
	cd latex; make