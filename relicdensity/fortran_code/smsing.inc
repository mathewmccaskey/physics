
	common/parms/msq,lambda,del1,del2,kap2,kap3,kap4,x0
	common/parmsmin/msqmin,lambdamin,del1min,del2min,
     .		kap2min,kap3min,kap4min,x0min
	common/parmsmax/msqmax,lambdamax,del1max,del2max,
     .		kap2max,kap3max,kap4max,x0max
	common/couplings/alphas,almz,g1,g2,sw,cw
	common/smmasses/mt,mc,mb,mmu,mw,mz,mtau,vsm,me,mu,md,ms
	common/singmass/mh,mdm,mhsm
	common/model/modelmin,modelmax,model
	common/iters/iters,itersmax,passcnt
	common/constraint/flag_fail
	common/pi/pi
	common/masssq/mhsq,mdmsq,m2h,rhiggs,thet
	common/numstate/nhiggs,ndm
	
	
	common/obliqueref/Tref,Sref,Uref
	real*8 Tref,Sref,Uref
	common/obliqueparm/Tparam,Sparam,Uparam
	real*8 Tparam,Sparam,Uparam
	common/precew/A0,gfermi,rwhat,rhohat,thetawhat,
     .	swhat,cwhat,mtref,mwref,mwval
	real*8 A0,gfermi,rwhat,rhohat,thetawhat,
     .	swhat,cwhat,mtref,mwref,mwval
	
	common/xsigma_cms/cmssig
	real*8 cmssig(2,8)
	
	common/sigscale/sigreduce
	real*8 sigreduce(2)	

        common/fermioncoup/le,lmu,ltau,lu,ld,ls,lb,lc,lt
        real*8 le,lmu,ltau,lu,ld,ls,lb,lc,lt           

	common/chisq/chisq
	real*8 chisq
	
	
	common/dmdetect/sip,sdp,oh2,oh2approx,xf,oh2unc,xfactor
	real*8 sip,sdp,oh2,oh2approx,xf,oh2unc,xfactor
		
	common/HS_lepdata/zzhmass,zzhcoup,zzh2
	common/files/allfile,nall
	common/random/idum
	
	real*8 msq,lambda,del1,del2,kap2,kap3,kap4,x0
	real*8 msqmin,lambdamin,del1min,del2min,
     .		kap2min,kap3min,kap4min,x0min
	real*8 msqmax,lambdamax,del1max,del2max,
     .		kap2max,kap3max,kap4max,x0max
	real*8 alphas,almz,g1,g2,sw,cw
	real*8 mt,mc,mb,mmu,mw,mz,mtau,vsm,me,mu,md,ms
	real*8 mh(2),mdm,mhsm

	real*8 Pi
	real*8 mhsq(2,2),mdmsq,m2h(2),rhiggs(2,2),thet
	
	real*8 zzhmass(217),zzhcoup(217),zzh2(217)

	character*80 allfile
	integer nall
	
	integer modelmin,modelmax,model
	integer iters,itersmax,passcnt
	integer flag_fail
	integer idum
	integer nhiggs,ndm


	common/v2vh/mhiv,xsectw,xsecterrw,xsect2w,
     .		xsectz,xsecterrz,xsect2z,
     .		xsectwh1,xsectwh2,xsectzh1,xsectzh2

	real*4 mhiv(500),xsectw(500),xsecterrw(500),xsect2w(500),
     .		xsectz(500),xsecterrz(500),xsect2z(500),
     .		xsectwh1,xsectwh2,xsectzh1,xsectzh2

	common/vbf/xsectvbf,xsecterrvbf,xsect2vbf,
     .		xsectvbfh1,xsectvbfh2

	real*4 xsectvbf(500),xsecterrvbf(500),xsect2vbf(500),
     .		xsectvbfh1,xsectvbfh2



