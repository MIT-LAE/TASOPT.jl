
     
      function airfun(cl,tau,Mach,
	idim,jdim,kdim,ldim,
	nAMa,nAcl,nAtau,nAfun,
	AMa,Acl,Atau,ARe,
	A,
	A_M,
	A_cl,
	A_tau,
	A_M_cl,
	A_M_tau,
	A_cl_tau,
	A_M_cl_tau,
	cdf, cdp, cdw, cm)
#-----------------------------------------------------------
#     Evaluates tri-cubic spline functions 
#
#      cdf(Mach,cl,tau) = A(...1)
#      cdp(Mach,cl,tau) = A(...2)
#       cm(Mach,cl,tau) = A(...3)
#
#     at specified paramater values Mach,cl,tau.
#-----------------------------------------------------------
#      implicit none
#
#      real cl,tau,Mach
#      integer idim,jdim,kdim,ldim
#      integer nAMa, nAcl, nAtau,nAfun
#      real AMa(idim), Acl(jdim), Atau(kdim), ARe
#      real       A(idim,jdim,kdim,ldim),
#	A_M(idim,jdim,kdim,ldim),
#	A_cl(idim,jdim,kdim,ldim),
#	A_tau(idim,jdim,kdim,ldim),
#	A_M_cl(idim,jdim,kdim,ldim),
#	A_M_tau(idim,jdim,kdim,ldim),
#	A_cl_tau(idim,jdim,kdim,ldim),
#	A_M_cl_tau(idim,jdim,kdim,ldim)
#      real cdf, cdp, cdw, cm
#
#      real Mamin, Mamax, dMa, tMa,
#	clmin, clmax, dcl, tcl,
#	taumin, taumax, dtau, ttau
##
#      integer ilow,imid,
#	jlow,jmid,
#	klow,kmid
#
#      integer i,j,k,l
#
#
#      integer io, im, id,
#	jo, jm, jd,
#	ko, km, kd
#      real fxm, fxo
#      real    Ai(2,2),
#	Ai_cl(2,2),
#	Ai_tau(2,2),
#	Ai_cl_tau(2,2)
#      real    Aij(2),
#	Aij_tau(2)
#
#      integer jx,kx,lx
#      parameter (lx=3)
#      real Aijk(lx)
#
#      if(nAfun>lx) 
#       write(*,*) 'AIRFUN: Local array overflow. Increase lx to', nAfun
#       stop
#      end
#
#---- find Ma interval i-1..i
      ilow = 1
      io = nAMa
#
   10 if(io-ilow <= 1) go to 11
#
      imid = (io+ilow)/2
      if(Mach < AMa(imid)) 
       io = imid
      else
       ilow = imid
      end
      go to 10
#
   11 continue
      im = io-1
      dMa = AMa(io) - AMa(im)
      tMa = (Mach - AMa(im)) / dMa


#---- find cl interval j-1..j
      jlow = 1
      jo = nAcl
#
   20 if(jo-jlow <= 1) go to 21
#
      jmid = (jo+jlow)/2
      if(cl < Acl(jmid)) 
       jo = jmid
      else
       jlow = jmid
      end
      go to 20
#
   21 continue
      jm = jo-1
      dcl = Acl(jo) - Acl(jm)
      tcl = (cl - Acl(jm)) / dcl


#---- find tau interval k-1..k
      klow = 1
      ko = nAtau
#
   30 if(ko-klow <= 1) go to 31
#
      kmid = (ko+klow)/2
      if(tau < Atau(kmid)) 
       ko = kmid
      else
       klow = kmid
      end
      go to 30
#
   31 continue
      km = ko-1
      dtau = Atau(ko) - Atau(km)
      ttau = (tau - Atau(km)) / dtau


      for l = 1: lx
        Aijk(l) = 0.
      end

#---- evaluate tri-cubic spline in i,j,k box
      for l = 1: nAfun

      for jd = 1: 2
      for kd = 1: 2
        j = jo + jd-2
        k = ko + kd-2

        fxm = dMa*A_M(im,j,k,l) - A(io,j,k,l) + A(im,j,k,l)
        fxo = dMa*A_M(io,j,k,l) - A(io,j,k,l) + A(im,j,k,l)
        Ai(jd,kd) = tMa *A(io,j,k,l)
	+ (1.0-tMa)*A(im,j,k,l) 
	+ tMa*(1.0-tMa)*((1.0-tMa)*fxm - tMa*fxo)

        fxm = dMa*A_M_cl(im,j,k,l)
	- A_cl(io,j,k,l)
	+ A_cl(im,j,k,l)
        fxo = dMa*A_M_cl(io,j,k,l)
	- A_cl(io,j,k,l)
	+ A_cl(im,j,k,l)
        Ai_cl(jd,kd) = tMa *A_cl(io,j,k,l)
	+ (1.0-tMa)*A_cl(im,j,k,l) 
	+ tMa*(1.0-tMa)*((1.0-tMa)*fxm - tMa*fxo)

        fxm = dMa*A_M_tau(im,j,k,l)
	- A_tau(io,j,k,l)
	+ A_tau(im,j,k,l)
        fxo = dMa*A_M_tau(io,j,k,l)
	- A_tau(io,j,k,l)
	+ A_tau(im,j,k,l)
        Ai_tau(jd,kd) = tMa *A_tau(io,j,k,l)
	+ (1.0-tMa)*A_tau(im,j,k,l) 
	+ tMa*(1.0-tMa)*((1.0-tMa)*fxm - tMa*fxo)

        fxm = dMa*A_M_cl_tau(im,j,k,l)
	- A_cl_tau(io,j,k,l)
	+ A_cl_tau(im,j,k,l)
        fxo = dMa*A_M_cl_tau(io,j,k,l)
	- A_cl_tau(io,j,k,l)
	+ A_cl_tau(im,j,k,l)
        Ai_cl_tau(jd,kd) = tMa *A_cl_tau(io,j,k,l)
	+ (1.0-tMa)*A_cl_tau(im,j,k,l) 
	+ tMa*(1.0-tMa)*((1.0-tMa)*fxm - tMa*fxo)
      end
      end

      for kd = 1: 2
        fxm = dcl*Ai_cl(1,kd) - Ai(2,kd) + Ai(1,kd)
        fxo = dcl*Ai_cl(2,kd) - Ai(2,kd) + Ai(1,kd)
        Aij(kd) = tcl *Ai(2,kd)
	+ (1.0-tcl)*Ai(1,kd) 
	+ tcl*(1.0-tcl)*((1.0-tcl)*fxm - tcl*fxo)
    
        fxm = dcl*Ai_cl_tau(1,kd) - Ai_tau(2,kd) + Ai_tau(1,kd)
        fxo = dcl*Ai_cl_tau(2,kd) - Ai_tau(2,kd) + Ai_tau(1,kd)
        Aij_tau(kd) = tcl *Ai_tau(2,kd)
	+ (1.0-tcl)*Ai_tau(1,kd) 
	+ tcl*(1.0-tcl)*((1.0-tcl)*fxm - tcl*fxo)
      end

      fxm = dtau*Aij_tau(1) - Aij(2) + Aij(1)
      fxo = dtau*Aij_tau(2) - Aij(2) + Aij(1)
      Aijk(l) =    ttau *Aij(2)
	+ (1.0-ttau)*Aij(1) 
	+ ttau*(1.0-ttau)*((1.0-ttau)*fxm - ttau*fxo)

      end

      cdf = Aijk(1)
      cdp = Aijk(2)
      cm  = Aijk(3)
      cdw = 0.
#     cdw = Aijk(4)

#---- add quadratic cd penalty for exceeding database's cl limits
      clmin = Acl(1)
      clmax = Acl(nAcl)
      if    (cl < clmin) 
       cdp = cdp + 1.0*(cl-clmin)^2
      elseif(cl > clmax) 
       cdp = cdp + 1.0*(cl-clmax)^2
      end

#---- add quadratic cd penalty for exceeding database's tau limits
      taumin = Atau(1)
      taumax = Atau(nAtau)
      if    (tau < taumin) 
       cdp = cdp + 25.0*(tau-taumin)^2
      elseif(tau > taumax) 
       cdp = cdp + 25.0*(tau-taumax)^2
      end

      return
      end # airfun

