      Program eqdsk_adj_sym_new
c     gfortran eqdsk_adj.f -o scaleIp
c     Bob Harvey,  Sept. '94, Initial work.

c-----------------------------------------------------------------------
c     Bob Harvey, March 1995, with several additions thereafter.
c     This program adjusts an eqdsk by redimensioning it from
c     nr x nz R,Z-points to nra x nza points (typically 33x65) 
c     It possible to crop out most of 
c     the portion of the plasma outside the last flux surface.  
c     (Leave noutside=3 points, as specified in data statement below).
c     (Choose noutside.ge.size_of_eqdsk, if want to ensure no
c      points are cropped.)
c
c     Code parameters below are to be adjusted with a text editor,
c     then recompile the code.
c
c     Parameters nr and nz, below, must be adjusted in accord with the
c     input eqdsk.
c
c     Toroidal and poloidal bfield can be scaled by resetting 
c     scalebt and scalebp (below).
c
c     Also, if isym (below).eq.1,
c     then eqdsk_adj symmetrizes the equilibrium
c     in the vertical direction about the height of the magnetic axis. 
c     The output eqdsk has the magnetic axis on the midplane, 
c     i.e., at zax=0.0
c
c     If isym.ne.1, then there are two possibilities, to vertically 
c     shift or not to vertically shift the magnetic axis:
c     ishift.eq.1: coords are shifted such that the
c     magnetic axis is at z=0, 
c     OR,
c     (ishift.ne.1), there is no shifting in the vertical dirn.
c     [This can be used to simply re-grid the equilibrium to 
c      a differently dimensioned eqdsk. BobH, 030624]
c
c     Or, this code can be used to simply reverse the sign of
c     of the psi quantities (except Ip), and write out the
c     eqdsk, if isgnonly=1.  (default: isgnonly=0). (Bob H, 980820).
c
c     BH020901: 
c     Added ability to read nonstandard eqdskin with nveqd(.ne.nxeqd)
c     length array of f,ff',p,p' versus psi.  This is activated if
c     nv_reset.eq.1.  The output eqdskout is in
c     standard form with nv=nra points in these arrays.
c     NOTE: Some eqdsk writers add text or numbers at the posn 
c           in the eqdsk referenced below by nveqd.  This will
c           cause the present code to fail, for nv_reset.eq.1.
c           If using nv_reset.eq.1, you could check the eqdsk 
c           to see there is no such confusion.
c     BH:  040225,  modified scale to scalebt,scalebp
c
c     To compile and load  
c
c     For Cray J-90:
c     f90 -o eqdsk_adj  -g -O0 -ev eqdsk_adj.f
c
c     For pgf77:
c     pgf77 -o eqdsk_adj  -g -r8 eqdsk_adj.f
c
c     For Linux:
c     gfortran -o eqdsk_adj  -g eqdsk_adj.f
c
c**** ===>>   REMEMBER: Check each of isym, scalebt. scalebp, 
c             and isgnonly, ishift
c             to make sure you are getting what you want.
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c     scalebt/scalebp = scaling factor to produce scaled output eqdsks,
c     which have toroidal/poloidal flux (respectively)adjusted by 
c     factor scalebt and scalebp.
c     If scalebp=scalebt, this keeps the safety factor unchanged, 
c     but toroidal field and
c     current are scaled by factor scalebt=scalebp.
c     If scalebt is negative, then the toroidal magnetic field terms
c     are reversed.
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c

c     Choose nr and nz to equal input eqdskin parameters.
c     Following choices are common:
       parameter (nr=257,nz=257)
c     The following parameters give output eqdsk parameters:
      parameter (nra=257,nza=257)

      dimension peqd(nr,nz),psieqd(nr),ffpeqd(nr),ppeqd(nr),er(nr),
     +  ez(nz),feqd(nr),preseqd(nr),qeqd(nr)
      dimension psieqd1(nr),temp1(nr)
      dimension peqd_ex(nr,nz),er_ex(nr),ez_ex(nz)
      dimension epsinew(nra,nza),psinew(nra),ffpnew(nra),ppnew(nra),
     +  ernew(nra),eznew(nza),fnew(nra),presnew(nra),qnew(nra)
      parameter (maxlimpt = 2000)
      dimension xlimiter(maxlimpt),ylimiter(maxlimpt)
      parameter (nconmax = 2000)
      dimension rcontr(nconmax), zcontr(nconmax)
      character*8 eqdskin,eqdskout,ntitle(5)
c     rakata
      data noutside/267/
      data scalebt/+1.0/,scalebp/+1.667/   
      data isym/0/
      data isgnonly/0/
      data ishift/0/
      data nv_reset/0/
c
      if (isym.eq.1 .and. ishift.ne.1) 
     +             STOP 'isym.eq.1 .and. ishift.ne.1 is not supported'
c
      write(*,*) 'eqdskin  to have nr= ',nr,'    nz=',nz
      write(*,*) 'eqdskout to have nra=',nra,'  nza=',nza
      write(*,*) 'Operation data: noutside,scalebt,scalebp,isym,',
     +           'isgnonly,ishift,nv_reset =',
     +            noutside,scalebt,scalebp,isym,isgnonly,ishift,nv_reset
c
      eqdskin='eqdskin'
      eqdskout='eqdskout'
c
      ncrt=6
      nin=20
      call myopen(nin,eqdskin,2,iread)
      if(iread.ne.1) then
        write(ncrt,8100)
 8100   format(' *** fatal error could not open eqdskin ****')
        stop 
      endif
c
cBH020901  Adding nveqd (.le. nxeqd) for different number of
cBH020901  flux surfaces on which p,f,p',f',q are tabulated.
cBH020901  The standard EQDSK does not incorporate this feature,
cBH020901  although is is available in cql3d.
cBH020901  Bonoli uses it for ACCOME eqdsk output.
cBH020901  Eqdskout will have nra points in these arrays.
      nveqd=0
      if (nv_reset.eq.1) then
c        *****Following line contains nveqd input****.
         read(nin,8190)(ntitle(i),i=1,5),dat,ipestg,nxeqd,nyeqd,nveqd
         if (nveqd.gt.nxeqd) stop 'nveqd.gt.nxeqd NOT ENABLED'
      else
         read(nin,8190)(ntitle(i),i=1,5),dat,ipestg,nxeqd,nyeqd
      endif
      if (nveqd.eq.0) nveqd=nxeqd

      if (nr.ne.nxeqd.or.nz.ne.nyeqd) stop 'check nr,nz'
      read(nin,8200) xdimeqd,ydimeqd,reqd,redeqd,ymideqd
      read(nin,8200) xma,zma,psimag,psilim,beqd
      read(nin,8200) toteqd,psimx1,psimx2,xax1,xax2
      read(nin,8200) zax1,zax2,psisep,xsep,zsep
      read(nin,8200) (feqd(i),i=1,nveqd)
      read(nin,8200) (preseqd(i),i=1,nveqd)
      read(nin,8200) (ffpeqd(i) ,i=1,nveqd)
      read(nin,8200) (ppeqd(i) ,i=1,nveqd)
      read(nin,8200) ((peqd(i,j),i=1,nxeqd),j=1,nyeqd)
      read(nin,8200) (qeqd(i),i=1,nveqd)
      read (nin, 8210, err=8221, end=8221)   ncontr, nlimiter
      write(*,*)'ncontr, nconmax ',ncontr, nlimiter
      if (ncontr.gt.0) then
        read (nin, 8200, err=8221, end=8221)  
     .        (rcontr(i), zcontr(i),i = 1,ncontr)
      endif
      if (nlimiter.gt.0) then
        read (nin, 8200, err=8221, end=8221)
     .       (xlimiter(i), ylimiter(i), i = 1,nlimiter)
      endif
 8190 format(6a8,4i4)
 8200 format(5e16.9)
 8210 format (2i5)

 8221 continue

c-----------------------------------------------------------------------
c     close guess file
c-----------------------------------------------------------------------
      close(nin)

c-----------------------------------------------------------------------
c     If (nveqd.ne.nxeqd),
c     interpolate single dimensioned arrays to nr length.
c-----------------------------------------------------------------------

      if (nveqd.ne.nxeqd) then

         do i=1,nveqd
            psieqd1(i)=psimag-(i-1)*(psimag-psilim)/(nveqd-1)
         enddo
         do i=1,nr
            psieqd(i)=psimag-(i-1)*(psimag-psilim)/(nr-1)
         enddo

         call intrp('free','free',psieqd1,feqd,nveqd,psieqd,temp1,nr)
         do i=1,nr
            feqd(i)=temp1(i)
         enddo

         call intrp('free','free',psieqd1,preseqd,nveqd,
     +                            psieqd,temp1,nr)
         do i=1,nr
            preseqd(i)=temp1(i)
         enddo

         call intrp('free','free',psieqd1,ffpeqd,nveqd,psieqd,temp1,nr)
         do i=1,nr
            ffpeqd(i)=temp1(i)
         enddo

         call intrp('free','free',psieqd1,ppeqd,nveqd,psieqd,temp1,nr)
         do i=1,nr
            ppeqd(i)=temp1(i)
         enddo

         call intrp('free','free',psieqd1,qeqd,nveqd,psieqd,temp1,nr)
         do i=1,nr
            qeqd(i)=temp1(i)
         enddo

         
      endif


c-----------------------------------------------------------------------
c     check for negative beqd. make appropriate changes if necessary.
c-----------------------------------------------------------------------
c     if(beqd.lt.0.) beqd=-beqd
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c     check/change sign of psi.
c     Write out results to eqdskout and exit, if isgnonly=1
c-----------------------------------------------------------------------
      if (psimag.gt.psilim) then
        write(*,*)
        write(*,*)'***********************************************'
        write(*,*)'init:  psimag.gt.psilim, so reversing psi signs'
        write(*,*)'***********************************************'
        write(*,*)
        psimag=-psimag
        psilim=-psilim
        psimx1=-psimx1
        psimx2=-psimx2
        psisep=-psisep
        do 150 i=1,nxeqd
          ffpeqd(i)=-ffpeqd(i)
          ppeqd(i)=-ppeqd(i)
          do 151 j=1,nyeqd
            peqd(i,j)=-peqd(i,j)
 151      continue
 150    continue
      endif
      if (isgnonly.eq.1) then

         nout=30
         call myopen(nout,eqdskout,2,iread)
         if(iread.ne.1) then
            write(ncrt,8101)
 8101       format(' *** fatal error could not open output file ****')
            stop 'Problem with output'
         endif
c
         write(nout,8190)(ntitle(i),i=1,5),dat,ipestg,nxeqd,nyeqd
         if (nr.ne.nxeqd.or.nz.ne.nyeqd) stop 'check nr,nz'
         write(nout,8200) xdimeqd,ydimeqd,reqd,redeqd,ymideqd
         write(nout,8200) xma,zma,psimag,psilim,beqd
         write(nout,8200) toteqd,psimx1,psimx2,xax1,xax2
         write(nout,8200) zax1,zax2,psisep,xsep,zsep
         write(nout,8200) (feqd(i),i=1,nxeqd)
         write(nout,8200) (preseqd(i),i=1,nxeqd)
         write(nout,8200) (ffpeqd(i) ,i=1,nxeqd)
         write(nout,8200) (ppeqd(i) ,i=1,nxeqd)
         write(nout,8200) ((peqd(i,j),i=1,nxeqd),j=1,nyeqd)
         write(nout,8200) (qeqd(i),i=1,nxeqd)
         write (nout, 8210)   ncontr, nlimiter
         if (ncontr .gt. nconmax)  go to 5000
         if (ncontr.gt.0) then
            write (nout, 8200)
     .           (rcontr(i), zcontr(i),i = 1,ncontr)
         endif
         if (nlimiter .gt. maxlimpt-2)  go to 5000
         if (nlimiter.gt.0) then
            write (nout, 8200)
     .           (xlimiter(i), ylimiter(i), i = 1,nlimiter)
         endif

         close(nout)

         go to 999
            
      endif


c-----------------------------------------------------------------------
c     check for consistency in sign between p' and j
c-----------------------------------------------------------------------
c     if(ppeqd(1)*toteqd.lt.0.) then
c     do 2150 i=1,nxeqd
c     ppeqd(i)=-ppeqd(i)
c     ffpeqd(i)=-ffpeqd(i)
c     2150 continue
c     endif
c-----------------------------------------------------------------------
c     use convention, totcur.ge.0
c-----------------------------------------------------------------------
c     if(toteqd.lt.0.) toteqd=-toteqd


c-----------------------------------------------------------------------
c     create er and ez arrays, as given by the eqdsk
c-----------------------------------------------------------------------
      drr=xdimeqd/(nxeqd-1)
      do 2200 ix=1,nxeqd
        er(ix)=(ix-1)*drr+redeqd
 2200 continue

c
      dzz=ydimeqd/(nyeqd-1)
      do 2220 iy=1,nyeqd
        ez(iy)=(iy-1)*dzz-0.5*ydimeqd+ymideqd
 2220 continue

      write(*,*)'ydimeqd,ymideqd,ez ',ydimeqd,ymideqd,ez

c-----------------------------------------------------------------------
c     Prepatory to symmetrizing step:
c     ishift.eq.1:
c     Expand vertical height of eqdsk grid by amount 2*(zma-ymideqd),
c     and shift making midplane of computational grid
c     at the magnetic axis.  This is both an expansion and a shift.
c     Interpolate psi onto new (expanded "_ex") grid.
c     ishift.ne.1:  Go through this section, but no vertical shift.
c-----------------------------------------------------------------------

      if (ishift.eq.1) then
         shift=zma-ymideqd
      else
         shift=0.0
      endif
      zdim_ex=ydimeqd+2.*abs(shift)
      dr_ex=drr
      dz_ex=zdim_ex/(nz-1)
      if (shift.ge.0.) then
        ez_ex(1)=ez(1)
      else
        ez_ex(1)=ez(1)-2.*abs(shift)
      endif
      do 310 ix=1,nr
 310  er_ex(ix)=er(ix)
      do 320 iz=2,nz
 320  ez_ex(iz)=ez_ex(iz-1)+dz_ex

c     Interpolating:
        kz1=1
        kz2=2
        do 350  j=1,nz
          zval=ez_ex(j)
 360      if (zval.le.ez(kz2).or.kz2.eq.nz) go to 365

          kz1=kz1+1
          kz2=kz1+1
          go to 360
c
 365      continue
          kr1=1
          kr2=2
          do 370  i=1,nr
            rval=er_ex(i)
 380        if (rval.le.er(kr2).or.kr2.eq.nr) go to 385

            kr1=kr1+1
            kr2=kr1+1
            go to 380
c
 385        continue
            f1=peqd(kr1,kz1)+(rval-er(kr1))*
     1        (peqd(kr2,kz1)-peqd(kr1,kz1))/(er(kr2)-er(kr1))
c
            f2=peqd(kr1,kz2)+(rval-er(kr1))*
     1        (peqd(kr2,kz2)-peqd(kr1,kz2))/(er(kr2)-er(kr1))
c
            val=f1+(zval-ez(kz1))*(f2-f1)/(ez(kz2)-ez(kz1))
            peqd_ex(i,j)=val
 370      continue
 350    continue

        xma_ex=xma
        rdim_ex=xdimeqd
        if (ishift.eq.1) then
           zmid_ex=0.0
           zma_ex=0.0
        else
           zmid_ex=ymideqd
           zma_ex=zma
        endif

c     Redefine the ez_ex-grid:
        ez_ex(1)=-zdim_ex/2.0 + zmid_ex
        do 390 i=2,nz
 390    ez_ex(i)=ez_ex(i-1)+dz_ex


      write(*,*)'zdim_ex,zmid_ex,ez_ex ',zdim_ex,zmid_ex,ez_ex




c......................................................................
c     Re-represent psi on nra x nza grid using bi-linear interpolation,
c     encompassing plasma plus noutside (e.g., 3) grid points.
c     Interpolate from expanded grid.
c......................................................................

c
c     Find region of grid containing the plasma by stepping
c     along rays issuing from the plasma magnetic axis, and
c     determining if (psi-psilim) has changed sign.


        rleft=er_ex(1)
        rright=er_ex(nr)
        zbot=ez_ex(1)
        ztop=ez_ex(nz)
        raxis=xma_ex
        zaxis=zma_ex

        dss=0.9*amin1(dr_ex,dz_ex)
        dthetap=8.*atan2(1.,1.)/201.
        thetap=0.0
        ismax=max((rright-rleft)/dss,(ztop-zbot)/dss)
        jmin=nz
        jmax=1
        iminn=nr
        imaxx=1

        do 200  it=1,200
          thetap=thetap+dthetap
          coss=cos(thetap)
          sinn=sin(thetap)
          s=0.0

          do 210  is=1,ismax
            s=s+dss
            rval=raxis+s*coss
            zval=zaxis+s*sinn

            if (rval.gt.rright.or.rval.lt.rleft) go to 200
            if (zval.gt.ztop.or.zval.lt.zbot) go to 200

            kr1=(rval-rleft)/dr_ex+1
            kz1=(zval-zbot)/dz_ex+1
            kr1=min0(kr1,nr-1)
            kz1=min0(kz1,nz-1)
            kr2=kr1+1
            kz2=kz1+1

            jmin=min0(jmin,kz1)
            jmax=max0(jmax,kz2)
            iminn=min0(iminn,kr1)
            imaxx=max0(imaxx,kr2)

            f1=peqd_ex(kr1,kz1)+(rval-er_ex(kr1))*
     1        (peqd_ex(kr2,kz1)-peqd_ex(kr1,kz1))
     +        /(er_ex(kr2)-er_ex(kr1))
c
            f2=peqd_ex(kr1,kz2)+(rval-er_ex(kr1))*
     1        (peqd_ex(kr2,kz2)-peqd_ex(kr1,kz2))
     +        /(er_ex(kr2)-er_ex(kr1))
c
            val=f1+(zval-ez_ex(kz1))*(f2-f1)/(ez_ex(kz2)-ez_ex(kz1))

            if ((val-psilim)*(psimag-psilim).lt.0.0) go to 200

 210      continue
 200    continue

      jmax=min0(jmax+noutside,nz)
      jmin=max0(jmin-noutside,1)
      jl=jmin-1
      ju=nz-jmax
      if (ju.le.jl) then
        jmin=nz-jmax+1
      else
        jmax=nz-jmin+1
      endif

        imaxx=min0(imaxx+noutside,nr)
        iminn=max0(1,iminn-noutside)

      if ((jmin.eq.1.and.jmax.eq.nz).and.
     1   (iminn.eq.1.and.imaxx.eq.nr)) then
          write(ncrt,8102)
 8102     format(' *** No change in eqdskin R and Z ranges ***')
        endif

        drr=(er_ex(imaxx)-er_ex(iminn))/(nra-1)
        dzz=(ez_ex(jmax)-ez_ex(jmin))/(nza-1)
        do 232  i=1,nra
 232    ernew(i)=er_ex(iminn)+(i-1)*drr
        do 233  j=1,nza
 233    eznew(j)=ez_ex(jmin)+(j-1)*dzz

        kz1=1
        kz2=2
        do 250  j=1,nza
          zval=eznew(j)
 260      if (zval.le.ez_ex(kz2)) go to 265

          kz1=kz1+1
          kz2=kz1+1
          go to 260
c
 265      continue
          kr1=1
          kr2=2
          do 270  i=1,nra
            rval=ernew(i)
 280        if (rval.le.er_ex(kr2)) go to 285

            kr1=kr1+1
            kr2=kr1+1
            go to 280
c
 285        continue
            f1=peqd_ex(kr1,kz1)+(rval-er_ex(kr1))*
     1        (peqd_ex(kr2,kz1)-peqd_ex(kr1,kz1))
     +        /(er_ex(kr2)-er_ex(kr1))
c
            f2=peqd_ex(kr1,kz2)+(rval-er_ex(kr1))*
     1        (peqd_ex(kr2,kz2)-peqd_ex(kr1,kz2))
     +        /(er_ex(kr2)-er_ex(kr1))
c
            val=f1+(zval-ez_ex(kz1))*(f2-f1)/(ez_ex(kz2)-ez_ex(kz1))
            epsinew(i,j)=val
 270      continue
 250    continue

c
        rdimnew=ernew(nra)-ernew(1)
        zdimnew=eznew(nza)-eznew(1)
        rednew=ernew(1)
        zmidnew=zmid_ex
        rmagnew=xma_ex
        zmagnew=zma_ex
        rnew=reqd
        zax1new=zma_ex
        zsepnew=0.0        !0.0, since probably don't need to 
                           !     figure it out.


      write(*,*)'zdimnew,zmidnew,eznew ',zdimnew,zmidnew,eznew


c     Scalebt,scalebp:
        scalebp2=scalebp**2
        do 300 i=1,nr
          psieqd(i)=(psimag-(i-1)*(psimag-psilim)/(nr-1))*abs(scalebp)
          feqd(i)=scalebt*feqd(i)
          preseqd(i)=scalebp2*preseqd(i)
          ffpeqd(i)=scalebt*ffpeqd(i)/scalebp
          ppeqd(i)=abs(scalebp)*ppeqd(i)
          qeqd(i)=qeqd(i)*scalebt/scalebp
 300    continue

        psimag=abs(scalebp)*psimag
        psilim=abs(scalebp)*psilim
        beqd=scalebt*beqd
cBH120713
        toteqd=scalebp*toteqd
        psimx1=abs(scalebp)*psimx1
        psimx2=abs(scalebp)*psimx2
        psisep=abs(scalebp)*psisep

        do 330 i=1,nra
          psinew(i)=psimag-(i-1)*(psimag-psilim)/(nra-1)
 330    continue

        call intrp('free','free',psieqd,feqd,nr,psinew,fnew,nra)
        call intrp('free','free',psieqd,preseqd,nr,psinew,presnew,nra)
        call intrp('free','free',psieqd,ffpeqd,nr,psinew,ffpnew,nra)
        call intrp('free','free',psieqd,ppeqd,nr,psinew,ppnew,nra)
        call intrp('free','free',psieqd,qeqd,nr,psinew,qnew,nra)

        if (scalebt.ne.1.0 .or. scalebp.ne.1.0) ntitle(4)='scaled'

        if (isym.eq.1) then
c..................................................................
c     Up-down symmetrize the eqdsk, about z=0. (this has no effect
c     if the equilibrium is initially up-down symmetric).
c..................................................................

        do 240  j=1,nza/2
          do 241  i=1,nra
            epsinew(i,j)=0.5*(epsinew(i,j)+epsinew(i,nza-(j-1)))
 241      continue
 240    continue
          do 242  j=nza-nza/2+1,nza
            do 243  i=1,nra
              epsinew(i,j)=epsinew(i,nza+1-j)
 243        continue
 242      continue

       endif
       
       do 290 j=1,nza
          do 291 i=1,nra
             epsinew(i,j)=abs(scalebp)*epsinew(i,j)
 291      continue
 290   continue
       
c..................................................................
c
       nout=30
       call myopen(nout,eqdskout,2,iread)
       if(iread.ne.1) then
          write(ncrt,8101)
          stop 'Problem with output'
       endif
c     
       write(nout,8190)(ntitle(i),i=1,5),dat,ipestg,nra,nza
       write(nout,8200) rdimnew,zdimnew,rnew,rednew,zmidnew
       write(nout,8200) rmagnew,zmagnew,psimag,psilim,beqd
       write(nout,8200) toteqd,psimx1,psimx2,xax1,xax2
       write(nout,8200) zax1new,zax2,psisep,xsep,zsepnew
       write(nout,8200) (fnew(i),i=1,nra)
       write(nout,8200) (presnew(i),i=1,nra)
       write(nout,8200) (ffpnew(i) ,i=1,nra)
       write(nout,8200) (ppnew(i) ,i=1,nra)
       write(nout,8200) ((epsinew(i,j),i=1,nra),j=1,nza)
       write(nout,8200) (qnew(i),i=1,nra)
       if (ncontr.gt.0 .or. nlimiter.gt.0) then
          write (nout, 8210) ncontr,nlimiter
       endif
       if (ncontr.gt.0) then
          write (nout,8200)
     .         (rcontr(i), zcontr(i)-shift,i = 1,ncontr)
       endif
       if (nlimiter.gt.0) then
          write(nout,8200)
     .         (xlimiter(i), ylimiter(i)-shift,i = 1,nlimiter)
       endif

c-----------------------------------------------------------------------
c     close output file
c-----------------------------------------------------------------------
       close(nout)
c------------------------------------------------------------------------
       
       go to 999

c
c --- fatal errors, stop code
c
 5000 write  (*, 5100)  nconmax,maxlimpt,ncontr,nlimtr
 5100 format (/ ' ERROR detected '          /
     .          '   parameter settings inconsistent with eqdsk' /
     .          '   nconmax, maxlimpt = ', 2(2x, i6)     /
     .          '   ncontr , nlimtr   = ', 2(2x, i6)     /
     .          '   cannot continue')

  999 continue
            end


            subroutine myopen(nunit,name,itype,iread)
            integer*4 nunit
            logical iopen
            character*8 name
            write(*,*)'myopen: nunit, name, itype',nunit, name, itype
            iread=0
            if (itype.eq.1)then
              open(nunit,file=name,status='unknown',form='unformatted')
            else
              open(nunit,file=name,status='unknown',form='formatted')
            endif
            inquire(file=name,opened=iopen)
            if(iopen) iread=1
            return
            end




            subroutine intrp(ilow,iup,x,y,nold,xnew,ynew,nnew)
c------------------------------------------------------------------
c     This subroutine uses the imsl cubic spline routines to calculate
c     ynew as a function of xnew given y as a function of x.  The
c     boundary conditions are specified as follows:
c     ilow : 'zero', set derivative at lower x value to zero
c     not eq'zero'or 'fixed', let derivative at lower x value be free
c     'fixed' set derivate according to bpar
c     iup  : 'zero', set derivative at upper x value to zero
c     not eq 'zero'or 'fixed', let derivative at upper x value be free
c     'fixed' set derivate according to bpar
c------------------------------------------------------------------
c
c
            dimension xnew(*),ynew(*),x(*),y(*)
            parameter (kjx=404)
            dimension bpar(4),c(kjx,3),xt(kjx),yt(kjx)
            character ilow*(*),iup*(*)
c
            if (nold.gt.kjx) stop 'nold.gt.kjx, in intrp'
c
c     must assure that x(1).le.x(2)....
c
            if(x(2).le.x(1)) go to 10
            iord=0
c
 60         continue
            do 4 i=1,4
              bpar(i)=0.
 4          continue
            if((ilow.eq.'fixed').or.(iup.eq.'fixed'))then
              stop 'intrp not set up for fixed boundary conditions'
            endif
            if(ilow.ne.'zero') go to 6
            bpar(1) = 1.
            bpar(2) = 6.*(y(2)-y(1))/(x(2)-x(1))**2
 6          if(iup.ne.'zero') go to 8
            bpar(3) = 1.
            bpar(4) = -6.*(y(nold)-y(nold-1))/(x(nold)-x(nold-1))**2
 8          continue
c
            call icsicu(x,y,nold,bpar,c,kjx,ier)
c
            if(ier.ne.0) go to 20
c
            call icsevu(x,y,nold,c,kjx,xnew,ynew,nnew,ier)
c
            if(iord.eq.0) return
c
            do 70 i=1,nold
              x(i)=xt(i)
              y(i)=yt(i)
 70         continue
c
            return
c
 10         continue
            iord=1
c
            do 30 i=1,nold
              xt(i)=x(i)
 30         yt(i)=y(i)
c
            do 50 i=1,nold
              x(i)=xt(nold+1-i)
 50         y(i)=yt(nold+1-i)
c
            go to 60
c
c
 20         write(66,40) ier,iord,nold
 40         format('1intrp    ier=',i3,'   iord=',i1,'   nold=',i3,/,
     .        5h  i  ,12h     x      ,12h     y       )
            write(66,80) (i,x(i),y(i),i=1,nold)
 80         format(i5,2e12.4)
c
            stop 'in subroutine INTRP'
c
            end


            subroutine icsicu (x,y,nx,bpar,c,ic,ier)
c     imsl routine name   - icsicu
c
c-----------------------------------------------------------------------
c
c     computer            - dec10/single
c
c     latest revision     - january 1, 1978
c
c     purpose             - interpolatory approximation by cubic splines
c     with arbitrary second derivative end
c     conditions.
c
c     usage               - call icsicu (x,y,nx,bpar,c,ic,ier)
c
c     arguments    x      - vector of length nx containing the abscissae
c     of the nx data points (x(i),y(i)) i=1,...,
c     nx. (input) x must be ordered so that
c     x(i) .lt. x(i+1).
c     y      - vector of length nx containing the ordinates
c     (or function values) of the nx data points.
c     (input)
c     nx     - number of elements in x and y. (input) nx
c     must be .ge. 2.
c     bpar   - vector of length 4 containing the end
c     condition parameters. (input)
c     2.0*spp(1)+bpar(1)*spp(2) = bpar(2),
c     bpar(3)*spp(nx-1)+2.0*spp(nx) = bpar(4),
c     where spp(i) = second derivative of the
c     cubic spline function s evaluated at x(i).
c     c      - spline coefficients. (output) c is an nx-1 by
c     3 matrix. the value of the spline
c     approximation at t is
c     s(t) = ((c(i,3)*d+c(i,2))*d+c(i,1))*d+y(i)
c     where x(i) .le. t .lt. x(i+1) and
c     d = t-x(i).
c     ic     - row dimension of matrix c exactly as
c     specified in the dimension statement in
c     the calling program. (input)
c     ier    - error parameter. (output)
c     terminal error
c     ier = 129, ic is less than nx-1.
c     ier = 130, nx is less than 2.
c     ier = 131, input abscissa are not ordered
c     so that x(1) .lt. x(2) ... .lt. x(nx).
c
c     precision/hardware  - single and double/h32
c     - single/h36,h48,h60
c
c     reqd. imsl routines - uertst,ugetio
c
c     notation            - information on special notation and
c     conventions is available in the manual
c     introduction or through imsl routine uhelp
c
c     copyright           - 1978 by imsl, inc. all rights reserved.
c
c     warranty            - imsl warrants only that imsl testing has been
c     applied to this code. no other warranty,
c     expressed or implied, is applicable.
c
c-----------------------------------------------------------------------
c
c     specifications for arguments
            integer            nx,ic,ier
            real               x(nx),y(nx),bpar(4),c(ic,3)
c     specifications for local variables
            integer            i,j,nxm1
            real               dx,dxj,dxjp1,dxp,dyj,dyjp1,half,one,pj,
     1        six,sixi,two,yppa,yppb,zero
            equivalence        (dxj,yppb),(pj,sixi),(dxjp1,yppa)
            data               zero/0.0/,half/0.5/,one/1.0/,
     1        two/2.0/,six/6.0/
c     first executable statement
            ier = 0
c     check error conditions
            nxm1 = nx-1
            if (ic .lt. nxm1) go to 30
            if (nx .lt. 2) go to 35
            if (nx .eq. 2) go to 10
c     compute coefficients and right
c     hand side of the tridiagonal
c     system defining the second
c     derivatives of the spline
c     interpolant for (x,y)
c     c(j,1) = lambda(j)
c     c(j,2) = mu(j)
c     c(j,3) = d(j)
            dxj = x(2)-x(1)
            if (dxj .le. zero) go to 40
            dyj = y(2)-y(1)
            do 5 j=2,nxm1
              dxjp1 = x(j+1)-x(j)
              if (dxjp1 .le. zero) go to 40
              dyjp1 = y(j+1)-y(j)
              dxp = dxj+dxjp1
              c(j,1) = dxjp1/dxp
              c(j,2) = one-c(j,1)
              c(j,3) = six*(dyjp1/dxjp1-dyj/dxj)/dxp
              dxj = dxjp1
              dyj = dyjp1
 5          continue
c     factor the tridiagonal matrix
c     and solve for u
c     c(j,2) = u(j)
c     c(j,1) = q(j)
c     bpar(1) = lambda(1)
c     bpar(2) = d(1)
c     bpar(3) = mu(nx)
c     bpar(4) = d(nx)
 10         c(1,1) = -bpar(1)*half
            c(1,2) = bpar(2)*half
            if (nx .eq. 2) go to 20
            do 15 j=2,nxm1
              pj = c(j,2)*c(j-1,1)+two
              c(j,1) = -c(j,1)/pj
              c(j,2) = (c(j,3)-c(j,2)*c(j-1,2))/pj
 15         continue
c     solve for cubic coefficients
c     of spline interpolant
c     c(j,1), c(j,2), and c(j,3)
 20         yppb = (bpar(4)-bpar(3)*c(nxm1,2))/(bpar(3)*c(nxm1,1)+two)
            sixi = one/six
            do 25 i=1,nxm1
              j = nx-i
              yppa = c(j,1)*yppb+c(j,2)
              dx = x(j+1)-x(j)
              c(j,3) = sixi*(yppb-yppa)/dx
              c(j,2) = half*yppa
              c(j,1) = (y(j+1)-y(j))/dx-(c(j,2)+c(j,3)*dx)*dx
              yppb = yppa
 25         continue
            go to 9005
 30         ier = 129
            go to 9000
 35         ier = 130
            go to 9000
 40         ier = 131
 9000       continue
cBH220415           call uertst(ier,6hicsicu)
 9005       return
            end


            subroutine icsevu  (x,y,nx,c,ic,u,s,m,ier)
c     imsl routine name   - icsevu
c
c-----------------------------------------------------------------------
c
c     computer            - dec10/single
c
c     latest revision     - january 1, 1978
c
c     purpose             - evaluation of a cubic spline
c
c     usage               - call icsevu(x,y,nx,c,ic,u,s,m,ier)
c
c     arguments    x      - vector of length nx containing the abscissae
c     of the nx data points (x(i),y(i)) i=1,...,
c     nx (input). x must be ordered so that
c     x(i) .lt. x(i+1).
c     y      - vector of length nx containing the ordinates
c     (or function values) of the nx data points
c     (input).
c     nx     - number of elements in x and y (input).
c     nx must be .ge. 2.
c     c      - spline coefficients (input). c is an nx-1 by
c     3 matrix.
c     ic     - row dimension of matrix c exactly as
c     specified in the dimension statement
c     in the calling program (input).
c     ic must be .ge. nx-1.
c     u      - vector of length m containing the abscissae
c     of the m points at which the cubic spline
c     is to be evaluated (input).
c     s      - vector of length m (output).
c     the value of the spline approximation at
c     u(i) is
c     s(i) = ((c(j,3)*d+c(j,2))*d+c(j,1))*d+y(j)
c     where x(j) .le. u(i) .lt. x(j+1) and
c     d = u(i)-x(j).
c     m      - number of elements in u and s (input).
c     ier    - error parameter (output).
c     warning error
c     ier = 33, u(i) is less than x(1).
c     ier = 34, u(i) is greater than x(nx).
c*****************************************************
c     output of warning errors has been suppressed
c*****************************************************
c
c     precision/hardware  - single and double/h32
c     - single/h36,h48,h60
c
c     reqd. imsl routines - uertst,ugetio
c
c     notation            - information on special notation and
c     conventions is available in the manual
c     introduction or through imsl routine uhelp
c
c     remarks  1.  the routine assumes that the abscissae of the nx
c     data points are ordered such that x(i) is less than
c     x(i+1) for i=1,...,nx-1. no check of this condition
c     is made in the routine. unordered abscissae will cause
c     the algorithm to produce incorrect results.
c     2.  the routine generates two warning errors. one error
c     occurs if u(i) is less than x(1), for some i in the
c     the interval (1,m) inclusively. the other error occurs
c     if u(i) is greater than x(nx), for some i in the
c     interval (1,m) inclusively.
c     3.  the ordinate y(nx) is not used by the routine. for
c     u(k) .gt. x(nx-1), the value of the spline, s(k), is
c     given by
c     s(k)=((c(nx-1,3)*d+c(nx-1,2))*d+c(nx-1,1))*d+y(nx-1)
c     where d=u(k)-x(nx-1).
c
c     copyright           - 1978 by imsl, inc. all rights reserved.
c
c     warranty            - imsl warrants only that imsl testing has been
c     applied to this code. no other warranty,
c     expressed or implied, is applicable.
c
c-----------------------------------------------------------------------
c
c     specifications for arguments
            integer            nx,ic,m,ier
            real               x(nx),y(nx),c(ic,3),u(m),s(m)
c     specifications for local variables
            integer            i,jer,ker,nxm1,k
            real               d,dd,zero
            data               i/1/,zero/0.0/
c     first executable statement
            jer = 0
            ker = 0
            if (m .le. 0) go to 9005
            nxm1 = nx-1
            if (i .gt. nxm1) i = 1
c     evaluate spline at m points
            do 40 k=1,m
c     find the proper interval
              d = u(k)-x(i)
              if (d) 5,25,15
 5            if (i .eq. 1) go to 30
              i = i-1
              d = u(k)-x(i)
              if (d) 5,25,20
 10           i = i+1
              d = dd
 15           if (i .ge. nx) go to 35
              dd = u(k)-x(i+1)
              if (dd .ge. zero) go to 10
              if (d .eq. zero) go to 25
c     perform evaluation
 20           s(k) = ((c(i,3)*d+c(i,2))*d+c(i,1))*d+y(i)
              go to 40
 25           s(k) = y(i)
              go to 40
c     u(k) .lt. x(1)
 30           jer = 33
              go to 20
c     u(k) .gt. x(nx)
 35           if (dd .gt. zero) ker = 34
              d = u(k)-x(nxm1)
              i = nxm1
              go to 20
 40         continue
            ier = max0(jer,ker)
c     if (jer .gt. 0) call uertst(jer,6hicsevu)
c     if (ker .gt. 0) call uertst(ker,6hicsevu)
 9005       return
            end



            subroutine uertst(ier,name)
c
c-uertst----------------library2--------------------------------------
c
c     function            - error message generation
c     usage               - call uertst(ier,name)
c     parameters   ier    - error parameter. type + n  where
c     type= 128 implies terminal error
c     64 implies warning with fix
c     32 implies warning
c     n   = error code relevant to calling routin
c     name   - input scalar (double precision on dec)
c     containing the name of the calling routine
c     as a 6-character literal string.
c     language            - fortran
c----------------------------------------------------------------------
c     latest revision     - october 1,1975
c     dec
c
            character*5 ityp
            dimension          ityp(4,4),ibit(4)
            double precision   name
            integer            warn,warf,term,printr
            equivalence        (ibit(1),warn),(ibit(2),warf),
     +        (ibit(3),term)
            data               ityp/'warni','ng   ','     ','     ',
     *        'warni','ng(wi','th fi','x)   ',
     *        'termi','nal  ','     ','     ',
     *        'non-d','efine','d    ','     '/,
     *        ibit/ 32,64,128,0/
            data               printr/ 6/
            ier2=ier
            if (ier2 .ge. warn) go to 5
c     non-defined
            ier1=4
            go to 20
 5          if (ier2 .lt. term) go to 10
c     terminal
            ier1=3
            go to 20
 10         if (ier2 .lt. warf) go to 15
c     warning(with fix)
            ier1=2
            go to 20
c     warning
 15         ier1=1
c     extract 'n'
 20         ier2=ier2-ibit(ier1)
c     print error message
            write (6,25) imslmd,(ityp(i,ier1),i=1,4),name,ier2,ier
 25         format(' from ',a8,': ',
     1        ' *** i m s l(uertst) ***  ',4a5,2x,a6,2x,i2,
     1        ' (ier = ',i3,')')
            return
            end


