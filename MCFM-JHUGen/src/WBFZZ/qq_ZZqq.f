      subroutine qq_ZZqq(p,msq)
      implicit none
c--- Author: R.K. Ellis, October 2014
c--- q(-p1)+q(-p2)->Z(p3,p4)+Z(p5,p6)+q(p7)+q(p8);
      include 'constants.f'
      include 'cmplxmass.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'runstring.f'
      include 'zprods_decl.f'
      include 'anom_higgs.f'
!      include 'first.f'
      include 'spinzerohiggs_anomcoupl.f'
      include 'interference.f'
      include 'pid_pdg.f'
      include 'WWbits.f'
      integer nmax,jmax
      parameter(jmax=12,nmax=10)
      integer j,k,l,
     & uqcq_uqcq,uquq_uquq,dqsq_dqsq,
     & dqdq_dqdq,uqbq_uqbq,dqcq_dqcq,
     & dquq_dquq,dqcq_uqsq,uqsq_dqcq
      parameter(
     & uqcq_uqcq=1,uquq_uquq=2,dqsq_dqsq=3,
     & dqdq_dqdq=4,uqbq_uqbq=5,dqcq_dqcq=6,
     & dquq_dquq=7,dqcq_uqsq=8,uqsq_dqcq=9)
      integer h1,h2,h3,h5
      double precision p(mxpart,4),msq(fn:nf,fn:nf),temp(fn:nf,fn:nf),
     & tempw(fn:nf,fn:nf),stat,spinavge,mult,
     & colfac34_56,ampsqfac
      double complex zab(mxpart,4,mxpart),zba(mxpart,4,mxpart),
     & amp(nmax,2,2,2,2),ampa(nmax,2,2,2,2),ampb(nmax,2,2,2,2),
     & amp_swap(nmax,2,2,2,2),
     & ampa_swap(nmax,2,2,2,2),ampb_swap(nmax,2,2,2,2)
      logical doHO,doBO,comb1278ok
      logical isALepton,isANeutrino
      parameter(spinavge=0.25d0,stat=0.5d0)
      integer,parameter:: j1(jmax)=(/1,2,8,8,7,2,7,1,1,7,2,7/)
      integer,parameter:: j2(jmax)=(/2,1,7,7,2,7,1,7,7,1,7,2/)
      integer,parameter:: j7(jmax)=(/7,7,2,1,1,8,2,8,2,8,1,8/)
      integer,parameter:: j8(jmax)=(/8,8,1,2,8,1,8,2,8,2,8,1/)
      save doHO,doBO,mult
!$omp threadprivate(doHO,doBO,mult)

      msq(:,:)=0d0
      ampsqfac = esq**6*spinavge

c--- This calculation uses the complex-mass scheme (c.f. arXiv:hep-ph/0605312)
c--- and the following lines set up the appropriate masses and sin^2(theta_w)
      cwmass2=dcmplx(wmass**2,0d0)
      czmass2=dcmplx(zmass**2,0d0)
      cxw=dcmplx(xw,0d0)

      doHO=.false.
      doBO=.false.
      if     (runstring(4:5) .eq. 'HO') then
        doHO=.true.
      elseif (runstring(4:5) .eq. 'BO') then
        doBO=.true.
      endif
      if (doHO) then
        Hbit=cone
        Bbit=czip
      elseif (doBO) then
        Hbit=czip
        Bbit=cone
      else
        Hbit=cone
        Bbit=cone
      endif

c--- rescaling factor for Higgs amplitudes, if anomalous Higgs width
       mult=1d0
       if (anom_Higgs) then
         mult=chi_higgs**2
       endif
       Hbit=mult*Hbit

C---call plabel/pdgid conversion
      call convertPLabelsToPDGIds()

C---setup spinors and spinorvector products
      call spinorcurr(8,p,za,zb,zab,zba)

c---setup Z/A couplings from PDG ids
      call couplzajk()

c---color factors for Z decays
      colfac34_56=1d0
      if (abs(pid_pdg(3)).ge.0 .and. abs(pid_pdg(3)).le.5) then
        colfac34_56=colfac34_56*xn
      endif
      if (abs(pid_pdg(5)).ge.0 .and. abs(pid_pdg(5)).le.5) then
        colfac34_56=colfac34_56*xn
      endif

      do j=1,jmax
      temp(:,:)=0d0
      tempw(:,:)=0d0
      amp(:,:,:,:,:)=czip
      ampa(:,:,:,:,:)=czip
      ampb(:,:,:,:,:)=czip
      amp_swap(:,:,:,:,:)=czip
      ampa_swap(:,:,:,:,:)=czip
      ampb_swap(:,:,:,:,:)=czip

C--   MARKUS: adding switches to remove VH or VBF contributions
      if( (vvhvvtoggle_vbfvh.eq.1) .and. (j.le.4) ) cycle
      ! U. Sarica: Test the combination
      call testWBFVVApartComb(j1(j),j2(j),j7(j),j8(j),comb1278ok)
      if (.not.comb1278ok) cycle
      if( (
     &    (isALepton(abs(pid_pdg(7))) .or. isANeutrino(abs(pid_pdg(7))))
     &    .or.
     &    (isALepton(abs(pid_pdg(8))) .or. isANeutrino(abs(pid_pdg(8))))
     &    ) .and. j.lt.9
     & ) then
         cycle
      endif

! U. SARICA TEST
!      write(6,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!      write(6,*) "Begin j=",j
!      write(6,*) " -> j1, j2, j7, j8=",j1(j),j2(j),j7(j),j8(j)
!      write(6,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

c--   Call the VVZZ amplitudes
      call getVVZZamps(amp,ampa,ampb,za,zb,zab,zba,
     & j1(j),j2(j),3,4,5,6,j7(j),j8(j),doHO,doBO)
      if (interference) then
        call getVVZZamps(amp_swap,ampa_swap,ampb_swap,za,zb,zab,zba,
     &   j1(j),j2(j),3,6,5,4,j7(j),j8(j),doHO,doBO)
      endif

      ! Kill amp or ampa in j=5,6,7,8 for VH
      ! Kill amp or ampa in j=9,10,11,12 for VBF
      if( (
     & (vvhvvtoggle_vbfvh.eq.1) .and. (j.le.8)
     & ) .or. (
     & (vvhvvtoggle_vbfvh.eq.0) .and. (j.ge.9)
     & )
     &  ) then
         amp(:,:,:,:,:)=czip
         amp_swap(:,:,:,:,:)=czip
         ampa(:,:,:,:,:)=czip
         ampa_swap(:,:,:,:,:)=czip
      ! Kill ampb in j=9,10,11,12 for VH
      ! Kill ampb in j=5,6,7,8 for VBF
      else if( (
     & (vvhvvtoggle_vbfvh.eq.1) .and. (j.ge.9)
     & ) .or. (
     & (vvhvvtoggle_vbfvh.eq.0) .and. (j.le.8 .and. j.ge.5)
     & )
     &  ) then
         ampb(:,:,:,:,:)=czip
         ampb_swap(:,:,:,:,:)=czip
      endif
      if(
     &    (isALepton(abs(pid_pdg(7))) .or. isANeutrino(abs(pid_pdg(7))))
     &    .or.
     &    (isALepton(abs(pid_pdg(8))) .or. isANeutrino(abs(pid_pdg(8))))
     & ) then
         ampa(:,:,:,:,:)=czip
         ampa_swap(:,:,:,:,:)=czip
         ampb(:,:,:,:,:)=czip
         ampb_swap(:,:,:,:,:)=czip
      endif


C-----setup for (uqbq_uqbq) (2,5)->(2,5)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2

      temp(2,5)=temp(2,5)+ampsqfac
     &   *dble(amp(uqbq_uqbq,h1,h2,h3,h5)
     & *dconjg(amp(uqbq_uqbq,h1,h2,h3,h5)))

      if (interference) then
      temp(2,5)=temp(2,5)+ampsqfac
     &   *dble(amp_swap(uqbq_uqbq,h1,h2,h3,h5)
     & *dconjg(amp_swap(uqbq_uqbq,h1,h2,h3,h5)))
         if(h3 .eq. h5) then
      temp(2,5)=temp(2,5)-2d0*ampsqfac/sqrt(colfac34_56)
     &   *dble(amp(uqbq_uqbq,h1,h2,h3,h5)
     & *dconjg(amp_swap(uqbq_uqbq,h1,h2,h3,h5)))
         endif
      endif

      enddo
      enddo
      enddo
      enddo
      temp(4,5)=temp(2,5)

C--------------------------------------------------------------------------
C-----setup for (uqcq_uqcq) (2,4)->(2,4)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2

      temp(2,4)=temp(2,4)+ampsqfac
     &   *dble(amp(uqcq_uqcq,h1,h2,h3,h5)
     & *dconjg(amp(uqcq_uqcq,h1,h2,h3,h5)))

      if (interference) then
      temp(2,4)=temp(2,4)+ampsqfac
     &   *dble(amp_swap(uqcq_uqcq,h1,h2,h3,h5)
     & *dconjg(amp_swap(uqcq_uqcq,h1,h2,h3,h5)))
         if(h3 .eq. h5) then
      temp(2,4)=temp(2,4)-2d0*ampsqfac/sqrt(colfac34_56)
     &   *dble(amp(uqcq_uqcq,h1,h2,h3,h5)
     & *dconjg(amp_swap(uqcq_uqcq,h1,h2,h3,h5)))
         endif
      endif

      enddo
      enddo
      enddo
      enddo


C-----setup for uqsq_dqcq W diagrams (2,3)->(1,4)
      do h1=1,1
      do h2=1,1
      do h3=1,2
      do h5=1,2

      tempw(2,3)=tempw(2,3)+ampsqfac
     &   *dble(amp(uqsq_dqcq,h1,h2,h3,h5)
     & *dconjg(amp(uqsq_dqcq,h1,h2,h3,h5)))

      if (interference) then
      tempw(2,3)=tempw(2,3)+ampsqfac
     &   *dble(amp_swap(uqsq_dqcq,h1,h2,h3,h5)
     & *dconjg(amp_swap(uqsq_dqcq,h1,h2,h3,h5)))
         if(h3 .eq. h5) then
      tempw(2,3)=tempw(2,3)-2d0*ampsqfac/sqrt(colfac34_56)
     &   *dble(amp(uqsq_dqcq,h1,h2,h3,h5)
     & *dconjg(amp_swap(uqsq_dqcq,h1,h2,h3,h5)))
         endif
      endif

      enddo
      enddo
      enddo
      enddo
      temp(2,3)=temp(2,5)

C-----setup for dqcq_uqsq (1,4)-->(2,3)
      do h1=1,1
      do h2=1,1
      do h3=1,2
      do h5=1,2

      tempw(1,4)=tempw(1,4)+ampsqfac
     &   *dble(amp(dqcq_uqsq,h1,h2,h3,h5)
     & *dconjg(amp(dqcq_uqsq,h1,h2,h3,h5)))

      if (interference) then
      tempw(1,4)=tempw(1,4)+ampsqfac
     &   *dble(amp_swap(dqcq_uqsq,h1,h2,h3,h5)
     & *dconjg(amp_swap(dqcq_uqsq,h1,h2,h3,h5)))
         if(h3 .eq. h5) then
      tempw(1,4)=tempw(1,4)-2d0*ampsqfac/sqrt(colfac34_56)
     &   *dble(amp(dqcq_uqsq,h1,h2,h3,h5)
     & *dconjg(amp_swap(dqcq_uqsq,h1,h2,h3,h5)))
         endif
      endif

      enddo
      enddo
      enddo
      enddo

C-----setup for (dqcq_dqcq) (1,4)-->(1,4)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2

      temp(1,4)=temp(1,4)+ampsqfac
     &   *dble(amp(dqcq_dqcq,h1,h2,h3,h5)
     & *dconjg(amp(dqcq_dqcq,h1,h2,h3,h5)))

      if (interference) then
      temp(1,4)=temp(1,4)+ampsqfac
     &   *dble(amp_swap(dqcq_dqcq,h1,h2,h3,h5)
     & *dconjg(amp_swap(dqcq_dqcq,h1,h2,h3,h5)))
         if(h3 .eq. h5) then
      temp(1,4)=temp(1,4)-2d0*ampsqfac/sqrt(colfac34_56)
     &   *dble(amp(dqcq_dqcq,h1,h2,h3,h5)
     & *dconjg(amp_swap(dqcq_dqcq,h1,h2,h3,h5)))
         endif
      endif

      enddo
      enddo
      enddo
      enddo

C--------------------------------------------------------
C-----setup for dquq_dquq W diagrams (1,2)-->(1,2)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2

      temp(1,2)=temp(1,2)+ampsqfac
     &   *dble(ampa(dquq_dquq,h1,h2,h3,h5)
     & *dconjg(ampa(dquq_dquq,h1,h2,h3,h5)))
      temp(1,2)=temp(1,2)+ampsqfac
     &   *dble(ampb(dquq_dquq,h1,h2,h3,h5)
     & *dconjg(ampb(dquq_dquq,h1,h2,h3,h5)))
      if ((h1 .eq. 1) .and. (h2.eq. 1)) then
      temp(1,2)=temp(1,2)-2d0/xn*ampsqfac
     &   *dble(ampa(dquq_dquq,h1,h2,h3,h5)
     & *dconjg(ampb(dquq_dquq,h1,h2,h3,h5)))
      endif

      if (interference) then
      temp(1,2)=temp(1,2)+ampsqfac
     &   *dble(ampa_swap(dquq_dquq,h1,h2,h3,h5)
     & *dconjg(ampa_swap(dquq_dquq,h1,h2,h3,h5)))
      temp(1,2)=temp(1,2)+ampsqfac
     &   *dble(ampb_swap(dquq_dquq,h1,h2,h3,h5)
     & *dconjg(ampb_swap(dquq_dquq,h1,h2,h3,h5)))
      if ((h1 .eq. 1) .and. (h2.eq. 1)) then
      temp(1,2)=temp(1,2)-2d0/xn*ampsqfac
     &   *dble(ampa_swap(dquq_dquq,h1,h2,h3,h5)
     & *dconjg(ampb_swap(dquq_dquq,h1,h2,h3,h5)))
      endif
         if(h3 .eq. h5) then
      temp(1,2)=temp(1,2)-2d0*ampsqfac/sqrt(colfac34_56)
     &   *dble(ampa(dquq_dquq,h1,h2,h3,h5)
     & *dconjg(ampa_swap(dquq_dquq,h1,h2,h3,h5)))
      temp(1,2)=temp(1,2)-2d0*ampsqfac/sqrt(colfac34_56)
     &   *dble(ampb(dquq_dquq,h1,h2,h3,h5)
     & *dconjg(ampb_swap(dquq_dquq,h1,h2,h3,h5)))
      if ((h1 .eq. 1) .and. (h2.eq. 1)) then
      temp(1,2)=temp(1,2)+2d0/xn*ampsqfac/sqrt(colfac34_56)
     &   *dble(ampa(dquq_dquq,h1,h2,h3,h5)
     & *dconjg(ampb_swap(dquq_dquq,h1,h2,h3,h5)))
      temp(1,2)=temp(1,2)+2d0/xn*ampsqfac/sqrt(colfac34_56)
     &   *dble(ampa_swap(dquq_dquq,h1,h2,h3,h5)
     & *dconjg(ampb(dquq_dquq,h1,h2,h3,h5)))
      endif
         endif
      endif

      enddo
      enddo
      enddo
      enddo
      temp(3,4)=temp(1,2)
C-----------------------------------------------------------------

C-----setup for (dqsq_dqsq) (1,3)-->(1,3)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2

      temp(1,3)=temp(1,3)+ampsqfac
     &   *dble(amp(dqsq_dqsq,h1,h2,h3,h5)
     & *dconjg(amp(dqsq_dqsq,h1,h2,h3,h5)))

      if (interference) then
      temp(1,3)=temp(1,3)+ampsqfac
     &   *dble(amp_swap(dqsq_dqsq,h1,h2,h3,h5)
     & *dconjg(amp_swap(dqsq_dqsq,h1,h2,h3,h5)))
         if(h3 .eq. h5) then
      temp(1,3)=temp(1,3)-2d0*ampsqfac/sqrt(colfac34_56)
     &   *dble(amp(dqsq_dqsq,h1,h2,h3,h5)
     & *dconjg(amp_swap(dqsq_dqsq,h1,h2,h3,h5)))
         endif
      endif

      enddo
      enddo
      enddo
      enddo
      temp(1,5)=temp(1,3)
      temp(3,5)=temp(1,3)

      if ((j.eq.2).or.(j.eq.4).or.(j.eq.6).or.(j.eq.8)
     & .or.(j.eq.10).or.(j.eq.12)) go to 100
C-----setup for ((uquq_uquq)  (2,2)-->(2,2)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2

      temp(2,2)=temp(2,2)+ampsqfac
     & *dble(ampa(uquq_uquq,h1,h2,h3,h5)
     & *dconjg(ampa(uquq_uquq,h1,h2,h3,h5)))
      temp(2,2)=temp(2,2)+ampsqfac
     & *dble(ampb(uquq_uquq,h1,h2,h3,h5)
     & *dconjg(ampb(uquq_uquq,h1,h2,h3,h5)))
      if (h1 .eq. h2) then
      temp(2,2)=temp(2,2)-2d0/xn*ampsqfac
     & *dble(ampa(uquq_uquq,h1,h2,h3,h5)
     & *dconjg(ampb(uquq_uquq,h1,h2,h3,h5)))
      endif

      if (interference) then
      temp(2,2)=temp(2,2)+ampsqfac
     & *dble(ampa_swap(uquq_uquq,h1,h2,h3,h5)
     & *dconjg(ampa_swap(uquq_uquq,h1,h2,h3,h5)))
      temp(2,2)=temp(2,2)+ampsqfac
     & *dble(ampb_swap(uquq_uquq,h1,h2,h3,h5)
     & *dconjg(ampb_swap(uquq_uquq,h1,h2,h3,h5)))
      if (h1 .eq. h2) then
      temp(2,2)=temp(2,2)-2d0/xn*ampsqfac
     & *dble(ampa_swap(uquq_uquq,h1,h2,h3,h5)
     & *dconjg(ampb_swap(uquq_uquq,h1,h2,h3,h5)))
      endif
         if(h3 .eq. h5) then
      temp(2,2)=temp(2,2)-2d0*ampsqfac/sqrt(colfac34_56)
     & *dble(ampa(uquq_uquq,h1,h2,h3,h5)
     & *dconjg(ampa_swap(uquq_uquq,h1,h2,h3,h5)))
      temp(2,2)=temp(2,2)-2d0*ampsqfac/sqrt(colfac34_56)
     & *dble(ampb(uquq_uquq,h1,h2,h3,h5)
     & *dconjg(ampb_swap(uquq_uquq,h1,h2,h3,h5)))
      if (h1 .eq. h2) then
      temp(2,2)=temp(2,2)+2d0/xn*ampsqfac/sqrt(colfac34_56)
     & *dble(ampa(uquq_uquq,h1,h2,h3,h5)
     & *dconjg(ampb_swap(uquq_uquq,h1,h2,h3,h5)))
      temp(2,2)=temp(2,2)+2d0/xn*ampsqfac/sqrt(colfac34_56)
     & *dble(ampa_swap(uquq_uquq,h1,h2,h3,h5)
     & *dconjg(ampb(uquq_uquq,h1,h2,h3,h5)))
      endif
         endif
      endif

      enddo
      enddo
      enddo
      enddo
      temp(4,4)=temp(2,2)


C-----setup for ((dqdq_dqdq)  (1,1)-->(1,1)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2

      temp(1,1)=temp(1,1)+ampsqfac
     & *dble(ampa(dqdq_dqdq,h1,h2,h3,h5)
     & *dconjg(ampa(dqdq_dqdq,h1,h2,h3,h5)))
      temp(1,1)=temp(1,1)+ampsqfac
     & *dble(ampb(dqdq_dqdq,h1,h2,h3,h5)
     & *dconjg(ampb(dqdq_dqdq,h1,h2,h3,h5)))
      if (h1 .eq. h2) then
      temp(1,1)=temp(1,1)-2d0/xn*ampsqfac
     & *dble(ampa(dqdq_dqdq,h1,h2,h3,h5)
     & *dconjg(ampb(dqdq_dqdq,h1,h2,h3,h5)))
      endif

      if (interference) then
      temp(1,1)=temp(1,1)+ampsqfac
     & *dble(ampa_swap(dqdq_dqdq,h1,h2,h3,h5)
     & *dconjg(ampa_swap(dqdq_dqdq,h1,h2,h3,h5)))
      temp(1,1)=temp(1,1)+ampsqfac
     & *dble(ampb_swap(dqdq_dqdq,h1,h2,h3,h5)
     & *dconjg(ampb_swap(dqdq_dqdq,h1,h2,h3,h5)))
      if (h1 .eq. h2) then
      temp(1,1)=temp(1,1)-2d0/xn*ampsqfac
     & *dble(ampa_swap(dqdq_dqdq,h1,h2,h3,h5)
     & *dconjg(ampb_swap(dqdq_dqdq,h1,h2,h3,h5)))
      endif
         if(h3 .eq. h5) then
      temp(1,1)=temp(1,1)-2d0*ampsqfac/sqrt(colfac34_56)
     & *dble(ampa(dqdq_dqdq,h1,h2,h3,h5)
     & *dconjg(ampa_swap(dqdq_dqdq,h1,h2,h3,h5)))
      temp(1,1)=temp(1,1)-2d0*ampsqfac/sqrt(colfac34_56)
     & *dble(ampb(dqdq_dqdq,h1,h2,h3,h5)
     & *dconjg(ampb_swap(dqdq_dqdq,h1,h2,h3,h5)))
      if (h1 .eq. h2) then
      temp(1,1)=temp(1,1)+2d0/xn*ampsqfac/sqrt(colfac34_56)
     & *dble(ampa(dqdq_dqdq,h1,h2,h3,h5)
     & *dconjg(ampb_swap(dqdq_dqdq,h1,h2,h3,h5)))
      temp(1,1)=temp(1,1)+2d0/xn*ampsqfac/sqrt(colfac34_56)
     & *dble(ampa_swap(dqdq_dqdq,h1,h2,h3,h5)
     & *dconjg(ampb(dqdq_dqdq,h1,h2,h3,h5)))
      endif
         endif
      endif

      enddo
      enddo
      enddo
      enddo
      temp(3,3)=temp(1,1)
      temp(5,5)=temp(1,1)
 100  continue

      do k=1,nf;do l=1,nf
      if (
     & .not.(
     &    (
     &    (pid_pdg(j1(j)).eq.0)
     &    .or. (
     &    j1(j).le.2 .and. (
     &    pid_pdg(j1(j)).eq.k
     &    .or.
     &    (
     &    modulo(k,2).eq.0 .and.
     &    (pid_pdg(j1(j)).eq.12
     &    .or. pid_pdg(j1(j)).eq.14
     &    .or. pid_pdg(j1(j)).eq.16)
     &    )
     &    .or.
     &    (
     &    modulo(k,2).eq.1 .and.
     &    (pid_pdg(j1(j)).eq.11
     &    .or. pid_pdg(j1(j)).eq.13
     &    .or. pid_pdg(j1(j)).eq.15)
     &    )
     &    )
     &    )
     &    .or. (
     &    j1(j).ge.7 .and. (
     &    pid_pdg(j1(j)).eq.-k
     &    .or.
     &    (
     &    modulo(k,2).eq.0 .and.
     &    (pid_pdg(j1(j)).eq.-12
     &    .or. pid_pdg(j1(j)).eq.-14
     &    .or. pid_pdg(j1(j)).eq.-16)
     &    )
     &    .or.
     &    (
     &    modulo(k,2).eq.1 .and.
     &    (pid_pdg(j1(j)).eq.-11
     &    .or. pid_pdg(j1(j)).eq.-13
     &    .or. pid_pdg(j1(j)).eq.-15)
     &    )
     &    )
     &    )
     &    ) .and. (
     &    (pid_pdg(j2(j)).eq.0)
     &    .or. (
     &    j2(j).le.2 .and. (
     &    pid_pdg(j2(j)).eq.l
     &    .or.
     &    (
     &    modulo(l,2).eq.0 .and.
     &    (pid_pdg(j2(j)).eq.12
     &    .or. pid_pdg(j2(j)).eq.14
     &    .or. pid_pdg(j2(j)).eq.16)
     &    )
     &    .or.
     &    (
     &    modulo(l,2).eq.1 .and.
     &    (pid_pdg(j2(j)).eq.11
     &    .or. pid_pdg(j2(j)).eq.13
     &    .or. pid_pdg(j2(j)).eq.15)
     &    )
     &    )
     &    )
     &    .or. (
     &    j2(j).ge.7 .and. (
     &    pid_pdg(j2(j)).eq.-l
     &    .or.
     &    (
     &    modulo(l,2).eq.0 .and.
     &    (pid_pdg(j2(j)).eq.-12
     &    .or. pid_pdg(j2(j)).eq.-14
     &    .or. pid_pdg(j2(j)).eq.-16)
     &    )
     &    .or.
     &    (
     &    modulo(l,2).eq.1 .and.
     &    (pid_pdg(j2(j)).eq.-11
     &    .or. pid_pdg(j2(j)).eq.-13
     &    .or. pid_pdg(j2(j)).eq.-15)
     &    )
     &    )
     &    )
     &    )
     & )
     & ) then
         temp(k,l)=zip
      endif
      enddo;enddo
      if (
     & isANeutrino(abs(pid_pdg(7))) .and.
     & abs(pid_pdg(8)).eq.abs(pid_pdg(7))
     & ) then
         tempw(:,:)=zip
         if(j.eq.10 .or. j.eq.12) then
            temp(4,5)=zip
         endif
      else if (
     & isALepton(abs(pid_pdg(7))) .and.
     & abs(pid_pdg(8)).eq.abs(pid_pdg(7))
     & ) then
         tempw(:,:)=zip
         if (j.eq.9 .or. j.eq.11) then
            temp(1:2,3)=zip
         else if (j.eq.10 .or. j.eq.12) then
            temp(1,3:5)=zip
         endif
      else if (
     & (isALepton(abs(pid_pdg(7))) .and.
     & isANeutrino(abs(pid_pdg(8)))) .or.
     & (isALepton(abs(pid_pdg(8))) .and.
     & isANeutrino(abs(pid_pdg(7))))
     & ) then
         temp(:,:)=zip
      endif

c---- Multiply by the decay color factor
      !temp(:,:) = temp(:,:)*colfac34_56*vsymfact
      !tempw(:,:) = tempw(:,:)*colfac34_56*vsymfact
      temp(:,:) = temp(:,:)*vsymfact
      tempw(:,:) = tempw(:,:)*vsymfact

c--- fill matrix elements
      if (j.eq.1) then
      do k=1,nf
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.k) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.k)
     & )
     & ) then
         call addtemptomsq(msq,temp,k,k,k,k,j,stat)
      endif
      do l=k+1,nf
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.k) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.l)
     & )
     & ) then
         call addtemptomsq(msq,temp,k,l,k,l,j,1d0)
      endif
      enddo
      enddo
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.4) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.1)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.1) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.4)
     & )
     & ) then
         call addtempwtomsq(msq,tempw,2,3,2,3,j,1d0)
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.3) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.2)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.2) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.3)
     & )
     & ) then
         call addtempwtomsq(msq,tempw,1,4,1,4,j,1d0)
      endif

      elseif (j.eq.2) then
      do k=1,nf
      do l=k+1,nf
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.k) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.l)
     & )
     & ) then
         call addtemptomsq(msq,temp,l,k,k,l,j,1d0)
      endif
      enddo
      enddo
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.4) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.1)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.1) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.4)
     & )
     & ) then
         call addtempwtomsq(msq,tempw,3,2,2,3,j,1d0)
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.3) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.2)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.2) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.3)
     & )
     & ) then
         call addtempwtomsq(msq,tempw,4,1,1,4,j,1d0)
      endif

      elseif (j.eq.3) then
      do k=-nf,-1
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.k) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.k)
     & )
     & ) then
         call addtemptomsq(msq,temp,k,k,-k,-k,j,stat)
      endif
      do l=k+1,-1
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.k) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.l)
     & )
     & ) then
         call addtemptomsq(msq,temp,k,l,-l,-k,j,1d0)
      endif
      enddo
      enddo
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-4) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-1)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-1) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-4)
     & )
     & ) then
         call addtempwtomsq(msq,tempw,-3,-2,1,4,j,1d0)
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-3) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-2)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-2) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-3)
     & )
     & ) then
         call addtempwtomsq(msq,tempw,-4,-1,2,3,j,1d0)
      endif

      elseif (j.eq.4) then
      do k=-nf,-1
      do l=k+1,-1
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.k) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.l)
     & )
     & ) then
         call addtemptomsq(msq,temp,l,k,-l,-k,j,1d0)
      endif
      enddo
      enddo
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-4) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-1)
     & )
     v .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-1) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-4)
     & )
     & ) then
         call addtempwtomsq(msq,tempw,-2,-3,1,4,j,1d0)
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-3) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-2)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-2) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-3)
     & )
     & ) then
         call addtempwtomsq(msq,tempw,-1,-4,2,3,j,1d0)
      endif

c--- qbar-q
      elseif (j.eq.5) then
      do k=-nf,-1
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.k) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-k)
     & ) .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-k) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.k)
     & )
     & ) then
         call addtemptomsq(msq,temp,k,-k,-k,-k,j,1d0)
      endif
      do l=1,nf
      if (abs(k) .lt. abs(l)) then
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.k) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.l)
     & )
     & ) then
         call addtemptomsq(msq,temp,k,l,-k,l,j,1d0)
      endif
      endif
      enddo
      enddo
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.4) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-2)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-2) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.4)
     & )
     & ) then
         call addtempwtomsq(msq,tempw,-1,3,2,3,j,1d0)
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.3) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-1)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-1) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.3)
     & )
     & ) then
         call addtempwtomsq(msq,tempw,-2,4,1,4,j,1d0)
      endif

c--- qbar-q
      elseif (j.eq.6) then
      do k=-nf,-1
      do l=1,nf
      if (abs(k) .gt. abs(l)) then
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.k) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.l)
     & )
     & ) then
         call addtemptomsq(msq,temp,k,l,l,-k,j,1d0)
      endif
      endif
      enddo
      enddo
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-4) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.2)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.2) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-4)
     & )
     & ) then
         call addtempwtomsq(msq,tempw,-3,1,1,4,j,1d0)
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-3) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.1)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.1) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-3)
     & )
     & ) then
         call addtempwtomsq(msq,tempw,-4,2,2,3,j,1d0)
      endif

c--- q-qbar
      elseif (j.eq.7) then
      do k=-nf,-1
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.k) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-k)
     & ) .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-k) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.k)
     & )
     & ) then
         call addtemptomsq(msq,temp,-k,k,-k,-k,j,1d0)
      endif
      do l=1,nf
      if (abs(k) .lt. abs(l)) then
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.k) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.l)
     & )
     & ) then
         call addtemptomsq(msq,temp,l,k,-k,l,j,1d0)
      endif
      endif
      enddo
      enddo
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.4) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-2)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-2) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.4)
     & )
     & ) then
         call addtempwtomsq(msq,tempw,3,-1,2,3,j,1d0)
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.3) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-1)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-1) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.3)
     & )
     & ) then
         call addtempwtomsq(msq,tempw,4,-2,1,4,j,1d0)
      endif

c--- q-qbar
      elseif (j.eq.8) then
      do k=-nf,-1
      do l=-nf,-1
      if (abs(k) .lt. abs(l)) then
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.l) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-k)
     & )
     & ) then
         call addtemptomsq(msq,temp,-k,l,-k,-l,j,1d0)
      endif
      endif
      enddo
      enddo
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-4) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.2)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.2) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-4)
     & )
     & ) then
         call addtempwtomsq(msq,tempw,1,-3,1,4,j,1d0)
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-3) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.1)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.1) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-3)
     & )
     & ) then
         call addtempwtomsq(msq,tempw,2,-4,2,3,j,1d0)
      endif

c--- q-qbar extra pieces
      elseif (j.eq.9) then
      do k=1,nf
      do l=1,nf
      if (k .lt. l) then
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.l) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-l)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-l) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.l)
     & )
     & .or.
     & (
     &    abs(pid_pdg(8)).eq.abs(pid_pdg(7)) .and. (
     &    (modulo(l,2).eq.0 .and. isANeutrino(abs(pid_pdg(7))))
     &    .or.
     &    (modulo(l,2).eq.1 .and. isALepton(abs(pid_pdg(7))))
     &    )
     & )
     & ) then
         call addtemptomsq(msq,temp,k,-k,k,l,j,1d0)
      endif
      endif
      enddo
      enddo
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-4
     &                  .or. isANeutrino(-pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.3
     &                  .or. isALepton(pid_pdg(8)))
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.3
     &                  .or. isALepton(pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-4
     &                  .or. isANeutrino(-pid_pdg(8)))
     & )
     & ) then
         call addtempwtomsq(msq,tempw,1,-2,1,4,j,1d0)
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-2
     &                  .or. isANeutrino(-pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.1
     &                  .or. isALepton(pid_pdg(8)))
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.1
     &                  .or. isALepton(pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-2
     &                  .or. isANeutrino(-pid_pdg(8)))
     & )
     & ) then
         call addtempwtomsq(msq,tempw,3,-4,1,4,j,1d0)
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.4
     &                  .or. isANeutrino(pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-3
     &                  .or. isALepton(-pid_pdg(8)))
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-3
     &                  .or. isALepton(-pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.4
     &                  .or. isANeutrino(pid_pdg(8)))
     & )
     & ) then
         call addtempwtomsq(msq,tempw,2,-1,2,3,j,1d0)
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.2
     &                  .or. isANeutrino(pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-1
     &                  .or. isALepton(-pid_pdg(8)))
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-1
     &                  .or. isALepton(-pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.2
     &                  .or. isANeutrino(pid_pdg(8)))
     & )
     & ) then
         call addtempwtomsq(msq,tempw,4,-3,2,3,j,1d0)
      endif

c--- q-qbar extra pieces
      elseif (j.eq.10) then
      do k=1,nf
      do l=1,nf
      if (k .gt. l) then
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.l) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-l)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-l) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.l)
     & )
     & .or.
     & (
     &    abs(pid_pdg(8)).eq.abs(pid_pdg(7)) .and. (
     &    (modulo(l,2).eq.0 .and. isANeutrino(abs(pid_pdg(7))))
     &    .or.
     &    (modulo(l,2).eq.1 .and. isALepton(abs(pid_pdg(7))))
     &    )
     & )
     & ) then
         call addtemptomsq(msq,temp,k,-k,l,k,j,1d0)
      endif
      endif
      enddo
      enddo

c--- qbar-q extra pieces
      elseif (j.eq.11) then
      do k=1,nf
      do l=1,nf
      if (k .lt. l) then
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.l) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-l)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-l) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.l)
     & )
     & .or.
     & (
     &    abs(pid_pdg(8)).eq.abs(pid_pdg(7)) .and. (
     &    (modulo(l,2).eq.0 .and. isANeutrino(abs(pid_pdg(7))))
     &    .or.
     &    (modulo(l,2).eq.1 .and. isALepton(abs(pid_pdg(7))))
     &    )
     & )
     & ) then
         call addtemptomsq(msq,temp,-k,k,k,l,j,1d0)
      endif
      endif
      enddo
      enddo
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-4
     &                  .or. isANeutrino(-pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.3
     &                  .or. isALepton(pid_pdg(8)))
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.3
     &                  .or. isALepton(pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-4
     &                  .or. isANeutrino(-pid_pdg(8)))
     & )
     & ) then
         call addtempwtomsq(msq,tempw,-2,1,1,4,j,1d0)
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-2
     &                  .or. isANeutrino(-pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.1
     &                  .or. isALepton(pid_pdg(8)))
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.1
     &                  .or. isALepton(pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-2
     &                  .or. isANeutrino(-pid_pdg(8)))
     & )
     & ) then
         call addtempwtomsq(msq,tempw,-4,3,1,4,j,1d0)
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.4
     &                  .or. isANeutrino(pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-3
     &                  .or. isALepton(-pid_pdg(8)))
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-3
     &                  .or. isALepton(-pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.4
     &                  .or. isANeutrino(pid_pdg(8)))
     & )
     & ) then
         call addtempwtomsq(msq,tempw,-1,2,2,3,j,1d0)
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.2
     &                  .or. isANeutrino(pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-1
     &                  .or. isALepton(-pid_pdg(8)))
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-1
     &                  .or. isALepton(-pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.2
     &                  .or. isANeutrino(pid_pdg(8)))
     & )
     & ) then
         call addtempwtomsq(msq,tempw,-3,4,2,3,j,1d0)
      endif

c--- qbar-q extra pieces
      elseif (j.eq.12) then
      do k=1,nf
      do l=1,nf
      if (k .gt. l) then
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.l) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-l)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-l) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.l)
     & )
     & .or.
     & (
     &    abs(pid_pdg(8)).eq.abs(pid_pdg(7)) .and. (
     &    (modulo(l,2).eq.0 .and. isANeutrino(abs(pid_pdg(7))))
     &    .or.
     &    (modulo(l,2).eq.1 .and. isALepton(abs(pid_pdg(7))))
     &    )
     & )
     & ) then
         call addtemptomsq(msq,temp,-k,k,l,k,j,1d0)
      endif
      endif
      enddo
      enddo

      endif

      enddo

      return

   77 format(' *      W-mass^2     (',f11.5,',',f11.5,')      *')
   78 format(' *      Z-mass^2     (',f11.5,',',f11.5,')      *')
   79 format(' *  sin^2(theta_w)   (',f11.5,',',f11.5,')      *')

      end
