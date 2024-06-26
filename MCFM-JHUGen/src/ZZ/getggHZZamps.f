c---
c--- MODIFICATION OF THE ORIGINAL MCFM SUBROUTINE TO ALLOW FOR ANOMALOUS H-Z-Z COUPLINGS  (e.g. nproc=128)
c--- SAME CHOICE OF CONVENTIONS AS IN JHUGEN
c---
      subroutine getggHZZamps(p,Mloop_bquark,Mloop_tquark)
c--- Returns a series of arrays representing the dressed amp[itudes
c--- for the process gg->Higgs->ZZ; there are:
c---        Mloop_bquark(h1,h2,h34,h56)   top quark mass=mt
c---        Mloop_tquark(h1,h2,h34,h56)   bottom quark mass=mb
c---
c--- The overall factor on the amplitude is:
c---
c---      4d0*esq*gsq/(16d0*pisq)*esq * delta(a,b)
c---
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'scale.f'
      include 'anom_higgs.f'
      include 'spinzerohiggs_anomcoupl.f'
      include 'plabel.f'
      include 'AnomZffCouplings.f'
      integer h1,h34,h56
      double precision p(mxpart,4),mb2,mt2,mtX2,mbX2
      double complex Mloop_bquark(2,2,2,2),Mloop_tquark(2,2,2,2),
     & ggHmq(2,2,2),prop12,prop34,prop56,
     & H4l(2,2),facHZZ,facHZA,facHAZ,facHAA,higgsprop
      double precision rescale
      double complex anomhzzamp,anomhzaamp,anomhaaamp

!==== for width studies rescale by appropriate factor
      if((keep_smhiggs_norm).and.(anom_higgs)) then
         rescale=chi_higgs**2
      else
         rescale=1d0
      endif

      Mloop_bquark(:,:,:,:)=czip
      Mloop_tquark(:,:,:,:)=czip
      if(hmass.lt.zip) then
         return
      endif
      ggHmq(:,:,:)=czip
      H4l(:,:)=czip

      call spinoru(6,p,za,zb)

c--- propagator factors
      prop12=higgsprop(s(1,2))
      prop34=cone/dcmplx(s(3,4)-zmass**2,zmass*zwidth)
      prop56=cone/dcmplx(s(5,6)-zmass**2,zmass*zwidth)

c--- Factor
      facHZZ=im*rescale*prop12*prop34*prop56/(2d0*xw*(1d0-xw))
      facHZA=-im*rescale*prop12*prop34/s(5,6)/(2d0*xw*(1d0-xw))
      facHAZ=-im*rescale*prop12/s(3,4)*prop56/(2d0*xw*(1d0-xw))
      facHAA=im*rescale*prop12/s(3,4)/s(5,6)/(2d0*xw*(1d0-xw))

c--- Amplitudes for production
      call anomhggvtxamp(1,2,1,za,zb,ggHmq)
      ! Overall factor=1
      !ggHmq(:,:,:) = ggHmq(:,:,:)

c--- Setting Anomalous Zff Couplings 
      if (AllowAnomalousZffCouplings .eq. 1) then
        if ((plabel(3) .eq. 'el') .or. (plabel(3) .eq. 'ml')
     &.or. (plabel(3) .eq. 'tl')) then
          l1 = leZ
          r1 = reZ 
        elseif (plabel(3) .eq. 'nl') then
          l1 = lnZ*dsqrt(3d0)
          r1 = rnZ*dsqrt(3d0) 
        elseif ((plabel(5) .eq. 'bq') .or. (plabel(5) .eq. 'sq')
     &.or. (plabel(5) .eq. 'dq')) then
          l1=lqdZ*dsqrt(3d0)
          r1=rqdZ*dsqrt(3d0)
        elseif ((plabel(5) .eq. 'uq') .or. (plabel(5) .eq. 'cq')) then
          l1=lquZ*dsqrt(3d0)
          r1=rquZ*dsqrt(3d0)
        endif 
      endif
      if (AllowAnomalousZffCouplings .eq. 1) then
        if ((plabel(5) .eq. 'el') .or. (plabel(5) .eq. 'ml')
     &.or. (plabel(5) .eq. 'tl')) then
          l2 = leZ
          r2 = reZ 
        elseif (plabel(5) .eq. 'nl') then
          l2 = lnZ*dsqrt(3d0)
          r2 = rnZ*dsqrt(3d0)
        elseif ((plabel(5) .eq. 'bq') .or. (plabel(5) .eq. 'sq')
     &.or. (plabel(5) .eq. 'dq')) then
            l2=lqdZ*dsqrt(3d0)
            r2=rqdZ*dsqrt(3d0)
        elseif ((plabel(5) .eq. 'uq') .or. (plabel(5) .eq. 'cq')) then
            l2=lquZ*dsqrt(3d0)
            r2=rquZ*dsqrt(3d0)  
        endif
      endif

c--- Amplitudes for decay
      H4l(1,1)=
     &  anomhzzamp(3,4,5,6,1,s(1,2),s(3,4),s(5,6),za,zb)*l1*l2*facHZZ
     & +anomhzaamp(3,4,5,6,1,s(1,2),s(3,4),s(5,6),za,zb)*l1*q2*facHZA
     & +anomhzaamp(5,6,3,4,1,s(1,2),s(5,6),s(3,4),za,zb)*q1*l2*facHAZ
     & +anomhaaamp(3,4,5,6,1,s(1,2),s(3,4),s(5,6),za,zb)*q1*q2*facHAA
      H4l(2,1)=
     &  anomhzzamp(4,3,5,6,1,s(1,2),s(3,4),s(5,6),za,zb)*r1*l2*facHZZ
     & +anomhzaamp(4,3,5,6,1,s(1,2),s(3,4),s(5,6),za,zb)*r1*q2*facHZA
     & +anomhzaamp(5,6,4,3,1,s(1,2),s(5,6),s(3,4),za,zb)*q1*l2*facHAZ
     & +anomhaaamp(4,3,5,6,1,s(1,2),s(3,4),s(5,6),za,zb)*q1*q2*facHAA
      H4l(1,2)=
     &  anomhzzamp(3,4,6,5,1,s(1,2),s(3,4),s(5,6),za,zb)*l1*r2*facHZZ
     & +anomhzaamp(3,4,6,5,1,s(1,2),s(3,4),s(5,6),za,zb)*l1*q2*facHZA
     & +anomhzaamp(6,5,3,4,1,s(1,2),s(5,6),s(3,4),za,zb)*q1*r2*facHAZ
     & +anomhaaamp(3,4,6,5,1,s(1,2),s(3,4),s(5,6),za,zb)*q1*q2*facHAA
      H4l(2,2)=
     &  anomhzzamp(4,3,6,5,1,s(1,2),s(3,4),s(5,6),za,zb)*r1*r2*facHZZ
     & +anomhzaamp(4,3,6,5,1,s(1,2),s(3,4),s(5,6),za,zb)*r1*q2*facHZA
     & +anomhzaamp(6,5,4,3,1,s(1,2),s(5,6),s(3,4),za,zb)*q1*r2*facHAZ
     & +anomhaaamp(4,3,6,5,1,s(1,2),s(3,4),s(5,6),za,zb)*q1*q2*facHAA

c--- Assemble
      do h1=1,2
      do h34=1,2
      do h56=1,2
      Mloop_bquark(h1,h1,h34,h56)=ggHmq(1,h1,h1)*H4l(h34,h56)
      Mloop_tquark(h1,h1,h34,h56)=ggHmq(2,h1,h1)*H4l(h34,h56)
      enddo
      enddo
      enddo

      return
      end





      subroutine getggH2ZZamps(p,Mloop_bquark,Mloop_tquark)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'scale.f'
      include 'anom_higgs.f'
      include 'spinzerohiggs_anomcoupl.f'
      include 'plabel.f'
      include 'AnomZffCouplings.f'
      integer h1,h34,h56
      double precision p(mxpart,4),mb2,mt2,mbX2,mtX2
      double complex Mloop_bquark(2,2,2,2),Mloop_tquark(2,2,2,2),
     & ggHmq(2,2,2),prop12,prop34,prop56,
     & H4l(2,2),facHZZ,facHZA,facHAZ,facHAA,higgs2prop
      double precision rescale
      double complex anomhzzamp,anomhzaamp,anomhaaamp

!==== for width studies rescale by appropriate factor
      if((keep_smhiggs_norm).and.(anom_higgs)) then
         rescale=chi_higgs**2
      else
         rescale=1d0
      endif

      Mloop_bquark(:,:,:,:)=czip
      Mloop_tquark(:,:,:,:)=czip
      if(h2mass.lt.zip) then
         return
      endif
      ggHmq(:,:,:)=czip
      H4l(:,:)=czip

      call spinoru(6,p,za,zb)

c--- propagator factors
      prop12=higgs2prop(s(1,2))
      prop34=cone/dcmplx(s(3,4)-zmass**2,zmass*zwidth)
      prop56=cone/dcmplx(s(5,6)-zmass**2,zmass*zwidth)

c--- Factor
      facHZZ=im*rescale*prop12*prop34*prop56/(2d0*xw*(1d0-xw))
      facHZA=-im*rescale*prop12*prop34/s(5,6)/(2d0*xw*(1d0-xw))
      facHAZ=-im*rescale*prop12/s(3,4)*prop56/(2d0*xw*(1d0-xw))
      facHAA=im*rescale*prop12/s(3,4)/s(5,6)/(2d0*xw*(1d0-xw))

c--- Amplitudes for production
      call anomhggvtxamp(1,2,2,za,zb,ggHmq)
      ! Overall factor=1
      !ggHmq(:,:,:) = ggHmq(:,:,:)

c--- Setting Anomalous Zff Couplings 
      if (AllowAnomalousZffCouplings .eq. 1) then
        if ((plabel(3) .eq. 'el') .or. (plabel(3) .eq. 'ml')
     &.or. (plabel(3) .eq. 'tl')) then
          l1 = leZ
          r1 = reZ 
        elseif (plabel(3) .eq. 'nl') then
          l1 = lnZ*dsqrt(3d0)
          r1 = rnZ*dsqrt(3d0) 
        elseif ((plabel(5) .eq. 'bq') .or. (plabel(5) .eq. 'sq')
     &.or. (plabel(5) .eq. 'dq')) then
          l1=lqdZ*dsqrt(3d0)
          r1=rqdZ*dsqrt(3d0)
        elseif ((plabel(5) .eq. 'uq') .or. (plabel(5) .eq. 'cq')) then
          l1=lquZ*dsqrt(3d0)
          r1=rquZ*dsqrt(3d0)
        endif 
      endif
      if (AllowAnomalousZffCouplings .eq. 1) then
        if ((plabel(5) .eq. 'el') .or. (plabel(5) .eq. 'ml')
     &.or. (plabel(5) .eq. 'tl')) then
          l2 = leZ
          r2 = reZ 
        elseif (plabel(5) .eq. 'nl') then
          l2 = lnZ*dsqrt(3d0)
          r2 = rnZ*dsqrt(3d0)
        elseif ((plabel(5) .eq. 'bq') .or. (plabel(5) .eq. 'sq')
     &.or. (plabel(5) .eq. 'dq')) then
            l2=lqdZ*dsqrt(3d0)
            r2=rqdZ*dsqrt(3d0)
        elseif ((plabel(5) .eq. 'uq') .or. (plabel(5) .eq. 'cq')) then
            l2=lquZ*dsqrt(3d0)
            r2=rquZ*dsqrt(3d0)  
        endif
      endif

c--- Amplitudes for decay
      H4l(1,1)=
     &  anomhzzamp(3,4,5,6,2,s(1,2),s(3,4),s(5,6),za,zb)*l1*l2*facHZZ
     & +anomhzaamp(3,4,5,6,2,s(1,2),s(3,4),s(5,6),za,zb)*l1*q2*facHZA
     & +anomhzaamp(5,6,3,4,2,s(1,2),s(5,6),s(3,4),za,zb)*q1*l2*facHAZ
     & +anomhaaamp(3,4,5,6,2,s(1,2),s(3,4),s(5,6),za,zb)*q1*q2*facHAA
      H4l(2,1)=
     &  anomhzzamp(4,3,5,6,2,s(1,2),s(3,4),s(5,6),za,zb)*r1*l2*facHZZ
     & +anomhzaamp(4,3,5,6,2,s(1,2),s(3,4),s(5,6),za,zb)*r1*q2*facHZA
     & +anomhzaamp(5,6,4,3,2,s(1,2),s(5,6),s(3,4),za,zb)*q1*l2*facHAZ
     & +anomhaaamp(4,3,5,6,2,s(1,2),s(3,4),s(5,6),za,zb)*q1*q2*facHAA
      H4l(1,2)=
     &  anomhzzamp(3,4,6,5,2,s(1,2),s(3,4),s(5,6),za,zb)*l1*r2*facHZZ
     & +anomhzaamp(3,4,6,5,2,s(1,2),s(3,4),s(5,6),za,zb)*l1*q2*facHZA
     & +anomhzaamp(6,5,3,4,2,s(1,2),s(5,6),s(3,4),za,zb)*q1*r2*facHAZ
     & +anomhaaamp(3,4,6,5,2,s(1,2),s(3,4),s(5,6),za,zb)*q1*q2*facHAA
      H4l(2,2)=
     &  anomhzzamp(4,3,6,5,2,s(1,2),s(3,4),s(5,6),za,zb)*r1*r2*facHZZ
     & +anomhzaamp(4,3,6,5,2,s(1,2),s(3,4),s(5,6),za,zb)*r1*q2*facHZA
     & +anomhzaamp(6,5,4,3,2,s(1,2),s(5,6),s(3,4),za,zb)*q1*r2*facHAZ
     & +anomhaaamp(4,3,6,5,2,s(1,2),s(3,4),s(5,6),za,zb)*q1*q2*facHAA

c--- Assemble
      do h1=1,2
      do h34=1,2
      do h56=1,2
      Mloop_bquark(h1,h1,h34,h56)=ggHmq(1,h1,h1)*H4l(h34,h56)
      Mloop_tquark(h1,h1,h34,h56)=ggHmq(2,h1,h1)*H4l(h34,h56)
      enddo
      enddo
      enddo

      return
      end

