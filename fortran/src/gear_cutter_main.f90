! --------------------------------------------------------------------------
!
! GearCutter 0.2
! Generation of involute gear tooth profiles
! Copyright (C) 2009-2024 Konstantinos Poulios
!
! --------------------------------------------------------------------------
!
! This file is part of GearCutter.
!
! GearCutter is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! --------------------------------------------------------------------------

!***************************************************
Module Konstanten

   !------ Konstanten ------------------------------------
   Real(8),Parameter    :: PI=3.14159265358979

End Module Konstanten
!***************************************************

!***************************************************
Subroutine get_io_unit(iunit)

   Implicit None

   Integer, Intent(out) :: iunit
   Integer  :: i, ios
   Logical  :: lopen

   iunit = 0
   Do i = 101, 999
      Inquire(UNIT=i,OPENED=lopen,IOSTAT=ios)
      If (ios==0) Then
         If (.not.lopen) Then
            iunit = i
            Return
         End if
      End if
   End do

End Subroutine get_io_unit
!***************************************************

!***************************************************
Module Strukturen

   !------ Strukturen ------------------------------------
   Type CSLage
      Real(8),Pointer   :: X(:),Y(:),fi(:)
   End Type CSLage

   Type Bogen
      Real(8)  :: xWZ_M,yWZ_M,rho,fiWZ_1,fiWZ_2    ,&
                  xSTCK_M,ySTCK_M,fiSTCK_1,fiSTCK_2,r_M,fi_M
   End Type Bogen

   Type Gerade
      Real(8)  :: xWZ_1,   yWZ_1,   xWZ_2,   yWZ_2,   fiWZ_12     ,&
                  xSTCK_1, ySTCK_1, xSTCK_2, ySTCK_2, fiSTCK_12   ,&
                  r1,r2
   End Type Gerade

   Type Kontur
      Integer  :: id
      Logical  :: flBog
      Type(Bogen),Pointer  :: Bogen
      Type(Gerade),Pointer :: Gerade
      Type(Kontur),Pointer :: next
   End Type Kontur

End Module Strukturen
!***************************************************

!***************************************************
Module Prozeduren

Contains
   !***************************************************
   Subroutine WZzuSTCK_CS(idx,WZ,STCK,xWZ,yWZ,xSTCK,ySTCK)

      Use Strukturen
      Implicit None
      !------------------------------------------------
      Integer,Intent(in)      :: idx
      Type(CSLage),Intent(in) :: WZ,STCK
      Real(8),Intent(in)      :: xWZ,yWZ
      Real(8),Intent(out)     :: xSTCK,ySTCK
      !------------------------------------------------
      Real(8)  :: xUrspr,yUrspr,myCos,mySin
      !------------------------------------------------

      myCos=dCos(WZ%fi(idx))
      mySin=dSin(WZ%fi(idx))

      xUrspr=xWZ*myCos-yWZ*mySin
      yUrspr=xWZ*mySin+yWZ*myCos

      xUrspr=xUrspr+WZ%X(idx)-STCK%X(idx)
      yUrspr=yUrspr+WZ%Y(idx)-STCK%Y(idx)

      myCos=dCos(-STCK%fi(idx))
      mySin=dSin(-STCK%fi(idx))

      xSTCK=xUrspr*myCos-yUrspr*mySin
      ySTCK=xUrspr*mySin+yUrspr*myCos

   End Subroutine WZzuSTCK_CS
   !***************************************************

   !***************************************************
   Subroutine Winkelkorr(Winkel)

      Use Konstanten
      Implicit None
      !------------------------------------------------
      Real(8)  :: Winkel
      !------------------------------------------------

      Do While (Winkel<0.d0)
         Winkel=Winkel+2*PI
      End Do
      Do While (Winkel>2*PI)
         Winkel=Winkel-2*PI
      End Do

   End Subroutine Winkelkorr
   !***************************************************

   !***************************************************
   Subroutine Winkelpaarkorr(Winkel1,Winkel2)

      Use Konstanten
      Implicit None
      !------------------------------------------------
      Real(8)  :: Winkel1,Winkel2
      !------------------------------------------------

      Do While (Winkel1<0.d0)
         Winkel1=Winkel1+2*PI
         Winkel2=Winkel2+2*PI
      End Do
      Do While (Winkel1>2*PI)
         Winkel1=Winkel1-2*PI
         Winkel2=Winkel2-2*PI
      End Do

   End Subroutine Winkelpaarkorr
   !***************************************************

   !***************************************************
   Subroutine WZParamInnenverzahnung(z_1,z_0,m_t,alfa_t,xm_n,r_a0,rho_a0,h_f1,s_P0,&
                                     r_w1,r_f1,r_Ff1,h_FaP0,a,x_M,y_M)
      Implicit None
      !------------------------------------------------
      Integer,Intent(in)  :: z_1,z_0
      Real(8),Intent(in)  :: m_t,alfa_t,xm_n,r_a0,rho_a0,h_f1,s_P0
      Real(8),Intent(out) :: r_w1,r_f1,r_Ff1,h_FaP0,a,x_M,y_M
      !------------------------------------------------
      Real(8)  :: AB,DE,AC,DF,EF,EG,AK,AD,BE,BC,EC,CG,HI
      Real(8)  :: cos_alfa_t,sin_alfa_t,tan_alfa_t,theta,phi,phi0
      !------------------------------------------------

      cos_alfa_t=dCos(alfa_t)
      sin_alfa_t=dSin(alfa_t)
      tan_alfa_t=dTan(alfa_t)

      AB=Abs(z_1)*m_t*cos_alfa_t/2
      DE=z_0*m_t*cos_alfa_t/2
      AC=Abs(z_1)*m_t/2
      r_w1=AC

      DF=r_a0-rho_a0
      EF=dSqrt(DF**2-DE**2)
      EG=EF+rho_a0
      AK=AC-xm_n+h_f1
      r_f1=AK

      AD=AK-r_a0
      a=AD

      BE=dSqrt(AD**2-(AB-DE)**2)
      BC=AC*sin_alfa_t
      EC=BC-BE
      CG=EG-EC
      h_FaP0=CG*sin_alfa_t+xm_n
      r_Ff1=dSqrt(AB**2+(BE+EG)**2)

      HI=CG*cos_alfa_t+h_FaP0*tan_alfa_t+s_P0/2

      theta=dACos((AB-DE)/AD)-alfa_t
      phi=dACos(DE/DF)-alfa_t-theta

      phi0=phi-(HI/r_w1-theta)*Abs(z_1)/Dble(z_0)
      x_M=DF*dSin(phi0)
      y_M=DF*dCos(phi0)

   End Subroutine WZParamInnenverzahnung
   !***************************************************

   !***************************************************
   Subroutine Fertigung(WZ_Kontur,WZ,STCK,N_sim,N_Profil,Profil,Profil_id,mode)

      Use Konstanten
      Use Strukturen
      Implicit None
      !------------------------------------------------
      Type(Kontur),Target,Intent(in) :: WZ_Kontur
      Type(CSLage),Intent(in)        :: WZ,STCK
      Integer,Intent(in)             :: N_sim,N_Profil
      Real(8),Intent(inout)          :: Profil(N_Profil,7)
      Integer,Intent(inout)          :: Profil_id(N_Profil,2)
      Integer,Intent(in)             :: mode
      !------------------------------------------------
      Type(Kontur),Pointer :: curSegm
      Type(Bogen),Pointer  :: curBog
      Type(Gerade),Pointer :: curGer
      Real(8) :: thetaMin,thetaMax,RCS,ROC,fiCO,thetaC,root   ,& !
                 myX,myY,myTheta,x1,y1,x2,y2,dx,dy,L,xP,yP,fiP,& !
                 R,Dfi,&
                 d,lambda,lambda1,lambda2,minlim,maxlim,xMin,xMax,yMin,yMax
      Integer :: it,it1,it2,idMin,idMax
      !------------------------------------------------

      Do it=1,N_sim
         Dfi=WZ%fi(it)-STCK%fi(it)
         curSegm=>WZ_Kontur
         Do While (Associated(curSegm))
            If (curSegm%flBog) Then
               curBog=>curSegm%Bogen
               Call WZzuSTCK_CS (it,WZ,STCK,curBog%xWZ_M,curBog%yWZ_M,curBog%xSTCK_M,curBog%ySTCK_M)
               curBog%r_M=dSqrt(curBog%xSTCK_M**2+curBog%ySTCK_M**2)
               curBog%fi_M=dAtan2(curBog%ySTCK_M,curBog%xSTCK_M)
               curBog%fiSTCK_1   =  curBog%fiWZ_1+Dfi
               curBog%fiSTCK_2   =  curBog%fiWZ_2+Dfi
               Call Winkelpaarkorr(curBog%fiSTCK_1,curBog%fiSTCK_2)
            Else
               curGer=>curSegm%Gerade
               Call WZzuSTCK_CS (it,WZ,STCK,curGer%xWZ_1,curGer%yWZ_1,curGer%xSTCK_1,curGer%ySTCK_1)
               Call WZzuSTCK_CS (it,WZ,STCK,curGer%xWZ_2,curGer%yWZ_2,curGer%xSTCK_2,curGer%ySTCK_2)
               curGer%r1=dSqrt(curGer%xSTCK_1**2+curGer%ySTCK_1**2)
               curGer%r2=dSqrt(curGer%xSTCK_2**2+curGer%ySTCK_2**2)
               curGer%fiSTCK_12  =  curGer%fiWZ_12+Dfi
               Call Winkelkorr(curGer%fiSTCK_12)
            End If
            curSegm=>curSegm%next
         End Do
         Do it1=1,N_profil
            R=Profil(it1,1)
            thetaMin=1.d10
            thetaMax=-1.d10
            curSegm=>WZ_Kontur
            Do While (Associated(curSegm))
               If (curSegm%flBog) Then !Kreisbogen
                  curBog=>curSegm%Bogen
                  RCS=curBog%rho
                  ROC=curBog%r_M
                  fiCO=curBog%fi_M+PI
                  thetaC=(RCS**2+ROC**2-R**2)/(2*ROC*RCS)
                  If (dAbs(thetaC)<=1.d0) Then
                     thetaC=dAcos(thetaC)
                     Do it2=-1,1,2
                        root=fiCO+it2*thetaC
                        Dfi=root-curBog%fiSTCK_1
                        Call Winkelkorr(Dfi)
                        If ((Dfi>0d0).AND.(Dfi<curBog%fiSTCK_2-curBog%fiSTCK_1)) Then
                           myX=curBog%xSTCK_M+RCS*dCos(root)
                           myY=curBog%ySTCK_M+RCS*dSin(root)
                           myTheta=dAtan2(myY,myX)
                           Call Winkelkorr(myTheta)
                           If (myTheta<thetaMin) Then
                              idMin=curSegm%id
                              thetaMin=myTheta
                              xMin=myX
                              yMin=myY
                           End If
                           If (myTheta>thetaMax) Then
                              idMax=curSegm%id
                              thetaMax=myTheta
                              xMax=myX
                              yMax=myY
                           End If
                        End If
                     End Do
                  End If
               Else   !Gerade
                  curGer=>curSegm%Gerade
                  x1=curGer%xSTCK_1
                  y1=curGer%ySTCK_1
                  x2=curGer%xSTCK_2
                  y2=curGer%ySTCK_2
                  dx=x2-x1
                  dy=y2-y1
                  L=dSqrt(dx**2+dy**2)

                  yP=(x1*dy-y1*dx)/L**2
                  xP=yP*dy
                  yP=-yP*dx
                  d=dSqrt(xP**2+yP**2)
                  fiP=dAtan2(yP,xP)

                  If (dAbs(dx)>dAbs(dy)) Then
                     lambda1=(x1-xP)/dx
                     lambda2=(x2-xP)/dx
                  Else
                     lambda1=(y1-yP)/dy
                     lambda2=(y2-yP)/dy
                  End If
                  lambda1=lambda1*L
                  lambda2=lambda2*L

                  minlim=dMin1(lambda1,lambda2)
                  maxlim=dMax1(lambda1,lambda2)
                 If (R>=d-1d-12) Then
                     lambda=dMax1(0d0,R**2-d**2) 
                     lambda=dSqrt(lambda)
                     Do it2=-1,1,2
                        root=it2*lambda
                        If ((root<maxlim).AND.(root>minlim)) Then
                           myX=xP+root*dx/L
                           myY=yP+root*dy/L
                           myTheta=dAtan2(myY,myX)
                           Call Winkelkorr(myTheta)
                           If (myTheta<thetaMin) Then
                              idMin=curSegm%id
                              thetaMin=myTheta
                              xMin=myX
                              yMin=myY
                           End If
                           If (myTheta>thetaMax) Then
                              idMax=curSegm%id
                              thetaMax=myTheta
                              xMax=myX
                              yMax=myY
                           End If
                        End If
                     End Do
                  End If
               End If
               curSegm=>curSegm%next
            End Do
            If ((mode==0).AND.(thetaMin<=thetaMax)) Then
               If (thetaMin<Profil(it1,2)) Then
                  Profil_id(it1,1)=idMin
                  Profil(it1,2)=thetaMin
                  Profil(it1,4)=xMin
                  Profil(it1,5)=yMin
               End If
               If (thetaMax>Profil(it1,3)) Then
                  Profil_id(it1,2)=idMax
                  Profil(it1,3)=thetaMax
                  Profil(it1,6)=xMax
                  Profil(it1,7)=yMax
               End If
            Else If ((mode==-1).AND.(thetaMin<=thetaMax)) Then
               If (thetaMax>Profil(it1,2)) Then
                  Profil_id(it1,1)=idMax
                  Profil(it1,2)=thetaMax
                  Profil(it1,4)=xMax
                  Profil(it1,5)=yMax
               End If
            Else If ((mode==1).AND.(thetaMin<=thetaMax)) Then
               If (thetaMin<Profil(it1,3)) Then
                  Profil_id(it1,2)=idMin
                  Profil(it1,3)=thetaMin
                  Profil(it1,6)=xMin
                  Profil(it1,7)=yMin
               End If
            End If
         End Do   !it1=1,N_profil
      End Do   !it=1,N_sim

   End Subroutine Fertigung
   !***************************************************

End Module Prozeduren
!***************************************************


!***************************************************
Program Verzahnen

   Use Konstanten
   Use Strukturen
   Use Prozeduren
   Implicit None

   !------ Variablen -------------------------------------
   Integer  :: uid_in,uid_kontur,uid_zahnstange
   Integer  :: N_sim,N_profil,N_start,N_kontur,N_fussgrund,N_zahnkopf
   Integer  :: z_1,it,it1,wz_id

   Real(8)  :: alfa_P0  ,  alfa_KP0 ,& !
               alfa_prP0,  m_n      ,& !
               h_fP0    ,  h_aP0    ,& !
               rho_aP0  ,  h_FaP0   ,& !
               h_prP0   ,  h_FprP0  ,& !
               s_P0     ,  x_1      ,& !
               beta     ,  m_t      ,& !
               d_a,r_a  ,  r_w      ,& !
               R,r_f,dTmp,myTheta,dx,dy
   !Benoetigt nur fuer Innenverzahnungen
   Integer  :: z_0
   Real(8)  :: d_a0,r_a0,a_schneidrad,x_M,y_M,r_Ff

   Real(8),Allocatable::   Profil(:,:),KonturXY(:,:)

   Integer,Allocatable::   Profil_id(:,:),Kontur_id(:)

   Character(200) :: arg1,arg2,input_file,output_file,output_file2
   Character(20)  :: cTmp

   Integer:: ioerror
   Real(8):: X,Y,X0,Y0

   Type(CSLage)   :: WZ,STCK

   Type(Bogen),Pointer  :: curBog
   Type(Gerade),Pointer :: curGer

   Type(Kontur),Target  :: myKontur,myKontur1,myKontur2
   Type(Kontur),Pointer :: curSegm,newSegm
   !------------------------------------------------------

   !Initialisierungen
   N_sim=1000
   N_profil=1200
   N_zahnkopf=7
   input_file='Verzahnen.txt'
   output_file='Kontur.dat'
   output_file2='Zahnstange.dat'

   !Kommandozeileparameter
   it = IARGC()
   Do While (it>0)
      Call GETARG(it, arg1)
      If (arg1(1:1)=='-') Then
         Call GETARG(it+1, arg2)
         Select Case (arg1(2:2))
            Case ('i')
               input_file = arg2
            Case ('o')
               output_file = arg2
               output_file2 = Trim(arg2)//'.1'
            Case ('N')
               If (arg1(2:6) == 'N_sim') Then
                  Read (arg2,*) N_sim
               Else If (arg1(2:9) == 'N_profil') Then
                  Read (arg2,*) N_profil
               Else If (arg1(2:11) == 'N_zahnkopf') Then
                  Read (arg2,*) N_zahnkopf
               End If
         End Select
      End If
      it = it - 1
   End Do

   !Matrizeninitialiserung
   Allocate(Profil(N_profil,7))     ; Profil=0.d0
   Allocate(Profil_id(N_profil,2))  ; Profil_id=0

   !------------------------------------------------------
   Call get_io_unit(uid_in)
   Open (uid_in,File=Adjustl(Trim(input_file)))
   Read (uid_in,*) wz_id

   !------ STCK-Data einlesen ----------------------------
   Read (uid_in,*) cTmp,z_1
   Read (uid_in,*) cTmp,d_a        ;  d_a=dAbs(d_a)  ;  r_a=d_a/2.d0
   Read (uid_in,*) cTmp,x_1
   Read (uid_in,*) cTmp,beta       ;  beta=beta*PI/180.d0

   If (wz_id==0) Then !Ohne Protuberanz
      !------ WZ-Data einlesen ------------------------------
      Read (uid_in,*) cTmp,alfa_P0    ;  alfa_P0=alfa_P0*PI/180.d0
      Read (uid_in,*) cTmp,m_n
      Read (uid_in,*) cTmp,h_fP0      ;  h_fP0    =h_fP0*m_n
      Read (uid_in,*) cTmp,h_aP0      ;  h_aP0    =h_aP0*m_n
      Read (uid_in,*) cTmp,rho_aP0    ;  rho_aP0  =rho_aP0*m_n

      !------ Schraegungswinkel beruecksichtigen ------------
      m_t = m_n/dCos(beta)
      r_w=Abs(z_1)*m_t/2.d0
      alfa_P0 = dATan(dTan(alfa_P0)/dCos(beta))
      !der Einfluss des Schraegungswinkels auf rho_aP0 wird vernachlaessigt

      !------ WZ-Hilfsvariablen berechnen -------------------
      s_P0=PI*m_t/2
      h_FaP0=h_aP0-(1.d0-dSin(alfa_P0))*rho_aP0

      !------ WZ-Struktur aufbauen --------------------------
      curSegm=>myKontur
         !1. Teil (Gerade)
         curSegm%id=1
         curSegm%flBog=.FALSE.
         curSegm%Bogen=>NULL()
         Allocate(curGer)
         curSegm%Gerade=>curGer
            curGer%xWZ_1   =  -s_P0/2.d0-h_fP0*dTan(alfa_P0)
            curGer%yWZ_1   =  2*h_fP0
            curGer%xWZ_2   =  -s_P0/2.d0-h_fP0*dTan(alfa_P0)
            curGer%yWZ_2   =  h_fP0
            curGer%fiWZ_12 =  dAtan2(curGer%yWZ_2-curGer%yWZ_1,curGer%xWZ_2-curGer%xWZ_1)
      Allocate(newSegm)
      curSegm%next=>newSegm
      curSegm=>curSegm%next
         !2. Teil (Gerade)
         curSegm%id=2
         curSegm%flBog=.FALSE.
         curSegm%Bogen=>NULL()
         Allocate(curGer)
         curSegm%Gerade=>curGer
            curGer%xWZ_1   =  -s_P0/2.d0-h_fP0*dTan(alfa_P0)
            curGer%yWZ_1   =  h_fP0
            curGer%xWZ_2   =  -s_P0/2.d0+h_FaP0*dTan(alfa_P0)
            curGer%yWZ_2   =  -h_FaP0
            curGer%fiWZ_12 =  dAtan2(curGer%yWZ_2-curGer%yWZ_1,curGer%xWZ_2-curGer%xWZ_1)
      Allocate(newSegm)
      curSegm%next=>newSegm
      curSegm=>curSegm%next
         !3. Teil (Bogen)
         curSegm%id=3
         curSegm%flBog=.TRUE.
         Allocate(curBog)
         curSegm%Bogen=>curBog
            curBog%xWZ_M   =  -s_P0/2.d0+h_FaP0*dTan(alfa_P0)+rho_aP0*dCos(alfa_P0)
            curBog%yWZ_M   =  -h_aP0+rho_aP0
            curBog%rho     =  rho_aP0
            curBog%fiWZ_1  =  PI+alfa_P0
            curBog%fiWZ_2  =  3.d0*PI/2.d0
         curSegm%Gerade=>NULL()
      Allocate(newSegm)
      curSegm%next=>newSegm
      curSegm=>curSegm%next
         !4. Teil (Gerade)
         curSegm%id=4
         curSegm%flBog=.FALSE.
         curSegm%Bogen=>NULL()
         Allocate(curGer)
         curSegm%Gerade=>curGer
            curGer%xWZ_1   =  -s_P0/2.d0+h_FaP0*dTan(alfa_P0)+rho_aP0*dCos(alfa_P0)
            curGer%yWZ_1   =  -h_aP0
            curGer%xWZ_2   =  -curGer%xWZ_1
            curGer%yWZ_2   =  curGer%yWZ_1
            curGer%fiWZ_12 =  0.d0
      Allocate(newSegm)
      curSegm%next=>newSegm
      curSegm=>curSegm%next
         !5. Teil (Bogen)
         curSegm%id=5
         curSegm%flBog=.TRUE.
         Allocate(curBog)
         curSegm%Bogen=>curBog
            curBog%xWZ_M   =  s_P0/2.d0-h_FaP0*dTan(alfa_P0)-rho_aP0*dCos(alfa_P0)
            curBog%yWZ_M   =  -h_aP0+rho_aP0
            curBog%rho     =  rho_aP0
            curBog%fiWZ_1  =  3.d0*PI/2.d0
            curBog%fiWZ_2  =  2.d0*PI-alfa_P0
         curSegm%Gerade=>NULL()
      Allocate(newSegm)
      curSegm%next=>newSegm
      curSegm=>curSegm%next
         !6. Teil (Gerade)
         curSegm%id=6
         curSegm%flBog=.FALSE.
         curSegm%Bogen=>NULL()
         Allocate(curGer)
         curSegm%Gerade=>curGer
            curGer%xWZ_1   =  s_P0/2.d0-h_FaP0*dTan(alfa_P0)
            curGer%yWZ_1   =  -h_FaP0
            curGer%xWZ_2   =  s_P0/2.d0+h_fP0*dTan(alfa_P0)
            curGer%yWZ_2   =  h_fP0
            curGer%fiWZ_12 =  dAtan2(curGer%yWZ_2-curGer%yWZ_1,curGer%xWZ_2-curGer%xWZ_1)
      Allocate(newSegm)
      curSegm%next=>newSegm
      curSegm=>curSegm%next
         !7. Teil (Gerade)
         curSegm%id=7
         curSegm%flBog=.FALSE.
         curSegm%Bogen=>NULL()
         Allocate(curGer)
         curSegm%Gerade=>curGer
            curGer%xWZ_1   =  s_P0/2.d0+h_fP0*dTan(alfa_P0)
            curGer%yWZ_1   =  h_fP0
            curGer%xWZ_2   =  s_P0/2.d0+h_fP0*dTan(alfa_P0)
            curGer%yWZ_2   =  2*h_fP0
            curGer%fiWZ_12 =  dAtan2(curGer%yWZ_2-curGer%yWZ_1,curGer%xWZ_2-curGer%xWZ_1)
      curSegm%next=>NULL()
      newSegm=>NULL()
   Else If (wz_id==1) Then !Mit Protuberanz
      !------ WZ-Data einlesen ------------------------------
      Read (uid_in,*) cTmp,alfa_P0    ;  alfa_P0=alfa_P0*PI/180.d0
      Read (uid_in,*) cTmp,alfa_KP0   ;  alfa_KP0=alfa_KP0*PI/180.d0
      Read (uid_in,*) cTmp,alfa_prP0  ;  alfa_prP0=alfa_prP0*PI/180.d0
      Read (uid_in,*) cTmp,m_n
      Read (uid_in,*) cTmp,h_fP0      ;  h_fP0    =h_fP0*m_n
      Read (uid_in,*) cTmp,h_FaP0     ;  h_FaP0   =h_FaP0*m_n
      Read (uid_in,*) cTmp,h_prP0     ;  h_prP0   =h_prP0*m_n
      Read (uid_in,*) cTmp,rho_aP0    ;  rho_aP0  =rho_aP0*m_n

      !------ Schraegungswinkel beruecksichtigen ------------
      m_t = m_n/dCos(beta)
      r_w=Abs(z_1)*m_t/2.d0
      alfa_P0 = dATan(dTan(alfa_P0)/dCos(beta))
      alfa_KP0 = dATan(dTan(alfa_KP0)/dCos(beta))
      alfa_prP0 = dATan(dTan(alfa_prP0)/dCos(beta))
      !der Einfluss der Schraegungswinkels auf rho_aP0 wird vernachlässigt

      !------ WZ-Hilfsvariablen berechnen -------------------
      s_P0=PI*m_t/2
      h_FprP0=h_prP0-(1.d0-dSin(alfa_prP0))*rho_aP0
      dTmp=-s_P0/2.d0+h_FaP0*dTan(alfa_P0)
      h_aP0 = h_prP0

      !------ WZ-Struktur aufbauen --------------------------
      curSegm=>myKontur
         !1. Teil (Gerade)
         curSegm%id=1
         curSegm%flBog=.FALSE.
         curSegm%Bogen=>NULL()
         Allocate(curGer)
         curSegm%Gerade=>curGer
            curGer%xWZ_1   =  -s_P0/2.d0-h_fP0*dTan(alfa_P0)
            curGer%yWZ_1   =  2*h_fP0
            curGer%xWZ_2   =  curGer%xWZ_1
            curGer%yWZ_2   =  h_fP0
            curGer%fiWZ_12 =  dAtan2(curGer%yWZ_2-curGer%yWZ_1,curGer%xWZ_2-curGer%xWZ_1)
      Allocate(newSegm)
      curSegm%next=>newSegm
      curSegm=>curSegm%next
         !2. Teil (Gerade)
         curSegm%id=2
         curSegm%flBog=.FALSE.
         curSegm%Bogen=>NULL()
         Allocate(curGer)
         curSegm%Gerade=>curGer
            curGer%xWZ_1   =  -s_P0/2.d0-h_fP0*dTan(alfa_P0)
            curGer%yWZ_1   =  h_fP0
            curGer%xWZ_2   =  -s_P0/2.d0+h_FaP0*dTan(alfa_P0)
            curGer%yWZ_2   =  -h_FaP0
            curGer%fiWZ_12 =  dAtan2(curGer%yWZ_2-curGer%yWZ_1,curGer%xWZ_2-curGer%xWZ_1)
      Allocate(newSegm)
      curSegm%next=>newSegm
      curSegm=>curSegm%next
         !3. Teil (Gerade)
         curSegm%id=3
         curSegm%flBog=.FALSE.
         curSegm%Bogen=>NULL()
         Allocate(curGer)
         curSegm%Gerade=>curGer
            curGer%xWZ_1   =  dTmp
            curGer%yWZ_1   =  -h_FaP0
            curGer%xWZ_2   =  dTmp+h_FprP0*dTan(alfa_prP0)
            curGer%yWZ_2   =  -h_FaP0-h_FprP0
            curGer%fiWZ_12 =  dAtan2(curGer%yWZ_2-curGer%yWZ_1,curGer%xWZ_2-curGer%xWZ_1)
      Allocate(newSegm)
      curSegm%next=>newSegm
      curSegm=>curSegm%next
         !4. Teil (Bogen)
         curSegm%id=4
         curSegm%flBog=.TRUE.
         Allocate(curBog)
         curSegm%Bogen=>curBog
            curBog%xWZ_M   =  dTmp+h_FprP0*dTan(alfa_prP0)+rho_aP0*dCos(alfa_prP0)
            curBog%yWZ_M   =  -h_FaP0-h_prP0+rho_aP0
            curBog%rho     =  rho_aP0
            curBog%fiWZ_1  =  PI+alfa_prP0
            curBog%fiWZ_2  =  3.d0*PI/2.d0
         curSegm%Gerade=>NULL()
      Allocate(newSegm)
      curSegm%next=>newSegm
      curSegm=>curSegm%next
         !5. Teil (Gerade)
         curSegm%id=5
         curSegm%flBog=.FALSE.
         curSegm%Bogen=>NULL()
         Allocate(curGer)
         curSegm%Gerade=>curGer
            curGer%xWZ_1   =  dTmp+h_FprP0*dTan(alfa_prP0)+rho_aP0*dCos(alfa_prP0)
            curGer%yWZ_1   =  -h_FaP0-h_prP0
            curGer%xWZ_2   =  -curGer%xWZ_1
            curGer%yWZ_2   =  curGer%yWZ_1
            curGer%fiWZ_12 =  0.d0
      Allocate(newSegm)
      curSegm%next=>newSegm
      curSegm=>curSegm%next
         !6. Teil (Bogen)
         curSegm%id=6
         curSegm%flBog=.TRUE.
         Allocate(curBog)
         curSegm%Bogen=>curBog
            curBog%xWZ_M   =  -dTmp-h_FprP0*dTan(alfa_prP0)-rho_aP0*dCos(alfa_prP0)
            curBog%yWZ_M   =  -h_FaP0-h_prP0+rho_aP0
            curBog%rho     =  rho_aP0
            curBog%fiWZ_1  =  3.d0*PI/2.d0
            curBog%fiWZ_2  =  2*PI-alfa_prP0
         curSegm%Gerade=>NULL()
      Allocate(newSegm)
      curSegm%next=>newSegm
      curSegm=>curSegm%next
         !7. Teil (Gerade)
         curSegm%id=7
         curSegm%flBog=.FALSE.
         curSegm%Bogen=>NULL()
         Allocate(curGer)
         curSegm%Gerade=>curGer
            curGer%xWZ_1   =  -dTmp
            curGer%yWZ_1   =  -h_FaP0
            curGer%xWZ_2   =  -dTmp-h_FprP0*dTan(alfa_prP0)
            curGer%yWZ_2   =  -h_FaP0-h_FprP0
            curGer%fiWZ_12 =  dAtan2(curGer%yWZ_2-curGer%yWZ_1,curGer%xWZ_2-curGer%xWZ_1)
      Allocate(newSegm)
      curSegm%next=>newSegm
      curSegm=>curSegm%next
         !8. Teil (Gerade)
         curSegm%id=8
         curSegm%flBog=.FALSE.
         curSegm%Bogen=>NULL()
         Allocate(curGer)
         curSegm%Gerade=>curGer
            curGer%xWZ_1   =  s_P0/2.d0+h_fP0*dTan(alfa_P0)
            curGer%yWZ_1   =  h_fP0
            curGer%xWZ_2   =  s_P0/2.d0-h_FaP0*dTan(alfa_P0)
            curGer%yWZ_2   =  -h_FaP0
            curGer%fiWZ_12 =  dAtan2(curGer%yWZ_2-curGer%yWZ_1,curGer%xWZ_2-curGer%xWZ_1)
      Allocate(newSegm)
      curSegm%next=>newSegm
      curSegm=>curSegm%next
         !9. Teil (Gerade)
         curSegm%id=9
         curSegm%flBog=.FALSE.
         curSegm%Bogen=>NULL()
         Allocate(curGer)
         curSegm%Gerade=>curGer
            curGer%xWZ_1   =  s_P0/2.d0+h_fP0*dTan(alfa_P0)
            curGer%yWZ_1   =  2*h_fP0
            curGer%xWZ_2   =  curGer%xWZ_1
            curGer%yWZ_2   =  h_fP0
            curGer%fiWZ_12 =  dAtan2(curGer%yWZ_2-curGer%yWZ_1,curGer%xWZ_2-curGer%xWZ_1)
      curSegm%next=>NULL()

      newSegm=>NULL()
   Else If ((wz_id==2).AND.(z_1<0.d0)) Then !Innenverzahnung
      !------ WZ-Data einlesen ------------------------------
      Read (uid_in,*) cTmp,alfa_P0    ;  alfa_P0=alfa_P0*PI/180.d0
      Read (uid_in,*) cTmp,m_n
      Read (uid_in,*) cTmp,h_aP0      ;  h_aP0    =h_aP0*m_n
      Read (uid_in,*) cTmp,rho_aP0    ;  rho_aP0  =rho_aP0*m_n
      Read (uid_in,*) cTmp,d_a0       ;  r_a0     =d_a0/2     !Scneidradaussendurchmesser bzw. -radius
      Read (uid_in,*) cTmp,z_0                                !Scneidradzaehnezahl

      !------ Schraegungswinkel beruecksichtigen ------------
      m_t = m_n/dCos(beta)
      alfa_P0 = dATan(dTan(alfa_P0)/dCos(beta))
      !der Einfluss des Schraegungswinkels auf rho_aP0 wird vernachlaessigt

      !------ WZ-Hilfsvariablen berechnen -------------------
      s_P0=PI*m_t/2
      Call WZParamInnenverzahnung(z_1,z_0,m_t,alfa_P0,x_1*m_n,r_a0,rho_aP0,h_aP0,s_P0,&
                                  r_w,r_f,r_Ff,h_FaP0,a_schneidrad,x_M,y_M)
      h_fP0=s_P0/(2*dTan(alfa_P0))

      !------ WZ-Struktur aufbauen --------------------------
      curSegm=>myKontur
         !1. Teil (Gerade)
         curSegm%id=1
         curSegm%flBog=.FALSE.
         curSegm%Bogen=>NULL()
         Allocate(curGer)
         curSegm%Gerade=>curGer
            curGer%xWZ_1   =  -s_P0/2.d0-2*h_aP0*dTan(alfa_P0)
            curGer%yWZ_1   =  2*h_aP0
            curGer%xWZ_2   =  0.d0
            curGer%yWZ_2   =  -h_fP0
            curGer%fiWZ_12 =  dAtan2(curGer%yWZ_2-curGer%yWZ_1,curGer%xWZ_2-curGer%xWZ_1)
      Allocate(newSegm)
      curSegm%next=>newSegm
      curSegm=>curSegm%next
         !2. Teil (Gerade)
         curSegm%id=2
         curSegm%flBog=.FALSE.
         curSegm%Bogen=>NULL()
         Allocate(curGer)
         curSegm%Gerade=>curGer
            curGer%xWZ_1   =  0.d0
            curGer%yWZ_1   =  -h_fP0
            curGer%xWZ_2   =  s_P0/2.d0+2*h_aP0*dTan(alfa_P0)
            curGer%yWZ_2   =  2*h_aP0
            curGer%fiWZ_12 =  dAtan2(curGer%yWZ_2-curGer%yWZ_1,curGer%xWZ_2-curGer%xWZ_1)
      curSegm%next=>NULL()

      curSegm=>myKontur1
         !3. Teil (Bogen)
         curSegm%id=3
         curSegm%flBog=.TRUE.
         Allocate(curBog)
         curSegm%Bogen=>curBog
            curBog%xWZ_M   =  x_M
            curBog%yWZ_M   =  y_M
            curBog%rho     =  rho_aP0
            curBog%fiWZ_1  =  0.d0
            curBog%fiWZ_2  =  PI/2 !dAtan2(y_M,x_M)
         curSegm%Gerade=>NULL()
      curSegm%next=>NULL()
      newSegm=>NULL()

      curSegm=>myKontur2
         !4. Teil (Bogen)
         curSegm%id=4
         curSegm%flBog=.TRUE.
         Allocate(curBog)
         curSegm%Bogen=>curBog
            curBog%xWZ_M   =  -x_M
            curBog%yWZ_M   =  y_M
            curBog%rho     =  rho_aP0
            curBog%fiWZ_1  =  PI/2 !dAtan2(y_M,-x_M)
            curBog%fiWZ_2  =  PI
         curSegm%Gerade=>NULL()
      curSegm%next=>NULL()
      newSegm=>NULL()
   Else If (wz_id==3) Then !Polygon
      !------ WZ-Data einlesen ------------------------------
      Read (uid_in,*) cTmp,m_n

      !------ Schraegungswinkel beruecksichtigen ------------
      m_t = m_n/dCos(beta)
      r_w=Abs(z_1)*m_t/2.d0

      !------ WZ-Struktur einlesen und aufbauen -------------
      curSegm=>myKontur
      Read (uid_in,*) X0,Y0
      Do it=1,100 ! Nicht mehr als 100 Punkte
         Read (uid_in,*,iostat=ioerror) X,Y
         If (ioerror /= 0) Exit
         If (it>1) Then
           Allocate(newSegm)
           curSegm%next=>newSegm
           curSegm=>curSegm%next
         End If
         curSegm%id=it
         curSegm%flBog=.FALSE.
         curSegm%Bogen=>NULL()
         Allocate(curGer)
         curSegm%Gerade=>curGer
            curGer%xWZ_1   =  X0
            curGer%yWZ_1   =  Y0/dCos(beta)
            curGer%xWZ_2   =  X
            curGer%yWZ_2   =  Y/dCos(beta)
            curGer%fiWZ_12 =  dAtan2(curGer%yWZ_2-curGer%yWZ_1,curGer%xWZ_2-curGer%xWZ_1)
         X0 = X
         Y0 = Y
         If (h_aP0 < -Y) h_aP0 = -Y
      End Do
      curSegm%next=>NULL()
      newSegm=>NULL()
   End If
   Close(uid_in)

   !------ Bahn des Werkstückes berechnen ----------------
   Allocate(STCK%X(N_sim))  ;  STCK%X=0.d0
   Allocate(STCK%Y(N_sim))  ;  STCK%Y=0.d0
   Allocate(STCK%fi(N_sim)) ;  STCK%fi=0.d0

   dTmp=PI/(2.d0*(N_sim-1))
   STCK%fi=(/(it*dTmp,it=0,N_sim-1)/)-PI/4.d0

   !------ Bahn des Werkzeuges berechnen -----------------
   Allocate(WZ%X(N_sim))    ;  WZ%X=0.d0
   Allocate(WZ%Y(N_sim))    ;  WZ%Y=0.d0
   Allocate(WZ%fi(N_sim))   ;  WZ%fi=0.d0

   WZ%X=-r_w*STCK%fi
   WZ%Y=r_w+Sign(1,z_1)*x_1*m_n

   !------ Profil berechnen ------------------------------
   r_f=r_w+Sign(1,z_1)*(x_1*m_n-h_aP0)
   dTmp=(r_a-r_f)/(N_profil-1)
   Profil(:,1)=(/(r_f+it1*dTmp,it1=0,N_profil-1)/)
   Profil(:,2)=1.d6*m_t
   Profil(:,3)=-1.d6*m_t

   Call Fertigung(myKontur,WZ,STCK,N_sim,N_Profil,Profil,Profil_id,0)

   If (wz_id==2) Then
      !------ Neue bahn des Werkzeuges berechnen ---------
      WZ%X=0.d0
      WZ%Y=a_schneidrad
      WZ%fi=(Abs(z_1)*STCK%fi)/z_0
      it=1
      Do While (Profil(it,1)>r_Ff)
         Profil(it,2)=-1.d6*m_t
         Profil(it,3)=1.d6*m_t
         Profil_id(it,1)=-1
         Profil_id(it,2)=-1
         it=it+1
      End Do
      Call Fertigung(myKontur1,WZ,STCK,N_sim,N_Profil,Profil,Profil_id,1)
      Call Fertigung(myKontur2,WZ,STCK,N_sim,N_Profil,Profil,Profil_id,-1)
   End If
   Deallocate(STCK%X)
   Deallocate(STCK%Y)
   Deallocate(STCK%fi)
   Deallocate(WZ%X)
   Deallocate(WZ%Y)
   Deallocate(WZ%fi)

   !Bestimmung des ersten Punktes des Profils (entspricht dem Zahnfussgrund)
   N_start=1
   If (z_1>0) Then
      Do While (Profil_id(N_start,2)==-1)
         N_start=N_start+1
      End Do
   Else
      Do While (Profil(N_start,2)<PI/2-PI/Abs(z_1).OR.(Profil_id(N_start,1)==-1))
         N_start=N_start+1
      End Do
   End If

   !Ergaenzung des Zahnfussgrundes mit einem Kreisbogen
   If (z_1>0) Then
      N_fussgrund=Ceiling((Profil(N_start  ,3)-PI/2            )/&
                          (Profil(N_start+1,3)-Profil(N_start,3)))
   Else
      N_fussgrund=Ceiling((Profil(N_start  ,2)-PI/2+PI/Abs(z_1))/&
                          (Profil(N_start+1,2)-Profil(N_start,2)))
   End If
   Allocate(KonturXY(N_fussgrund+N_zahnkopf+N_profil-N_start+1,2))
   KonturXY=0.d0
   Allocate(Kontur_id(N_fussgrund+N_zahnkopf+N_profil-N_start+1))
   Kontur_id=0

   N_Kontur=0
   If (N_fussgrund>0) Then
      If (z_1>0) Then
         dTmp=(Profil(N_start,3)-PI/2)/N_fussgrund
      Else
         dTmp=(Profil(N_start,2)-PI/2+PI/Abs(z_1))/N_fussgrund
      End If
      Do it=1,N_fussgrund
         Kontur_id(it)=-1
         myTheta=PI/2-PI/Abs(z_1)+dTmp*(it-1)
         KonturXY(it,1)=r_f*dCos(myTheta)
         KonturXY(it,2)=r_f*dSin(myTheta)
      End Do
      N_Kontur=N_fussgrund
   End If

   !Uebertragen des Vektors Profil in den Vektor Kontur
   !(mit Verdehung um eine halbe Teilung bei Aussenverzahnungen)
   If (z_1>0) Then
      dx=dSin(PI/z_1)
      dy=dCos(PI/z_1)
      Do it=N_start,N_profil
         If (Profil(it,2)<Profil(it,3)) Then
            N_kontur=N_kontur+1
            Kontur_id(N_kontur)=Profil_id(it,2)
            KonturXY(N_kontur,1)= Profil(it,6)*dy+Profil(it,7)*dx
            KonturXY(N_kontur,2)=-Profil(it,6)*dx+Profil(it,7)*dy
         Else
            Write(*,'("Fehler in der erechneten Kontur bei r=",F10.5)') Profil(it,1)
         End If
      End Do
   Else
      Do it=N_start,N_profil
         If (Profil(it,2)<Profil(it,3)) Then
            N_kontur=N_kontur+1
            Kontur_id(N_kontur)=Profil_id(it,1)
            KonturXY(N_kontur,1)= Profil(it,4)
            KonturXY(N_kontur,2)= Profil(it,5)
         Else
            Write(*,'("Fehler in der erechneten Kontur bei r=",F10.5)') Profil(it,1)
         End If
      End Do
   End If

   !Ergaenzung des Zahnkopfes mit einem Kreisbogen
   R=Profil(N_profil,1)
   If (z_1>0) Then
      myTheta=Profil(N_profil,3)-PI/z_1
   Else
      myTheta=Profil(N_profil,2)
   End If
   dTmp=(PI/2-myTheta)/N_zahnkopf
   Do it=1,N_zahnkopf
      N_kontur=N_kontur+1
      Kontur_id(N_kontur)=-2
      KonturXY(N_kontur,1)=R*dCos(myTheta+dTmp*it)
      KonturXY(N_kontur,2)=R*dSin(myTheta+dTmp*it)
   End Do
   Deallocate(Profil)
   Deallocate(Profil_id)

   Call get_io_unit(uid_kontur)
   Open(uid_kontur,File=Adjustl(Trim(output_file)))
   Write(uid_kontur,'("#",T15,"x",T31,"y",T34,"id")')
   Do it=1,N_kontur
      Write(uid_kontur,'(F15.10,",",F15.10,",",I3)') KonturXY(it,1:2),Kontur_id(it)
   End Do
   Close(uid_kontur)
   Deallocate(KonturXY)
   Deallocate(Kontur_id)

   Call get_io_unit(uid_zahnstange)
   Open(uid_zahnstange,File=Adjustl(Trim(output_file2)))
      curSegm=>myKontur
      Do While (Associated(curSegm))
         If (curSegm%flBog) Then
            curBog=>curSegm%Bogen
            Write (uid_zahnstange,*) 'Bogen'
            Write (uid_zahnstange,'(3(A4,F12.5,3X))') 'x_M=',curBog%xWZ_M  ,&
                                                      'y_M=',curBog%yWZ_M  ,&
                                                      'rho=',curBog%rho
            Deallocate(curBog)
         Else
            curGer=>curSegm%Gerade
            Write (uid_zahnstange,*) 'Gerade'
            Write (uid_zahnstange,'(4(A4,F12.5,3X))') 'x_A=',curGer%xWZ_1  ,&
                                                      'y_A=',curGer%yWZ_1  ,&
                                                      'x_E=',curGer%xWZ_2  ,&
                                                      'y_E=',curGer%yWZ_2
            Deallocate(curGer)
         End If
         curSegm=>curSegm%next
      End Do
      newSegm=>myKontur%next
      Do While (Associated(newSegm))
         curSegm=>newSegm
         newSegm=>newSegm%next
         Deallocate(curSegm)
      End Do
   Close(uid_zahnstange)

End Program Verzahnen
!***************************************************
