module soilhydrologymod

use arveparams, only : dp

implicit none

public :: soilwaterflux
public :: newsnowflux
public :: snowdynamics

private :: comb  !snow layer combination utility calculator

contains

!-------------------
subroutine soilwaterflux(j)

use arveparams, only : lb,Pmax,grav,daysec,kd,Timp,pliq,pice,Lf,Tfreeze,Smin,&
                        ns,nl,npft,peatlim,alpha,aq_thick,psimax
use statevars, only : sv,doy,dtime
use soilstate_vars, only : surf
use tridiagonal
use pft_state_variables, only : veg
use pftparametersmod,    only : prm
use iovariables, only         : soil

use metvarsmod,                 only : met,dm  !only for testing

implicit none

!arguments
integer, intent(in) :: j                        !index value of the current grid cell

!pointers
integer, pointer :: snl                         !index value of top snow layer (negative is more layers)
integer, pointer :: jint                        !index of the soil layer directly above the water table
integer, pointer  :: gnl                        !index value of lowest 'soil' layer
logical, pointer, dimension(:) :: present       !true if pft is present
real(dp), pointer, dimension(:)   :: Aliq       !daily storage of soil liquid water content (fraction)
real(dp), pointer, dimension(:,:) :: rootfracl  !root fraction in that soil layer
real(dp), pointer, dimension(:) :: dz           !thickness of the soil layers (m)
real(dp), pointer, dimension(:) :: zpos         !z coordinate of the middle of the soil layer relative to soil surface (m), positive downwards
real(dp), pointer, dimension(:) :: dzmm         !thickness of the soil layers (mm)
real(dp), pointer, dimension(:) :: zipos        !z coordinate of the middle of the soil layer relative to soil surface (m), positive downwards
real(dp), pointer, dimension(:) :: zposmm       !z coordinate of the middle of the soil layer relative to soil surface (mm), positive downwards
real(dp), pointer, dimension(:) :: sorg         !soil organic matter content (percent)
real(dp), pointer, dimension(:) :: Ksat         !soil water saturated conductivity at layer midpoint (mm s-1)
real(dp), pointer, dimension(:) :: Tsat         !soil water volumetric water content at saturation (fraction)
real(dp), pointer, dimension(:) :: Bexp         !soil water B exponent used in the Brooks & Corey Pedotransfer Functions
real(dp), pointer, dimension(:) :: Psat         !soil water matric potential at saturation (mm)
real(dp), pointer, dimension(:) :: Wliq         !soil liquid water content at layer midpoint (mm)
real(dp), pointer, dimension(:) :: Wice         !soil ice content at layer midpoint (mm)
real(dp), pointer, dimension(:) :: Psi          !soil water potential at layer midpoint (mm)
real(dp), pointer, dimension(:) :: Tice         !soil ice content (fraction)
real(dp), pointer, dimension(:) :: Tliq         !timestep soil liquid water content (fraction)
real(dp), pointer, dimension(:) :: Tpor         !soil volumetric porosity (liquid minus ice content (or stones?)) (fraction)
real(dp), pointer, dimension(:) :: daet         !daily actual evapotranspiration (mm)
real(dp), pointer, dimension(:) :: Ffrz         !fractional impermeable area as a function of soil ice content at a layer
real(dp), pointer, dimension(:) :: betaT        !soil moisture function limiting transpiration (fraction)
real(dp), pointer, dimension(:) :: Psi_eq       ! restriction for min of soil potential (mm)
real(dp), pointer, dimension(:) :: Ku           !soil water instantaneous (unsaturated) conductivity across layer boundary (mm s-1)
real(dp), pointer, dimension(:) :: rPsi         !Psi calc used in root extraction calculations
real(dp), pointer, dimension(:) :: Tsoil        !soil temperature (K)
real(dp), pointer :: Wsno                       !snow water equivalent of the snowpack (mm)
real(dp), pointer :: zw                         !mean water table depth (m)
real(dp), pointer, dimension(:) :: ald          !active layer depth (m) 
real(dp), pointer, dimension(:) :: wetsum       !Plant wilting factor (CLM 4.0 eq 8.18) weighted by rooting zone
real(dp), pointer, dimension(:) :: psi_c        !soil water potential (mm) whtn stomata are fully closed
real(dp), pointer, dimension(:) :: psi_o        !soil water potential (mm) whtn stomata are fully open
real(dp), pointer               :: fmax         !gridcell maximum saturated fraction
real(dp), pointer               :: qseva_tot    !timestep soil evaporation (mm)
real(dp), pointer               :: qsubl_tot    !timestep soil sublimation (mm)
real(dp), pointer               :: qsdew_tot    !timestep dew (mm)
real(dp), pointer               :: qfrost_tot   !timestep frost (mm)
real(dp), pointer               :: qliq         !liquid precipitation (kg m-2 sec-1)
real(dp), pointer               :: qover        !liquid water surface runoff (kg m-2 sec-1)

!parameters
real(dp), parameter :: fdec = 0.5_dp            !decay factor (m-1) for 0.5 degree grid, see CLM 4.0 pg. 135-136
real(dp), parameter :: minV = 1.e-10            !min value for betaT (CLM 3.0 eq. 8.10)
real(dp), parameter :: Tthres = Tfreeze - 2._dp  !soil temperature below which soil water is assumed to be unavailable to plants (deg C)
real(dp), parameter :: ealpha = exp(-alpha)     !exponent of negative alpha

!local variables
real(dp) :: qsurf                                  !estimate of water flux at surface (net of rain and evaporation rates)
real(dp) :: qinfl                                  !liquid water infiltrating the surface soil layer (kg m-2 s-1)
real(dp) :: qliq0                                  !liquid water reaching soil surface (kg m-2 sec-1)
real(dp) :: s1 = 0._dp                          !term used in unsaturated conductivity calculations
real(dp) :: s2 = 0._dp                          !term used in unsaturated conductivity calculations
real(dp) :: nterm                               !numerator term used in flux calculations
real(dp) :: dterm                               !denominator term used in flux calculations
real(dp) :: Fin                                 !water flux into the soil layer
real(dp) :: Fout                                !water flux out of the soil layer
real(dp) :: ddFinTliq0                          !derivative of the water flux into the layer with respect to theta up
real(dp) :: ddFinTliq1                          !derivative of the water flux into the layer with respect to theta down
real(dp) :: ddFoutTliq1                         !derivative of the water flux out of the layer with respect to theta up
real(dp) :: ddFoutTliq2                         !derivative of the water flux out of the layer with respect to theta down
real(dp), dimension(ns:nl+1) :: ddKuTliq        !derivative of unsaturated conductivity with respect to theta
real(dp), dimension(ns:nl+1) :: ddPsiTliq       !derivative of matric potential with respect to theta
real(dp) :: ddPsiTliq1
real(dp), dimension(ns:nl+1) :: avect           !vectors for the tridiagonal solver
real(dp), dimension(ns:nl+1) :: bvect
real(dp), dimension(ns:nl+1) :: cvect
real(dp), dimension(ns:nl+1) :: rvect
real(dp), dimension(ns:nl+1) :: dTliq           !change in liquid water content (mm)
real(dp) :: fsat                                !saturated fraction (used in infiltration calculation)
real(dp) :: qinflmax                            !maximum soil infiltration capacity (kg m-2 s-1)
real(dp) :: Seffpor                             !liquid water top soil layer relative to effective porosity adjusted for saturated frac
real(dp) :: varV                                !variable V from B5 in CLM 3.5 documentation
real(dp), dimension(nl) :: Sr                   !soil wetness
integer :: k,pft,l                              !loop counters
integer :: dl                                   !number of iterations for daily timestep
real(dp) :: dft                                 !timestep for fast subdaily timestep for water flux calculations (sec)
real(dp), dimension(ns:0) :: snopor
real(dp) :: qliqsnoin
real(dp) :: qliqsnout
real(dp), dimension(nl)   :: wetf               !Plant wilting factor (CLM 4.0 eq 8.18)
real(dp), dimension(nl)   :: rootfeff           !water extraction by plants per soil level (mm s-1)
real(dp), dimension(1:nl) :: mid_term           !temporary array used for rPsi calc.
real(dp), dimension(nl,npft) :: demandperlevel
real(dp) :: ffilled_pores                       !CLM 4.0 eq. 7.65 fraction of water filled pores
real(dp) :: funsat                              !unsaturated fraction CLM 4.0 eq. 7.65
real(dp) :: zwmm                                !water table depth in mm
real(dp), dimension(nl+1)  :: vol_eq            !equilibrium volumetric water content
real(dp) :: tempi                               !temp var.
real(dp) :: temp0                               !temp var.
real(dp) :: voleq1                              !equilibrium volumetric water content
real(dp), dimension(0:nl) :: ziposmm            !soil layer interface depth (mm)
real(dp) :: dPsi_eq                             !
real(dp) :: Psi1                                !
real(dp), dimension(nl+1) :: T_wat              ! total soil water (liq + ice) (percent)
real(dp) :: dzmm_aq                             !aquifer layer thickness (mm)
real(dp) :: qseva_r                             !soil evaporation flux (mm s-1)
real(dp) :: qsdew_r                             !soil dew flux (mm s-1)
integer :: ald_layer                            !index of the lowest layer that is still frozen

!assign the pointers
zw        => sv(j)%zw
jint      => sv(j)%jint
gnl       => sv(j)%gnl
snl       => sv(j)%snl
dz        => sv(j)%dz
zpos      => sv(j)%zpos
zipos     => sv(j)%zipos
dzmm      => sv(j)%dzmm
zposmm    => sv(j)%zposmm
sorg      => sv(j)%sorg
Ksat      => sv(j)%Ksat
Tsat      => sv(j)%Tsat
Bexp      => sv(j)%Bexp
Psat      => sv(j)%Psat
Wice      => sv(j)%Wice
Wliq      => sv(j)%Wliq
Wsno      => sv(j)%Wsno
Tliq      => sv(j)%Tliq
Aliq      => surf%Aliq
Tice      => sv(j)%Tice
Tpor      => sv(j)%Tpor
Psi       => sv(j)%Psi
rootfracl => sv(j)%rootfracl
present   => sv(j)%present
betaT     => sv(j)%betaT
Ffrz      => sv(j)%Ffrz
Psi_eq    => sv(j)%Psi_eq
Ku        => sv(j)%Ku
daet      => veg%daet
rPsi      => sv(j)%rPsi
Tsoil     => sv(j)%Tsoil
ald       => surf%ald
wetsum    => sv(j)%wetsum
psi_c     => prm%psi_c
psi_o     => prm%psi_o
fmax      => soil(j)%fmax
qseva_tot => surf%qseva_tot
qsubl_tot => surf%qsubl_tot
qsdew_tot => surf%qsdew_tot
qfrost_tot => surf%qfrost_tot
qliq      => surf%qliq
qover     => surf%qover

!------

ziposmm(0:nl) = zipos(0:nl) * 1.e3        !interface depth in mm
qseva_r = qseva_tot / dtime               !evaporative flux in mm/s
qsdew_r = qsdew_tot / dtime               !dew flux in mm/s

! Update the surface ice content, snow water content, and surface infiltration

    ! water flux between snow layers
    if (-snl > 0) then

      ! update top layer water content for liquid precipitation, dew, frost and sublimation
      Wice(snl+1) = Wice(snl+1) + qfrost_tot - qsubl_tot

      if (Wice(snl+1) < 0._dp) then
        Wliq(snl+1) = Wliq(snl+1) + Wice(snl+1)
        Wice(snl+1) = 0._dp
      end if

      Wliq(snl+1) = max(0._dp,Wliq(snl+1) + qliq * dtime + qsdew_tot - qseva_tot)

      ! test whc of snow and transport water down snow column
      qliqsnoin = 0._dp

      do l = snl+1,0
        Tice(l) = min(Wice(l) / (dz(l) * pice), 1._dp)
        snopor(l)  = 1._dp - Tice(l)
        Tliq(l) = min(Wliq(l) / (dz(l) * pliq), snopor(l))
      end do

      do l = snl+1,0

        Wliq(l) = Wliq(l) + qliqsnoin

        if (l <= -1) then

          if (snopor(l) < Timp .or. snopor(l+1) < Timp) then                 !no flow
            qliqsnout = 0._dp
          else                                                               !water flows out of this layer
            qliqsnout = max(0._dp, (Tliq(l) - Smin * snopor(l)) * dz(l))
            qliqsnout = min(qliqsnout, (1._dp - Tice(l+1) - Tliq(l+1)) * dz(l+1))
          end if
        else
          qliqsnout = max(0._dp, (Tliq(l) - Smin * snopor(l)) * dz(l))
        end if

        qliqsnout = qliqsnout * pliq
        Wliq(l) = Wliq(l) - qliqsnout
        qliqsnoin = qliqsnout

      end do

      qliq = qliqsnout / dtime

      ! the total water flux to the soil cannot exceed the total water in the snowpack
      qliq = min(qliq,sum(Wliq(snl+1:0)+Wice(snl+1:0))/dtime)

    end if

qliq0 = qliq   ! amount of liquid water reaching the soil surface

!------------------------------------------

 ! Water extraction by plant roots. This extraction is the actual evapotranspiration
    do pft = 1,npft
      demandperlevel(1:gnl,pft) = daet(pft) * rootfracl(1:gnl,pft)
    end do

    do l = 1,gnl
      rootfeff(l) = sum(demandperlevel(l,1:npft)) / dtime  !mm s-1
    end do

!------------------------------------------

! Calculate the number of iterations and timestep for water flux as
! the fraction of the infiltration rate over the porosity of the top soil layer
         Tice(1) = min(Tsat(1), Wice(1) / (dz(1) * pice))
         Tpor(1) = Tsat(1) - Tice(1)

         ! estimate the surface water flux as the net of rain and evaporation
         ! FLAG try adding dew as it was left out originally JM 26.06.2011  
         qsurf = abs(qliq0 * dtime - qseva_tot + qsdew_tot)

         ! Calculate the wetness of the top soil layer to select an appropriate timestep for the hydrology calculations
         Sr(1) = min(1._dp,(Tice(1) + Tliq(1)) / Tsat(1))

         ! The sub-daily timestep is given by the net surface flux over the top layer soil wetness
               if (qsurf > 0._dp) then
                 if (Sr(1) > 1.e-4) then
                   dl = min(48,1 + nint(qsurf / Sr(1)))
                 else
                   dl = 48
                 end if
               else
                 dl = 1
               end if

               dft = dtime / real(dl)  !sec

!------------------------------------------

! Begin the loop for water flux.

dayloop : do k = 1,dl

  ! Determine which soil layer corresponds to the one immediately above the present water table height.
    jint = gnl
    do l = 2,gnl
        if(zw <= zipos(l)) then
           jint = l-1
           exit
        end if
    end do

  zwmm = zw * 1.e3  !water table depth in mm

  ! Calc infiltration and surface runoff based on current days rainfall and/or snowmelt
  ! instantaneous effective porosity as the partial volume of ice and liquid

  Tice(1:gnl) = min(Tsat(1:gnl), Wice(1:gnl) / (dz(1:gnl) * pice))        !CLM 3.0 eq. 7.121
  Tpor(1:gnl) = max(0.01_dp,Tsat(1:gnl) - Tice(1:gnl))                  !changed to be consistent with CLM4 code, was 0._dp. JM 07.04.2011
  Tliq(1:gnl) = min(Tpor(1:gnl), Wliq(1:gnl) / (dz(1:gnl) * pliq))      !CLM 3.0 eq. 7.122
  T_wat(1:gnl) = Tice(1:gnl) + Tliq(1:gnl)

  ! Fractional permeability from Niu & Yang, Hydrometeorology 2006,v. 7,p. 937
  do l = 1,gnl
      Ffrz(l) = max(0._dp,(exp(-alpha * (1._dp - Wice(l) / (Wice(l) + Wliq(l)))) - ealpha)/(1._dp - ealpha))
  end do

  ! Saturated fraction (with fractional permeability added)
          ! NOTE fmax is the percent of pixels in a grid cell whose topographic index is larger than or equal to the grid cell
          ! mean topographic index. It is read-in in arve_dgvm.f90 into the iovariables soil structure.
          ! water table depth is needed here in meters.
  fsat = (1._dp - Ffrz(1)) * fmax * exp(-0.5 * zw * fdec) + Ffrz(1)

  ! Maximum soil infiltration capacity 
    ffilled_pores = max(0.01_dp, (Tliq(1) / (max(Timp,(Tsat(1) - Tice(1))))))   !CLM 4.0 eq. 7.65
    funsat = max(0.01_dp, (1._dp - fsat))                                       !CLM 4.0 eq. 7.65
    Seffpor = max(0._dp, ((ffilled_pores - fsat) / funsat))                     !CLM 4.0 eq. 7.65
    
    varV = Bexp(1) * Psat(1) * 1._dp / (0.5 * dzmm(1))        !CLM 4.0 eq. 7.66 layer thickness in mm!

    qinflmax = Ksat(1) * (1._dp + varV * (Seffpor - 1._dp)) !CLM 4.0 eq. 7.64

  ! runoff and infiltration
   if (Tpor(1) < Timp) then    !qover is liquid water surface runoff, qliq0 is liquid water reaching soil surface.
    qover = qliq0
   else
    qover = fsat * qliq0 + (1._dp - fsat) * max(0._dp, qliq0 - qinflmax)
   end if

  ! Calc how much water inflitrates the soil (if snow covered or not)
   if (-snl == 0) then
    qinfl = qliq0 - qover - qseva_r + qsdew_r  !FLAG! JM 26.06.2011 
   else
    qinfl = qliq0 - qover              !CLM Eqn. 7.62
   end if
   
!write(*,*)
!write(*,'(a,5es12.4)')'infl',qinfl,qliq0,qover,qseva_r,qsdew_r

    ! Calculate the equilibrium water content based on the water table depth (CLM 4 Eqns 7.120 - 7.122)
    do l = 1,gnl
          if (zwmm < ziposmm(l-1)) then   !fully saturated when zw is less than the layer top

             vol_eq(l) = Tsat(l)

          ! use the weighted average from the saturated part (depth > zw) and the equilibrium solution for the
          ! rest of the layer
          else if ((zwmm < ziposmm(l)) .and. (zwmm > ziposmm(l-1))) then !CLM 4 Eqn 7.121

             !Find the equilbrium volumetric water content for the unsaturated part of the layer (Eqn. 7.122)
             tempi = 1._dp
             temp0 = ((Psat(l) - zwmm + ziposmm(l-1)) / Psat(l))**(1._dp - 1._dp/Bexp(l))
             voleq1 = Psat(l) * Tsat(l) / (1._dp - 1._dp / Bexp(l)) / (zwmm - ziposmm(l-1)) * (tempi - temp0)

             !Find the equilbrium volumetric water content for the total layer (Eqn. 7.121)
             vol_eq(l) = (voleq1 * (zwmm - ziposmm(l-1)) + Tsat(l) * (ziposmm(l) - zwmm)) / (ziposmm(l) - ziposmm(l-1))
             vol_eq(l) = min(Tsat(l), vol_eq(l))
             vol_eq(l) = max(vol_eq(l),0._dp)

          else  !layers fully above the water table (zw) (CLM 4 Eqn 7.120)

             tempi = ((Psat(l) - zwmm + ziposmm(l)) / Psat(l))**(1._dp - 1._dp / Bexp(l))
             temp0 = ((Psat(l) - zwmm + ziposmm(l-1)) / Psat(l))**(1._dp - 1._dp / Bexp(l))
             vol_eq(l) = Psat(l) * Tsat(l) / (1._dp - 1._dp / Bexp(l)) / (ziposmm(l) - ziposmm(l-1)) * (tempi - temp0)
             vol_eq(l) = max(vol_eq(l), 0._dp)
             vol_eq(l) = min(Tsat(l), vol_eq(l))

          end if

          Psi_eq(l) = Psat(l) * (max(vol_eq(l) / Tsat(l),0.01_dp))**(-Bexp(l))  !CLM 4 Eqn 7.125
          Psi_eq(l) = max(Pmax, Psi_eq(l))

    end do

    ! If water table is below soil column calculate Psi_eq for the gnl+1 layer (CLM 4 eqn 7.123)
    l = gnl
       if(jint == gnl .and. zwmm > ziposmm(l)) then

          tempi = 1._dp
          temp0 = ((Psat(l) - zwmm + ziposmm(l)) / Psat(l))**(1._dp - 1._dp / Bexp(l))
          vol_eq(l+1) = Psat(l) * Tsat(l) / (1._dp - 1._dp / Bexp(l)) / (zwmm - ziposmm(l)) * (tempi - temp0)
          vol_eq(l+1) = max(vol_eq(l+1), 0._dp)
          vol_eq(l+1) = min(Tsat(l), vol_eq(l+1))

          Psi_eq(l+1) = Psat(l) * (max(vol_eq(l+1) / Tsat(l), 0.01_dp))**(-Bexp(l))
          Psi_eq(l+1) = max(Pmax, Psi_eq(l+1))

       end if

    ! Hydraulic conductivity and soil matric potential and their derivatives
    ! Based upon Campbell 1974
    do l = 1, gnl

          s1 = 0.5_dp * (T_wat(l) + T_wat(min(gnl,l+1))) / (0.5_dp * (Tsat(l) + Tsat(min(gnl,l+1))))
          s1 = min(1._dp, s1)
          s2 = Ksat(l) * s1**(2._dp * Bexp(l) + 2._dp)  
          
          !CLM 4.0 Eqn. 7.80
          Ku(l) = (1._dp - 0.5_dp * (Ffrz(l) + Ffrz(min(gnl,l+1)))) * s2  * s1   

          !CLM 4.0 Eqn. 7.115
          ddKuTliq(l) = (1._dp - 0.5_dp * (Ffrz(l) + Ffrz(min(gnl,l+1)))) * (2._dp * Bexp(l) + 3._dp)&
                                 * s2  * 0.5_dp / Tsat(l)

         ! Calc the soil wetness
          Sr(l) = min(1._dp,((Tliq(l) + Tice(l)) / Tsat(l)))
          Sr(l) = max(0.01_dp,Sr(l))

          Psi(l) = Psat(l) * Sr(l)**(-Bexp(l))
          Psi(l) = max(Pmax, Psi(l))

          ddPsiTliq(l) = -Bexp(l) * Psi(l) / (Sr(l) * Tsat(l))

    end do

    ! aquifer (gnl+1) layer

       zposmm(gnl+1) = 0.5_dp * (zwmm + zposmm(gnl))

       if(jint < gnl) then
         dzmm_aq = dzmm(gnl)
       else
         dzmm_aq = zwmm - zposmm(gnl)
       end if

    !-------------------------------------------------
    ! Set up r, a, b, and c vectors for tridiagonal solution

    ! Node l=1 (top)

    l = 1

       Fin = qinfl
       dterm = zposmm(l+1) - zposmm(l)
       dPsi_eq = Psi_eq(l+1) - Psi_eq(l)
       nterm = (Psi(l+1) - Psi(l)) - dPsi_eq
       Fout = -Ku(l) * nterm / dterm
       ddFoutTliq1 = -(-Ku(l) * ddPsiTliq(l) + nterm * ddKuTliq(l)) / dterm
       ddFoutTliq2 = -( Ku(l) * ddPsiTliq(l+1) + nterm * ddKuTliq(l)) / dterm
       rvect(l) =  Fin - Fout - rootfeff(l)
       avect(l) =  0._dp
       bvect(l) =  dzmm(l) * (1._dp/dft) + ddFoutTliq1
       cvect(l) =  ddFoutTliq2

    !   !middle soil layers (l = 2 to gnl-1)

    do l = 2,gnl-1

          dterm = zposmm(l) - zposmm(l-1)
          dPsi_eq = Psi_eq(l) - Psi_eq(l-1)
          nterm = (Psi(l) - Psi(l-1)) - dPsi_eq
          Fin    = -Ku(l-1) * nterm / dterm
          ddFinTliq0 = -(-Ku(l-1) * ddPsiTliq(l-1) + nterm * ddKuTliq(l-1)) / dterm
          
          ddFinTliq1 = -( Ku(l-1) * ddPsiTliq(l)   + nterm * ddKuTliq(l-1)) / dterm
          dterm = zposmm(l+1) - zposmm(l)
          dPsi_eq = Psi_eq(l+1) - Psi_eq(l)
          nterm = (Psi(l+1) - Psi(l)) - dPsi_eq
          Fout = -Ku(l) * nterm / dterm
          ddFoutTliq1 = -(-Ku(l) * ddPsiTliq(l) + nterm * ddKuTliq(l)) / dterm
          ddFoutTliq2 = -( Ku(l) * ddPsiTliq(l+1) + nterm * ddKuTliq(l)) / dterm
         
          rvect(l) =  Fin - Fout - rootfeff(l)
          avect(l) = -ddFinTliq0
          bvect(l) =  dzmm(l) / dft - ddFinTliq1 + ddFoutTliq1
          cvect(l) =  ddFoutTliq2

    end do

    ! Node l=gnl (bottom)

    l = gnl

       if(l > jint) then !water table is in soil column

         dterm = zposmm(l) - zposmm(l-1)
         dPsi_eq = Psi_eq(l) - Psi_eq(l-1)
         nterm = (Psi(l) - Psi(l-1)) - dPsi_eq
         Fin = -Ku(l-1) * nterm / dterm
         ddFinTliq0 = -(-Ku(l-1) * ddPsiTliq(l-1) + nterm * ddKuTliq(l-1)) / dterm
         ddFinTliq1 = -( Ku(l-1) * ddPsiTliq(l)   + nterm * ddKuTliq(l-1)) / dterm
         Fout =  0._dp
         ddFoutTliq1 = 0._dp
         rvect(l) =  Fin - Fout - rootfeff(l)
         avect(l) = -ddFinTliq0
         bvect(l) =  dzmm(l) / dft - ddFinTliq1 + ddFoutTliq1
         cvect(l) =  0._dp

         !Next set up aquifer layer; hydrologically inactive
         rvect(l+1) = 0._dp
         avect(l+1) = 0._dp
         bvect(l+1) = dzmm_aq / dft
         cvect(l+1) = 0._dp

       else ! water table is below soil column

         !Compute aquifer soil moisture as average of layer gnl and saturation
         Sr(l) = max(0.5_dp * (1._dp + T_wat(l) / Tsat(l)), 0.01_dp)
         Sr(l) = min(1.0_dp, Sr(l))

         !Compute Psi for aquifer layer
         Psi1 = Psat(l) * Sr(l)**(-Bexp(l))
         Psi1 = max(Pmax, Psi1)

         !Compute ddPsiTliq for aquifer layer
         ddPsiTliq1 = -Bexp(l) * Psi1 / (Sr(l) * Tsat(l))

         !First set up bottom layer of soil column
         dterm = zposmm(l) - zposmm(l-1)
         dPsi_eq = Psi_eq(l) - Psi_eq(l-1)
         nterm = (Psi(l) - Psi(l-1)) - dPsi_eq
         Fin = -Ku(l-1) * nterm / dterm
         ddFinTliq0 = -(-Ku(l-1) * ddPsiTliq(l-1) + nterm * ddKuTliq(l-1)) / dterm
         ddFinTliq1 = -( Ku(l-1) * ddPsiTliq(l)   + nterm * ddKuTliq(l-1)) / dterm
         dterm = zposmm(l+1) - zposmm(l)
         dPsi_eq = Psi_eq(l+1) - Psi_eq(l)
         nterm = (Psi1 - Psi(l)) - dPsi_eq
         Fout = -Ku(l) * nterm / dterm
         ddFoutTliq1 = -(-Ku(l) * ddPsiTliq(l) + nterm * ddKuTliq(l)) / dterm
         ddFoutTliq2 = -( Ku(l) * ddPsiTliq1 + nterm * ddKuTliq(l)) / dterm

         rvect(l) =  Fin - Fout - rootfeff(l)
         avect(l) = -ddFinTliq0
         bvect(l) =  dzmm(l) / dft - ddFinTliq1 + ddFoutTliq1
         cvect(l) =  ddFoutTliq2

         !scs: next set up aquifer layer; dterm/nterm unchanged, qin=qout
         Fin = Fout  !Fout from layer gnl
         ddFinTliq0 = -(-Ku(l) * ddPsiTliq(l) + nterm * ddKuTliq(l)) / dterm
         ddFinTliq1 = -( Ku(l) * ddPsiTliq1   + nterm * ddKuTliq(l)) / dterm
         Fout = 0._dp                  ! zero-flow bottom boundary condition
         ddFoutTliq1 = 0._dp          ! zero-flow bottom boundary condition
         
         rvect(l+1) =  Fin - Fout
         avect(l+1) = -ddFinTliq0
         bvect(l+1) =  dzmm_aq / dft - ddFinTliq1 + ddFoutTliq1
         cvect(l+1) =  0._dp

       end if

  !tridiagonal solver

  call tridiag(avect(1:gnl+1),bvect(1:gnl+1),cvect(1:gnl+1),rvect(1:gnl+1),dTliq(1:gnl+1))
  
  !----------------------------
  
! recalculate the water contents of the soil layers (the update to water levels for this
! timestep occurs in aquifer.also accounting for deficit and excess also water table depth and aquifer state

  call aquifer(j,dft,dTliq(1:gnl+1),dzmm_aq)

end do dayloop

!update the daily layer 1 soil moisture (for the fire routine) NOTE: changed from Wliq to Tliq and made
                                                                !avg of top 5 soil layers (0.26m)- JM July 10th 08
Aliq(doy) = Aliq(doy) + (sum(Tliq(1:5) / 5._dp)) * (dtime / daysec)

!calc betaT (soil wetness function)

   !calculate rPsi (which is slightly different than in rest of model) 
   do l = 1,gnl
     if (Tpor(l) > 0._dp) then
       mid_term(l) = max(0.01_dp, Tliq(l) / Tpor(l))  !Eq. 8.20
     else
       mid_term(l) = 0.01_dp
     end if
   end do
 
   !Plant wilting factor (CLM 4.0 eq 8.18)
   do pft = 1,npft
     do l = 1,gnl
     
       !CLM 4.0 Eq. 8.19
       rPsi(l) = max(psi_c(pft), Psat(l) * mid_term(l)**(-Bexp(l)))

       if (Tsoil(l) > Tthres) then
         wetf(l) = min(1._dp, (psi_c(pft) - rPsi(l)) / (psi_c(pft) - psi_o(pft)) * (Tpor(l) / Tsat(l)))
       else
         wetf(l) = 0._dp
       end if  
     
     end do !soil layers

     !soil moisture function limiting transpiration (fraction) CLM eq. 8.10
     wetsum(pft) = sum(wetf(1:gnl) * rootfracl(1:gnl,pft))

     !and now the soil wetness function
     betaT(pft) = max(minV,wetsum(pft))

   end do !pft loop

! Find the active layer depth (only needs to be done 1 time per day).
! Assume that the active layer depth can be found by linearly interpolating between layers

       !find depth of frozen lowest layer
           ald_layer = nl
           do l = nl,1,-1
               if(Tsoil(l) <= Tfreeze) then
                  ald_layer = l
               else
                  exit
               end if
           end do

        !linearly interpolate a depth depending on the temperature difference between
        ! the ald_layer temperature and the layer above it.
           if (ald_layer /= 1 .and. Tsoil(ald_layer) <= Tfreeze) then
        
              ald(doy) = (Tfreeze - Tsoil(ald_layer-1))*(zpos(ald_layer+1)-zpos(ald_layer))/&
                         (Tsoil(ald_layer+1) - Tsoil(ald_layer-1)) + zpos(ald_layer-1)
   
           else if (ald_layer == 1) then !frozen all the way to the surface

           ald(doy) = 0.d0  
           
           else !no frozen ground
           
           ald(doy) = 999.d0

           end if

end subroutine soilwaterflux

!--------------------------------------------------------------------------------------

subroutine newsnowflux(j,i)

! Subroutine to calculate current day new snow flux, called before soil temperature routine

use arveparams, only : Tfreeze,Tc,pdrysnow,pwetsnow,daysec,npft,zsno_max,z0mg
use soilstate_vars, only : surf
use statevars, only : sv,dtime
use metvarsmod, only : met
use canopywater, only : canopy_interception
use pft_state_variables, only : veg

implicit none

!arguments
integer, intent(in) :: j    !index of the current gridcell
integer, intent(in) :: i    !intra-day iteration count

!pointers to the state variables
logical, pointer, dimension(:) :: present               !true if pft is present in gridcell
integer,  pointer :: snl                                !index value of top snow layer (negative is more layers)
real(dp), pointer :: Wsno                               !snow water equivalent of the snowpack (mm)
real(dp), pointer :: zsno                               !total thickness of the snowpack (m)
real(dp), pointer :: tausno                             !non-dimensional snow age
real(dp), pointer :: temp                               !air temperature (c)
real(dp), pointer :: dprec                              !total daily precipitation (mm)
real(dp), pointer, dimension(:) :: dz                   !thickness of the soil layers (m)
real(dp), pointer, dimension(:) :: zpos                 !z coordinate of the middle of the soil layer relative to soil surface (m), positive downwards
real(dp), pointer, dimension(:) :: zipos                !z coordinate of the interface between soil layers relative to soil surface (m), positive downwards
real(dp), pointer, dimension(:) :: Wliq                 !soil liquid water content at layer midpoint (mm)
real(dp), pointer, dimension(:) :: Wice                 !soil ice content at layer midpoint (mm)
real(dp), pointer, dimension(:) :: Tsoil                !soil temperature (K)
real(dp), pointer, dimension(:)  :: Tsoiln              !soil temperature (K) previous time step
real(dp), pointer  :: fpc_sum                           !instantaneous cumulative fpc_grid over all pfts for gridcell
real(dp), pointer, dimension(:)  :: fpc_grid            !instantaneous cumulative fpc_grid per pft
real(dp), pointer, dimension(:)  :: xs_can_drip         !excess canopy water (mm)
real(dp), pointer                :: fsnow               !fraction of the gridcell covered by snow (fraction)
real(dp), pointer               :: qliq                 !liquid precipitation (kg m-2 sec-1)
real(dp), pointer               :: qsno                 !snowfall (kg m-2 sec -1)
real(dp), pointer, dimension(:) :: qgrnd_l               !total rate of liquid precip reaching the ground under canopy(mm s-1)
real(dp), pointer, dimension(:) :: qgrnd_s               !total rate of solid precip reaching the ground under canopy(mm s-1)

!parameters
real(dp), parameter :: psno_new        = 100._dp        ! density of new snow (kg m-3)
integer, parameter :: m = 1                             ! exponent for fsnow calc. Suggested to be 1 for global applications

!local variables
real(dp) :: TairK               !air temperature (K)
real(dp) :: Prate               !top of the canopy precipitation rate (kg m-2 sec-1)
real(dp) :: Fpliq               !rain to snow transfer scalar (fraction)
real(dp) :: pnsno = 0._dp       !density of newly fallen snow (kg m-3)
real(dp) :: dzsno = 0._dp       !change in snow depth (m)
integer :: pft                  !counter

!pointers to the state variables
present                 => sv(j)%present
snl                     => sv(j)%snl
dz                      => sv(j)%dz
zpos                    => sv(j)%zpos                  
zipos                   => sv(j)%zipos                 
temp                    => met(j,0)%temp(i)            
Wice                    => sv(j)%Wice                  
Wliq                    => sv(j)%Wliq                  
Wsno                    => sv(j)%Wsno                  
zsno                    => sv(j)%zsno                  
tausno                  => sv(j)%tausno                
Tsoil                   => sv(j)%Tsoil                 
Tsoiln                  => sv(j)%Tsoiln
dprec                   => met(j,0)%prec      
fpc_sum                 => veg%fpc_sum
fpc_grid                => veg%fpc_grid
xs_can_drip             => veg%xs_can_drip
fsnow                   => surf%fsnow
qliq                    => surf%qliq
qsno                    => surf%qsno
qgrnd_l                 => surf%qgrnd_l
qgrnd_s                 => surf%qgrnd_s

!-----------------------------
fsnow = 0._dp

! Determine state of precipitation

! put temp in Kelvin
TairK = temp + Tfreeze

! find top of canopy precipitation rate
Prate = dprec / daysec    !mm s-1

if (TairK <= Tfreeze) then
  Fpliq = 0._dp   !then the rain to snow transfer scalar is set to 0
else if (TairK <= Tfreeze + 2._dp) then  !if it is above freezing but below the Tc then some melting
  Fpliq = max(0._dp,-54.632_dp + 0.2 * TairK)
else if (TairK <= Tfreeze + Tc) then  !!Tc is critical temperature threshold seperating rain from snow,Tc = 2.5
  Fpliq = 0.4_dp
else  !The temperature is too high for snow so the Fpliq is set to 1
  Fpliq = 1._dp
end if

!partition precip into snow or rain
qliq = Prate * Fpliq
qsno = Prate * (1._dp - Fpliq)

!--------
 !calculate canopy water (interception, canopy drip and throughfall)
 !this will determine the amount of snow/rain to reach the ground
 do pft = 1,npft
   if (present(pft)) then
      call canopy_interception(j,pft)
   else
      qgrnd_l(pft) = 0._dp
      qgrnd_s(pft) = 0._dp
   end if
 end do

!Sum the throughfall and canopy drip across all pfts here 
qliq = (1._dp - fpc_sum) * qliq + sum((qgrnd_l(:) + xs_can_drip(:) * Fpliq / dtime)* fpc_grid(:))
qsno = (1._dp - fpc_sum) * qsno + sum((qgrnd_s(:) + xs_can_drip(:) * Fpliq / dtime) * fpc_grid(:))

!-------
!if snow then determine if the snow is wet or dry
if (qsno > 0._dp) then

 !calc the density of the NEW snow based upon air temperature
  if (TairK > Tfreeze + 2._dp) then
    pnsno = pwetsnow
  else if (TairK <= Tfreeze + 2._dp .and. TairK > Tfreeze - 15._dp) then
    pnsno = pdrysnow + 1.7_dp * (TairK - Tfreeze + 15._dp)**1.5_dp
  else
    pnsno = pdrysnow
  end if

!update the change in snow depth, total thickness of the snowpack and water content
  dzsno = (qsno * dtime) / pnsno
  zsno  = zsno + dzsno

  ! Some areas have climate that does not allow summer melting of snow which results in the
  ! the build up of snow (i.e. glacier building). The snow can not presently be shifted away by wind
  ! so we just set a maximum snow depth. Any snow above that is assumed to be blown away (i.e. it is not
  ! added to soil). JM 21.04.2011
  if (zsno > zsno_max) then
     zsno = zsno_max
     qsno = 0._dp
     dzsno = 0._dp
  end if
  
  Wsno  = Wsno + qsno * dtime

  if (zsno >= 0.01_dp) then
    if (snl+1 == 1) then
    
      !initialize a snow layer
      snl           = -1
      dz(snl+1)     =  zsno
      zpos(snl+1)   = -0.5 * dz(snl+1)
      zipos(snl)    = -dz(snl+1)
      tausno        =  0._dp
      Tsoil(snl+1)  =  min(Tfreeze,TairK)
      Tsoiln(snl+1) =  min(Tfreeze,TairK)  !this is for radflux (which requires a previous timestep value)
      Wice(snl+1)   =  Wsno
      Wliq(snl+1)   =  0._dp

    else
    
      !update existing top snow layer
      Wice(snl+1) = Wice(snl+1) + qsno * dtime
      dz(snl+1) = dz(snl+1) + dzsno

      zpos(snl+1)= zipos(snl+1) - 0.5 * dz(snl+1)
      zipos(snl) = zipos(snl+1) - dz(snl+1)

    end if !snow layer exists?
  end if !need to initialize a layer?
end if !new snow today?

if (Wsno > 0._dp) then
  ! fractional area of the gridcell covered by snow (Niu and Yang 2007)  
   fsnow = tanh(zsno / (2.5_dp * z0mg * ((min(Wsno/zsno, 800._dp)) / psno_new)**m))
end if

end subroutine newsnowflux

!--------------------------------------------------------------------------------------

subroutine snowdynamics(j)

use arveparams, only : Tfreeze,Cice,Cliq,Lf,snomax,pwetsnow,ns,nl,pliq,pice,z0mg
use statevars, only : sv,dtime
use soilstate_vars, only : ithaw,surf

implicit none

!arguments
integer, intent(in) :: j                     !index of the current gridcell

!pointers to the state variables
integer , pointer :: snl                     !index value of top snow layer (negative is more layers)
real(dp), pointer :: Wsno                    !snow water equivalent of the snowpack (mm)
real(dp), pointer :: Wsno_old                !snow water equivalent of the snowpack (mm) from previous timestep
real(dp), pointer :: zsno                    !total thickness of the snowpack (m)
real(dp), pointer :: tausno                  !non-dimensional snow age
real(dp), pointer :: tausno_old              !non-dimensional snow age of previous timestep
real(dp), pointer, dimension(:) :: dz        !thickness of the soil layers (m)
real(dp), pointer, dimension(:) :: zpos      !z coordinate of the middle of the soil layer relative to soil surface (m), positive downwards
real(dp), pointer, dimension(:) :: zipos     !z coordinate of the interface between soil layers relative to soil surface (m), positive downwards
real(dp), pointer, dimension(:) :: Wliq      !soil liquid water content at layer midpoint (mm)
real(dp), pointer, dimension(:) :: Wice      !soil ice content at layer midpoint (mm)
real(dp), pointer, dimension(:) :: Tsoil     !soil temperature (K)
real(dp), pointer, dimension(:) :: Tsoiln    !soil temperature (K) from previous timestep
real(dp), pointer, dimension(:) :: fice0     !layer ice fraction, previous timestep
real(dp), pointer               :: fsnow     !fraction of the gridcell covered by snow (fraction)

!local variables
integer :: l
logical :: combdivflag
real(dp), dimension(-snomax:0) :: Cr          !snow fractional compaction rate (s-1)
real(dp) :: Cr1                               !compaction rate due to destructive metamorphism (s-1)
real(dp) :: Cr2                               !compaction rate due to overburden (s-1)
real(dp) :: Cr3                               !compaction rate due to melting (s-1)
real(dp) :: c1                                !coefficients for calculating the compaction rates
real(dp) :: c2
real(dp), parameter :: c3 = 2.777e-6          !(s-1)
real(dp), parameter :: c4 = 0.04_dp           !(K-1)
real(dp), parameter :: c5 = 0.08_dp           !(K-1)
real(dp), parameter :: c6 = 0.023_dp          !(m3 kg-1)
real(dp), parameter :: n0 = 9.e5              !(kg s m-2)
real(dp), parameter :: tao_0 = 1.e-6          !(s-1)
real(dp) :: Ps                                !snow load pressure (?)
real(dp) :: Nu                                !viscosity coefficient for overburden compaction (kg s m-2)
real(dp) :: r1                                !for soil age calculations (grain growth
real(dp) :: r2                                !effect of melt
real(dp) :: r3                                !, soot and dirt)
real(dp) :: del_tausno                        !increase in snow age this timestep
real(dp) :: snonew                            !increase in snow amount since last timestep (mm)
integer :: nsno                               !current number of snow layers
integer :: ll                                 !index of the lower snow layer used in division calculations
integer :: ul                                 !index of the upper snow layer used in division calculations
integer :: sl                                 !negative index value for snow layer division
integer :: ls                                 !index for layer shift down
real(dp), dimension(snomax) :: dzsno          !snow layer thickness (m)
real(dp), dimension(snomax) :: dzsno_old      !snow layer thickness from previous timestep (m)
real(dp), dimension(snomax) :: wlsno          !Wliq for snow layers (m3 m-3)
real(dp), dimension(snomax) :: wisno          !Wice for snow layers (m3 m-3)
real(dp), dimension(snomax) :: tsnow          !snow layer temperature (K)
real(dp) :: xice
real(dp) :: xliq
real(dp) :: xfra
real(dp) :: xsno
real(dp) :: Wtot                              !total water mass of the snow layer (kg m-2)
real(dp) :: Tpor                              !air filled porosity of the snow layer (fraction)
real(dp), dimension(ns:nl) :: fice            !layer ice fraction, current timestep
real(dp), dimension(-snomax:0):: psno            !layer snow density (kg m-3)

!snow layer minimum and maximum thickness limits (m)
!note there is no limit to the thickness of the bottom layer
!layers get thinner towards the top, as for soil

real(dp), dimension(snomax),  parameter :: dzmin = [ 0.010_dp, 0.015_dp, 0.025_dp, 0.055_dp, 0.115_dp ] !min
real(dp), parameter :: psno_new        = 100._dp        ! density of new snow (kg m-3)
integer, parameter :: m = 1                             ! exponent for fsnow calc. Suggested to be 1 for global applications

!pointers to the state variables
snl   => sv(j)%snl
dz    => sv(j)%dz
zpos  => sv(j)%zpos
zipos => sv(j)%zipos
Wice  => sv(j)%Wice
Wliq  => sv(j)%Wliq
Wsno  => sv(j)%Wsno
Wsno_old => sv(j)%Wsno_old
zsno  => sv(j)%zsno
tausno=> sv(j)%tausno
tausno_old => sv(j)%tausno_old
Tsoil => sv(j)%Tsoil
Tsoiln => sv(j)%Tsoiln
fice0 => sv(j)%fice0
fsnow => surf%fsnow

!--------------------------------------------
fsnow = 0._dp

! Begin calculations

! Snow layer compaction

        !write(*,*)'compact',dz(0),wice(0)

Ps = 0._dp
combdivflag = .false.

do l = snl+1,0

  Wtot = Wice(l) + Wliq(l)
  Tpor = 1._dp - (Wice(l) / pice + Wliq(l) / pliq) / dz(l)
  
        !write(*,'(a,i4,5f12.4)')'start snow dyn',l,Tpor,(Wice(l) + Wliq(l)) / dz(l),Wice(l),Wliq(l),dz(l)

  !to account for the current layer add half of the ice and water contents of the layer being
  !compacted (CLM 4 Eqn 7.44 15.06.2010 JM)
  Ps = Ps + Wtot * 0.5

  !allow compaction only for non-saturated node and higher ice lens node.
  if (Tpor > 0.001_dp .and. Wice(l) > 0.1_dp) then

    !compaction from metamorphism

    if (Wice(l) / dz(l) <= 100._dp) then  !kg m-3
      c1 = 1._dp
    else
      c1 = exp(-0.046_dp * (Wice(l) / dz(l) - 100._dp))
    end if

    if (Wliq(l) / dz(l) > 0.01_dp) then
      c2 = 2._dp
    else
      c2 = 1._dp
    end if

    Cr1 = -c3 * c1 * c2 * exp(-c4 * (Tfreeze - Tsoil(l)))

    !compaction from overburden

    Nu = n0 * exp(c5 * (Tfreeze - Tsoil(l)) + c6 * (Wice(l) / dz(l)))

    Cr2 = -Ps / Nu

    !compaction from melting

    fice(l) = Wice(l) / (Wice(l) + Wliq(l))

    if (ithaw(l) == 1) then
      Cr3 = -1._dp / dtime * max(0._dp, (fice0(l) - fice(l)) / fice0(l))
    else
      Cr3 = 0._dp
    end if

    !total compaction rate (s-1)

    Cr(l) = Cr1 + Cr2 + Cr3

    !update snow layer thickness

    dz(l) = dz(l) * (1._dp + Cr(l) * dtime)
    
                !write(*,'(a,i4,2f12.4)')'new layer thickness',l,dz(l),(Wice(l)+Wliq(l))/dz(l)

    !check if this makes the layer denser than ice
    psno(l) = (Wice(l) + Wliq(l)) / dz(l)
    
    if (psno(l) > pice) then !if so then set the thickness so that it is an ice layer.
        
        dz(l) = (Wice(l) + Wliq(l))/ pice
        
    end if
    
      !if the snow layer thickness is less than 0 then remove the layer.
      if (dz(l) < 0._dp) then
        dz(l) = 0._dp
        Wice(l) = 0._dp
        Wliq(l) = 0._dp
      end if

  else

  end if

  Ps = Ps + Wtot !sum(Wice(snl+1:l-1) + Wliq(snl+1:l-1))  !amount of ice and water in the layers above current layer (kg m-2)

end do

!-------------------------------
!snow layer combination

!first pass to determine if any layer is nearly melted
do l = snl+1,0

  if (Wice(l) <= 0.01_dp) then
  
!write(*,'(a,i5,5f15.3)')'start compaction:layer combo ',l,Wice(l),Wliq(l),dz(l),(Wice(l) + Wliq(l)) / dz(l),Tsoil(l)-Tfreeze

    Wliq(l+1) = Wliq(l+1) + Wliq(l)
    Wice(l+1) = Wice(l+1) + Wice(l)

    Tsoil(l) = Tsoil(l+1)

    Wice(l)  = Wice(l+1)
    Wliq(l)  = Wliq(l+1)
    dz(l)    = dz(l+1)
    snl      = snl + 1
    
   if (dz(l) <= 0._dp) cycle
   
   !check if this makes either layer denser than ice
    psno(l) = (Wice(l) + Wliq(l)) / dz(l)

     if (psno(l) > pice) then !if so then set the thickness so that it is an ice layer.
      dz(l) = (Wice(l) + Wliq(l))/ pice  
    end if

        !write(*,'(a,i5,5f15.3)')'end compaction:layer combo ',l,Wice(l),Wliq(l),dz(l),(Wice(l) + Wliq(l)) / dz(l),Tsoil(l)-Tfreeze
   
  end if
end do

if (snl == 0) then
  !no snow layers remain
  Wsno = 0._dp
  zsno = 0._dp
  return

else

  Wsno = sum(Wice(snl+1:0) + Wliq(snl+1:0))
  zsno = sum(dz(snl+1:0))

  if (zsno < 0.01_dp) then       !total snowpack is too thin

    Wsno = sum(Wice(snl+1:0))
    Wliq(1) = Wliq(1) + sum(Wliq(snl+1:0))
    snl = 0
    return

  else
  
    !check that all remaining layers are at minimum thickness, or combine
    do l = snl+1,0

      sl = l - snl

               ! write(*,'(a,4i5,4f12.4)')'thickness',l,sl,snl,snl+1,dz(l),dzmin(sl),Tsoil(l)-Tfreeze,(Wice(l)+Wliq(l))/dz(l)

      if (dz(l) < dzmin(sl)) then      !layer too thin, combine it...
        if (l == snl+1) then          !surface layer, combine with layer below
          ul = l
          ll = l+1
        else if (l == 0) then        !bottom layer, combine with layer above
          ul = l-1
          ll = l
        else                          !combine with the thinnest neighboring layer
          if (dz(l+1) < dz(l-1)) then
            ul = l
            ll = l+1
          else
            ul = l-1
            ll = l
          end if
        end if

        call comb(dz(ll),Wliq(ll),Wice(ll),Tsoil(ll),dz(ul),Wliq(ul),Wice(ul),Tsoil(ul))
         combdivflag = .true.

        if (ll-1 > snl+1) then
          do ls = ll-1,snl+2,-1
            Tsoil(ls) = Tsoil(ls-1)
            Wice(ls)  = Wice(ls-1)
            Wliq(ls)  = Wliq(ls-1)
            dz(ls)    = dz(ls-1)
          end do
        end if

        snl = snl + 1
        if (snl >= -1) exit

      end if
    end do
  end if
end if

!---------------------------------------------------
!snow layer division

!There are four state variables that describe the snowpack:
!thickness (dz), temperature (Tsoil), water content (Wliq), ice content (Wice)

!The top snow layer is always the thinnest layer, thickening with depth.
!The snow layer against the ground (bottom layer) has index 0, with layers numbered negatively upwards.

!When performing layer division, we use temporary variables to avoid needing to renumber
!the layers each time a new one is formed.
!At the end of the routine, we translate all of the temporary variables back to the state variables.

!New temperature scheme from CLM 4.0 implimented. Without this temp scheme it was noticed that ARVE had
!temp spikes after a division. This new scheme will correct that problem. See CLM 4.0 pg. 134-135.

do l = snl+1,0
  sl = l - snl
  dzsno(sl) = dz(l)
  dzsno_old(sl) = dz(l)
  wlsno(sl) = Wliq(l)
  wisno(sl) = Wice(l)
  tsnow(sl) = Tsoil(l)
      
                !write(*,'(i5,a18,2i5,3f8.3)')-snl,'start division: layers, assign:',l,sl,dzsno(sl),Tsoil(l)-Tfreeze,(Wice(l) + Wliq(l)) / dz(l)

end do

nsno = abs(snl)

if (nsno == 1) then
  if (dzsno(1) > 0.03_dp) then          !subdivide layer 2 from layer 1 (1,2)

    nsno = 2
    ul = 1
    ll = 2

    dzsno(1) = 0.5_dp * dzsno(1)
    wisno(1) = 0.5_dp * wisno(1)
    wlsno(1) = 0.5_dp * wlsno(1)

    dzsno(2) = dzsno(1)
    wisno(2) = wisno(1)
    wlsno(2) = wlsno(1)
    tsnow(2) = tsnow(1)
      
  end if

else if (nsno > 1) then

  if (dzsno(1) > 0.02_dp) then

    xsno = dzsno(1) - 0.02_dp  !drr

    xfra = xsno / dzsno(1)    !propor
    xice = xfra * wisno(1)    !zwice
    xliq = xfra * wlsno(1)    !zwliq

    xfra = 0.02_dp / dzsno(1)
    wisno(1) = xfra * wisno(1)
    wlsno(1) = xfra * wlsno(1)
    dzsno(1) = 0.02_dp

    call comb(dzsno(2),wlsno(2),wisno(2),tsnow(2),xsno,xliq,xice,tsnow(1))
    combdivflag = .true.

    if (nsno <= 2 .and. dzsno(2) > 0.07_dp) then !subdivide layer 3 from layer 2 (1,2,3)

      nsno = 3

      dzsno(2) = 0.5_dp * dzsno(2)
      wisno(2) = 0.5_dp * wisno(2)
      wlsno(2) = 0.5_dp * wlsno(2)

      dzsno(3) = dzsno(2)
      wisno(3) = wisno(2)
      wlsno(3) = wlsno(2)
     
      !adjust the new layer temperature as in CLM 4.0 Eqn. 7.59
      if (dzsno_old(1) + dzsno_old(2) > 0._dp) then
        tsnow(3) = tsnow(2) - ((tsnow(1) - tsnow(2)) / ((dzsno_old(1) + dzsno_old(2)) * 0.5_dp)) * (dzsno(2) * 0.5_dp)
      else
        tsnow(3) = tsnow(2)
      end if
        
      if (tsnow(3) >= Tfreeze) then  !Eqn 7.60
       tsnow(3) = tsnow(2)
      else
        if (dzsno_old(1) + dzsno_old(2) > 0._dp) then
          tsnow(2) = tsnow(2) - ((tsnow(1) - tsnow(2)) / ((dzsno_old(1) + dzsno_old(2)) * 0.5_dp)) * (dzsno(2) * 0.5_dp)
        end if
      end if
      
    end if
  end if
end if

if (nsno > 2) then
  if (dzsno(2) > 0.05_dp) then

    xsno = dzsno(2) - 0.05_dp
    xfra = xsno / dzsno(2)
    xice = xfra * wisno(2)
    xliq = xfra * wlsno(2)

    xfra = 0.05_dp / dzsno(2)
    wisno(2) = xfra * wisno(2)
    wlsno(2) = xfra * wlsno(2)
    dzsno(2) = 0.05_dp

    call comb(dzsno(3),wlsno(3),wisno(3),tsnow(3),xsno,xliq,xice,tsnow(2))
    combdivflag = .true.

    if (nsno <= 3 .and. dzsno(3) > 0.18_dp) then !subdivide layer 4 from layer 3 (1,2,3,4)

      nsno = 4

      dzsno(3) = 0.5_dp * dzsno(3)
      wisno(3) = 0.5_dp * wisno(3)
      wlsno(3) = 0.5_dp * wlsno(3)

      dzsno(4) = dzsno(3)
      wisno(4) = wisno(3)
      wlsno(4) = wlsno(3)
      
      !adjust the new layer temperature as in CLM 4.0 Eqn. 7.59
      if (dzsno_old(2) + dzsno_old(3) > 0._dp) then
        tsnow(4) = tsnow(3) - ((tsnow(2) - tsnow(3)) / ((dzsno_old(2) + dzsno_old(3)) * 0.5_dp)) * (dzsno(3) * 0.5_dp)
      else 
        tsnow(4) = tsnow(3)
      end if  
 
        if (tsnow(4) >= Tfreeze) then !Eqn 7.60
         tsnow(4) = tsnow(3)
        else

           if (dzsno_old(2) + dzsno_old(3) > 0._dp) then
             tsnow(3) = tsnow(3) - ((tsnow(2) - tsnow(3)) / ((dzsno_old(2) + dzsno_old(3)) * 0.5_dp)) * (dzsno(3) * 0.5_dp)
           end if  

        end if
       
    end if
  end if
end if

if (nsno > 3) then
  if (dzsno(3) > 0.11_dp) then

    xsno = dzsno(3) - 0.11_dp
    xfra = xsno / dzsno(3)
    xice = xfra * wisno(3)
    xliq = xfra * wlsno(3)

    xfra = 0.11_dp / dzsno(3)
    wisno(3) = xfra * wisno(3)
    wlsno(3) = xfra * wlsno(3)
    dzsno(3) = 0.11_dp

    call comb(dzsno(4),wlsno(4),wisno(4),tsnow(4),xsno,xliq,xice,tsnow(3))
    combdivflag = .true.

    if (nsno <= 4 .and. dzsno(4) > 0.41_dp) then !subdivide layer 5 from layer 4 (1,2,3,4,5)

      nsno = 5

      dzsno(4) = 0.5_dp * dzsno(4)
      wisno(4) = 0.5_dp * wisno(4)
      wlsno(4) = 0.5_dp * wlsno(4)

      dzsno(5) = dzsno(4)
      wisno(5) = wisno(4)
      wlsno(5) = wlsno(4)
     
      !adjust the new layer temperature as in CLM 4.0 Eqn. 7.59 (modified for cases where the snow jumps to 
      !five layers (result of using weather generator)
      if (dzsno_old(3) + dzsno_old(4) > 0._dp) then
        tsnow(5) = tsnow(4) - ((tsnow(3) - tsnow(4)) / ((dzsno_old(3) + dzsno_old(4)) * 0.5_dp)) * (dzsno(4) * 0.5_dp)
      else
        tsnow(5) = tsnow(4)
      end if  
      
        if (tsnow(5) >= Tfreeze) then !Eqn 7.60
         tsnow(5) = tsnow(4)
        else
         
         if (dzsno_old(3) + dzsno_old(4) > 0._dp) then
           tsnow(4) = tsnow(4) - ((tsnow(3) - tsnow(4)) / ((dzsno_old(3) + dzsno_old(4)) * 0.5_dp)) * (dzsno(4) * 0.5_dp)
         end if  
      
        end if

    end if
  end if
end if

if (nsno > 4) then
  if (dzsno(4) > 0.23_dp) then

    xsno = dzsno(4) - 0.23_dp
    xfra = xsno / dzsno(4)
    xice = xfra * wisno(4)
    xliq = xfra * wlsno(4)

    xfra = 0.23_dp / dzsno(4)
    wisno(4) = xfra * wisno(4)
    wlsno(4) = xfra * wlsno(4)
    dzsno(4) = 0.23_dp

    call comb(dzsno(5),wlsno(5),wisno(5),tsnow(5),xsno,xliq,xice,tsnow(4))
    combdivflag = .true.

  end if
end if

!---------------------------------
!update state variables

        snl = -nsno

        do l = snl+1,0
          sl = l - snl
          dz(l)    = dzsno(sl)
          Wice(l)  = wisno(sl)
          Wliq(l)  = wlsno(sl)
          Tsoil(l) = tsnow(sl)

          if (Tsoiln(snl+1) == 0._dp) then
            Tsoiln(snl+1) = min(Tfreeze,Tsoil(snl+1))
          end if

          !check if any layer is denser than ice
          psno(l) = (Wice(l) + Wliq(l)) / dz(l)

          if (psno(l) > pice) then !if so then set the thickness so that it is an ice layer.

              dz(l) = (Wice(l) + Wliq(l))/ pice
          end if
          
        end do
        
        zsno = sum(dz(snl+1:0))
        
!recalculate layer depths and interfaces from the soil surface upwards
        do l = 0,snl+1,-1
          zpos(l) = zipos(l) - 0.5 * dz(l)
          zipos(l-1) = zipos(l) - dz(l)
        end do

!update snow age -added May 5 08 JM
        snonew = max(0._dp,(Wsno - Wsno_old))

        if (snonew > 10.) then !whole snow surface is renewed

                   tausno = 0._dp

            else !snow is aging

                !effect of grain growth due to vapour diffusion
                r1 = 5.e3 * (1. / Tfreeze - 1. / Tsoil(snl+1))
                
                !effect of near and at freezing of melt water
                r2 = min(0._dp,r1 * 10._dp)

                r1 = exp(r1)
                r2 = exp(r2)
                
                !effect of dirt and soot    
                 r3 = 0.3                   

                !find the snow aging for this timestep
                if (Wsno > 0._dp .and. Wsno <= 800.) then
                  del_tausno = tao_0 * (r1 + r2 + r3) * dtime
                else
                  del_tausno = 0.
                end if

                tausno = max(0._dp,((tausno_old + del_tausno)*(1. - 0.1 * max(0._dp,(Wsno - Wsno_old)))))
                
        end if
        
        if (Wsno > 0._dp) then
          ! fractional area of the gridcell covered by snow (Niu and Yang 2007)   
           fsnow = tanh(zsno / (2.5_dp * z0mg * ((min(Wsno/zsno, 800._dp)) / psno_new)**m))
        end if


        !if the snow layers are combined or subdivided and the number of snow layers is less than the max
           !set taosno to zero
        if (combdivflag .and. snl < snomax) tausno = 0._dp

        !set the previous timestep variables for next timestep
         tausno_old = tausno
          Wsno_old = Wsno


end subroutine snowdynamics

!--------------------------------------------------------------------------------------

subroutine comb(dz,Wliq,Wice,t,dz2,Wliq2,Wice2,t2)

use arveparams, only : Tfreeze, Lf, Cice, Cliq,pice

implicit none

!Combines characteristics of adjoining snow layers, operating on state vars: dz, t, wliq, wice.
!The combined temperature is based on the equation:
! the sum of the enthalpies of the two elements = that of the combined element.
!from CLM3

!arguments

real(dp), intent(inout) :: dz    ! nodal thickness of no. 1 element being combined [m]
real(dp), intent(inout) :: Wliq  ! liquid water of element 1
real(dp), intent(inout) :: Wice  ! ice of element 1 [kg/m2]
real(dp), intent(inout) :: t     ! nodal temperature of elment 1 [K]

real(dp), intent(in)    :: dz2   ! nodal thickness of no. 2 element being combined [m]
real(dp), intent(in)    :: Wliq2 ! liquid water of element 2 [kg/m2]
real(dp), intent(in)    :: Wice2 ! ice of element 2 [kg/m2]
real(dp), intent(in)    :: t2    ! nodal temperature of element 2 [K]

!local variables

real(dp) :: dzc   ! Total thickness of nodes 1 and 2 (dzc=dz+dz2).
real(dp) :: wliqc ! Combined liquid water [kg/m2]
real(dp) :: wicec ! Combined ice [kg/m2]
real(dp) :: Tc    ! Combined node temperature [K]
real(dp) :: h     ! enthalpy of element 1 [J/m2]
real(dp) :: h2    ! enthalpy of element 2 [J/m2]
real(dp) :: hc    ! temporary

!--------------

dzc = dz + dz2
Wicec = Wice + Wice2
Wliqc = Wliq + Wliq2

h  = (Cice * Wice  + Cliq * Wliq)  * (t  - Tfreeze) + Lf * Wliq
h2 = (Cice * Wice2 + Cliq * Wliq2) * (t2 - Tfreeze) + Lf * Wliq2

hc = h + h2

!This follows CLM 4.0 Eqn 7.55 to correct an error and ensure that enthalpy is always
!conserved during combination.
  Tc = Tfreeze + (hc - Lf * Wliqc) / (Cice * Wicec + Cliq * Wliqc)

dz   = dzc
Wice = Wicec
Wliq = Wliqc
t = Tc

end subroutine comb

!--------------------------------------------------------------------------------------

end module soilhydrologymod

