module soiltemperaturemod

! Calculates snow and soil temperatures excluding phase change (soilphasechg.f90)

implicit none

public :: soilthermalprop
public :: surfaceheatflux
public :: soiltemperature

contains

!-----------------------------------------------------------------------
subroutine soilthermalprop(j)

! This has been updated based upon Ballard & Arp 2005, J. Environ Eng. Sci 4:549-588 to
! account for organic matter and coarse fragments in the soil
! May 14 08 removed sections to soilinit (Kdry, Ksolid, Csolid) to make more
! efficient -JM

use arveparams, only : dp,nl,Tfreeze,Kice,Kair,Cair,Cice,Kliq,Cliq,Lf,pice,pwetsnow
use statevars, only : sv
use soilstate_vars, only : Cl,Kh,Kl,surf

implicit none

!arguments
integer, intent(in) :: j                !index of the current gridcell

!pointers
integer,  pointer :: snl                !index value of top snow layer (negative is more layers)
integer,  pointer :: gnl                !index value of lowest 'soil' layer
real(dp), pointer :: Wsno               !snow water equivalent of the snowpack (mm)
real(dp), pointer, dimension(:) :: dz      !thickness of the soil layers (m)
real(dp), pointer, dimension(:) :: zpos    !z coordinate of the middle of the soil layer relative to soil surface (m), positive downwards
real(dp), pointer, dimension(:) :: zipos   !z coordinate of the interface between soil layers relative to soil surface (m), positive downwards
real(dp), pointer, dimension(:) :: sand    !soil sand content (percent)
real(dp), pointer, dimension(:) :: clay    !soil clay content (percent)
real(dp), pointer, dimension(:) :: sorg    !soil organic matter content (percent)
real(dp), pointer, dimension(:) :: rock    !soil rock content (MASS fraction)
real(dp), pointer, dimension(:) :: bulk    !soil bulk density (kg m-3)
real(dp), pointer, dimension(:) :: Tsat    !soil water volumetric water content at saturation (fraction)
real(dp), pointer, dimension(:) :: Bexp    !soil water B exponent used in the Brooks & Corey Pedotransfer Functions
real(dp), pointer, dimension(:) :: Psat    !soil water matric potential at saturation (mm)
real(dp), pointer, dimension(:) :: Wliq    !soil liquid water content at layer midpoint (mm)
real(dp), pointer, dimension(:) :: Wice    !soil ice content at layer midpoint (mm)
real(dp), pointer, dimension(:) :: Tsoil   !soil temperature (K)
real(dp), pointer, dimension(:) :: Tliq    !soil liquid water content (fraction)
real(dp), pointer, dimension(:) :: Tice    !soil ice content (fraction)
real(dp), pointer, dimension(:) :: Vcf     !FRACTION of coarse fragments (rocks) by volume
real(dp), pointer, dimension(:) :: Vom     !volumetric fraction of soil organic matter
real(dp), pointer, dimension(:) :: Vsand   !volumetric fraction of sand
real(dp), pointer, dimension(:) :: Kdry    !soil thermal conductivity of dry natural soil (W m-1 K-1)
real(dp), pointer, dimension(:) :: Ksolid  !soil mineral thermal conductivity at layer midpoint (W m-1 K-1)
real(dp), pointer, dimension(:) :: Csolid  !soil solids volumetric heat capacity at layer midpoint (J m-3 K-1)

!parameters
real(dp), parameter :: alpha =  0.24_dp  ! +/- 0.04 constant for Knum calc.from Balland and Arp
real(dp), parameter :: beta  = 18.3_dp   ! +/- 1.1 constant for Knum calc.from Balland and Arp

!local variables
integer :: l
real(dp) :: Knum          !Kersten number (CLM technical note eqn 6.63, pg. 94)
real(dp) :: Ktsat         !soil saturated thermal conductivity (W m-1 K-1)
real(dp) :: psno          !snow bulk density
real(dp) :: Sr            !soil wetness with respect to saturation
real(dp) :: Knum1         !interm variables.
real(dp) :: Knum2
real(dp) :: term_a        !Balland and Arp snow thermal conductivity adjustable term

!point to pointers
gnl   => sv(j)%gnl
snl   => sv(j)%snl
dz    => sv(j)%dz
zpos  => sv(j)%zpos
zipos => sv(j)%zipos
sand  => sv(j)%sand
clay  => sv(j)%clay
sorg  => sv(j)%sorg
bulk  => sv(j)%bulk
rock  => sv(j)%rock
Tsat  => sv(j)%Tsat
Bexp  => sv(j)%Bexp
Psat  => sv(j)%Psat
Wice  => sv(j)%Wice
Wliq  => sv(j)%Wliq
Wsno  => sv(j)%Wsno
Tsoil => sv(j)%Tsoil
Tliq  => sv(j)%Tliq
Tice  => sv(j)%Tice
Vcf   => sv(j)%Vcf
Vom   => sv(j)%Vom
Vsand => sv(j)%Vsand
Kdry  => sv(j)%Kdry
Csolid  => sv(j)%Csolid
Ksolid  => sv(j)%Ksolid

!------------------->>

! While the bottom layers (gnl+1,nl) are not hydrologically active, they do freeze/thaw
! so soil thermal properties are calculated for all soil levels.

!thermal conductivity and heat capacity of soil solids (from Hillel, Env. Soil Phys Book) and snow
do l = snl+1,nl 

  !snow layers instantaneous conductivity and heat capacity
  if (l < 1) then

            !snow bulk density (kg m-3) eqn 6.66 CLM
            psno = (Wice(l) + Wliq(l)) / dz(l)

            ! snow thermal conductivity CLM original formulation
!            Kh(l) = Kair + (7.75e-5 * psno + 1.105e-6 * psno * psno) * (Kice - Kair) !(W m-1 K-1) eqn 6.65 CLM

            ! snow thermal conductivity from Balland and Arp 2005.
            term_a = 0.3
            Kh(l) = ((term_a * Kice - Kair) * psno * 1.e3 + Kair * pice * 1.e3) / (pice * 1.e3 - (1._dp - term_a) * psno * 1.e3)

            !snow volumetric heat capacity
            Cl(l) = (Wice(l) / dz(l)) * Cice + (Wliq(l) / dz(l)) * Cliq  !(J m-3 K-1) eqn 6.69 CLM


  else     !all soil layers

          !saturated thermal conductivity
          if (Tsoil(l) > Tfreeze) then !non-frozen soils
             Ktsat = Ksolid(l)**(1._dp - Tsat(l)) * Kliq**Tsat(l) !eqn 12 in Balland and Arp
          else  !frozen soils
             Ktsat = Ksolid(l)**(1._dp - Tsat(l)) * Kliq**Tliq(l) * Kice**(Tsat(l) - Tliq(l)) ! eqn. 13 in Balland & Arp
          end if

            !wetness with respect to saturation
            Sr = min((Tliq(l) + Tice(l)) / Tsat(l), 1._dp)  !eqn. 6

            !Kersten number
            if (Tsoil(l) > Tfreeze) then    !thawed soils (eqn 17 balland & arp)

              Knum1 = (0.5 * (1._dp + Vom(l) - (alpha * Vsand(l)) - Vcf(l)))
              Knum2 = (((1. / (1. + exp(-beta * Sr)))**3.) - ((1. - Sr) / 2.)**3.)**(1. - Vom(l))

              Knum = Sr**Knum1 * Knum2

            else  !frozen or partially frozen soils   (eqn 18 balland & arp)

              Knum = Sr**(1. + Vom(l))

            end if

            !instantaneous soil thermal conductivity at layer midpoint
            Kh(l) = Knum * (Ktsat - Kdry(l)) + Kdry(l)    !(eqn 5 balland & arp)

            !soil heat capacity at layer midpoint
            !NOTE: new formulation - Apr 18 08 JM
            Cl(l) = ((Csolid(l) * (1. - Tsat(l))) + (Tice(l) * 1000._dp * Cice) + (Tliq(l) * 1000._dp * Cliq)) !+ &
                !(1000._dp * Cair * (Tsat(l) - Tliq(l) - Tice(l))) ) !NOTE the last part is to account for the heat capacity of the air
            !within the pore space, ignored at present.

  end if

end do

!find soil/snow thermal conductivity across layer boundary (W m-1 K-1)
do l = snl+1,nl  
   if (l < nl) then
        Kl(l) = (Kh(l) * Kh(l+1) * (zpos(l+1) - zpos(l))) /  &
                (Kh(l) * (zpos(l+1) - zipos(l)) + Kh(l+1) * (zipos(l) - zpos(l)))
   else
        Kl(l) = 0._dp  !no flux bottom condition
   end if
end do


!set the volumetric heat capacity of a snow layer if it is the inital one.
      if (snl+1 == 1 .and. Wsno > 0._dp .and. dz(snl+1) /= 0._dp) then

        Cl(snl+1) = (Cl(snl+1) + Cice * Wsno / dz(snl+1))

        !FLAG not sure if these last 2 are necessary, I am adding them to cover my bases- JM 15.06.2010
        !I put this back to snl+1, it was l, from cvs version 1.9 though it does not appear to make sense as
        ! l. JM 25.10.10

        !snow bulk density (kg m-3) eqn 6.66 CLM
        psno = (Wice(snl+1) + Wliq(snl+1)) / dz(snl+1)

        !snow thermal conductivity
        Kh(snl+1) = Kair + (7.75e-5 * psno + 1.105e-6 * psno * psno) * (Kice - Kair) !(W m-1 K-1) eqn 6.65 CLM


      end if

end subroutine soilthermalprop

!------------------------------------------------------------------------------------------------------------

subroutine surfaceheatflux(j)

! Calculate the heat flux from the atmosphere to the land surface

use arveparams,     only : dp,npft
use statevars,      only : sv,doy !doy only needed for testing.
use soilstate_vars, only : Eg_v,Eg_b,Eg,dHgdt,dHgdt_v,dHgdt_b,Hg,Hg_v,Hg_b,dhsdt,hs,dEgdt,&
                           dEgdt_v,dEgdt_b,dLgdt_b,Swg_v,Swg_b,Lwg_v,Lwg_b,dLgdt_v,dLgdt_b, &
                           Lwvd,lam,TairK
use latentsensiblefluxes, only : latentsens
use pft_state_variables, only : veg

implicit none

!arguments
integer, intent(in) :: j                        !index of the current gridcell

!pointer
real(dp), pointer :: fpc_sum                    !instantaneous cumulative fpc_grid over all pfts for gridcell
real(dp), pointer, dimension(:) :: fpc_grid    !instantaneous gridcell foliar projective cover (FPC) for gridcell (m2 m-2)

!local variables
real(dp) :: Hg_v_tot            !sensible heat flux into the ground (vegetated surfaces) (W m-2)
real(dp) :: Eg_v_tot            !latent heat flux into the ground (vegetated surfaces) (W m-2)
real(dp) :: Swg_tot             !shortwave flux into the ground (vegetated surfaces) (W m-2)
real(dp) :: Lwg_tot             !longwave flux into the ground (vegetated surfaces) (W m-2)
real(dp) :: dHgdt_v_tot         !derivative of the sensible heat flux into the ground w.r.t. temp (vegetated surfaces) (W m-2)
real(dp) :: dEgdt_v_tot         !derivative of the latent heat flux into the ground w.r.t. temp (vegetated surfaces) (W m-2)
real(dp) :: dLgdt_v_tot         !derivative of the longwave flux into the ground w.r.t. temp (vegetated surfaces) (W m-2)
real(dp) :: dLgdt               !derivative of the longwave flux into the ground w.r.t. temp (all surfaces) (W m-2)

!point pointers
fpc_sum         => veg%fpc_sum
fpc_grid        => veg%fpc_grid

!-------
! Begin calcs

! Use the vegetation/soil temperature to calculate the sensible and latent heat fluxes (mm -s) for the timestep
        call latentsens(j)

! Find the weighted Sw,Lw,Hg and Eg values.
  Eg_v_tot         = sum(Eg_v(1:npft) * fpc_grid(1:npft))
  Hg_v_tot         = sum(Hg_v(1:npft) * fpc_grid(1:npft))
  dHgdt_v_tot         = sum(dHgdt_v(1:npft) * fpc_grid(1:npft))
  dEgdt_v_tot         = sum(dEgdt_v(1:npft) * fpc_grid(1:npft))
  Swg_tot         = sum(Swg_v(1:npft) * fpc_grid(1:npft))
  Lwg_tot         = sum(Lwg_v(1:npft) * fpc_grid(1:npft))
  dLgdt_v_tot        = sum(dLgdt_v(1:npft) * fpc_grid(1:npft))

! Add the bare grounds fluxes (scaled by the amount of ground that is bare)
  Swg_tot = Swg_tot + Swg_b * (1._dp - fpc_sum)
  Lwg_tot = Lwg_tot + Lwg_b * (1._dp - fpc_sum)
  Hg = Hg_v_tot + Hg_b * (1._dp - fpc_sum)
  Eg = Eg_v_tot + Eg_b * (1._dp - fpc_sum)
  dHgdt = dHgdt_v_tot + dHgdt_b * (1._dp - fpc_sum)
  dEgdt = dEgdt_v_tot + dEgdt_b * (1._dp - fpc_sum)
  dLgdt = dLgdt_v_tot + dLgdt_b * (1._dp - fpc_sum)

! The total ground heat flux is given by
  hs = Swg_tot - Lwg_tot - Hg - lam * Eg  !5.129
  dhsdT =  - dLgdT - dHgdT - lam * dEgdT

!write(0,'(a,i6,a,3f12.4)')'DAY: ',doy,'   ',Swg_tot,Swg_b,fpc_sum
!write(0,'(a,8f15.2)')'hs',hs,dhsdt,Swg_tot,Lwg_tot,Hg,Eg*lam,sv(j)%Tsoil(1),TairK
!write(0,'(a,6f15.2)')'dhs',dHgdt_v_tot,dEgdt_v_tot ,dLgdt_v_tot,dHgdt,dEgdt,dLgdt 
!write(*,*)

end subroutine surfaceheatflux

!====================================================================================================

subroutine soiltemperature(j)

use arveparams,     only : dp,ns,nl,Lf,Tstune,CNfac,Tfreeze,npft
use statevars,      only : sv,dt
use soilstate_vars, only : hs,dhsdt,Eg_b,Eg_v,dEgdT_v,dEgdT_b,Lwg_b,Lwg_v,Lwv,Latm,Hg_v,Hg_b, &
                           Kh,Kl,Cl,fact,Fhti,Swg_b,                 &
                           dHgdT_v,Eg,dHgdT_b,tairk,lam,surf
use metvarsmod,     only : met
use increment,      only : fast_increment
use tridiagonal,    only : tridiag
use pft_state_variables, only : veg

implicit none

!arguments
integer, intent(in) :: j                        !index of the current gridcell

!pointers
integer,  pointer :: snl                        !index value of top snow layer (negative is more layers)
real(dp), pointer :: fpc_sum                    !instantaneous cumulative fpc_grid over all pfts for gridcell

real(dp), pointer, dimension(:)   :: dz         !thickness of the soil layers (m)
real(dp), pointer, dimension(:)   :: zpos       !z coordinate of the middle of the soil layer relative to soil surface (m), positive downwards
real(dp), pointer, dimension(:)   :: zipos      !z coordinate of the interface between soil layers relative to soil surface (m), positive downwards
real(dp), pointer, dimension(:)   :: Wliq       !soil liquid water content at layer midpoint (mm)
real(dp), pointer, dimension(:)   :: Wice       !soil ice content at layer midpoint (mm)
real(dp), pointer, dimension(:)   :: Tsoil      !soil temperature (K)
real(dp), pointer, dimension(:)   :: Tnsoi      !updated soil temperature at timestep n+1
real(dp), pointer, dimension(:)   :: Tsoiln     !soil temperature (K) from previous timestep
real(dp), pointer, dimension(:)   :: fpc_grid  !instantaneous gridcell foliar projective cover (FPC) for gridcell (m2 m-2)
real(dp), pointer, dimension(:,:) :: temp_veg   !temperature of the vegetation (K) (for the previous (1) and present (2) timesteps)
real(dp), pointer, dimension(:,:) :: Tveg_ave   !timestep average vegetation temperature (K)(for the previous (1) and present (2) timesteps)
real(dp), pointer :: qseva                      !soil evaporation flux (mm s-1)     
real(dp), pointer :: qsubl                      !soil sublimation flux (mm s-1)     
real(dp), pointer :: qsdew                      !soil surface dew (mm s-1)          
real(dp), pointer :: qfrost                     !soil surface frost (mm s-1)        

!local variables
real(dp) :: dzm                                 !distance between layer nodes upwards (m)
real(dp) :: dzp                                 !distance between layer nodes downwards (m)
real(dp) :: fevap                               !soil evaporation flux?
real(dp) :: Egn_b                               !updated surface latent heat flux for bare ground
real(dp) :: evapmax                             !max. evaporation which soil can provide at one time step
real(dp) :: soil_evap                           !total evaporative flux
real(dp) :: Eg_v_tot                            !latent heat flux across all pfts (W m-2)
real(dp) :: delZtop                             !recalc top layer thickness for more accurate surface temp (m)

real(dp), dimension(ns:nl) :: avect             !vectors for tridiagonal solver
real(dp), dimension(ns:nl) :: bvect
real(dp), dimension(ns:nl) :: cvect
real(dp), dimension(ns:nl) :: rvect

real(dp), dimension(npft) :: Egn_v              !updated surface latent heat flux for vegetated areas

integer :: a                                    !upper and lower vector boundaries
integer :: b                                    !for tridiagonal solver
integer :: l,pft                                !counters

!point pointers

Tnsoi     => sv(j)%Tnsoi
Tsoil     => sv(j)%Tsoil
Tsoiln    => sv(j)%Tsoiln
Tveg_ave  => sv(j)%Tveg_ave
Wice      => sv(j)%Wice
Wliq      => sv(j)%Wliq
dz        => sv(j)%dz
snl       => sv(j)%snl
temp_veg  => sv(j)%temp_veg
zipos     => sv(j)%zipos
zpos      => sv(j)%zpos
fpc_sum   => veg%fpc_sum
fpc_grid  => veg%fpc_grid
qseva     => surf%qseva
qsubl     => surf%qsubl
qsdew     => surf%qsdew
qfrost    => surf%qfrost

!----------------------------------------------
! Conductivity and heat flux at layer interfaces

!top snow/soil layer
l = snl+1

! Top snow/soil layer is the layer-averaged temp, so to be more accurate
! we adjust the heat capacity of the top layer. This will give a more realistic top soil layer
! temperature. See P.107 of CLM 4.0. Eqn 6.29

delZtop = 0.5_dp * (zpos(l) - zipos(l-1) + Tstune * (zpos(l+1) - zipos(l-1)))

fact(l) = dt / Cl(l) * dz(l) / delZtop
Fhti(l) = Kl(l) * (Tsoil(l+1) - Tsoil(l)) / (zpos(l+1) - zpos(l))

! snow/soil layers except bottom

do l = snl+2,nl-1
  fact(l) = dt / Cl(l)
  Fhti(l) =-Kl(l) * (Tsoil(l) - Tsoil(l+1)) / (zpos(l+1) - zpos(l))
end do

! bottom soil layer

fact(nl) = dt / Cl(nl)
Fhti(nl) = 0._dp

! set up the coefficients for the tridiagonal matrix solver

! top snow/soil layer
dzp = zpos(snl+2) - zpos(snl+1)

avect(snl+1) = 0._dp
bvect(snl+1) = 1._dp + (1._dp - CNfac) * fact(snl+1) * (Kl(snl+1) / dzp) - fact(snl+1) * dhsdT
cvect(snl+1) =       - (1._dp - CNfac) * fact(snl+1) *  Kl(snl+1) / dzp
rvect(snl+1) = Tsoil(snl+1) + fact(snl+1) * (hs - dhsdT * Tsoil(snl+1) + CNfac * Fhti(snl+1))

! snow/soil layers except bottom

do l = snl+2,nl-1

  dzm = zpos(l) - zpos(l-1)
  dzp = zpos(l+1) - zpos(l)

  avect(l) =       - (1._dp - CNfac) * fact(l) *  Kl(l-1) / dzm
  bvect(l) = 1._dp + (1._dp - CNfac) * fact(l) * (Kl(l)   / dzp + Kl(l-1) / dzm)
  cvect(l) =       - (1._dp - CNfac) * fact(l) *  Kl(l)   / dzp
  rvect(l) = Tsoil(l) + CNfac * fact(l) * (Fhti(l) - Fhti(l-1))

end do

! bottom soil layer
l = nl
dzm = zpos(l) - zpos(l-1)

avect(l) =       - (1._dp - CNfac) * fact(l) *  Kl(l-1) / dzm
bvect(l) = 1._dp + (1._dp - CNfac) * fact(l) *  Kl(l-1) / dzm
cvect(l) = 0._dp
rvect(l) = Tsoil(l) - CNfac * fact(l) * Fhti(l-1)

! Call tridiagonal solver

a = snl+1
b = nl

call tridiag(avect(a:b),bvect(a:b),cvect(a:b),rvect(a:b),Tnsoi(a:b))

! Perform phase change

call soilphasechg(j)

! Update the sensible and latent heat surface fluxes and calculate surface evaporation
! CLM Tech note 3.0 section 5.4

!ground under vegetation

Hg_v(:) = Hg_v(:) + (Tnsoi(snl+1) - Tsoil(snl+1)) * dHgdT_v(:)   !Eqn. 5.120
Eg_v(:) = Eg_v(:) + (Tnsoi(snl+1) - Tsoil(snl+1)) * dEgdT_v(:)   !Eqn. 5.121

!bare ground (no canopy)

Hg_b = Hg_b + (Tnsoi(snl+1) - Tsoil(snl+1)) * dHgdT_b   !Eqn. 5.120
Eg_b = Eg_b + (Tnsoi(snl+1) - Tsoil(snl+1)) * dEgdT_b   !Eqn. 5.121


! Adjust the ground evaporation for the VEGETATED areas

soil_evap = sum(max(0._dp,Eg_v(1:npft) * fpc_grid(1:npft)))  !sum only the areas where it is evaporating

if (soil_evap == 0._dp) then !if no evaporation

  Egn_v(:) = Eg_v(:)

else !evaporation is occuring

  !find maximum soil evaporation for this timestep

  evapmax = (Wice(snl+1) + Wliq(snl+1)) / dt

  !NOTE: This was previously set up wrong. It allowed fevap (which is a fraction) to be always 1 or greater!
  ! it is now changed to coincide with tech notes. JM 08.06.2010

  if (soil_evap > 0._dp) then
    fevap = min(evapmax / soil_evap, 1._dp)   !Eqn. 5.122   
  else
    fevap =  1._dp
  end if

  if (fevap < 1._dp) then
    do pft = 1,npft
      Egn_v(pft) = fevap * Eg_v(pft)                    !Eqn. 5.123

      !Give the energy defecit to sensible heat
      Hg_v(pft) = Hg_v(pft) + lam * (Eg_v(pft) - Egn_v(pft))          !Eqn. 5.124
    end do
  else
    Egn_v(:) = Eg_v(:)
  end if

end if

!adjust the ground evaporation for the BARE ground areas

soil_evap = Eg_b * (1._dp - fpc_sum)

if (soil_evap == 0._dp) then !if no evaporation

  Egn_b = Eg_b

else !evaporation is occuring

  !find maximum soil evaporation for this timestep
  evapmax = (Wice(snl+1) + Wliq(snl+1)) / dt

  !NOTE: This was set up wrong. It allowed fevap (which is a fraction) to be always 1 or greater!
  !it is now changed to coincide with tech notes. JM 08.06.2010

  if (soil_evap > 0._dp) then
    fevap = min(evapmax / soil_evap, 1._dp)   !Eqn. 5.122
  else
   fevap =  1._dp
  end if

  if (fevap < 1._dp) then
    Egn_b = fevap * Eg_b                !Eqn. 5.123
    Hg_b = Hg_b + lam * (Eg_b - Egn_b)  !Eqn. 5.124
  else
    Egn_b = Eg_b
  end if

end if

!Now find the gridcell Eg and Hg to determine the water fluxes

Eg_v_tot = sum(Egn_v(1:npft) * fpc_grid(1:npft))

Eg = Egn_b * (1._dp - fpc_sum) + Eg_v_tot

!Update the evaporation, sublimation and frost fluxes for soil hydrology
if (Eg > 0._dp) then !evap.

  if (Wliq(snl+1) + Wice(snl+1) > 0._dp) then
    qseva  = max(Eg * Wliq(snl+1) / (Wliq(snl+1) + Wice(snl+1)) ,0._dp)  !Eqn. 5.125
  else
    qseva  = 0._dp  !Eqn. 5.125
  end if

  qsubl  = max(Eg - qseva, 0._dp)                    !Eqn. 5.126

  qsdew  = 0._dp
  qfrost = 0._dp

else !condens.

  if (Tsoil(snl+1) >= Tfreeze) then
    qsdew = abs(Eg)                      !Eqn. 5.127
    qfrost = 0._dp
  else
    qfrost = abs(Eg)                     !Eqn. 5.128
    qsdew = 0._dp
  end if

  qseva = 0._dp
  qsubl = 0._dp

end if

!write(0,'(5es12.4)')Eg,qseva,qsubl,qsdew,qfrost
!save updated temps

Tsoiln(snl+1:nl) = Tsoil(snl+1:nl) !save previous temp for radflux
Tsoil(snl+1:nl) = Tnsoi(snl+1:nl)  !update to new temp.
Tveg_ave(2,:) = Tveg_ave(2,:) + temp_veg(2,:)  !increment the average veg tracker

end subroutine soiltemperature

end module soiltemperaturemod
