module soilinit

! Set up the soil layers and depth. Allocates soil hydrology variables. Also assigns
! the soil properties based upon the soil physical characteristics (sand, silt, etc.)

implicit none

public :: soillayers
public :: soilprep

contains

!-----------------------

subroutine soillayers(j)

! Set up the soil layers and find the bottom of the soil column. This assigns gnl.
! Created by JM Aug 2010

use arveparams, only : dp,ns,nl
use statevars, only : sv
use iovariables, only : soil
use soilstate_vars, only : surf

implicit none

!arguments
integer, intent(in) :: j  !gridcell index

!pointers
integer, pointer  :: gnl                   !index value of lowest 'soil' layer
real(dp), pointer :: soildepth             !depth to bedrock (m)
real(dp), pointer, dimension(:) :: dz      !thickness of the soil layers (m)
real(dp), pointer, dimension(:) :: zpos    !z coordinate of the middle of the soil layer relative to soil surface (m), positive downwards
real(dp), pointer, dimension(:) :: zipos   !z coordinate of the interface between soil layers relative to soil surface (m), positive downwards
real(dp), pointer, dimension(:) :: dzmm    !thickness of the soil layers (mm)
real(dp), pointer, dimension(:) :: zposmm  !z coordinate of the middle of the soil layer relative to soil surface (mm), positive downwards

!parameters
real(dp), parameter :: fs = 0.0225_dp                                              !scaling factor for soil depth calculation
real(dp), dimension(5), parameter :: depth_lookup = [ 0.2, 0.4, 0.6, 0.8, 1.0 ]    ! FAO soil depths

!local variables
integer :: l
integer :: nln                              !normal calculated layers
integer :: depth_int                        !layer of the lowest dataset soil properties

!point pointers
dz     => sv(j)%dz
zpos   => sv(j)%zpos
zipos  => sv(j)%zipos
soildepth => sv(j)%soildepth
dzmm   => sv(j)%dzmm
zposmm => sv(j)%zposmm
gnl    => sv(j)%gnl

!----------------------------------------

!set the initial soil layer thicknesses (layer 1 is the depth of the snow)
!node (layer center depths, 0 is soil surface, positive downwards)

        !layer centre depths

        !nl = 17
        !nl(1:15) normal
        !nl(16) 30m slab
        !nl(17) 1m bottom layer

        nln = nl - 2

        forall (l = 1:nln)
          zpos(l) = fs * (exp(0.5_dp * (l - 0.5_dp)) - 1._dp)
        end forall

        !layer thicknesses

        dz(1) = 0.5_dp * (zpos(1) + zpos(2))

        forall (l = 2:nln-1)
          dz(l) = 0.5_dp * (zpos(l+1) - zpos(l-1))
        end forall

        dz(nln) = zpos(nln) - zpos(nln-1)

        !interface depths

        zipos(0) = 0._dp
        forall (l = 1:nln-1)
          zipos(l) = 0.5_dp * (zpos(l) + zpos(l+1))
        end forall

        zipos(nln) = zpos(nln) + 0.5_dp * dz(nln)

        !lowest 2 layers 
        !thickness is fixed and everything else is calculated from that

        !thickness
        dz(nl-1) = 30._dp !m
        dz(nl)   =  1._dp !m

        !layer midpoints

        zpos(nl-1) = zipos(nl-2) + 0.5_dp * dz(nl-1)
        zpos(nl) = zpos(nl-1) + 0.5_dp * (dz(nl-1) + dz(nl))

        !layer interfaces

        zipos(nl-1) = zpos(nl-1) + 0.5_dp * dz(nl-1)
        zipos(nl)   = zpos(nl)   + 0.5_dp   * dz(nl)

        dzmm   = 1000._dp * dz(1:nl)
        zposmm(1:nl) = 1000._dp * zpos(1:nl)
        zposmm(nl+1) = zposmm(nl)  !set the 'aquifer' thickness to be 0, ie. the depth of the bottom layer.

! Assign soil depth.

   depth_int = 5  !this is set to the lowest level in the FAO soil dataset

   do l = 1,5
    if (soil(j)%sdto(l) < 0._dp) then !if the dataset has negative value then no soil in this location
       depth_int = l-1
       exit
     end if
   end do

   if (depth_int < 4) then  !check if the soil depth is lower than the lowest soil dataset point
   soildepth = depth_lookup(depth_int)
    else
   soildepth = max(depth_lookup(depth_int), soil(j)%maxd * 0.01)
   end if

! For areas with more than 250cm of soil (the deepest that is known by the FAO
! soil depth dataset) we assign a correction that assumes the soil is as deep as our total soil column.
! the properties of the bedrock will be set in the next section.

!FLAG, this might not be appropriate, JM 21.04.2011
!   if (soildepth > 1._dp) then  !double the depth
!     soildepth = soil(j)%maxd * 0.01 * 2.
!   else if (soildepth > 2.0) then !make it as deep as the total soil column
!     soildepth = zipos(nl)
!   end if

 ! if (soildepth > 2.0) then !make it as deep as the total soil column
 !    soildepth = zipos(nl)
 ! end if

! Set the index of the ground level that will be the bottom of 'soil'. This will be used in hydrology calculations
! while the total ground column will be used for temperature calculations.
    gnl = nl
    do l = 2,nl
        if(soildepth <= zipos(l)) then !zipos is the interface depth of the bottom of the soil column.
           gnl = l
           exit
        end if
    end do

end subroutine soillayers

!--------------------------------------------------------------

subroutine soilprep(j)

!coded from a CLM 3.0 base code with relatively major changes by Jed Kaplan and Joe Melton, 2007

use arveparams, only : dp,ns,nl,Tfreeze,pliq,pice,OMorgC,peatlim,Kair,omd,sandd,siltd,clayd,rockd, &
                        sandquartz,rockquartz,bedrock_tc,Corg,Csand,Cclay,Csilt,Crock,Kom,Kquartz,Kmineral,Kliq, &
                        Cice,Cliq,Kom_dry
use statevars, only : sv
use soilparameters, only : fKsat, fTsat, fBexp, fPsat, fbulk, fRock!fVcf fKsato, fTsato, fBexpo, fPsato
use iovariables, only : soil
use organicsoilprop, only : fKsato, fTsato, fBexpo, fPsato
use soilstate_vars, only : surf

implicit none

!arguments
integer, intent(in) :: j  !gridcell index

!pointers
integer, pointer :: snl                 !index value of top snow layer (negative is more layers)
integer, pointer  :: gnl                !index value of lowest 'soil' layer
real(dp), pointer :: Wsno               !snow water equivalent of the snowpack (mm)
real(dp), pointer :: zsno               !snow depth (m)
real(dp), pointer :: soildepth          !depth to bedrock (m)
real(dp), pointer, dimension(:) :: dz      !thickness of the soil layers (m)
real(dp), pointer, dimension(:) :: zpos    !z coordinate of the middle of the soil layer relative to soil surface (m), positive downwards
real(dp), pointer, dimension(:) :: zipos   !z coordinate of the interface between soil layers relative to soil surface (m), positive downwards
real(dp), pointer, dimension(:) :: dzmm    !thickness of the soil layers (mm)
real(dp), pointer, dimension(:) :: zposmm  !z coordinate of the middle of the soil layer relative to soil surface (mm), positive downwards
real(dp), pointer, dimension(:) :: sand    !soil sand content (mass percent)
real(dp), pointer, dimension(:) :: silt    !soil silt content (mass percent)
real(dp), pointer, dimension(:) :: clay    !soil clay content (mass percent)
real(dp), pointer, dimension(:) :: sorg    !soil organic matter content (mass percent)
real(dp), pointer, dimension(:) :: bulk    !soil bulk density (kg m-3)
real(dp), pointer, dimension(:) :: rock    !mass fraction of coarse fragments in soil
real(dp), pointer, dimension(:) :: Ksat    !soil water saturated conductivity at layer midpoint (mm s-1)
real(dp), pointer, dimension(:) :: Tsat    !soil water volumetric water content at saturation (fraction)
real(dp), pointer, dimension(:) :: Bexp    !soil water B exponent used in the Brooks & Corey Pedotransfer Functions
real(dp), pointer, dimension(:) :: Psat    !soil water matric potential at saturation (mm)
real(dp), pointer, dimension(:) :: Wliq    !soil liquid water content at layer midpoint (mm)
real(dp), pointer, dimension(:) :: Wice    !soil ice content at layer midpoint (mm)
real(dp), pointer, dimension(:) :: Tsoil   !soil temperature (K)
real(dp), pointer, dimension(:) :: Tsoiln  !soil temperature for previous timestep (K)
real(dp), pointer, dimension(:) :: Tliq    !fractional soil water content
real(dp), pointer, dimension(:) :: Tice    !fractional soil ice content
real(dp), pointer, dimension(:) :: Psi     !soil water potential
real(dp), pointer, dimension(:) :: Vcf     !FRACTION of coarse fragments (rocks) by volume
real(dp), pointer, dimension(:) :: Vom     !volumetric fraction of soil organic matter
real(dp), pointer, dimension(:) :: Vsand   !volumetric fraction of sand
real(dp), pointer, dimension(:) :: Vclay   !volumetric fraction of clay
real(dp), pointer, dimension(:) :: Vsilt   !volumetric fraction of silt
real(dp), pointer, dimension(:) :: Kdry    !soil thermal conductivity of dry natural soil (W m-1 K-1)
real(dp), pointer, dimension(:) :: Ksolid  !soil mineral thermal conductivity at layer midpoint (W m-1 K-1)
real(dp), pointer, dimension(:) :: Csolid  !soil solids volumetric heat capacity at layer midpoint (J m-3 K-1)
real(dp), pointer :: Waquif_a              !unconfined aquifer water content (mm)
real(dp), pointer :: Waquif_t              !total aquifer water content (mm)
real(dp), pointer :: zw                    !height of water table (m)
real(dp), pointer :: theta_fc              !the volumetric water content at field capacity (fraction)
logical, pointer :: peat                   !peatland flag (true if peat)

!parameters
real(dp), parameter :: a = 0.053_dp        ! From Balland and Arp for dry soil thermal conductivity calc
real(dp), parameter :: alpha = 0.24_dp     ! +/- 0.04 constant for Knum calc.from Balland and Arp
real(dp), parameter :: beta = 18.3_dp      ! +/- 1.1 constant for Knum calc.from Balland and Arp

!local variables
integer :: l,c,x
real(dp) :: densp                         !particle density
real(dp) :: Vtotal                        !total volume
real(dp), dimension(1:nl) :: Kdrysolid    !soil dry mineral thermal conductivity at layer midpoint (W m-1 K-1)

!point pointers
snl                    => sv(j)%snl
gnl                    => sv(j)%gnl
dz                     => sv(j)%dz
peat                   => sv(j)%peat
zpos                   => sv(j)%zpos
zipos                  => sv(j)%zipos
soildepth              => sv(j)%soildepth
dzmm                   => sv(j)%dzmm
zposmm                 => sv(j)%zposmm
sand                   => sv(j)%sand
silt                   => sv(j)%silt
clay                   => sv(j)%clay
sorg                   => sv(j)%sorg
bulk                   => sv(j)%bulk
rock                   => sv(j)%rock
Ksat                   => sv(j)%Ksat
Tsat                   => sv(j)%Tsat
Bexp                   => sv(j)%Bexp
Psat                   => sv(j)%Psat
Wice                   => sv(j)%Wice
Wliq                   => sv(j)%Wliq
Wsno                   => sv(j)%Wsno
zsno                   => sv(j)%zsno
Tsoil                  => sv(j)%Tsoil
Tsoiln                 => sv(j)%Tsoiln
Tliq                   => sv(j)%Tliq
Tice                   => sv(j)%Tice
Psi                    => sv(j)%Psi
Waquif_a               => sv(j)%Waquif_a
Waquif_t               => sv(j)%Waquif_t
zw                     => sv(j)%zw     
Vcf                    => sv(j)%Vcf    
Vom                    => sv(j)%Vom    
Vsand                  => sv(j)%Vsand
Vclay                  => sv(j)%Vclay
Vsilt                  => sv(j)%Vsilt
Kdry                   => sv(j)%Kdry
Csolid                 => sv(j)%Csolid
Ksolid                 => sv(j)%Ksolid
theta_fc               => sv(j)%theta_fc

!----------------

! Assign the FAO soil data into the appropriate ARVE soil layers. If using a
! different soil dataset, this will require changes!

!soil layers FAO soil data are for depth midpoints, 10, 30, 50, 70, 90 cm

!FAO 0 - 20cm, soil model to 0 - 14.9cm

!do for soil layers 1 to 4, as long as 4 is less than gnl, otherwise 1 to gnl
x = min(4,gnl)

do l = 1,x

    sand(l) = soil(j)%sdto(1)
    silt(l) = soil(j)%stpc(1)
    clay(l) = soil(j)%clpc(1)
    sorg(l) = soil(j)%totc(1)
    sorg(l) =(sorg(l) / OMorgC) * 0.1
    sorg(l) = min(100._dp,sorg(l))
    Vcf(l) = soil(j)%cfrag(1) * 0.01
    bulk(l) = soil(j)%bulk(1) * 1000.

end do

  c = 2 !now for FAO level 2.

!do for soil layers 5 to 8, as long as 4 is less than gnl, otherwise 5 to gnl
x = min(8,gnl)

if (x >= 5) then

do l = 5,x

  !FAO(2) 20 - 40cm, soil model(5) 14.9 - 26cm
  !FAO(3) 40 - 60cm, soil model(6) 26 -44.4cm
  !FAO(4) 60 - 80cm, soil model(7) 44.4 - 74.6cm
  !FAO(5) 80 - 100cm, soil model(8) 74.6 - 124.45

  if (soil(j)%sdto(c) > 0._dp) then
    sand(l) = soil(j)%sdto(c)
    silt(l) = soil(j)%stpc(c)
    clay(l) = soil(j)%clpc(c)
    sorg(l) = soil(j)%totc(c)
    sorg(l) = (sorg(l) / OMorgC) * 0.1
    sorg(l) = min(100._dp,sorg(l))
    Vcf(l) = soil(j)%cfrag(c) * 0.01
    bulk(l) = soil(j)%bulk(c) * 1000.
  else !no soil values at this depth, use the layer above
    sand(l) = sand(l - 1)
    silt(l) = silt(l - 1)
    clay(l) = clay(l - 1)
    sorg(l) = sorg(l - 1)
    Vcf(l) = Vcf(l - 1)
    bulk(l) = bulk(l - 1)
  end if

  c = c + 1  !increment the FAO depth level.

 end do
end if

!using FAO 80 - 100cm, soil model 124.5 - bottom
! NOTE: This assumes that the soil does not change composition with depth
! beyond that available in the FAO datset.

!do for soil layers 9 to nl, as long as nl is less than gnl, otherwise 9 to gnl 
x = min(nl,gnl)

if (x >= 9) then
 do l = 9,x

  sand(l) = soil(j)%sdto(5)
  silt(l) = soil(j)%stpc(5)
  clay(l) = soil(j)%clpc(5)
  sorg(l) = soil(j)%totc(5)
  sorg(l) = (sorg(l) / OMorgC) * 0.1
  sorg(l) = min(100._dp,sorg(l))
  Vcf(l) = soil(j)%cfrag(5) * 0.01
  bulk(l) = soil(j)%bulk(5) * 1000.

 end do
end if

!----

! Water holding capacity and saturated conductivity
! soil moisture to 30% of saturation, no ice, no snow

do l = 1,gnl

  bulk(l) = fbulk(zipos(l),sorg(l),Vcf(l),clay(l),silt(l))

  rock(l) = fRock(Vcf(l),bulk(l)) !find the mass fraction of coarse fragments Skirvin correction to WEPP model

  !NOTE: Reduce the sand,silt,clay mass percent by the coarse fragment mass percent (while retaining the relative weight percents
  ! since sand,clay,silt add up to 100 regardless of Vcf or sorg in the FAO dataset. -JM Apr 04 08
  !FLAG is this correct to change the values for coarse frags???????!
 ! sand(l) = sand(l) * (1. - rock(l))
 ! silt(l) = silt(l) * (1. - rock(l))
 ! clay(l) = clay(l) * (1. - rock(l))

  if (rock(l) /= 0._dp) then !if there is rock in the soil then redo the bulk (since sand,clay,silt have changed)-JM Apr 20 08
    bulk(l) = fbulk(zipos(l),sorg(l),Vcf(l),clay(l),silt(l))
  end if

    if (sorg(l) < peatlim) then ! peatlim is the percent organic cut-off for org soil (peat) in Hillier

        !for mostly mineral soil
          Ksat(l) = fKsat(sand(l),sorg(l))
          Tsat(l) = fTsat(sand(l),clay(l),sorg(l),Vcf(l),bulk(l),zpos(l)) !note zpos since it needs the mid point
          Bexp(l) = fBexp(clay(l),silt(l),sand(l),rock(l),sorg(l))
          Psat(l) = fPsat(sand(l),rock(l))

          if (Vcf(l) > 0.) then !from Brakensiek & Rawls, Catena 1994

           !Changed from Peck and Watson's formulation, as it overpredicts conductivities, esp at higher gravel fractions
           !see Bouwer & Rice 1984 p. 698. -JM Apr 16 08
             !old P & W formula: Ksat(l) = Ksat(l) * (2 - 2 * Vcf(l)) * (1 / (2 + (Vcf(l))))

            !Brakensiek & Rawls, Catena 1994
            Ksat(l) = Ksat(l) * ((1._dp - rock(l)) / (1._dp - (rock(l) * 0.25)))

          end if

      else !'peat soil'

         ! flag the soil as peatland
         peat = .true.

         ! Ksat(l) = fKsato(bulk(l))
         ! Tsat(l) = fTsat(sand(l),clay(l),sorg(l),Vcf(l),bulk(l),zpos(l),rock(l)) !note zpos since it needs the mid point
         ! Bexp(l) = fBexpo(bulk(l))
         ! Psat(l) = fPsato(bulk(l),sand(l),rock(l),sorg(l),Ksat(l))

        !New way from CLASS (Letts 2000) (Coded by JOK 13.07.2010) FLAG CHECK!
          Ksat(l) = fKsato(zpos(l))
          Tsat(l) = fTsato(zpos(l))
          Bexp(l) = fBexpo(zpos(l))
          Psat(l) = fPsato(zpos(l))

    end if
end do

   ! If soil column does not extend completely to nl, set lower layers to be bedrock.
   ! if gnl == nl, then assume that the soil column extends at least beyond zipos(nl), otherwise
   ! set the rest of the ground column to bedrock (deeper than gnl)
           if (gnl < nl) then
               do l = min(gnl+1,nl),nl   !layers between gnl+1 and nl

                   Ksat(l) =    1.e-3  !(mm s-1)Value taken from Katsura et al 2006 Vadose Zone Journal for weathered granitic
                   Tsat(l) =    0.073  !taken from paper above (average of the data listed)
                   Bexp(l) =   -0.1_dp  !value estimated from curves in Kew and Gilkes 2004 Supersoil
                   Psat(l) = -398._dp   !mm value taken from Katsura et al 2006

                   !for all layers below the bedrock depth set to 100% coarse fragments and no sand, silt, clay or org matter
                   sand(l) = 0._dp
                   silt(l) = 0._dp
                   clay(l) = 0._dp
                   Vcf(l) = 1._dp
                   rock(l) = 1._dp
                   sorg(l) = 0._dp
                   bulk(l) = 2650._dp

               end do
           end if

 !hydrology initial pre-sets for first run.
    do l = 1,nl 

      Tliq(l) = 0.3_dp * Tsat(l)
      Wliq(l) = Tliq(l) * dz(l) * pliq
      Tice(l) = 0._dp ! assume no ice.
      Wice(l) = Tice(l) * dz(l) * pice

      Psi(l) = Psat(l) * (Tliq(l) / Tsat(l))**(-Bexp(l))

      if (Psi(l) < -1.e8) then
        Psi(l) = -1.e8
      end if

    end do

!calculate dry soil thermal conductivity and heat capacity.
!also initial soil temp presets.

do l = 1,nl

           !inital presets
           Tsoil(l) = Tfreeze + 5._dp !set to 5 deg. C.
           Tsoiln(l) = Tfreeze + 5.05_dp

           !find the volume fractions of the soil consituents.
           Vom(l)   = sorg(l) * 0.01 / omd
           Vsand(l) = sand(l) * 0.01 / sandd
           Vclay(l) = clay(l) * 0.01 / clayd
           Vsilt(l) = silt(l) * 0.01 / siltd

           Vtotal = Vom(l) + Vsand(l) + Vclay(l) + Vcf(l) + Vsilt(l)

           Vom(l)   = Vom(l)   / Vtotal
           Vsand(l) = Vsand(l) / Vtotal
           Vclay(l) = Vclay(l) / Vtotal
           Vsilt(l) = Vsilt(l) / Vtotal

           Csolid(l) = (Vom(l) * (Corg) + Vsand(l) * Csand + Vclay(l) * Cclay + Vcf(l)&
                             * Crock + Vsilt(l) * Csilt)   !(units: Jm-3K-1)


           Ksolid(l) = (Kom**Vom(l)) * (Kquartz**((sandquartz * Vsand(l)) + (rockquartz * Vcf(l)))) * &
                     (Kmineral**(1 - Vom(l) - (sandquartz * Vsand(l)) - (rockquartz * Vcf(l)))) !(units: WM-1K-1) B&A Eqn 15
                     
           Kdrysolid(l) = (Kom_dry**Vom(l)) * (Kquartz**((sandquartz * Vsand(l)) + (rockquartz * Vcf(l)))) * &
                     (Kmineral**(1 - Vom(l) - (sandquartz * Vsand(l)) - (rockquartz * Vcf(l)))) !(units: WM-1K-1) B&A Eqn 15
 
           if (l <= gnl+1) then  !for soil column     

            !dry soil thermal conductivity  

            densp = 1. / (Vom(l) + (1. - sorg(l) * 0.01) / clayd) !eqn 3

            Kdry(l) = ((a * Kdrysolid(l) - Kair) * (bulk(l) * 0.001) + Kair * densp) &
             / (densp - (1. - a) * (bulk(l) * 0.001))   !eqn 16

           else  !set bedrock thermal conductivity
           	
            Kdry(l) = bedrock_tc 
           
           end if

end do 

        ! the volumetric water content at field capacity is derived by assuming a hydraulic
        ! conductivity of 0.1mm/day and inverting the hydraulic conductivity function (CLM 4 eq. 5.70)
        theta_fc = Tsat(1) * (0.1 / (86400. * Ksat(1)))**(1._dp / (2. * Bexp(1) + 3.))

!More initial pre-sets for first run

        !no snow or ice.
        snl = 0
        Wsno = 0._dp
        zsno = 0._dp

        !set initial water table depth 
        Waquif_a = 5100._dp!4800._dp
        Waquif_t = 5100._dp!4800._dp
        zw = zipos(gnl) + 1._dp

end subroutine soilprep

end module soilinit
