subroutine resistance(j)

! Calculation of resistance to sensible and latent heat transfer for
! vegetated and non-vegetated surfaces. Based upon the CLM 3.1 formulation with multiple modifications
! NOTE: there is presently no difference in the calculation of sensible vs. latent resistances!
! Added Oct 08 08 Joe Melton 2008

use arveparams, only : dp,npft,minplant,vkarm,visco,z0mg,grav,Tfreeze
use statevars, only : sv
use soilstate_vars, only : TairK,surf
use pftparametersmod, only : prm
use metvarsmod, only : atm

implicit none

!ARGUMENTS
integer, intent(in) :: j        !gridcell index

!pointers
real(dp), pointer, dimension(:) :: lai_sno              !instantaneous leaf area index
real(dp), pointer, dimension(:) :: stem_sno             !instantaneous stem area index
real(dp), pointer, dimension(:) :: H                    !height of the vegetation per gridcell (m)
real(dp), pointer :: winds                              !instantaneous 'wind speed' diurnally varying (m s-1)
real(dp), pointer, dimension(:)  :: rawp                !aerodynamic resistance for latent heat transfer btwn canopy air & ground (s m-1)
real(dp), pointer, dimension(:)  :: raw                 !aerodynamic resistance for latent heat transfer btwn canopy air & atmosphere (s m-1)
real(dp), pointer, dimension(:)  :: rah                 !aerodynamic resistance for sensible heat transfer btwn canopy air & atmosphere (s m-1)
real(dp), pointer  :: raw_b                             !aerodynamic resistance for latent heat transfer btwn bare ground & atmosphere (s m-1)
real(dp), pointer  :: rah_b                             !aerodynamic resistance for sensible heat transfer btwn bare ground & atmosphere (s m-1)
real(dp), pointer, dimension(:)  :: rahp                !aerodynamic resistance for sensible heat transfer btwn canopy air & ground (s m-1)
real(dp), pointer, dimension(:) :: rb                   !leaf boundary layer resistance (s m-1)
real(dp), pointer, dimension(:) :: Rzmom                !Ratio of momentum roughness length to canopy top height
real(dp), pointer, dimension(:) :: Tsoil                !soil temperature (K)
integer, pointer :: snl                                 !index value of top snow layer (negative is more layers)
real(dp), pointer :: ustar_b                            !friction velocity estimated from scalar wind speed (m s-1)
real(dp), pointer :: fsnow                              !fraction of the gridcell covered by snow (fraction)

!local variables
integer :: pft                  !counter
real(dp) :: zpdisp,zpdisp_b     !zero plane displacement (cm) (vegetated,bare)
real(dp) :: roughl,roughl_b     !roughness length  (vegetated,bare)
real(dp) :: mroughl             !momentum calculation roughness length (cm)
real(dp) :: slroughl            !sensible and latent roughness length for vegetated surface (cm)
real(dp) :: ram                 !aerodynamic resistance to momentum transfer (s m-1)
real(dp) :: ustar               !friction velocity estimated from scalar wind speed (m s-1)
real(dp) :: Cs                  !turbulent transfer coefficient between the underlying soil and the canopy air (m s-1/2)
real(dp) :: W                   !Weight (determinant of how dense the canopy is)
real(dp) :: Cs_bare             !turbulent transfer coefficient for bare soil (m s-1/2)
real(dp) :: z,z_b               !height above surface (cm) (vegetated,bare)
real(dp) :: x                   !temp var.
real(dp) :: rah0                !
real(dp) :: P,ned               !variables for eqns 13 and 14 Lhomme et al. 1994

!parameters
real(dp), parameter :: refH = 5._dp             !reference height (m)
real(dp), parameter :: a = 0.13_dp              !
real(dp), parameter :: b = vkarm / a            !
real(dp), parameter :: Cs_dense = 0.004_dp      !turbulent transfer coefficient for dense canopy (Dickinson et al 1993)(m s-1/2)
real(dp), parameter :: Cv = 0.01_dp             !turbulent transfer coefficient between canopy surface and canopy air (m s-1/2)
real(dp), parameter :: d_leaf = 0.04_dp         !Characteristic dimension of the leaves in the direction of wind flow (m)
real(dp), parameter :: mrsnow = 0.0024_dp       !momentum roughness length for snow covered surfaces

!point pointers
ustar_b         => surf%ustar_b
Tsoil           => sv(j)%Tsoil
H               => sv(j)%height
winds           => atm%winds
stem_sno        => sv(j)%stem_sno
lai_sno         => sv(j)%lai_sno
rawp            => surf%rawp
raw             => surf%raw
rah             => surf%rah
raw_b           => surf%raw_b
rah_b           => surf%rah_b
rahp            => surf%rahp
rb              => surf%rb
snl             => sv(j)%snl
Rzmom           => prm%Rzmom
fsnow           => surf%fsnow

!--------

!Determine the momentum roughness length for bare soil (m)

  if (-snl > 0) then
    mroughl = fsnow * mrsnow + (1._dp - fsnow) * z0mg   !snow covered
  else
    mroughl = z0mg      !bare soil
  end if

do pft = 1,npft
   if ((lai_sno(pft) + stem_sno(pft)) > minplant) then !vegetated surface

      !Computation of zero plane displacement and roughness length as a function
      !of LAI and height based on Shaw & Pereira (1982) (m)
      z = (H(pft) + refH) !5 meters above the vegetation height is chosen (as in Hydrall Model)

      !find the zero plane displacement (cm) (commonly comes out to be ~2/3 the vegetation height)
      x = 0.2_dp * lai_sno(pft)
      zpdisp = H(pft) * (log(1._dp + x**0.166) + 0.03_dp * log(1._dp + x**6))
      if (zpdisp >= H(pft)) zpdisp = 0.99 * H(pft)

      !find the roughness length (cm)
      if (x < 0.2_dp) then
        roughl = (0.01_dp + 0.28*sqrt(x) * H(pft))
      else
        roughl = 0.3_dp * H(pft) * (1._dp - (0.01_dp * zpdisp)/H(pft))
      end if

    !determine the resistances for transfers between the ground and the canopy

      !the aerodynamic resistance changes according to conditions that influence turbulence. The main variables are wind speed and
      !the roughness length. This approach is taken from Szeicz et al. 1969 Water Resources Res. (5) 380.
      !the aerodynamic resistance for momentum transfer between the canopy air and the atmosphere (s m-1)
      !NOTE: Equation is valid only in a neutral atmosphere,and when the vertical mixing of the turbulent boundary layer
      ! is enhanced or suppressed by buoyancy a correction term has to be introduced,depending on the sign and
      !size of the vertical temperature gradient.In climatological and hydrological investigations restricted to
      ! temperate climates, buoyancy corrections can safely be neglected.
      ram = (log((z - zpdisp) / roughl))**2 / (vkarm**2 * winds) !s m-1

      !magnitude of wind velocity incident on the leaves (equivalent here to friction velocity)
      ustar = winds * (1._dp / (ram * winds))**0.5  !CLM 5.99

      !Weight, determinant of how dense the canopy is
      W = exp(-(lai_sno(pft) + stem_sno(pft))) !CLM 5.101

      !turbulent transfer coefficient for bare soil
      Cs_bare = b * (mroughl * ustar / visco)**(-0.45_dp)

      !turbulent transfer coefficient between the underlying soil and the canopy air
      Cs = Cs_bare * W + Cs_dense * (1._dp - W)

      !resistance for latent and sensible heat transfer between the ground and canopy
      rahp(pft) = 1._dp / (Cs * ustar)
      rawp(pft) = rahp(pft)

    !find the resistances for transfers between the canopy air and the atmosphere

      !the roughness lengths for vegetated surfaces are asuumed the same for sensible and latent heat (m)
      slroughl = H(pft) * Rzmom(pft)

      !resistances for transfers between the canopy air and the atmosphere (Szeicz et al.)
      rah0 = (log((z - zpdisp) / slroughl))**2 / (vkarm**2 * winds) !s m-1
      raw(pft) = rah0

      !Apply an adjustment from Lhomme et al 1994
      ned = 5._dp * (z - zpdisp) * grav * (Tsoil(snl+1) - tairK) / (tairK * winds)

      !assume if (1 + ned) is less than zero that is stable conditions (use P =2) otherwise unstable.
      if ((1._dp + ned) > 0._dp) then
        P = 0.75
       else
        P = 2
      end if

      rah(pft) = rah0 / (max(0.01_dp,(1._dp + ned)))**P !min constraint added for numerical stability JM 26.10.10

      !leaf boundary layer resistance
      rb(pft) = 1._dp / Cv * (ustar / d_leaf)**(-0.5)

   end if
end do

!Find the resistances for bare ground

  z_b = refH   !2m in m
  zpdisp_b = 0._dp !zero plane displacement (m)

  ustar_b = 0.14 * winds !friction velocity estimated from the scalar wind speed (Weber 1999 Boundary Layer Meteor, 93, 197)

  !convert from momentum roughness length to sensible and latent roughness length. This is because the transfer
  !of momentum is affected by pressure fluctuations in the turbulent waves behind the roughness elements, while for
  !heat and water vapour transfer no such dynamical mechanism exists.(m)

  roughl_b = mroughl * exp(-a * (ustar_b * mroughl / visco)**0.45)  !CLM 5.67

  !the aerodynamic resistance for momentum transfer between the ground and the atmosphere (s m-1)
  !from Szeicz et al. 1969 Water Resources Res. (5) 380. No difference is accounted for between
  !sensible and latent heat transfers
  rah_b = (log((z_b - zpdisp_b) / roughl_b))**2 / (vkarm**2 * winds) !s m-1
  raw_b = rah_b

!write(*,'(a,4es12.4,2f10.3)')'resist',raw,rah,rawp,rahp,rah_b,raw_b

end subroutine resistance
