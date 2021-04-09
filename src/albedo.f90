subroutine albedo(j)

!calculates the soil, vegetation and surface albedos. Also does the solar flux absorbed by
!the canopy for diffuse and direct fluxes. Coded by JM May 8 08
!called after new snow flux but before soil heat flux routines.

use arveparams,  only : dp,Tfreeze,npft,d2r,pi,z0mg,band
use pftparametersmod, only : prm
use statevars,  only : sv
use soilstate_vars,  only : asurf_g,surf
use pft_state_variables, only : veg

implicit none

!arguments
integer, intent(in)             :: j            !index of the current gridcell

!parameters
real(dp), parameter, dimension(band)  :: scatcf_s =  [ 0.8 , 0.4 ]       !snow scattering coeffic (vis, nir)
real(dp), parameter, dimension(band)  :: upsc0_s =  [ 0.5 , 0.5 ]        !snow upscattering coefficient (direct)
real(dp), parameter, dimension(band)  :: upsc_s = [ 0.5 , 0.5 ]           !snow upscattering coefficient (diffuse)
real(dp), parameter, dimension(band)  :: asoil_sat = [ 0.08 , 0.17 ]     !soil albedo saturated
real(dp), parameter, dimension(band)  :: asoil_dry = [ 0.17 , 0.34 ]     !soil albedo dry
real(dp), parameter, dimension(band)  :: bigC = [ 0.2 , 0.5 ]                   !empirical constant
real(dp), parameter, dimension(band)  :: newsnoalb = [ 0.95 , 0.65 ]     !albedo of new snow for sol zen < 60 deg
real(dp), parameter :: bpar = 2.0                                           !controls the solar zenith angle dependance

!pointers
logical, pointer, dimension(:)  :: virtual_leaf !true if the leaf being calculated is a virtual one to determine phenology.
real(dp), pointer               :: tausno       !non-dimensional snow age
real(dp), pointer               :: fsnow        !fraction of the gridcell covered by snow (fraction)
real(dp), pointer, dimension(:) :: f_wet        !fraction of the canopy that is wet
real(dp), pointer               :: mu           !cosine of the zenith angle (rad) of the incident beam
real(dp), pointer, dimension(:,:) :: Tveg_ave   !temperature of the vegetation (average time step value for the previous day)(K)
real(dp), pointer, dimension(:) :: stem_sno     !stem area index (m2 stem m-2 ground area) accounting for snow burial
real(dp), pointer, dimension(:) :: lai_sno      !individual leaf area index (accounting for snow burial)
real(dp), pointer, dimension(:) :: lai_vir      !virtual individual leaf area index (for phenology)
real(dp), pointer, dimension(:) :: lfang        !leaf angle for pft
real(dp), pointer, dimension(:) :: fi1          !relates leaf angles to projected area
real(dp), pointer, dimension(:) :: fi2          !
real(dp), pointer, dimension(:) :: alpha_s1     !stem reflectances (visible)
real(dp), pointer, dimension(:) :: alpha_s2     !stem reflectances (nir)
real(dp), pointer, dimension(:) :: alpha_l1     !leaf reflectances (visible)
real(dp), pointer, dimension(:) :: alpha_l2     !leaf reflectances (nir)
real(dp), pointer, dimension(:) :: tau_l1       !leaf transmittance (visible)
real(dp), pointer, dimension(:) :: tau_l2       !leaf transmittance (nir)
real(dp), pointer, dimension(:) :: tau_s1       !stem transmittance (visible)
real(dp), pointer, dimension(:) :: tau_s2       !stem transmittance (nir)
real(dp), pointer, dimension(:) :: Tliq         !timestep soil liquid water content (fraction)
real(dp), pointer, dimension(:,:) :: ftid       !downward diffuse fluxes per unit incident direct beam
real(dp), pointer, dimension(:,:) :: ftii       !downward diffuse fluxes per unit incident diffuse beam
real(dp), pointer, dimension(:,:) :: scatcf_tot !fraction of intercepted radiation that is scattered (0 to 1)
real(dp), pointer, dimension(:,:) :: scatcf_tot_vir !fraction of intercepted radiation that is scattered (0 to 1) (using virtual leaves)
real(dp), pointer, dimension(:,:) :: ftdd       !down direct flux below veg per unit dir flx
real(dp), pointer, dimension(:,:) :: fabd       !flux absorbed by vegetation per unit direct flux
real(dp), pointer, dimension(:,:) :: fabi       !flux absorbed by veg per unit diffuse flux
real(dp), pointer, dimension(:,:) :: fabd_vir   !flux absorbed by vegetation per unit direct flux (virtual leaves)
real(dp), pointer, dimension(:,:) :: fabi_vir   !flux absorbed by veg per unit diffuse flux (virtual leaves)
real(dp), pointer, dimension(:) :: bigK         !the optical depth of direct beam per unit leaf and stem area
real(dp), pointer               :: Wsno         !snow water equivalent of the snowpack (mm)
real(dp), pointer, dimension(:) :: mu_bar       !average inverse optical depth per unit leaf and stem area
logical, pointer, dimension(:)  :: tree
logical, pointer, dimension(:)  :: present

!local variables
integer :: t,k,v           !counter integer
integer :: pft              !counter
integer :: lim_loop             
real(dp) :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,temp0,temp1,temp2  !temp var.
real(dp) :: s2          !down direct flux below veg per unit dir flx
real(dp) :: Fage           !transformed snow age for solar zenith angle of less than 60 deg.
real(dp) :: del           !measure of soil wetness (used in soil albedo)
real(dp) :: Fmu                !function giving the increase in snow albedo due to the solar zenith angle exceeding 60 deg
real(dp) :: b,c,dterm,f,h,sig,u1,u2,u3,s1,p1,p2,p3,p4,d1,d2,h1,h2,h3,h4,h5,h6,h7,h8,h9,h10  !interm variables
real(dp), dimension(band,npft) :: asurfdir_u    !upward diffuse fluxes per unit incident direct beam  (surface albedos)
real(dp), dimension(band,npft) :: asurfdif_u    !upward diffuse fluxes per unit incident diffuse beam  (surface albedos)
real(dp), dimension(band) :: alsoil_dif    !soil albedo (diffuse)
real(dp), dimension(band) :: alsoil_dir    !soil albedo (direct)
real(dp), dimension(band) :: asnow_dif     !snow albedo (diffuse)
real(dp), dimension(band) :: asnow_dir     !snow albedo (direct)
real(dp), dimension(band) :: alphar        !weighted combination of the leaf and stem reflectances
real(dp), dimension(band) :: tau           !weighted combination of the leaf and stem transmittances
real(dp), dimension(band) :: scatupsc_v    !upscatter for diffuse vegetation
real(dp), dimension(band) :: scatupsc0_v   !upscatter for direct beam radiation
real(dp), dimension(band) :: ssalb         !single scattering albedo
real(dp), dimension(band) :: scatupsc      !
real(dp), dimension(band) :: scatupsc0     !
real(dp), dimension(band,npft) :: scatcf_v   !scattering coefficient for vegetation
real(dp), dimension(npft) :: Gmu          !relative projected area of leaves and stems in the direction cos-1 (mu)
real(dp) :: frac_leaf  !fraction of total leaf and stem that is leaf
real(dp) :: frac_stem  !fraction of total leaf and stem that is stem
real(dp) :: lai_temp
real(dp) :: albs                    ! temporary vis snow albedo
real(dp) :: albl                    ! temporary nir snow albedo
real(dp) :: czf                     ! solar zenith correction for new snow albedo [-]
real(dp) :: coszeni                 ! temporary cos of zenith angle variable
    
!point pointers
Tliq            => sv(j)%Tliq
alpha_l1        => prm%alpha_l1
alpha_l2        => prm%alpha_l2
alpha_s1        => prm%alpha_s1
alpha_s2        => prm%alpha_s2
ftii            => veg%ftii
ftid            => veg%ftid
f_wet           => veg%f_wet
lai_sno         => sv(j)%lai_sno
lai_vir         => sv(j)%lai_vir
lfang           => prm%lfang
!mu              => surf%avg_cZ         !for the daily timestep
mu              => surf%zen            !for the fast timestep
fsnow           => surf%fsnow
stem_sno        => sv(j)%stem_sno
tau_l1          => prm%tau_l1
tau_l2          => prm%tau_l2
tau_s1          => prm%tau_s1
tau_s2          => prm%tau_s2
fi1             => prm%fi1
fi2             => prm%fi2
tausno          => sv(j)%tausno
Tveg_ave        => sv(j)%Tveg_ave       !NOTE: this uses the value from the day/night before present timestep
scatcf_tot      => veg%scatcf_tot
scatcf_tot_vir  => veg%scatcf_tot_vir
ftdd            => veg%ftdd
fabi            => veg%fabi
fabd            => veg%fabd
fabi_vir        => veg%fabi_vir
fabd_vir        => veg%fabd_vir
bigK            => veg%bigK
Wsno            => sv(j)%Wsno
virtual_leaf    => sv(j)%virtual_leaf
mu_bar          => prm%mu_bar
tree            => prm%tree
present         => sv(j)%present

!--------

! Initial presets
ftdd(1:band,1:npft) = 1._dp
ftid(1:band,1:npft) = 0._dp  
ftii(1:band,1:npft) = 1._dp  
fabd(1:band,1:npft) = 0._dp
fabi(1:band,1:npft) = 0._dp
scatcf_tot(1:band,1:npft) = 0._dp
scatcf_v(1:band,1:npft) = 0._dp
Fmu = 0._dp
asnow_dif(1:band) = 0._dp
asnow_dir(1:band) = 0._dp
asurf_g(1:2,1:band) = 1._dp

if (mu > 0._dp) then  !sun must be above the horizon.

! Snow albedo calculations 

   if (Wsno > 0._dp) then
   
      ! calc. fractional area of the gridcell covered by snow (Niu and Yang 2007) 
      ! ->moved to newsnowflux  

      Fage = 1._dp - 1._dp / (1._dp + tausno) !transformed snow age (solar zenith angle less than 60 deg) eq. 3.55
      albs = newsnoalb(1) * (1._dp - bigC(1) * Fage)
      albl = newsnoalb(2) * (1._dp - bigC(2) * Fage)

      !for diffuse
      asnow_dif(1) = albs
      asnow_dif(2) = albl

      !for direct
      !the function Fmu is a factor between 0 and 1 giving the increase in snow albedo dur to the 
      !solar zenith angle exceeeding 60 deg 
      Fmu = max(0._dp, (1._dp + 1._dp/ bpar) / (1._dp + max(0.001_dp,mu) * 2._dp * bpar) - 1._dp / bpar )
      czf  = 0.4_dp * Fmu * (1._dp - albs)
      asnow_dir(1) = albs + czf
      czf  = 0.4_dp * Fmu * (1._dp - albl)
      asnow_dir(2) = albl + czf

   end if

!Ground albedo calculations
 
  do k = 1,band  !1 is vis, 2 is nir      

    !soil albedo (adjusted to soil wetness)
    del = max(0._dp,(0.11_dp - 0.4_dp * Tliq(1)))

    alsoil_dif(k) = min((asoil_sat(k) + del),asoil_dry(k))   !eq. 3.51  (both diffuse and direct)
    alsoil_dir(k) = alsoil_dif(k)                            !eq. 3.51  (both diffuse and direct)

    !overall direct (1) and diffuse (2) ground albedos weighted by snow and soil
    asurf_g(1,k) = alsoil_dir(k) * (1._dp - fsnow) + asnow_dir(k) * fsnow  !eq. 3.47
    asurf_g(2,k) = alsoil_dif(k) * (1._dp - fsnow) + asnow_dif(k) * fsnow  !eq. 3.48

  end do

  coszeni = max(0.001,mu)  !values < 0 already handled by initial if, so now acting on values between 0 and 0.001 

! Canopy albedos and fluxes
  do pft = 1,npft
  
    if (present(pft)) then !if the LAI and SAI are above zero
      
       !find the relative projected area of leaves and stems in the direction cos-1(mu)

       !precalc'd and placed in pftparametersmod
       !fi1(pft) = 0.5_dp - 0.633_dp * lfang(pft) - 0.33_dp * lfang(pft) * lfang(pft)
       !fi2(pft) = 0.877_dp * (1._dp - 2._dp * fi1(pft))
       !mu_bar = (1._dp - fi1(pft) / fi2(pft) * log((fi1(pft) + fi2(pft)) / fi1(pft)) ) / fi2(pft)  !eq. 3.4
       !NOTE: this is only valid for fi1(pft) and fi2(pft) both > 0. For the case where they are 0, check Appendix 1
       !of Dai et al. 2004 J. of Climate

      Gmu(pft) = fi1(pft) + fi2(pft) * coszeni   !eq. 3.3

      !find the optical depth of direct beam per unit leaf and stem area
      bigK(pft) = Gmu(pft) / coszeni  !CLM 4.0 Eqn 4.8
  
      tmp2 = mu_bar(pft) * bigK(pft)

        if (virtual_leaf(pft) .and. tree(pft)) then  !do virtual leaves first then redo calc for real leaf state.
                                                     !grass only enter once for the virtual leaves.
          lim_loop = 2 !run it once for the real leaf state and once for virtual leaves
        else
          lim_loop = 1 ! only run once for real leaf state.  
        end if

      do v = 1,lim_loop  !1 is the virtual leaf, 2 is the real leaf state.

         if (v == 1 .and. virtual_leaf(pft)) then
          lai_temp = lai_vir(pft)  !perform calc for real leaves
         else if (v == 1 .and. .not. virtual_leaf(pft)) then
          lai_temp = lai_sno(pft)  !perform calc for real leaves
         else !v == 2
          lai_temp = lai_sno(pft) 
         end if  

          if (lai_temp + stem_sno(pft) == 0.) cycle !if no lai or sai, then no influence on albedo, cycle

          frac_leaf = lai_temp / (lai_temp + stem_sno(pft))
          frac_stem = stem_sno(pft) / (lai_temp + stem_sno(pft))
         
          !vis weighted leaf and stem reflectances
          alphar(1) = alpha_l1(pft) * frac_leaf + alpha_s1(pft) * frac_stem  !eq. 3.11
          tau(1)   = tau_l1(pft) * frac_leaf + tau_s1(pft) * frac_stem       !eq. 3.12
          
          !nir weighted leaf and stem reflectances
          alphar(2)= alpha_l2(pft) * frac_leaf + alpha_s2(pft) * frac_stem   !eq. 3.11
          tau(2)  = tau_l2(pft) * frac_leaf + tau_s2(pft) * frac_stem        !eq. 3.12

         do k = 1,band  !1 is vis, 2 is nir

            !scattering coefficient for vegetation
            scatcf_v(k,pft) = alphar(k) + tau(k)

            !upscatter for direct radiation
              !find the single scattering albedo (eq. 3.15)

               temp0 = Gmu(pft) + fi2(pft) * coszeni
               temp1 = fi1(pft)*coszeni
               temp2 = 1._dp - temp1 / temp0 * log((temp1 + temp0) / temp1)

            ssalb(k) = scatcf_v(k,pft) * 0.5_dp * Gmu(pft) / temp0 * temp2

            !nir can become more than 1. when there is not direct beam, so set to 1 in that case
            scatupsc0_v(k) = (1._dp + tmp2) / (scatcf_v(k,pft) * tmp2) * ssalb(k)  !eq. 3.14 taken from CLM code

            !upscatter for diffuse radiation
            scatupsc_v(k) = 0.5_dp * (scatcf_v(k,pft) + (alphar(k) - tau(k)) * ((1._dp + lfang(pft)) * 0.5_dp)**2) / scatcf_v(k,pft) !eq. 3,13 from code

            !find the combined terms

             if (Tveg_ave(1,pft) > Tfreeze) then  !no snow

                  scatcf_tot(k,pft)= scatcf_v(k,pft)                        !eq. 3.8
                  scatupsc(k)      = scatupsc_v(k)                          !eq. 3.9
                  scatupsc0(k)     = scatupsc0_v(k)                         !eq.3.10

              else !snow in the canopy

                  scatcf_tot(k,pft)= scatcf_v(k,pft) * (1._dp - f_wet(pft))                    + scatcf_s(k) * f_wet(pft)                      !eq. 3.5
                  scatupsc(k)      = (scatupsc_v(k)   * scatcf_v(k,pft) * (1._dp - f_wet(pft)) + scatcf_s(k) * upsc_s(k) * f_wet(pft) )  / scatcf_tot(k,pft)  !eq. 3.6
                  scatupsc0(k)     = (scatupsc0_v(k)  * scatcf_v(k,pft) * (1._dp - f_wet(pft)) + scatcf_s(k) * upsc0_s(k) * f_wet(pft))  / scatcf_tot(k,pft)  !eq. 3.7

             end if

            !find upward and downward diffuse fluxes per unit incident direct beam and diffuse flux (surface albedos)
            !intial equations (CLM tech note 3.20 - 3.46)

            ! Absorbed, reflected, transmitted fluxes per unit incoming radiation
                 b = 1._dp - scatcf_tot(k,pft) + scatupsc(k) * scatcf_tot(k,pft)
                 c = scatupsc(k) * scatcf_tot(k,pft)
                 dterm = scatupsc0(k) * tmp2 * scatupsc0(k)
                 f = scatcf_tot(k,pft) * tmp2  * (1._dp - scatupsc0(k))
                 tmp1 = b * b - c * c

                 h = sqrt(tmp1) / mu_bar(pft)
                 sig = tmp2 * tmp2 - tmp1
                 s1 = exp(-(min(h * (lai_temp + stem_sno(pft)),40._dp)))  !updated to CLM 4.0 Eqn 3.30
                 s2 = exp(-(min(bigK(pft) * (lai_temp + stem_sno(pft)), 40._dp)))

                 p1 = b + mu_bar(pft) * h
                 p2 = b - mu_bar(pft) * h
                 p3 = b + tmp2
                 p4 = b - tmp2

                 ! Determine fluxes for vegetated pft for unit incoming direct
                 ! Loop over incoming direct and incoming diffuse
                 ! 1=incoming direct; 2=incoming diffuse

           do t = 1,2  !1 is direct, 2 is diffuse
                 u1 = b - c / asurf_g(t,k)
                 u2 = b - c * asurf_g(t,k)
                 u3 = f + c * asurf_g(t,k)

                 tmp3 = u1 - mu_bar(pft) * h
                 tmp4 = u1 + mu_bar(pft) * h

                 d1 = p1 * tmp3 / s1 - p2 * tmp4 * s1

                 tmp5 = u2 - mu_bar(pft) * h
                 tmp6 = u2 + mu_bar(pft) * h

                 d2 = tmp6 / s1 - tmp5 * s1
                 h1 = -dterm * p4 - c * f

                 tmp7 = dterm - h1 * p3 / sig
                 tmp8 = (dterm - c - h1 / sig * (u1 + tmp2)) * s2

                 h2 = (tmp7 * tmp3 / s1 - p2 * tmp8) / d1
                 h3 = -(tmp7 * tmp4 * s1 - p1 * tmp8) / d1
                 h4 = -f * p3 - c * dterm

                 tmp9 = h4 / sig
                 tmp10 = (u3 - tmp9 * (u2 - tmp2)) * s2

                 h5 = -(tmp9 * tmp6 / s1 + tmp10) / d2
                 h6 = (tmp9 * tmp5 * s1 + tmp10) / d2

                 h7 = (c * tmp3) / (d1 * s1)
                 h8 = (-c * tmp4 * s1) / d1
                 h9 = tmp6 / (d2 * s1)
                 h10 = (-s1 * tmp5) / d2

             if (t == 1) then !direct

                !Downward diffuse fluxes per direct beam below vegetation
                 ftid(k,pft) = tmp9 * s2 + h5 * s1 + h6 / s1 !eq. 3.18

                 !Downward direct flux below veg per unit dir flx (used in radflux)
                 ftdd(k,pft) = s2 

                 !Flux reflected by vegetation (direct albedo)
                 asurfdir_u(k,pft) = h1 / sig + h2 + h3 !eq. 3.16

                 !Flux absorbed by vegetation
                 fabd(k,pft) = 1._dp - asurfdir_u(k,pft) &
                      - (1._dp - asurf_g(1,k)) * ftdd(k,pft) - (1._dp - asurf_g(2,k))*ftid(k,pft)

             else  !diffuse

                 !Downward  diffuse fluxes per diffuse beam below vegetation
                 ftii(k,pft) = h9 * s1 + h10 / s1     !eq. 3.19

                 !Flux reflected by vegetation
                 asurfdif_u(k,pft) = h7 + h8                   !eq. 3.17

                 !Flux absorbed by vegetation
                 fabi(k,pft) = 1._dp - asurfdif_u(k,pft) - (1._dp - asurf_g(2,k)) * ftii(k,pft)

             end if !diffuse,direct if
           end do !diffuse/ direct loop
         end do !wavelength loop
    
         if (virtual_leaf(pft)) then !copy to virtual variables.
               scatcf_tot_vir(1:band,pft) = scatcf_tot(1:band,pft)
               fabd_vir(1:band,pft) = fabd(1:band,pft) 
               fabi_vir(1:band,pft) = fabi(1:band,pft)
         end if
         
         if (virtual_leaf(pft) .and. .not. tree(pft)) then  !when grass have no lai_sno, treat as no veg (no stems)
            ftdd(1:band,pft) = 1._dp
            ftid(1:band,pft) = 0._dp  
            ftii(1:band,pft) = 1._dp  
            fabd(1:band,pft) = 0._dp
            fabi(1:band,pft) = 0._dp
            scatcf_tot(1:band,pft) = 0._dp
            scatcf_v(1:band,pft) = 0._dp          
         end if

    end do !virtual, real    
  end if !min. leaf and stem areas.
 end do !pft loop

else 

         asurfdir_u(1:band,1:npft) = asurf_g(1,k)
         asurfdif_u(1:band,1:npft) = asurf_g(2,k)
                   
end if !mu > 0 

        !increment monthly albedo for output
      !  write(*,*)asurfdir_u(1,5),asurfdif_u(1,5),asurf_g(1,1),asurf_g(2,1)
!        albedo_mo_tot = albedo_mo_tot +  ) * dt/daysec


end subroutine albedo
