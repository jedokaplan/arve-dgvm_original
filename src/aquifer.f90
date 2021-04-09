subroutine aquifer(j,dft,dTliq,dzmm_aq)

!Unconfined aquifer is added. Based on Niu et. al 2007 but revamped in CLM 4.0
!Code adapted from CLM 4.0 by JM (22.04.2010)

use arveparams,only : dp,nl,fdecay,Wpond,Tfreeze,alpha,Pmax
use statevars, only : sv,doy
use soilstate_vars, only : surf

implicit none

!arguments
integer, intent(in) :: j                           !gridcell index                                  
real(dp), intent(in) :: dft                        !fast timestep(s)                                
real(dp), intent(in), dimension(1:nl+1) :: dTliq   !water flowing into unconfined aquifer (mm)   
real(dp), intent(in) :: dzmm_aq                    !aquifer layer thickness (mm)   

!parameters
real(dp), parameter :: Sy = 0.2                    !the fraction of water volume that can be drained by gravity in an unconfined aquifer
real(dp), parameter :: aquif_max = 5000._dp        !max water in aquifer below soil profile (mm)
real(dp), parameter :: qdraimax = 5.5e-3           !This new version is from CLM 4.0 tech note p.155
                                                   !old=4.5e-4 !maximum subsurface drainage when the grid avg water table is 0 (kg m-1 s-1)
real(dp), parameter :: Wlmin = 0.01_dp             !min water content per layer (mm)

real(dp), parameter :: ealpha        = exp(-alpha)                          !for computational efficiency
real(dp), parameter :: oneoverealpha = 1._dp / (1._dp - exp(-alpha))        !for computational efficiency

!variables
real(dp) :: qrecharge                      !recharge to the aquifer (mm s-1) positive when water enters aquifer
real(dp) :: qdrai                          !liquid water draining from the soil (mm s-1)
real(dp) :: fimp                           !fraction of impermeable area
integer :: l,i                             !counter
real(dp) :: Wexc                           !total soil column water surplus (mm s-1)
real(dp) :: wh_zt                          !water head at the water table depth
real(dp), dimension(nl) :: eff_porosity    !measure of how much pore space is available, adjusted for ice (fraction)
real(dp) :: dzmmsum                        !soil thickness between water table and bottom of soil column (mm)
real(dp) :: icefracsum                     !fraction of the soil between water table and bottom of soil column that contains ice 
real(dp) :: Sr,s1                          !
real(dp) :: ka                             !
real(dp) :: smp1                           !
real(dp) :: wh                             !
real(dp) :: ws                             !water used to fill soil air pores regardless of water content
real(dp) :: Waquif_tsub                    !
real(dp) :: xsi                            !
real(dp) :: xs1,xs                         !
real(dp) :: available_Wliq                 !
real(dp) :: zw_old                         !previous timestep's water level

!pointers
integer, pointer :: jint                   !index of the soil layer directly above the water table 
integer, pointer  :: gnl                   !index value of lowest 'soil' layer                
real(dp), pointer  :: zw                   !mean water table depth (m)                             
real(dp), pointer  :: Waquif_a             !water in unconfined aquifer below soil (mm)            
real(dp), pointer  :: Waquif_t             !water in aquifer within soil column (mm).              
real(dp), pointer, dimension(:) :: Ksat    !soil water saturated conductivity at layer midpoint (mm s-1)
real(dp), pointer, dimension(:) :: Psi     !soil water potential at layer midpoint (mm)
real(dp), pointer, dimension(:) :: zpos    !z coordinate of the middle of the soil layer relative to soil surface (m),positive downwards
real(dp), pointer, dimension(:) :: Psat    !soil water matric potential at saturation (mm)
real(dp), pointer, dimension(:) :: dzmm    !thickness of the soil layers (mm)
real(dp), pointer, dimension(:) :: zipos   !z coordinate of the interface between soil layers relative to soil surface (m), positive downwards
real(dp), pointer, dimension(:) :: Tice    !soil ice content (fraction)
real(dp), pointer, dimension(:) :: Tsat    !soil water volumetric water content at saturation (fraction)
real(dp), pointer, dimension(:) :: Tpor    !soil volumetric porosity (liquid minus ice content (or stones?)) (fraction)
real(dp), pointer, dimension(:) :: Wliq    !soil liquid water content at layer midpoint (mm)
real(dp), pointer, dimension(:) :: Wice    !soil ice content at layer midpoint (mm)
real(dp), pointer, dimension(:) :: Bexp    !
real(dp), pointer, dimension(:) :: Psi_eq  !
real(dp), pointer, dimension(:) :: Ku      !soil water instantaneous (unsaturated) conductivity across layer boundary (mm s-1)
real(dp), pointer               :: qover   !liquid water surface runoff (kg m-2 sec-1)

!point pointers
gnl             => sv(j)%gnl
Ku              => sv(j)%Ku
jint            => sv(j)%jint
zw              => sv(j)%zw
Waquif_a        => sv(j)%Waquif_a
Waquif_t        => sv(j)%Waquif_t
Ksat            => sv(j)%Ksat
Psi             => sv(j)%Psi
zpos            => sv(j)%zpos
zipos           => sv(j)%zipos
Psat            => sv(j)%Psat
dzmm            => sv(j)%dzmm
Tice            => sv(j)%Tice
Tsat            => sv(j)%Tsat
Tpor            => sv(j)%Tpor
Wliq            => sv(j)%Wliq
Wice            => sv(j)%Wice
Bexp            => sv(j)%Bexp
Psi_eq          => sv(j)%Psi_eq
qover           => surf%qover

!--------------------

!Initial presets.
  Wexc = 0._dp
  zw_old = zw
  
!write(*,'(i5,f12.4)')doy,surf%qliq*dft
!write(*,'(17f12.4)')Wliq(1:nl)
!write(*,'(17f12.4)')dTliq(1:nl)*dzmm(1:nl)

 ! Renew the mass of liquid water (add in change from this timestep)
   Wliq(1:gnl) = Wliq(1:gnl) + dTliq(1:gnl) * dzmm(1:gnl)

       ! calculate qrecharge 
       if(jint < gnl) then !if the water table is within the soil column

          !water head at the water table depth
          wh_zt = 0._dp   !defined as zero.

          Sr = max(Wliq(jint)/dzmm(jint)/Tsat(jint), 0.01_dp)
          Sr = min(1._dp, Sr)

          !use average moisture between water table and layer Waquif_t (which has to be 1)
          s1 = 0.5_dp * (1._dp + Sr)
          s1 = min(1._dp, s1)

          !this is the expression for unsaturated hydraulic conductivity
          ka = Ksat(jint) * s1**(2._dp * Bexp(jint) + 3._dp)

          ! Recharge rate qrecharge to groundwater (positive to aquifer)
          smp1 = Psat(jint) * Sr**(-Bexp(jint))  
          smp1 = max(Pmax, Psi(jint))
          wh      = smp1 - Psi_eq(jint)
          qrecharge = -ka * (wh_zt - wh)  /((zw - zpos(jint)) * 1000._dp)

          ! To limit qrecharge  (for the first several timesteps)
          qrecharge = max(-10.0_dp * (nint(dft/1200)) / dft,qrecharge)
          qrecharge = min( 10.0_dp * (nint(dft/1200)) / dft,qrecharge)

       else

          !if water table is below soil column, compute qrecharge from dTliq(gnl+1)
          qrecharge = dTliq(gnl+1) * dzmm_aq / dft
    
       end if

    eff_porosity(1:gnl) = max(0.01_dp,Tsat(1:gnl) - Tice(1:gnl))

    ! Topographic runoff
    
       dzmmsum = sum(dzmm(jint:gnl)) 

       icefracsum = sum((Wice(jint:gnl) / (Wice(jint:gnl) + Wliq(jint:gnl)) * dzmm(jint:gnl)))

       fimp = max(0._dp,(exp(-alpha * (1._dp - (icefracsum / dzmmsum)))- ealpha) * oneoverealpha)
      
       qdrai = (1._dp - fimp) * qdraimax * exp(-fdecay * zw)

    ! Water table calculation

       ! Water storage in aquifer + soil
       Waquif_t = Waquif_t + (qrecharge - qdrai) * dft

       if (jint == gnl) then             ! water table is below the soil column

          Waquif_a  = Waquif_a + (qrecharge - qdrai) * dft
          Waquif_t  = Waquif_a
          zw = (zipos(gnl) + 25._dp) - Waquif_a / 1000._dp / Sy
          Wliq(gnl) = Wliq(gnl) + max(0._dp,(Waquif_a - aquif_max))
          Waquif_a  = min(Waquif_a, aquif_max)

       else   ! water table within soil layers

          if (jint == gnl-1) then       ! water table within bottom soil layer

             zw = zipos(gnl)- (Waquif_t - Sy * 1000._dp * 25._dp) / eff_porosity(gnl) / 1000._dp

          else  ! water table within soil layers 1-gnl-1

             ws = 0._dp   ! water used to fill soil air pores regardless of water content

             do l = jint+2,gnl
               ws = ws + eff_porosity(l) * dzmm(l)
             end do

             zw = zipos(jint+1) - (Waquif_t - Sy * 1000_dp * 25._dp - ws) / eff_porosity(jint+1) / 1000._dp

          end if

          Waquif_tsub = 0._dp
          do l = jint+1, gnl
             Waquif_tsub = Waquif_tsub + Ku(l) * dzmm(l)
          end do

    ! Remove subsurface runoff
	  
	  ! NOTE! CLM subsurface runoff does not make sense. As the lateral subsruface runoff is removed
	  ! it simply disappears from the model as there is no lateral transfer of water between gridcells. In ARVE
	  ! we then assume that any water that should laterally transfer is compensated by lateral transfer
	  ! from the other cells adjoining (i.e. it is a zero-sum transfer). JM 25.10.10.
          ! Further note (20.06.2011 JM): This actually can not be turned off, it prevents small fluctuations in 
          ! the water table depth. Don't see an easy way around this. It might however lead to too dry soils...
	  
          do l = jint+1, gnl
            if (Waquif_tsub > 0._dp) then !FLAG this was causing a floating inexact without
               ! this if statement. The if statement is somehow not required in CLM code. It causes problems in areas with almost pure
               ! sand soils. Essentially the subsurface runoff keeps pulling water from the soil even if it is bone dry.
               ! so I made it so that it can not pull water if the conductivity is zero, which makes sense. JM 22.06.2011
                Wliq(l) = Wliq(l) - qdrai * dft * Ku(l) * dzmm(l) / Waquif_tsub
             end if
          end do
	  
       end if

        zw = max(-2.0_dp,zw)  !NOTE: changed to allow standing water in wetlands! JM 27.10.10
        zw = min(80._dp,zw)

    ! Handle any excess water

         !  excessive water above saturation added to the above unsaturated layer like a bucket
         !  if column fully saturated, excess water goes to runoff

         do l = gnl,2,-1
               xsi = max(Wliq(l) - eff_porosity(l) * dzmm(l),0._dp)
               Wliq(l) = min(eff_porosity(l) * dzmm(l), Wliq(l))
               Wliq(l-1) = Wliq(l-1) + xsi
         end do

         ! top layer
            xs1 = max(max(Wliq(1),0._dp) - max(0._dp,(Wpond + Tsat(1) * dzmm(1) - Wice(1))),0._dp)
            Wliq(1) = min(max(0._dp, Wpond + Tsat(1) * dzmm(1) - Wice(1)), Wliq(1))
            Wexc     = xs1 / dft

            ! add the excess water to the over land run-off
            ! NOTE: this is a change from CLM which adds it to qdrai. However, that does not make
            ! sense that the excess water at the surface should be added to subsurface drainage since it is not
            ! below the surface. We add it to over land runoff instead. JM 27.10.10
           ! qover = qover + Wexc
            qdrai = qdrai + Wexc
  
    ! Handle any excess dryness

           ! Limit Wliq to be greater than or equal to Wlmin.
           ! Get water needed to bring Wliq equal Wlmin from lower layer.
           ! If insufficient water in soil layers, get from aquifer water

           do l = 1,gnl-1
                 if (Wliq(l) < Wlmin) then
                    xs = Wlmin - Wliq(l)
                 else
                    xs = 0._dp
                 end if
                 Wliq(l) = Wliq(l) + xs
                 Wliq(l+1) = Wliq(l+1) - xs
           end do

       ! Get water for bottom layer from layers above if possible
           l = gnl
              if (Wliq(l) < Wlmin) then
                 xs = Wlmin - Wliq(l)
                 searchforwater: do i = gnl-1,1,-1

                    available_Wliq = max(Wliq(i) - Wlmin - xs,0._dp)

                    if (available_Wliq >= xs) then
                      Wliq(l) = Wliq(l) + xs
                      Wliq(i) = Wliq(i) - xs
                      xs = 0._dp
                      exit searchforwater
                    else
                      Wliq(l) = Wliq(l) + available_Wliq
                      Wliq(i) = Wliq(i) - available_Wliq
                      xs = xs - available_Wliq
                    end if
                 end do searchforwater
              else
                 xs = 0._dp
              end if

       ! Needed in case there is no water to be found
              Wliq(l) = Wliq(l) + xs
              Waquif_t = Waquif_t - xs

!write(*,'(17f12.4)')Wliq(1:nl)
!write(*,*)

! if (Wexc > 0.) then
! write(*,*)'high Wexc'
! write(0,'(f12.4,3i5,f12.4)')Wexc,sv(j)%snl,gnl,jint,zw
! write(0,'(17f12.4)')sv(j)%Tliq(1:nl)
! write(*,*)
! write(0,'(17f12.4)')sv(j)%Wliq(1:nl)
! write(*,*)
! write(0,'(17f12.4)')sv(j)%Tice(1:nl)
! write(*,*)
! write(0,'(17f12.4)')Tsat(1:nl)
! write(*,*)
! write(0,'(17f12.4)')eff_porosity(1:nl) * dzmm(1:nl)
! read(*,*)
! end if

       ! Instead of removing water from aquifer, take it back out of
       ! drainage. !FLAG Should this be part of an if statement then? JM 07.04.2011
             ! qdrai = qdrai - xs / dft

        !Check for large jumps in water table and limit
        if (ABS(zw_old - zw) > 0.1) then

          if (zw_old - zw > 0.d0) then
              zw = zw_old - 0.1
          else 
                zw = zw_old + 0.1
          end if

        end if
        
end subroutine aquifer
