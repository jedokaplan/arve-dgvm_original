module leafarea

! Calculate the instantanous leaf and stem area due to snow burial. 
! calculate the sunlight and shaded leaf fractions and areas
! Also calculate the variant SLA for each PFT

use arveparams,          only : dp,npft,d2r,minsno
use statevars,           only : sv
use pftparametersmod,    only : prm
use pft_state_variables, only : veg
use iovariables,         only : do_allocation_daily

implicit none

public :: leafburial
public :: sla_twoleaf

contains

!--

subroutine  leafburial(j)

! Accounts for burial by snow of the vegetation. Also sums the
! leaf and stem areas. Created June 6 08 by JM

implicit none

!arguments
integer, intent(in) :: j   !gridcell index

!pointers
logical, pointer, dimension(:) :: present
logical,  pointer, dimension(:) :: tree         !true if pft is a tree
real(dp), pointer :: zsno                       !total thickness of the snowpack (m)
real(dp), pointer :: fpc_sum                    !instantaneous total foliar protective cover
real(dp), pointer, dimension(:) :: fpc_grid    !instantaneous foliar protective cover per gridcell (m2 m-2)
real(dp), pointer, dimension(:) :: lai_ind      !individual leaf area index
real(dp), pointer, dimension(:) :: lai_sno      !LAI with snow burial
real(dp), pointer, dimension(:) :: lai_vir      !virtual LAI with snow burial
real(dp), pointer, dimension(:) :: lai_vir_no_burial      !virtual LAI with no snow burial
real(dp), pointer, dimension(:) :: stema        !stem area index
real(dp), pointer, dimension(:) :: stem_sno     !stema with snow burial
real(dp), pointer, dimension(:) :: nind         !gridcell individual density (indiv/m2)
real(dp), pointer, dimension(:) :: fpc_ind     !instantaneous foliar protective cover per pft individual (m2 m-2)
real(dp), pointer, dimension(:) :: crownarea    !crown area(m2)
integer, pointer, dimension(:) :: pstate       !nstantaneous phenological state for deciduous plants 
                                                        ! (0=no leaves, 1=full leaf out, 2=leafing out, 3=senescing)
real(dp), pointer, dimension(:) :: lai_max
real(dp), pointer, dimension(:) :: height
logical, pointer, dimension(:)  :: virtual_leaf !true if the leaf being calculated is a virtual one to determine phenology.
real(dp), pointer, dimension(:) :: phfrac               !instantaneous phenological state fraction 
                                                                ! (0=no leaves, 1=full leaf out)
                                                                
!local variables
integer :: pft
real(dp) :: stem_temp
real(dp) :: fsnoveg      !vertical fraction of vegetation covered by snow

!local parameters
real(dp), parameter :: exstem = 0.1     !minimum fraction of stems that are exposed

!point pointers
present         => sv(j)%present
fpc_sum         => veg%fpc_sum
fpc_grid        => veg%fpc_grid
lai_ind         => sv(j)%lai_ind
lai_sno         => sv(j)%lai_sno
lai_vir         => sv(j)%lai_vir
lai_vir_no_burial => sv(j)%lai_vir_no_burial
zsno            => sv(j)%zsno
stema           => sv(j)%stema
stem_sno        => sv(j)%stem_sno
nind            => sv(j)%nind
fpc_ind         => veg%fpc_ind
pstate          => sv(j)%pstate
crownarea       => sv(j)%crownarea
tree            => prm%tree
lai_max         => sv(j)%lai_max
height          => sv(j)%height
virtual_leaf    => sv(j)%virtual_leaf
phfrac          => sv(j)%phfrac

!---------------

fpc_sum = 0._dp

if (.not. do_allocation_daily)  lai_ind(:) = lai_max(:) * phfrac(:)  
           
! Find the leaf area index sum and foliar projective cover sum for all pfts accounting for snow burial

 do pft = 1,npft
        
  if (present(pft)) then
  
      if (tree(pft)) then   ! Stem area as related to maximal leaf area. This is from CLM 3.5 dynamic ecosystem module.
       stema(pft) = lai_max(pft) * 0.25_dp
      else
       stema(pft) = lai_max(pft) * 0.05_dp
      end if

      if (zsno > minsno) then !calculate vegetation burial by snow
              
              if (tree(pft)) then !trees have neglible burying so add as normal.

                lai_sno(pft)  = lai_ind(pft)
                stem_sno(pft) = stema(pft) * (max(exstem,1._dp - (lai_ind(pft)/lai_max(pft)))) 

              else !grass, will be buried to different extents dependent upon the phen state.

              if (do_allocation_daily) then
              
                ! in the daily C allocation scheme, the grass height is not constant so use that to
                ! determine burial.
              
                if (height(pft) /= minsno) then
                  fsnoveg = min(1._dp,(zsno - minsno) / (height(pft) - minsno))
                  fsnoveg = max(0._dp,fsnoveg)
                else
                  fsnoveg = 1._dp
                end if
              
              else 
              
               ! annual C allocation. grass height is constant at 1m. So for leaf burial, use the
               ! phfrac to scale the plant height
                if (height(pft) * phfrac(pft) /= minsno) then
                  fsnoveg = min(1._dp,(zsno - minsno) / (height(pft) * phfrac(pft) - minsno))
                  fsnoveg = max(0._dp,fsnoveg)
                else
                  fsnoveg = 1._dp
                end if
             
              end if !annual/daily C allocation
               
                stem_temp= stema(pft) * (1._dp - fsnoveg) * (max(exstem,1._dp - (lai_ind(pft)/lai_max(pft))))
                stem_sno(pft) = stem_temp

                if (.not. virtual_leaf(pft)) then  !burial of 'real' leaves
                
                   lai_sno(pft) = lai_ind(pft) * (1._dp - fsnoveg)

                else !burial of 'virtual' leaves
                  
                  ! the lai_vir needs to account for snow burial too, otherwise the model will
                  ! allow for leafout due to lai_vir NPP but once the leaves are put on, snow
                  ! burial makes lai_sno = 0 so they have to retreat. Accounting for burial 
                  ! avoids this.
                  
                   lai_vir(pft) = lai_vir_no_burial(pft) * (1._dp - fsnoveg)

                end if

              end if
 
      else ! no burial

            lai_vir(pft) = lai_vir_no_burial(pft)
            lai_sno(pft) = lai_ind(pft)
            stem_sno(pft) = stema(pft) * (max(exstem,1._dp - (lai_ind(pft)/lai_max(pft))))

      end if !snow/no snow      
      
      ! Calculate instantaneous FPC
      
      ! FPC (foliar projective cover individual) is the area of ground covered directly by foliage above it
      ! original LPJ formulation of FPAR. NOTE this is not used to determine the amount of PAR absorbed as in LPJ, 
      ! it is used only in FPC calcs and in daily_allocation to help determine how the plant allocates C
      fpc_ind(pft) = 1._dp - exp(-0.5_dp * lai_sno(pft)) 
      
      !instantaneous FPC_grid values (FPC as a fraction of the total grid cell)
      fpc_grid(pft) = crownarea(pft) * nind(pft) * fpc_ind(pft)
          
 else !not present
 
    fpc_grid(pft) = 0._dp
    fpc_ind(pft) = 0._dp
    
 end if !present 
end do !pft loop

   !find fpc_sum
   fpc_sum = sum(fpc_grid(1:npft))
   
   ! Can not allow fpc_sum to be greater than 1, FLAG JM 19.08.2010
   fpc_sum = min(1._dp, fpc_sum)
   

end subroutine leafburial

!----------------------------------------------------------------

subroutine sla_twoleaf(j,pft)

! Calculates the sunlight and shaded LAI. Created Feb 05 09 by JM from CLM base.

implicit none

!arguments
integer, intent(in) :: j   !gridcell index
integer, intent(in) :: pft !plant functional type

!pointers
real(dp), pointer :: lai_sun      !instantaneous sunlit leaf area (m2 m-2)
real(dp), pointer :: lai_sha      !instantaneous shaded leaf area (m2 m-2)
real(dp), pointer :: bigK         !the optical depth of direct beam per unit leaf and stem area
real(dp), pointer :: sla_sun      !specific leaf area sunlight (m2 gC-1)
real(dp), pointer :: sla_sha      !specific leaf area shaded (m2 gC-1)
real(dp), pointer :: m            !linear slope coefficient (dSLA/dLAI, projected area basis [m^2/gC])
real(dp), pointer :: sla_top      !SLA at the top of the canopy (m2 / gC)
real(dp), pointer :: stem_sno     !stema with snow burial
real(dp), pointer :: lai_sno      !LAI with snow burial
real(dp), pointer, dimension(:) :: fsun         !sunlit fraction of canopy
real(dp), pointer, dimension(:) :: fsha         !shaded fraction of canopy

!local variables
real(dp) :: c1
real(dp) :: temp

!point pointers
lai_sun  => veg%lai_sun(pft)
lai_sha  => veg%lai_sha(pft)
bigK     => veg%bigK(pft)
sla_sun  => veg%sla_sun(pft)
sla_sha  => veg%sla_sha(pft)
m        => prm(pft)%m
sla_top  => prm(pft)%sla_top
stem_sno => sv(j)%stem_sno(pft)
lai_sno  => sv(j)%lai_sno(pft)
fsun     => veg%fsun
fsha     => veg%fsha

!--------------
!begin calculations

  if (lai_sno >= 0.01_dp .and. bigK /= 0._dp) then

        ! Calculate the one-sided leaf area (NOTE: special CLM forulation)
        temp = min(40._dp,bigK * lai_sno)

        c1 = exp(-temp)

        !sunlight fraction of the canopy
        fsun(pft) = (1._dp - c1) / temp  !CLM 4.0 Eqn. 4.7

        !shaded fraction of the canopy
        fsha(pft) = 1._dp - fsun(pft)

        ! Find the instantaneous sunlit leaf area index (m2 m-2) per pft
        lai_sun = fsun(pft) * lai_sno

        ! Find the instantaneous shaded leaf area index (m2 m-2) per pft
        lai_sha = lai_sno - lai_sun

        ! SLA for the sunlight part of the canopy
        sla_sun = -(c1 * m * bigK * lai_sno + &
                        c1 * m + &
                        c1 * sla_top * bigK - &
                        m - &
                        sla_top * bigK) / &
                        (bigK * bigK * lai_sun)
        
        ! SLA for the shaded
        sla_sha = ((sla_top + &
                       m * lai_sno / 2._dp) * lai_sno - &
                       lai_sun * sla_sun) / &
                       lai_sha

  else

        lai_sun = lai_sno
        lai_sha = 0._dp
        fsun(pft) = 1._dp
        fsha(pft) = 0._dp
        sla_sun = sla_top
        sla_sha = 0._dp

  end if

end subroutine sla_twoleaf

end module leafarea
