!SUBROUTINE decay
!Calculation of radioactive 14C decay
!recoded from f to f.90 by Joe Melton 2007 from original LPJ code

subroutine decay(j)

use arveparams,only : dp,npft,ncvar
use statevars, only : sv

implicit none

!arguments
integer, intent(in) :: j        !gridcell index

!PARAMETERS
real(dp), parameter :: lambda = 8267._dp

!Pointers
real(dp), pointer, dimension(:,:) :: lm_ind                !individual leaf mass (gC)
real(dp), pointer, dimension(:,:) :: sm_ind                !individual sapwood mass (gC)
real(dp), pointer, dimension(:,:) :: hm_ind                !individual heartwood mass (gC)
real(dp), pointer, dimension(:,:) :: rm_ind                 !individual fine root mass (gC)
real(dp), pointer, dimension(:,:) :: litter_ag_fast        !gridcell above-ground litter (gC/m2)
real(dp), pointer, dimension(:,:) :: litter_ag_slow        !gridcell above-ground litter (gC/m2)
real(dp), pointer, dimension(:,:) :: litter_bg                !gridcell below-ground litter (gC/m2)
real(dp), pointer, dimension(:)   :: cpool_fast          !fast-decomposing soil C pool (gC/m2)
real(dp), pointer, dimension(:)   :: cpool_slow                !slow-decomposing soil C pool (gC/m2)
logical, pointer, dimension(:)    :: present                !true if pft is present in that gridcell

!LOCAL VARIABLES
integer  :: pft,c        !counters
real(dp) :: dfac

!point pointers
lm_ind                  => sv(j)%lm_ind
sm_ind                => sv(j)%sm_ind
hm_ind                => sv(j)%hm_ind
rm_ind                => sv(j)%rm_ind
litter_ag_fast          => sv(j)%litter_ag_fast
litter_ag_slow          => sv(j)%litter_ag_slow
litter_bg             => sv(j)%litter_bg
cpool_fast          => sv(j)%cpool_fast
cpool_slow            => sv(j)%cpool_slow
present           => sv(j)%present

!-------------------------------------------------------------------------------
!begin calculations

dfac = exp(-1._dp / lambda)
c = 3

cpool_slow(c) = (((cpool_slow(c) / 1000._dp + 1._dp) * dfac) - 1._dp) * 1000._dp
cpool_fast(c) = (((cpool_fast(c) / 1000._dp + 1._dp) * dfac) - 1._dp) * 1000._dp



do pft = 1,npft
   lm_ind(pft,c) = (((lm_ind(pft,c) / 1000._dp + 1._dp) * dfac) - 1._dp) * 1000._dp
   sm_ind(pft,c) = (((sm_ind(pft,c) / 1000._dp + 1._dp) * dfac) - 1._dp) * 1000._dp
   hm_ind(pft,c) = (((hm_ind(pft,c) / 1000._dp + 1._dp) * dfac) - 1._dp) * 1000._dp
   rm_ind(pft,c) = (((rm_ind(pft,c) / 1000._dp + 1._dp) * dfac) - 1._dp) * 1000._dp
   litter_ag_fast(pft,c) = (((litter_ag_fast(pft,c) / 1000._dp + 1._dp) * dfac) - 1._dp) * 1000._dp
   litter_ag_slow(pft,c) = (((litter_ag_slow(pft,c) / 1000._dp + 1._dp) * dfac) - 1._dp) * 1000._dp
   litter_bg(pft,c) = (((litter_bg(pft,c) / 1000._dp + 1._dp) * dfac) - 1._dp) * 1000._dp
end do

return

end subroutine decay
