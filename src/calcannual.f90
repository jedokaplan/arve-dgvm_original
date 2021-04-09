module calcannual

! Central module for all annual timestep calculations INCLUDING annual and daily C 
! allocation. 

use arveparams,                 only : dp,npft,ncvar
use statevars,                  only : sv, ov
use metvarsmod,                 only : atm
use increment,                  only : annual_increment
use iovariables,                only : do_allocation_daily
use pft_state_variables,        only : veg
use soilstate_vars,             only : surf       

implicit none

!Annual subroutines, not all used by each C allocation scheme.
public  :: annual_calcs
public  :: bioclimatic
public  :: reproduction
public  :: turnover
public  :: killplant
public  :: allocation
public  :: light
public  :: mortality
public  :: fire
public  :: establishment
public  :: roots

!------

contains

subroutine annual_calcs(j,yr)

!Central annual driver

implicit none

!arguments
integer, intent(in) :: j   !gridcell index
integer, intent(in) :: yr  !run year

!pointers
logical, pointer, dimension(:)          :: present              !true if pft is present in gridell
real(dp), pointer, dimension(:)         :: aleafdays            !count of days with foliage this year
real(dp), pointer, dimension(:)         :: awscal               !annual water stress factor
real(dp), pointer, dimension(:,:)       :: annual_npp           !annual net primary productivity (g C m-2 yr-1)
real(dp), pointer, dimension(:,:)       :: annual_resp          !annual autotrophic respiration (gC m-2 yr-1)
real(dp), pointer, dimension(:,:)       :: bm_inc               !annual biomass increment (gC/m2)
real(dp), pointer, dimension(:)         :: arh                  !annual heterotrophic respiration (gC m-2 yr-1)
real(dp), pointer, dimension(:)         :: Aliq                 !daily soil water content for layer 1
real(dp), pointer, dimension(:,:)       :: annual_gpp           !annual gross primary productivity (g C m-2 yr-1)
real(dp), pointer                       :: absrad_annual        !annual total radiation absorbed (J m-2)
real(dp), pointer                       :: zw                   !groundwater depth (m)
real(dp), pointer, dimension(:)         :: reprod               !reproduction cost (gC/m2)  
real(dp), pointer                       :: twm                  !maximum monthly temperature for the year being simulated
real(dp), pointer                       :: tminm                !minimum monthly temperature for the year being simulated
real(dp), pointer                       :: aprec                !annual precipitation (mm)
real(dp), pointer, dimension(:)         :: aaet                 !annual actual evapotranspiration (mm/yr)
        
!local variables
integer :: l            !counter

!point pointers
present         => sv(j)%present
aleafdays       => sv(j)%aleafdays
awscal          => veg%awscal
bm_inc          => veg%bm_inc
annual_npp      => veg%annual_npp
annual_resp     => veg%annual_resp
arh             => veg%arh
annual_gpp      => veg%annual_gpp
Aliq            => surf%Aliq
absrad_annual   => surf%absrad_annual
zw              => sv(j)%zw
reprod          => veg%reprod
twm             => atm%twm
tminm           => atm%tminm
aprec           => atm%aprec
aaet            => veg%aaet

!---------------

! Begin calculations

! Make annual averages of daily incremented variables that need averages or sums
  call annual_increment(j)

!----------------------------Print some output
goto 200
    !write out some vegetation output to screen
        if (yr > 1 .and. any(present)) then
          write(*,'(a,i5,2f10.2)')'end of year totals for YEAR:',yr,sv(j)%lon,sv(j)%lat
           write(*,4)'PFT','AGPP','autoresp','ANPP','Rh','LAImax','Height','nind','fpc_gmax','leafon days','ZW'
        do l = 1,npft
          if (sv(j)%present(l)) then
           write(*,5)l,annual_gpp(l,1),annual_resp(l,1),annual_npp(l,1),arh(1),sv(j)%lai_max(l),sv(j)%height(l),&
                        sv(j)%nind(l),sv(j)%fpc_gmax(l),aleafdays(l),zw
          end if
        end do
            write(*,2)'   ','crown',' rootdep',' lm_ind','rm_ind',' awscal',' fastlit',' smind',' hm-ind','fireC'
        do l = 1,npft
          if (sv(j)%present(l)) then
 write(*,3)l,sv(j)%crownarea(l),sv(j)%rootdepth(l),sv(j)%lm_ind(l,1),sv(j)%rm_ind(l,1),veg%awscal(l),sv(j)%litter_ag_fast(l,1),sv(j)%sm_ind(l,1)&
                     ,sv(j)%hm_ind(l,1),ov(j)%acflux_fire(1)
           end if
        end do
       end if
2 format(a5,10a10)
3 format(i5,5f10.3,f10.3,4f10.3)
4 format(a5,10a10)
5 format(i5,4f10.3,6f10.3)
200 continue
!--------------------------

! Implement PFT bioclimatic limits
call bioclimatic(yr,j)

! Allocation to reproduction (annual C allocation scheme)
if (.not. do_allocation_daily) call reproduction(j)
                
! Calculation of leaf, sapwood, and fine-root turnover (annual C allocation scheme)
if (.not. do_allocation_daily) call turnover(j)
                
! Removal of PFTs with negative C increment this year
call killplant(j)
                
! Allocation of annual carbon increment to leaf, stem and fine
! root compartments (annual C allocation scheme)
if (.not. do_allocation_daily) call allocation(j)
                
! Implement light competition between trees and grasses
call light(j)
                
! Implement tree background and stress mortality
! (including heat damage and due to lower limit of npp for boreal trees)
call mortality(j)
                
! Calculation of biomass destruction by fire disturbance
call fire(j)
                
! Establishment of new individuals (saplings) of woody PFTs, grass establishment,
! removal of PFTs not adapted to current climate,update of individual structure and FPC.
call establishment(j)
                
! Root distribution, rooting depth, and rooting fraction for depth 
! (annual C allocation scheme, but used by daily too)
if (.not. do_allocation_daily) call roots(j)

! Radioactive decay of 14C
! call decay()

!----------Print some output
goto 300
    !write out some vegetation output to screen
        if (any(present)) then
          write(*,'(a,i5,2f10.2)')'Calc annual end: YEAR:',yr,sv(j)%lon,sv(j)%lat
           write(*,4)'PFT','AGPP','tot_resp','ANPP','Rh','LAImax','Height','nind','fpc_gmax','leafon days','ZW'
        do l = 1,npft
          if (sv(j)%present(l)) then
           write(*,5)l,annual_gpp(l,1),annual_resp(l,1),annual_npp(l,1),arh(1),sv(j)%lai_max(l),sv(j)%height(l),&
                        sv(j)%nind(l),sv(j)%fpc_gmax(l),aleafdays(l),zw
          end if
        end do
            write(*,2)'   ','crown',' rootdep',' lm_ind','rm_ind',' awscal',' fastlit',' smind',' hm-ind','fireC',' lm_max'
        do l = 1,npft
          if (sv(j)%present(l)) then 
          write(*,3)l,sv(j)%crownarea(l),sv(j)%rootdepth(l),sv(j)%lm_ind(l,1),sv(j)%rm_ind(l,1),veg%awscal(l),&
          sv(j)%litter_ag_fast(l,1),sv(j)%sm_ind(l,1),sv(j)%hm_ind(l,1),ov(j)%acflux_fire(1),sv(j)%lm_max(l,1)
           end if
        end do
       end if
300 continue
                     
!--------------------------

! assign values to some annual output variables
  ov(j)%grid_npp = real(sum(annual_npp(:,1)*sv(j)%fpc_gmax,mask=sv(j)%present))
  ov(j)%grid_gpp = real(sum(annual_gpp(:,1)*sv(j)%fpc_gmax,mask=sv(j)%present))
  ov(j)%grid_autoresp = real(sum(annual_resp(:,1)*sv(j)%fpc_gmax,mask=sv(j)%present))
  ov(j)%plant_carbon = real(sum(sv(j)%nind(:) * (sv(j)%lm_ind(:,1) + sv(j)%sm_ind(:,1) + sv(j)%hm_ind(:,1) + sv(j)%rm_ind(:,1))))  
  ov(j)%annual_hetresp=real(arh(1))
  ov(j)%leafondays = real(aleafdays)
  ov(j)%abs_rad_annual = real(absrad_annual / 1.e9)
  ov(j)%grndwat_depth = real(zw)
  ov(j)%aaet_out = real(aaet)


! End of year resets
      veg%gdd_bioclim = 0._dp
      aaet  = 0._dp
      arh(1)   = 0._dp
      annual_npp(:,1) = 0._dp
      annual_gpp(:,1) = 0._dp
      annual_resp(:,1)= 0._dp
      awscal(:) = 0._dp
   !   twm = -100._dp
   !   tminm = 100._dp
      aprec = 0._dp
      aleafdays = 0._dp
      Aliq = 0._dp
      absrad_annual = 0._dp
      reprod = 0._dp

end subroutine annual_calcs


!-------------------------------------------------------------------------------------------------------------

subroutine bioclimatic(yr,j)

! Apply bioclimatic limits on PFT survival and establishment
! Has been changed to be able to operate without a priori knowledge of the upcoming
! years climate data. Originally based upon LPJ's approach.
! coded by Joe Melton 2008

use arveparams,                 only : dp,npft,climbuf
use pftparametersmod,           only : prm
use statevars,                  only : sv
use metvarsmod,                 only : atm
use pft_state_variables,        only : veg

implicit none

! arguments
integer, intent(in)  :: yr           !year of the simlulation
integer,  intent(in) :: j            !grid index

! Pointers
real(dp), pointer, dimension(:) :: gdd                !growing degree days at or above 5 deg C base   
real(dp), pointer, dimension(:) :: mtemp_min_buf      !array of coldest month temperature values for climbuf years of data          
real(dp), pointer, dimension(:)  :: mintcm            !PFT-specific minimum coldest-month temperature
real(dp), pointer, dimension(:)  :: maxtcm            !PFT-specific maximum coldest-month temperature
real(dp), pointer, dimension(:)  :: gddmin            !PFT-specific minimum GDD
real(dp), pointer, dimension(:)  :: maxtwm            !upper limit of temperature of the warmest month
logical, pointer, dimension(:)  :: estab              !true if pft can establish
logical, pointer, dimension(:)  :: survive            !true if pft can survive in gridcell
logical, pointer, dimension(:)  :: present            !true if pft is present in gridcell          
real(dp), pointer               :: twm                !maximum monthly temperature for the year being simulated
real(dp), pointer               :: tminm              !minimum monthly temperature for the year being simulated

! LOCAL VARIABLES
real(dp) :: yrs                 !years in climate buffer.
real(dp) :: mtemp_min20         !mean minimum coldest month air temperature over last 20 years

! point pointers
mtemp_min_buf           => sv(j)%mtemp_min_buf
gdd                     => veg%gdd_bioclim
mintcm                  => prm%mintcm         
maxtcm                  => prm%maxtcm         
gddmin                  => prm%gddmin         
maxtwm                  => prm%maxtwm  
twm                     => atm%twm       
estab                   => veg%pft_estab    
survive                 => veg%pft_survive  
present                 => sv(j)%present  
tminm                   => atm%tminm    

!---------

! begin calculations

if (yr > climbuf) then
  yrs = climbuf
else
  yrs = yr
end if

! The survival and establishment criteria only needs to be calculated once per year
! (since the information is used only in annual calculations)
 
mtemp_min_buf = eoshift(mtemp_min_buf,1,tminm)
mtemp_min20 = sum(mtemp_min_buf) / yrs

!  Limits based on 20-year running averages of coldest-month mean
!  temperature and growing degree days (5 degree base).
!  For SURVIVAL, coldest month temperature and GDD should be
!  at least as high as PFT-specific limits.
!  For REGENERATION, PFT must be able to survive AND coldest month
!  temperature should be no higher than a PFT-specific limit.

survive = .false.
estab   = .false.

where (mtemp_min20 >= mintcm) survive = .true.
  
where (gdd >= gddmin .and. mtemp_min20 <= maxtcm .and. twm <= maxtwm) estab = .true.

end subroutine bioclimatic

!-------------------------------------------------------------------------------------------------------------

subroutine reproduction(j)

! SUBROUTINE REPRODUCTION (annual C allocation scheme)
! Deduction of reproduction costs from annual biomass increment
! recoded from f to f.90 by Joe Melton, 2007 from orginal LPJ code

use arveparams,         only : dp,npft,ncvar,reprod_cost
use statevars,          only : sv,co2
use pftparametersmod,   only : prm

implicit none

!ARGUMENTS
integer, intent(in) :: j        !gridcell index

real(dp), pointer, dimension(:,:)        :: litter_ag_fast      !above ground litter pool fast
real(dp), pointer, dimension(:,:)        :: litter_ag_slow      !above ground litter pool slow
logical, pointer, dimension(:)           :: present             !true if pft is present
logical, pointer, dimension(:)           :: tree                !true if pft is tree
logical, pointer, dimension(:)           :: c4                  !true if pft is c4
real(dp), pointer, dimension(:,:)        :: lm_sapl             !initial (sapling) leaf mass (gC/m2)
real(dp), pointer, dimension(:,:)        :: hm_sapl             !initial (sapling) heartwood mass (gC/m2)
real(dp), pointer, dimension(:,:)        :: rm_sapl             !initial (sapling) fine root mass (gC/m2)
real(dp), pointer, dimension(:,:)        :: sm_sapl             !initial (sapling) sapwood mass (gC/m2)
real(dp), pointer, dimension(:,:)        :: bm_inc              !annual biomass increment (gC/m2)
real(dp), pointer, dimension(:)          :: reprod              !allocation to reproduction (gC/m2)

!LOCAL VARIABLES
integer  :: pft,c             !pft counter
real(dp) :: temp              !temporary variable

!point pointers
litter_ag_fast => sv(j)%litter_ag_fast
litter_ag_slow => sv(j)%litter_ag_slow
present        => sv(j)%present
lm_sapl        => veg%lm_sapl
hm_sapl        => veg%hm_sapl
rm_sapl        => veg%rm_sapl
sm_sapl        => veg%sm_sapl
bm_inc         => veg%bm_inc
tree           => prm%tree
c4             => prm%c4
reprod         => veg%reprod

!------------------------

! begin calculations

do pft = 1,npft 

  if (present(pft)) then  

    ! Calculate allocation to reproduction
    ! Reproduction costs taken simply as a constant fraction of annual NPP

    reprod = max(bm_inc(pft,1) * reprod_cost,0._dp)

    ! assume the costs go to reproductive structures which will
    ! eventually enter the litter pool

    temp = litter_ag_fast(pft,1)
    litter_ag_fast(pft,1) = litter_ag_fast(pft,1) + reprod(pft)

    if (litter_ag_fast(pft,1) > 1._dp) then !JM 01.03.2011
      do c = 2,ncvar
        litter_ag_fast(pft,c) = (litter_ag_fast(pft,c) * temp + bm_inc(pft,c) * reprod(pft)) / (litter_ag_fast(pft,1))
      end do
    else
        litter_ag_fast(pft,2:ncvar) = 0._dp  
    end if

    ! Reduce biomass increment by reproductive cost

    bm_inc(pft,1) = bm_inc(pft,1) - reprod(pft)

    if (tree(pft)) then
             if (bm_inc(pft,1) > 0._dp) then

                lm_sapl(pft,2:ncvar) = bm_inc(pft,2:ncvar)
                sm_sapl(pft,2:ncvar) = bm_inc(pft,2:ncvar)
                hm_sapl(pft,2:ncvar) = bm_inc(pft,2:ncvar)
                rm_sapl(pft,2:ncvar) = bm_inc(pft,2:ncvar)

             else

                lm_sapl(pft,2) = 17.8 - co2(2) !C3 value
                sm_sapl(pft,2) = 17.8 - co2(2) !from lloyd & farquhar,1994
                hm_sapl(pft,2) = 17.8 - co2(2)
                rm_sapl(pft,2) = 17.8 - co2(2)
                lm_sapl(pft,3) = co2(3)
                sm_sapl(pft,3) = co2(3)
                hm_sapl(pft,3) = co2(3)
                rm_sapl(pft,3) = co2(3)

             end if
             
    else  ! grass
    
              if (bm_inc(pft,1) > 0._dp) then

                lm_sapl(pft,2:ncvar) = bm_inc(pft,2:ncvar)
                rm_sapl(pft,2:ncvar) = bm_inc(pft,2:ncvar)

              else
                 if (c4(pft)) then

                  lm_sapl(pft,2) = 3.6 - co2(2)  !C4 value
                  rm_sapl(pft,2) = 3.6 - co2(2)  !from lloyd & farquhar,1994
                  lm_sapl(pft,3) = co2(3)
                  rm_sapl(pft,3) = co2(3)

                else

                  lm_sapl(pft,2) = 17.8 - co2(2) !C3 value
                  rm_sapl(pft,2) = 17.8 - co2(2) !from lloyd & farquhar,1994
                  lm_sapl(pft,3) = co2(3)
                  rm_sapl(pft,3) = co2(3)

                   end if
                    end if
    end if !tree/grass
  end if !present
end do !pft do loop

return

end subroutine reproduction

!-------------------------------------------------------------------------------------------------------------

subroutine turnover(j)

! SUBROUTINE TURNOVER (annual C allocation scheme)
! Turnover of PFT-specific fraction from each living C pool
! Leaf and root C transferred to litter, sapwood C to heartwood
! recoded from f to f.90 by Joe Melton 2007 from original LPJ code
! NOTE: new way of calculating fine root turnover

use arveparams,         only : dp,npft,ncvar
use statevars,          only : sv
use pftparametersmod,   only : prm

implicit none

! ARGUMENTS:
integer, intent(in) :: j

logical, pointer, dimension(:)        :: present
real(dp), pointer, dimension(:,:)     :: lm_ind                 !individual leaf mass (gC)
real(dp), pointer, dimension(:,:)     :: sm_ind                 !individual sapwood mass (gC)
real(dp), pointer, dimension(:,:)     :: hm_ind                 !individual heartwood mass (gC)
real(dp), pointer, dimension(:,:)     :: rm_ind                 !individual fine root mass (gC)
real(dp), pointer, dimension(:,:)     :: litter_ag_fast         !gridcell above-ground litter (gC/m2)
real(dp), pointer, dimension(:,:)     :: litter_ag_slow         !gridcell above-ground litter (gC/m2)
real(dp), pointer, dimension(:,:)     :: litter_bg              !gridcell below-ground litter (gC/m2)
real(dp), pointer, dimension(:)       :: nind                   !gridcell individual density (indiv/m2)
real(dp), pointer, dimension(:)       :: leafturnover           !leaf turnover period (years)
real(dp), pointer, dimension(:)       :: sapturnover            !sapwood turnover period (sapwood converted to heartwood) (years)
real(dp), pointer, dimension(:)       :: frootQ10               !fine root turnover Q10 exponent (Gill and Jackson 2000)
real(dp), pointer, dimension(:)       :: turnover_ind           !total C turnover per individual (gC)

! LOCAL VARIABLES:
integer  :: pft                 !counter
real(dp) :: l_torate            !leaf turnover rate (yr-1)
real(dp) :: s_torate            !sapwood turnover rate(yr-1)
real(dp) :: r_torate            !fine root turnover rate(yr-1)
real(dp) :: lm_turn             !leaf mass turnover in one year (gC)
real(dp) :: sm_turn             !sapwood mass turnover in one year (gC)
real(dp) :: rm_turn             !fine root mass turnover in one year (gC)
real(dp) :: temp                !temporary variable

! parameters
real(dp), parameter, dimension(npft) :: rootturnover = [ 2., 1., 2., 1., 1., 2., 1., 2., 2. ] 

! point pointers
present                 => sv(j)%present
lm_ind                  => sv(j)%lm_ind
sm_ind                  => sv(j)%sm_ind
hm_ind                  => sv(j)%hm_ind
rm_ind                  => sv(j)%rm_ind
litter_ag_fast          => sv(j)%litter_ag_fast
litter_ag_slow          => sv(j)%litter_ag_slow
litter_bg               => sv(j)%litter_bg
nind                    => sv(j)%nind
turnover_ind            => veg%turnover_ind
leafturnover            => prm%leafturnover
sapturnover             => prm%sapturnover
frootQ10                => prm%frootQ10

!--------------------

! calculation begin

do pft = 1,npft
    if (present(pft)) then

    ! Turnover rates are reciprocals of tissue longevity
    l_torate = 1._dp / leafturnover(pft)
    s_torate = 1._dp / sapturnover(pft)
    r_torate = 1._dp / rootturnover(pft)

    ! Calculate the biomass turnover in this year
    lm_turn = lm_ind(pft,1) * l_torate
    sm_turn = sm_ind(pft,1) * s_torate
    rm_turn = rm_ind(pft,1) * r_torate

    ! Update the pools
    lm_ind(pft,1) = lm_ind(pft,1) - lm_turn
    sm_ind(pft,1) = sm_ind(pft,1) - sm_turn
    rm_ind(pft,1) = rm_ind(pft,1) - rm_turn

    ! Convert sapwood to heartwood
    temp = hm_ind(pft,1)
    hm_ind(pft,1) = hm_ind(pft,1) + sm_turn

    if (hm_ind(pft,1) > 0._dp) then

          hm_ind(pft,2:ncvar) = (hm_ind(pft,2:ncvar) * temp + sm_ind(pft,2:ncvar) * sm_turn) / hm_ind(pft,1)
     end if


    ! Transfer to litter pools
    temp = litter_ag_fast(pft,1)
    litter_ag_fast(pft,1) = litter_ag_fast(pft,1) + lm_turn * nind(pft)

    if (litter_ag_fast(pft,1) > 1._dp) then !JM 01.03.2011

          litter_ag_fast(pft,2:ncvar) = (litter_ag_fast(pft,2:ncvar) * temp &
                                        + lm_ind(pft,2:ncvar) * lm_turn * nind(pft)) / litter_ag_fast(pft,1)
    end if

    temp = litter_bg(pft,1)
    litter_bg(pft,1) = litter_bg(pft,1) + rm_turn * nind(pft)


    if (litter_bg(pft,1) > 1._dp) then  !JM 01.03.2011

          litter_bg(pft,2:ncvar) = (litter_bg(pft,2:ncvar) * temp + rm_ind(pft,2:ncvar) * rm_turn * nind(pft)) / litter_bg(pft,1)

    end if

    ! Record total turnover
    turnover_ind(pft) = lm_turn + sm_turn + rm_turn

    ! Set isotopes to zero if pools are empty
       if (lm_ind(pft,1) <= 0._dp) lm_ind(pft,1:ncvar) = 0._dp
       if (sm_ind(pft,1) <= 0._dp) sm_ind(pft,1:ncvar) = 0._dp
       if (rm_ind(pft,1) <= 0._dp) rm_ind(pft,1:ncvar) = 0._dp

  end if !present
end do  !pft

return

end subroutine turnover

!-------------------------------------------------------------------------------------------------------------

subroutine killplant(j)

! SUBROUTINE killplant
! Removal of PFTs with negative annual C increment or with an awscal of zero.
! recoded from f to f.90 by Joe Melton 2007 from original LPJ code
! Note: Killing of PFTs newly beyond their bioclimatic limits is done in subroutine establishment
! awscal requirement added by JM (17.11.2010) as bioclimatic no longer has the 
! knowledge of the entire years climate before the year begins so can allow plants to start to
! grow based upon a good season, but that good season would have passed and the next could be 
! too dry for growth.

use arveparams,         only : dp,npft,ncvar
use pftparametersmod,   only : prm
use statevars,          only : sv

implicit none

! ARGUMENTS
integer, intent(in) :: j                        !grid cell

! pointers
real(dp), pointer, dimension(:,:) :: lm_ind                !individual leaf mass (gC)
real(dp), pointer, dimension(:,:) :: sm_ind                !individual sapwood mass (gC)
real(dp), pointer, dimension(:,:) :: hm_ind                !individual heartwood mass (gC)
real(dp), pointer, dimension(:,:) :: rm_ind                !individual fine root mass (gC)
real(dp), pointer, dimension(:,:) :: litter_ag_fast        !gridcell above-ground litter (gC/m2)
real(dp), pointer, dimension(:,:) :: litter_ag_slow        !gridcell above-ground litter (gC/m2)
real(dp), pointer, dimension(:,:) :: litter_bg             !gridcell below-ground litter (gC/m2)
real(dp), pointer, dimension(:)   :: nind                  !gridcell individual density (indiv/m2)
logical, pointer, dimension(:)    :: tree                  !true if pft is tree
logical, pointer, dimension(:)    :: present               !true if pft is present in that gridcell
real(dp), pointer, dimension(:,:) :: bm_inc                !annual biomass increment (gC/m2)
real(dp), pointer, dimension(:)   :: awscal                !annual fractional water scalar

! LOCAL VARIABLES
integer  :: pft,c                       !counters
real(dp) :: temp                        !temporary variable
real(dp) :: litter_inc                  !temporary variable, litter increase

! point pointers
present                 => sv(j)%present
lm_ind                  => sv(j)%lm_ind
sm_ind                  => sv(j)%sm_ind
hm_ind                  => sv(j)%hm_ind
rm_ind                  => sv(j)%rm_ind
litter_ag_fast          => sv(j)%litter_ag_fast
litter_ag_slow          => sv(j)%litter_ag_slow
litter_bg               => sv(j)%litter_bg
nind                    => sv(j)%nind
bm_inc                  => veg%bm_inc
awscal                  => veg%awscal
tree                    => prm%tree

!------------------------

! begin calculations

do pft = 1,npft

  if (present(pft)) then

    if (bm_inc(pft,1) < 0._dp .or. awscal(pft) == 0.d0) then  ! negative C increment this year
       
       present(pft) = .false. ! remove PFT
       
       ! Transfer killed biomass to litter

      if (tree(pft)) then

        temp = litter_ag_fast(pft,1)
        litter_ag_fast(pft,1) = litter_ag_fast(pft,1) + nind(pft) * lm_ind(pft,1)

!        if (lm_ind(pft,1) > 0._dp .and. litter_ag_fast(pft,1) > 1._dp) then  
        if (litter_ag_fast(pft,1) > 1._dp) then  !JM 04.03.2011
          do c = 2,ncvar
             litter_inc = lm_ind(pft,1) * lm_ind(pft,c)
             litter_ag_fast(pft,c) = (litter_ag_fast(pft,c) * temp + litter_inc * nind(pft)) / litter_ag_fast(pft,1)
          end do
        end if 

        temp = litter_ag_slow(pft,1)
        litter_ag_slow(pft,1) = litter_ag_slow(pft,1) + nind(pft) * (sm_ind(pft,1) + hm_ind(pft,1))

        if (litter_ag_slow(pft,1) > 1._dp) then  !JM 01.03.2011
        do c = 2,ncvar
           litter_inc = (sm_ind(pft,1) * sm_ind(pft,c) + hm_ind(pft,1) * hm_ind(pft,c))
           litter_ag_slow(pft,c) = (litter_ag_slow(pft,c) * temp + litter_inc * nind(pft)) / litter_ag_slow(pft,1)
        end do
        end if

      else  !grasses

        temp = litter_ag_fast(pft,1)
        litter_inc = lm_ind(pft,1) * lm_ind(pft,2)
        litter_ag_fast(pft,1) = litter_ag_fast(pft,1) + lm_ind(pft,1) * nind(pft)

        if (litter_ag_fast(pft,1) > 1._dp) then  !JM 01.03.2011
        do c = 2,ncvar
           litter_inc = lm_ind(pft,1) * lm_ind(pft,c)
           litter_ag_fast(pft,c) = (litter_ag_fast(pft,c) * temp + litter_inc * nind(pft)) / litter_ag_fast(pft,1)
        end do
        end if

      end if !tree/grass loop

      temp = litter_bg(pft,1) !fine roots
      litter_bg(pft,1) = litter_bg(pft,1) + rm_ind(pft,1) * nind(pft)
       
      if (litter_bg(pft,1) > 1._dp) then  !JM 01.03.2011
        do c = 2,ncvar
           litter_bg(pft,c) = (litter_bg(pft,c) * temp + rm_ind(pft,c) * nind(pft)) / litter_bg(pft,1)
        end do
      end if  

    end if  !bm_inc
  end if !present
end do !pft

return

end subroutine killplant

!-------------------------------------------------------------------------------------------------------------

subroutine allocation(j)

! SUBROUTINE ALLOCATION (annual C allocation scheme)
! Allocation of annual C increment to leaf, stem and fine root
! compartments, update of individual structure
! recoded from f to f.90 by Joe Melton 2007 from original LPJ code

use arveparams,         only : dp,npft,ncvar,pi,reinickerp,allom2,allom3
use pftparametersmod,   only : prm
use statevars,          only : sv

implicit none

! ARGUMENTS
integer, intent(in) :: j

! pointers
real(dp), pointer, dimension(:,:)       :: lm_ind                !individual leaf mass (gC)
real(dp), pointer, dimension(:,:)       :: sm_ind                !individual sapwood mass (gC)
real(dp), pointer, dimension(:,:)       :: hm_ind                !individual heartwood mass (gC)
real(dp), pointer, dimension(:,:)       :: rm_ind                !individual fine root mass (gC)
real(dp), pointer, dimension(:,:)       :: litter_ag_fast        !gridcell above-ground litter (gC/m2)
real(dp), pointer, dimension(:,:)       :: litter_ag_slow        !gridcell above-ground litter (gC/m2)
real(dp), pointer, dimension(:,:)       :: litter_bg             !gridcell below-ground litter (gC/m2)
real(dp), pointer, dimension(:)         :: nind                  !gridcell individual density (indiv/m2)
real(dp), pointer, dimension(:)         :: crownarea             !crown area (m2)
real(dp), pointer, dimension(:)         :: lai_max               !individual leaf area index (max)
real(dp), pointer, dimension(:)         :: height                !tree height (m)
logical, pointer, dimension(:)          :: tree                  !true if pft is tree
logical, pointer,  dimension(:)         :: present               !true if pft is present in that gridcell
real(dp), pointer, dimension(:)         :: maxcrowna             !tree maximum crown area (m2)
real(dp), pointer, dimension(:)         :: sla                   !PFT specific leaf area (m2/gC)
real(dp), pointer, dimension(:,:)       :: bm_inc                !annual biomass increment (gC/m2)
real(dp), pointer, dimension(:)         :: awscal                !mean daily water scalar (among leaf-on days) (0-1 scale)
real(dp), pointer, dimension(:)         :: latosa                !leaf area to sapwood area (m2 m-2)
real(dp), pointer, dimension(:)         :: allom1                !Packing constrain from Reinicke's rule exponent
real(dp), pointer, dimension(:)         :: fpc_gmax              !max gridcell foliar projective cover (m2 m-2)
real(dp), pointer, dimension(:)         :: wooddens              !density of wood (sap + heart) (g m-3)

! PARAMETERS
real(dp), parameter :: xacc = 0.1     !threshold x-axis precision of allocation soln
real(dp), parameter :: yacc = 1.0e-10 !threshold y-axis precision of allocation soln
integer, parameter  :: nseg = 20      !number of equal segments
integer, parameter  :: jmax = 40      !max number of iterations in search for alloc soln

real(dp), parameter, dimension(npft) :: leafrootratio = [ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.75, 0.75 ]!leaf to root ratio under non-water stressed conditions

! LOCAL VARIABLES
integer  :: pft,c                       !counters
real(dp) :: lminc_ind                   !individual leafmass increment this year
real(dp) :: sminc_ind                   !individual sapmass increment this year
real(dp) :: rminc_ind                   !individual fineroot mass increment this year
real(dp) :: bm_inc_ind                  !individual total biomass increment this year
real(dp) :: sap_xsa                     !cross sectional area of sapwood
real(dp) :: stem_diam                   !stem diameter
real(dp) :: x1                          !working vars in bisection
real(dp) :: x2
real(dp) :: rtbis
real(dp) :: dx
real(dp) :: xmid
real(dp) :: sign
real(dp) :: fx1
real(dp) :: fmid
real(dp) :: lminc_ind_min               !min leafmass increment to maintain current sapwood
real(dp) :: rminc_ind_min               !min rootmass increment to support new leafmass
integer  :: i                           !counter in bisection loop
real(dp) :: lmtorm                      !ratio of leafmass to fine rootmass
real(dp) :: crownarea_max               !maximum crown area (m2)
real(dp) :: temp                        !temporary variable
real(dp) :: lm_temp                     !temporary variable leaf
real(dp) :: rm_temp                     !temporary variable root
real(dp) :: sm_temp                     !temporary variable sapwood
logical  :: normal                      !type of allocation flag (normal = positive increment to all living C compartments)

! point pointers
present                 => sv(j)%present
lm_ind                  => sv(j)%lm_ind
sm_ind                  => sv(j)%sm_ind
hm_ind                  => sv(j)%hm_ind
rm_ind                  => sv(j)%rm_ind
litter_ag_fast          => sv(j)%litter_ag_fast
litter_ag_slow          => sv(j)%litter_ag_slow
litter_bg               => sv(j)%litter_bg
nind                    => sv(j)%nind
crownarea               => sv(j)%crownarea
lai_max                 => sv(j)%lai_max
height                  => sv(j)%height
sla                     => veg%sla
bm_inc                  => veg%bm_inc
awscal                  => veg%awscal
latosa                  => prm%latosa
allom1                  => prm%allom1
tree                    => prm%tree
maxcrowna               => prm%maxcrowna
fpc_gmax                => sv(j)%fpc_gmax
wooddens                => prm%wooddens

!-----------------

! calculations begin

do pft = 1,npft

 if (present(pft)) then

    bm_inc_ind = bm_inc(pft,1) / nind(pft)

    !calculate this year's leaf to fine root mass ratio from mean annual
    !water scalar and pft specific parameter

    lmtorm = leafrootratio(pft) * awscal(pft)

       if (tree(pft)) then

         normal = .false.
         crownarea_max = maxcrowna(pft)

         !TREE ALLOCATION

         !Allocation of this year's biomass increment (bm_inc_ind) to the
         !three living carbon pools, such that the basic allometric
         !relationships (A-C below) are always satisfied.

           !Eqn A below is being related by a relationship from Hickler's hydraulic routine.
         !(A) (leaf area) = latosa * (sapwood xs area)
           !    (Pipe Model, Shinozaki et al. 1964a,b; Waring et al 1982)
         !(B) (leaf mass) = lmtorm * (root mass)
         !(C) height = allom2 * (stem diameter)**allom3
         !    (eg. Huang et al. 1992 )
         !(D) (crown area) = min (allom1 * (stem diameter)**reinickerp,
         !                              crownarea_max)

         !Mathematical derivation:

           !(1) bm_inc_ind = lminc_ind + sminc_ind + rminc_ind
           !(2) leaf_area_new = latosa * sap_xsa_new   [from (A)]
           !(3) leaf_area_new = (lm_ind + lminc_ind) * sla
         !from (2) & (3),
           !(4) (lm_ind + lminc_ind) * sla = latosa * sap_xsa_new
         !from (4),
           !(5) sap_xsa_new = (lm_ind + lminc_ind) * sla / latosa
           !(6) (lm_ind + lminc_ind) = lmtorm * (rm_ind + rminc_ind)
           !      [from (B)]
           !(7) height_new = allom2 * stem_diam_new**allom3  [from (C)]
         !from (1),
           !(8) sminc_ind = bm_inc_ind - lminc_ind - rminc_ind
         !from (6),
           !(9) rminc_ind=((lm_ind + lminc_ind) / lmtorm) - rm_ind
         !from (8) & (9),
          !(10) sminc_ind = bm_inc_ind - lminc_ind
           !      - ((lm_ind + lminc_ind)  / lmtorm) + rm_ind
          !(11) wooddens = (sm_ind + sminc_ind + hm_ind) / stemvolume_new
          !(12) stemvolume_new = height_new * pi * stem_diam_new**2 / 4
         !from (10), (11) & (12)
          !(13) stem_diam_new = [ ((sm_ind + bm_inc_ind - lminc_ind
           !      - ((lm_ind + lminc_ind) / lmtorm) + rm_ind + hm_ind)
           !      / wooddens) / (height_new * pi / 4) ]**(1/2)
         !combining (7) and (13),
          !(14) height_new = allom2 * [ ((sm_ind + bm_inc_ind - lminc_ind
           !      - ((lm_ind + lminc_ind) / lmtorm) + rm_ind + hm_ind)
           !     / wooddens) / (height_new * pi / 4) ]**(1/2 * allom3)
         !from (14),
          !(15) height_new**(1 + 2 / allom3) = allom2**(2 / allom3)
           !      * ((sm_ind + bm_inc_ind - lminc_ind - ((lm_ind + lminc_ind)
            !     / lmtorm) + rm_ind + hm_ind) / wooddens) / (pi / 4)
          !(16) wooddens = (sm_ind + sminc_ind) / sapvolume_new
         !from (10) and (16),
          !(17) wooddens = (sm_ind + bm_inc_ind - lminc_ind
           !      - ((lm_ind + lminc_ind) / lmtorm) + rm_ind) / sapvolume_new
          !(18) sapvolume_new = height_new * sap_xsa_new
         !from (17) and (18),
          !(19) sap_xsa_new = (sm_ind + bm_inc_ind - lminc_ind
           !      - ((lm_ind + lminc_ind) / lmtorm) + rm_ind)
           !      / (height_new * wooddens)
         !from (19),
          !(20) height_new = (sm_ind + bm_inc_ind - lminc_ind
           !      - ((lm_ind + lminc_ind) / lmtorm) + rm_ind )
           !      / (sap_xsa_new * wooddens)
         !from (5) and (20),
          !(21) height_new**(1 + 2 / allom3) = [ (sm_ind + bm_inc_ind
           !      - lminc_ind - ((lm_ind + lminc_ind) / lmtorm) + rm_ind )
           !     / ((lm_ind + lminc_ind) * sla * wooddens / latosa) ]
           !      **(1 + 2 / allom3)
        ! -------------------------------------------------------------------
          !(15) and (21) are two alternative expressions for
           !    height_new**(1 + 2 / allom3). Combining these,

          !(22) allom2**(2 / allom3) * ((sm_ind + bm_inc_ind - lminc_ind
           !      - ((lm_ind + lminc_ind) / lmtorm) + rm_ind + hm_ind)
           !      / wooddens) / (pi / 4) - [ (sm_ind + bm_inc_ind - lminc_ind
           !      - ((lm_ind + lminc_ind) / lmtorm) + rm_ind )
           !      / ((lm_ind + lminc_ind) * sla * wooddens / latosa) ]
           !      **(1 + 2 / allom3)
           !      = 0

         !Equation (22) can be expressed in the form f(lminc_ind)=0.

         !Numerical methods are used to solve the equation for the
         !unknown lminc_ind.
         !-------------------------------------------------------------------

         !Work out minimum leaf production to maintain current sapmass

         ! (23) sap_xsa = sm_ind / wooddens / height
         !from (A) and (23),
         ! (24) leaf_mass * sla = latosa * sap_mass / wooddens / height
         !from (24),
         ! (25) leaf_mass = latosa * sap_mass / (wooddens * height * sla)
         !from (25), assuming sminc_ind=0,
         ! (26) lm_ind + lminc_ind_min = latosa * sm_ind
         !      / (wooddens * height * sla)
         !from (26),
         ! (27) lminc_ind_min = latosa * sm_ind / (wooddens * height * sla)
         !              - lm_ind

              !write(0,*)'alloc0hak',wooddens,height(pft),sla(pft)

         lminc_ind_min = latosa(pft) * sm_ind(pft,1) / (wooddens(pft) * height(pft) * sla(pft)) - lm_ind(pft,1)  !eqn (27)

         !Work out minimum root production to support this leaf mass
         !(i.e. lm_ind + lminc_ind_min)
         !May be negative following a reduction in soil water limitation
         !(increase in lmtorm) relative to last year.

         !from (B) and (25),
         ! (28) root_mass = latosa * sap_mass / (wooddens * height * sla)
         !              / lmtorm
         !from (28), assuming sminc_ind=0,
         ! (29) rm_ind + rminc_ind_min = latosa * sm_ind
         !      / (wooddens * height * sla * lmtorm)
         !from (29),
         ! (30) rminc_ind_min = latosa * sm_ind
         !              / (wooddens * height * sla * lmtorm) - rm_ind

         rminc_ind_min = latosa(pft) * sm_ind(pft,1) / (wooddens(pft) * height(pft) * sla(pft) * lmtorm) - rm_ind(pft,1)!eqn (30)

         if (rminc_ind_min > 0._dp .and. lminc_ind_min > 0._dp .and. rminc_ind_min + lminc_ind_min <= bm_inc_ind) then

           normal = .true.

           !Normal allocation (positive increment to all living C
           !compartments)

           !Calculation of leaf mass increment (lminc_ind) that satisfies
           !Eqn (22) using Bisection Method (Press et al 1986, p 346)

           !Seeking a root for non-negative lminc_ind, rminc_ind and
           !sminc_ind.  There should be exactly one (no proof presented, but
           !Steve has managed one) and it should lie between x1=0 and
           !x2=(bm_inc_ind-(lm_ind/lmtorm-rm_ind))/(1+1/lmtorm).

           x1 = 0._dp
           x2 = (bm_inc_ind - (lm_ind(pft,1) / lmtorm - rm_ind(pft,1))) / (1._dp + 1._dp / lmtorm)
           dx = (x2 - x1) / real(nseg)

           if (lm_ind(pft,1) == 0._dp) x1 = x1 + dx  !to avoid division by zero

           !evaluate f(x1)=LHS of eqn (22) at x1

           fx1=allom2**(2.0/allom3)*((sm_ind(pft,1)+bm_inc_ind-x1 &
               -((lm_ind(pft,1)+x1)/lmtorm)+rm_ind(pft,1) &
               +hm_ind(pft,1))/wooddens(pft))/(pi/4.0)-((sm_ind(pft,1) &
               +bm_inc_ind-x1-((lm_ind(pft,1)+x1)/lmtorm) &
               +rm_ind(pft,1))/((lm_ind(pft,1)+x1)*sla(pft) &
               *wooddens(pft)/latosa(pft)))**(1.0+2.0/allom3)

           !Find approximate location of leftmost root on the interval
           !(x1,x2).  Subdivide (x1,x2) into nseg equal segments seeking
           !change in sign of f(xmid) relative to f(x1).

           fmid = fx1
           xmid = x1

           do while (fmid * fx1 > 0._dp .and. xmid < x2)

             xmid = xmid + dx
                 fmid=allom2**(2.0/allom3)*((sm_ind(pft,1) &
                 +bm_inc_ind-xmid-((lm_ind(pft,1)+xmid)/lmtorm) &
                 +rm_ind(pft,1)+hm_ind(pft,1)) &
                 /wooddens(pft))/(pi/4.0)-((sm_ind(pft,1)+bm_inc_ind-xmid &
                 -((lm_ind(pft,1)+xmid)/lmtorm)+rm_ind(pft,1)) &
                 /((lm_ind(pft,1)+xmid)*sla(pft)*wooddens(pft)/latosa(pft))) &
                 **(1.0+2.0/allom3)
           end do

           x1 = xmid - dx
           x2 = xmid

           !Apply bisection method to find root on the new interval (x1,x2)

                fx1=allom2**(2.0/allom3)*((sm_ind(pft,1)+bm_inc_ind-x1 &
               -((lm_ind(pft,1)+x1)/lmtorm)+rm_ind(pft,1) &
               +hm_ind(pft,1))/wooddens(pft))/(pi/4.0)-((sm_ind(pft,1) &
               +bm_inc_ind-x1-((lm_ind(pft,1)+x1)/lmtorm) &
               +rm_ind(pft,1))/((lm_ind(pft,1)+x1)*sla(pft)*wooddens(pft) &
               /latosa(pft)))**(1.0+2.0/allom3)

           if (fx1 >= 0._dp) then
             sign = -1._dp
           else
             sign = 1._dp
           end if

           rtbis = x1
           dx = x2 - x1

           !Bisection loop
           !Search iterates on value of xmid until xmid lies within
           !xacc of the root, i.e. until |xmid-x|<xacc where f(x)=0

           fmid = 1._dp  !dummy value to guarantee entry to loop
           i = 0._dp          !number of iterations so far (maximum tries=jmax)

           do while (dx >= xacc .and. abs(fmid) > yacc)

             dx = dx * 0.5
             xmid = rtbis + dx

             !calculate fmid=f(xmid) [eqn (22)]

              fmid=allom2**(2.0/allom3)*((sm_ind(pft,1)+bm_inc_ind&
                 -xmid-((lm_ind(pft,1)+xmid)/lmtorm)+rm_ind(pft,1)&
                 +hm_ind(pft,1))/wooddens(pft))/(pi/4.0)-((sm_ind(pft,1)&
                 +bm_inc_ind-xmid-((lm_ind(pft,1)+xmid)/lmtorm)&
                 +rm_ind(pft,1))/((lm_ind(pft,1)+xmid)*sla(pft)&
                 *wooddens(pft)/latosa(pft)))**(1.0+2.0/allom3)

             if (fmid * sign <= 0._dp) rtbis = xmid
             i = i + 1

           end do

           !Now rtbis contains numerical solution for lminc_ind given
           !eqn (22)

           lminc_ind = rtbis

           !Calculate increments in other compartments using allometry
           !relationships

           rminc_ind = (lm_ind(pft,1) + lminc_ind) / lmtorm - rm_ind(pft,1)!eqn (9)
           sminc_ind = bm_inc_ind - lminc_ind - rminc_ind    !eqn (1)

         else

           !Abnormal allocation: reduction in some C compartment(s)
           !to satisfy allometry

           !Attempt to distribute this year's production among leaves and
           !roots only
           !1 (31) bm_inc_ind = lminc_ind + rminc_ind
           !from (31) and (9),
           ! (32) bm_inc_ind = lminc_ind + ((lm_ind + lminc_ind) / lmtorm)
           !        - rm_ind
           !from (32)
           ! (33) lminc_ind = (bm_inc_ind - lmind / lmtorm + rm_ind) /
           !        (1 + 1 / lmtorm)

           lminc_ind = (bm_inc_ind - lm_ind(pft,1) / lmtorm + rm_ind(pft,1)) / (1._dp + 1._dp / lmtorm)!eqn (33)

           if (lminc_ind >= 0._dp) then

            ! Positive allocation to leafmass

             rminc_ind = bm_inc_ind - lminc_ind  !eqn (31)

             !Add killed roots (if any) to below-ground litter

                     if (rminc_ind < 0._dp) then

                     !Trying a different approach because negative rminc_ind does not seem to be reasonable.
                     !Idea: Max lminc_ind cannot be bigger than bm_inc_ind

                      !   lminc_ind=bm_inc_ind
                      !   rminc_ind=0.0

                     !THE FOLLOWING WAS ORIGINAL LPJ CODE

                        lminc_ind = bm_inc_ind
                        rminc_ind = ((lm_ind(pft,1) + lminc_ind) / lmtorm) - rm_ind(pft,1)

                      !UNTIL HERE

                        if (rminc_ind < 0._dp) then
                           temp = litter_bg(pft,1)
                           litter_bg(pft,1) = litter_bg(pft,1) + (-rminc_ind) * nind(pft)

                           if (litter_bg(pft,1) > 1._dp) then !JM 01.03.2011
                           do c = 2,ncvar
                              litter_bg(pft,c) = (litter_bg(pft,c) * temp + (-rminc_ind) * rm_ind(pft,c) * nind(pft))/ &
                                                  litter_bg(pft,1)
                           end do
                           end if

                        elseif (rminc_ind >= 0._dp) then
                           temp = litter_bg(pft,1)
                           litter_bg(pft,1) = litter_bg(pft,1) + (rminc_ind) * nind(pft)
        
                           if (litter_bg(pft,1) > 1._dp) then  !JM 01.03.2011
                           do c = 2,ncvar
                              litter_bg(pft,c) = (litter_bg(pft,c) * temp + (rminc_ind) * rm_ind(pft,c) * nind(pft)) / &
                                                   litter_bg(pft,1)

                           end do
                           end if
                          end if
                       end if
                   else

                     !Negative allocation to leaf mass

                     rminc_ind = bm_inc_ind
                     lminc_ind =(rm_ind(pft,1) + rminc_ind) * lmtorm - lm_ind(pft,1)!from eqn (9)

                     !Add killed leaves to litter

                     temp = litter_ag_fast(pft,1)
                     litter_ag_fast(pft,1) = litter_ag_fast(pft,1) + (-lminc_ind) * nind(pft)

                     if (litter_ag_fast(pft,1) > 1._dp) then  !JM 01.03.2011
                     do c = 2,ncvar
                        litter_ag_fast(pft,c) = (litter_ag_fast(pft,c) * temp + (-lminc_ind) * lm_ind(pft,c) * nind(pft))&
                                                 / litter_ag_fast(pft,1)
                     end do
                     end if

                   end if  !rminc_ind

           !Calculate sminc_ind (must be negative)

           !from (25),
           ! (34) lm_ind + lminc_ind = latosa * (sm_ind + sminc_ind)
           !        / (wooddens * height * sla)
           !from (34),
           ! (35) sminc_ind = (lm_ind + lminc_ind) * wooddens * height * sla
           !        / latosa - sm_ind

           sminc_ind = (lm_ind(pft,1) + lminc_ind) * wooddens(pft) * height(pft) * sla(pft) / latosa(pft) - sm_ind(pft,1)!eqn (35)

           !Convert killed sapwood to heartwood

           temp = hm_ind(pft,1)
           hm_ind(pft,1) = hm_ind(pft,1) + (-sminc_ind)

           do c = 2,ncvar
              hm_ind(pft,c) = (hm_ind(pft,c) * temp + (-sminc_ind) * sm_ind(pft,c)) / hm_ind(pft,1)
           end do

         end if  !lminc_ind

         !Increment C compartments

         lm_temp = lm_ind(pft,1)
         rm_temp = rm_ind(pft,1)
         sm_temp = sm_ind(pft,1)

         lm_ind(pft,1) = lm_ind(pft,1) + lminc_ind
         rm_ind(pft,1) = rm_ind(pft,1) + rminc_ind
         sm_ind(pft,1) = sm_ind(pft,1) + sminc_ind

                 if (normal) then

                    do c = 2,ncvar
                       lm_ind(pft,c) = (lm_ind(pft,c) * lm_temp + lminc_ind * bm_inc(pft,c)) / lm_ind(pft,1)
                       rm_ind(pft,c) = (rm_ind(pft,c) * rm_temp + rminc_ind * bm_inc(pft,c)) / rm_ind(pft,1)
                       sm_ind(pft,c) = (sm_ind(pft,c) * sm_temp + sminc_ind * bm_inc(pft,c)) / sm_ind(pft,1)
                    end do

                 else

                   if (lminc_ind > 0._dp) then
                      do c = 2,ncvar
                         lm_ind(pft,c) = (lm_ind(pft,c) * lm_temp + lminc_ind * bm_inc(pft,c)) / lm_ind(pft,1)
                      end do
                   end if

                   if (rminc_ind > 0._dp) then
                      do c = 2,ncvar
                         rm_ind(pft,c) = (rm_ind(pft,c) * rm_temp + rminc_ind * bm_inc(pft,c)) / rm_ind(pft,1)
                      end do
                   end if

                 end if

         !Calculate new height, diameter and crown area

         sap_xsa = lm_ind(pft,1) * sla(pft) / latosa(pft)  !eqn (5)
         height(pft) = sm_ind(pft,1) / sap_xsa / wooddens(pft)
         stem_diam = (height(pft) / allom2)**(1._dp / allom3)  !eqn (C)
         crownarea(pft) = min(allom1(pft) * stem_diam**reinickerp,crownarea_max)  !eqn (D)

       else !grass

         !GRASS ALLOCATION
         !Distribute this year's production among leaves and fine roots
         !according to leaf to rootmass ratio [eqn (33)]
         !Relocation of C from one compartment to the other not allowed:
         !negative increment in either compartment transferred to litter

         lminc_ind = (bm_inc_ind - lm_ind(pft,1) / lmtorm + rm_ind(pft,1)) / (1._dp + 1._dp / lmtorm)
         rminc_ind = bm_inc_ind - lminc_ind

                 if (lminc_ind >= 0._dp) then

                   !Add killed roots (if any) to below-ground litter

                   if (rminc_ind < 0._dp) then
                     temp = litter_bg(pft,1)
                     litter_bg(pft,1) = litter_bg(pft,1) + (-rminc_ind) *nind(pft)

                     if (litter_bg(pft,1) > 1._dp) then !JM 01.03.2011
                     do c = 2,ncvar
                       litter_bg(pft,c) = (litter_bg(pft,c) * temp + (-rminc_ind) * rm_ind(pft,c) * nind(pft)) / litter_bg(pft,1)
                     end do
                     end if
                   end if

                 else

                   !Negative allocation to leaf mass

                   rminc_ind = bm_inc_ind
                   lminc_ind =(rm_ind(pft,1) + rminc_ind) * lmtorm - lm_ind(pft,1)!from eqn (9)

                   !Add killed leaves to litter

                   temp = litter_ag_fast(pft,1)
                   litter_ag_fast(pft,1) = litter_ag_fast(pft,1) + (-lminc_ind) * nind(pft)
                   
                   if (litter_ag_fast(pft,1) > 1._dp) then !JM 01.03.2011
                    do c = 2,ncvar
                     litter_ag_fast(pft,c) = (litter_ag_fast(pft,c) * temp + (-lminc_ind) * lm_ind(pft,c) * nind(pft)) / litter_ag_fast(pft,1)
                    end do
                   end if

                 end if !lmin_ind

      !Increment C compartments
      lm_temp = lm_ind(pft,1)
      rm_temp = rm_ind(pft,1)

      lm_ind(pft,1) = lm_ind(pft,1) + lminc_ind
      rm_ind(pft,1) = rm_ind(pft,1) + rminc_ind

              if (lminc_ind >= 0._dp) then
                 do c = 2,ncvar
                    lm_ind(pft,c) = (lm_ind(pft,c) * lm_temp + lminc_ind * bm_inc(pft,c)) / lm_ind(pft,1)
                 end do
              end if

              if (rminc_ind >= 0._dp) then
                 do c = 2,ncvar
                    rm_ind(pft,c) = (rm_ind(pft,c) * rm_temp + rminc_ind * bm_inc(pft,c)) / rm_ind(pft,1)
                 end do
              end if

    end if !grass vs/ tree


    !Update LAI and FPC

     if (crownarea(pft) > 0._dp) then
       lai_max(pft) = (lm_ind(pft,1) * sla(pft)) / crownarea(pft)
       fpc_gmax(pft) = crownarea(pft) * nind(pft) * (1._dp - exp(-0.5 * lai_max(pft)))
     else
       lai_max(pft) = 0._dp             
       fpc_gmax(pft) = 0._dp
     end if

  end if !pft present


end do !pft

return

end subroutine allocation

!-------------------------------------------------------------------------------------------------------------

!SUBROUTINE LIGHT
!Competition for light among PFTs
!recoded from f to f.90 by Joe Melton 2007 from original LPJ code
!changes to tree competition to fix bug identified by J.Kaplan

subroutine light(j)

use arveparams,only : dp,npft,ncvar
use pftparametersmod, only : prm
use statevars, only : sv

implicit none

!PARAMETERS
real(dp), parameter :: fpc_tree_max = 0.95                     !maximum total tree FPC

!ARGUMENTS
integer, intent(in) :: j

real(dp), pointer, dimension(:,:) :: lm_ind                !individual leaf mass (gC)
real(dp), pointer, dimension(:,:) :: lm_max
real(dp), pointer, dimension(:,:) :: sm_ind                !individual sapwood mass (gC)
real(dp), pointer, dimension(:,:) :: hm_ind                !individual heartwood mass (gC)
real(dp), pointer, dimension(:,:) :: rm_ind                !individual fine root mass (gC)
real(dp), pointer, dimension(:,:) :: litter_ag_fast        !gridcell above-ground litter (gC/m2)
real(dp), pointer, dimension(:,:) :: litter_ag_slow        !gridcell above-ground litter (gC/m2)
real(dp), pointer, dimension(:,:) :: litter_bg             !gridcell below-ground litter (gC/m2)
real(dp), pointer, dimension(:)   :: nind                  !gridcell individual density (indiv/m2)
real(dp), pointer, dimension(:)   :: crownarea             !crown area (m2)
real(dp), pointer, dimension(:)   :: fpc_gmax              !maximal gridcell foliar projective cover (m2 m-2)
real(dp), pointer, dimension(:)   :: lai_max
logical, pointer, dimension(:)    :: tree                  !true if pft is tree
logical, pointer, dimension(:)    :: present               !true if pft is present in that gridcell
real(dp), pointer, dimension(:)   :: sla                   !PFT specific leaf area (m2/gC)

!LOCAL VARIABLES
integer  :: pft,c                                !counters
integer  :: ntree                                !no of tree PFTs currently present
real(dp) :: fpc_tree_total                       !total grid FPC for tree PFTs
real(dp) :: fpc_grass_total                      !total grid FPC for grass PFTs
real(dp) :: grasscover                           !grass PFT proportional cover ("crown area")
real(dp) :: excess                               !tree FPC or grass cover to be reduced
real(dp) :: nind_kill                            !reduction in individual density to reduce tree FPC to permitted maximum (indiv/m2)
real(dp) :: lm_kill                              !reduction in grass PFT leaf mass to reduce grass cover to permitted maximum (gC)
real(dp) :: rm_kill                              !reduction in grass PFT root mass to reduce grass cover to permitted maximum (gC)
real(dp) :: lm_old                               !
integer :: ngrass                                !
real(dp) :: temp                                 !temporary variable
real(dp) :: litter_inc                           !
real(dp), dimension(npft) :: pft_excess          !
real(dp) :: fpc_grass_max                        !max allowed grass FPC given the current tree cover
real(dp) :: lost_lm_inst                         !leaf matter that is passed to litter pools

!point pointers
lm_ind                => sv(j)%lm_ind
lm_max                => sv(j)%lm_max
sm_ind                => sv(j)%sm_ind
hm_ind                => sv(j)%hm_ind
rm_ind                => sv(j)%rm_ind
litter_ag_fast        => sv(j)%litter_ag_fast
litter_ag_slow        => sv(j)%litter_ag_slow
litter_bg             => sv(j)%litter_bg
nind                  => sv(j)%nind
crownarea             => sv(j)%crownarea
lai_max               => sv(j)%lai_max
present               => sv(j)%present
sla                   => veg%sla
fpc_gmax              => sv(j)%fpc_gmax
tree                  => prm%tree

!-------------------------------------------------------------------------------
! begin calculations

! initial presets
pft_excess = 0._dp
lm_old = 0._dp

if (do_allocation_daily) then

    ! Update the maximal FPC and LAI for plant for growth in the previous year

        where (crownarea > 0._dp)

          lai_max = lm_max(:,1) * sla / crownarea
          fpc_gmax = crownarea * nind * (1._dp - exp(-0.5 * lai_max))

        elsewhere

          lai_max = 0._dp
          fpc_gmax = 0._dp

        end where
        
else  ! annual C allocation scheme

       where (crownarea > 0.) 

         lai_max  = lm_ind(:,1) * sla / crownarea
         fpc_gmax = crownarea * nind * (1._dp - exp(-0.5 * lai_max))
        
       elsewhere

         lai_max  = 0.
         fpc_gmax = 0.
        
       end where
       
end if

        ! Calculate total woody FPC, FPC increment and grass cover (= crown area) !FLAG check on this.
        ntree           = count(present .and. tree)
        fpc_tree_total  = sum(fpc_gmax, mask = present .and. tree)
        ngrass          = count(present .and. .not. tree)
        grasscover      = sum(crownarea, mask = present .and. .not. tree)
        fpc_grass_total = sum(fpc_gmax,  mask = present .and. .not. tree)

!LIGHT COMPETITION

    if (fpc_tree_total > fpc_tree_max) then    ! reduce tree cover

    excess = fpc_tree_total - fpc_tree_max

        do pft = 1,npft

          if (present(pft) .and. tree(pft) .and. fpc_gmax(pft) > 0.) then

            !this formulation ensures equal competition (precludes total dominance by one PFT)

            pft_excess(pft) = min(fpc_gmax(pft),excess *  fpc_gmax(pft) / fpc_tree_total)

            !original LPJ formulation allows one PFT to become dominant if it has no fpc_inc (so the others are reduced)

            !if (fpc_inc_tree > 0.) then
            !  pft_excess(pft) = min(fpc_grid(pft),excess * (fpc_inc(pft) / fpc_inc_tree))
            !else
            !  pft_excess(pft) = min(fpc_grid(pft),excess / real(ntree))
            !end if

          else
            pft_excess(pft) = 0._dp
          end if

          if (pft_excess(pft) > 0._dp) then

            !use individual density (and thereby gridcell-level biomass)
            !that total tree FPC reduced to 'fpc_tree_max'

             nind_kill = nind(pft) * pft_excess(pft) / fpc_gmax(pft)
             nind(pft) = nind(pft) - nind_kill

            !Transfer lost biomass to litter

             temp = litter_ag_fast(pft,1) !leaves
             litter_ag_fast(pft,1) = litter_ag_fast(pft,1) + nind_kill * lm_ind(pft,1)  !FLAG for do_alloca_daily should this be changed to reduce lm_max???

             if (litter_ag_fast(pft,1) > 1._dp) then  !JM 01.03.2011
             do c = 2,ncvar
                litter_inc = lm_ind(pft,1) * lm_ind(pft,c)
                litter_ag_fast(pft,c) = (litter_ag_fast(pft,c) * temp + litter_inc * nind_kill) / litter_ag_fast(pft,1)
             end do
             end if

             temp = litter_ag_slow(pft,1)  !stems
             litter_ag_slow(pft,1) = litter_ag_slow(pft,1) + nind_kill * (sm_ind(pft,1) + hm_ind(pft,1))

             if (litter_ag_slow(pft,1) > 1._dp) then  !JM 01.03.2011
             do c = 2,ncvar
                litter_inc = (sm_ind(pft,1) * sm_ind(pft,c) + hm_ind(pft,1) * hm_ind(pft,c))
                litter_ag_slow(pft,c) = (litter_ag_slow(pft,c) * temp + litter_inc * nind_kill) / litter_ag_slow(pft,1)
             end do
             end if 

            temp = litter_bg(pft,1)  !fine roots
            litter_bg(pft,1) = litter_bg(pft,1) + nind_kill * rm_ind(pft,1)

            if (litter_bg(pft,1) > 1._dp) then  !JM 01.03.2011
            do c = 2,ncvar
               litter_inc = rm_ind(pft,1) * rm_ind(pft,c)
               litter_bg(pft,c) = (litter_bg(pft,c) * temp+nind_kill * litter_inc) / litter_bg(pft,1)
            end do
            end if

          end if  !pft_excess
        end do  !pft loop

   end if

    ! grass competition
    
   fpc_grass_max = 1. - min(fpc_tree_total,fpc_tree_max)

   if (fpc_grass_total > fpc_grass_max) then  !reduce grass cover

      excess = fpc_grass_total - fpc_grass_max

      if (do_allocation_daily) then  !daily C allocation scheme

        do pft = 1,npft

          if (present(pft) .and. .not. tree(pft) .and. lm_max(pft,1) > 0.) then
          
               excess = (min(fpc_tree_total,fpc_tree_max) + fpc_grass_total - 1._dp) * (fpc_gmax(pft) / fpc_grass_total)
            
               lm_old = lm_max(pft,1)
               lm_max(pft,1) = max(-2._dp * LOG(1._dp - (fpc_gmax(pft) - excess)) / sla(pft),0._dp)

               lm_kill = lm_old - lm_max(pft,1)
               rm_kill = rm_ind(pft,1) * (lm_kill / lm_old)
               rm_ind(pft,1) = rm_ind(pft,1) - rm_kill

              ! Transfer lost biomass to litter

              temp = litter_ag_fast(pft,1)  !leaves
              lost_lm_inst = max(0._dp,lm_ind(pft,1) - lm_max(pft,1))  ! if the present lm is larger than the new lm_max, this 
                                                                       ! leaf mass can be passed to litter pool
              litter_ag_fast(pft,1) = litter_ag_fast(pft,1) + lost_lm_inst
              
              do c = 2,ncvar
                litter_ag_fast(pft,c) = (litter_ag_fast(pft,c) * temp + lm_ind(pft,c) * lost_lm_inst) / litter_ag_fast(pft,1)
              end do
            
              temp = litter_bg(pft,1)  !roots
              litter_bg(pft,1) = litter_bg(pft,1) + rm_kill

              do c = 2,ncvar
                 litter_bg(pft,c) = (litter_bg(pft,c) * temp + rm_ind(pft,c) * rm_kill) / litter_bg(pft,1)
              end do
              
          end if !reduce
         
        end do  !pft
      
      else  !annual C allocation scheme
       
        do pft = 1,npft

          if (present(pft) .and. .not. tree(pft) .and. lm_ind(pft,1) > 0.) then
          
               excess = (min(fpc_tree_total,fpc_tree_max) + fpc_grass_total - 1._dp) * (fpc_gmax(pft) / fpc_grass_total)

               lm_old = lm_ind(pft,1)
               lm_ind(pft,1) = max(-2._dp * LOG(1._dp - (fpc_gmax(pft) - excess)) / sla(pft),0._dp)

               lm_kill = lm_old - lm_ind(pft,1)
               rm_kill = rm_ind(pft,1) * (lm_kill / lm_old)
               rm_ind(pft,1) = rm_ind(pft,1) - rm_kill

              ! Transfer lost biomass to litter

              temp = litter_ag_fast(pft,1)  !leaves
              litter_ag_fast(pft,1) = litter_ag_fast(pft,1) + lm_kill
              do c = 2,ncvar
                litter_ag_fast(pft,c) = (litter_ag_fast(pft,c) * temp + lm_ind(pft,c) * lm_kill) / litter_ag_fast(pft,1)
              end do
              
              temp = litter_bg(pft,1)  !roots
              litter_bg(pft,1) = litter_bg(pft,1) + rm_kill

              do c = 2,ncvar
                 litter_bg(pft,c) = (litter_bg(pft,c) * temp + rm_ind(pft,c) * rm_kill) / litter_bg(pft,1)
              end do
              
          end if !reduce
         
        end do  !pft
    
      end if  !daily vs. annual
              
   end if !reduction of grasses

if (do_allocation_daily) then

    ! Update FPC for changes due to light.
    where (crownarea > 0._dp)
      fpc_gmax = crownarea * nind * (1._dp - exp(-0.5 * lai_max))
    elsewhere
      fpc_gmax = 0._dp
    end where

else !annual C alloc scheme

     ! update lai,fpc
    where (crownarea > 0._dp)

      lai_max = lm_ind(:,1) * sla / crownarea
      fpc_gmax = (1._dp - exp(-0.5 * lai_max)) * nind * crownarea
      
    elsewhere
      
      present = .false.
      lm_ind(:,1) = 0._dp
      lai_max  = 0._dp
      fpc_gmax = 0._dp
      
    end where
    
end if

return

end subroutine light

!-------------------------------------------------------------------------------------------------------------

subroutine mortality(j)

! SUBROUTINE MORTALITY
! Tree background and stress mortality
! recoded from f to f.90 by Joe Melton 2007 from original LPJ code

use arveparams,only : dp,npft,ncvar
use pftparametersmod, only : prm
use statevars, only : sv
use metvarsmod, only : atm

implicit none

! PARAMETERS
real(dp), parameter :: mort_max = 0.01          !asymptotic maximum mortality rate [year)
real(dp), parameter :: k_mort = 0.3             !coefficient of growth efficiency in mortality equation
real(dp), parameter :: ramp_gddtw = 300.0       !ramp for heat damage function.
                                                !Above 200 growing degree days above the upper limit
                                                !tw, establishment is zero and mortality 100%

! ARGUMENTS
integer, intent(in) :: j

! pointers
real(dp), pointer, dimension(:,:) :: lm_max                !annual maximal individual leaf mass (gC)
real(dp), pointer, dimension(:,:) :: lm_ind                !individual leaf mass (gC)
real(dp), pointer, dimension(:,:) :: sm_ind                !individual sapwood mass (gC)
real(dp), pointer, dimension(:,:) :: hm_ind                !individual heartwood mass (gC)
real(dp), pointer, dimension(:,:) :: rm_ind                !individual fine root mass (gC)
real(dp), pointer, dimension(:,:) :: litter_ag_fast        !gridcell above-ground litter (gC/m2)
real(dp), pointer, dimension(:,:) :: litter_ag_slow        !gridcell above-ground litter (gC/m2)
real(dp), pointer, dimension(:,:) :: litter_bg             !gridcell below-ground litter (gC/m2)
real(dp), pointer, dimension(:)   :: nind                  !gridcell individual density (indiv/m2)
logical, pointer, dimension(:)    :: tree                  !true if pft is tree
logical, pointer, dimension(:)    :: present               !true if pft is present in that gridcell
real(dp), pointer, dimension(:)   :: maxtwm                !upper limit of temperature of the warmest month
real(dp), pointer, dimension(:)   :: sla                   !PFT specific leaf area (m2/gC)
real(dp), pointer, dimension(:)   :: turnover_ind          !total turnover of living biomass per individual (gC)
real(dp), pointer, dimension(:)   :: heatstress            !reduction in individual density (& establishment) due to heat induced mortality  (indiv/m2)
real(dp), pointer, dimension(:,:) :: bm_inc                !annual biomass increment (gC/m2)
real(dp), pointer, dimension(:)   :: daytempv              !vector of daytime averaged temperatures
real(dp), pointer                 :: twm                   !maximum monthly temperature for the year being simulated

!LOCAL VARIABLES
integer  :: pft,c                        !counters
real(dp) :: bm_delta                     !net individual living biomass increment (incorporating loss through leaf, root and sapwood turnover) (gC)
real(dp) :: mort                         !tree mortality rate
real(dp) :: nind_kill                    !reduction in individual density due to mortality (indiv/m2)
integer  :: d                            !counter
real(dp) :: temp                         !temporary variable
real(dp) :: litter_inc                   !
real(dp) :: greffic                      !
real(dp) :: gddtw                        !

!point pointers
lm_max                          => sv(j)%lm_max
lm_ind                          => sv(j)%lm_ind
sm_ind                          => sv(j)%sm_ind
hm_ind                          => sv(j)%hm_ind
rm_ind                          => sv(j)%rm_ind
litter_ag_fast                  => sv(j)%litter_ag_fast
litter_ag_slow                  => sv(j)%litter_ag_slow
litter_bg                       => sv(j)%litter_bg
nind                            => sv(j)%nind
present                         => sv(j)%present
sla                             => veg%sla
turnover_ind                    => veg%turnover_ind
heatstress                      => veg%heatstress
bm_inc                          => veg%bm_inc
tree                            => prm%tree
maxtwm                          => prm%maxtwm  
daytempv                        => atm%daytempv 
twm                             => atm%twm

!------------------------------------------------------
!begin calculations

 !initialisation
 heatstress = 0._dp

do pft = 1,npft

  if (present(pft) .and. tree(pft)) then

    !Calculate net individual living biomass increment
    bm_delta = max(0._dp,bm_inc(pft,1) / nind(pft) - turnover_ind(pft))

    !Calculate growth efficiency (net biomass increment per unit leaf area)
    if (do_allocation_daily) then
       greffic = bm_delta / (lm_max(pft,1) * sla(pft))
    else
       greffic = bm_delta / (lm_ind(pft,1) * sla(pft))
    end if
    
    !Mortality rate inversely related to growth efficiency (Prentice et al 1993)
    mort = mort_max / (1._dp + k_mort * greffic)

    !heat damage mortality in boreal trees

    if (twm > maxtwm(pft)) then  ! heat damage

      !calculate growing degree days above maxtwm
      gddtw = 0._dp

      do d = 1,365
        gddtw = gddtw + max(0._dp,daytempv(d) - maxtwm(pft))  !NOTE this uses daytempv, which is not retained year over year.
                                                        !if mortality will be called more frequently it should be moved to a
                                                        !location where it is retained.
      end do

      heatstress(pft) = min(1._dp,gddtw / ramp_gddtw)

    else
      heatstress(pft) = 0._dp
    end if

    !Reduce individual density (and thereby gridcell-level biomass) by mortality rate

    mort = min(1._dp,mort + heatstress(pft))
    nind_kill = nind(pft) * mort
    nind(pft) = nind(pft) - nind_kill

    !Transfer lost biomass to litter

    temp = litter_ag_fast(pft,1)
    litter_ag_fast(pft,1) = litter_ag_fast(pft,1) + nind_kill * lm_ind(pft,1)

    if (litter_ag_fast(pft,1) > 1._dp) then  !JM 01.03.2011
    do c = 2,ncvar
      litter_inc = lm_ind(pft,1) * lm_ind(pft,c)
      litter_ag_fast(pft,c) = (litter_ag_fast(pft,c) * temp + litter_inc * nind_kill) / &
                                                litter_ag_fast(pft,1)
    end do
    end if

    temp = litter_ag_slow(pft,1)
    litter_ag_slow(pft,1) = litter_ag_slow(pft,1) + nind_kill * (sm_ind(pft,1) + hm_ind(pft,1))
    
    if (litter_ag_slow(pft,1) > 1._dp) then  !JM 01.03.2011
    do c = 2,ncvar
      litter_inc = (sm_ind(pft,1) * sm_ind(pft,c) + hm_ind(pft,1) * hm_ind(pft,c))
      litter_ag_slow(pft,c) = (litter_ag_slow(pft,c) * temp + litter_inc * nind_kill) / &
                                              litter_ag_slow(pft,1)
    end do
    end if

    temp = litter_bg(pft,1)  !fine roots
    litter_bg(pft,1) = litter_bg(pft,1) + nind_kill * rm_ind(pft,1)

    if (litter_bg(pft,1) > 1._dp) then  !JM 01.03.2011
    do c = 2,ncvar
       litter_bg(pft,c) = (litter_bg(pft,c) * temp + nind_kill * rm_ind(pft,1) * rm_ind(pft,c))&
                           / litter_bg(pft,1)
    end do
    end if

  end if

  if (nind(pft) < 0._dp) present(pft) = .false.

end do

end subroutine mortality

!-------------------------------------------------------------------------------------------------------------

!SUBROUTINE FIRE
!Biomass destruction through disturbance by fire
!recoded from f to f.90 by Joe Melton 2007 from original LPJ code
!methane production and isotopic value added by JM Jan 2008

subroutine fire(j)

use arveparams,only : dp,npft,ncvar,pi,sp
use pftparametersmod, only : prm
use statevars, only : sv,ov
use metvarsmod, only : atm

implicit none

!PARAMETERS
real(dp), parameter :: minfuel = 200.0  !fuel threshold to carry a fire (gC/m2)

!ARGUMENTS
integer, intent(in) :: j

!pointers
real(dp), pointer, dimension(:)   :: dw1                !daily soil layer 1 water content
real(dp), pointer, dimension(:,:) :: lm_ind             !individual leaf mass (gC)
real(dp), pointer, dimension(:,:) :: sm_ind             !individual sapwood mass (gC)
real(dp), pointer, dimension(:,:) :: hm_ind             !individual heartwood mass (gC)
real(dp), pointer, dimension(:,:) :: rm_ind             !individual fine root mass (gC)
real(dp), pointer, dimension(:,:) :: litter_ag_fast     !gridcell above-ground litter (gC/m2)
real(dp), pointer, dimension(:,:) :: litter_ag_slow     !gridcell above-ground litter (gC/m2)
real(dp), pointer, dimension(:)   :: nind               !gridcell individual density (indiv/m2)
real(dp), pointer, dimension(:)   :: fpc_gmax           !gridcell foliar projective cover (FPC) max
logical, pointer, dimension(:)    :: tree               !true if pft is tree
logical, pointer, dimension(:)    :: present            !true if pft is present in that gridcell
real(dp), pointer, dimension(:)   :: flam               !flammability threshold
real(dp), pointer, dimension(:)   :: fireres            !fire resistance index
real(dp), pointer, dimension(:)   :: CH4ef              !Emission factor for CH4 production from
                                                        !fire (g CH4/ kg C per PFT) (Andreae & Merlet GBC 2001)
real(dp), pointer, dimension(:)   :: COef               !Emission factor for CO production from fire (g / kg C per PFT) (Andreae & Merlet GBC 2001)
real(dp), pointer, dimension(:)   :: VOCef              !Emission factor for VOC production from fire (g/ kg C per PFT) (Andreae & Merlet GBC 2001)
real(dp), pointer, dimension(:)   :: TPMef              !Emission factor for total particulate matter production from fire (g/ kgC)(Andreae & Merlet GBC 2001)
real(dp), pointer, dimension(:)   :: NOxef              !Emission factor for NOx production from fire (g/ kg C per PFT) (Andreae & Merlet GBC 2001)
real(dp), pointer, dimension(:)   :: daytempv           !vector of daytime averaged temperatures

!output (sp)
real(sp), pointer                 :: afire_frac         !fraction of gridcell burnt this year
real(sp), pointer                 :: COfire_flux        !CO produced by fire (g/m2)        
real(sp), pointer                 :: VOCfire_flux       !VOC produced by fire (g/m2)       
real(sp), pointer                 :: TPMfire_flux       !TPM produced by fire (g/m2)      
real(sp), pointer                 :: NOxfire_flux       !NOx produced by fire (g/m2)   
real(sp), pointer, dimension(:)   :: acflux_fire        !C flux to atmosphere due to fire (gC/m2)
real(sp), pointer, dimension(:)   :: CH4flux_fire       !CH4 flux to atmosphere due to fire (gCH4/m2)

!LOCAL VARIABLES
integer  :: pft,d,c                                     !
real(dp) :: fire_length                                 !
real(dp), dimension(365) :: fire_prob                   !probability of fire today
real(dp) :: fire_index                                  !
real(dp) :: disturb                                     !
real(dp) :: moistfactor                                 !litter moisture weighting factor
real(dp) :: litter_ag_total                             !total aboveground litter per gridcell (gC m-2)
real(dp) :: fire_term                                   !
real(dp) :: temp                                        !
real(dp) :: all_ind                                     !
real(dp) :: litter_inc                !
real(dp), dimension(npft,ncvar) :: acflux_fire_int      !
real(dp), dimension(npft,ncvar) :: CH4flux_fire_int     !
real(dp), dimension(npft) :: COfire_flux_int            !
real(dp), dimension(npft) :: VOCfire_flux_int           !
real(dp), dimension(npft) :: TPMfire_flux_int           !
real(dp), dimension(npft) :: NOxfire_flux_int           !

!point pointers
lm_ind                  => sv(j)%lm_ind
sm_ind                  => sv(j)%sm_ind
hm_ind                  => sv(j)%hm_ind
rm_ind                  => sv(j)%rm_ind
litter_ag_fast          => sv(j)%litter_ag_fast
litter_ag_slow          => sv(j)%litter_ag_slow
nind                    => sv(j)%nind
present                 => sv(j)%present
dw1                     => surf%Aliq
flam                    => prm%flam
fireres                 => prm%fireres
tree                    => prm%tree
fpc_gmax                => sv(j)%fpc_gmax
CH4ef                   => prm%CH4ef
COef                    => prm%COef
VOCef                   => prm%VOCef
TPMef                   => prm%TPMef
NOxef                   => prm%NOxef
daytempv                => atm%daytempv
acflux_fire             => ov(j)%acflux_fire
CH4flux_fire            => ov(j)%CH4flux_fire
COfire_flux             => ov(j)%COfire_flux
VOCfire_flux            => ov(j)%VOCfire_flux
TPMfire_flux            => ov(j)%TPMfire_flux
NOxfire_flux            => ov(j)%NOxfire_flux
afire_frac              => ov(j)%afire_frac

!------------------------------------
!begin calculations
acflux_fire_int(1:npft,1:ncvar) = 0._dp
acflux_fire(1:ncvar) = 0._dp

!Calculate total above-ground litter per gridcell
litter_ag_total = sum(litter_ag_fast(1:npft,1) + litter_ag_slow(1:npft,1))

!Calculate litter moisture weighting factor
moistfactor  = 0._dp

  if (litter_ag_total > 0._dp) then
    do pft = 1,npft
        moistfactor = moistfactor + ((litter_ag_fast(pft,1) + litter_ag_slow(pft,1)) / litter_ag_total) * flam(pft)
    end do
  end if

!First calculate the annual fraction of the grid cell affected by fire

!Calculate the length of the fire season (units=days)
if (moistfactor > 0._dp) then

  do d = 1,365

   !Calculate today's fire probability, fire_prob. Assume fire is only possible when temperature is above zero
   ! this uses daytempv, which is not retained year over year. This is not a problem if fire is called only once per year
   ! if fire is to be called more freqently it should be moved into a different structure so it is retained.
   if (daytempv(d) > 0._dp) then
     fire_prob(d) = exp((-pi / 4._dp) * (max(0._dp,dw1(d) / moistfactor)**2._dp))
   else
     fire_prob(d) = 0._dp
   end if

  end do

  fire_length = sum(fire_prob(1:365))

  !Calculate annual fire index
  fire_index = fire_length / 365._dp

else

  fire_length = 0._dp
  fire_index = 0._dp

end if

!Calculate the fraction of the grid cell affected by fire
!Assign a minimum fire fraction (for presentational purposes)

 !afire_frac = 1._dp - (exp(-0.2 * fire_index**1.4))**1.1   !NOTE: commented out in orig LPJ code...JM

 fire_term = fire_index - 1._dp
 
 !convert to single precision for netcdf output.
 afire_frac = real(max(fire_index * exp(fire_term / (-0.13 * fire_term**3._dp &
                + 0.6 * fire_term**2._dp + 0.8 * fire_term + 0.45)),0.001))


do pft = 1,npft

  if (present(pft)) then

  !Calculate the available fuel (above-ground litter) to carry the fire

  !Reduce fraction of grid cell affected by fire when fuel
  !becomes limiting (reduced carrying capacity)

                if (litter_ag_total < minfuel * fpc_gmax(pft)) afire_frac = 0.001_dp

        
        !afire_frac = afire_frac * (1._dp / (minfuel**2._dp)) * fuel**2._dp  !NOTE: commented out in orig LPJ code...JM

!Implement the effect of the fire on vegetation structure and litter
!in the disturbed fraction.

        !Each PFT is assigned a resistance to fire, representing the fraction of
        !the PFT which survives a fire. Grasses assumed already to have completed
        !their life cycle and thus are not affected by fire, giving them
        !a competitive advantage against woody PFTs.

  if (tree(pft)) then

    !Calculate the fraction of individuals in grid cell which die

     disturb = (1._dp - fireres(pft)) * afire_frac

     !Calculate carbon flux to atmosphere (gC/m2) due to burnt biomass

     all_ind = lm_ind(pft,1) + sm_ind(pft,1) + hm_ind(pft,1) + rm_ind(pft,1) !FLAG should this be lm_max for daily alloc?

     acflux_fire_int(pft,1) = disturb * (nind(pft) * all_ind)

     if (all_ind > 0._dp) then
          acflux_fire_int(pft,2:ncvar) = (lm_ind(pft,2:ncvar) * lm_ind(pft,1) + sm_ind(pft,2:ncvar) * sm_ind(pft,1)&
                                  + hm_ind(pft,2:ncvar) * hm_ind(pft,1)&
                                  + rm_ind(pft,2:ncvar) * rm_ind(pft,1)) / all_ind
     end if                                  

     !Update the individual density
     nind(pft) = nind(pft) * (1._dp - disturb)

  end if

  !Add combusted litter to carbon flux to atmosphere term

  temp = acflux_fire_int(pft,1)
  acflux_fire_int(pft,1) = acflux_fire_int(pft,1) + (afire_frac * (litter_ag_fast(pft,1) + litter_ag_slow(pft,1)))

  if (acflux_fire_int(pft,1) > 0._dp) then
  
        do c = 2,ncvar
          litter_inc = litter_ag_fast(pft,c) * litter_ag_fast(pft,1) + &
                                litter_ag_slow(pft,c) * litter_ag_slow(pft,1)
          acflux_fire_int(pft,c) = (acflux_fire_int(pft,c) * temp + afire_frac * litter_inc)&
                                 / acflux_fire_int(pft,1)
                                
        end do
        
  else

        acflux_fire_int(pft,2:ncvar) = 0._dp

  end if

  !---
  !Find amount of CH4 and trace gases produced due to fire, and its isotope value.
   CH4flux_fire_int(pft,1) = (acflux_fire_int(pft,1) / 1000._dp) * CH4ef(pft)
   COfire_flux_int(pft)  =  (acflux_fire_int(pft,1) / 1000._dp) * COef(pft)
   VOCfire_flux_int(pft) =  (acflux_fire_int(pft,1) / 1000._dp) * VOCef(pft)
   TPMfire_flux_int(pft) =  (acflux_fire_int(pft,1) / 1000._dp) * TPMef(pft)
   NOxfire_flux_int(pft) =  (acflux_fire_int(pft,1) / 1000._dp) * NOxef(pft)

  !Update the above ground litter term
  litter_ag_fast(pft,1) = (1._dp - afire_frac) * litter_ag_fast(pft,1)
  litter_ag_slow(pft,1) = (1._dp - afire_frac) * litter_ag_slow(pft,1)

  !Calculate carbon flux summed over all pft's
  !also convert to single precision for output
  temp = acflux_fire(1)  
  acflux_fire(1) = real(acflux_fire(1) + acflux_fire_int(pft,1))

  if (acflux_fire(1) > 0._dp) then
       do c = 2,ncvar
         acflux_fire(c) = real((acflux_fire(c) * temp + acflux_fire_int(pft,c) * &
                                acflux_fire_int(pft,1)) / acflux_fire(1))
       end do
  else
       acflux_fire(2:ncvar) = 0.      
  end if

             

 end if !present
end do !pft loop

     !total amount CH4 (g CH4/ m2) and trace gases this year
     ! convert to single precision for output to netcdf (ov)
      CH4flux_fire(1) = real(sum(CH4flux_fire_int(1:npft,1)))
      COfire_flux  = real(sum(COfire_flux_int(1:npft)))
      VOCfire_flux = real(sum(VOCfire_flux_int(1:npft)))
      TPMfire_flux = real(sum(TPMfire_flux_int(1:npft)))
      NOxfire_flux = real(sum(NOxfire_flux_int(1:npft)))

!NOTE This next block of code was causing a crash earlier, keep an eye on it.
!   !find isotopic value
!FLAG JM 01.03.2011, I think it is coded wrong anyway, redo!
    if (CH4flux_fire(1) > 0._dp) then
      do c = 2,ncvar
        CH4flux_fire(c) = real((CH4flux_fire(c) * CH4flux_fire(1) + sum(CH4flux_fire_int(1:npft,c)) * &
                          sum(CH4flux_fire_int(1:npft,1))) / CH4flux_fire(1))
      end do
    end if
!to here.        

!write(*,'(a,4f10.4)')'fire:',afire_frac,fire_length,litter_ag_total,disturb
!write(*,*)acflux_fire(1)

return

end subroutine fire

!-------------------------------------------------------------------------------------------------------------

!SUBROUTINE ESTABLISHMENT
!Establishment of new individuals (saplings) of woody PFTs,
!grass establishment, removal of PFTs not adapted to current climate,
!update of individual structure and FPC.
!recoded from f to f.90 by Joe Melton 2007 from original LPJ code

subroutine establishment(j)

use arveparams,only : dp,npft,ncvar,pi,reinickerp, &
                        allom2,allom3,k5,k6
use pftparametersmod, only : prm
use iovariables, only : dovegetation
use statevars, only : sv,co2
use metvarsmod, only : atm
use iovariables, only : dovegetation

implicit none

!ARGUMENTS
integer, intent(in) :: j        !grid cell index
!integer, intent(in) :: yr        !NOT strictly needed, useful for debugging

!pointers
real(dp), pointer, dimension(:,:) :: lm_ind                     !individual leaf mass (gC)
real(dp), pointer, dimension(:,:) :: sm_ind                     !individual sapwood mass (gC)
real(dp), pointer, dimension(:,:) :: hm_ind                     !individual heartwood mass (gC)
real(dp), pointer, dimension(:,:) :: rm_ind                     !individual fine root mass (gC)
real(dp), pointer, dimension(:,:) :: lm_sapl                    !initial (sapling) leaf mass (gC)
real(dp), pointer, dimension(:,:) :: sm_sapl                    !initial (sapling) sapwood mass (gC)
real(dp), pointer, dimension(:,:) :: hm_sapl                    !initial (sapling) heartwood mass (gC)
real(dp), pointer, dimension(:,:) :: rm_sapl                    !initial (sapling) root mass (gC)
real(dp), pointer, dimension(:,:) :: litter_ag_fast             !gridcell above-ground litter (gC/m2)
real(dp), pointer, dimension(:,:) :: litter_ag_slow             !gridcell above-ground litter (gC/m2)
real(dp), pointer, dimension(:,:) :: litter_bg                  !gridcell below-ground litter (gC/m2)
real(dp), pointer, dimension(:)   :: nind                       !gridcell individual density (indiv/m2)
real(dp), pointer, dimension(:)   :: crownarea                  !crown area (m2)
real(dp), pointer, dimension(:)   :: fpc_gmax                   !maximal gridcell foliar projective cover (FPC)
real(dp), pointer, dimension(:)   :: lai_ind                    !individual leaf area index
real(dp), pointer, dimension(:)   :: lai_max                    !maximumal present individual leaf area index
real(dp), pointer, dimension(:)   :: height                     !tree height (m)
real(dp), pointer, dimension(:,:) :: acflux_estab               !annual biomass increment due to establishment (gC/m2)
logical , pointer, dimension(:)   :: present                    !true if pft is present in that gridcell
real(dp), pointer, dimension(:)   :: sla                        !PFT specific leaf area (m2/gC)
real(dp), pointer, dimension(:)   :: maxcrowna                  !tree maximum crown area (m2)
logical,  pointer, dimension(:)   :: tree                       !true if pft is tree
logical,  pointer, dimension(:)   :: estab                      !whether PFT within bioclimatic limits for establishment
logical,  pointer, dimension(:)   :: survive                    !whether PFT within bioclimatic limits for survival
real(dp), pointer, dimension(:)   :: fpc_ind                    !foliar protective cover per pft individual (m2 m-2)
real(dp), pointer, dimension(:)   :: latosa                     !leaf area to sapwood area (m2 m-2)
real(dp), pointer, dimension(:)   :: allom1                     !Packing constrain from Reinicke's rule exponent
real(dp), pointer, dimension(:,:) :: bm_inc                     !annual biomass increment (gC/m2)    
logical, pointer, dimension(:)    :: c4                         !true if pft is c4
real(dp), pointer, dimension(:,:) :: lm_max                     !
integer, pointer, dimension(:)    :: phentype                   !
integer, pointer, dimension(:)    :: pstate                     !
real(dp), pointer, dimension(:)   :: crownarea_max              !tree maximum crown area (m2) 
integer, pointer, dimension(:)    :: leafonday                  !day of year of leaf onset
real(dp), pointer, dimension(:,:) :: NSC_storage                !individual plant non-structural carbon stores (gC)
real(dp), pointer, dimension(:,:) :: NSC_sapl                   !initial (sapling) plant non-strucutral carbon store (gC)
real(dp), pointer, dimension(:)   :: NSC_limit                  !allometric plant non-structural carbon stores size (gC/indiv)
real(dp), pointer, dimension(:)   :: NSC_opt                    !allometric plant non-structural carbon multiple of maximal leaf mass 
real(dp), pointer, dimension(:)   :: wooddens                   !density of wood (sap + heart) g m-3
real(dp), pointer                 :: aprec                      !annual precipitation (mm)

!LOCAL VARIABLES
integer  :: pft,c                       !counters
integer  :: npft_estab                  !number of regenerating tree PFTs
real(dp) :: fpc_tree_total              !total grid FPC for tree PFTs
real(dp) :: estab_rate                  !sapling establishment rate over area available for establishment (indiv/m2)
real(dp) :: estab_grid                  !grid-level establishment rate (indiv/m2)
real(dp) :: nind_old                    !number of individuals /m2 before establishment
real(dp) :: stem_diam                    !stem diameter (m)
real(dp) :: fpc_total                   !total grid FPC
real(dp) :: bare                        !gridcell bare ground fraction
integer  :: ngrass
real(dp) :: fpc_grass_total
real(dp) :: temp,litter_inc
real(dp) :: all
real(dp) :: sm_ind_temp

!PARAMETERS
real(dp), parameter :: aprec_min_estab = 100.0          !minimum annual precipitation for establishment (mm)
real(dp), parameter :: estab_max = 0.24                 !maximum sapling establishment rate (indiv/m2) (old -estab_max=0.12)
real(dp), parameter :: nind_min = 1.0e-10               !minimum individual density for persistence of PFT (indiv/m2)
real(dp), parameter :: eps = 1.e-6                      !

!point pointers
lm_ind                => sv(j)%lm_ind
sm_ind                => sv(j)%sm_ind
hm_ind                => sv(j)%hm_ind
rm_ind                => sv(j)%rm_ind
lm_sapl               => veg%lm_sapl
sm_sapl               => veg%sm_sapl
hm_sapl               => veg%hm_sapl
rm_sapl               => veg%rm_sapl
litter_ag_fast        => sv(j)%litter_ag_fast      
litter_ag_slow        => sv(j)%litter_ag_slow      
litter_bg             => sv(j)%litter_bg
nind                  => sv(j)%nind
crownarea             => sv(j)%crownarea
fpc_gmax              => sv(j)%fpc_gmax
lai_ind               => sv(j)%lai_ind
lai_max               => sv(j)%lai_max
present               => sv(j)%present
height                => sv(j)%height
acflux_estab          => veg%acflux_estab
sla                   => veg%sla
fpc_ind               => veg%fpc_ind
maxcrowna             => prm%maxcrowna
tree                  => prm%tree
estab                 => veg%pft_estab
survive               => veg%pft_survive
pstate                => sv(j)%pstate
latosa                => prm%latosa
allom1                => prm%allom1
bm_inc                => veg%bm_inc
c4                    => prm%c4
lm_max                => sv(j)%lm_max
phentype              => prm%phentype
crownarea_max         => prm%maxcrowna
leafonday             => sv(j)%leafonday
NSC_storage           => sv(j)%NSC_storage 
NSC_sapl              => veg%NSC_sapl  
NSC_limit             => veg%NSC_limit  
NSC_opt               => prm%NSC_opt
wooddens              => prm%wooddens
aprec                 => atm%aprec

!-------------------------------------------------------------------------------
!begin calculations
!begin calculations

do pft = 1,npft

if (do_allocation_daily) then
 if (present(pft)) then
  if (tree(pft)) then
  
  ! Reproduction.
  ! isotope values for new establishment saplings (NOTE! NSC added in without testing or checking for different
  ! fractionation values! JM 07.01.2011)

             if (bm_inc(pft,1) > 0._dp) then

                lm_sapl(pft,2:ncvar) = bm_inc(pft,2:ncvar)
                sm_sapl(pft,2:ncvar) = bm_inc(pft,2:ncvar)
                hm_sapl(pft,2:ncvar) = bm_inc(pft,2:ncvar)
                rm_sapl(pft,2:ncvar) = bm_inc(pft,2:ncvar)
                NSC_sapl(pft,2:ncvar) = bm_inc(pft,2:ncvar)
                
             else

                lm_sapl(pft,2) = 17.8 - co2(2) !C3 value
                sm_sapl(pft,2) = 17.8 - co2(2) !from lloyd & farquhar,1994
                hm_sapl(pft,2) = 17.8 - co2(2)
                rm_sapl(pft,2) = 17.8 - co2(2)
                NSC_sapl(pft,2) = 17.8 - co2(2)
                lm_sapl(pft,3) = co2(3)
                sm_sapl(pft,3) = co2(3)
                hm_sapl(pft,3) = co2(3)
                rm_sapl(pft,3) = co2(3)
                NSC_sapl(pft,3) = co2(3)

             end if
    else  !grass
              if (bm_inc(pft,1) > 0._dp) then

                lm_sapl(pft,2:ncvar) = bm_inc(pft,2:ncvar)
                rm_sapl(pft,2:ncvar) = bm_inc(pft,2:ncvar)
                NSC_sapl(pft,2:ncvar) = bm_inc(pft,2:ncvar)

              else
                 if (c4(pft)) then

                  lm_sapl(pft,2) = 3.6 - co2(2)  !C4 value
                  rm_sapl(pft,2) = 3.6 - co2(2)  !from lloyd & farquhar,1994
                  NSC_sapl(pft,2) = 3.6 - co2(2)
                  lm_sapl(pft,3) = co2(3)
                  rm_sapl(pft,3) = co2(3)
                  NSC_sapl(pft,3) = co2(3)
                  
                else

                  lm_sapl(pft,2) = 17.8 - co2(2) !C3 value
                  rm_sapl(pft,2) = 17.8 - co2(2) !from lloyd & farquhar,1994
                  NSC_sapl(pft,2) = 17.8 - co2(2)
                  lm_sapl(pft,3) = co2(3)
                  rm_sapl(pft,3) = co2(3)
                  NSC_sapl(pft,3) = co2(3)

                end if
              end if
    end if
  end if
end if

!---------------------- 

!Kill PFTs not adapted to current climate, introduce newly "adapted" PFTs

  if (present(pft) .and. (.not. survive(pft) .or. nind(pft) < nind_min)) then !kill PFT

    present(pft) = .false.

    !Add killed biomass to litter

    if (tree(pft)) then

      temp = litter_ag_fast(pft,1)
      litter_ag_fast(pft,1) = litter_ag_fast(pft,1) + nind(pft) * lm_ind(pft,1)

      if (litter_ag_fast(pft,1) > 1._dp) then !JM 01.03.2011
      do c = 2,ncvar
        litter_inc = lm_ind(pft,1) * lm_ind(pft,c)
        litter_ag_fast(pft,c) = (litter_ag_fast(pft,c) * temp + litter_inc * nind(pft)) / litter_ag_fast(pft,1)
      end do
      end if
      
      temp = litter_ag_slow(pft,1)
      litter_ag_slow(pft,1) = litter_ag_slow(pft,1) + nind(pft) * (sm_ind(pft,1) + hm_ind(pft,1))

      if (litter_ag_slow(pft,1) > 1._dp) then !JM 01.03.2011
      do c = 2,ncvar
        litter_inc = (sm_ind(pft,1) * sm_ind(pft,c) + hm_ind(pft,1) * hm_ind(pft,c))
        litter_ag_slow(pft,c) = (litter_ag_slow(pft,c) * temp + litter_inc * nind(pft)) / litter_ag_slow(pft,1)
      end do
      end if

    else                  !grasses

      temp = litter_ag_fast(pft,1)
      litter_ag_fast(pft,1) = litter_ag_fast(pft,1) + lm_ind(pft,1) * nind(pft)

      if (litter_ag_fast(pft,1) > 1._dp) then !JM 01.03.2011
      do c = 2,ncvar
        litter_inc = lm_ind(pft,1) * lm_ind(pft,c)
        litter_ag_fast(pft,c) = (litter_ag_fast(pft,c) * temp + litter_inc * nind(pft)) / litter_ag_fast(pft,1)
      end do
      end if

    end if

    temp = litter_bg(pft,1)
    litter_bg(pft,1) = litter_bg(pft,1) + rm_ind(pft,1) * nind(pft)

    if (litter_bg(pft,1) > 1._dp) then !JM 01.03.2011
     litter_bg(pft,2:ncvar) = (litter_bg(pft,2:ncvar) * temp+rm_ind(pft,2:ncvar) * rm_ind(pft,1) * nind(pft)) / litter_bg(pft,1)
    end if

    fpc_gmax(pft) = 0._dp

  else if (.not. present(pft) .and. survive(pft) .and. estab(pft) .and. aprec >= aprec_min_estab) then

    !Introduce PFT if conditions suitable for establishment

    if (dovegetation) then
    
    	present(pft) = .true.
        
        !for annual C allocation, introduce pft with no leaves on.
        pstate(pft) = 0
        
        if (do_allocation_daily) then
          
          if (phentype(pft) /= 1) then
           pstate(pft) = 3         !set pstate to 3 so that a lai_vir is assigned.
          else
           pstate(pft) = 1         !evergreens are always in pstate 1 (full leafout)
          end if
        
          leafonday(pft) = 1      !set to 1 so that needed_sap is defined.
    
        end if 
    
    else  !in joboptions file it is set to false so the simulation is bare ground, no vegetation!
    
    	present(pft) = .false. 
	   
    end if
    
    
    ! initial values prior to the sapling establishment  
    if (tree(pft)) then
      nind(pft) = 0._dp
    else
      nind(pft) = 1._dp    !each grass PFT = 1 "individual"
      crownarea(pft) = 1._dp !assumed as a 'green carpet'
    end if

    lm_ind(pft,1:ncvar) = 0._dp
    sm_ind(pft,1:ncvar) = 0._dp
    rm_ind(pft,1:ncvar) = 0._dp
    hm_ind(pft,1:ncvar) = 0._dp
    fpc_gmax(pft) = 0._dp

  end if

end do

!----------
!SAPLING AND GRASS ESTABLISHMENT

!Calculate total woody FPC and number of woody PFTs present and able to establish (do NOT replace below)
 
    fpc_total       = sum(fpc_gmax, mask = present)
    fpc_tree_total  = sum(fpc_gmax, mask = present .and. tree)
    ngrass          = count(present .and. .not. tree)
    npft_estab      = count(present .and. tree .and. estab)
    fpc_grass_total = sum(fpc_gmax,  mask = present .and. .not. tree)

!Prohibit establishment under extreme temperature or water stress.

do pft = 1,npft

  if (aprec >= aprec_min_estab .and. npft_estab > 0) then

    !Calculate establishment rate over available space, per tree PFT
    !Maximum establishment rate reduced by shading as tree FPC approaches 1
    !Total establishment rate partitioned equally among regenerating woody PFTs

    estab_rate = estab_max * (1._dp - exp(5._dp * (fpc_tree_total - 1._dp))) / real(npft_estab)

    !Calculate grid-level establishment rate per woody PFT
    !Space available for woody PFT establishment is proportion of grid cell
    !not currently occupied by woody PFTs.

    estab_grid = estab_rate * (1._dp - fpc_tree_total)

  else                   !unsuitable climate for establishment

    estab_grid = 0._dp

  end if
  
  
!if (yr > 1) then !FLAG!!
! if (pft == 1) write(*,*)'ESTABLISHMENT FOR FURTHER YEARS IS TURNED OFF'
!if (sv(j)%yrpres(pft) > 1) then
!estab_rate = 0._dp
!estab_grid = 0._dp
!estab(pft) = .false.
!end if
!end if !FLAG!

  if (present(pft) .and. tree(pft) .and. estab(pft)) then

    !Add new saplings to current population

    nind_old = nind(pft)
    nind(pft) = nind_old + estab_grid

    if (nind(pft) > 0._dp) then
    
      if (do_allocation_daily) then
    
        !Add leafmass to lm_ind if it is an evergreen
         if (phentype(pft) == 1) then  

                temp = lm_ind(pft,1)
                lm_ind(pft,1) = (lm_ind(pft,1) * nind_old + lm_sapl(pft,1) * estab_grid) / nind(pft)

                if (lm_ind(pft,1) > 0._dp) then
                  lm_ind(pft,2:ncvar) = (lm_ind(pft,2:ncvar) * temp * nind_old + lm_sapl(pft,2:ncvar) &
                  * lm_sapl(pft,1) * estab_grid) / (lm_ind(pft,1) * nind(pft))
                end if

                lm_max(pft,1) = max(lm_max(pft,1), lm_ind(pft,1))
                
                temp = rm_ind(pft,1)
                rm_ind(pft,1) = (rm_ind(pft,1) * nind_old + rm_sapl(pft,1) * estab_grid) / nind(pft)

                if (rm_ind(pft,1) > 0.d0) then
                  rm_ind(pft,2:ncvar) = (rm_ind(pft,2:ncvar) * temp * nind_old + rm_sapl(pft,2:ncvar) * rm_sapl(pft,1) * estab_grid) / &
                                        (rm_ind(pft,1) * nind(pft))
                end if
                
                !rm_max(pft,1) = max(rm_max(pft,1), rm_ind(pft,1))
                
          else !deciduous
          
                !The carbon to be added is saved for the next leaf out.
          !FLAG this is not completed! JM 28.03.2011
                temp = lm_ind(pft,1)
                lm_ind(pft,1) = (lm_ind(pft,1) * nind_old + lm_sapl(pft,1) * estab_grid) / nind(pft)

                if (lm_ind(pft,1) > 0._dp) then
                  lm_ind(pft,2:ncvar) = (lm_ind(pft,2:ncvar) * temp * nind_old + lm_sapl(pft,2:ncvar) &
                  * lm_sapl(pft,1) * estab_grid) / (lm_ind(pft,1) * nind(pft))
                end if

                lm_max(pft,1) = max(lm_max(pft,1), lm_ind(pft,1))
                
                temp = rm_ind(pft,1)
                rm_ind(pft,1) = (rm_ind(pft,1) * nind_old + rm_sapl(pft,1) * estab_grid) / nind(pft)

                if (rm_ind(pft,1) > 0.d0) then
                  rm_ind(pft,2:ncvar) = (rm_ind(pft,2:ncvar) * temp * nind_old + rm_sapl(pft,2:ncvar) * rm_sapl(pft,1) * estab_grid) / &
                                        (rm_ind(pft,1) * nind(pft))
                end if
                
               ! rm_max(pft,1) = max(rm_max(pft,1), rm_ind(pft,1))                

          end if
                  
            temp = sm_ind(pft,1)
            sm_ind(pft,1) = (sm_ind(pft,1) * nind_old + sm_sapl(pft,1) * estab_grid) / nind(pft)

            if (sm_ind(pft,1) > 0._dp) then
              sm_ind(pft,2:ncvar) = (sm_ind(pft,2:ncvar) * temp * nind_old + sm_sapl(pft,2:ncvar) * sm_sapl(pft,1) * estab_grid) / &
                                    (sm_ind(pft,1) * nind(pft))
            end if

            temp = NSC_storage(pft,1)
            NSC_storage(pft,1) = (NSC_storage(pft,1) * nind_old + NSC_sapl(pft,1) * estab_grid) / nind(pft)

            if (NSC_storage(pft,1) > 0._dp) then
              NSC_storage(pft,2:ncvar) = (NSC_storage(pft,2:ncvar) * temp * nind_old + NSC_sapl(pft,2:ncvar) * NSC_sapl(pft,1) &
              * estab_grid) / (NSC_storage(pft,1) * nind(pft))
            end if
           
           ! NOTE: for daily C allocation, there is assumed to be no heartwood in the seedlings
           ! so this is an adjustment only for the new nind.
           temp = hm_ind(pft,1)      
           hm_ind(pft,1) = (hm_ind(pft,1) * nind_old + hm_sapl(pft,1) * estab_grid) / nind(pft)

           if (hm_ind(pft,1) > 0._dp) then
             hm_ind(pft,2:ncvar) = (hm_ind(pft,2:ncvar) * temp * nind_old + hm_sapl(pft,2:ncvar) * hm_sapl(pft,1) * estab_grid) / &
                                   (hm_ind(pft,1) * nind(pft))
           end if

            
        else  !annual C allocation scheme        

          temp = lm_ind(pft,1)
          lm_ind(pft,1) = (lm_ind(pft,1) * nind_old + lm_sapl(pft,1) * estab_grid) / nind(pft)

          if (lm_ind(pft,1) > 0._dp) then
            lm_ind(pft,2:ncvar) = (lm_ind(pft,2:ncvar) * temp * nind_old + lm_sapl(pft,2:ncvar) * lm_sapl(pft,1) * estab_grid) / &
                                  (lm_ind(pft,1) * nind(pft))
          end if
          
          temp = sm_ind(pft,1)
          sm_ind_temp = (sm_ind(pft,1) * nind_old + sm_sapl(pft,1) * estab_grid) / nind(pft)

          if (sm_ind_temp > 0._dp) then
            sm_ind(pft,2:ncvar) = (sm_ind(pft,2:ncvar) * temp * nind_old + sm_sapl(pft,2:ncvar) * sm_sapl(pft,1) * estab_grid) / &
                                  (sm_ind_temp * nind(pft))
          end if

         temp = rm_ind(pft,1)
         rm_ind(pft,1) = (rm_ind(pft,1) * nind_old + rm_sapl(pft,1) * estab_grid) / nind(pft)

         if (rm_ind(pft,1) > 0.d0) then
           rm_ind(pft,2:ncvar) = (rm_ind(pft,2:ncvar) * temp * nind_old + rm_sapl(pft,2:ncvar) * rm_sapl(pft,1) * estab_grid) / &
                                 (rm_ind(pft,1) * nind(pft))
         end if

         temp = hm_ind(pft,1)      
         hm_ind(pft,1) = (hm_ind(pft,1) * nind_old + hm_sapl(pft,1) * estab_grid) / nind(pft)

         if (hm_ind(pft,1) > 0._dp) then
           hm_ind(pft,2:ncvar) = (hm_ind(pft,2:ncvar) * temp * nind_old + hm_sapl(pft,2:ncvar) * hm_sapl(pft,1) * estab_grid) / &
                                 (hm_ind(pft,1) * nind(pft))
         end if

          !annual does not track NSC
           NSC_storage(:,:) = 0._dp  

        end if !daily/annual
      

    else  !no individuals
      lm_ind(pft,:) = 0._dp
      sm_ind(pft,:) = 0._dp
      hm_ind(pft,:) = 0._dp
      rm_ind(pft,:) = 0._dp
      NSC_storage(pft,:) = 0._dp
    end if

    !Accumulate biomass increment due to sapling establishment (NOTE: Added NSC_sapl to this. JM 07.01.2011)

    temp = acflux_estab(pft,1)
    all = lm_sapl(pft,1) + sm_sapl(pft,1) + hm_sapl(pft,1) + rm_sapl(pft,1) + NSC_sapl(pft,1)

    if (all * estab_grid > eps) then
      acflux_estab(pft,1) = acflux_estab(pft,1) + all * estab_grid

      if (acflux_estab(pft,1) > 0._dp .and. temp >= 0._dp) then
        acflux_estab(pft,2:ncvar) = (acflux_estab(pft,2:ncvar) * temp + (lm_sapl(pft,1) * lm_sapl(pft,2:ncvar) + sm_sapl(pft,1)&
                        * sm_sapl(pft,2:ncvar) + hm_sapl(pft,1) * hm_sapl(pft,2:ncvar) + rm_sapl(pft,1)&
                        * rm_sapl(pft,2:ncvar) + NSC_sapl(pft,1) * NSC_sapl(pft,2:ncvar)) * estab_grid) / acflux_estab(pft,1)

      else if (acflux_estab(pft,1) > 0._dp .and. temp < 0._dp) then
        acflux_estab(pft,2:ncvar) = (lm_sapl(pft,1) * lm_sapl(pft,2:ncvar) + sm_sapl(pft,1) * sm_sapl(pft,2:ncvar) + hm_sapl(pft,1)&
                              * hm_sapl(pft,c) + rm_sapl(pft,1) * rm_sapl(pft,2:ncvar) + NSC_sapl(pft,1) * NSC_sapl(pft,2:ncvar)) / all
      end if
    end if


    if (do_allocation_daily) then
    
             ! ---NEW WAY
         ! Calculate height, diameter and crown area for new average
         ! individual such that the basic allometric relationships (A-E below)
         ! are satisfied. NOTE: We do not include NSC pools in the plant structural
         ! calculations
    
  	 ! Eqn A below is being related by a relationship from Hickler's hydraulic routine.
	 ! (A) (leaf area) = latosa * (sapwood xs area)
	 !    (Pipe Model, Shinozaki et al. 1964a,b; Waring et al 1982)
    
	 ! (B) (leaf mass) = (betallom(pft) * root mass)**alphallom(pft)
         !      (relations from Enquist & Niklas Science 2002, Niklas 2005 Ann of Botany)
         
	 ! (C) height = k5 * (stem diameter)**2.3 - k6
	 !    (Niklas & Spatz PNAS 2004)
  
	 ! (D) (crown area) = min (allom1 * (stem diameter)**allom5,
	 !			      crownarea_max)
         !      (Inversion of fig 2d, Enquist & Niklas Nature 2001)
        
         ! definition:
             !(0) leaf_area = (lm_ind) * sla    
         ! From (A) and (0)
           !(1) stem_xsa = pi * stem_diam**2 / 4
           !(2) wooddens = (sm_ind + hm_ind) / stemvolume
           !(3) stemvolume = stem_xsa * height = (sm_ind + hm_ind) / wooddens
         ! From (1), (2) & (3),
           !(4) stem_xsa = (sm_ind + hm_ind) / wooddens / height

         ! From (5) & (4) & (1),
         ! (6) height = (sm_ind / wooddens) / ((lm_ind * sla) / latosa)
         ! stem_diam = FLAG finish this derivation for the notes.

! JM 21.02.2011
! Try new stem diam calculation from Niklas and Spatz 2004          
!         stem_diam = (4.d0 * (hm_ind(pft,1) + sm_ind(pft,1)) * lm_max(pft,1) * sla(pft) / (pi * sm_ind(pft,1) * latosa(pft)))**0.5 

          stem_diam = ((sm_ind(pft,1) + hm_ind(pft,1)) * 1.e-3 /(202.3 * k5))**(0.375)
!--

         height(pft) = max(0.05,k5 * stem_diam**(2. / 3.) - k6) 
         crownarea(pft) = min(allom1(pft) * stem_diam**reinickerp,crownarea_max(pft))  !LPJ relation
         
    else !LPJ relations!
    
        !Calculate height, diameter and crown area for new average
        !individual such that the basic allometric relationships (A-C below)
        !are satisfied.

        !(A) (leaf area) = latosa * (sapwood xs area)
            !(Pipe Model, Shinozaki et al. 1964a,b; Waring et al 1982)
        !(B) (leaf mass) = lmtorm * (root mass)
        !(C) height = allom2 * (stem diameter)**allom3    !   (source?)
        !(D) (crown area) = min (allom1 * (stem diameter)**reinickerp, maxcrowna(pft))

        !From (A),
          !(1) sap_xsa = lm_ind * sla / latosa
          !(2) wooddens = (sm_ind + hm_ind) / stemvolume
          !(3) stemvolume = stem_xsa * height
        !From (1), (2) & (3),
          !(4) stem_xsa = (sm_ind + hm_ind) / wooddens / height
          !(5) stem_xsa = pi * (stem_diam**2) / 4
        !From (5),
          !(6) stem_diam = ( 4 * stem_xsa / pi )**0.5
        !From (4) & (6),
          !(7) stem_diam = ( 4 * (sm_ind + hm_ind) / wooddens / height / pi )**0.5
        !From (C) & (7),
          !(8) stem_diam = ( 4 * (sm_ind + hm_ind) / wooddens / ( allom2 * stem_diam**allom3 ) / pi )**0.5
        !From (8),
          !(9) stem_diam = ( 4 * (sm_ind + hm_ind ) / wooddens / pi / allom2 )**( 1 / (2 + allom3) )

        stem_diam = (4._dp * (sm_ind_temp + hm_ind(pft,1)) / wooddens(pft) / pi / &
                    allom2)**(1._dp / (2._dp + allom3)) !Eqn 9
        height(pft) = allom2 * stem_diam**allom3 !Eqn C
        crownarea(pft) = min(maxcrowna(pft),allom1(pft) * stem_diam**reinickerp) !Eqn D

        !Recalculate sapwood mass, transferring excess sapwood to heartwood compartment, if necessary to satisfy Eqn A

        sm_ind(pft,1) = lm_ind(pft,1) * height(pft) * wooddens(pft) * sla(pft) / latosa(pft)

        temp = hm_ind(pft,1)
        hm_ind(pft,1) = hm_ind(pft,1) + (sm_ind_temp - sm_ind(pft,1))

        if (hm_ind(pft,1) > 0._dp) then
          hm_ind(pft,2:ncvar) = (hm_ind(pft,2:ncvar) * temp + (sm_ind_temp - sm_ind(pft,1)) * sm_ind(pft,2:ncvar)) / hm_ind(pft,1)
        end if

    end if !end daily/annual allocation section

  else if (present(pft) .and. .not. tree(pft)) then

    if (estab(pft)) then

      !Grasses can establish in non-vegetated areas
      if (ngrass > 0) then
        bare = (1._dp - fpc_total) / real(ngrass)
      else
        bare = 0._dp
      end if
      
      if (do_allocation_daily) then  !daily alloc routine
      
          !Add increment to lm_max FLAG!....Is this okay? What about when their are grasses already there?
          !is there a better way for this? JM 23.02.2011
           temp = lm_max(pft,1)
           lm_max(pft,1) = lm_max(pft,1) + bare * lm_sapl(pft,1)

           if (lm_max(pft,1) /= 0._dp) then
             lm_max(pft,2:ncvar) = (lm_max(pft,2:ncvar) * temp + bare * lm_sapl(pft,1) * lm_sapl(pft,2:ncvar)) / lm_max(pft,1)
           end if

           if (phentype(pft) == 1) then  !evergreens are given the lm immediately
           lm_ind(pft,1) = lm_max(pft,1)
           end if
           
            temp = NSC_storage(pft,1)
           NSC_storage(pft,1) = NSC_storage(pft,1) + bare * NSC_sapl(pft,1)

           if (NSC_storage(pft,1) /= 0._dp) then
             NSC_storage(pft,2:ncvar) = (NSC_storage(pft,2:ncvar) * temp + bare * NSC_sapl(pft,1)&
                                         * NSC_sapl(pft,2:ncvar)) / NSC_storage(pft,1)
           end if

           
      else  !annual alloc routine

          temp = lm_ind(pft,1)
          lm_ind(pft,1) = lm_ind(pft,1) + bare * lm_sapl(pft,1)

          if (lm_ind(pft,1) /= 0._dp) then
            lm_ind(pft,2:ncvar) = (lm_ind(pft,2:ncvar) * temp + bare * lm_sapl(pft,1) * lm_sapl(pft,2:ncvar)) / lm_ind(pft,1)
          end if
          
      end if

      temp = rm_ind(pft,1)
      rm_ind(pft,1) = rm_ind(pft,1) + bare * rm_sapl(pft,1)

      if (rm_ind(pft,1) /= 0._dp) then
        rm_ind(pft,2:ncvar) = (rm_ind(pft,2:ncvar) * temp + bare * rm_sapl(pft,1) * rm_sapl(pft,2:ncvar)) / rm_ind(pft,1)
      end if

      !Accumulate biomass increment due to grass establishment

      temp = acflux_estab(pft,1)
      all = bare * (lm_sapl(pft,1) + rm_sapl(pft,1) + NSC_sapl(pft,1))  !NSC_sapl is 0 for annual allocation

      if (all * crownarea(pft) > eps) then
        acflux_estab(pft,1) = acflux_estab(pft,1) + all * crownarea(pft)

     
             if (acflux_estab(pft,1) > 0._dp .and. temp >= 0._dp) then
               acflux_estab(pft,2:ncvar) = (acflux_estab(pft,2:ncvar) * temp + bare * (lm_sapl(pft,1) * &
                                             lm_sapl(pft,2:ncvar) + rm_sapl(pft,1) * rm_sapl(pft,2:ncvar)&
                                              + NSC_sapl(pft,1) * NSC_sapl(pft,2:ncvar)) * crownarea(pft)) / &
                                             acflux_estab(pft,1)
             else if (acflux_estab(pft,1) > 0._dp .and. temp < 0._dp) then
               acflux_estab(pft,2:ncvar) = bare * (lm_sapl(pft,1) * lm_sapl(pft,2:ncvar) + rm_sapl(pft,1) * &
               rm_sapl(pft,2:ncvar) + NSC_sapl(pft,1) * NSC_sapl(pft,2:ncvar)) / all
             end if
             
      end if  

    end if  !estab
    
    if (do_allocation_daily) then
    
        ! CHANGED to rm, from lm.
        if (rm_ind(pft,1) <= 0._dp) then
        
          ! no rm_ind, so kill plant      
          present(pft) = .false.

        else !plant has rm_ind

          ! use CTEM (Arora and Boer GBC 2005) grass height (m) (expects leaf mass in kg)
          ! min height is 0.05m for resistance calcs.
          height(pft) = max(0.05,3.5_dp * (lm_ind(pft,1) * 0.001_dp)**0.5_dp) 

        end if !rm_ind
        
    else !annual C allocation 

       ! if no lm_ind then kill the plant off 
       if (lm_ind(pft,1) <= 0._dp) present(pft) = .false.

    end if  !annual/daily alloc.   
  end if  !grass/tree

  if (present(pft)) then
  
  if (do_allocation_daily) then

   !Update LAI and FPC for changes due to establishment.

    if (crownarea(pft) > 0._dp) then
    
      !find the optimal size of the non-structural carbon pool for this maximal leaf mass
      NSC_limit(pft) = lm_max(pft,1) * NSC_opt(pft)

      lai_ind(pft) = lm_ind(pft,1) * sla(pft) / crownarea(pft)
      lai_max(pft) = lm_max(pft,1) * sla(pft) / crownarea(pft)

      fpc_ind(pft) = 1._dp - exp(-0.5 * lai_ind(pft))  !FLAG right?
      
      fpc_gmax(pft) = crownarea(pft) * nind(pft) * (1._dp - exp(-0.5 * lai_max(pft)))

    else
    
      lai_ind(pft) = 0._dp
      lai_max(pft) = 0._dp
      fpc_ind(pft) = 0._dp

      fpc_gmax(pft) = 0._dp

    end if   !crownarea
  
  else  !annual

    !Update LAI and FPC

    if (crownarea(pft) > 0._dp) then
      lai_max(pft) = lm_ind(pft,1) * sla(pft) / crownarea(pft)
      fpc_gmax(pft) = crownarea(pft) * nind(pft) * (1._dp - exp(-0.5 * lai_max(pft)))

    else
      lai_max(pft) = 0._dp  
      fpc_gmax(pft) = 0._dp
    end if
   end if !annual/daily
 end if  !present
end do   !pft loop

end subroutine establishment

!-------------------------------------------------------------------------------------------------------------

subroutine roots(j)
!Rooting depth and lateral spread is based around Arora & Boer, Earth Interactions, Vol 7 paper no. 6. 2003
!coded and adapted by Joe Melton Nov 24 2007

use arveparams,          only : dp,npft,nl,OMorgC
use statevars,          only : sv
use pftparametersmod, only : prm
use soilstate_vars, only : surf

implicit none

!------------------------

!arguments
integer, intent(in) :: j                        !gridcell index

!local variables
integer :: l,pft
real(dp) :: Broot                               !root biomass in kg organic/m2
real(dp) :: rootgro                             !parameter for root growth direction (0 -1 )
real(dp) :: Bunit
real(dp) :: a                                   !vegetation dependent inverse e-folding length scale
real(dp) :: b
real(dp), dimension(nl)      :: rootf           !individual cumulative distribution profiles
!real(dp), dimension(nl,npft) :: rootdens        !root density per pft and per layer (g/m2)
real(dp) :: trm_ind                             !individual total root mass (gC)
real(dp) :: term_a                              !temporary variable
real(dp) :: term_b                              !temporary variable

!input pointers
integer, pointer  :: gnl                        !index value of lowest 'soil' layer
real(dp), pointer, dimension(:)  :: abar        !parameter representing mean root distribution profile
real(dp), pointer, dimension(:)  :: Bbar        !average standing root biomass (kg m-2)
real(dp), pointer, dimension(:)  :: littleb     !parameter representing variable root distribution profile (=abar * Bbar exp alpha)
real(dp), pointer :: soildepth                  !depth to bedock (m)
real(dp), pointer, dimension(:,:) :: rm_ind     !individual fine root mass (gC)
real(dp), pointer, dimension(:,:) :: rm_max     !maximal individual fine root mass (gC)
real(dp), pointer, dimension(:,:) :: lm_ind     !individual leaf mass (gC)
real(dp), pointer, dimension(:)   :: zipos      !z coordinate of the interface between soil layers relative to soil surface (m), positive downwards
real(dp), pointer, dimension(:,:) :: rootfracl   !root fraction in that soil layer
real(dp), pointer, dimension(:)   :: rootdepth  !root depth (m)
real(dp), pointer, dimension(:)   :: nind       !gridcell individual density (indiv/m2)
logical, pointer, dimension(:)    :: present    !if PFT is present in grid cell
real(dp), pointer, dimension(:)   :: crownarea  !m2 crown per individual
logical, pointer, dimension(:)    :: tree       !true if tree
real(dp), pointer, dimension(:)   :: frac       !fraction of sapwood and heartwood that is below ground (in coarse roots)(gC)
real(dp), pointer, dimension(:,:) :: sm_ind     !individual sapwood mass (gC)
real(dp), pointer, dimension(:,:) :: hm_ind     !individual heartwood mass (gC)

!point pointers
gnl               => sv(j)%gnl
rm_ind            => sv(j)%rm_ind
rm_max            => sv(j)%rm_max
zipos             => sv(j)%zipos
rootfracl         => sv(j)%rootfracl
soildepth         => sv(j)%soildepth
rootdepth         => sv(j)%rootdepth
present           => sv(j)%present
nind              => sv(j)%nind
abar              => prm%abar
Bbar              => prm%Bbar
littleb           => prm%littleb
lm_ind            => sv(j)%lm_ind
tree              => prm%tree
crownarea         => sv(j)%crownarea
frac              => veg%frac
sm_ind            => sv(j)%sm_ind
hm_ind            => sv(j)%hm_ind

!---------------------------
!calculations begin

do pft = 1,npft

  if (present(pft)) then

    !this subroutine is expecting 'total' root mass (i.e. both coarse (which we do not explicitly
    !simulate at present) and fine). So simulate coarse roots as a multiple of fine. Basically, some
    !amount of support is required for all the fine roots to be attached to.
    
      
      if (do_allocation_daily) then 
        if (tree(pft)) then
         !frac is used here to approximate the coarse root amounts
          trm_ind = rm_max(pft,1) + frac(pft) * (sm_ind(pft,1) + hm_ind(pft,1)) 
        else
          trm_ind = rm_max(pft,1)
        end if
        write(*,*)'roots',pft,trm_ind
        
      else !annual C allocation
        trm_ind = rm_ind(pft,1) * 5._dp !fine roots are usually ~20% of the total root mass
      end if
      

    !Find the rooting depth

    !intialize root growth direction (0.80 value is from Arora & Boer)
      !if we wish to add in an ability for the allocation of resources to change dependent upon
      !something like PET, climate, or something else, we would change rootgro here. It is the
      !determinant of whether the roots go preferentially vertical or horizontally. 1 = vertical, 0 = horizontal
      !if try to put in climate, grasses are heavily swayed by mean annual precip, trees are not. See
      !Schenk & Jackson J Ecol. 2002 v90 480-494. However with the daily C allocation, this is effectively 
      !occuring.
    rootgro = 0.80_dp

    !put root biomass (gC/indiv) into kg organic/m2 of vegetated surface
    Broot = trm_ind * 1.e-3  / crownarea(pft) / OMorgC

    !find the rooting depth
    rootdepth(pft) = 3._dp * Broot**rootgro / littleb(pft)

        if (rootdepth(pft) > soildepth) then
        
          !if the rooting depth is deeper than the soil depth the root distribution profile is not allowed to deepen
          !with increasing root biomass. This is equivalent to reducing rootgro since to keep rootdepth constant at soil
          !depth it requires that soildepth = 3/abar (Broot/Bbar)exp (rootgro) hence rootgro =
          !ln(soildepth*abar/3)/ln(Broot/Bbar)
          
          !NOTE: the parameterization from the paper does not make much sense. It says that the rootgro value should 
          !decrease as the roots get to the max soil depth but the value actually increases, exponentially. It is easy 
          !to plot and check. Instead of that relation I am just setting alpha to 0 when the rootdepth is greater than
          !the soil depth. This makes all future growth be horizontal as you would expect (this assumes a constant soil
          !depth across a grid cell and no cracks in the bedrock for roots to exploit

            !rootgro  = log(soildepth * abar(pft) / 3._dp) / log(Broot / Bbar(pft)) Arora pub relation
            rootgro  = 0

            !now set the rooting depth to be the soil depth
            rootdepth(pft) = soildepth

        end if

     b = abar(pft) * Bbar(pft)**rootgro

     a = b / (Broot**rootgro)

    do l = 1,gnl

         ! These calculations can produce underflow errors, as they exponentially decline. 
         ! this is limited by these max statements below. JM 10.11.2010
         ! Update: This was likely a result of the problem with root depth once below soil
         ! depth, so likely not needed anymore. JM 17.03.2011
         
         !calculate the individual cumulative distribution profiles  (eqn 16)
         term_a = max(-50._dp,(-abar(pft) * (Bbar(pft) / Broot)**rootgro * zipos(l)))
         rootf(l) = 1._dp - exp(term_a)

         !find the root density at this soil depth  !Eqn 9
         term_b = max(-50._dp,(-a * zipos(l)))
         Bunit = exp(term_b)

         !Root density, can be useful for diagnostics but not strictly needed.
         !rootdens(l,pft) = b * Broot**(1._dp - rootgro) * Bunit !* 1000._dp !convert to g/m2

    end do !layers loop

        !make rootfraction per layer, not cumulative
         rootfracl(1,pft) = rootf(1)
         
         do l = 2,gnl
          rootfracl(l,pft) = rootf(l) - rootf(l-1)
         end do
         
        if (rootfracl(1,pft) == 1._dp) then  !don't allow all of root mass to be in the upper layer. Spread
                                                !it between the top 2
            rootfracl(1,pft) = 0.5
            rootfracl(2,pft) = 0.5
        end if    
                                                

  end if !present loop

end do !pft loop

return

end subroutine roots

!-------------------------------------------------------------------------------------------------------------


end module calcannual
