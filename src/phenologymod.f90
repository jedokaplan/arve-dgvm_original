module phenologymod

! Phenology subroutine for the daily C allocation scheme.
! Based upon Arora and Boer Global Change Biol. 2005
! Coded by JM 2010

use arveparams,                 only : dp,daysec,nts,Tfreeze,nts,npft,minsno
use statevars,                  only : sv,doy,dayl
use pftparametersmod,           only : prm
use metvarsmod,                 only : met,dm
use pft_state_variables,        only : veg

implicit none

public  :: phenology
private :: rainphen
private :: summerphen

contains

!---------------

subroutine phenology(j)

implicit none

! arguments
integer, intent (in) :: j       !gridcell index

! pointers
integer,  pointer, dimension(:)         :: phentype     !phenological type: evergreen (1), summergreen (2), raingreen (3), any type (4)
logical, pointer, dimension(:)          :: present      !true if PFT is present in gridcell.
real(dp), pointer, dimension(:,:)       :: npp_ra       !array of npp for last 7 days
integer,  pointer, dimension(:)         :: pstate       !instantaneous phenological state flag 
                                                        !(0=no leaves, 1=full leaf out,2=leafing out,3=senescing)
real(dp), pointer, dimension(:)         :: NPP7         !mean NPP over the last 7 days (gC m-2)
real(dp), pointer, dimension(:,:)       :: dnpp         !net primary productivity (gC m-2 timestep-1)
real(dp), pointer, dimension(:)         :: Tsoil        !soil temperature (K)
real(dp), pointer, dimension(:,:)       :: rootfracl    !root fraction in that soil layer
real(dp), pointer, dimension(:)         :: meanT_ra     !array of size nts holding mean air temp (C)
real(dp), pointer, dimension(:,:)       :: soilT_ra     !array of size nts holding mean soil temp in rooting zone (C)
real(dp), pointer                       :: tmin         !24 hour minimum temperature (C)
integer, pointer                        :: gnl          !index value of lowest 'soil' layer
real(dp), pointer, dimension(:)         :: gddbase      !base temperature for growing degree days
real(dp), pointer, dimension(:)         :: gdd_bioclim  !for use in bioclimatic
real(dp), pointer                       :: tmean        !mean temperature of 24 hour period (C)
       
! variables
real(dp) :: meanT_av             !average of nts for mean AIR temp.
integer  :: pft                  !plant functional type
real(dp) :: soilT_av             !average of nts for soil temp.
real(dp), dimension(npft) :: rootzone             !mean soil temperature in the rooting soil layers (K)

! point pointers
phentype        => prm%phentype          
present         => sv(j)%present         
npp_ra          => sv(j)%npp_ra          
NPP7            => sv(j)%NPP7            
meanT_ra        => sv(j)%meanT_ra       
rootfracl       => sv(j)%rootfracl
soilT_ra        => sv(j)%soilT_ra
Tsoil           => sv(j)%Tsoil 
tmin            => met(j,0)%tmin  !FLAG need mean?
gnl             => sv(j)%gnl
gddbase         => prm%gddbase
gdd_bioclim     => veg%gdd_bioclim
tmean           => met(j,0)%temp24
pstate          => sv(j)%pstate
dnpp            => veg%dnpp

!--------------

! begin calculations
  
  ! 24 hr mean surface air temperature
  meanT_ra(:) = eoshift(meanT_ra(:),1,tmin)
  meanT_av = sum(meanT_ra(:)) / real(nts)     

  do pft = 1,npft
        
        ! accumulate the gdd_bioclim for this pft
        gdd_bioclim(pft) = gdd_bioclim(pft) + max(tmean - gddbase(pft),0._dp)

       if (present(pft)) then
    
         ! add today's NPP to array
         npp_ra(pft,:) = eoshift(npp_ra(pft,:),1,dnpp(pft,1))
!write(0,*)doy,pft,dnpp(1)
         ! find running average across nts days
         NPP7(pft) = sum(npp_ra(pft,:)) / nts

         ! find average soil temperature in rooting zone
         rootzone(pft) = sum(Tsoil(1:gnl) * rootfracl(1:gnl,pft))
         soilT_ra(pft,:) = eoshift(soilT_ra(pft,:),1,rootzone(pft) - Tfreeze)
         
         soilT_av = sum(soilT_ra(pft,:)) / real(nts)

           !select type of phenology and calculate
           select case (phentype(pft))

             case(1)  !evergreen phenology
               pstate(pft) = 1

             case(2) !summergreen phenology
               call summerphen(j,pft,meanT_av,soilT_av)

             case(3) !raingreen phenology (treated as evergreens but with drought sensitivity)
               call rainphen(j,pft)

             case(4) !herbaceous PFTs act as raingreens
               call rainphen(j,pft)

           end select

        end if !if present
  end do !pft loop

end subroutine phenology

!----------------------------------------------------------------------------------------------------------------------

subroutine rainphen(j,pft)

! Raingreen phenology plants will act like evergreens when conditions are such that the plant keeps
! leaves on (i.e. leaf shedding due to drought or cold are not severe enough to defoliate). If the plant
! loses leaves to a lai of 0, then it behaves like a summerphen by putting all growth to leaves initially

implicit none

! arguments
integer, intent (in) :: j                       !gridcell index
integer, intent (in) :: pft                     !plant functional type

! pointers
logical, pointer                        :: virtual_leaf                 !true if the leaf being calculated is virtual
logical, pointer                        :: tree                         !true if tree
real(dp), pointer, dimension(:,:)       :: sm_ind                       !individual sapwood mass (gC)
integer,  pointer                       :: pstate                       !instantaneous phenological state flag
                                                                        ! (0=no leaves, 1=full leaf out,2=leafing out,3=senescing)
integer, pointer                        :: leafonday                    !day of year of leaf onset
integer, pointer                        :: leafoffday                   !day of year senescence complete
real(dp), pointer, dimension(:)         :: NPP7                         !
real(dp), pointer, dimension(:,:)       :: dnpp                         !net primary productivity (gC m-2 timestep-1)
real(dp), pointer, dimension(:,:)       :: dgpp                         !gross primary productivity (gC m-2 timestep-1)
real(dp), pointer                       :: lai_sno                      !
real(dp), pointer                       :: lai_vir                      !
real(dp), pointer                       :: lai_vir_no_burial            !
real(dp), pointer, dimension(:)         :: sla                          !
real(dp), pointer                       :: lai_max                      !
real(dp), pointer                       :: crownarea                    !
real(dp), pointer, dimension(:,:)         :: dresp                        !maintenance respiration (gC m-2 timestep-1)
real(dp), pointer                       :: lai_ind                      !
real(dp), pointer                       :: lm_ind                       !
real(dp), pointer                       :: needed_sap                   !the minimum amount of sapwood for the plant
                                                                        ! (at start of dormant season) gC indiv
real(dp), pointer                       :: sapturnover                  !sapwood turnover period (sapwood converted to heartwood) (years)
real(dp), pointer                       :: latosa                       !leaf area to sapwood area (m2 m-2)                       
real(dp), pointer, dimension(:,:)       :: lm_max                       !maximum individual leaf mass (gC)
real(dp), pointer                       :: height                       !plant height (m)     
real(dp), pointer                       :: leaf_thres                   !fraction of maximal leaf mass at which 
                                                                        ! the plant switches to full leafout allocation pattern
real(dp), pointer, dimension(:,:)       :: NSC_storage                  !individual plant non-structural carbon stores (gC)
real(dp), pointer, dimension(:)         :: wooddens                     !density of wood (sap + heart) (g m-3)

! local variables
real(dp) :: s_torate           !sapwood turnover rate(d-1)
real(dp) :: sm_turn            !sapwood mass turnover  (gC)
integer  :: pstate_new         !updated pstate after checking phenological conditions

! point pointers
leafonday         => sv(j)%leafonday(pft)
leafoffday        => sv(j)%leafoffday(pft)
pstate            => sv(j)%pstate(pft)
sm_ind            => sv(j)%sm_ind
NPP7              => sv(j)%NPP7
dnpp              => veg%dnpp
dgpp              => veg%dgpp
lai_sno           => sv(j)%lai_sno(pft)
lai_vir           => sv(j)%lai_vir(pft)
lai_vir_no_burial => sv(j)%lai_vir_no_burial(pft)
sla               => veg%sla
lai_max           => sv(j)%lai_max(pft)
crownarea         => sv(j)%crownarea(pft)
dresp             => veg%dresp
lai_ind           => sv(j)%lai_ind(pft)
lm_ind            => sv(j)%lm_ind(pft,1)
needed_sap        => sv(j)%needed_sap(pft)
sapturnover       => prm(pft)%sapturnover    
latosa            => prm(pft)%latosa
tree              => prm(pft)%tree
lm_max            => sv(j)%lm_max
height            => sv(j)%height(pft)
leaf_thres        => prm(pft)%leaf_thres
NSC_storage       => sv(j)%NSC_storage
virtual_leaf      => sv(j)%virtual_leaf(pft)
wooddens          => prm%wooddens
       
!---------------

! begin calculations

! initialize pstate_new to the present value
pstate_new = pstate

select case (pstate)

case(0)  ! no leaves
  
        ! check if NPP allows leaf out (FLAG I think I will need other moisture elements?)
        if (NPP7(pft) > 0._dp) then 

              ! can grow real leaves 
              virtual_leaf = .false.
              
              ! assign lai_vir to lai_sno and  lai_vir_no_burial to lai_ind
              lai_sno = lai_vir            ! set snow burial lai to lai_vir
              lai_ind = lai_vir_no_burial  ! set instant lai to lai_vir_no_burial 
              lai_vir = 0._dp              ! reset lai_vir and lai_vir_no_burial
              lai_vir_no_burial = 0._dp
              pstate_new = 2                   ! leafing out
              leafonday = doy                ! start leaf on day counter
              
              ! the leaf mass for this lai comes from the NSC_storage
              lm_ind = lai_sno * crownarea / sla(pft) 
              NSC_storage(pft,1) = NSC_storage(pft,1) - lm_ind  

               write(*,'(a,i4,5f12.5)')'grass, allowed to leaf out',doy,NPP7(pft),lai_sno,lai_ind,dgpp(pft,1),dnpp(pft,1)
               
        else !continue with no real leaves

            ! Conditions are not suitable for leaf out, keep no leaves state
              dgpp(pft,1) = 0._dp                   ! dgpp is 0 as only virtual leaves
              dnpp(pft,1) = dgpp(pft,1) - dresp(pft,1)      ! reassign npp assuming no gpp (no leaves)

        end if 

        
case(1)  ! full leaf
  
      ! if the leaves have been lost to an lai_ind of 0 then move to pstate 0 (no leaves).
      if (lai_ind <= 0._dp) then
      
      write(*,'(a,i5,3f12.4)')'grass lost back to 0',doy,lai_sno,lai_ind,sv(j)%zsno

         ! assign the virtual leaves for the raingreen trees/grasses
         virtual_leaf = .true.
           
         if (tree) then
           lai_vir_no_burial = min(0.3_dp, 0.075 * lai_max)
         else
           lai_vir_no_burial = min(0.2_dp, 0.075 * lai_max)
         end if  
                
         ! set pstate and assign the day of year to leafoffday              
         pstate_new = 0  !no leaves 
         leafoffday = doy
                  
         if (tree) then
           
          ! Set up the conversion of excess sapwood to heartwood for this coming 
          ! dormant season !FLAG or should they be treated like a evergreen????? JM 

            needed_sap = wooddens(pft) * height * lm_max(pft,1) * sla(pft) / latosa

          ! increment needed_sap for the amount of turnover that will occur during the dormant season 
          ! This assumes that dormant season will be as long as the previous year and that the tree
          ! will plan for sapwood turnover when setting how much sapwood to convert to heartwood.
            s_torate = 1._dp / (365._dp * sapturnover)  ! sapturnover is a tissue longevity value
            sm_turn = sm_ind(pft,1) * s_torate          ! daily turnover

            needed_sap = needed_sap + sm_turn * (leafonday + 365 - leafoffday)
            
         end if !tree
        
        !else continue at full leaf out
         
        end if     

case(2)  ! leafing out
  
        ! if the leaves have been gained back to leaf_thres of lai_max, move to pstate 1
        if (lai_ind > leaf_thres * lai_max) then 
        
        write(*,'(a,i4,4es12.4)')'grass moving to full leaf out',doy,NPP7(pft),lai_sno,lai_ind,lai_max
         
           pstate_new = 1 !full leaf out

        ! if leaves have been lost back to zero, restart leaf out when conditions permit.
        else if (lai_ind == 0._dp .or. NPP7(pft) < 0._dp) then  
         
           pstate_new = 0
           virtual_leaf = .true.
           
           if (tree) then   
              lai_vir_no_burial = min(0.3_dp, 0.075 * lai_max)
           else
              lai_vir_no_burial = min(0.2_dp, 0.075 * lai_max)  
           end if  !tree 
        
        !else continue leaf out
            
        end if  !lai>thres
        
case(3)  ! senescing

  ! There is assumed no senescence for rainphen. Assume that loss to drought and cold are how senescence occurs.

  ! plants are initally given a pstate of 3 to assign lai_vir_no_burial values
  ! now reset it to 0. (no leaves state)  
   pstate_new = 0
   
   !assign the virtual leaves
   if (tree) then
     lai_vir_no_burial = min(0.3_dp, 0.075 * lai_max)
   else
     lai_vir_no_burial = min(0.2_dp, 0.075 * lai_max)
   end if

   virtual_leaf = .true.
   leafonday = doy
   
end select

! update the pstate for changes
pstate = pstate_new

end subroutine rainphen

!----------------------------------------------------------------------------------------------------------------------

subroutine summerphen(j,pft,meanT_av,soilT_av)

implicit none

!arguments
integer, intent (in) :: j                               !gridcell index
integer, intent (in) :: pft                             !plant functional type
real(dp), intent(in) :: meanT_av                        !average of nts for min temp.
real(dp), intent(in) :: soilT_av                        !average soil temperature in the root zone over nts days

!pointers
logical, pointer                        :: virtual_leaf !true if the leaf being calculated is a virtual one to determine phenology.
real(dp), pointer, dimension(:,:)       :: sm_ind       !individual sapwood mass (gC)
real(dp), pointer, dimension(:,:)       :: rm_ind       !individual fine root mass (gC)
integer,  pointer                       :: pstate       !instantaneous phenological state flag 
                                                        !(0=no leaves, 1=full leaf out,2=leafing out,3=senescing)
integer, pointer                        :: leafonday    !day of year of leaf onset
integer, pointer                        :: leafoffday   !day of year senescence complete
real(dp), pointer, dimension(:)         :: NPP7         !
real(dp), pointer, dimension(:,:)       :: dnpp         !net primary productivity (gC m-2 timestep-1)
integer, pointer                        :: leaftype     !
real(dp), pointer                       :: lai_sno      !
real(dp), pointer                       :: lai_vir      !
real(dp), pointer                       :: lai_vir_no_burial      !
real(dp), pointer, dimension(:)         :: sla          !
real(dp), pointer                       :: lai_max      !
real(dp), pointer                       :: crownarea    !
real(dp), pointer, dimension(:,:)       :: NSC_storage  !individual plant non-structural carbon stores (gC)
real(dp), pointer, dimension(:,:)         :: dresp        !maintenance respiration (gC m-2 timestep-1)
real(dp), pointer                       :: lai_ind      !present individual leaf area index based on present lm_ind (m2 m-2)                  
real(dp), pointer                       :: lm_ind       !                                                                                               
real(dp), pointer                       :: needed_sap   !the minimum amount of sapwood for the plant (at start of dormant season) gC indiv              
real(dp), pointer                       :: sapturnover  !sapwood turnover period (sapwood converted to heartwood) (years)                               
real(dp), pointer                       :: latosa       !leaf area to sapwood area (m2 m-2)                                                             
real(dp), pointer, dimension(:,:)       :: lm_max       !maximum individual leaf mass (gC)
real(dp), pointer                       :: height       !plant height (m)  
real(dp), pointer, dimension(:,:)       :: dgpp         !gross primary productivity (gC m-2 timestep-1)   
real(dp), pointer                       :: leaf_thres   !fraction of maximal leaf mass at which the plant switches to full leafout allocation pattern
real(dp), pointer, dimension(:)         :: wooddens     !density of wood (sap + heart) (g m-3)

!local variables
real(dp) :: s_torate           !sapwood turnover rate(d-1)
real(dp) :: sm_turn            !sapwood mass turnover  (gC)
integer  :: pstate_new         !updated pstate after checking phenological conditions

!point pointers
leafonday         => sv(j)%leafonday(pft)
leafoffday        => sv(j)%leafoffday(pft)
pstate            => sv(j)%pstate(pft)
sm_ind            => sv(j)%sm_ind
rm_ind            => sv(j)%rm_ind
NPP7              => sv(j)%NPP7
dnpp              => veg%dnpp
leaftype          => prm(pft)%leaftype  
lai_sno           => sv(j)%lai_sno(pft)
lai_vir           => sv(j)%lai_vir(pft)
lai_vir_no_burial => sv(j)%lai_vir_no_burial(pft)
sla               => veg%sla
lai_max           => sv(j)%lai_max(pft)
crownarea         => sv(j)%crownarea(pft)
NSC_storage       => sv(j)%NSC_storage
dresp             => veg%dresp
lai_ind           => sv(j)%lai_ind(pft)
lm_ind            => sv(j)%lm_ind(pft,1)
needed_sap        => sv(j)%needed_sap(pft)
sapturnover       => prm(pft)%sapturnover    
latosa            => prm(pft)%latosa
lm_max            => sv(j)%lm_max
height            => sv(j)%height(pft)
dgpp              => veg%dgpp
leaf_thres        => prm(pft)%leaf_thres
virtual_leaf      => sv(j)%virtual_leaf(pft)
wooddens          => prm%wooddens
    
!--------------

! Calculations begin

! initialize pstate_new to the present value
pstate_new = pstate

select case (pstate)

case(0)  !no leaves
   
            ! Need a NPP7 > 0, daylength and above 0 air temps to be allowed to leaf out
            if (NPP7(pft) > 0._dp .and. dayl(1) > 39600._dp .and. meanT_av > 0._dp) then 
            
              virtual_leaf = .false.            !change flag to real eaves
              lai_sno = lai_vir                 !set snow burial lai to lai_vir
              lai_ind = lai_vir_no_burial       !set instant (no burial) lai to lai_vir_no_burial
              lai_vir = 0._dp                   !reset lai_vir and lai_vir_no_burial
              lai_vir_no_burial = 0._dp
              pstate_new = 2                    !leafing out
              leafonday = doy                     !set leafonday
              
              ! From the non-structural carbon store allow plants to develop leaves at leaf out at the same lai
              ! as lai_vir. The lm_ind allowed at the start of the season is calculated from the lai_vir
              ! value.        
       
              lm_ind = lai_sno * crownarea / sla(pft)   !since this is for trees, fine to leave as lai_sno
              
              ! this initial leaf carbon comes from the NSC storage
              NSC_storage(pft,1) = NSC_storage(pft,1) - lm_ind
              
                write(*,'(a20,2i5,4f12.4)')'tree-leafing begins',doy,pft,lai_sno,lm_ind,NPP7(pft),sv(j)%rm_ind(pft,1)

            else  !not ready to grow leaves
             
              !keep no leave state (reset gpp to zero as it was for virtual leaves)
              dgpp(pft,1) = 0._dp
              dnpp(pft,1) = dgpp(pft,1) - dresp(pft,1) !assign npp assuming no gpp (no leaves)
             
            end if
            
    
case(1)  !full leaf (conditions depend on leaf type)

       ! Needle leaves will senesce if avg min temp goes below -5.
       if (leaftype == 2 .and. meanT_av < -5._dp) then
                
                write(*,'(a20,2i5,f12.4)')'2snes needle senes-minT',doy,pft,meanT_av
                
         pstate_new = 3 !senescing

       ! Broad leaves are sensitive to light levels and soil temperature  !FLAG ADD IN NPP7<0 constraint?
       else if (leaftype == 1 .and. ((dayl(1) < 39600._dp .and. soilT_av < 11._dp) .or. (soilT_av < 2._dp .and. meanT_av < 0._dp))) then 
                
                write(*,'(a20,2i5,3f12.4)')'2sens broad-light,soil',doy,pft,dayl(1),soilT_av,meanT_av
                
         pstate_new = 3 !senescing

       !else continue with full leaf out 
       
       end if

case(2)  !leafing out

         ! if present LAI is leaf_thres of the max leafout then move to 'normal' allocation stage.
         if (lai_ind > leaf_thres * lai_max) then 
         
                write(*,'(a20,2i5,2f12.4)')'tree in full leaf',doy,pft,soilT_av,NPP7(pft)

           pstate_new = 1 !full leaf out

        ! if leaves have been lost back to zero, restart leaf out when conditions permit.
         else if (lai_ind == 0._dp) then !.or. NPP7(pft) < 0._dp) then  
         
           pstate_new = 0
           virtual_leaf = .true.

           !reassign the virtual leaves
              lai_vir_no_burial = min(0.3_dp, 0.075 * lai_max)
            
                write(*,'(a20,2i5,3f12.4)')'tree back2leafout',doy,pft,NPP7(pft),lai_sno,sv(j)%lm_ind(pft,1)
!               read(*,*)
         !else continue to leaf out      
                
         end if

case(3)  !senescing. 
     
        ! Allow for a return to normal (full leaf) state if the conditions improve while the leaves are still there
        ! this allows the plant to shed some leaves under stress but recover.
        
        ! Needleleaves are allowed to recover from a very cold weather period by putting needles back on.
        if (leaftype == 2 .and. meanT_av > -5._dp) then
                
                write(*,'(a20,2i5,1f12.4)')'tree broad-return',doy,pft,meanT_av 
                 
         pstate_new = 1 !back to full leaf out

        ! Broad leaves can recover if the daylength is long enough and temperature warms back up
        else if (leaftype == 1 .and. dayl(1) > 39600._dp .and. soilT_av > 2._dp .and. meanT_av > 0._dp) then 
       
                write(*,'(a20,2i5,3f12.4)')'tree broad-return',doy,pft,dayl(1),soilT_av,meanT_av

         pstate_new = 1 !back to full leaf out
         
        !else plant will continue sensesence 

        end if

        ! if the leaves have senseced to an lai of 0 then move to pstate 0 (no leaves).
        if (lai_ind <= 0._dp) then

           !assign the virtual leaves
           lai_vir_no_burial = min(0.3_dp, 0.075 * lai_max)
        
                write(*,'(a20,2i5,1f12.4)')'tree no leaves',doy,pft,lai_vir

           !set pstate and assign the day of year to leafoffday              
           pstate_new = 0  !no leaves (done senescence)
           virtual_leaf = .true.
           leafoffday = doy

           ! Set up the conversion of excess sapwood to heartwood for this coming 
           ! dormant season

           needed_sap = wooddens(pft) * height * lm_max(pft,1) * sla(pft) / latosa

           ! increment needed_sap for the amount of turnover that will occur during the dormant season 
           ! This assumes that dormant season will be as long as the previous year and that the tree
           ! will plan for sapwood turnover when setting how much sapwood to convert to heartwood.
            s_torate = 1._dp / (365._dp * sapturnover)  ! sapturnover is a tissue longevity value
            sm_turn = sm_ind(pft,1) * s_torate          ! daily turnover

            needed_sap = needed_sap + sm_turn * (leafonday + 365 - leafoffday)  
        
        !else plant continues to lose leaves
        
        end if  !laisno
        
   end select
      
! update the pstate for changes
pstate = pstate_new

end subroutine summerphen

end module phenologymod
