subroutine daily_allocation(j)

!DAILY Vegetation carbon dynamics

use arveparams,only : dp,npft,ncvar,bet_com,alph_com,k5,k6,pi,&
                        OMorgC,frac_cr_tot,reinickerp,reprod_cost
use pftparametersmod, only : prm
use statevars, only : sv,doy
use metvarsmod, only : atm,met
use pft_state_variables, only : veg

implicit none

!arguments
integer, intent(in) :: j

!PARAMETERS  !FLAG! These will need sensitivity testing!
real(dp), parameter, dimension(npft) :: epsilon_s = [ 0.05, 0.05, 0.05, 0.05, 0.10, 0.05, 0.10, 0.00, 0.00 ]
real(dp), parameter, dimension(npft) :: epsilon_l = [ 0.40, 0.40, 0.40, 0.40, 0.35, 0.10, 0.35, 0.45, 0.45 ]
real(dp), parameter, dimension(npft) :: epsilon_r = [ 0.55, 0.55, 0.55, 0.55, 0.55, 0.85, 0.55, 0.55, 0.55 ]

real(dp), parameter, dimension(npft) :: par1 = [ 0.80, 0.80, 0.80, 0.80, 0.80, 0.80, 0.80, 0.80, 0.80 ] 
real(dp), parameter, dimension(npft) :: beta_lr = [ 0.76, 0.30, 0.76, 0.76, 0.30, 0.76, 0.30, 1., 1. ] !check these against those in pftparams
real(dp), parameter, dimension(npft) :: alpha_lr = [ 0.86, 0.76, 0.86, 0.86, 0.76, 0.86, 0.76, 1., 1. ] !check these against those in pftparams

real(dp), parameter, dimension(npft) :: Tcold = [ 5._dp, 5._dp, 5._dp, 5._dp, 5._dp, -25._dp, 0._dp, 0._dp, 0._dp ] 
real(dp), parameter, dimension(npft) :: cold_stress_max = [ 0.3, 0.3, 0.15, 0.15, 0.3, 0.15, 0.3, 0.3, 0.15 ] 
real(dp), parameter, dimension(npft) :: moist_stress_max = [ 0.005, 0.025, 0.005, 0.005, 0.025, 0.005, 0.005, 0.005, 0.005 ] 

real(dp), parameter :: b_W = 3._dp                            !
real(dp), parameter :: b_T = 3._dp                            !
real(dp), parameter :: sapkill_pc = 0.01_dp                   !percent of sapwood mass that can be converted to heartwood in 1 day
real(dp), parameter :: NSCboost_frac = 2.0              !fractional increase in bm_inc during leafout from NSC !FLAG

!pointers
real(dp), pointer, dimension(:,:) :: lm_ind                   !individual leaf mass (gC)
real(dp), pointer, dimension(:,:) :: lm_max                   !maximum individual leaf mass (gC)
real(dp), pointer, dimension(:,:) :: sm_ind                   !individual sapwood mass (gC)
real(dp), pointer, dimension(:,:) :: hm_ind                   !individual heartwood mass (gC)
real(dp), pointer, dimension(:,:) :: rm_ind                   !individual fine root mass (gC)
real(dp), pointer, dimension(:,:) :: rm_max                   !maximum individual fine root mass (gC)
real(dp), pointer, dimension(:,:) :: litter_ag_fast           !gridcell above-ground litter (gC/m2)
real(dp), pointer, dimension(:,:) :: litter_ag_slow           !gridcell above-ground litter (gC/m2)
real(dp), pointer, dimension(:,:) :: litter_bg                !gridcell below-ground litter (gC/m2)
real(dp), pointer, dimension(:)   :: nind                     !gridcell individual density (indiv/m2)  
real(dp), pointer, dimension(:)   :: leafturnover             !leaf turnover period (years)
real(dp), pointer, dimension(:)   :: sapturnover              !sapwood turnover period (sapwood converted to heartwood) (years)
real(dp), pointer, dimension(:)   :: frootQ10                 !fine root turnover Q10 exponent (Gill and Jackson 2000)
real(dp), pointer                 :: tmean                    !mean temperature of 24 hour period (C)
real(dp), pointer, dimension(:)   :: lai_ind                  !individual leaf area index                               
real(dp), pointer, dimension(:)   :: lai_sno                  !
logical, pointer , dimension(:)   :: tree                     !true if pft is tree                                      
logical, pointer,  dimension(:)   :: present                  !true if pft is present in that gridcell                  
integer, pointer, dimension(:)    :: pstate                   !
integer, pointer, dimension(:)    :: leaftype                 !type of leaf for plant
real(dp), pointer, dimension(:)   :: latosa                   !leaf area to sapwood area (m2 m-2)                       
real(dp), pointer, dimension(:)   :: crownarea                !crown area (m2)                                             
integer, pointer, dimension(:)    :: phentype                 !tree maximum crown area (m2)                             
real(dp), pointer, dimension(:)   :: sla                      !PFT specific leaf area (m2/gC)                           
real(dp), pointer, dimension(:)   :: fpc_ind                  !
real(dp), pointer, dimension(:)   :: height                   !tree height (m)                                          
real(dp), pointer, dimension(:,:) :: dnpp                     !net primary productivity (gC m-2 timestep-1)
real(dp), pointer, dimension(:)   :: lai_max                  !
real(dp), pointer, dimension(:)   :: dreprod                  !reproduction cost (gC/m2)
real(dp), pointer, dimension(:)   :: dturnover_ind            !daily total turnover of living biomass per individual (gC)
real(dp), pointer, dimension(:)   :: wetsum                   !Plant wilting factor (CLM 4.0 eq 8.18) weighted by root zone
real(dp), pointer, dimension(:)   :: lm_tot_max               !running max leaf mass for the year (gC/indiv)
real(dp), pointer, dimension(:)   :: rm_tot_max               !running max fine root mass for the year (gC/indiv)
real(dp), pointer, dimension(:,:) :: NSC_storage              !plant non-structural carbon stores (gC/indiv)
real(dp), pointer, dimension(:)   :: crownarea_max            !tree maximum crown area (m2) 
real(dp), pointer, dimension(:)   :: allom1                   !Packing constrain from Reinicke's rule exponent 
real(dp), pointer, dimension(:)   :: needed_sap               !the minimum amount of sapwood for the plant (at start of dormant season) gC indiv
real(dp), pointer, dimension(:)   :: frac                     !fraction of sapwood and heartwood that is below ground (in coarse roots)(gC)
real(dp), pointer, dimension(:)   :: NSC_limit                !allometric plant non-structural carbon stores size (gC/indiv)
real(dp), pointer, dimension(:)   :: NSC_opt                  !allometric plant non-structural carbon multiple of maximal leaf mass 
real(dp), pointer, dimension(:)   :: wooddens                 !density of wood (sap + heart) (g m-3)
real(dp), pointer                 :: MAT                      !mean annual temperature (C)
        
!LOCAL VARIABLES
integer  :: pft                 !counters
logical  :: needroots           !plant needs more root mass to satisfy allometric relationships
real(dp) :: a_leaves            !allocation fraction to leaves
real(dp) :: a_sapw              !allocation fraction to sapwood
real(dp) :: a_roots             !allocation fraction to roots
real(dp) :: a_leaves_new        !update to allocation fraction to leaves
real(dp) :: a_sapw_new          !update to allocation fraction to sapwood
real(dp) :: a_roots_new         !update to allocation fraction to roots
real(dp) :: bm_inc_ind          ! individual total biomass increment this year
real(dp) :: light_lev           !scalar index for measure of light availability
real(dp) :: soil_mois           !plant wilting factor (scalar index)
real(dp) :: roots2gro           !root mass (gC/ind) to grow to satisfy allometric relationships
real(dp) :: woodmass            !plant woody mass (kg dry OM / ind) includes rm,sm, and hm
real(dp) :: diff_root           !
real(dp) :: sap2gro             !
real(dp) :: sapwood_a           !
real(dp) :: cold_st             !
real(dp) :: temp_meas           !
real(dp) :: lm_loss             !
real(dp) :: moist_st            !
real(dp) :: temp                !temporary variable
real(dp) :: lm_temp             !temporary variable leaf
real(dp) :: rm_temp             !temporary variable root
real(dp) :: sm_temp             !temporary variable sapwood
real(dp) :: fine_roott          !
real(dp) :: l_torate            !leaf turnover rate (d-1)
real(dp) :: s_torate            !sapwood turnover rate(d-1)
real(dp) :: r_torate            !fine root turnover rate(d-1)
real(dp) :: sm_turn             !sapwood mass turnover  (gC)
real(dp) :: rm_turn             !fine root mass turnover  (gC)
real(dp) :: stem_diam           !plant stem diameter (m)
real(dp) :: a_sapw_temp         !
real(dp) :: a_roots_temp        !
real(dp) :: sapkill             !sapwood converted to heartwood this timestep (gC / ind)
real(dp) :: npp_st              !
real(dp) :: a_roots_inc         !
real(dp) :: tot_root            !total root mass (coarse and fine) (kg / ind)
real(dp) :: tot_leav            !maximum leaf mass for PFT this year (gC/ind)
logical  :: needwood            !plant needs wood production to satisfy allometric relationships
real(dp) :: needed_sap_i        !daily check if the minimum amount of sapwood is present (gC)
real(dp) :: allom_tot_leaves    !allometric relation for total leaf mass (in kg C / indiv) 
real(dp) :: NSC_need            !amount of carbon the NSC_storage is in deficit to the optimal size (gC)
real(dp) :: nsc_fill            !amount of carbon added to the nsc_storage for this day (gC)
real(dp) :: bm_inctemp          !temp variable to store bm_inc_ind

!point pointers
present           => sv(j)%present
lm_max            => sv(j)%lm_max
lm_ind            => sv(j)%lm_ind
sm_ind            => sv(j)%sm_ind
hm_ind            => sv(j)%hm_ind
rm_ind            => sv(j)%rm_ind
rm_max            => sv(j)%rm_max
litter_ag_fast    => sv(j)%litter_ag_fast
litter_ag_slow    => sv(j)%litter_ag_slow
litter_bg         => sv(j)%litter_bg
nind              => sv(j)%nind
leafturnover      => prm%leafturnover
sapturnover       => prm%sapturnover    
frootQ10          => prm%frootQ10
tmean             => met(j,0)%temp24
tree              => prm%tree
pstate            => sv(j)%pstate
crownarea         => sv(j)%crownarea
lai_ind           => sv(j)%lai_ind
lai_sno           => sv(j)%lai_sno
height            => sv(j)%height
sla               => veg%sla
latosa            => prm%latosa
phentype          => prm%phentype
lai_max           => sv(j)%lai_max
wetsum            => sv(j)%wetsum
NSC_storage       => sv(j)%NSC_storage
crownarea_max     => prm%maxcrowna
allom1            => prm%allom1
needed_sap        => sv(j)%needed_sap
lm_tot_max        => veg%lm_tot_max  
rm_tot_max        => veg%rm_tot_max  
frac              => veg%frac
NSC_opt           => prm%NSC_opt
NSC_limit         => veg%NSC_limit
fpc_ind           => veg%fpc_ind
leaftype          => prm%leaftype    
wooddens          => prm%wooddens
dreprod           => veg%dreprod    
dnpp              => veg%dnpp    
dturnover_ind     => veg%dturnover_ind
MAT               => atm%MAT

!---------------------------------

!calculations begin

do pft = 1,npft

 if (present(pft)) then
 
        !initial presets
        needroots = .false.
        needwood = .false.
        sapkill = 0._dp
        lm_loss = 0._dp  
        sm_turn = 0._dp  
        sapkill = 0._dp
        
!REPRODUCTION ----------------------------------------------------------------------------------------------

             ! Remove reproduction costs, which are taken simply as a constant fraction of NPP
               dreprod(pft) = max(dnpp(pft,1) * reprod_cost,0._dp)

                 ! Assume the costs go to reproductive structures which will
                 ! eventually enter the litter pool so add now.
     
                 temp = litter_ag_fast(pft,1)
                 litter_ag_fast(pft,1) = litter_ag_fast(pft,1) + dreprod(pft)

                 if (dreprod(pft) > 0._dp) then
                   litter_ag_fast(pft,2:ncvar) = (litter_ag_fast(pft,2:ncvar) * temp + dnpp(pft,2:ncvar)&
                                                 * dreprod(pft)) / (litter_ag_fast(pft,1))
                 end if

             ! Find biomass increment per individual. Reduce dnpp by reproduction costs and divide by indiv. density. 
               bm_inc_ind = (dnpp(pft,1) - dreprod(pft)) / nind(pft) 
               
             ! if the plant is in the max leafout stage, assist the leafout with C from NSC.
             ! the amount of C is a multiple (NSCboost_frac) of the bm_inc gained by photosynthesis.
             if (pstate(pft) == 2 .and. NSC_storage(pft,1) > 0.) then
             bm_inctemp = bm_inc_ind
             bm_inc_ind = bm_inctemp * NSCboost_frac
             NSC_storage(pft,1) = NSC_storage(pft,1) - (bm_inctemp * NSCboost_frac - bm_inctemp) 
             
             end if 
               
!write(*,'(a12,i5)')'=======  day',doy               
!write(*,'(a12,i5,es12.4)')'enter d a',pft,bm_inc_ind 
        
!NON-STRUCTURAL CARBON STORAGE ----------------------------------------------------------------------------

        ! Non-structural carbon stores either accumulate or deplete depending on the carbon
        ! source / sink status of the plant. If 1) the NSC_storage is less than the optimal amount
        ! and 2)the plant is not in full leaf out, and 3) the NPP is positive allow for contributions
        ! to the NSC pool at a rate that will replenish the pool after 30 days (rate from Thornton and Zimmerman 2007).
        ! if NPP is negative, the carbon can be depleted from the NSC pool. The optimal amount of NSC for each plant is 
        ! calculated as a proportion of the maximum leaf mass.
        ! the NSC pool is not considered to exist in any specific plant compartment and is assumed to be 
        ! available for relocation as needed. 
        
        ! Determine the present limit 
      
        NSC_limit(pft) = lm_max(pft,1) * NSC_opt(pft)
     
            if (NSC_storage(pft,1) < NSC_limit(pft) .and. bm_inc_ind > 0._dp .and. pstate(pft) /= 2) then  
            
            ! if NSC_storage is less than NSC_opt, use a portion of bm_inc_ind to refill at a rate
            ! that will fill it in 30 days, the remainder is used for growth
             
                NSC_need = NSC_limit(pft) - NSC_storage(pft,1)
                nsc_fill = NSC_need / 30._dp
                bm_inctemp = bm_inc_ind - nsc_fill
                
                if (bm_inctemp > 0._dp) then
                  bm_inc_ind = bm_inctemp
                else 
                  !NSC stores are quite depleted, all C must go to refilling the NSC stores
                  nsc_fill = bm_inc_ind
                  bm_inc_ind = 0._dp
                end if
                 
                NSC_storage(pft,1) = NSC_storage(pft,1) + nsc_fill
                
            else if (NSC_storage(pft,1) > NSC_limit(pft) .and. pstate(pft) /= 0) then
            
            ! NSC_storage is overfilled, use some of that excess carbon to grow. The amount is added to 
            ! bm_inc_ind at the same rate that was used to fill it up
              NSC_need = NSC_storage(pft,1) - NSC_limit(pft)
              nsc_fill = NSC_need !/ 3._dp !FLAG I made it so it is immediately returned to growth JM 22.03.2011
               write(*,*)doy,pft,'top up',bm_inc_ind,bm_inc_ind+nsc_fill,lm_max(pft,1)
              bm_inc_ind = bm_inc_ind + nsc_fill
              
              NSC_storage(pft,1) = NSC_storage(pft,1) - nsc_fill

                
            else if (bm_inc_ind < 0._dp) then

              !the negative bm_inc_ind value can be taken from NSC stores
              NSC_storage(pft,1) = NSC_storage(pft,1) + bm_inc_ind
              bm_inc_ind = 0._dp
                            
            end if 
            
!write(*,'(a12,3es12.4)')'after NSC',bm_inc_ind,NSC_storage(pft,1),NSC_limit(pft)

!write(*,'(a12,i5)')'enter al',pstate(pft)
!write(*,'(i5,a5,f12.4,a5,f12.4,a5,f12.4,a5,f12.4)')pft,'rm',rm_ind(pft,1),'sm',sm_ind(pft,1),'hm',hm_ind(pft,1),'lm',lm_ind(pft,1)
!write(*,*)

!ALLOCATION ----------------------------------------------------------------------------------------------

   ! determine each PFTs light and soil moisture constraints
       ! following Arora and Boer 2005 GCB
         if (tree(pft)) then
           light_lev = fpc_ind(pft)   
         else
           light_lev = max(0._dp, 1._dp - lai_sno(pft) / 4.5_dp)
         end if
   
         soil_mois = wetsum(pft) !FLAG Does this one do a good job?
                
    if (bm_inc_ind /= 0._dp) then  
    
      if (tree(pft)) then

            ! Using the day's bm_inc_ind, the allocation pattern depend on pstate.
            
            ! if the tree is in max leaf out stage, all growth to leaves
            if (pstate(pft) == 2) then
            
            !FLAG I might have to give some C to roots here, they appear to start growth about the same time as the 
            !leaves, or maybe lower the threshold to begin pstate 1, however I don't want shoot growth too early and 
            !usually fine roots grow before shoots. JM 24.03.2011.
               a_leaves = 0.5_dp !FLAG!! 4.4.2011
               a_sapw   = 0._dp
               a_roots  = 0.5_dp  !FLAG!! JM 4.4.2011
               
            ! if the tree is in normal growth stage (evergreens are always in this state), allocation is 
            ! split between the different compartments accounting for light and moisture constraints
            else if (pstate(pft) == 1) then

               a_sapw   = (epsilon_s(pft) + par1(pft) * (1._dp - light_lev)) / (1._dp + par1(pft) * (2._dp - light_lev - soil_mois))
               a_roots  = (epsilon_r(pft) + par1(pft) * (1._dp - soil_mois)) / (1._dp + par1(pft) * (2._dp - light_lev - soil_mois))
               a_leaves = 1._dp - a_sapw - a_roots
                          
               ! First Interm allocation
                 lm_temp = max(0._dp,lm_ind(pft,1) + (bm_inc_ind * a_leaves))
                 sm_temp = max(0._dp,sm_ind(pft,1) + (bm_inc_ind * a_sapw))
                 rm_temp = max(0._dp,rm_ind(pft,1) + (bm_inc_ind * a_roots)) 
                  
!write(*,*)'1st',pft,bm_inc_ind,a_leaves,a_roots,a_sapw
!write(*,'(a,a5,f12.4,a5,f12.4,a5,f12.4,a5,f12.4)')'first allocation','rm',rm_temp,'sm',sm_temp,'hm',hm_ind(pft,1),'lm',lm_temp
                  
               ! First, check that the present woody biomass is enough to support the leaf mass
               ! this includes all support material such as roots, sapwood and heartwood.
               ! units needed are kg dry organic matter / plant so convert from g C /plant
               ! woody mass to leaf mass relationship is from Niklas & Enquist, PNAS 2001 pg. 2923, fig 3A
               woodmass = (rm_temp + sm_temp + hm_ind(pft,1)) * 0.001_dp / OMorgC

               if (woodmass < ((lm_temp * 0.001_dp / OMorgC / bet_com)**(1./ alph_com))) then
                 
                  ! woody biomass is then not sufficient to support leaf mass so
                  ! shift allocation away from leaves to woody biomass. This is done by
                  ! first checking which needs more allocation (either sapwood or fine roots or both)
                  ! and then changing the allocation, which is done in the next 2 checks. 

                   needwood = .true.
                   
               end if
               
                    ! Now check that the fine root biomass is sufficient to supply this leaf and stem mass
                    ! if not, then increase the allocation to fine roots at expense of leaves (if needwood) or
                    ! possibly stems.

                    ! This will use Enquist and Niklas Science 2002 Table 1 relationship which is for total
                    ! root mass. rm_ind is only fine root mass. In ARVE (and LPJ), coarse root mass is implicitly considered part
                    ! of sapwood and heartwood compartments. So first need to determine how much of the sapwood 
                    ! and heartwood are underground (i.e. coarse roots) and how much is above ground (i.e. stems)

                    ! Since the plant structure will be structure for periods of full leafout, compare the plant
                    ! structure for the frac calculation with maximum leaf mass (not instantaneous)
                    tot_leav = max(lm_max(pft,1), lm_temp)
                    
                    ! the fraction of sapw and heartw below ground is based on the assumption that coarse roots 
                    ! are 20% of total aboveground biomass (see Ryan, Binkley & Formes Adv.Ecol. Res. 1997)
                    frac(pft) = (frac_cr_tot * (tot_leav + hm_ind(pft,1) + sm_temp)) / (hm_ind(pft,1) + sm_temp)
    
                    ! thus total root mass (kg dry wt. / indiv) is:
                    tot_root = 0.001 / OMorgC * (rm_temp + frac(pft) * (sm_temp + hm_ind(pft,1)))
                    
                    ! the allometric relation for total leaf mass (in kg dry wt. / indiv) is then:
                    allom_tot_leaves = (lm_temp * 0.001 / OMorgC / beta_lr(pft))**(1./ alpha_lr(pft))
                                     
                    if (tot_root < allom_tot_leaves) then

                         needroots = .true.
                        
                        ! determine how much TOTAL roots to grow (gC /ind) based upon the plant compartments 
                        ! after the first interm allocation, take it principally from leaves IF needwood = false
                         roots2gro = (allom_tot_leaves - tot_root) / 0.001 * OMorgC
                         
                         a_roots_inc = min(1._dp, ABS(roots2gro / bm_inc_ind)) 
                                                   
                         ! a_roots_inc can then be divided up into the fine and coarse roots
                         ! FLAG this needs to be fixed! JM 22.02.2011                                                  
                         a_roots_new  = min(1._dp, a_roots_inc + a_roots)
                         diff_root    = a_roots_new - a_roots 
                           
                           if (needwood) then  !the C comes from the leaves allotment

                            a_leaves_new = a_leaves - diff_root

                               if (a_leaves_new < 0._dp) then  ! if the amount allocated to leaves is less than 0, 
                                                               ! reduce the allocation to sapw                         
                                   a_sapw_new   = max(a_sapw + a_leaves_new, 0._dp)
                                   a_sapw       = a_sapw_new
                                   a_leaves_new = max(0._dp, a_leaves_new)

                               end if

                                a_roots =  a_roots_new
                                a_leaves = a_leaves_new

                           else  !comes from the sapwood allotment

                            a_sapw_new = a_sapw - diff_root  

                               if (a_sapw_new < 0._dp) then  ! if the amount allocated to leaves is less than 0, reduce the allocation
                                                             ! to sapw
                                   a_leaves_new   = max(a_leaves + a_sapw_new, 0._dp)
                                   a_leaves       = a_leaves_new
                                   a_sapw_new = max(0._dp, a_sapw_new)

                               end if
                               
                                a_roots =  a_roots_new
                                a_sapw = a_sapw_new
                                
                            end if
                       
                              !Second Interm allocation
                              lm_temp = max(0._dp,lm_ind(pft,1) + (bm_inc_ind * a_leaves))
                              sm_temp = max(0._dp,sm_ind(pft,1) + (bm_inc_ind * a_sapw))
                              rm_temp = max(0._dp,rm_ind(pft,1) + (bm_inc_ind * a_roots))

                    end if  !root check

                    ! Now check if sapwood mass is sufficient
                      
                    needed_sap_i = wooddens(pft) * height(pft) * lm_temp * sla(pft) / latosa(pft)

                    if (sm_temp < needed_sap_i) then  !tree needs more sapwood mass

                        if (needroots) then   ! the plant already needs more roots so 1) make a_sapw = 0, and if
                                              ! that is not enough, 2) kill some sapwood to heartwood
                          
                              ! First try, from present leaf and sapwood mass, determine how much more sapwood is required
                              sap2gro = needed_sap_i - sm_temp
                              sapwood_a = min(1._dp, ABS(sap2gro / bm_inc_ind))  
                              
                              ! if the needed sapwood allocation can come from a_leaves then do that
                              if (sapwood_a < (a_sapw + a_leaves)) then

                                a_sapw_new = min(1._dp, a_sapw + sapwood_a)
                                a_leaves = max(0._dp,a_leaves - (a_sapw_new - a_sapw))
                                a_sapw = a_sapw_new
                                
                              else !need to kill some sapwood to heartwood.
                              
                                ! assume that a max percent (sapkill_pc) of the sapwood can be converted to heartwood in one time step
                                ! the amount killed is the deficit of wood to grow minus that grown this round (the 
                                ! sapwood and leaves allocation).
                                sapkill  = min(sapkill_pc, (sap2gro - (a_sapw + a_leaves) * ABS(bm_inc_ind)) / sm_ind(pft,1))

                                a_leaves = 0._dp
                                a_sapw = max(0._dp,1._dp - a_roots)

                              end if
                              
                        else  ! has enough roots, needed C can come from roots and leaves
                        
                              ! From present leaf and sapwood mass, determine how much more sapwood is required                      
                              sap2gro = needed_sap_i - sm_temp
                             
                              sapwood_a = min(1._dp, ABS(sap2gro / bm_inc_ind))
                              a_sapw_new = min(1._dp, a_sapw + sapwood_a)
                              
                              ! reduce a_roots and a_leaves proportionally
                              a_roots = max(0._dp,a_roots - ((a_sapw_new - a_sapw) * a_roots / (a_roots + a_leaves)))
                              a_leaves = max(0._dp, 1._dp - a_roots - a_sapw_new)
                              a_sapw = a_sapw_new
                              
                        end if  !needroots check within sapwood check

                    !Third Interm allocation
                    lm_temp = max(0._dp,lm_ind(pft,1) + (bm_inc_ind * a_leaves))
                    sm_temp = max(0._dp,sm_ind(pft,1) + (bm_inc_ind * a_sapw) - (sapkill * sm_ind(pft,1)))
                    rm_temp = max(0._dp,rm_ind(pft,1) + (bm_inc_ind * a_roots))

                    end if   !sapwood check

            ! if the tree is senescing, all growth only to stems or roots
            else if (pstate(pft) == 3) then
            
            !FLAG I wonder if this should only be to shoots? JM 24.03.2011

                a_sapw_temp = (epsilon_s(pft) + par1(pft) * (1._dp - light_lev)) / (1._dp + par1(pft) * (2._dp - light_lev - soil_mois))
                a_roots_temp = (epsilon_r(pft) + par1(pft) * (1._dp - soil_mois)) / (1._dp + par1(pft) * (2._dp - light_lev - soil_mois))
                
                a_sapw = a_sapw_temp / (a_sapw_temp + a_roots_temp)
                a_roots = a_roots_temp / (a_sapw_temp + a_roots_temp)
                a_leaves = 0._dp
               
            !no leaves, no allocation
            else if (pstate(pft) == 0) then

               a_leaves = 0._dp
               a_sapw = 0._dp
               a_roots = 0._dp

            end if  !pstate

            !Perform final allocation
             lm_ind(pft,1) = max(0._dp,lm_ind(pft,1) + (bm_inc_ind * a_leaves))
             sm_ind(pft,1) = max(0._dp,sm_ind(pft,1) + (bm_inc_ind * a_sapw) - (sapkill * sm_ind(pft,1)))
             rm_ind(pft,1) = max(0._dp,rm_ind(pft,1) + (bm_inc_ind * a_roots))
             hm_ind(pft,1) = hm_ind(pft,1) + (sapkill * sm_ind(pft,1))
      
      else  !grass allocation

             ! if the grass is in max leaf out stage, all growth to leaves
            if (pstate(pft) == 2) then

        !FLAG same questions as for the trees, should fine roots be given more C earlier? JM 24.03.2011
               a_leaves = 1._dp
               a_roots = 0._dp

            ! if the grass is in normal growth stage
            else if (pstate(pft) == 1) then

               a_roots = (epsilon_r(pft) + par1(pft) * (1._dp - soil_mois)) / (1._dp + par1(pft) * (1._dp + light_lev - soil_mois))
               a_leaves = (epsilon_l(pft) + par1(pft) * light_lev) / (1._dp + par1(pft) * (1._dp + light_lev - soil_mois))

            ! if the grass is senescing, all growth only to roots
            else if (pstate(pft) == 3) then

        !FLAG would grass do this? Or would they put growth to making a big NSC store? JM 24.03.2011
               a_roots = 1._dp
               a_leaves = 0._dp
               
            !no leaves, no allocation
            else if (pstate(pft) == 0) then

               a_leaves = 0._dp
               a_roots = 0._dp

            end if
            
            ! Peform final grass allocation
             lm_ind(pft,1) = max(0._dp,lm_ind(pft,1) + (bm_inc_ind * a_leaves))
             rm_ind(pft,1) = max(0._dp,rm_ind(pft,1) + (bm_inc_ind * a_roots))

      end if !tree/grass
    end if !bm_inc_ind/=0
      
!write(*,'(a8,a5,f12.4,a5,f12.4,a5,f12.4,a5,f12.4)')'post al','rm',rm_ind(pft,1),'sm',sm_ind(pft,1),'hm',hm_ind(pft,1),'lm',lm_ind(pft,1)

!LEAF LITTER ----------------------------------------------------------------------------------------------   

     ! Perform generation of leaf litter through normal turnover, moisture/drought stress, cold stress
     if (lm_ind(pft,1) > 0._dp .or. pstate(pft) /= 2) then  !if there is some leaf mass and it is not in full leaf out
     
     cold_st = 0._dp !inital presets
     npp_st = 0._dp
     moist_st = 0._dp
 
       if (phentype(pft) /= 1) then ! evergreens are subject to only normal leaf turnover. Deciduous have also
                                    ! cold, moisture, and an npp stress component
                                    
         ! Normal leaf turnover
         ! decididuous assumes that normal leaf turnover is accomplished by moisture, cold and NPP stress.
         !FLAG should grass have a normal leaf turnover?? JM 22.03.2011
         l_torate = 0._dp
                           
         ! Moisture stress
         moist_st = moist_stress_max(pft) * (1._dp - soil_mois)**b_W

        ! Grasses are assumed only sensitive to moisture stress. Trees are also sensitive to 
        ! cold and NPP
          if (leaftype(pft) /= 3) then  !trees

                ! Cold stress
                ! using 24 hr mean air temperature (C)
                    if (tmean >= Tcold(pft)) then

                      temp_meas = 1._dp

                    else if (Tcold(pft) > tmean .and. tmean > (Tcold(pft) - 5._dp)) then

                      temp_meas = (tmean - (Tcold(pft) - 5._dp)) / 5._dp

                    else

                      temp_meas = 0._dp

                    end if      

                cold_st = cold_stress_max(pft) * (1._dp - temp_meas)**b_T

                ! If the dnpp(1) is negative this is assumed to dramatically increase leaf loss in deciduous 
                ! PFTs. Evergreens are assumed to remain only sensitive to drought and cold stress.
                if (dnpp(pft,1) < 0._dp) then   
                    npp_st = (-dnpp(pft,1)) ** 0.2_dp
                end if    
         
           end if !tree vs grass
 
         else !evergreens

           ! Normal leaf turnover 
           ! turnover rates (per day) are reciprocals of tissue longevity (can be pre-calc'd elsewhere)
           l_torate = 1._dp / (365._dp * leafturnover(pft))

         end if !evergreen/deciduous
 
           ! Sum leaf loss rate due to turnover, cold, negative NPP, and drought (grass only turnover and drought!)
           lm_loss = lm_ind(pft,1) * (l_torate + moist_st + cold_st + npp_st)

           !update the leaf C 
           lm_ind(pft,1) = lm_ind(pft,1) - lm_loss

           if (lm_ind(pft,1) < 1.e-4) then

              ! since the leaf loss term takes a proportion of the lm_ind, it will never go to 0 unless
              ! forced. Therefore set to 0 with lm_ind below 1.e-4 and give to litter 
              lm_loss = lm_loss + lm_ind(pft,1)
              lm_ind(pft,1) = 0._dp

           end if

           !Transfer killed leaves to litter pools
           temp = litter_ag_fast(pft,1)
           litter_ag_fast(pft,1) = litter_ag_fast(pft,1) + lm_loss * nind(pft)

           if (litter_ag_fast(pft,1) > 0._dp) then

                litter_ag_fast(pft,2:ncvar) = (litter_ag_fast(pft,2:ncvar) * temp &
                                              + lm_ind(pft,2:ncvar) * lm_loss * nind(pft)) / litter_ag_fast(pft,1)
           end if
           
        end if !leaf mass and not in full leaf out.
        
!ROOTS TURNOVER  ----------------------------------------------------------------------------------------------
       
          !Normal fine root turnover.
          
                ! 23.03.2011 JM. Okay new way of looking at fine root turnover. I am looking at Norby et al. 2004 PNAS, it is one of the few
                ! fine root turnover references that shows the seasonal cycle of fine root production and mortality. Most show only annual
                ! values. Of course this data is limited by only being for the FACE locations (temperate). Assume that fine root turnover is dependent upon
                ! the leaf status. Since leaves require water due to transpiration, leaf mortality is small when the leaves
                ! are actively growing. Once they are established the turnover can increase and should respond to temperature
                ! since the respiration from the fine roots is strongly sensitive to temperature. Essentially with at least temperate
                ! deciduous fine roots they should get quite depleted by the late winter. This is to reduce the costs associate with the fine roots
                ! while the plant is not actively photosynthesizing. 
                
                !if needroots, then no enhanced turnover, simple baseline turnover
                if (needroots .or. pstate(pft) == 2 .or. soil_mois < 0.8) then  !FLAG what about when water supplies are low???
             
               ! New way of calculating fine root turnover based upon Gill and Jackson 2000 New Phytologist, 147 13 -31)
               ! Based upon a global database (with some noted holes), they found 'Q10' relations for fine root turnover for grasslands, shrublands,
               ! and forest fine roots. The relation is dependent upon mean annual temperature (C). JM Dec 1 08
               fine_roott = frootQ10(pft) ** (0.1 * (MAT - 25._dp))  

               ! Turnover rates (per day) are reciprocals of tissue longevity (these can be pre-calc'd)
               r_torate = 1._dp / 365._dp * fine_roott          ! fine_roott is a tissue turnover value

              else !enhanced turnover
               
               ! turnover is enhanced when the soil temperatures are elevated or when the cost to plant is excessive (leaves are gone or
               ! seneseced. 
               
               if (pstate(pft) == 0 .or. pstate(pft) == 3) then
               r_torate = 0.5
               
               else
                ! New way of calculating fine root turnover based upon Gill and Jackson 2000 New Phytologist, 147 13 -31)
               ! Based upon a global database (with some noted holes), they found 'Q10' relations for fine root turnover for grasslands, shrublands,
               ! and forest fine roots. The relation is dependent upon mean annual temperature (C). JM Dec 1 08
               fine_roott = frootQ10(pft) ** (0.1 * (MAT - 25._dp))  

               ! Turnover rates (per day) are reciprocals of tissue longevity (these can be pre-calc'd)
               r_torate = 1._dp / 365._dp * fine_roott          ! fine_roott is a tissue turnover value
               
               end if
               
               
               end if
             
    !  write(*,'(a,2i4,l4,i4,3f12.4)')'root turnover',doy,pft,needroots,pstate(pft),soil_mois,r_torate,rm_ind(pft,1) 
       
               !Calculate the biomass turnover in this day
               rm_turn = rm_ind(pft,1) * r_torate

               !Update the pools
               rm_ind(pft,1) = rm_ind(pft,1) - rm_turn

               !Transfer root turnover to litter pools
               temp = litter_bg(pft,1)
               litter_bg(pft,1) = litter_bg(pft,1) + rm_turn * nind(pft)

               if (litter_bg(pft,1) > 0._dp) then

                     litter_bg(pft,2:ncvar) = (litter_bg(pft,2:ncvar) * temp + rm_ind(pft,2:ncvar) * rm_turn * nind(pft)) / litter_bg(pft,1)

               end if

! NORMAL and EXCEPTIONAL SAPWOOD TURNOVER---------------------------------------------------------------------------------------

        ! Assume: Evergreens have constant sapwood turnover, deciduous do not.
        
        ! if the plant has more sapwood than required for peak leaf mass, then kill some of the sapwood to
        ! heartwood. This can occur past peak leafout and continues at the rate of sapkill until sapwood mass
        ! has reached an appropriate amount for max leaf mass
        
        ! Deciduous:
        ! Taylor et al. 2002 Wood and Fibre Science suggest from a literature review that the majority of
        ! sapwood to heartwood conversion occurs during the dormant season. This will apply to 
        ! deciduous, for evergreens, assume that conversion will occur at a constant rate depending on how
        ! much excess sapwood exists.  
 
          ! Evergreens have constant sapwood turnover
          if (tree(pft) .and. phentype(pft) == 1) then

               s_torate = 1._dp / (365._dp * sapturnover(pft))  ! sapturnover is a tissue longevity value

               sm_turn = sm_ind(pft,1) * s_torate

               sm_ind(pft,1) = sm_ind(pft,1) - sm_turn

               ! Convert turned over sapwood to heartwood
               temp = hm_ind(pft,1)
               hm_ind(pft,1) = hm_ind(pft,1) + sm_turn

               if (hm_ind(pft,1) > 0._dp) then
                     hm_ind(pft,2:ncvar) = (hm_ind(pft,2:ncvar) * temp + sm_ind(pft,2:ncvar) * sm_turn) / hm_ind(pft,1)
               end if

          end if

         ! Deciduous in the dormant season or evergreen anytime can have exceptional sapwood turnover due to excess sapwood
         ! FLAG should it only be dormant season? JM 23.03.2011
          if (tree(pft) .and. ((pstate(pft) == 0 .and. phentype(pft) /= 1) .or. phentype(pft) == 1)) then  

               ! Find amount of sap needed.
               if (phentype(pft) == 1) then  !for evergreens
                 needed_sap(pft) = sm_ind(pft,1)!wooddens(pft) * height(pft) * lm_max(pft,1) * sla(pft) / latosa(pft)  !FLAG!
               else  !non-evergreen  
                 needed_sap(pft) = wooddens(pft) * height(pft) * lm_max(pft,1) * sla(pft) / latosa(pft)
               end if  

                 if (sm_ind(pft,1) > needed_sap(pft)) then  

                    ! assume that a max percent (sapkill_pc) of the sapwood can be converted to heartwood in one time step
                    ! the amount killed is then the smaller of  sapkill_pc*sm_ind and the difference between the amount 
                    ! required and present sm_ind
                    sapkill  = min(sapkill_pc * sm_ind(pft,1), sm_ind(pft,1) - needed_sap(pft))

                    hm_ind(pft,1) = hm_ind(pft,1) + sapkill
                    sm_ind(pft,1) = sm_ind(pft,1) - sapkill

                 end if  
          end if    

         !Record total turnover (used in mortality)
            dturnover_ind(pft) = lm_loss + rm_turn + sm_turn + sapkill  !FLAG check on this use of sapkill

         !Set isotopes to zero if pools are empty
            if (lm_ind(pft,1) <= 0._dp) lm_ind(pft,1:ncvar) = 0._dp
            if (sm_ind(pft,1) <= 0._dp) sm_ind(pft,1:ncvar) = 0._dp
            if (rm_ind(pft,1) <= 0._dp) rm_ind(pft,1:ncvar) = 0._dp

! READJUST LAI, HEIGHT, CROWNAREA,FPC for changes---------------------------------------------------------------------------------------

              !adjust the running max lm and rm
              lm_tot_max(pft) = max(lm_tot_max(pft),lm_ind(pft,1))  
              rm_tot_max(pft) = max(rm_tot_max(pft),rm_ind(pft,1)) 
                            
              !if lm/rm_tot_max is now larger than lm/rm_max, increase lm/rm_max accordingly
              lm_max(pft,1) = max(lm_tot_max(pft),lm_max(pft,1))
              rm_max(pft,1) = max(rm_tot_max(pft),rm_max(pft,1))
              
              !adjust for changes to lai
              lai_max(pft) = lm_max(pft,1) * sla(pft) / crownarea(pft)
              lai_ind(pft) = (lm_ind(pft,1) * sla(pft)) / crownarea(pft)
   
       if (tree(pft) .and. sm_ind(pft,1) > 0._dp) then
       
! JM 21.02.2011
! Try new stem diam calculation from Niklas and Spatz 2004          
!         stem_diam = (4.d0 * (hm_ind(pft,1) + sm_ind(pft,1)) * lm_max(pft,1) * sla(pft) / (pi * sm_ind(pft,1) * latosa(pft)))**0.5 
!
              stem_diam = ((sm_ind(pft,1) + hm_ind(pft,1)) * 1.e-3 /(202.3 * k5))**(0.375)
!--                 
              height(pft) = max(0.05,k5 * stem_diam**(2. / 3.) - k6) 
              crownarea(pft) = min(allom1(pft) * stem_diam**reinickerp,crownarea_max(pft))  !LPJ relation
              
       else !grass   
            
               !use CTEM (Arora and Boer GBC 2005) grass height (m) (expects leaf mass in kg)
               ! set a minimum height of 0.05 m, for resistance roughness length.
               height(pft) = max(0.05, 3.5_dp * (lm_ind(pft,1) * 0.001_dp)**0.5_dp) 

       end if !tree/grass

!write(*,'(a,2i5,l5,i5,6f12.4)')'exit dail all',doy,pft,tree(pft),pstate(pft),NSC_limit(pft),NSC_storage(pft,1),lai_sno(pft),bm_inc_ind,rm_max(pft,1),rm_ind(pft,1)
!write(*,'(i5,a5,f12.4,a5,f12.4,a5,f12.4,a5,f12.4,a5,f12.4)')pft,'rm',rm_ind(pft,1),'sm',sm_ind(pft,1),'hm',hm_ind(pft,1),'lm',lm_ind(pft,1),'NSC',NSC_storage(pft,1)
!write(*,*)
!if (doy == 185) read(*,*)  
  
  end if !pft present
end do !pft

return

end subroutine daily_allocation
