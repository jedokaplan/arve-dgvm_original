subroutine auto_resp(j,i,pft)

!SUBROUTINE auto_resp
!Calculation of maintenance and growth respiration and NPP
!recoded from f to f.90 by Joe Melton 2007 from original LPJ code (npp.f) and converted from monthly to twice-daily call.
!Also includes correction by Annette Wolf ETH to respcoeff and implements virtual leaves.

use arveparams,         only : dp,npft,ncvar,c12mass,n14mass,Tfreeze,daysec,minplant
use statevars,          only : sv,co2,dtime
use pftparametersmod,   only : prm
use pft_state_variables,only : veg
use metvarsmod,         only : met
use iovariables,        only : do_allocation_daily

implicit none

!ARGUMENTS
integer, intent(in)                     :: j                    !gridcell index
integer, intent(in)                     :: i                    !subdaily timestep
integer, intent(in)                     :: pft                  !pft

!pointers
logical, pointer, dimension(:)          :: virtual_leaf         !true if the leaf being calculated is a virtual one
integer, pointer                        :: gnl                  !index value of lowest 'soil' layer
real(dp), pointer, dimension(:,:)       :: dnpp                 !net primary productivity (gC m-2 timestep-1)
real(dp), pointer, dimension(:,:)       :: dgpp                 !gross primary productivity (gC m-2 timestep-1)
real(dp), pointer, dimension(:,:)       :: dresp                !maintenance respiration (gC m-2 timestep-1)
real(dp), pointer, dimension(:,:)       :: gresp                !growth respiration (gC m-2 timestep-1)
real(dp), pointer                       :: nind                 !individual density (indiv m-2)               
real(dp), pointer                       :: stem_sno             !stem area index (m2 m-2)                     
real(dp), pointer                       :: lai_sno              !leaf area index (m2 m-2)                     
real(dp), pointer                       :: Tveg_ave             !timestep average vegetation temperature (K)        
real(dp), pointer                       :: lm_ind               !individual leaf mass (gC)                          
real(dp), pointer                       :: sm_ind               !individual sapwood mass (gC)                       
real(dp), pointer                       :: rm_ind               !individual fine root mass (gC)                     
real(dp), pointer                       :: leafcn               !leaf C:N mass ratio                                
real(dp), pointer                       :: sapwoodcn            !sapwood C:N mass ratio                             
real(dp), pointer                       :: rootcn               !root C:N mass ratio                                
real(dp), pointer                       :: respcoef             !maintenance respiration coefficient                                        
real(dp), pointer                       :: temp                 !air temp (C)                                      
logical, pointer                        :: tree                 !true if tree                                      
real(dp), pointer, dimension(:,:)       :: rootfracl            !rootfraction per layer
real(dp), pointer, dimension(:)         :: Tsoil_ave            !day/night average soil temp (K)
real(dp), pointer                       :: lai_vir              !lai of virtual leaf (m2 m-2)
real(dp), pointer, dimension(:)         :: sla                  !PFT specific leaf area (m2/gC) 
real(dp), pointer                       :: fpc_grid             !instantaneous gridcell foliar projective cover 
                                                                ! (FPC) per pft for gridcell (m2 m-2)
real(dp), pointer                       :: phfrac               !instantaneous phenological state fraction 
                                                                ! (0=no leaves, 1=full leaf out)
real(dp), pointer                       :: fpc_gmax             !max instantaneous gridcell foliar projective cover 
                                                                ! (FPC) per pft for gridcell (m2 m-2)

!LOCAL VARIABLES
real(dp) :: lresp                       !leaf respiration 
real(dp) :: sresp                       !sapwood respiration
real(dp) :: rresp                       !fine root respiration
real(dp) :: gtemp_air                   ! value of temperature response function given air temperature
real(dp) :: gtemp_soil                  ! value of temperature response function given soil temperature
real(dp) :: k                           !
real(dp) :: rztC                        !temperature in the root zone average (C)
real(dp) :: leafT                       !temporary leaf temp.
real(dp) :: lm_temp                     !temp variable for leaf mass (virtual or real)

!parameter
real(dp), parameter :: st1 = 1._dp / 56.02_dp   !temp var.

!point pointers
gnl          => sv(j)%gnl
dgpp         => veg%dgpp
dnpp         => veg%dnpp
dresp        => veg%dresp
gresp        => veg%gresp
leafcn       => prm(pft)%leafcn
lm_ind       => sv(j)%lm_ind(pft,1)
nind         => sv(j)%nind(pft)
respcoef     => prm(pft)%respcoef
rm_ind       => sv(j)%rm_ind(pft,1)
rootcn       => prm(pft)%rootcn
sapwoodcn    => prm(pft)%sapwoodcn
sm_ind       => sv(j)%sm_ind(pft,1)
temp         => met(j,0)%temp(i)
tree         => prm(pft)%tree
stem_sno     => sv(j)%stem_sno(pft)
lai_sno      => sv(j)%lai_sno(pft)
Tveg_ave     => sv(j)%Tveg_ave(2,pft)
Tsoil_ave    => sv(j)%Tsoil_ave
rootfracl    => sv(j)%rootfracl
lai_vir      => sv(j)%lai_vir(pft)
sla          => veg%sla
fpc_grid     => veg%fpc_grid(pft)
virtual_leaf => sv(j)%virtual_leaf
phfrac       => sv(j)%phfrac(pft)
fpc_gmax     => sv(j)%fpc_gmax(pft)

!-----------------

! Formulae for calculation of maintenance respiration components in living tissue.
! DERIVATION appears below calculations.

gtemp_air = 0._dp
gtemp_soil = 0._dp

        ! set leaf temperature based on lai_sno
        if ((lai_sno + stem_sno) > minplant) then
          leafT = Tveg_ave - Tfreeze
        else
          leafT = temp !set to air temp
        end if
                       
        ! assign the leaf mass temp values 
        if (.not. virtual_leaf(pft)) then
            lm_temp = lm_ind
        else !virtual leaves
            lm_temp =  lai_vir * sv(j)%crownarea(pft) / sla(pft)  
        end if

        rztC = (sum(Tsoil_ave(1:gnl) * rootfracl(1:gnl,pft))) - Tfreeze
          
        if (rztC > -40._dp) gtemp_soil = exp(308.56_dp * (st1 - 1._dp / (rztC + 46.02_dp)))

        if (leafT >= -40._dp) gtemp_air  = exp(308.56_dp * (st1 - 1._dp / (leafT + 46.02_dp)))

        k = 1.28e-6 * C12mass / N14mass * dtime  !NOTE: this is new value! Annette Wolf @ ETHZ documents this

        ! Calculate tissue maintenance respiration values today [Eqn (7)] 
        ! Fluxes are expressed in units g m-2, pools are g ind-1. Multiply pool sizes by nind() to convert units
        if (tree) then
          sresp = respcoef * k * sm_ind / sapwoodcn * gtemp_air
        else   !grass
          sresp = 0._dp ! no sapwood 
        end if
        
        ! find the leaf and root respiration values
        lresp = respcoef * k * lm_temp / leafcn * gtemp_air 
        rresp = respcoef * k * rm_ind  / rootcn * gtemp_soil    
          
        if (do_allocation_daily) then
          
          if (tree) then
                  
            dresp(pft,1) = (lresp + rresp + sresp) * nind  
      
          else
          
          ! Grasses scale by fpc_grid. Grasses are modelled with a nind of 1 ('a green carpet')
          ! therefore to scale properly the resp, use the fpc_gmax (since the grass phenology
          ! is already taken care of by the amount of lm and rm. JM 03.03.2011
          dresp(pft,1) = (lresp + rresp) * nind * fpc_gmax
          
          end if
          
          gresp(pft,1) = 0._dp !FLAG test JM 22.03.2011
                      
        else !annual C allocation scheme
        
         if (tree) then
        
            ! annual C allocation scheme requires to account for phfrac of the veg
            ! for the roots and leaves
            dresp(pft,1) = ((lresp + rresp) * phfrac + sresp) * nind

          else
          
            ! grasses scale by fpc_grid. Grasses are modelled with a nind of 1 ('a green carpet')
            ! therefore to scale properly the resp, use the fpc_grid. JM 03.03.2011
            dresp(pft,1) = (lresp + rresp) * nind * fpc_grid  !FLAG need phfrac? JM 04.03.2011

          end if
          
        ! Incorporate growth respiration = 25% of (GPP - maintenance respiration)
        ! in daily NPP calculation (growth resp always non-negative)
        gresp(pft,1) = max(((dgpp(pft,1) - dresp(pft,1)) * 0.25_dp),0._dp)

        end if !daily/annual alloc.
        
        ! total carbon fixed (NPP) is then:
        dnpp(pft,1) = dgpp(pft,1) - dresp(pft,1) - gresp(pft,1)
        
!write(*,'(i4,6es12.4)')pft,dnpp(1),dgpp(1),dresp(1),gresp(1),lresp,rresp
        
        dresp(pft,2) = dgpp(pft,2)  !isotope ratio of maintenance respiration
        gresp(pft,2) = dgpp(pft,2)  !isotope ratio of growth resp.
        
        dnpp(pft,2) = dgpp(pft,2)
        dnpp(pft,3) = dgpp(pft,3)

        if (virtual_leaf(pft)) then !virtual leaves

            !reset the gresp values (no real growth...) 
            gresp(pft,1) = 0._dp
            
            !reset dresp, now assuming no contribution from leaves
            lresp = 0._dp
            
            if (tree) then                                    
               dresp(pft,1) = (lresp + sresp + rresp) * nind 
            else
               dresp(pft,1) = (lresp + rresp) * nind * fpc_gmax
            end if
            
            !total carbon fixed (NPP) is for virtual leaves is then:
            dnpp(pft,1) = dgpp(pft,1) - dresp(pft,1) - gresp(pft,1)

        end if

!Based on the relations
!(A) Tissue respiration response to temperature
!    (Sprugel et al 1995, Eqn 7)

!    (A1) Rm = 7.4e-7 * N * f(T)
!    (A2) f(T) = EXP (beta * T)

!      where Rm   = tissue maintenance respiration rate in mol C/sec
!            N         = tissue nitrogen in mol N
!            f(T) = temperature response function
!            beta = ln Q10 / 10
!            Q10  = change in respiration rate with a 10 K change
!                   in temperature
!            T         = tissue absolute temperature in K

!(B) Temperature response of soil respiration across ecosystems
!    incorporating damping of Q10 response due to temperature acclimation
!    (Lloyd & Taylor 1994, Eqn 11)

!    (B1) R = R10 * g(T)
!    (B2) g(T) = EXP [308.56 * (1 / 56.02 - 1 / (T - 227.13))]

!      where R         = respiration rate
!            R10  = respiration rate at 10 deg C
!            g(T) = temperature response function at 10 deg C
!            T         = soil absolute temperature in K

!Mathematical derivation:

!For a tissue with C:N mass ratio cton, and C mass, c_mass, N concentration
!in mol given by
!  (1) N = c_mass / cton / atomic_mass_N
!Tissue respiration in gC/day given by
!  (2) R = Rm * atomic_mass_C * seconds_per_day
!From (A1), (1) and (2),
!  (3) R = 7.4e-7 * c_mass / cton / atomic_mass_N * atomic_mass_C
!        * seconds_per_day * f(T)
!Let
!  (4) k = 7.4e-7 * atomic_mass_C / atomic_mass_N * seconds_per_day
!        = 0.0548
!from (3), (4)
!  (5) R = k * c_mass / cton * f(T)
!substituting ecosystem temperature response function g(T) for f(T)
!(Eqn B2),
!  (6) R = k * c_mass / cton * EXP [308.56 * (1 / 56.02 - 1 /
!        (T - 227.13))]
!incorporate PFT-specific respiration coefficient to model acclimation
!of respiration rates to average (temperature) conditions for PFT (Ryan
!1991)
!  (7) R_pft = respcoef_pft * k * c_mass / cton
!        * EXP [308.56 * (1 / 56.02 - 1 / (T - 227.13))]

end subroutine auto_resp
