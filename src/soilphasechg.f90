subroutine soilphasechg(j)
!calculates the energy budget of phase change for the soil water constituents

use arveparams, only : dp,ns,nl,Lf,Tfreeze,grav,oneminusCNfac,CNfac,pliq,pice,z0mg
use statevars, only : sv,dt
use soilstate_vars, only : hs,dhsdt,Kl,qsnomelt,ithaw,fact,Fhti,surf

implicit none

!arguments
integer, intent(in) :: j    !index of the current gridcell

!implicit arguments
integer, pointer :: snl     !index value of top snow layer (negative is more layers)
real(dp), pointer :: Wsno   !snow water equivalent of the snowpack (mm)
real(dp), pointer :: zsno   !total thickness of the snowpack (m)
real(dp), pointer, dimension(:) :: Wliq    !soil liquid water content at layer midpoint (mm)
real(dp), pointer, dimension(:) :: Wice    !soil ice content at layer midpoint (mm)
real(dp), pointer, dimension(:) :: Tsoil   !soil temperature (K)
real(dp), pointer, dimension(:) :: Tnsoi   !updated soil temperature at timestep n+1
real(dp), pointer, dimension(:) :: Tsat    !soil water volumetric water content at saturation (fraction)
real(dp), pointer, dimension(:) :: Bexp    !soil water B exponent used in the Brooks & Corey Pedotransfer Functions
real(dp), pointer, dimension(:) :: Psat    !soil water matric potential at saturation (mm)
real(dp), pointer, dimension(:) :: zpos    !z coordinate of the middle of the soil layer relative to soil surface (m), positive downwards
real(dp), pointer, dimension(:) :: dz      !thickness of the soil layers (m)
real(dp), pointer, dimension(:) :: dzmm    !thickness of the soil layers (mm)
real(dp), pointer, dimension(:) :: Tice           !soil ice content (fraction)
real(dp), pointer, dimension(:) :: Tpor    !soil volumetric porosity (fraction)
real(dp), pointer, dimension(:) :: Tliq    !soil liquid water content (fraction)
real(dp), pointer, dimension(:) :: fice0   !layer ice fraction, previous timestep
real(dp), pointer               :: fsnow   !fraction of the gridcell covered by snow (fraction)

!parameters
real(dp), parameter :: psno_new        = 100._dp        ! density of new snow (kg m-3)
integer, parameter :: m = 1                             ! exponent for fsnow calc. Suggested to be 1 for global applications

!local variables
real(dp), dimension(ns:nl) :: Hi     !excess heat before phase change calculation (W m-2)
real(dp), dimension(ns:nl) :: Hin    !excess heat after phase change calculation (W m-2)
real(dp), dimension(ns:nl) :: Wliqn  !quantity of liquid water in the soil at time t+1, after phase change calc. (mm)
real(dp), dimension(ns:nl) :: Wicen  !quantity of ice in the soil at time t+1, after phase change calc. (mm)
real(dp), dimension(ns:nl) :: Wtot   !total mass of liquid and ice in the soil (mm)
real(dp), dimension(ns:nl) :: thaw   !quantity of water thawed, negative means frozen (mm)
real(dp), dimension(ns:nl) :: Fnhti  !heat flux across soil layer boundary at time t+1 (W m-2)
real(dp), dimension(ns:nl) :: Bi     !temporary variable for phase change calculations
real(dp) :: Ephase    !total latent heat of phase change W m-2
real(dp) :: dTemp     !change in temperature of the soil layer between timesteps
real(dp) :: frac
real(dp) :: tmp1
real(dp) :: Tliqmax
real(dp), dimension(nl) :: Wliqmax   !maximum liquid water fraction when the soil temp is below freezing   (mm)
integer :: l            !counters

!Point pointers
snl   => sv(j)%snl
Tsoil => sv(j)%Tsoil
Wice  => sv(j)%Wice
Wliq  => sv(j)%Wliq
Wsno  => sv(j)%Wsno
zsno  => sv(j)%zsno
Tsat  => sv(j)%Tsat
Bexp  => sv(j)%Bexp
Psat  => sv(j)%Psat
zpos  => sv(j)%zpos
Tnsoi => sv(j)%Tnsoi
dz    => sv(j)%dz
dzmm  => sv(j)%dzmm
Tice  => sv(j)%Tice
Tpor  => sv(j)%Tpor
Tliq  => sv(j)%Tliq
fice0 => sv(j)%fice0
fsnow => surf%fsnow

!----------------------------
!calculate soil energy excess or deficit vs. water phase change

!fn1 in CLM (6.24)
do l = snl+1,nl-1
    Fnhti(l) =  Kl(l) * (Tnsoi(l+1) - Tnsoi(l)) / (zpos(l+1) - zpos(l))   !CLM 6.24
end do
                
Fnhti(nl) = 0._dp  !no heat flux in bottom layer

!Brr in CLM  (middle two terms of CLM 6.42)
   !top snow/soil layer
    Bi(snl+1) = CNfac * Fhti(snl+1) + oneminusCNfac * Fnhti(snl+1)

    !other layers
do l = snl+2,nl
    Bi(l) = CNfac * (Fhti(l) - Fhti(l-1)) + oneminusCNfac * (Fnhti(l) - Fnhti(l-1))
end do

        !Frozen soil method by Niu & Yang, 2006, J Hydrometeorology (7) 937-952.
        !coded by Joe Melton Dec 14 2007

do l = 1,nl !soil layers only
 if (Tnsoi(l) < Tfreeze) then

  !(Niu & Yang (eqn 3)

    Tliqmax = Tsat(l) * (1.e3 * Lf * (Tnsoi(l) - Tfreeze) / (grav * Psat(l) * Tnsoi(l)))**(-1._dp / Bexp(l))

    if (Tliqmax > Tsat(l)) Tliqmax = Tsat(l)

    Wliqmax(l) = Tliqmax * dzmm(l)

 end if
end do

!--

Wliqn(:) = Wliq(:)
Wicen(:) = Wice(:)

Ephase   = 0._dp

!--

do l = snl+1,nl  !FLAG, not sure about whether this should be only gnl or not...

    ithaw(l) = 0
    Hi(l)    = 0._dp
    thaw(l)  = 0._dp
    Wtot(l)  = Wice(l) + Wliq(l)

           if (Wtot(l) > 0._dp) then
             fice0(l) = Wice(l) / Wtot(l)
           else  !no water/ice so skip
             fice0(l) = 0._dp
             cycle
           end if

  !assess conditions
          if (l >= 1) then ! soil layers

           if (Tnsoi(l) > Tfreeze .and. Wice(l) > 0._dp) then       !thawing conditions
               ithaw(l) = 1
               Tnsoi(l) = Tfreeze
             else if (Tnsoi(l) < Tfreeze .and. Wliq(l) > Wliqmax(l)) then  !freezing conditions
               ithaw(l) = 2
               Tnsoi(l) = Tfreeze
             else                                                 !neither freeze nor thaw
               ithaw(l) = 0
           end if

           if (snl+1 == 1 .and. Wsno > 0._dp .and. l == 1) then     !conditions where there is snow but no explict snow layer
             if (Tsoil(l) > Tfreeze) then
               ithaw = 1  !thaw
               Tnsoi(l) = Tfreeze
             end if
           end if

         else !snow layers
           if (Tnsoi(l) > Tfreeze .and. Wice(l) > 0._dp) then       !thawing conditions
               ithaw(l) = 1
               Tnsoi(l) = Tfreeze
             else if (Tnsoi(l) < Tfreeze .and. Wliq(l) > 0._dp) then  !freezing conditions
               ithaw(l) = 2
               Tnsoi(l) = Tfreeze
             else                                                 !neither freeze nor thaw
               ithaw(l) = 0
           end if

           if (snl+1 == 1 .and. Wsno > 0._dp .and. l == 1) then     !conditions where there is snow but no explict snow layer
             if (Tsoil(l) > Tfreeze) then
               ithaw = 1
               Tnsoi(l) = Tfreeze
             end if
           end if

       
         end if

  !calculate amount of water frozen or thawed
           if (ithaw(l) /= 0) then   !if there is freeze or thaw

                 dTemp = Tnsoi(l) - Tsoil(l)

                 !calculate heat excess (deficit)   (CLM 6.42)
                 if (l > snl+1) then
                  !all layers except surface
                   Hi(l) = Bi(l) - dTemp / fact(l)
                 else
                   !surface layer
                   Hi(l) = hs + dhsdT * dTemp + Bi(l) - dTemp / fact(l)
                 end if

           end if

           !if there are thawing (freezing) conditions but the heat excess (deficit) is negative (positive),
           !then there is no thaw (freeze)

           !Notes: Hi(l) = CLM hm, thaw(l) = CLM xm, xmf = Ephase

           if ((ithaw(l) == 1 .and. Hi(l) < 0._dp) .or. (ithaw(l) == 2 .and. Hi(l) > 0._dp)) then
             Hi(l) = 0._dp
             ithaw(l) = 0
           end if

           if (ithaw(l) /= 0 .and. abs(Hi(l)) > 0._dp) then


                thaw(l) = Hi(l) * dt / Lf           !amount of water thawed (frozen) kg m-2

                !----
                !if snow exists but there is no explicit layer

                if (l == 1) then
                  if (snl+1 == 1 .and. Wsno > 0._dp .and. thaw(l) > 0._dp) then

                    tmp1 = Wsno
                    Wsno = max(0._dp,tmp1 - thaw(l))
                    frac = Wsno / tmp1
                    zsno = frac * zsno
                    
                    ! fractional area of the gridcell covered by snow (Niu and Yang 2007)   
                    if (Wsno > 0._dp) then
                      fsnow = tanh(zsno / (2.5_dp * z0mg * ((min(Wsno/zsno, 800._dp)) / psno_new)**m))
                    else
                      fsnow = 0._dp
                    end if  
                      
                    Hin(l) = Hi(l) - Lf * (tmp1 - Wsno) / dt

                    if (Hin(l) > 0._dp) then
                      thaw(l) = Hin(l) * dt / Lf
                      Hi(l)   = Hin(l)
                    else
                      thaw(l) = 0._dp
                      Hi(l)   = 0._dp
                    end if

                    qsnomelt = max(0._dp, (tmp1 - Wsno))/ dt  !kg m-2 s-1
                    Ephase = Lf * qsnomelt
                  end if
                end if

    !calculate additional heat excess (deficit) after thawing (freezing)
             if (l > 0) then !soil layers only

                  Hin(l) = 0._dp

                      if (thaw(l) > 0._dp) then  !thawing (CLM 6.43)

                             Wicen(l) = max(0._dp,(Wice(l) - thaw(l)))

                     Hin(l) = Hi(l) - Lf * (Wice(l) - Wicen(l)) / dt


                   else if (thaw(l) < 0._dp) then !freezing  (CLM 6.43)

                            Wicen(l) = min(Wtot(l),Wice(l) - thaw(l))

                     Hin(l) = Hi(l) - Lf * (Wice(l) - Wicen(l)) / dt

                   end if


            else  !snow!

                  Hin(l) = 0._dp

                  if (thaw(l) > 0._dp) then  !thawing

                    Wicen(l) = max(0._dp, Wice(l) - thaw(l))
                    Hin(l) = Hi(l) - Lf * (Wice(l) - Wicen(l)) / dt

                  else if (thaw(l) < 0._dp) then  !freezing

                    Wicen(l) = min(Wtot(l), Wice(l) - thaw(l))
                    Hin(l) = Hi(l) - Lf * (Wice(l) - Wicen(l)) / dt

                  end if

                end if

      Wliqn(l) = max(0._dp, Wtot(l) - Wicen(l))

        !----------------------------------------------------------------

                if (abs(Hin(l)) > 0._dp) then

                  !if there is excess (deficit) heat, update the soil temperature

                  if (l > snl+1) then
                    Tnsoi(l) = Tnsoi(l) + fact(l) * Hin(l)

                  else
                    Tnsoi(l) = Tnsoi(l) + fact(l) * Hin(l) / (1._dp - fact(l) * dhsdT)

                  end if

                  !if there are both ice and water in the soil the temperature must be 0 deg C
                   !=---NOTE: This constraint is now removed with the Niu & Yang correction implemented above! -JM Dec 14 2007
                  !if (Wliqn(l) * Wicen(l) > 0._dp) Tnsoi(l) = Tfreeze

                end if

                Ephase = Ephase + Lf * (Wice(l) - Wicen(l)) / dt

                if (l < 1 .and. ithaw(l) == 1) then
                  qsnomelt = max(0._dp, (Wice(l) - Wicen(l))) / dt  !FLAG does this get used anywhere?
                end if

           end if


end do

!-----
Wliq(:) = Wliqn(:)
Wice(:) = Wicen(:)

  !update the soil volumetric ice content, porosity and water content.
  Tice(1:nl) = min(Tsat(1:nl), Wice(1:nl) / (dz(1:nl) * pice))
  Tpor(1:nl) = Tsat(1:nl) - Tice(1:nl)
  Tpor(1:nl) = max(0._dp,Tpor(1:nl))
  Tliq(1:nl) = min(Tpor(1:nl), Wliq(1:nl) / (dz(1:nl) * pliq))

return

end subroutine soilphasechg
