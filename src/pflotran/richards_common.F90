module Richards_Common_module

  use Richards_Aux_module
  use Global_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "finclude/petscsys.h"
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscsnes.h"
#include "finclude/petscviewer.h"
#include "finclude/petsclog.h"


! Cutoff parameters
  PetscReal, parameter :: eps       = 1.D-8
  PetscReal, parameter :: floweps   = 1.D-24
  PetscReal, parameter :: perturbation_tolerance = 1.d-6
  
  public RichardsAccumulation, &
         RichardsAccumDerivative, &
         RichardsFlux, &
         RichardsFluxDerivative, &
         RichardsBCFlux, &
         RichardsBCFluxDerivative
  
contains

! ************************************************************************** !

subroutine RichardsAccumDerivative(rich_auxvar,global_auxvar, &
                                   material_auxvar, &
                                   option,sat_func,J)
  ! 
  ! Computes derivatives of the accumulation
  ! term for the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/13/07
  ! 

  use Option_module
  use Saturation_Function_module
  use Material_module, only : MaterialCompressSoil
  use Material_Aux_class
  
  implicit none

  type(richards_auxvar_type) :: rich_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  type(saturation_function_type) :: sat_func
  PetscReal :: J(option%nflowdof,option%nflowdof)
     
  PetscInt :: ispec 
  PetscReal :: vol_over_dt
  PetscReal :: por
  PetscReal :: tempreal
  PetscReal :: derivative
  PetscReal :: compressed_porosity, dcompressed_porosity_dp

  PetscInt :: iphase, ideriv
  type(richards_auxvar_type) :: rich_auxvar_pert
  type(global_auxvar_type) :: global_auxvar_pert
  ! leave as type
  type(material_auxvar_type) :: material_auxvar_pert
  PetscReal :: x(1), x_pert(1), pert, res(1), res_pert(1), J_pert(1,1)

  vol_over_dt = material_auxvar%volume/option%flow_dt
  por = material_auxvar%porosity
  
  derivative = 0.d0
  if (soil_compressibility_index > 0) then
    tempreal = global_auxvar%sat(1)*global_auxvar%den(1)
    call MaterialCompressSoil(material_auxvar,global_auxvar%pres(1), &
                                 compressed_porosity,dcompressed_porosity_dp)
    por = compressed_porosity
    derivative = derivative + dcompressed_porosity_dp * tempreal
  endif

  ! accumulation term units = dkmol/dp
  derivative = derivative + (global_auxvar%sat(1)*rich_auxvar%dden_dp + &
                             rich_auxvar%dsat_dp*global_auxvar%den(1)) * por

  J(1,1) = derivative * vol_over_dt

  if (option%numerical_derivatives_flow) then
    call GlobalAuxVarInit(global_auxvar_pert,option)  
    call MaterialAuxVarInit(material_auxvar_pert,option)  
    call RichardsAuxVarCopy(rich_auxvar,rich_auxvar_pert,option)
    call GlobalAuxVarCopy(global_auxvar,global_auxvar_pert,option)
    call MaterialAuxVarCopy(material_auxvar,material_auxvar_pert,option)
    x(1) = global_auxvar%pres(1)
    call RichardsAccumulation(rich_auxvar,global_auxvar,material_auxvar, &
                              option,res)
    ideriv = 1
    pert = max(dabs(x(ideriv)*perturbation_tolerance),0.1d0)
    x_pert = x
    if (x_pert(ideriv) < option%reference_pressure) pert = -1.d0*pert
    x_pert(ideriv) = x_pert(ideriv) + pert
    
    call RichardsAuxVarCompute(x_pert(1),rich_auxvar_pert,global_auxvar_pert, &
                               material_auxvar_pert,sat_func,option)
#if 0      
      select case(ideriv)
        case(1)
!         print *, 'dvis_dp:', auxvar%dvis_dp, (auxvar_pert%vis-auxvar%vis)/pert(ideriv)
!         print *, 'dkr_dp:', auxvar%dkr_dp, (auxvar_pert%kr-auxvar%kr)/pert(ideriv)
          print *, 'dsat_dp:', auxvar%dsat_dp, (global_auxvar_pert%sat-global_auxvar%sat)/pert
          print *, 'dden_dp:', auxvar%dden_dp, (global_auxvar_pert%den-global_auxvar%den)/pert
!          print *, 'dkvr_dp:', auxvar%dkvr_dp, (rich_auxvar_pert%kvr-rich_auxvar%kvr)/pert
          print *, 'dkvr_x_dp:', auxvar%dkvr_x_dp, (rich_auxvar_pert%kvr_x-rich_auxvar%kvr_x)/pert
          print *, 'dkvr_y_dp:', auxvar%dkvr_y_dp, (rich_auxvar_pert%kvr_y-rich_auxvar%kvr_y)/pert
          print *, 'dkvr_z_dp:', auxvar%dkvr_z_dp, (rich_auxvar_pert%kvr_z-rich_auxvar%kvr_z)/pert
      end select     
#endif     
    call RichardsAccumulation(rich_auxvar_pert,global_auxvar_pert, &
                              material_auxvar_pert, &
                              option,res_pert)
    J_pert(1,1) = (res_pert(1)-res(1))/pert
    J = J_pert
    call GlobalAuxVarStrip(global_auxvar_pert)  
    call MaterialAuxVarStrip(material_auxvar_pert)  
  endif
   
end subroutine RichardsAccumDerivative

! ************************************************************************** !

subroutine RichardsAccumulation(rich_auxvar,global_auxvar, &
                                material_auxvar, &
                                option,Res)
  ! 
  ! Computes the non-fixed portion of the accumulation
  ! term for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/13/07
  ! 

  use Option_module
  use Material_module
  use Material_Aux_class
  
  implicit none

  type(richards_auxvar_type) :: rich_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  PetscReal :: Res(1:option%nflowdof)
  
  PetscReal :: por, por1, vol_over_dt
  PetscReal :: compressed_porosity, dcompressed_porosity_dp
       
  vol_over_dt = material_auxvar%volume/option%flow_dt
  por = material_auxvar%porosity
  
  if (soil_compressibility_index > 0) then
    call MaterialCompressSoil(material_auxvar,global_auxvar%pres(1), &
                              compressed_porosity,dcompressed_porosity_dp)
    por = compressed_porosity
  endif
    
  ! accumulation term units = kmol/s
  Res(1) = global_auxvar%sat(1) * global_auxvar%den(1) * por * &
           vol_over_dt
  
end subroutine RichardsAccumulation

! ************************************************************************** !

subroutine RichardsFluxDerivative(rich_auxvar_up,global_auxvar_up, &
                                  material_auxvar_up,sir_up, & 
                                  rich_auxvar_dn,global_auxvar_dn, &
                                  material_auxvar_dn,sir_dn, &
                                  area, dist, &
                                  option,sat_func_up,sat_func_dn,Jup,Jdn)
  ! 
  ! Computes the derivatives of the internal flux terms
  ! for the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/13/07
  ! 
  use Option_module 
  use Saturation_Function_module 
  use Material_Aux_class
  use Connection_module
  
  implicit none
  
  type(richards_auxvar_type) :: rich_auxvar_up, rich_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_up, sir_dn
  PetscReal :: v_darcy, area, dist(-1:3)
  type(saturation_function_type) :: sat_func_up, sat_func_dn
  PetscReal :: Jup(option%nflowdof,option%nflowdof)
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)
     
  PetscReal :: q
  PetscReal :: ukvr,Dq
  PetscReal :: dd_up, dd_dn, perm_up, perm_dn
  PetscReal :: dist_gravity
  PetscReal :: upweight,density_ave,cond,gravity,dphi
  
  PetscReal :: dden_ave_dp_up, dden_ave_dp_dn
  PetscReal :: dgravity_dden_up, dgravity_dden_dn
  PetscReal :: dphi_dp_up, dphi_dp_dn
  PetscReal :: dukvr_dp_up, dukvr_dp_dn
!  PetscReal :: dukvr_x_dp_up, dukvr_x_dp_dn
!  PetscReal :: dukvr_y_dp_up, dukvr_y_dp_dn
!  PetscReal :: dukvr_z_dp_up, dukvr_z_dp_dn
  PetscReal :: dq_dp_up, dq_dp_dn

  PetscInt :: iphase, ideriv
  type(richards_auxvar_type) :: rich_auxvar_pert_up, rich_auxvar_pert_dn
  type(global_auxvar_type) :: global_auxvar_pert_up, global_auxvar_pert_dn
  ! leave as type
  type(material_auxvar_type) :: material_auxvar_pert_up, material_auxvar_pert_dn
  PetscReal :: x_up(1), x_dn(1), x_pert_up(1), x_pert_dn(1), pert_up, pert_dn, &
            res(1), res_pert_up(1), res_pert_dn(1), J_pert_up(1,1), J_pert_dn(1,1)
  
  v_darcy = 0.D0  
  ukvr = 0.d0
  
  Jup = 0.d0
  Jdn = 0.d0 
  
  dden_ave_dp_up = 0.d0
  dden_ave_dp_dn = 0.d0
  dgravity_dden_up = 0.d0
  dgravity_dden_dn = 0.d0
  dphi_dp_up = 0.d0
  dphi_dp_dn = 0.d0
  dukvr_dp_up = 0.d0
  dukvr_dp_dn = 0.d0
  dq_dp_up = 0.d0
  dq_dp_dn = 0.d0
  
  call ConnectionCalculateDistances(dist,option%gravity,dd_up,dd_dn, &
                                    dist_gravity,upweight)
  call material_auxvar_up%PermeabilityTensorToScalar(dist,perm_up)
  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)
  
  Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)
  
! Flow term
  if (global_auxvar_up%sat(1) > sir_up .or. global_auxvar_dn%sat(1) > sir_dn) then
    if (global_auxvar_up%sat(1) <eps) then 
      upweight=0.d0
    else if (global_auxvar_dn%sat(1) <eps) then 
      upweight=1.d0
    endif
    density_ave = upweight*global_auxvar_up%den(1)+ &
                  (1.D0-upweight)*global_auxvar_dn%den(1)
    dden_ave_dp_up = upweight*rich_auxvar_up%dden_dp
    dden_ave_dp_dn = (1.D0-upweight)*rich_auxvar_dn%dden_dp

    gravity = (upweight*global_auxvar_up%den(1) + &
              (1.D0-upweight)*global_auxvar_dn%den(1)) &
              * FMWH2O * dist_gravity
    dgravity_dden_up = upweight*FMWH2O*dist_gravity
    dgravity_dden_dn = (1.d0-upweight)*FMWH2O*dist_gravity

    dphi = global_auxvar_up%pres(1) - global_auxvar_dn%pres(1)  + gravity
    dphi_dp_up = 1.d0 + dgravity_dden_up*rich_auxvar_up%dden_dp
    dphi_dp_dn = -1.d0 + dgravity_dden_dn*rich_auxvar_dn%dden_dp

    if (dphi>=0.D0) then
#ifdef USE_ANISOTROPIC_MOBILITY
      if (dabs(dist(1))==1) then
        ukvr = rich_auxvar_up%kvr_x
        dukvr_dp_up = rich_auxvar_up%dkvr_x_dp
      else if (dabs(dist(2))==1) then
        ukvr = rich_auxvar_up%kvr_y
        dukvr_dp_up = rich_auxvar_up%dkvr_y_dp
      else if (dabs(dist(3))==1) then
        ukvr = rich_auxvar_up%kvr_z
        dukvr_dp_up = rich_auxvar_up%dkvr_z_dp
      end if
#else
      ukvr = rich_auxvar_up%kvr
      dukvr_dp_up = rich_auxvar_up%dkvr_dp
#endif
    else
#ifdef USE_ANISOTROPIC_MOBILITY    
      if (dabs(dist(1))==1) then
        ukvr = rich_auxvar_dn%kvr_x
        dukvr_dp_dn = rich_auxvar_dn%dkvr_x_dp
      else if (dabs(dist(2))==1) then
        ukvr = rich_auxvar_dn%kvr_y
        dukvr_dp_dn = rich_auxvar_dn%dkvr_y_dp
      else if (dabs(dist(3))==1) then
        ukvr = rich_auxvar_dn%kvr_z
        dukvr_dp_dn = rich_auxvar_dn%dkvr_z_dp
      end if
#else
      ukvr = rich_auxvar_dn%kvr
      dukvr_dp_dn = rich_auxvar_dn%dkvr_dp
#endif
    endif      
   
    if (ukvr>floweps) then
      v_darcy= Dq * ukvr * dphi
   
      q = v_darcy * area
      dq_dp_up = Dq*(dukvr_dp_up*dphi+ukvr*dphi_dp_up)*area
      dq_dp_dn = Dq*(dukvr_dp_dn*dphi+ukvr*dphi_dp_dn)*area
      
      Jup(1,1) = (dq_dp_up*density_ave+q*dden_ave_dp_up)
      Jdn(1,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)
    endif
  endif

 ! note: Res is the flux contribution, for node up J = J + Jup
 !                                              dn J = J - Jdn  

  if (option%numerical_derivatives_flow) then
    call GlobalAuxVarInit(global_auxvar_pert_up,option)
    call GlobalAuxVarInit(global_auxvar_pert_dn,option)  
    call MaterialAuxVarInit(material_auxvar_pert_up,option)
    call MaterialAuxVarInit(material_auxvar_pert_dn,option)  
    call RichardsAuxVarCopy(rich_auxvar_up,rich_auxvar_pert_up,option)
    call RichardsAuxVarCopy(rich_auxvar_dn,rich_auxvar_pert_dn,option)
    call GlobalAuxVarCopy(global_auxvar_up,global_auxvar_pert_up,option)
    call GlobalAuxVarCopy(global_auxvar_dn,global_auxvar_pert_dn,option)
    call MaterialAuxVarCopy(material_auxvar_up,material_auxvar_pert_up,option)
    call MaterialAuxVarCopy(material_auxvar_dn,material_auxvar_pert_dn,option)
    x_up(1) = global_auxvar_up%pres(1)
    x_dn(1) = global_auxvar_dn%pres(1)
    call RichardsFlux(rich_auxvar_up,global_auxvar_up,material_auxvar_up, &
                      sir_up, &
                      rich_auxvar_dn,global_auxvar_dn,material_auxvar_dn, &
                      sir_dn, &
                      area, dist, &
                      option,v_darcy,res)
    ideriv = 1
!    pert_up = x_up(ideriv)*perturbation_tolerance
    pert_up = max(dabs(x_up(ideriv)*perturbation_tolerance),0.1d0)
    if (x_up(ideriv) < option%reference_pressure) pert_up = -1.d0*pert_up
!    pert_dn = x_dn(ideriv)*perturbation_tolerance
    pert_dn = max(dabs(x_dn(ideriv)*perturbation_tolerance),0.1d0)
    if (x_dn(ideriv) < option%reference_pressure) pert_dn = -1.d0*pert_dn
    x_pert_up = x_up
    x_pert_dn = x_dn
    x_pert_up(ideriv) = x_pert_up(ideriv) + pert_up
    x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
    call RichardsAuxVarCompute(x_pert_up(1),rich_auxvar_pert_up, &
                               global_auxvar_pert_up, &
                               material_auxvar_pert_up,sat_func_up, &
                               option)
    call RichardsAuxVarCompute(x_pert_dn(1),rich_auxvar_pert_dn, &
                               global_auxvar_pert_dn, &
                               material_auxvar_pert_dn,sat_func_dn, &
                               option)
    call RichardsFlux(rich_auxvar_pert_up,global_auxvar_pert_up, &
                      material_auxvar_pert_up,sir_up, &
                      rich_auxvar_dn,global_auxvar_dn, &
                      material_auxvar_dn,sir_dn, &
                      area, dist, &
                      option,v_darcy,res_pert_up)
    call RichardsFlux(rich_auxvar_up,global_auxvar_up, &
                      material_auxvar_up,sir_up, &
                      rich_auxvar_pert_dn,global_auxvar_pert_dn, &
                      material_auxvar_pert_dn,sir_dn, &
                      area, dist, &
                      option,v_darcy,res_pert_dn)
    J_pert_up(1,ideriv) = (res_pert_up(1)-res(1))/pert_up
    J_pert_dn(1,ideriv) = (res_pert_dn(1)-res(1))/pert_dn
    Jup = J_pert_up
    Jdn = J_pert_dn
    call GlobalAuxVarStrip(global_auxvar_pert_up)
    call GlobalAuxVarStrip(global_auxvar_pert_dn)    
    call MaterialAuxVarStrip(material_auxvar_pert_up)
    call MaterialAuxVarStrip(material_auxvar_pert_dn)    
  endif

end subroutine RichardsFluxDerivative

! ************************************************************************** !

subroutine RichardsFlux(rich_auxvar_up,global_auxvar_up, &
                        material_auxvar_up,sir_up, &
                        rich_auxvar_dn,global_auxvar_dn, &
                        material_auxvar_dn,sir_dn, &
                        area, dist, &
                        option,v_darcy,Res)
  ! 
  ! Computes the internal flux terms for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/13/07
  ! 
  use Option_module
  use Material_Aux_class
  use Connection_module
  
  implicit none
  
  type(richards_auxvar_type) :: rich_auxvar_up, rich_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_up, sir_dn
  PetscReal :: v_darcy, area, dist(-1:3)
  PetscReal :: Res(1:option%nflowdof) 
     
  PetscInt :: ispec
  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: dd_up, dd_dn, perm_up, perm_dn
  PetscReal :: fluxm, q
  PetscReal :: ukvr,Dq
  PetscReal :: upweight, density_ave, cond, gravity, dphi
  
  fluxm = 0.d0
  v_darcy = 0.D0  
  ukvr = 0.d0
  
  call ConnectionCalculateDistances(dist,option%gravity,dd_up,dd_dn, &
                                    dist_gravity,upweight)
  call material_auxvar_up%PermeabilityTensorToScalar(dist,perm_up)
  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)

  Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)
  
! Flow term
  if (global_auxvar_up%sat(1) > sir_up .or. global_auxvar_dn%sat(1) > sir_dn) then
    if (global_auxvar_up%sat(1) <eps) then 
      upweight=0.d0
    else if (global_auxvar_dn%sat(1) <eps) then 
      upweight=1.d0
    endif
    

    density_ave = upweight*global_auxvar_up%den(1)+ &
                  (1.D0-upweight)*global_auxvar_dn%den(1)

    gravity = (upweight*global_auxvar_up%den(1) + &
              (1.D0-upweight)*global_auxvar_dn%den(1)) &
              * FMWH2O * dist_gravity

    dphi = global_auxvar_up%pres(1) - global_auxvar_dn%pres(1)  + gravity

    if (dphi>=0.D0) then
#ifdef USE_ANISOTROPIC_MOBILITY       
      if (dabs(dabs(dist(1))-1) < 1e-6) then
        ukvr = rich_auxvar_up%kvr_x
      else if (dabs(dabs(norm(2))-1) < 1e-6) then
        ukvr = rich_auxvar_up%kvr_y
      else if (dabs(dabs(norm(3))-1) < 1e-6) then
        ukvr = rich_auxvar_up%kvr_z
      end if
#else
      ukvr = rich_auxvar_up%kvr
#endif
    else
#ifdef USE_ANISOTROPIC_MOBILITY       
      if (dabs(dabs(norm(1))-1) < 1e-6) then
        ukvr = rich_auxvar_dn%kvr_x
      else if (dabs(dabs(norm(2))-1) < 1e-6) then
        ukvr = rich_auxvar_dn%kvr_y
      else if (dabs(dabs(norm(3))-1) < 1e-6) then
        ukvr = rich_auxvar_dn%kvr_z
      end if
#else
      ukvr = rich_auxvar_dn%kvr
#endif
    endif      

    if (ukvr>floweps) then
      v_darcy= Dq * ukvr * dphi
   
      q = v_darcy * area

      fluxm = q*density_ave       
    endif
  endif 

  Res(1) = fluxm
 ! note: Res is the flux contribution, for node 1 R = R + Res_FL
 !                                              2 R = R - Res_FL  

end subroutine RichardsFlux

! ************************************************************************** !

subroutine RichardsBCFluxDerivative(ibndtype,auxvars, &
                                    rich_auxvar_up,global_auxvar_up, &
                                    rich_auxvar_dn,global_auxvar_dn, &
                                    material_auxvar_dn, &
                                    sir_dn, &
                                    area,dist,option, &
                                    sat_func_dn,Jdn)
  ! 
  ! Computes the derivatives of the boundary flux
  ! terms for the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/13/07
  ! 
  use Option_module
  use Saturation_Function_module
  use Material_Aux_class
  use EOS_Water_module
  use Utility_module
 
  implicit none
  
  PetscInt :: ibndtype(:)
  type(richards_auxvar_type) :: rich_auxvar_up, rich_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_dn
  PetscReal :: auxvars(:) ! from aux_real_var array in boundary condition
  PetscReal :: area
  ! dist(-1) = fraction_upwind
  ! dist(0) = magnitude
  ! dist(1:3) = unit vector
  ! dist(0)*dist(1:3) = vector
  PetscReal :: dist(-1:3)
  type(saturation_function_type) :: sat_func_dn  
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)
  
  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: perm_dn
  
  PetscReal :: v_darcy
  PetscReal :: q,density_ave
  PetscReal :: ukvr,diffdp,Dq
  PetscReal :: upweight,cond,gravity,dphi

  PetscReal :: dden_ave_dp_dn
  PetscReal :: dgravity_dden_dn
  PetscReal :: dphi_dp_dn
  PetscReal :: dukvr_dp_dn
  PetscReal :: dq_dp_dn
  PetscInt :: pressure_bc_type
  PetscReal :: dphi_x_dp_dn,dphi_y_dp_dn,dphi_z_dp_dn
  PetscReal :: dHdn_x_dp_dn, dHdn_y_dp_dn, dHdn_z_dp_dn

  PetscInt :: iphase, ideriv
  type(richards_auxvar_type) :: rich_auxvar_pert_dn, rich_auxvar_pert_up
  type(global_auxvar_type) :: global_auxvar_pert_dn, global_auxvar_pert_up
  class(material_auxvar_type), allocatable :: material_auxvar_pert_dn, &
                                              material_auxvar_pert_up
  PetscReal :: perturbation
  PetscReal :: x_dn(1), x_up(1), x_pert_dn(1), x_pert_up(1), pert_dn, res(1), &
            res_pert_dn(1), J_pert_dn(1,1)
  PetscReal :: rho
  PetscReal :: dum1
  PetscReal :: dq_lin,dP_lin
  PetscReal :: q_approx, dq_approx
  PetscErrorCode :: ierr

  v_darcy = 0.d0
  ukvr = 0.d0
  density_ave = 0.d0
  q = 0.d0

  Jdn = 0.d0 
  
  dden_ave_dp_dn = 0.d0
  dgravity_dden_dn = 0.d0
  dphi_dp_dn = 0.d0
  dukvr_dp_dn = 0.d0
  dq_dp_dn = 0.d0

  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)
  
  ! Flow
  pressure_bc_type = ibndtype(RICHARDS_PRESSURE_DOF)
  select case(pressure_bc_type)
    ! figure out the direction of flow
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC,HET_SURF_SEEPAGE_BC)

      ! dist(0) = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3) = vector(3) - unit vector
      dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))

      if (pressure_bc_type == CONDUCTANCE_BC) then
        Dq = auxvars(RICHARDS_CONDUCTANCE_DOF)
      else
        Dq = perm_dn / dist(0)
      endif
      ! Flow term
      if (global_auxvar_up%sat(1) > sir_dn .or. global_auxvar_dn%sat(1) > sir_dn) then
        upweight=1.D0
        if (global_auxvar_up%sat(1) < eps) then 
          upweight=0.d0
        else if (global_auxvar_dn%sat(1) < eps) then 
          upweight=1.d0
        endif
        
        density_ave = upweight*global_auxvar_up%den(1) + (1.D0-upweight)*global_auxvar_dn%den(1)
        dden_ave_dp_dn = (1.D0-upweight)*rich_auxvar_dn%dden_dp

        gravity = (upweight*global_auxvar_up%den(1) + &
                  (1.D0-upweight)*global_auxvar_dn%den(1)) &
                  * FMWH2O * dist_gravity
        dgravity_dden_dn = (1.d0-upweight)*FMWH2O*dist_gravity

        dphi = global_auxvar_up%pres(1) - global_auxvar_dn%pres(1) + gravity
        dphi_dp_dn = -1.d0 + dgravity_dden_dn*rich_auxvar_dn%dden_dp

        if (pressure_bc_type == SEEPAGE_BC .or. &
            pressure_bc_type == CONDUCTANCE_BC .or. &
            pressure_bc_type == HET_SURF_SEEPAGE_BC) then
              ! flow in         ! boundary cell is <= pref
          if (dphi > 0.d0 .and. global_auxvar_up%pres(1)-option%reference_pressure < eps) then
            dphi = 0.d0
            dphi_dp_dn = 0.d0
          endif
        endif
        
        if (dphi>=0.D0) then
#ifdef USE_ANISOTROPIC_MOBILITY  
          if (dabs(dabs(dist(1))-1) < 1e-6) then
            ukvr = rich_auxvar_up%kvr_x
          else if (dabs(dabs(dist(2))-1) < 1e-6) then
            ukvr = rich_auxvar_up%kvr_y
          else if (dabs(dabs(dist(3))-1) < 1e-6) then
            ukvr = rich_auxvar_up%kvr_z
          end if
#else
          ukvr = rich_auxvar_up%kvr
#endif
        else
#ifdef USE_ANISOTROPIC_MOBILITY
          if (dabs(dabs(dist(1))-1) < 1e-6) then
            ukvr = rich_auxvar_dn%kvr_x
            dukvr_dp_dn = rich_auxvar_dn%dkvr_x_dp
          else if (dabs(dabs(dist(2))-1) < 1e-6) then
            ukvr = rich_auxvar_dn%kvr_y
            dukvr_dp_dn = rich_auxvar_dn%dkvr_y_dp
          else if (dabs(dabs(dist(3))-1) < 1e-6) then
            ukvr = rich_auxvar_dn%kvr_z
            dukvr_dp_dn = rich_auxvar_dn%dkvr_z_dp
          end if
#else
          ukvr = rich_auxvar_dn%kvr
          dukvr_dp_dn = rich_auxvar_dn%dkvr_dp
#endif
        endif      
     
        if (ukvr*Dq>floweps) then
          v_darcy = Dq * ukvr * dphi
          q = v_darcy * area
          dq_dp_dn = Dq*(dukvr_dp_dn*dphi + ukvr*dphi_dp_dn)*area

          ! If running with surface-flow model, ensure (darcy_velocity*dt) does
          ! not exceed depth of standing water.
          if (.not. rich_auxvar_dn%bcflux_default_scheme) then
            if (pressure_bc_type == HET_SURF_SEEPAGE_BC .and. option%surf_flow_on) then
              call EOSWaterdensity(option%reference_temperature, &
                                   option%reference_pressure,rho,dum1,ierr)

              if (global_auxvar_dn%pres(1) <= rich_auxvar_dn%P_min) then
              
                ! Linear approximation
                call Interpolate(rich_auxvar_dn%range_for_linear_approx(2), &
                                 rich_auxvar_dn%range_for_linear_approx(1), &
                                 global_auxvar_dn%pres(1), &
                                 rich_auxvar_dn%range_for_linear_approx(4), &
                                 rich_auxvar_dn%range_for_linear_approx(3), &
                                 q_approx)
                v_darcy = q_approx/area
                q       = q_approx

                dP_lin = rich_auxvar_dn%range_for_linear_approx(2) - &
                         rich_auxvar_dn%range_for_linear_approx(1)
                dq_lin = rich_auxvar_dn%range_for_linear_approx(4) - &
                         rich_auxvar_dn%range_for_linear_approx(3)
                dq_dp_dn = dq_lin/dP_lin

              else
                if (global_auxvar_dn%pres(1) <= rich_auxvar_dn%P_max) then

                  ! Cubic approximation
                  call CubicPolynomialEvaluate(rich_auxvar_dn%coeff_for_cubic_approx, &
                                               global_auxvar_dn%pres(1) - option%reference_pressure, &
                                               q_approx, dq_approx)
                  v_darcy = q_approx/area
                  q = q_approx
                  dq_dp_dn = dq_approx
                endif
              endif
            endif
          endif

        endif

      endif

    case(NEUMANN_BC)
      if (dabs(auxvars(RICHARDS_PRESSURE_DOF)) > floweps) then
        v_darcy = auxvars(RICHARDS_PRESSURE_DOF)
        if (v_darcy > 0.d0) then 
          density_ave = global_auxvar_up%den(1)
        else 
          density_ave = global_auxvar_dn%den(1)
          dden_ave_dp_dn = rich_auxvar_dn%dden_dp
        endif 
        q = v_darcy * area
      endif

    case(UNIT_GRADIENT_BC)

      if (global_auxvar_dn%sat(1) > sir_dn) then

        dphi = dot_product(option%gravity,dist(1:3))*global_auxvar_dn%den(1)*FMWH2O      
        density_ave = global_auxvar_dn%den(1)

        dphi_dp_dn = dot_product(option%gravity,dist(1:3))*rich_auxvar_dn%dden_dp*FMWH2O
        dden_ave_dp_dn = rich_auxvar_dn%dden_dp

        ! since boundary auxvar is meaningless (no pressure specified there), only use cell
#ifdef USE_ANISOTROPIC_MOBILITY
        if (dabs(dabs(dist(1))-1) < 1e-6) then
          ukvr = rich_auxvar_dn%kvr_x
          dukvr_dp_dn = rich_auxvar_dn%dkvr_x_dp
        else if (dabs(dabs(dist(2))-1) < 1e-6) then
          ukvr = rich_auxvar_dn%kvr_y
          dukvr_dp_dn = rich_auxvar_dn%dkvr_y_dp
        else if (dabs(dabs(dist(3))-1) < 1e-6) then
          ukvr = rich_auxvar_dn%kvr_z
          dukvr_dp_dn = rich_auxvar_dn%dkvr_z_dp
        end if
#else
        ukvr = rich_auxvar_dn%kvr
        dukvr_dp_dn = rich_auxvar_dn%dkvr_dp
#endif
     
        if (ukvr*perm_dn>floweps) then
          v_darcy = perm_dn * ukvr * dphi
          q = v_darcy*area
          dq_dp_dn = perm_dn*(dukvr_dp_dn*dphi+ukvr*dphi_dp_dn)*area
        endif 
     
      endif

  end select

  Jdn(1,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)

  if (option%numerical_derivatives_flow) then
    call GlobalAuxVarInit(global_auxvar_pert_up,option)
    call GlobalAuxVarInit(global_auxvar_pert_dn,option)  
    allocate(material_auxvar_pert_up,material_auxvar_pert_dn)
    call MaterialAuxVarInit(material_auxvar_pert_up,option)  
    call MaterialAuxVarInit(material_auxvar_pert_dn,option)  
    call RichardsAuxVarCopy(rich_auxvar_up,rich_auxvar_pert_up,option)
    call RichardsAuxVarCopy(rich_auxvar_dn,rich_auxvar_pert_dn,option)
    call GlobalAuxVarCopy(global_auxvar_up,global_auxvar_pert_up,option)
    call GlobalAuxVarCopy(global_auxvar_dn,global_auxvar_pert_dn,option)
    call MaterialAuxVarCopy(material_auxvar_dn,material_auxvar_pert_up, &
                            option)
    call MaterialAuxVarCopy(material_auxvar_dn,material_auxvar_pert_dn, &
                            option)
    x_up(1) = global_auxvar_up%pres(1)
    x_dn(1) = global_auxvar_dn%pres(1)
    ideriv = 1
    if (ibndtype(ideriv) == ZERO_GRADIENT_BC) then
      x_up(ideriv) = x_dn(ideriv)
    endif
    call RichardsBCFlux(ibndtype,auxvars, &
                        rich_auxvar_up,global_auxvar_up, &
                        rich_auxvar_dn,global_auxvar_dn, &
                        material_auxvar_dn, &
                        sir_dn, &
                        area,dist,option,v_darcy,res)
    if (pressure_bc_type == ZERO_GRADIENT_BC) then
      x_pert_up = x_up
    endif
    ideriv = 1
!    pert_dn = x_dn(ideriv)*perturbation_tolerance    
    pert_dn = max(dabs(x_dn(ideriv)*perturbation_tolerance),0.1d0)
    if (x_dn(ideriv) < option%reference_pressure) pert_dn = -1.d0*pert_dn
    x_pert_dn = x_dn
    x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
    x_pert_up = x_up
    if (ibndtype(ideriv) == ZERO_GRADIENT_BC) then
      x_pert_up(ideriv) = x_pert_dn(ideriv)
    endif   
    call RichardsAuxVarCompute(x_pert_dn(1),rich_auxvar_pert_dn, &
                               global_auxvar_pert_dn, &
                               material_auxvar_pert_dn,sat_func_dn, &
                               option)
    call RichardsAuxVarCompute(x_pert_up(1),rich_auxvar_pert_up, &
                               global_auxvar_pert_up, &
                               material_auxvar_pert_up,sat_func_dn, &
                               option)
    call RichardsBCFlux(ibndtype,auxvars, &
                        rich_auxvar_pert_up,global_auxvar_pert_up, &
                        rich_auxvar_pert_dn,global_auxvar_pert_dn, &
                        material_auxvar_pert_dn, &
                        sir_dn, &
                        area,dist,option,v_darcy,res_pert_dn)
    J_pert_dn(1,ideriv) = (res_pert_dn(1)-res(1))/pert_dn
    Jdn = J_pert_dn
    call GlobalAuxVarStrip(global_auxvar_pert_up)
    call GlobalAuxVarStrip(global_auxvar_pert_dn)   
    call MaterialAuxVarStrip(material_auxvar_pert_up)
    call MaterialAuxVarStrip(material_auxvar_pert_dn)
    deallocate(material_auxvar_pert_up,material_auxvar_pert_dn)
  endif

end subroutine RichardsBCFluxDerivative

! ************************************************************************** !

subroutine RichardsBCFlux(ibndtype,auxvars, &
                          rich_auxvar_up, global_auxvar_up, &
                          rich_auxvar_dn, global_auxvar_dn, &
                          material_auxvar_dn, &
                          sir_dn, &
                          area, dist, option,v_darcy,Res)
  ! 
  ! Computes the  boundary flux terms for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/13/07
  ! 
  use Option_module
  use Material_Aux_class
  use EOS_Water_module
  use Utility_module
 
  implicit none
  
  PetscInt :: ibndtype(:)
  type(richards_auxvar_type) :: rich_auxvar_up, rich_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_dn
  PetscReal :: auxvars(:) ! from aux_real_var array
  PetscReal :: v_darcy, area
  ! dist(-1) = fraction_upwind - not applicable here
  ! dist(0) = magnitude
  ! dist(1:3) = unit vector
  ! dist(0)*dist(1:3) = vector
  PetscReal :: dist(-1:3)
  PetscReal :: Res(1:option%nflowdof) 
  
  PetscReal :: dist_gravity  ! distance along gravity vector
          
  PetscInt :: ispec
  PetscReal :: perm_dn
  PetscReal :: fluxm,q,density_ave
  PetscReal :: ukvr,diffdp,Dq
  PetscReal :: upweight,cond,gravity,dphi
  PetscInt :: pressure_bc_type
  PetscReal :: dphi_x,dphi_y,dphi_z
  PetscReal :: rho
  PetscReal :: dum1
  PetscReal :: q_approx, dq_approx
  PetscErrorCode :: ierr
  
  fluxm = 0.d0
  v_darcy = 0.d0
  density_ave = 0.d0
  q = 0.d0
  ukvr = 0.d0

  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)  

  ! Flow  
  pressure_bc_type = ibndtype(RICHARDS_PRESSURE_DOF)
  select case(pressure_bc_type)
    ! figure out the direction of flow
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC,HET_SURF_SEEPAGE_BC)

      ! dist(0) = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3) = vector(3) - unit vector
      dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))

      if (pressure_bc_type == CONDUCTANCE_BC) then
        Dq = auxvars(RICHARDS_CONDUCTANCE_DOF)
      else
        Dq = perm_dn / dist(0)
      endif
      ! Flow term
      if (global_auxvar_up%sat(1) > sir_dn .or. global_auxvar_dn%sat(1) > sir_dn) then
        upweight=1.D0
        if (global_auxvar_up%sat(1) < eps) then 
          upweight=0.d0
        else if (global_auxvar_dn%sat(1) < eps) then 
          upweight=1.d0
        endif
        density_ave = upweight*global_auxvar_up%den(1)+(1.D0-upweight)*global_auxvar_dn%den(1)
   
        gravity = (upweight*global_auxvar_up%den(1) + &
                  (1.D0-upweight)*global_auxvar_dn%den(1)) &
                  * FMWH2O * dist_gravity
       
        dphi = global_auxvar_up%pres(1) - global_auxvar_dn%pres(1) + gravity

        if (pressure_bc_type == SEEPAGE_BC .or. &
            pressure_bc_type == CONDUCTANCE_BC .or. &
            pressure_bc_type == HET_SURF_SEEPAGE_BC) then
              ! flow in         ! boundary cell is <= pref
          if (dphi > 0.d0 .and. global_auxvar_up%pres(1)-option%reference_pressure < eps) then
            dphi = 0.d0
          endif
        endif
   
       if (dphi>=0.D0) then
#ifdef USE_ANISOTROPIC_MOBILITY       
         if (dabs(dabs(dist(1))-1) < 1e-6) then
           ukvr = rich_auxvar_up%kvr_x
         else if (dabs(dabs(dist(2))-1) < 1e-6) then
           ukvr = rich_auxvar_up%kvr_y
         else if (dabs(dabs(dist(3))-1) < 1e-6) then
           ukvr = rich_auxvar_up%kvr_z
         end if
#else
         ukvr = rich_auxvar_up%kvr
#endif
       else
#ifdef USE_ANISOTROPIC_MOBILITY
         if (dabs(dabs(dist(1))-1) < 1e-6) then
           ukvr = rich_auxvar_dn%kvr_x
         else if (dabs(dabs(dist(2))-1) < 1e-6) then
           ukvr = rich_auxvar_dn%kvr_y
         else if (dabs(dabs(dist(3))-1) < 1e-6) then
           ukvr = rich_auxvar_dn%kvr_z
         end if
#else
         ukvr = rich_auxvar_dn%kvr
#endif
       endif
        
       if (ukvr*Dq>floweps) then
        v_darcy = Dq * ukvr * dphi

        ! If running with surface-flow model, ensure (darcy_velocity*dt) does
        ! not exceed depth of standing water.
        if (pressure_bc_type == HET_SURF_SEEPAGE_BC .and. option%surf_flow_on) then
          call EOSWaterdensity(option%reference_temperature, &
                               option%reference_pressure,rho,dum1,ierr)

          if (.not. rich_auxvar_dn%bcflux_default_scheme) then
            if (global_auxvar_dn%pres(1) <= rich_auxvar_dn%P_min) then

              if (rich_auxvar_dn%range_for_linear_approx(1) == -99999.d0) then
                call printErrMsg(option,'Coeffs for linear approx for darcy flux not set')
              endif

              ! Linear approximation
              call Interpolate(rich_auxvar_dn%range_for_linear_approx(2), &
                               rich_auxvar_dn%range_for_linear_approx(1), &
                               global_auxvar_dn%pres(1), &
                               rich_auxvar_dn%range_for_linear_approx(4), &
                               rich_auxvar_dn%range_for_linear_approx(3), &
                               q_approx)
              v_darcy = q_approx/area

            else if (global_auxvar_dn%pres(1) <= rich_auxvar_dn%P_max) then

              if (rich_auxvar_dn%coeff_for_cubic_approx(1) == -99999.d0) then
                call printErrMsg(option,'Coeffs for cubic approx for darcy flux not set')
              endif

              ! Cubic approximation
              call CubicPolynomialEvaluate(rich_auxvar_dn%coeff_for_cubic_approx, &
                                           global_auxvar_dn%pres(1) - option%reference_pressure, &
                                           q_approx, dq_approx)
              v_darcy = q_approx/area
            endif
          endif

        endif
       endif
      endif 

    case(NEUMANN_BC)
      if (dabs(auxvars(RICHARDS_PRESSURE_DOF)) > floweps) then
        v_darcy = auxvars(RICHARDS_PRESSURE_DOF)

        if (v_darcy > 0.d0) then 
          density_ave = global_auxvar_up%den(1)
        else 
          density_ave = global_auxvar_dn%den(1)
        endif 
      endif

    case(UNIT_GRADIENT_BC)

      dphi = dot_product(option%gravity,dist(1:3))*global_auxvar_dn%den(1)*FMWH2O
      density_ave = global_auxvar_dn%den(1)

      ! since boundary auxvar is meaningless (no pressure specified there), only use cell
#ifdef USE_ANISOTROPIC_MOBILITY
      if (dabs(dabs(dist(1))-1) < 1e-6) then
        ukvr = rich_auxvar_dn%kvr_x
      else if (dabs(dabs(dist(2))-1) < 1e-6) then
        ukvr = rich_auxvar_dn%kvr_y
      else if (dabs(dabs(dist(3))-1) < 1e-6) then
        ukvr = rich_auxvar_dn%kvr_z
      end if
#else
      ukvr = rich_auxvar_dn%kvr
#endif
     
      if (ukvr*perm_dn>floweps) then
        v_darcy = perm_dn * ukvr * dphi
      endif 

  end select

  q = v_darcy * area

  fluxm = q*density_ave

  Res(1)=fluxm

end subroutine RichardsBCFlux

end module Richards_Common_module
