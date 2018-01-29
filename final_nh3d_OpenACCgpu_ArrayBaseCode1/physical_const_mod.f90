!!===================================================================
!! Basic theromodynamics constants used for NH model

 Module  physical_const_mod 
 Use basic_mod 

 Implicit None
 Public 

  !! NH Physical Constant -------------------------------------------
  real (kind=real_kind), Parameter :: grv =  9.81D0 !gravity 
  real (kind=real_kind), Parameter :: C_p = 1004.0D0
  real (kind=real_kind), Parameter :: C_v =  717.0D0
  real (kind=real_kind), Parameter :: R_d =  287.0D0
  real (kind=real_kind), Parameter :: Cp_Cv = C_p / C_v
  real (kind=real_kind), Parameter :: Rd_Cp = R_d / C_p
  real (kind=real_kind), Parameter :: Cp_Rd = C_p / R_d
  real (kind=real_kind), Parameter :: C_0 = 27.5629410929719D0
  !real (kind=real_kind), Parameter :: C_0 = rd**gama * p0**(-rd/cv) 
  real (kind=real_kind), Parameter :: T_0 = 300.0D0
  !real (kind=real_kind), Parameter :: T_0 = 250.0D0
  real (kind=real_kind), Parameter :: N_0 = 0.01D0     !B-V frequency
  real (kind=real_kind), Parameter :: P_0 = 100000.0D0 !Surface Prs  
  real (kind=real_kind), Parameter :: k_nu =   75.0D0  !Diffusion Coeff 
  real (kind=real_kind), Parameter :: omg_earth =   7.292D-5 !Earth rotation rate 
  real (kind=real_kind), Parameter :: r_earth = 6371.229D3  !Easrth mean radius 
  !! NH Physical Constant -------------------------------------------

 End  Module  physical_const_mod
