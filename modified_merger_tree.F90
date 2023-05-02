module Modified_Merger_Tree
  implicit none
  save
  ! Parameters used to modify the merger rate used in split.F90
  ! to be slightly different to the standard Press-Schechter formula
  !
  !
  ! Best fit to Millennium simulation: gamma_1=0.38, gamma_2=-0.01, gamma_3=0.0, G0=0.57
  ! Original Galform behaviour       : gamma_1=0.00, gamma_2= 0.00, gamma_3=0.0, G0=1.00
  ! Andrew Benson's 2020 fit         : gamma_1=-0.158, gamma_2=0.0488, gamma_3=0.202, G0=0.943 (from AJB email 13-Nov-2020)
  !
  ! The modify factor is
  ! G0 [sigma(m1)/sigma(m2)]^gamma_1 [w/sigma(m2)]^gamma_2 [1-sigma^2(m2)/sigma^2(m1)]^gamma_3
  !
  !real, parameter :: gamma_1=0.02,gamma_2=0.1,G0=0.82
  real :: gamma_1,gamma_2,gamma_3,gamma_4,gamma_5,G0=0.0
end module Modified_Merger_Tree
