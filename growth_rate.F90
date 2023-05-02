module growth_rate_module

  use Cosmological_Parameters

  implicit none
  private

  public :: growth_rate_at_z

contains
    
  double precision function growth_rate(omega_m_z, omega_lambda_z)

    implicit none
    real, intent(in) :: omega_m_z, omega_lambda_z

    growth_rate = omega_m_z**(4./7.) + (1.0+omega_m_z/2.)*omega_lambda_z/70.0

    return
  end function growth_rate


  double precision function growth_rate_at_z(z)
    !
    ! Compute the growth rate at redshift z given omega_m and omega_lambda
    ! at redshift zero (taken from the cosmology module)
    !
    implicit none
    real, intent(in) :: z
    double precision :: a
    real :: omega_m_z, omega_lambda_z, omega_k_z
    double precision :: h_of_a

    a = 1.0/(1.0+z)
    h_of_a = sqrt(omega0*a**(-3) + (1.-omega0-lambda0)*a**(-2) + lambda0)
    omega_m_z = omega0 * a**(-3) / h_of_a**2
    omega_lambda_z = 1.0 - omega_m_z
    omega_k_z = 0.0 ! Assume flat cosmology
    
    growth_rate_at_z = growth_rate(omega_m_z, omega_lambda_z)

    return
  end function growth_rate_at_z

end module growth_rate_module
