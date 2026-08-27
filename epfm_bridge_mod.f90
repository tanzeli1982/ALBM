module epfm_bridge_mod
   use shr_kind_mod,       only : r8
   use shr_param_mod,      only : NPARAM, sa_params
   use shr_param_mod,      only : Param_Feta, Param_Hscale, Param_Dscale
   use shr_param_mod,      only : Param_TinDiff
   use shr_ctrl_mod,       only : WATER_LAYER
   use math_utilities_mod, only : gaussian_randn 
   use data_buffer_mod,    only : m_waterTemp, m_mixTopIndex, m_Zw
   use thermal_mod,        only : EnforceThermalConsistency 
   use epfm_da_mod,        only : DefaultObsOperator, EPFM_Config
   use epfm_da_mod,        only : PerturbPositiveParameter 
   use epfm_da_mod,        only : PerturbRealParameter

   implicit none
   private
   public :: ObsOperator
   public :: PerturbTemperatureProfile
   public :: EvolveParameters 

contains
   subroutine PerturbTemperatureProfile(cfg, partId)
      type(EPFM_Config), intent(in) :: cfg
      integer, intent(in) :: partId
      ! local variables
      integer :: ii, kk, jj
      real(r8) :: weight_sum, weight
      real(r8) :: sigma_k, depth
      real(r8) :: noise(WATER_LAYER+1)
      real(r8) :: smooth_noise(WATER_LAYER+1)

      ! keep one particle exactly equal to the original ALBM configuration
      if (partId==1) then
         return
      end if

      ! Independent standard-normal values
      do ii = 1, WATER_LAYER+1, 1
         noise(ii) = gaussian_randn(0._r8, 1._r8)
      end do

      ! Gaussian vertical smoothing.
      do kk = 1, WATER_LAYER+1, 1
         smooth_noise(kk) = 0.0_r8
         weight_sum = 0.0_r8

         do jj = 1, WATER_LAYER+1, 1
            weight = exp(-0.5_r8 * ((m_Zw(jj) - m_Zw(kk)) / cfg%corr_len)**2)

            smooth_noise(kk) = smooth_noise(kk) + weight * noise(jj)
            weight_sum = weight_sum + weight
         end do

         smooth_noise(kk) = smooth_noise(kk) / max(weight_sum, 1.0e-12_r8)
      end do

      ! Normalize the smoothed profile to unit RMS.
      smooth_noise = smooth_noise / &
         max(sqrt(sum(smooth_noise**2) / DBLE(WATER_LAYER+1)), 1.0e-12_r8)

      ! Apply depth-dependent uncertainty.
      do kk = 1, WATER_LAYER+1, 1
         depth = m_Zw(1) - m_Zw(kk)
         sigma_k = cfg%sigma_deep + (cfg%sigma_lwst - cfg%sigma_deep) * &
            exp( -depth/max(cfg%corr_len, 1.0e-6_r8))

         m_waterTemp(kk) = m_waterTemp(kk) + sigma_k * smooth_noise(kk)
      end do

      ! Enforce broad physical bounds.
      m_waterTemp = min(max(m_waterTemp, cfg%temp_min), cfg%temp_max)
      call EnforceThermalConsistency()
   end subroutine

   subroutine EvolveParameters(cfg, month)
      type(EPFM_Config), intent(in) :: cfg
      integer, intent(in) :: month
      ! local variables
      real(r8) :: base_value

      base_value = sa_params(Param_Feta)
      sa_params(Param_Feta) = PerturbPositiveParameter(base_value, &
         cfg%feta_log_sd, cfg%feta_bounds(1), cfg%feta_bounds(2))

      base_value = sa_params(Param_Hscale)
      sa_params(Param_Hscale) = PerturbPositiveParameter(base_value, &
         cfg%hscale_log_sd, cfg%hscale_bounds(1), cfg%hscale_bounds(2))

      base_value = sa_params(Param_Dscale)
      sa_params(Param_Dscale) = PerturbPositiveParameter(base_value, &
         cfg%dscale_log_sd, cfg%dscale_bounds(1), cfg%dscale_bounds(2))

      if (month>=5 .and. month<=9) then
         base_value = sa_params(Param_TinDiff)
         sa_params(Param_TinDiff) = PerturbRealParameter(base_value, &
            cfg%tindiff_sd, cfg%tindiff_bounds(1), cfg%tindiff_bounds(2))
      else
         sa_params(Param_TinDiff) = 0._r8
      end if
   end subroutine

   subroutine ObsOperator(tw, params, sim_lwst, sim_secchi)
      real(r8), intent(in)  :: tw(:)
      real(r8), intent(in)  :: params(:)
      real(r8), intent(out) :: sim_lwst
      real(r8), intent(out) :: sim_secchi

      ! Replace this with your native optics/Secchi diagnostic if available.
      call DefaultObsOperator(tw, params, sim_lwst, sim_secchi)
   end subroutine

end module epfm_bridge_mod
