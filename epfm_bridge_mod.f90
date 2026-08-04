module epfm_bridge_mod
   use shr_kind_mod,       only : r8
   use shr_param_mod,      only : NPARAM, sa_params
   use shr_param_mod,      only : Param_Feta, Param_Hscale, Param_Dscale
   use shr_ctrl_mod,       only : WATER_LAYER
   use math_utilities_mod, only : Randn
   use data_buffer_mod,    only : m_waterTemp, m_mixTopIndex, m_Zw
   use thermal_mod,        only : EnforceThermalConsistency 
   use epfm_da_mod,        only : DefaultObsOperator, epfm_cfg 
   use epfm_da_mod,        only : PerturbPositiveParameter 

   implicit none
   private
   public :: ObsOperator
   public :: RerunWindow
   public :: PerturbInitialTemperatureProfile
   public :: perturbinitialsensparameters 

contains
   subroutine PerturbInitialTemperatureProfile(partId)
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
         noise(ii) = Randn()
      end do

      ! Gaussian vertical smoothing.
      do kk = 1, WATER_LAYER+1, 1
         smooth_noise(kk) = 0.0_r8
         weight_sum = 0.0_r8

         do jj = 1, WATER_LAYER+1, 1
            weight = exp(-0.5_r8 * ((m_Zw(jj) - m_Zw(kk)) / &
               epfm_cfg%corr_len)**2)

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
         sigma_k = epfm_cfg%sigma_deep + (epfm_cfg%sigma_lwst - &
            epfm_cfg%sigma_deep) * exp( -depth/max(epfm_cfg%corr_len, &
            1.0e-6_r8))

         m_waterTemp(kk) = m_waterTemp(kk) + sigma_k * smooth_noise(kk)
      end do

      ! Enforce broad physical bounds.
      m_waterTemp = min(max(m_waterTemp, epfm_cfg%temp_min), epfm_cfg%temp_max)
      call EnforceThermalConsistency()
   end subroutine

   subroutine PerturbInitialSensParameters(partId)
      integer, intent(in) :: partId
      ! local variables
      real(r8) :: base_value

      ! keep one particle exactly equal to the original ALBM configuration
      if (partId==1) then
         return
      end if 
      
      base_value = sa_params(Param_Feta)
      sa_params(Param_Feta) = PerturbPositiveParameter(base_value, &
         epfm_cfg%feta_log_sd, epfm_cfg%feta_bounds(1), &
         epfm_cfg%feta_bounds(2))

      base_value = sa_params(Param_Hscale)
      sa_params(Param_Hscale) = PerturbPositiveParameter(base_value, &
         epfm_cfg%hscale_log_sd, epfm_cfg%hscale_bounds(1), &
         epfm_cfg%hscale_bounds(2))

      base_value = sa_params(Param_Dscale)
      sa_params(Param_Dscale) = PerturbPositiveParameter(base_value, &
         epfm_cfg%dscale_log_sd, epfm_cfg%dscale_bounds(1), &
         epfm_cfg%dscale_bounds(2))
   end subroutine

   subroutine ObsOperator(tw, params, sim_lwst, sim_secchi)
      real(r8), intent(in)  :: tw(:)
      real(r8), intent(in)  :: params(:)
      real(r8), intent(out) :: sim_lwst
      real(r8), intent(out) :: sim_secchi

      ! Replace this with your native optics/Secchi diagnostic if available.
      call DefaultObsOperator(tw, params, sim_lwst, sim_secchi)
   end subroutine

   subroutine RerunWindow(tw_prev, params, tw_new, sim_lwst, sim_secchi)
      real(r8), intent(in)  :: tw_prev(:)
      real(r8), intent(in)  :: params(:)
      real(r8), intent(out) :: tw_new(:)
      real(r8), intent(out) :: sim_lwst
      real(r8), intent(out) :: sim_secchi

      real(r8) :: params_save(NPARAM)

      !------------------------------------------------------------
      ! This is the only ALBM-specific rerun hook you need to fill in.
      !
      ! Required behavior:
      !   1) restore previous analysis state tw_prev into model state
      !   2) set sa_params(Param_Feta/Wstr/Hscale) = proposed values
      !   3) rerun ALBM from previous analysis time to current obs time
      !      using the same forcing window already used for the forecast
      !   4) return updated tw_new and simulated observation equivalents
      !
      ! You likely already have a one-day integrator or a DA forecast
      ! wrapper in your DA branch. Put that call here.
      !------------------------------------------------------------

      params_save = sa_params

      sa_params(Param_Feta)   = params(Param_Feta)
      sa_params(Param_Hscale) = params(Param_Hscale)
      sa_params(Param_Dscale) = params(Param_Dscale)

      ! Example skeleton:
      !
      ! call RestorePreviousAnalysisState(tw_prev)
      ! call RunModelOneAssimilationWindow()
      ! tw_new    = m_waterTemp(1:nz)
      ! sim_lwst  = m_waterTemp(1)
      ! sim_secchi = DiagnoseSecchiFromALBM()
      !
      ! Temporary fallback:
      tw_new     = tw_prev
      sim_lwst   = tw_new(1)
      sim_secchi = 1.7_r8 / max(1.0e-6_r8, 0.5_r8 * sa_params(Param_Feta))

      sa_params = params_save
   end subroutine

end module epfm_bridge_mod
