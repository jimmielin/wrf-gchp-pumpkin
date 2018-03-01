module module_chem_utilities
   use module_domain
   use module_model_constants
   use module_state_description
   use module_configure

contains
   subroutine chem_prep(config_flags,
                        u, v, p, pb, alt, ph,
                        phb, t, moist, n_moist,
                        rho, p_phy,
                        u_phy, v_phy, p8w, t_phy, t8w,
                        z, z_at_w, dz8w, rh,
                        fzm, fzp,
                        ids, ide, jds, jde, kds, kde, &
                        ims, ime, jms, jme, kms, kme, &
                        its, ite, jts, jte, kts, kte)
   implicit none

   type(grid_config_rec_type), intent(in) :: config_flags
   integer, intent(in) :: ids, ide, jds, jde, kds, kde, &
                          ims, ime, jms, jme, kms, kme, &
                          its, ite, jts, jte, kts, kte
   integer, intent(in) :: n_moist

   real, dimension(ims:ime, kms:kme, jms:jme, n_moist), intent(in) :: moist

   real, dimension(ims:ime, kms:kme, jms:jme), &
      intent(out) :: u_phy, v_phy, p_phy, p8w, &
                     t_phy, t8w, rho, z, dz8w, &
                     rh, z_at_w

   real, dimension(ims:ime, kms:kme, jms:jme), &
      intent(in) :: pb, &
                    p, &
                    u, &
                    v, &
                    alt, &
                    ph, &
                    phb, &
                    t

   real, dimension(kms:kme), intent(in) :: fzm, &
                                           fzp

   integer :: i_start, i_end, j_start, j_end, k_start, k_end
   integer :: i, j, k
   real    :: w1, w2, z0, z1, z2

!-----------------------------------------------------------------------
!  set up loop bounds for this grid's boundary conditions

   i_start = its
   i_end = min(ite, ide - 1)
   j_start = jts
   j_end = min(jte, jde - 1)

   k_start = kts
   k_end = min(kte, kde - 1)

!  compute thermodynamics and velocities at pressure points
   do j = j_start, j_end
   do k = k_start, k_end
   do i = i_start, i_end

      p_phy(i, k, j) = p(i, k, j) + pb(i, k, j)
      t_phy(i, k, j) = (t(i, k, j) + t0)*(p_phy(i, k, j)/p1000mb)**rcp
      rho(i, k, j) = 1./alt(i, k, j)*(1.+moist(i, k, j, P_QV))
      u_phy(i, k, j) = 0.5*(u(i, k, j) + u(i + 1, k, j))
      v_phy(i, k, j) = 0.5*(v(i, k, j) + v(i, k, j + 1))

   enddo
   enddo
   enddo

!  wig: added to make sure there is no junk in the top level even
!       though it should not be used
   do j = j_start, j_end
   do i = i_start, i_end
      p_phy(i, kte, j) = p_phy(i, k_end, j)
      t_phy(i, kte, j) = t_phy(i, k_end, j)
      rho(i, kte, j) = rho(i, k_end, j)
      u_phy(i, kte, j) = u_phy(i, k_end, j)
      v_phy(i, kte, j) = v_phy(i, k_end, j)
   enddo
   enddo

!  compute z at w points
   do j = j_start, j_end
   do k = k_start, kte
   do i = i_start, i_end
      z_at_w(i, k, j) = (phb(i, k, j) + ph(i, k, j))/g
   enddo
   enddo
   enddo

   do j = j_start, j_end
   do k = k_start, kte - 1
   do i = i_start, i_end
      dz8w(i, k, j) = z_at_w(i, k + 1, j) - z_at_w(i, k, j)
   enddo
   enddo
   enddo

   do j = j_start, j_end
   do i = i_start, i_end
      dz8w(i, kte, j) = 0.
   enddo
   enddo

!  compute z at p points (average of z at w points)
   do j = j_start, j_end
   do k = k_start, k_end
   do i = i_start, i_end
      z(i, k, j) = 0.5*(z_at_w(i, k, j) + z_at_w(i, k + 1, j))
      rh(i, k, j) = max(.1, MIN(.95, moist(i, k, j, p_qv)/ &
                                (3.80*exp(17.27*(t_phy(i, k, j) - 273.)/ &
                                          (t_phy(i, k, j) - 36.))/(.01*p_phy(i, k, j)))))
   enddo
   enddo
   enddo

!  interp t and p at w points

   do j = j_start, j_end
   do k = 2, k_end
   do i = i_start, i_end
      p8w(i, k, j) = fzm(k)*p_phy(i, k, j) + fzp(k)*p_phy(i, k - 1, j)
      t8w(i, k, j) = fzm(k)*t_phy(i, k, j) + fzp(k)*t_phy(i, k - 1, j)
   enddo
   enddo
   enddo

!  extrapolate p and t to surface and top.
!  we'll use an extrapolation in z for now

   do j = j_start, j_end
   do i = i_start, i_end

! bottom
      z0 = z_at_w(i, 1, j)
      z1 = z(i, 1, j)
      z2 = z(i, 2, j)
      w1 = (z0 - z2)/(z1 - z2)
      w2 = 1.-w1
      p8w(i, 1, j) = w1*p_phy(i, 1, j) + w2*p_phy(i, 2, j)
      t8w(i, 1, j) = w1*t_phy(i, 1, j) + w2*t_phy(i, 2, j)

! top
      z0 = z_at_w(i, kte, j)
      z1 = z(i, k_end, j)
      z2 = z(i, k_end - 1, j)
      w1 = (z0 - z2)/(z1 - z2)
      w2 = 1.-w1

!    p8w(i,kde,j) = w1*p_phy(i,kde-1,j)+w2*p_phy(i,kde-2,j)
!!!  bug fix      extrapolate ln(p) so p is positive definite
      p8w(i, kde, j) = exp(w1*log(p_phy(i, kde - 1, j)) + w2*log(p_phy(i, kde - 2, j)))
      t8w(i, kde, j) = w1*t_phy(i, kde - 1, j) + w2*t_phy(i, kde - 2, j)

   enddo
   enddo
END SUBROUTINE chem_prep
END MODULE module_chem_utilities
