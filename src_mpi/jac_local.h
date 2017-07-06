//subroutine loop_Jacobian(r_n, r_ca, r_c, Jac)
void loop_Jacobian(double r_n[3][3], double r_ca[3][3], double r_c[3][3], double *Jac);
double det3(double j1[3], double j2[3], double j3[3]);
void loop_Jacobian(double r_n[3][3], double r_ca[3][3], double r_c[3][3], double *Jac)
 {
//  real(dp), intent(in) :: r_n(3,3), r_ca(3,3), r_c(3,3) 
//  real(dp), intent(out) :: Jac
//  integer :: n, i, c1, c2, m
  int n, i, c1, c2, m;
//  real(dp) :: axis(3,6), pivot(3,6), r_ca3(3), r_cac3(3), j(6,6), det
  double axis[6][3], pivot[6][3], r_ca3[3], r_cac3[3], j[6][6], det;
//  real(dp) :: ran_axis(3), ran_q(4), ran_U(3,3)
  double ran_axis[3], ran_q[4], ran_U[3][3];
//  real(dp) :: rn(3,3), rca(3,3), rc(3,3), x, y, z, ran_angle 
  double rn[3][3], rca[3][3], rc[3][3], x, y, z, ran_angle;
  double arg_a[3], arg_b[3];
  int k;
  double va_1[3], va_2[3], va_3[3];
  double vb_1[3], vb_2[3], vb_3[3];
  double vc_1[3], vc_2[3], vc_3[3];
  double vd_1[3], vd_2[3], vd_3[3];
  double pro_rn[3], pro_rca[3], pro_rc[3];
  double normal_value;

//  m = 0
  m = 0;
  for(i=0;i<6;i++)
    for(k=0;k<6;k++)
      j[i][k] = 0.;
                                                                                                                           
//  rn(:,:) = r_n(:,:)
//  rca(:,:) = r_ca(:,:)
//  rc(:,:) = r_c(:,:)
  for(i=0;i<3;i++)
    for(k=0;k<3;k++)
     {
      rn[i][k] = r_n[i][k];
      rca[i][k] = r_ca[i][k];
      rc[i][k] = r_c[i][k];
     }

//111 continue
here:
//  m = m + 1
  m++;

//  ! rotation axis and pivot for phi/psi
//  n = 0
  n = 0;
//  do i = 1, 3
  for(i=0;i<3;i++)
   {
//     ! phi
//     n = n + 1
    n++;
//     axis(:,n) = rca(:,i) - rn(:,i)
    axis[n-1][0] = rca[i][0] - rn[i][0];
    axis[n-1][1] = rca[i][1] - rn[i][1];
    axis[n-1][2] = rca[i][2] - rn[i][2];
    normal_value = sqrt(dot_product(axis[n-1], axis[n-1]));
//     axis(:,n) = axis(:,n)/sqrt(dot_product(axis(:,n), axis(:,n)))
    axis[n-1][0] = axis[n-1][0]/normal_value;
    axis[n-1][1] = axis[n-1][1]/normal_value;
    axis[n-1][2] = axis[n-1][2]/normal_value;
//     pivot(:,n) = rca(:,i)
    pivot[n-1][0] = rca[i][0];
    pivot[n-1][1] = rca[i][1];
    pivot[n-1][2] = rca[i][2];
//     ! psi
//     n = n + 1
    n++;
//     axis(:,n) = rc(:,i) - rca(:,i)
    axis[n-1][0] = rc[i][0] - rca[i][0];
    axis[n-1][1] = rc[i][1] - rca[i][1];
    axis[n-1][2] = rc[i][2] - rca[i][2];
    normal_value = sqrt(dot_product(axis[n-1], axis[n-1]));
//     axis(:,n) = axis(:,n)/sqrt(dot_product(axis(:,n), axis(:,n)))
    axis[n-1][0] = axis[n-1][0]/normal_value;
    axis[n-1][1] = axis[n-1][1]/normal_value;
    axis[n-1][2] = axis[n-1][2]/normal_value;
//     pivot(:,n) = rc(:,i)
    pivot[n-1][0] = rc[i][0];
    pivot[n-1][1] = rc[i][1];
    pivot[n-1][2] = rc[i][2];
//  end do
   }

//  r_ca3(:) = rca(:,3)
  r_ca3[0] = rca[2][0];
  r_ca3[1] = rca[2][1];
  r_ca3[2] = rca[2][2];
//  r_cac3(:) = axis(:,6)
  r_cac3[0] = axis[5][0];
  r_cac3[1] = axis[5][1];
  r_cac3[2] = axis[5][2];

//  if (abs(r_cac3(3)) < 1.0d-10) then 
  if (fabs(r_cac3[2]) < 1.0e-10)
   {
//     ! x and y components of cross(axis(:,n), r_cac3(:)) become dependent.
//     ! so take x, and z components instead.
//     c1 = 4
    c1 = 3;
//     c2 = 6
    c2 = 5;
   }
//  else 
  else
   {
//     c1 = 4
    c1 = 3;
//     c2 = 5
    c2 = 4;
   }
//  end if

//  ! Jacobian elements (j11 - j14), (j21 - j24), (j31 - j34) 
//  do n = 1, 4
//     call cross_local(axis(:,n), r_ca3(:) - pivot(:,n), j(1:3,n))
//  end do
  for(n=0;n<4;n++)
   {
    arg_a[0] = r_ca3[0] - pivot[n][0];
    arg_a[1] = r_ca3[1] - pivot[n][1];
    arg_a[2] = r_ca3[2] - pivot[n][2];
    cross(axis[n], arg_a, arg_b);
    j[n][0] = arg_b[0];
    j[n][1] = arg_b[1];
    j[n][2] = arg_b[2];
   }
//  ! Jacobian elements (j41 - j45), (j51 - j55)
//  do n = 1, 5
//     call cross_local(axis(:,n), r_cac3(:), j(4:6,n))
//  end do
  for(n=0;n<5;n++)
   {
    cross(axis[n], r_cac3, arg_b);
    j[n][3] = arg_b[0];
    j[n][4] = arg_b[1];
    j[n][5] = arg_b[2];
   }

//  ! compute determinant
//  det = -(j(c2,5)*j(c1,1) - j(c1,5)*j(c2,1)) * det3(j(1:3,2), j(1:3,3), j(1:3,4)) &
//       +(j(c2,5)*j(c1,2) - j(c1,5)*j(c2,2)) * det3(j(1:3,1), j(1:3,3), j(1:3,4)) &
//       -(j(c2,5)*j(c1,3) - j(c1,5)*j(c2,3)) * det3(j(1:3,1), j(1:3,2), j(1:3,4)) &
//       +(j(c2,5)*j(c1,4) - j(c1,5)*j(c2,4)) * det3(j(1:3,1), j(1:3,2), j(1:3,3)) 
  for(i=0;i<3;i++)
   {
    va_1[i] = j[1][i];
    va_2[i] = j[2][i];
    va_3[i] = j[3][i];
    vb_1[i] = j[0][i];
    vb_2[i] = j[2][i];
    vb_3[i] = j[3][i];
    vc_1[i] = j[0][i];
    vc_2[i] = j[1][i];
    vc_3[i] = j[3][i];
    vd_1[i] = j[0][i];
    vd_2[i] = j[1][i];
    vd_3[i] = j[2][i];
   }
  det = -(j[4][c2]*j[0][c1] - j[4][c1]*j[0][c2]) * det3(va_1, va_2, va_3)
        +(j[4][c2]*j[1][c1] - j[4][c1]*j[1][c2]) * det3(vb_1, vb_2, vb_3)  
        -(j[4][c2]*j[2][c1] - j[4][c1]*j[2][c2]) * det3(vc_1, vc_2, vc_3) 
        +(j[4][c2]*j[3][c1] - j[4][c1]*j[3][c2]) * det3(vd_1, vd_2, vd_3);

//  if (abs(det) < 1.0d-10) then
  if (fabs(det) < 1.0e-10)
   {
//!     print*, 'randomly rotating:'
//     ! rotate randomly : rn(3,3), rca(3,3), rc(3,3)
//     ! random rotation axis:
//     x = random()
    x = drand48();
//     y = sqrt(1.0d0 - x*x)*random()
    y = sqrt(1.0e0 - x*x)*drand48();
//     z = sqrt(1.0d0 - x*x - y*y)
    z = sqrt(1.0e0 - x*x - y*y);
//     ran_axis(1:3) = (/ x, y, z /)
    ran_axis[0] = x;
    ran_axis[1] = y;
    ran_axis[2] = z;
//     ran_angle = two_pi*random()
    ran_angle = 2.*PI*drand48();
//     call quaternion_local(ran_axis, 0.25d0*ran_angle, ran_q)
    quaternion(ran_axis, 0.25*ran_angle, ran_q);
//     call rotation_matrix_local(ran_q, ran_U)
    rotation_matrix(ran_q, ran_U);
//     do i = 1, 3
//        rn(:,i) = matmul(ran_U, rn(:,i))
//        rca(:,i) = matmul(ran_U, rca(:,i))
//        rc(:,i) = matmul(ran_U, rc(:,i))
//     end do
   for(i=0;i<3;i++)
    {
     matmul(ran_U, rn[i], pro_rn);
     matmul(ran_U, rca[i], pro_rca);
     matmul(ran_U, rc[i], pro_rc);
     for(k=0;k<3;k++)
      {
       rn[i][k] = pro_rn[k];
       rca[i][k] = pro_rca[k];
       rc[i][k] = pro_rc[k];
      }
    }

//     if (m < 2) goto 111
   if(m<2)
     goto here;

//     print*, 'det', det
//     print*, rc(:,3)
//     print*, rca(:,3)

//     print*, -(j(5,5)*j(4,1) - j(4,5)*j(5,1)) * det3(j(1:3,2), j(1:3,3), j(1:3,4))
//     print*, (j(5,5)*j(4,2) - j(4,5)*j(5,2)) * det3(j(1:3,1), j(1:3,3), j(1:3,4))
//     print*,-(j(5,5)*j(4,3) - j(4,5)*j(5,3)) * det3(j(1:3,1), j(1:3,2), j(1:3,4))
//     print*,(j(5,5)*j(4,4) - j(4,5)*j(5,4)) * det3(j(1:3,1), j(1:3,2), j(1:3,3))
//     stop
    *Jac = -1.;
    return;
//    printf("ERROR: Jacobian is ZERO.\n");
//    exit(1);   
   }
//  end if
//  Jac = 1.0d0 / abs(det) 
  *Jac = 1.0e0 / fabs(det);

//end subroutine loop_Jacobian
  return;
 }
//!-----------------------------------------------------------------------
//function det3(j1, j2, j3)
double det3(double j1[3], double j2[3], double j3[3])
 {
//  ! determinant of a 3X3 matrix
//  implicit none
//  real(dp), intent(in) :: j1(3), j2(3), j3(3)
//  real(dp) :: det3
  double deter;
//  real(dp) :: j11, j12, j13, j21, j22, j23, j31, j32, j33
  double j11,j12,j13,j21,j22,j23,j31,j32,j33;

//  j11 = j1(1); j12 = j1(2); j13 = j1(3)
//  j21 = j2(1); j22 = j2(2); j23 = j2(3)
//  j31 = j3(1); j32 = j3(2); j33 = j3(3)
  j11 = j1[0]; j12 = j2[0]; j13 = j3[0];
  j21 = j1[1]; j22 = j2[1]; j23 = j3[1];
  j31 = j1[2]; j32 = j2[2]; j33 = j3[2];

//  det3 = j11 * (j22*j33 - j23*j32) - j12 * (j21*j33 - j23*j31) + j13 * (j21*j32 - j22*j31)
  deter = j11 * (j22*j33 - j23*j32) - j12 * (j21*j33 - j23*j31) + j13 * (j21*j32 - j22*j31);

//end function det3
  return deter;
 }
//!-----------------------------------------------------------------------


