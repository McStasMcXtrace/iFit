!!-------------------------------------------------------
!!---- Crystallographic Fortran Modules Library (CrysFML)
!!-------------------------------------------------------
!!---- The CrysFML project is distributed under LGPL. In agreement with the
!!---- Intergovernmental Convention of the ILL, this software cannot be used
!!---- in military applications.
!!----
!!---- Copyright (C) 1999-2012  Institut Laue-Langevin (ILL), Grenoble, FRANCE
!!----                          Universidad de La Laguna (ULL), Tenerife, SPAIN
!!----                          Laboratoire Leon Brillouin(LLB), Saclay, FRANCE
!!----
!!---- Authors: Juan Rodriguez-Carvajal (ILL)
!!----          Javier Gonzalez-Platas  (ULL)
!!----
!!---- Contributors: Laurent Chapon     (ILL)
!!----               Marc Janoschek     (Los Alamos National Laboratory, USA)
!!----               Oksana Zaharko     (Paul Scherrer Institute, Switzerland)
!!----               Tierry Roisnel     (CDIFX,Rennes France)
!!----               Eric Pellegrini    (ILL)
!!----
!!---- This library is free software; you can redistribute it and/or
!!---- modify it under the terms of the GNU Lesser General Public
!!---- License as published by the Free Software Foundation; either
!!---- version 3.0 of the License, or (at your option) any later version.
!!----
!!---- This library is distributed in the hope that it will be useful,
!!---- but WITHOUT ANY WARRANTY; without even the implied warranty of
!!---- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!!---- Lesser General Public License for more details.
!!----
!!---- You should have received a copy of the GNU Lesser General Public
!!---- License along with this library; if not, see <http://www.gnu.org/licenses/>.
!!----
!!----
!!---- MODULE: CFML_Math_3D
!!----   INFO: Simple mathematics general utilities for 3D Systems
!!----
!!---- HISTORY
!!----    Update: 04/03/2011
!!----
!!---- DEPENDENCIES
!!--++    Use CFML_GlobalDeps,   only: cp, sp, dp, pi, to_rad, to_deg
!!--++    Use CFML_Math_General, only: cosd, sind
!!----
!!---- VARIABLES
!!--++    EPS
!!----    ERR_Math3D
!!----    ERR_Math3D_Mess
!!----
!!---- PROCEDURES
!!----    Functions:
!!----       CROSS_PRODUCT
!!--++       CROSS_PRODUCT_CMPL_dp     [Overloaded]
!!--++       CROSS_PRODUCT_CMPL_sp     [Overloaded]
!!--++       CROSS_PRODUCT_dp          [Overloaded]
!!--++       CROSS_PRODUCT_in          [Overloaded]
!!--++       CROSS_PRODUCT_sp          [Overloaded]
!!----       DETERM_A
!!--++       DETERM_A_I                [Overloaded]
!!--++       DETERM_A_R                [Overloaded]
!!----       DETERM_V
!!--++       DETERM_V_I                [Overloaded]
!!--++       DETERM_V_R                [Overloaded]
!!----       INVERT_A
!!--++       INVERT_DP                 [Overloaded]
!!--++       INVERT_SP                 [Overloaded]
!!----       MAT_CROSS
!!--++       MAT_CROSS_CMPL_dp     [Overloaded]
!!--++       MAT_CROSS_CMPL_sp     [Overloaded]
!!--++       MAT_CROSS_dp          [Overloaded]
!!--++       MAT_CROSS_in          [Overloaded]
!!--++       MAT_CROSS_sp          [Overloaded]
!!----       POLYHEDRON_VOLUME
!!----       ROTATE_OX
!!----       ROTATE_OY
!!----       ROTATE_OZ
!!----       TENSOR_PRODUCT
!!--++       TENSOR_PRODUCT_CMPL_dp     [Overloaded]
!!--++       TENSOR_PRODUCT_CMPL_sp     [Overloaded]
!!--++       TENSOR_PRODUCT_dp          [Overloaded]
!!--++       TENSOR_PRODUCT_in          [Overloaded]
!!--++       TENSOR_PRODUCT_sp          [Overloaded]
!!----       VECLENGTH
!!----
!!----    Subroutines:
!--..
!!--..    Init Routine
!!----       INIT_ERR_MATH3D
!!----       SET_EPS
!!----       SET_EPS_DEFAULT
!--..
!!--..    Matrix and Vectors Subroutines
!!----       GET_CART_FROM_CYLIN
!!--++       GET_CART_FROM_CYLIN_DP    [Overloaded]
!!--++       GET_CART_FROM_CYLIN_SP    [Overloaded]
!!----       GET_CENTROID_COORD
!!----       GET_CYLINDR_COORD
!!--++       GET_CYLINDR_COORD_DP      [Overloaded]
!!--++       GET_CYLINDR_COORD_SP      [Overloaded]
!!----       GET_CART_FROM_SPHER
!!--++       GET_CART_FROM_SPHER_DP    [Overloaded]
!!--++       GET_CART_FROM_SPHER_SP    [Overloaded]
!!----       GET_PLANE_FROM_POINTS
!!----       GET_SPHERIC_COORD
!!--++       GET_SPHERIC_COORD_DP      [Overloaded]
!!--++       GET_SPHERIC_COORD_SP      [Overloaded]
!!----       MATRIX_DIAGEIGEN
!!----       MATRIX_INVERSE
!!----       RESOLV_SIST_1X2
!!----       RESOLV_SIST_1X3
!!----       RESOLV_SIST_2X2
!!----       RESOLV_SIST_2X3
!!----       RESOLV_SIST_3X3
!!----
!!
 Module CFML_Math_3D
    !---- Use Modules ----!
    Use CFML_GlobalDeps,   only: cp, sp, dp, pi, to_rad, to_deg
    Use CFML_Math_General, only: cosd, sind, euclidean_norm

    implicit none

    private

    !---- List of public functions ----!
    public :: Polyhedron_Volume, Rotate_OX, Rotate_OY, Rotate_OZ, Veclength

    !---- List of public overloaded procedures: functions ----!
    public :: Cross_Product, Determ_A, Determ_V, Invert_A, Mat_Cross, Tensor_Product

    !---- List of public subroutines ----!
    public :: Init_Err_Math3D, Set_Eps, Set_Eps_Default, Matrix_DiagEigen, Matrix_Inverse, &
              Resolv_Sist_1X2, Resolv_Sist_1X3, Resolv_Sist_2X2, Resolv_Sist_2X3,          &
              Resolv_Sist_3X3, Get_Plane_from_Points, Get_Centroid_Coord

    !---- List of public overloaded procedures: subroutines ----!
    public :: Get_Cart_From_Cylin, Get_Cylindr_Coord, Get_Cart_From_Spher, Get_Spheric_Coord

    !----  Make private the overloaded procedures ----!
    private :: Cross_Product_dp, Cross_Product_sp, Determ_A_I, Determ_A_R, Determ_V_I,    &
               Determ_V_R, Invert_dp, Invert_sp, Get_Cart_From_Cylin_dp,                  &
               Get_Cart_From_Cylin_sp, Get_Cylindr_Coord_dp, Get_Cylindr_Coord_sp,        &
               Get_Cart_From_Spher_dp, Get_Cart_From_Spher_sp, Get_Spheric_Coord_dp,      &
               Get_Spheric_Coord_sp, Cross_Product_cmpl_dp, Cross_Product_cmpl_sp,        &
               Mat_Cross_dp,Mat_Cross_sp,Mat_Cross_in,Mat_Cross_cmpl_dp,Mat_Cross_cmpl_sp,&
               Tensor_Product_dp,Tensor_Product_sp,Tensor_Product_in,                     &
               Tensor_Product_cmpl_dp,Tensor_Product_cmpl_sp

    !---- Definitions ----!
    !!--++
    !!--++  EPS
    !!--++     real(kind=cp), private ::  eps=0.00001_cp
    !!--++
    !!--++  (PRIVATE)
    !!--++     Epsilon value
    !!--++
    !!--++  Update: February - 2005
    !!
    real(kind=cp),  private  ::  eps=0.00001_cp

    !!----
    !!---- ERR_Math3D
    !!----    logical :: ERR_Math3D
    !!----
    !!----    Logical Variable indicating an error in CFML_Math_3D module
    !!----
    !!---- Update: February - 2005
    !!
    logical, public  :: ERR_Math3D

    !!----
    !!---- ERR_Math3D_Mess
    !!----    character(len=150) :: ERR_Math3D_Mess
    !!----
    !!----    String containing information about the last error
    !!----
    !!---- Update: February - 2005
    !!
    character(len=150), public :: ERR_Math3D_Mess

    !---- Interfaces - Overlapp ----!
    Interface  Cross_Product
       Module Procedure Cross_product_sp
       Module Procedure Cross_product_dp
       Module Procedure Cross_product_in
       Module Procedure Cross_product_cmpl_sp
       Module Procedure Cross_product_cmpl_dp
    End Interface

    Interface  Determ_A
       Module Procedure Determ_A_I
       Module Procedure Determ_A_R
    End Interface

    Interface  Determ_V
       Module Procedure Determ_V_I
       Module Procedure Determ_V_R
    End Interface

    Interface  Invert_A
       Module Procedure Invert_sp
       Module Procedure Invert_dp
    End Interface

    Interface  Get_Cart_from_Cylin
       Module Procedure Get_Cart_from_Cylin_dp
       Module Procedure Get_Cart_from_Cylin_sp
    End Interface

    Interface  Get_Cylindr_Coord
       Module Procedure Get_Cylindr_Coord_dp
       Module Procedure Get_Cylindr_Coord_sp
    End Interface

    Interface  Get_Cart_from_Spher
       Module Procedure Get_Cart_from_Spher_dp
       Module Procedure Get_Cart_from_Spher_sp
    End Interface

    Interface  Get_Spheric_Coord
       Module Procedure Get_Spheric_Coord_dp
       Module Procedure Get_Spheric_Coord_sp
    End Interface

    Interface  Mat_Cross
       Module Procedure Mat_Cross_sp
       Module Procedure Mat_Cross_dp
       Module Procedure Mat_Cross_in
       Module Procedure Mat_Cross_cmpl_sp
       Module Procedure Mat_Cross_cmpl_dp
    End Interface

    Interface  Tensor_Product
       Module Procedure Tensor_product_sp
       Module Procedure Tensor_product_dp
       Module Procedure Tensor_product_in
       Module Procedure Tensor_product_cmpl_sp
       Module Procedure Tensor_product_cmpl_dp
    End Interface

 Contains

    !!----
    !!---- Function  Cross_Product(U,V) Result(W)
    !!----    real(kind=sp/dp), dimension(3), intent( in) :: u   !  In -> Vector 1
    !!----    real(kind=sp/dp), dimension(3), intent( in) :: v   !  In -> Vector 2
    !!----    real(kind=sp/dp), dimension(3)              :: w   ! Out -> Vector 1 x vector 2
    !!----
    !!----    Calculates the cross product of vectors u and v
    !!----    Vectors, w= u x v, are given in cartesian components.
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Function  Cross_Product_cmpl_dp(U,V) Result(W)
    !!--++    complex(kind=dp/sp), dimension(3), intent( in) :: u   !  In -> Vector 1
    !!--++    complex(kind=dp/sp), dimension(3), intent( in) :: v   !  In -> Vector 2
    !!--++    complex(kind=dp/sp), dimension(3)              :: w   ! Out -> Vector 1 x vector 2
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the cross product of the complex vectors u and v
    !!--++    Vectors, w = u x v, are given in cartesian components.
    !!--++
    !!--++ Update: June - 2012
    !!
    Function Cross_Product_cmpl_dp(u,v) Result(w)
       !---- Argument ----!
       complex(kind=dp), dimension(3), intent( in) :: u,v
       complex(kind=dp), dimension(3)              :: w

       w(1)=u(2)*v(3)-u(3)*v(2)
       w(2)=u(3)*v(1)-u(1)*v(3)
       w(3)=u(1)*v(2)-u(2)*v(1)

       return
    End Function Cross_Product_cmpl_dp

    Function Cross_Product_cmpl_sp(u,v) Result(w)
       !---- Argument ----!
       complex(kind=sp), dimension(3), intent( in) :: u,v
       complex(kind=sp), dimension(3)              :: w

       w(1)=u(2)*v(3)-u(3)*v(2)
       w(2)=u(3)*v(1)-u(1)*v(3)
       w(3)=u(1)*v(2)-u(2)*v(1)

       return
    End Function Cross_Product_cmpl_sp

    !!--++
    !!--++ Function  Cross_Product_dp(U,V) Result(W)
    !!--++    real(kind=dp), dimension(3), intent( in) :: u   !  In -> Vector 1
    !!--++    real(kind=dp), dimension(3), intent( in) :: v   !  In -> Vector 2
    !!--++    real(kind=dp), dimension(3)              :: w   ! Out -> Vector 1 x vector 2
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the cross product of vectors u and v
    !!--++    Vectors, w= u x v, are given in cartesian components.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Cross_Product_dp(u,v) Result(w)
       !---- Argument ----!
       real(kind=dp), dimension(3), intent( in) :: u,v
       real(kind=dp), dimension(3)              :: w

       w(1)=u(2)*v(3)-u(3)*v(2)
       w(2)=u(3)*v(1)-u(1)*v(3)
       w(3)=u(1)*v(2)-u(2)*v(1)

       return
    End Function Cross_Product_dp

    !!--++
    !!--++ Function  Cross_Product_in(U,V) Result(W)
    !!--++    integer, dimension(3), intent( in) :: u   !  In -> Vector 1
    !!--++    integer, dimension(3), intent( in) :: v   !  In -> Vector 2
    !!--++    integer, dimension(3)              :: w   ! Out -> Vector 1 x vector 2
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the cross product of integer vectors u and v
    !!--++    In the indices are givent w.r.t the direct lattice, the cross product
    !!--++    are indices w.r.t. reciprocal lattice and viceversa.
    !!--++
    !!--++ Update: November - 2008
    !!
    Function Cross_Product_in(u,v) Result(w)
       !---- Argument ----!
       integer, dimension(3), intent( in) :: u,v
       integer, dimension(3)              :: w

       w(1)=u(2)*v(3)-u(3)*v(2)  ! i  j   k !
       w(2)=u(3)*v(1)-u(1)*v(3)  !u1  u2  u3! = (u2.v3 - u3.v2)i + (v1.u3 - u1.v3)j + (u1.v2-u2.v1)k
       w(3)=u(1)*v(2)-u(2)*v(1)  !v1  v2  v3!

       return
    End Function Cross_Product_in

    !!--++
    !!--++ Function  Cross_Product_sp(U,V) Result(W)
    !!--++    real(kind=sp), dimension(3), intent( in) :: u   !  In -> Vector 1
    !!--++    real(kind=sp), dimension(3), intent( in) :: v   !  In -> Vector 2
    !!--++    real(kind=sp), dimension(3)              :: w   ! Out -> Vector 1 x vector 2
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the cross product of vectors u and v
    !!--++    Vectors, w= u x v, are given in cartesian components.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Cross_Product_sp(u,v) Result(w)
       !---- Argument ----!
       real(kind=sp), dimension(3), intent( in) :: u,v
       real(kind=sp), dimension(3)              :: w

       w(1)=u(2)*v(3)-u(3)*v(2)  ! i  j   k !
       w(2)=u(3)*v(1)-u(1)*v(3)  !u1  u2  u3! = (u2.v3 - u3.v2)i + (v1.u3 - u1.v3)j + (u1.v2-u2.v1)k
       w(3)=u(1)*v(2)-u(2)*v(1)  !v1  v2  v3!

       return
    End Function Cross_Product_sp

    !!----
    !!---- Function Determ_A(A)
    !!----    integer/real(kind=cp), dimension(3,3), intent(in)  :: a
    !!----
    !!----    Calculates the determinant of an integer/real 3x3 matrix
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Function Determ_A_I(A)
    !!--++    integer, dimension(3,3), intent(in)  :: a
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the determinant of an integer 3x3 matrix
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Determ_A_I(A) Result(determ)
       !---- Argument ----!
       integer, dimension(3,3), intent(in) :: A
       integer                             :: determ

       determ=A(1,1)*A(2,2)*A(3,3)+A(2,1)*A(3,2)*A(1,3)+A(1,2)*A(2,3)*A(3,1) &
             -A(1,3)*A(2,2)*A(3,1)-A(1,1)*A(3,2)*A(2,3)-A(1,2)*A(2,1)*A(3,3)

       return
    End Function Determ_A_I

    !!--++
    !!--++ Function Determ_A_R(A)
    !!--++    real(kind=cp), dimension(3,3), intent(in)  :: a
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the determinant of a real 3x3 matrix
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Determ_A_R(A) Result (determ)
       !---- Argument ----!
       real(kind=cp), dimension(3,3), intent(in) :: A
       real(kind=cp)                             :: determ

       determ=A(1,1)*A(2,2)*A(3,3)+A(2,1)*A(3,2)*A(1,3)+A(1,2)*A(2,3)*A(3,1) &
             -A(1,3)*A(2,2)*A(3,1)-A(1,1)*A(3,2)*A(2,3)-A(1,2)*A(2,1)*A(3,3)

       return
    End Function Determ_A_R

    !!----
    !!---- Function  Determ_V(a,b,c)
    !!----    integer/real(kind=cp), dimension(3), intent(in) :: a,b,c
    !!----
    !!----    Calculates the determinant of the components of three vectors
    !!----
    !!----  Update: February - 2005
    !!

    !!--++
    !!--++ Function Determ_V_I(A,B,C)
    !!--++    integer, dimension(3), intent(in) :: a,b,c
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the determinant of the components of three vectors
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Determ_V_I(a,b,c) Result(det)
       !---- Arguments ----!
       integer, dimension(3), intent(in) :: a,b,c
       integer                           :: det

       !---- Local variables ----!
       integer :: i,j,k

       det = 0
       do i = 1,3
          j = i+1
          if (j == 4) j = 1
          k = 6-i-j
          det = det+a(i)*(b(j)*c(k)-b(k)*c(j))
       end do

       return
    End Function Determ_V_I

    !!--++
    !!--++ Function Determ_V_R(A,B,C)
    !!--++    real(kin=cp), dimension(3), intent(in) :: a,b,c
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the determinant of the components of three vectors
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Determ_V_R(a,b,c) Result(det)
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent(in) :: a,b,c
       real(kind=cp)                           :: det

       !---- Local variables ----!
       integer :: i,j,k

       det = 0.0
       do i = 1,3
          j = i+1
          if (j == 4) j = 1
          k = 6-i-j
          det = det+a(i)*(b(j)*c(k)-b(k)*c(j))
       end do

       return
    End Function Determ_V_R

    !!----
    !!---- Funcion Invert_A(A) Result(b)
    !!----    real(kind=sp/dp), dimension(3,3), intent(in) :: a
    !!----    real(Kind=sp/dp), dimension(3,3)             :: b
    !!----
    !!----    Calculate de inverse of a real 3x3 matrix. If the routine fails,
    !!----    then a 0.0 matrix is returned.
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Funcion Invert_Dp(A) Result(b)
    !!--++    real(kind=dp), dimension(3,3), intent(in) :: a
    !!--++    real(Kind=dp), dimension(3,3)             :: b
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate de inverse of a real 3x3 matrix
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Invert_Dp(a) Result(b)
       !---- Arguments ----!
       real(kind=dp),dimension(3,3), intent(in) :: a
       real(kind=dp),dimension(3,3)             :: b

       !---- Local variables ----!
       real(kind=dp)  :: dmat

       b(1,1) =   a(2,2)*a(3,3)-a(2,3)*a(3,2)
       b(2,1) = -(a(2,1)*a(3,3)-a(2,3)*a(3,1))
       b(3,1) =   a(2,1)*a(3,2)-a(2,2)*a(3,1)
       b(1,2) = -(a(1,2)*a(3,3)-a(1,3)*a(3,2))
       b(2,2) =   a(1,1)*a(3,3)-a(1,3)*a(3,1)
       b(3,2) = -(a(1,1)*a(3,2)-a(1,2)*a(3,1))
       b(1,3) =   a(1,2)*a(2,3)-a(1,3)*a(2,2)
       b(2,3) = -(a(1,1)*a(2,3)-a(1,3)*a(2,1))
       b(3,3) =   a(1,1)*a(2,2)-a(1,2)*a(2,1)
       dmat = a(1,1)*b(1,1)+a(1,2)*b(2,1)+a(1,3)*b(3,1) !determinant of A

       if (abs(dmat) > tiny(dmat)) then
          b= b/dmat
       else
          b=0.0_dp
       end if

       return
    End Function Invert_Dp

    !!--++
    !!--++ Funcion Invert_Sp(A) Result(b)
    !!--++    real(kind=sp), dimension(3,3), intent(in) :: a
    !!--++    real(Kind=sp), dimension(3,3)             :: b
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate de inverse of a real 3x3 matrix
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Invert_Sp(a) Result(b)
       !---- Arguments ----!
       real(kind=sp),dimension(3,3), intent(in) :: a
       real(kind=sp),dimension(3,3)             :: b

       !---- Local variables ----!
       real(kind=sp)  :: dmat

       b(1,1) =   a(2,2)*a(3,3)-a(2,3)*a(3,2)
       b(2,1) = -(a(2,1)*a(3,3)-a(2,3)*a(3,1))
       b(3,1) =   a(2,1)*a(3,2)-a(2,2)*a(3,1)
       b(1,2) = -(a(1,2)*a(3,3)-a(1,3)*a(3,2))
       b(2,2) =   a(1,1)*a(3,3)-a(1,3)*a(3,1)
       b(3,2) = -(a(1,1)*a(3,2)-a(1,2)*a(3,1))
       b(1,3) =   a(1,2)*a(2,3)-a(1,3)*a(2,2)
       b(2,3) = -(a(1,1)*a(2,3)-a(1,3)*a(2,1))
       b(3,3) =   a(1,1)*a(2,2)-a(1,2)*a(2,1)
       dmat = a(1,1)*b(1,1)+a(1,2)*b(2,1)+a(1,3)*b(3,1) !determinant of A

       if (abs(dmat) > tiny(dmat)) then
          b= b/dmat
       else
          b=0.0
       end if

       return
    End Function Invert_Sp

    !!----
    !!---- Function  Mat_Cross(U) Result(M)
    !!----    real/complex(kind=sp/dp)/integer, dimension(3), intent( in) :: u   !  In -> Vector 1
    !!----    real/complex(kind=sp/dp)/integer, dimension(3,3)            :: M   ! Out -> Matrix [u]cross
    !!----
    !!----    Calculates the matrix corresponding to the operator u x
    !!----    Antisymmetric matrix of the form:
    !!----                /  0   -u(3)  u(2)\
    !!----    M=[u]cross=|  u(3)   0   -u(1) |
    !!----                \-u(2)  u(1)   0  /
    !!----
    !!----  Updated: June - 2012
    !!
    Function Mat_Cross_cmpl_dp(u) Result(M)
       !---- Argument ----!
       complex(kind=dp), dimension(3), intent( in) :: u
       complex(kind=dp), dimension(3,3)            :: M

       M = reshape( (/  (0.0_dp,0.0_dp),   -u(3),         u(2),  &
                            u(3),   (0.0_dp,0.0_dp),     -u(1),  &
                           -u(2),           u(1),   (0.0_dp,0.0_dp)/),(/3,3/))
       return
    End Function Mat_Cross_cmpl_dp

    Function Mat_Cross_cmpl_sp(u) Result(M)
       !---- Argument ----!
       complex(kind=sp), dimension(3), intent( in) :: u
       complex(kind=sp), dimension(3,3)            :: M

       M = reshape( (/  (0.0_sp,0.0_sp),   -u(3),         u(2),  &
                            u(3),   (0.0_sp,0.0_sp),     -u(1),  &
                           -u(2),           u(1),   (0.0_sp,0.0_sp)/),(/3,3/))
       return
    End Function Mat_Cross_cmpl_sp

    Function Mat_Cross_dp(u) Result(M)
       !---- Argument ----!
       real(kind=dp), dimension(3), intent( in) :: u
       real(kind=dp), dimension(3,3)            :: M

       M = reshape( (/ 0.0_dp,   -u(3),     u(2),  &
                        u(3),    0.0_dp,   -u(1),  &
                       -u(2),     u(1),    0.0_dp/),(/3,3/))
       return
    End Function Mat_Cross_dp

    Function Mat_Cross_in(u) Result(M)
       !---- Argument ----!
       integer, dimension(3), intent( in) :: u
       integer, dimension(3,3)            :: M

       M = reshape( (/   0,    -u(3),    u(2),  &
                        u(3),    0,     -u(1),  &
                       -u(2),   u(1),     0 /),(/3,3/))
       return
    End Function Mat_Cross_in

    Function Mat_Cross_sp(u) Result(M)
       !---- Argument ----!
       real(kind=sp), dimension(3), intent( in) :: u
       real(kind=sp), dimension(3,3)            :: M

       M = reshape( (/ 0.0_sp, -u(3),    u(2),  &
                        u(3),  0.0_sp,  -u(1),  &
                       -u(2),   u(1),   0.0_sp/),(/3,3/))
       return
    End Function Mat_Cross_sp

    !!----
    !!---- Function Polyhedron_Volume(Nv,Vert,Cent) Result(vol)
    !!----    integer,                       intent(in) :: Nv       ! Vertices Number
    !!----    real(kind=cp), dimension(:,:), intent(in) :: Vert     ! Cartesian coordinates of vertices
    !!----    real(kind=cp), dimension(3),   intent(in) :: Cent     ! Cartesian coordinates a central point
    !!----
    !!---- This procedure calculate the volume of polyhedral with Nv vertices.
    !!---- It is based on volcal program of L. W. FINGER.
    !!---- Adapted by Javier Gonzalez Platas
    !!----
    !!---- Update: February - 2010
    !!
    Function Polyhedron_Volume(NV,Vert,Cent) Result(vol)
       !---- Arguments ----!
       integer,                       intent(in) :: Nv       ! Number of Vertices
       real(kind=cp), dimension(:,:), intent(in) :: Vert     ! Cartesian coordinates of atoms
       real(kind=cp), dimension(3),   intent(in) :: Cent     ! Cartesian coordinates of Central atom
       real(kind=cp)                             :: vol
       !---- Local Variables ----!
       integer                       :: i,j,k,l,i1,j1
       real(kind=cp)                 :: z,z0,area,factor
       real(kind=cp),dimension(6)    :: vxyz
       real(kind=cp),dimension(3)    :: d
       real(kind=cp),dimension(3,Nv) :: Atm_cart

       vol=0.0
       call init_err_Math3d()

       if (nv <= 3) then
          ERR_Math3D=.true.
          ERR_Math3D_Mess='The number of vertices for polyhedron volume is less than 4'
          return
       end if

       do i=1,nv
          Atm_cart(:,i)=Vert(:,i)- Cent
       end do

       do i=1,nv-2
          i1=i+1
          do j=i1,nv-1
             j1=j+1
             vxyz(1:3)=Atm_cart(:,j)-Atm_cart(:,i)
        loop:do k=j1,nv
                vxyz(4:6)=Atm_cart(:,k)-Atm_cart(:,i)
                d(1)=vxyz(2)*vxyz(6)-vxyz(5)*vxyz(3)
                d(2)=vxyz(4)*vxyz(3)-vxyz(1)*vxyz(6)
                d(3)=vxyz(1)*vxyz(5)-vxyz(4)*vxyz(2)
                area=0.5*sqrt(d(1)*d(1)+d(2)*d(2)+d(3)*d(3))
                z0=0.5*(Atm_cart(1,i)*d(1)+Atm_cart(2,i)*d(2)+Atm_cart(3,i)*d(3))/area

                ! check for and avoid plane through origin
                if (abs(z0) < 1.0e-5) cycle
                factor = 3.0
                do l=1,nv
                   if(l==i .or. l==j .or. l==k) cycle

                   ! calculate distance of point l from plane of ijk
                   z=0.5*((Atm_cart(1,i)-Atm_cart(1,l))*d(1)+ &
                          (Atm_cart(2,i)-Atm_cart(2,l))*d(2)+ &
                          (Atm_cart(3,i)-Atm_cart(3,l))*d(3))/area

                   ! z and z0 must have the same sign
                   if (z * z0 < -0.001) cycle loop
                   if (abs(z * z0) < 0.001)then
                      ! if more than 3 corners on this face, the area will be counted twice.
                      ! change factor to handle this case.
                     factor = 6.0
                   end if
                end do

                ! all points on same side,  thus ijk are face
                ! Direction Cosines Of Plane Normal
                d=d/(2.0*area)

                vol=vol+area*abs(z0)/factor

             end do loop
          end do
       end do

       return
    End Function Polyhedron_Volume

    !!----
    !!---- Function Rotate_OX(X,Angle) Result (Vec)
    !!----    real(kind=cp), dimension(3), intent(in) :: x       !  In -> Vector
    !!----    real(kind=cp),               intent(in) :: angle   !  In -> Angle (Degrees)
    !!----    real(kind=cp), dimension(3)             :: vec     ! Out -> Vector
    !!----
    !!----    X Rotation. Positive rotation is counter-clockwise
    !!----
    !!---- Update: February - 2005
    !!
    Function Rotate_OX(X,Angle) Result(vec)
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent(in) :: x
       real(kind=cp),               intent(in) :: angle
       real(kind=cp), dimension(3)             :: vec

       !---- Variables locales ----!
       real(kind=cp), dimension(3,3)           :: rot

       rot(1,1)=  1.0
       rot(2,1)=  0.0_cp
       rot(3,1)=  0.0_cp

       rot(1,2)=  0.0_cp
       rot(2,2)=  cosd(angle)
       rot(3,2)=  sind(angle)

       rot(1,3)=  0.0_cp
       rot(2,3)=  -sind(angle)
       rot(3,3)=  cosd(angle)

       vec=matmul(rot,x)

       return
    End Function Rotate_OX

    !!----
    !!---- Function Rotate_OY(Y,Angle) Result (Vec)
    !!----    real(kind=cp), dimension(3), intent(in) :: y       !  In -> Vector
    !!----    real(kind=cp),               intent(in) :: angle   !  In -> Angle (Degrees)
    !!----    real(kind=cp), dimension(3)             :: vec     ! Out -> Vector
    !!----
    !!----    Y Rotation.
    !!----
    !!---- Update: February - 2005
    !!
    Function Rotate_OY(Y,Angle) Result(vec)
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent(in) :: y
       real(kind=cp),               intent(in) :: angle     ! Angle in degrees
       real(kind=cp), dimension(3)             :: vec

       !---- Variables locales ----!
       real(kind=cp), dimension(3,3)           :: rot

       rot(1,1)=  cosd(angle)
       rot(2,1)=  0.0_cp
       rot(3,1)=  -sind(angle)

       rot(1,2)=  0.0_cp
       rot(2,2)=  1.0_cp
       rot(3,2)=  0.0_cp

       rot(1,3)= sind(angle)
       rot(2,3)= 0.0_cp
       rot(3,3)= cosd(angle)

      vec=matmul(rot,y)

       return
    End Function Rotate_OY

    !!----
    !!---- Function Rotate_OZ(Z,Angle) Result (Vec)
    !!----    real(kind=cp), dimension(3), intent(in) :: z       !  In -> Vector
    !!----    real(kind=cp),               intent(in) :: angle   !  In -> Angle (Degrees)
    !!----    real(kind=cp), dimension(3)             :: vec     ! Out -> Vector
    !!----
    !!----    Z Rotation
    !!----
    !!---- Update: February - 2005
    !!
    Function Rotate_OZ(Z,Angle) Result(vec)
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent(in) :: z
       real(kind=cp),               intent(in) :: angle
       real(kind=cp), dimension(3)             :: vec

       !---- Variables locales ----!
       real(kind=cp), dimension(3,3)           :: rot

       rot(1,1)=  cosd(angle)
       rot(2,1)=  sind(angle)
       rot(3,1)=  0.0_cp

       rot(1,2)=  -sind(angle)
       rot(2,2)=  cosd(angle)
       rot(3,2)=  0.0_cp

       rot(1,3)=  0.0_cp
       rot(2,3)=  0.0_cp
       rot(3,3)=  1.0_cp

       vec=matmul(rot,z)

       return
    End Function Rotate_OZ

    !!----
    !!---- Function  Tensor_Product(U,V) Result(W)
    !!----    complex/real(kind=sp/dp)/integer, dimension(3), intent( in) :: u   !  In -> Vector 1
    !!----    complex/real(kind=sp/dp)/integer, dimension(3), intent( in) :: v   !  In -> Vector 2
    !!----    complex/real(kind=sp/dp)/integer, dimension(3,3)            :: w   ! Out -> Tensor product Vector1 (o) Vector2
    !!----
    !!----    Calculates the tensor product of vectors u and v
    !!----
    !!---- Updated: June - 2012
    !!
    Function Tensor_Product_cmpl_dp(u,v) Result(w)
       !---- Argument ----!
       complex(kind=dp), dimension(3), intent( in) :: u,v
       complex(kind=dp), dimension(3,3)            :: w
       !
       complex(kind=dp), dimension(3,3)            :: mu,mv
       mu=0.0_dp;  mv=0.0_dp
       mu(:,1)=u
       mv(1,:)=v
       w=matmul(mu,mv)
       return
    End Function Tensor_Product_cmpl_dp

    Function Tensor_Product_cmpl_sp(u,v) Result(w)
       !---- Argument ----!
       complex(kind=sp), dimension(3), intent( in) :: u,v
       complex(kind=sp), dimension(3,3)            :: w
       !
       complex(kind=sp), dimension(3,3)            :: mu,mv
       mu=0.0_sp;  mv=0.0_sp
       mu(:,1)=u
       mv(1,:)=v
       w=matmul(mu,mv)
       return
    End Function Tensor_Product_cmpl_sp

    Function Tensor_Product_dp(u,v) Result(w)
       !---- Argument ----!
       real(kind=dp), dimension(3), intent( in) :: u,v
       real(kind=dp), dimension(3,3)            :: w
       !
       real(kind=dp), dimension(3,3)            :: mu,mv
       mu=0.0_dp;  mv=0.0_dp
       mu(:,1)=u
       mv(1,:)=v
       w=matmul(mu,mv)
       return
    End Function Tensor_Product_dp

    Function Tensor_Product_in(u,v) Result(w)
       !---- Argument ----!
       integer, dimension(3), intent( in) :: u,v
       integer, dimension(3,3)            :: w
       !
       integer, dimension(3,3)            :: mu,mv
       mu=0;  mv=0
       mu(:,1)=u
       mv(1,:)=v
       w=matmul(mu,mv)
       return
    End Function Tensor_Product_in

    Function Tensor_Product_sp(u,v) Result(w)
       !---- Argument ----!
       real(kind=sp), dimension(3), intent( in) :: u,v
       real(kind=sp), dimension(3,3)            :: w
       !
       real(kind=sp), dimension(3,3)            :: mu,mv
       mu=0.0_sp;  mv=0.0_sp
       mu(:,1)=u
       mv(1,:)=v
       w=matmul(mu,mv)
       return
    End Function Tensor_Product_sp
    !!----
    !!---- Function Veclength(A,B) Result(c)
    !!----    real(kind=cp), dimension(3,3), intent(in)  :: a
    !!----    real(kind=cp), dimension(3),   intent(in)  :: b
    !!----    real(kind=cp),                             :: c
    !!----
    !!----    Length of vector B when A is the Crystallographic
    !!----    to orthogonal matrix length=c
    !!----
    !!---- Update: February - 2005
    !!
    Function Veclength(a,b) Result(c)
       !---- Arguments ----!
       real(kind=cp), intent(in)  , dimension(3,3)       :: a
       real(kind=cp), intent(in)  , dimension(3  )       :: b
       real(kind=cp)                                     :: c

       !---- Local variables ----!
       integer                     :: i,j
       real(kind=cp), dimension(3) :: v

       v=0.0
       do i = 1,3
          do j = 1,3
             v(i) = v(i)+a(i,j)*b(j)
          end do
       end do

       c = sqrt(v(1)**2+v(2)**2+v(3)**2)

       return
    End Function Veclength

    !---------------------!
    !---- Subroutines ----!
    !---------------------!

    !!----
    !!---- Subroutine Init_Err_Math3D()
    !!----
    !!----    Initialize the errors flags in CFML_Math_3D
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Init_Err_Math3D()

       ERR_Math3D=.false.
       ERR_Math3D_Mess=" "

       return
    End Subroutine Init_Err_Math3D

    !!----
    !!---- Subroutine Set_Eps(Neweps)
    !!----    real(kind=cp), intent( in) :: neweps
    !!----
    !!----    Sets global EPS to the value "neweps"
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Set_Eps(Neweps)
       !---- Arguments ----!
       real(kind=cp), intent( in) :: neweps

       eps=neweps

       return
    End Subroutine Set_Eps

    !!----
    !!---- Subroutine Set_Eps_Default()
    !!----
    !!----    Sets global EPS to the default value: eps=0.00001
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Set_Eps_Default()

       eps=0.00001

       return
    End Subroutine Set_Eps_Default

    !!----
    !!---- Subroutine Get_Cart_from_Cylin(rho,Phi,zeta,Xo,Mode)
    !!----    real(kind=sp/dp),              intent( in)           :: rho
    !!----    real(kind=sp/dp),              intent( in)           :: phi
    !!----    real(kind=sp/dp),              intent( in)           :: zeta
    !!----    real(kind=sp/dp), dimension(3),intent(out)           :: xo
    !!----    character(len=*),              intent( in), optional :: mode
    !!----
    !!----    Determine the Cartesian coordinates from cylindrical coordinates.
    !!----    If Mode='D' the angle phi is provided in Degrees.
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine  Get_Cart_from_Cylin_dp(rho,Phi,zeta,Xo,Mode)
    !!--++    real(kind=dp),              intent( in)           ::  rho
    !!--++    real(kind=dp),              intent( in)           ::  phi
    !!--++    real(kind=dp),              intent( in)           ::  zeta
    !!--++    real(kind=dp), dimension(3),intent(out)           ::  xo
    !!--++    character(len=*),           intent( in), optional ::  mode
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the Cartesian coordinates from cylindrical coordinates.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Get_Cart_from_Cylin_dp(rho,Phi,zeta,Xo,Mode)
       !---- Arguments ----!
       real(kind=dp),              intent( in)           ::  rho
       real(kind=dp),              intent( in)           ::  phi
       real(kind=dp),              intent( in)           ::  zeta
       real(kind=dp), dimension(3),intent(out)           ::  xo
       character(len=*),           intent( in), optional ::  mode

       !---- Local Variables ----!
       real(kind=dp) :: ph

       ph=phi
       if (present(mode)) then
          if (mode(1:1) == "D" .or. mode(1:1) == "d") ph=phi*to_rad
       end if
       xo(1)=rho*cos(ph)
       xo(2)=rho*sin(ph)
       xo(3)=zeta

       return
    End Subroutine Get_Cart_from_Cylin_dp

    !!--++
    !!--++ Subroutine  Get_Cart_from_Cylin_sp(rho,Phi,zeta,Xo,Mode)
    !!--++    real(kind=sp),              intent( in)           ::  rho
    !!--++    real(kind=sp),              intent( in)           ::  phi
    !!--++    real(kind=sp),              intent( in)           ::  zeta
    !!--++    real(kind=sp), dimension(3),intent(out)           ::  xo
    !!--++    character(len=*),           intent( in), optional ::  mode
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the Cartesian coordinates from cylindrical coordinates.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Get_Cart_from_Cylin_sp(rho,Phi,zeta,Xo,Mode)
       real(kind=sp),              intent( in)           ::  rho
       real(kind=sp),              intent( in)           ::  phi
       real(kind=sp),              intent( in)           ::  zeta
       real(kind=sp), dimension(3),intent(out)           ::  xo
       character(len=*),           intent( in), optional ::  mode

       !---- Local Variables ----!
       real(kind=sp) :: ph

       ph=phi
       if (present(mode)) then
          if (mode(1:1) == "D" .or. mode(1:1) == "d") ph=phi*to_rad
       end if
       xo(1)=rho*cos(ph)
       xo(2)=rho*sin(ph)
       xo(3)=zeta

       return
    End Subroutine Get_Cart_from_Cylin_sp

    !!----
    !!---- Subroutine Get_Centroid_Coord(Cn,Atm_Cart,Centroid,Baricenter)
    !!----    integer,                       intent(in) :: Cn          ! Coordination Number
    !!----    real(kind=cp), dimension(:,:), intent(in) :: Atm_Cart    ! Cartesian coordinates of atoms
    !!----    real(kind=cp), dimension(3),   intent(out):: Centroid    ! Centroid
    !!----    real(kind=cp), dimension(3),   intent(out):: Baricenter  ! Baricenter
    !!----
    !!---- Procedure to calculate Centroid and BariCenter of a pPolyhedron according to
    !!---- Tonci Balic-Zunic (Acta Cryst. B52, 1996, 78-81; Acta Cryst. B54, 1998, 766-773)
    !!---- Centroid is here different from Baricentre and it is defined in the above reference.
    !!----
    !!---- Update: February - 2010
    !!
    Subroutine Get_Centroid_Coord(Cn,Atm_Cart,Centroid,Baricenter)
       !---- Arguments ----!
       integer,                       intent(in) :: Cn          ! Coordination Number
       real(kind=cp), dimension(:,:), intent(in) :: Atm_Cart    ! Cartesian coordinates of atoms, gathered as: (1:3,1:Cn)
       real(kind=cp), dimension(3),   intent(out):: Centroid    ! Centroid
       real(kind=cp), dimension(3),   intent(out):: Baricenter  ! Baricenter

       !---- Local variables ----!
       real(kind=cp), dimension(4)   :: plane1,plane2,plane3
       real(kind=cp), dimension(3)   :: p0,p1,p2,p3,u,v,r,t
       real(kind=cp), dimension(3,3) :: w, w1
       real(kind=cp)                 :: d,umod,vmod,rmod,d1
       real(kind=cp)                 :: sx, sy, sz, sx2, sy2, sz2, sx3, sy3, sz3
       real(kind=cp)                 :: sxy, sxz, syz, sxy2, sxz2
       real(kind=cp)                 :: sx2y, sx2z, syz2, sy2z
       integer                       :: i

       call init_err_math3d()
       centroid=0.0
       baricenter=0.0

       p1=Atm_Cart(1:3,1)
       p2=Atm_Cart(1:3,2)
       p3=Atm_Cart(1:3,3)

       select case (cn)
          case (:2)
             err_Math3D=.true.
             err_Math3D_Mess='Centroid calculation needs 3 vertives as minimum'
             return

          case (3)
             !---- Plane 1: Defined with those 3 Points ----!
             call Get_Plane_From_Points(p1, p2, p3, &
                                        plane1(1), plane1(2), plane1(3), plane1(4))
             r=plane1(1:3)
             rmod=euclidean_norm(3,r)
             if (abs(rmod) <= 0.0001) then
                err_Math3D=.true.
                err_Math3D_Mess='Imposible to define a Plane with the three given points '
                return
             end if
             r=r/rmod

             !---- Vectors ----!
             u=p2-p1
             umod=euclidean_norm(3,u)
             if (abs(umod) <= 0.0001) then
                err_Math3D=.true.
                err_Math3D_Mess='Check your points! Seems that two of them are equal'
                return
             end if

             v=p3-p1
             vmod=euclidean_norm(3,v)
             if (abs(vmod) <= 0.0001) then
                err_Math3D=.true.
                err_Math3D_Mess='Check your points! Seems that two of them are equal'
                return
             end if

             !---- Plane 2 ----!
             p0=p1+0.5*u
             u=u/umod
             plane2(1:3)=u
             plane2(4)=-( plane2(1)*p0(1)+plane2(2)*p0(2)+plane2(3)*p0(3) )

             !---- Plane 3 ----!
             p0=p1+0.5*v
             v=v/vmod
             plane3(1:3)=v
             plane3(4)=-( plane3(1)*p0(1)+plane3(2)*p0(2)+plane3(3)*p0(3) )

             !---- Centroid ----!
             w(1,1:3)=plane1(1:3)
             w(2,1:3)=plane2(1:3)
             w(3,1:3)=plane3(1:3)
             d=determ_a(w)

             if (abs(d) <= 0.0001) then
                err_Math3D=.true.
                err_Math3D_Mess='Determinant is singular to calculate Centroid point'
                return
             end if

             w(1:3,1)=(/-plane1(4),-plane2(4), -plane3(4)/)
             d1=determ_a(w)
             centroid(1)=d1/d

             w(1,1:3)=plane1(1:3)
             w(2,1:3)=plane2(1:3)
             w(3,1:3)=plane3(1:3)
             w(1:3,2)=(/-plane1(4),-plane2(4), -plane3(4)/)
             d1=determ_a(w)
             centroid(2)=d1/d

             w(1,1:3)=plane1(1:3)
             w(2,1:3)=plane2(1:3)
             w(3,1:3)=plane3(1:3)
             w(1:3,3)=(/-plane1(4),-plane2(4), -plane3(4)/)
             d1=determ_a(w)
             centroid(3)=d1/d

             sx =0.0; sy =0.0; sz =0.0
             do i=1,3
                sx=sx+Atm_Cart(1,i)
                sy=sy+Atm_Cart(2,i)
                sz=sz+Atm_Cart(3,i)
             end do

          case (4:)
             sx =0.0; sy =0.0; sz =0.0
             sx2=0.0; sy2=0.0; sz2=0.0
             sx3=0.0; sy3=0.0; sz3=0.0
             sxy=0.0; sxz=0.0; syz=0.0
             sxy2=0.0; sxz2=0.0
             sx2y=0.0; sx2z=0.0
             syz2=0.0; sy2z=0.0
             do i=1,cn
                sx=sx+Atm_Cart(1,i)
                sy=sy+Atm_Cart(2,i)
                sz=sz+Atm_Cart(3,i)

                sx2=sx2+Atm_Cart(1,i)*Atm_Cart(1,i)
                sy2=sy2+Atm_Cart(2,i)*Atm_Cart(2,i)
                sz2=sz2+Atm_Cart(3,i)*Atm_Cart(3,i)

                sx3=sx3+Atm_Cart(1,i)*Atm_Cart(1,i)*Atm_Cart(1,i)
                sy3=sy3+Atm_Cart(2,i)*Atm_Cart(2,i)*Atm_Cart(2,i)
                sz3=sz3+Atm_Cart(3,i)*Atm_Cart(3,i)*Atm_Cart(3,i)

                sxy=sxy+Atm_Cart(1,i)*Atm_Cart(2,i)
                sxz=sxz+Atm_Cart(1,i)*Atm_Cart(3,i)
                syz=syz+Atm_Cart(2,i)*Atm_Cart(3,i)

                sxy2=sxy2+Atm_Cart(1,i)*Atm_Cart(2,i)*Atm_Cart(2,i)
                sxz2=sxz2+Atm_Cart(1,i)*Atm_Cart(3,i)*Atm_Cart(3,i)

                sx2y=sx2y+Atm_Cart(2,i)*Atm_Cart(1,i)*Atm_Cart(1,i)
                sx2z=sx2z+Atm_Cart(3,i)*Atm_Cart(1,i)*Atm_Cart(1,i)

                syz2=syz2+Atm_Cart(2,i)*Atm_Cart(3,i)*Atm_Cart(3,i)
                sy2z=sy2z+Atm_Cart(3,i)*Atm_Cart(2,i)*Atm_Cart(2,i)
             end do

             w(1,1)=sx2 - (sx**2)/real(cn)
             w(1,2)=sxy - (sx*sy)/real(cn)
             w(1,3)=sxz - (sx*sz)/real(cn)
             t(1)=0.5*(sx3 + sxy2 + sxz2 - ((sx2*sx + sy2*sx + sz2*sx)/real(cn)))

             w(2,1)=sxy - (sx*sy)/real(cn)
             w(2,2)=sy2 - (sy**2)/real(cn)
             w(2,3)=syz - (sy*sz)/real(cn)
             t(2)=0.5*(sx2y + sy3 + syz2 - ((sx2*sy + sy2*sy + sz2*sy)/real(cn)))

             w(3,1)=sxz - (sx*sz)/real(cn)
             w(3,2)=syz - (sy*sz)/real(cn)
             w(3,3)=sz2 - (sz**2)/real(cn)
             t(3)=0.5*(sx2z + sy2z + sz3 - ((sx2*sz + sy2*sz + sz2*sz)/real(cn)))

             d=determ_a(w)
             if (abs(d) <= 0.0001) then
                err_Math3D=.true.
                err_Math3D_Mess='Determinant is singular to calculate Centroid point'
                return
             end if

             w1=w
             w1(:,1)=t
             d1=determ_a(w1)
             centroid(1)=d1/d

             w1=w
             w1(:,2)=t
             d1=determ_a(w1)
             centroid(2)=d1/d

             w1=w
             w1(:,3)=t
             d1=determ_a(w1)
             centroid(3)=d1/d
       end select

       baricenter=(/ sx/real(cn), sy/real(cn), sz/real(cn) /)

       return
    End Subroutine Get_Centroid_Coord

    !!----
    !!---- Subroutine Get_Cylindr_Coord(Xo,rho,Phi,zeta,Mode)
    !!----    real(kind=sp/dp), dimension(3),intent( in)           :: xo
    !!----    real(kind=sp/dp),              intent(out)           :: rho
    !!----    real(kind=sp/dp),              intent(out)           :: phi
    !!----    real(kind=sp/dp),              intent(out)           :: zeta
    !!----    character(len=*),              intent( in), optional :: mode
    !!----
    !!----    Determine the cylindrical coordinates from Cartesian coordinates.
    !!----    If Mode='D' the angle phi is provided in Degrees.
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine  Get_Cylindr_Coord_dp(Xo,rho,Phi,zeta,Mode)
    !!--++    real(kind=dp), dimension(3),intent( in)           ::  xo
    !!--++    real(kind=dp),              intent(out)           ::  rho
    !!--++    real(kind=dp),              intent(out)           ::  phi
    !!--++    real(kind=dp),              intent(out)           ::  zeta
    !!--++    character(len=*),           intent( in), optional ::  mode
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the cylindrical coordinates from Cartesian coordinates.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Get_Cylindr_Coord_dp(Xo,rho,Phi,zeta,Mode)
       !---- Arguments ----!
       real(kind=dp), dimension(3),intent( in)           ::  xo
       real(kind=dp),              intent(out)           ::  rho
       real(kind=dp),              intent(out)           ::  phi
       real(kind=dp),              intent(out)           ::  zeta
       character(len=*),           intent( in), optional ::  mode

       !---- Local Variables ----!
       integer :: j

       zeta=xo(3)
       if( abs(xo(2)) > eps .or. abs(xo(1)) > eps) then
          phi=atan2(xo(2),xo(1))
       else
          phi= 0.0_dp
       end if
       rho=0.0_dp
       do j=1,2
          rho=rho+xo(j)*xo(j)
       end do
       rho=sqrt(rho)

       if (present(mode)) then
          if (mode(1:1) == "D" .or. mode(1:1) == "d") phi=phi*to_deg
       end if

       return
    End Subroutine Get_Cylindr_Coord_dp

    !!--++
    !!--++ Subroutine  Get_Cylindr_Coord_sp(Xo,rho,Phi,zeta,Mode)
    !!--++    real(kind=sp), dimension(3),intent( in)           ::  xo
    !!--++    real(kind=sp),              intent(out)           ::  rho
    !!--++    real(kind=sp),              intent(out)           ::  phi
    !!--++    real(kind=sp),              intent(out)           ::  zeta
    !!--++    character(len=*),           intent( in), optional ::  mode
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the cylindrical coordinates from Cartesian coordinates.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Get_Cylindr_Coord_sp(Xo,rho,Phi,zeta,Mode)
       !---- Arguments ----!
       real(kind=sp), dimension(3),intent( in)           ::  xo
       real(kind=sp),              intent(out)           ::  rho
       real(kind=sp),              intent(out)           ::  phi
       real(kind=sp),              intent(out)           ::  zeta
       character(len=*),           intent( in), optional ::  mode

       !---- Local Variables ----!
       integer :: j

       zeta=xo(3)
       if( abs(xo(2)) > eps .or. abs(xo(1)) > eps) then
          phi=atan2(xo(2),xo(1))
       else
          phi= 0.0_sp
       end if
       rho=0.0_sp
       do j=1,2
          rho=rho+xo(j)*xo(j)
       end do
       rho=sqrt(rho)

       if (present(mode)) then
          if (mode(1:1) == "D" .or. mode(1:1) == "d") phi=phi*to_deg
       end if

       return
    End Subroutine Get_Cylindr_Coord_sp

    !!----
    !!---- Subroutine Get_Cart_from_Spher(r,Theta,Phi,Xo,Mode)
    !!----    real(kind=sp/dp),              intent( in)           :: r
    !!----    real(kind=sp/dp),              intent( in)           :: Theta
    !!----    real(kind=sp/dp),              intent( in)           :: Phi
    !!----    real(kind=sp/dp), dimension(3),intent(out)           :: xo
    !!----    character(len=*),              intent( in), optional :: mode
    !!----
    !!----    Determine the Cartesian coordinates from spherical coordinates.
    !!----    If Mode='D' the angle phi is provided in Degrees.
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine Get_Cart_from_Spher_dp(r,Theta,Phi,Xo,Mode)
    !!--++    real(kind=dp),              intent( in)           :: r
    !!--++    real(kind=dp),              intent( in)           :: Theta
    !!--++    real(kind=dp),              intent( in)           :: Phi
    !!--++    real(kind=dp), dimension(3),intent(out)           :: xo
    !!--++    character(len=*),           intent( in), optional :: mode
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the Cartesian coordinates from spherical coordinates.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Get_Cart_from_Spher_dp(r,Theta,Phi,Xo,Mode)
       !---- Arguments ----!
       real(kind=dp),              intent( in)           :: r
       real(kind=dp),              intent( in)           :: Theta
       real(kind=dp),              intent( in)           :: phi
       real(kind=dp), dimension(3),intent(out)           :: xo
       character(len=*),           intent( in), optional :: mode

       !---- Local Variables ----!
       real(kind=dp) :: ph,th

       ph=Phi
       th=Theta
       if (present(mode)) then
          if (mode(1:1) == "D" .or. mode(1:1) == "d") then
             ph=Phi*to_rad
             th=Theta*to_rad
          end if
       end if
       xo(1)=r*cos(ph)*sin(th)
       xo(2)=r*sin(ph)*sin(th)
       xo(3)=r*cos(th)

       return
    End Subroutine Get_Cart_from_Spher_dp

    !!--++
    !!--++ Subroutine Get_Cart_from_Spher_sp(r,Theta,Phi,Xo,Mode)
    !!--++    real(kind=sp),              intent( in)           :: r
    !!--++    real(kind=sp),              intent( in)           :: Theta
    !!--++    real(kind=sp),              intent( in)           :: Phi
    !!--++    real(kind=sp), dimension(3),intent(out)           :: xo
    !!--++    character(len=*),           intent( in), optional :: mode
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the Cartesian coordinates from spherical coordinates.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Get_Cart_from_Spher_sp(r,Theta,Phi,Xo,Mode)
       !---- Arguments ----!
       real(kind=sp),              intent( in)           :: r
       real(kind=sp),              intent( in)           :: Theta
       real(kind=sp),              intent( in)           :: phi
       real(kind=sp), dimension(3),intent(out)           :: xo
       character(len=*),           intent( in), optional :: mode

       !---- Local Variables ----!
       real(kind=sp) :: ph,th

       ph=Phi
       th=Theta
       if (present(mode)) then
          if (mode(1:1) == "D" .or. mode(1:1) == "d") then
             ph=Phi*to_rad
             th=Theta*to_rad
          end if
       end if
       xo(1)=r*cos(ph)*sin(th)
       xo(2)=r*sin(ph)*sin(th)
       xo(3)=r*cos(th)

       return
    End Subroutine Get_Cart_from_Spher_sp

    !!----
    !!---- Subroutine Get_Plane_from_Points(P1,P2,P3,A,B,C,D)
    !!----    real(kind=cp), dimension(3), intent(in) :: P1
    !!----    real(kind=cp), dimension(3), intent(in) :: P2
    !!----    real(kind=cp), dimension(3), intent(in) :: P3
    !!----    real(kind=cp),               intent(out):: A
    !!----    real(kind=cp),               intent(out):: B
    !!----    real(kind=cp),               intent(out):: C
    !!----    real(kind=cp),               intent(out):: D
    !!----
    !!----    Caculate the implicit form of a Plane in 3D as
    !!----    A * X + B * Y + C * Z + D = 0
    !!----
    !!---- Update: July - 2005
    !!
    Subroutine Get_Plane_from_Points(P1, P2, P3, A, B, C, D)
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent(in) :: P1
       real(kind=cp), dimension(3), intent(in) :: P2
       real(kind=cp), dimension(3), intent(in) :: P3
       real(kind=cp),               intent(out):: A
       real(kind=cp),               intent(out):: B
       real(kind=cp),               intent(out):: C
       real(kind=cp),               intent(out):: D

       a = ( p2(2) - p1(2) ) * ( p3(3) - p1(3) ) &
           - ( p2(3) - p1(3) ) * ( p3(2) - p1(2) )

       b = ( p2(3) - p1(3) ) * ( p3(1) - p1(1) ) &
           - ( p2(1) - p1(1) ) * ( p3(3) - p1(3) )

       c = ( p2(1) - p1(1) ) * ( p3(2) - p1(2) ) &
           - ( p2(2) - p1(2) ) * ( p3(1) - p1(1) )

       d = - p2(1) * a - p2(2) * b - p2(3) * c

       return
    End Subroutine Get_Plane_from_Points

    !!----
    !!---- Subroutine Get_Spheric_Coord(Xo,Ss,Theta,Phi,Mode)
    !!----    real(kind=sp/dp), dimension(3),intent( in)           :: xo
    !!----    real(kind=sp/dp),              intent(out)           :: ss
    !!----    real(kind=sp/dp),              intent(out)           :: theta
    !!----    real(kind=sp/dp),              intent(out)           :: phi
    !!----    character(len=*),              intent( in), optional :: mode
    !!----
    !!----    Determine the spheric coordinates from rectangular coordinates.
    !!----    If Mode='D' the angles will be done in Degrees.
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine Get_Spheric_Coord_dp(Xo,Ss,Theta,Phi,Mode)
    !!--++    real(kind=dp), dimension(3),intent( in)           :: xo
    !!--++    real(kind=dp),              intent(out)           :: ss
    !!--++    real(kind=dp),              intent(out)           :: theta
    !!--++    real(kind=dp),              intent(out)           :: phi
    !!--++    character(len=*),           intent( in), optional :: mode
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the spheric coordinates from rectangular coordinates
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Get_Spheric_Coord_dp(xo,ss,theta,phi,mode)
       !---- Arguments ----!
       real(kind=dp), intent( in), dimension(3)   :: xo
       real(kind=dp), intent(out)                 :: ss
       real(kind=dp), intent(out)                 :: theta
       real(kind=dp), intent(out)                 :: phi
       character(len=*), intent(in), optional     :: mode

       !---- Local Variables ----!
       integer :: j

       ss=0.0_dp
       do j=1,3
          ss=ss+xo(j)*xo(j)
       end do
       ss=sqrt(ss)
       if (ss > 0.0_dp) then
          theta=xo(3)/ss
          if (abs(theta) > 1.0_dp) then
             theta=sign(1.0_dp,theta)
          end if
          theta=acos(theta)
          if (abs(theta) < eps .or. abs(theta-pi) < eps) then
             phi=0.0_dp
          else
             phi=atan2(xo(2),xo(1))
          end if
       else
          theta=0.0_dp
          phi=0.0_dp
       end if
       if (present(mode)) then
          if (mode(1:1) == "D" .or. mode(1:1) == "d") then
             theta=theta*to_deg
             phi=phi*to_deg
          end if
       end if

       return
    End Subroutine Get_Spheric_Coord_dp

    !!--++
    !!--++ Subroutine Get_Spheric_Coord_sp(Xo,Ss,Theta,Phi,Mode)
    !!--++    real(kind=sp), dimension(3),intent( in)           :: xo
    !!--++    real(kind=sp),              intent(out)           :: ss
    !!--++    real(kind=sp),              intent(out)           :: theta
    !!--++    real(kind=sp),              intent(out)           :: phi
    !!--++    character(len=*),           intent( in), optional :: mode
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the spheric coordinates from rectangular coordinates
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Get_Spheric_Coord_sp(xo,ss,theta,phi,mode)
       !---- Arguments ----!
       real(kind=sp), intent( in), dimension(3)   :: xo
       real(kind=sp), intent(out)                 :: ss
       real(kind=sp), intent(out)                 :: theta
       real(kind=sp), intent(out)                 :: phi
       character(len=*), intent(in), optional     :: mode

       !---- Local Variables ----!
       integer :: j

       ss=0.0_sp
       do j=1,3
          ss=ss+xo(j)*xo(j)
       end do
       ss=sqrt(ss)
       if (ss > 0.0_sp) then
          theta=xo(3)/ss
          if (abs(theta) > 1.0_sp) then
             theta=sign(1.0_sp,theta)
          end if
          theta=acos(theta)
          if (abs(theta) < eps .or. abs(theta-pi) < eps) then
             phi=0.0_sp
          else
             phi=atan2(xo(2),xo(1))
          end if
       else
          theta=0.0_sp
          phi=0.0_sp
       end if
       if (present(mode)) then
          if (mode(1:1) == "D" .or. mode(1:1) == "d") then
             theta=theta*to_deg
             phi=phi*to_deg
          end if
       end if

       return
    End Subroutine Get_Spheric_Coord_sp

    !!----
    !!---- Subroutine Matrix_DiagEigen(A, V, C)
    !!----    real(kind=cp), dimension(3,3), intent(in)  :: a
    !!----    real(kind=cp), dimension(3),   intent(out) :: v
    !!----    real(kind=cp), dimension(3,3), intent(out) :: c
    !!----
    !!----    Diagonalize the matrix A, put eigenvalues in V and
    !!----    eigenvectors in C
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Matrix_DiagEigen(a,v,c)
       !---- Arguments ----!
       real(kind=cp), intent(in)  , dimension(3,3)    :: a
       real(kind=cp), intent(out) , dimension(3)      :: v
       real(kind=cp), intent(out) , dimension(3,3)    :: c

       !---- Local Variables ----!
       integer, parameter            :: n=3
       integer                       :: i, j, k, itmax, nm1, ip1, iter
       real(kind=cp), dimension(3)   :: u
       real(kind=cp), dimension(3,3) :: e
       real(kind=cp), parameter      :: eps1=1.e-7 , eps2=1.e-7 , eps3=1.e-7
       real(kind=cp)                 :: sigma1, offdsq, p, q, spq, csa, sna
       real(kind=cp)                 :: holdik, holdki, sigma2

       call init_err_math3d()
       nm1=n-1
       itmax=50
       do i=1,n
          do j=1,n
             e(i,j)=a(i,j)
             c(i,j)=0.0
             if (j < i) e(i,j)=0.0
          end do
       end do
       sigma1=0.0
       offdsq=0.0

       do i=1,n
          sigma1=sigma1+e(i,i)**2
          c(i,i)=1.0
          ip1=i+1
          if (i >= n) exit
          do j=ip1,n
             offdsq=offdsq+e(i,j)**2
          end do
       end do

       do iter=1,itmax
          do i=1,nm1
             ip1=i+1
             do j=ip1,n
                q=abs(e(i,i)-e(j,j))
                if (q <= eps1) then
                   csa=1.0/sqrt(2.0)
                   sna=csa
                else
                   if (abs(e(i,j)) <= eps2) then
                      e(i,j)=0.0
                      cycle
                   end if
                   p=2.0*e(i,j)*q/(e(i,i)-e(j,j))
                   spq=sqrt(p*p+q*q)
                   csa=sqrt((1.0+q/spq)/2.0)
                   sna=p/(2.0*csa*spq)
                end if
                do k=1,n
                   holdki=c(k,i)
                   c(k,i)=holdki*csa+c(k,j)*sna
                   c(k,j)=holdki*sna-c(k,j)*csa
                end do
                do k=i,n
                   if (k > j) then
                      holdik=e(i,k)
                      e(i,k)=csa*holdik+sna*e(j,k)
                      e(j,k)=sna*holdik-csa*e(j,k)
                   else
                      u(k)=e(i,k)
                      e(i,k)=csa*u(k)+sna*e(k,j)
                      if (k /= j) cycle
                      e(j,k)=sna*u(k)-csa*e(j,k)
                   end if
                end do
                u(j)=sna*u(i)-csa*u(j)
                do k=1,j
                   if (k <= i)  then
                      holdki=e(k,i)
                      e(k,i)=csa*holdki+sna*e(k,j)
                      e(k,j)=sna*holdki-csa*e(k,j)
                   else
                      e(k,j)=sna*u(k)-csa*e(k,j)
                   end if
                end do
                e(i,j)=0.0
             end do
          end do
          sigma2=0.0
          do i=1,n
             v(i)=e(i,i)
             sigma2=sigma2+v(i)*v(i)
          end do
          if (1.0-sigma1/sigma2 <= eps3) return
          sigma1=sigma2
       end do

       ERR_Math3D =.true.
       ERR_Math3D_Mess=" Convergence not reached in diagonalization "

       return
    End Subroutine Matrix_DiagEigen

    !!----
    !!---- Subroutine Matrix_Inverse(A, B, Ifail)
    !!----    real(kind=cp), dimension(3,3), intent(in)  :: a
    !!----    real(kind=cp), dimension(3,3), intent(out) :: b
    !!----    integer                      , intent(out) :: ifail
    !!----                                                  0 = OK; 1 = Fail
    !!----
    !!----    Inverts a 3x3 Matrix
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Matrix_Inverse(a,b,ifail)
       !---- Argument ----!
       real(kind=cp), dimension(3,3), intent(in)  :: a
       real(kind=cp), dimension(3,3), intent(out) :: b
       integer                      , intent(out) :: ifail

       !---- Local variables ----!
       real(kind=cp), parameter :: epso=1.0e-20
       real(kind=cp)            :: dmat

       ifail=0
       call init_err_math3d()

       b(1,1) = a(2,2)*a(3,3)-a(2,3)*a(3,2)
       b(2,1) = -(a(2,1)*a(3,3)-a(2,3)*a(3,1))
       b(3,1) = a(2,1)*a(3,2)-a(2,2)*a(3,1)
       b(1,2) = -(a(1,2)*a(3,3)-a(1,3)*a(3,2))
       b(2,2) = a(1,1)*a(3,3)-a(1,3)*a(3,1)
       b(3,2) = -(a(1,1)*a(3,2)-a(1,2)*a(3,1))
       b(1,3) = a(1,2)*a(2,3)-a(1,3)*a(2,2)
       b(2,3) = -(a(1,1)*a(2,3)-a(1,3)*a(2,1))
       b(3,3) = a(1,1)*a(2,2)-a(1,2)*a(2,1)
       dmat = a(1,1)*b(1,1)+a(1,2)*b(2,1)+a(1,3)*b(3,1)

       if (abs(dmat) < epso) then
          ifail=1
          ERR_Math3D =.true.
          ERR_Math3D_Mess="Singular Matrix: inversion imposible"
          return
       end if

       b = b/dmat

       return
    End Subroutine Matrix_Inverse

    !!----
    !!---- Subroutine Resolv_Sist_1X2(W,T,Ts,X,Ix)
    !!----    integer,       dimension(2),      intent(in) :: w     !  In -> Input vector
    !!----    real(kind=cp),                    intent(in) :: t     !  In -> Input value
    !!----    real(kind=cp), dimension(2),      intent(out):: ts    ! Out -> Fixed value of solution
    !!----    real(kind=cp), dimension(2),      intent(out):: x     ! Out -> Fixed value for x,y
    !!----    integer, dimension(2),            intent(out):: ix    ! Out -> determine if solution
    !!----                                                                   1: x, 2: y, 3: z
    !!--<<
    !!----              w11 x1 + w12 x2  = t1
    !!----              x_sol(i)= ts(i) + x(i) ix(i)
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Resolv_Sist_1x2(w,t,ts,x,ix)
       !---- Arguments ----!
       integer,dimension(2), intent( in) :: w
       real(kind=cp),                 intent( in) :: t
       real(kind=cp), dimension(2),   intent(out) :: ts
       real(kind=cp), dimension(2),   intent(out) :: x
       integer,dimension(2), intent(out) :: ix

       !---- Initialize ----!
       ts = 0.0
       x  = 1.0
       ix = 0
       call init_err_math3d()

       !---- Both are zeros ----!
       if ( all(w == 0)) then
          if (abs(t) < eps) then
             ix(1)=1
             ix(2)=2
          else
             ERR_Math3D=.true.
             ERR_Math3D_Mess="Inconsistent solution (1x2)"
          end if
          return
       end if

       !---- Any is zero ----!
       if (any(w == 0)) then
          if ( w(1) == 0 ) then
             ix(1)=1
             ts(2)=t/real(w(2))
              x(2)=0.0
          else
             ts(1)=t/real(w(1))
              x(1)=0.0
             ix(2)=2
          end if
       else
          ix(1)=1
          ts(2)=t/real(w(2))
           x(2)=-real(w(1))/real(w(2))
          ix(2)=1
       end if

       return
    End Subroutine Resolv_Sist_1x2

    !!----
    !!---- Subroutine Resolv_Sist_1X3(W,T,Ts,X,Ix)
    !!----    integer, dimension(3),            intent(in) :: w     !  In -> Input vector
    !!----    real(kind=cp),                    intent(in) :: t     !  In -> Input value
    !!----    real(kind=cp), dimension(3),      intent(out):: ts    ! Out -> Fixed value of solution
    !!----    real(kind=cp), dimension(3),      intent(out):: x     ! Out -> Fixed value for x,y,z
    !!----    integer, dimension(3),            intent(out):: ix    ! Out -> determine if solution
    !!----                                                                   1: x, 2: y, 3: z
    !!--<<
    !!----               w11 x1 + w12 x2 + w13 x3 = t1
    !!----               x_sol(i)= ts(i) + x(i) ix(i)
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Resolv_Sist_1x3(w,t,ts,x,ix)
       !---- Arguments ----!
       integer,dimension(3), intent( in) :: w
       real(kind=cp),                 intent( in) :: t
       real(kind=cp), dimension(3),   intent(out) :: ts
       real(kind=cp), dimension(3),   intent(out) :: x
       integer,dimension(3), intent(out) :: ix

       !---- Local Variables ----!
       integer               :: i, zeros
       integer, dimension(2) :: w1
       integer, dimension(2) :: ix1
       real(kind=cp), dimension(2)    :: ts1
       real(kind=cp), dimension(2)    :: x1

       !---- Initialize ----!
       ts = 0.0
       x  = 1.0
       ix = 0
       call init_err_math3d()

       !---- Are there zeros? ----!
       zeros=0
       do i=1,3
          if (w(i) == 0) zeros=zeros+1
       end do
       select case (zeros)
          case (3)
             if (abs(t) < eps) then
                do i=1,3
                   ix(i)=i
                end do
             else
                ERR_Math3D=.true.
                ERR_Math3D_Mess="Inconsistent solution (1 x 3)"
             end if

          case (2)
             do i=1,3
                if (w(i) /= 0) then
                   ts(i)=t/real(w(i))
                   x(i) =0.0
                else
                   ix(i)=i
                end if
             end do

          case (1)
             do i=1,3
                if (w(i) == 0) exit
             end do
             select case (i)
                case (1)
                   w1=w(2:3)

                case (2)
                   w1(1)=w(1)
                   w1(2)=w(3)

                case (3)
                   w1=w(1:2)
             end select
             call resolv_sist_1x2(w1,t,ts1,x1,ix1)
             select case (i)
                case (1)
                   ix(1)  = 1
                   ts(2:3)= ts1
                   x(2:3) = x1
                   if (ix1(1)==1) ix(2)=2
                   if (ix1(1)==2) ix(2)=3
                   if (ix1(2)==1) ix(3)=2
                   if (ix1(2)==2) ix(3)=3

                  case (2)
                     ix(2)= 2
                     ts(1)= ts1(1)
                     ts(3)= ts1(2)
                     x(1) = x1(1)
                     x(3) = x1(2)
                     if (ix1(1)==1) ix(1)=1
                     if (ix1(1)==2) ix(1)=3
                     if (ix1(2)==1) ix(3)=1
                     if (ix1(2)==2) ix(3)=3

                  case (3)
                     ix(3)  = 3
                     ts(1:2)= ts1
                     x(1:2) = x1
                     ix(1:2)= ix1
               end select

          case (0)
             ERR_Math3D=.true.
             ERR_Math3D_Mess="Inconsistent case ax+by+cz=t (1x3)"
       end select

       return
    End Subroutine Resolv_Sist_1x3

    !!----
    !!---- Subroutine Resolv_Sist_2X2(W,T,Ts,X,Ix)
    !!----    integer, dimension(2,2),          intent(in) :: w     !  In -> Input vector
    !!----    real(kind=cp), dimension(2),      intent(in) :: t     !  In -> Input value
    !!----    real(kind=cp), dimension(2),      intent(out):: ts    ! Out -> Fixed value of solution
    !!----    real(kind=cp), dimension(2),      intent(out):: x     ! Out -> Fixed value for x,y
    !!----    integer, dimension(2),            intent(out):: ix    ! Out -> determine if solution
    !!----                                                                   1: x, 2: y, 3: z
    !!--<<
    !!----                 w11 x1 + w12 x2  = t1
    !!----                 w21 x1 + w22 x2  = t2
    !!----                 x_sol(i)= ts(i) + x(i) ix(i)
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Resolv_Sist_2x2(w,t,ts,x,ix)
       !---- Arguments ----!
       integer,dimension(2,2), intent( in) :: w
       real(kind=cp),dimension(2),      intent( in) :: t
       real(kind=cp),dimension(2),      intent(out) :: ts
       real(kind=cp),dimension(2),      intent(out) :: x
       integer,dimension(2),   intent(out) :: ix

       !---- Local Variables ----!
       integer                 :: i,deter
       integer, dimension(2)   :: zeros,colum
       real(kind=cp)           :: rden, rnum

       !---- Initialize ----!
       ts    = 0.0
       x     = 1.0
       ix    = 0
       call init_err_math3d()

       deter = w(1,1)*w(2,2) - w(1,2)*w(2,1)
       rden=real(deter)
       if (deter /= 0) then
          !---- X(1) ----!
          rnum=t(1)*w(2,2) - w(1,2)*t(2)
          ts(1)=rnum/rden

          !---- X(2) ----!
          rnum=w(1,1)*t(2) - t(1)*w(2,1)
          ts(2)=rnum/rden

          x =0.0

       else                        ! Singular Matrix
          !---- Are there zero rows? ----!
          zeros=0
          do i=1,2
             if (w(i,1) == 0 .and. w(i,2) == 0 )  zeros(i)=1
          end do
          select case (sum(zeros))
             case (2)
                if (abs(t(1)) <= eps .and. abs(t(2)) <= eps) then
                   ix(1)=1
                   ix(2)=2
                else
                   ERR_Math3D=.true.
                   ERR_Math3D_Mess="Inconsistent solution (2x2)"
                end if

             case (1)
                do i=1,2
                   if (zeros(i) == 0) exit
                end do
                call resolv_sist_1x2(w(i,:),t(i),ts,x,ix)

             case (0)
                !---- Are there zero columns? ----!
                colum=0
                do i=1,2
                   if (w(1,i) == 0 .and. w(2,i) == 0 ) colum(i)=1
                end do
                select case (sum(colum))
                   case (1)
                      do i=1,2
                         if (colum(i) == 0) exit
                      end do
                      if (w(1,i) /= 0) then
                         ts(i)=t(1)/real(w(1,i))
                      else
                         ts(i)=t(2)/real(w(2,i))
                      end if
                      x(i)=0.0
                      if (i == 1) then
                         ix(2)=2
                      else
                         ix(1)=1
                      end if

                   case (0)
                      call resolv_sist_1x2(w(1,:),t(1),ts,x,ix)

                end select
          end select
       end if

       return
    End Subroutine Resolv_Sist_2x2

    !!----
    !!---- Subroutine Resolv_Sist_2X3(W,T,Ts,X,Ix)
    !!----    integer, dimension(2,3),          intent(in) :: w      !  In -> Input vector
    !!----    real(kind=cp), dimension(2),      intent(in) :: t      !  In -> Input value
    !!----    real(kind=cp), dimension(3),      intent(out):: ts     ! Out -> Fixed value of solution
    !!----    real(kind=cp), dimension(3),      intent(out):: x      ! Out -> Fixed value for x,y
    !!----    integer, dimension(3),            intent(out):: ix     ! Out -> determine if solution
    !!----                                                                    1: x, 2: y, 3: z
    !!----               w11 x1 + w12 x2 + w13 x3 = t1
    !!----               w21 x1 + w22 x2 + w23 x3 = t2
    !!----               x_sol(i)= ts(i) + x(i) ix(i)
    !!----
    !!----   Update: February - 2005
    !!
    Subroutine Resolv_Sist_2x3(w,t,ts,x,ix)
       !---- Arguments ----!
       integer,dimension(2,3),          intent( in) :: w
       real(kind=cp),dimension(2),      intent( in) :: t
       real(kind=cp),dimension(3),      intent(out) :: ts
       real(kind=cp),dimension(3),      intent(out) :: x
       integer,dimension(3),            intent(out) :: ix

       !---- Local Variables ----!
       integer                 :: i, j
       integer, dimension(2)   :: fila
       integer, dimension(2)   :: ix1
       integer, dimension(3)   :: colum
       integer, dimension(2,2) :: w1
       integer, dimension(2,3) :: wm
       integer, dimension(2)   :: wc
       real(kind=cp)                    :: tc
       real(kind=cp), dimension(2)      :: tm
       real(kind=cp), dimension(2)      :: ts1, x1

       !---- Initialize ----!
       ts    = 0.0
       x     = 1.0
       ix    = 0
       call init_err_math3d()

       !---- Are there zero columns? ----!
       colum=0
       do i=1,3
            if (all(w(:,i) == 0)) colum(i)=1
       end do
       select case (sum(colum))
          case (3)
             if (abs(t(1)) <= eps .and. abs(t(2)) <= eps) then
                do i=1,3
                   ix(i)=i
                end do
             else
                ERR_Math3D=.true.
                ERR_Math3D_Mess="Inconsistent solution in (2x3)"
             end if

          case (2)
             do i=1,3
                if (colum(i) == 0) exit
             end do
             if (w(1,i) /= 0) then
                ts(i)=t(1)/real(w(1,i))
             else
                ts(i)=t(2)/real(w(2,i))
             end if
             x(i)=0.0
             select case (i)
                case (1)
                   ix(2)=2
                   ix(3)=3

                case (2)
                   ix(1)=1
                   ix(3)=3

                case (3)
                   ix(1)=1
                   ix(2)=2
             end select

          case (1)
             do i=1,3
                if (colum(i) == 1) exit
             end do
             select case (i)
                case (1)
                   w1=w(:,2:3)

                case (2)
                   w1(1,1)=w(1,1)
                   w1(1,2)=w(1,3)
                   w1(2,1)=w(2,1)
                   w1(2,2)=w(2,3)

                case (3)
                   w1=w(:,1:2)
             end select
             call resolv_sist_2x2(w1,t,ts1,x1,ix1)
             select case (i)
                case (1)
                   ix(1)  = 1
                   ts(2:3)= ts1
                   x (2:3)= x1
                   if (ix1(1) == 1) ix(2)=2
                   if (ix1(1) == 2) ix(2)=3
                   if (ix1(2) == 1) ix(3)=2
                   if (ix1(2) == 2) ix(3)=3

                case (2)
                   ix(2)=2
                   ts(1)=ts1(1)
                   ts(3)=ts1(2)
                   x(1) = x1(1)
                   x(3) = x1(2)
                   if (ix1(1) == 1) ix(1)=1
                   if (ix1(1) == 2) ix(1)=3
                   if (ix1(2) == 1) ix(3)=1
                   if (ix1(2) == 2) ix(3)=3

                case (3)
                   ix(3)  = 3
                   ts(1:2)= ts1
                   x (1:2)= x1
                   ix(1:2)= ix1
             end select

          case (0)
             !---- Are there zeros in any element of rows? ----!
             fila = 0
             do i=1,2
                if (all(w(i,:)==0)) fila(i)=1
             end do
             select case (sum(fila))
                case (1)
                   if (w(1,1) /= 0) then
                      call resolv_sist_1x3(w(1,:),t(1),ts,x,ix)
                   else
                      call resolv_sist_1x3(w(2,:),t(2),ts,x,ix)
                   end if

                case (0)
                   fila = 0
                   wm   = w
                   tm   = t
                   !---- Are there zeros in any element of rows? ----!
                   do i=1,2
                      do j=1,3
                         if (w(i,j)==0) fila(i)=fila(i)+1
                      end do
                   end do
                   if ( fila(2) > fila(1) ) then
                      wm(1,:)=w(2,:)
                      wm(2,:)=w(1,:)
                      tm(1)  =t(2)
                      tm(2)  =t(1)
                          j  =fila(1)
                      fila(1)=fila(2)
                      fila(2)=j
                   end if
                   select case (fila(1))
                      case (2)
                         do i=1,3
                            if (wm(1,i) /= 0) exit
                         end do
                         ts(i)=tm(1)/real(wm(1,i))
                         x(i)=0.0
                         select case (i)
                            case (1)
                               wc(1)=wm(2,2)
                               wc(2)=wm(2,3)
                               tc=tm(2)-(wm(2,1)*ts(i))

                            case (2)
                               wc(1)=wm(2,1)
                               wc(2)=wm(2,3)
                               tc=tm(2)-(wm(2,2)*ts(i))

                            case (3)
                               wc(1)=wm(2,1)
                               wc(2)=wm(2,2)
                               tc=tm(2)-(wm(2,3)*ts(i))
                         end select
                         call resolv_sist_1x2(wc,tc,ts1,x1,ix1)
                         select case(i)
                            case (1)
                               ts(2:3)=ts1
                                x(2:3)=x1
                                if (ix1(1)==1) ix(2)=2
                                if (ix1(1)==2) ix(2)=3
                                if (ix1(2)==1) ix(3)=2
                                if (ix1(2)==2) ix(3)=3

                            case (2)
                               ts(1)=ts1(1)
                               ts(3)=ts1(2)
                                x(1)=x1(1)
                                x(3)=x1(2)
                                if (ix1(1)==1) ix(1)=1
                                if (ix1(1)==2) ix(1)=3
                                if (ix1(2)==1) ix(3)=1
                                if (ix1(2)==2) ix(3)=3

                            case (3)
                               ts(1:2)=ts1
                                x(1:2)=x1
                               ix(1:2)=ix1
                         end select

                      case (1)
                         do i=1,3
                            if (wm(1,i) == 0) exit
                         end do
                         select case (fila(2))
                            case (1)
                               do j=1,3
                                  if (wm(2,j) == 0) exit
                               end do
                               select case (i)
                                  case (1)             ! 0 en w(1,1)
                                     select case (j)
                                        case (2)
                                           wc(1)=-wm(2,1)/wm(2,3)
                                           wc(2)= wm(1,2)/wm(1,3)
                                           tc=tm(1)/real(wm(1,3)) - tm(2)/real(wm(2,3))
                                           call resolv_sist_1x2(wc,tc,ts1,x1,ix1)
                                           ts(1:2)=ts1
                                           x(1:2) =x1
                                           ix(1:2)=ix1
                                           if (ix(1) == 0) then
                                              ts(3)=tm(2)/real(wm(2,3)) - ts(1)*wm(2,1)/real(wm(2,3))
                                              x(3)=0.0
                                           else
                                              if (ix(2) == 0) then
                                                 ts(3)=tm(1)/real(wm(1,3)) - ts(2)*wm(1,2)/real(wm(1,3))
                                                 x(3)=0.0
                                              else
                                                 ts(3)=tm(2)/real(wm(2,3))
                                                 x(3)=-real(wm(2,1))/real(wm(2,3))
                                                 ix(3)=1

                                                 ts(2)=tc/real(wc(2))
                                                 x(2) =-real(wc(1))/real(wc(2))
                                                 ix(2)=1
                                              end if
                                           end if

                                        case (3)
                                           wc(1)=-wm(2,1)/wm(2,2)
                                           wc(2)= wm(1,3)/wm(1,2)
                                           tc=tm(1)/real(wm(1,2)) - tm(2)/real(wm(2,2))
                                           call resolv_sist_1x2(wc,tc,ts1,x1,ix1)
                                           ts(1)=ts1(1)
                                           ts(3)=ts1(2)
                                           x(1) =x1(1)
                                           x(3) =x1(2)
                                           if (ix1(1) == 1) ix(1)=1
                                           if (ix1(1) == 2) ix(1)=3
                                           if (ix1(2) == 1) ix(3)=1
                                           if (ix1(2) == 2) ix(3)=3
                                           if (ix(1) == 0) then
                                              ts(2)=tm(2)/real(wm(2,2)) - ts(1)*wm(2,1)/real(wm(2,2))
                                              x(2)=0.0
                                           else
                                              if (ix(3) == 0) then
                                                 ts(2)=tm(1)/real(wm(1,2)) - ts(3)*wm(1,3)/real(wm(1,2))
                                                 x(2)=0.0
                                              else
                                                 ts(2)=tm(2)/real(wm(2,2))
                                                 x(3)=-real(wm(2,1))/real(wm(2,2))
                                                 ix(2)=1

                                                 ts(3)=tc/real(wc(2))
                                                 x(3) =-real(wc(1))/real(wc(2))
                                                 ix(3)=1
                                              end if
                                           end if
                                     end select

                                  case (2)             ! 0 en w(1,2)
                                     select case (j)
                                        case (1)
                                           wc(1)= wm(1,1)/wm(1,3)
                                           wc(2)=-wm(2,2)/wm(2,3)
                                           tc=tm(1)/real(wm(1,3)) - tm(2)/real(wm(2,3))
                                           call resolv_sist_1x2(wc,tc,ts1,x1,ix1)
                                           ts(1:2)=ts1
                                           x(1:2) =x1
                                           ix(1:2)=ix1
                                           if (ix(1) == 0) then
                                              ts(3)=tm(1)/real(wm(1,3)) - ts(1)*wm(1,1)/real(wm(1,3))
                                              x(3)=0.0
                                           else
                                              if (ix(2) == 0) then
                                                 ts(3)=tm(2)/real(wm(2,3)) - ts(2)*wm(2,2)/real(wm(2,3))
                                                 x(3)=0.0
                                              else
                                                 ts(3)=tm(1)/real(wm(1,3))
                                                 x(3)=-real(wm(1,1))/real(wm(1,3))
                                                 ix(3)=1

                                                 ts(2)=tc/real(wc(2))
                                                 x(2) = -real(wc(1))/real(wc(2))
                                                 ix(2)= 1
                                              end if
                                           end if

                                        case (3)
                                           wc(1)=-wm(2,2)/wm(2,1)
                                           wc(2)= wm(1,3)/wm(1,1)
                                           tc=tm(1)/real(wm(1,1)) - tm(2)/real(wm(2,1))
                                           call resolv_sist_1x2(wc,tc,ts1,x1,ix1)
                                           ts(2:3)=ts1
                                           x(2:3) =x1
                                           if (ix1(1) == 1) ix(2)=2
                                           if (ix1(1) == 2) ix(2)=3
                                           if (ix1(2) == 1) ix(3)=2
                                           if (ix1(2) == 2) ix(3)=3
                                           if (ix(2) == 0) then
                                              ts(1)=tm(2)/real(wm(2,1)) - ts(2)*wm(2,2)/real(wm(2,1))
                                              x(1)=0.0
                                           else
                                              if (ix(3) == 0) then
                                                 ts(1)=tm(1)/real(wm(1,1)) - ts(3)*wm(1,3)/real(wm(1,1))
                                                 x(1)=0.0
                                              else
                                                 ix(1)=1

                                                 ts(2)=tm(2)/real(wm(2,2))
                                                 x(2) =-real(wm(2,1))/real(wm(2,2))
                                                 ix(2)=1

                                                 ts(3)=tm(1)/real(wm(1,3))
                                                 x(3) =-real(wm(1,1))/real(wm(1,3))
                                                 ix(3)=1
                                              end if
                                           end if
                                     end select

                                  case (3)             ! 0 en w(1,3)
                                     select case (j)
                                        case (1)
                                           wc(1)= wm(1,1)/wm(1,2)
                                           wc(2)=-wm(2,3)/wm(2,2)
                                           tc=tm(1)/real(wm(1,2)) - tm(2)/real(wm(2,2))
                                           call resolv_sist_1x2(wc,tc,ts1,x1,ix1)
                                           ts(1)=ts1(1)
                                           ts(3)=ts1(2)
                                           x(1) =x1(1)
                                           x(3) =x1(2)
                                           if (ix1(1) == 1) ix(1)=1
                                           if (ix1(1) == 2) ix(1)=3
                                           if (ix1(2) == 1) ix(3)=1
                                           if (ix1(2) == 2) ix(3)=3
                                           if (ix(1) == 0) then
                                              ts(2)=tm(1)/real(wm(1,2)) - ts(1)*wm(1,1)/real(wm(1,2))
                                              x(2)=0.0
                                           else
                                              if (ix(3) == 0) then
                                                 ts(2)=tm(2)/real(wm(2,2)) - ts(3)*wm(2,3)/real(wm(2,2))
                                                 x(2)=0.0
                                              else
                                                 ts(2)=tm(1)/real(wm(1,2))
                                                 x(2) =-real(wm(1,1))/real(wm(1,2))
                                                 ix(2)=1

                                                 ts(3)=tc/real(wc(2))
                                                 x(3) =-real(wc(1))/real(wc(2))
                                                 ix(3)=1
                                              end if
                                           end if

                                        case (2)
                                           wc(1)= wm(1,2)/wm(1,1)
                                           wc(2)=-wm(2,3)/wm(2,1)
                                           tc=tm(1)/real(wm(1,1)) - tm(2)/real(wm(2,1))
                                           call resolv_sist_1x2(wc,tc,ts1,x1,ix1)
                                           ts(2:3)=ts1
                                           x(2:3) =x1
                                           if (ix1(1) == 1) ix(2)=2
                                           if (ix1(1) == 2) ix(2)=3
                                           if (ix1(2) == 1) ix(3)=2
                                           if (ix1(2) == 2) ix(3)=3
                                           if (ix(2) == 0) then
                                              ts(1)=tm(1)/real(wm(1,1)) - ts(2)*wm(1,2)/real(wm(1,1))
                                              x(1)=0.0
                                           else
                                              if (ix(3) == 0) then
                                                 ts(1)=tm(2)/real(wm(2,1)) - ts(3)*wm(2,3)/real(wm(2,1))
                                                 x(1)=0.0
                                              else
                                                 ix(1)=1

                                                 ts(2)=tm(1)/real(wm(1,2))
                                                 x(2) =-real(wm(1,1))/real(wm(1,2))
                                                 ix(2)=1

                                                 ts(3)=tm(2)/real(wm(2,3))
                                                 x(3) =-real(wm(2,1))/real(wm(2,3))
                                                 ix(3)=1
                                              end if
                                           end if
                                     end select
                               end select

                            case (0)
                               select case (i)
                                  case (1)
                                     wc(1)=wm(2,1)
                                     wc(2)=wm(2,2)- wm(2,3)*wm(1,2)/wm(1,3)
                                     tc=tm(2)-tm(1)*wm(2,3)/real(wm(1,3))
                                     call resolv_sist_1x2(wc,tc,ts1,x1,ix1)
                                     ts(1:2)=ts1
                                     x(1:2)=x1
                                     ix(1:2)=ix1
                                     if (ix(2) == 0) then
                                        ts(3)=tm(1)/real(wm(1,3)) - ts(2)*real(wm(1,2))/real(wm(1,3))
                                        x(3)=0.0
                                     else
                                        ix(1)=1

                                        ts(2)=(tm(2) - tm(1)*wm(2,3)/real(wm(1,3))) / &
                                              (real(wm(2,2)) - real(wm(2,3)*wm(1,2))/real(wm(1,3)) )
                                        x(2) =-real(wm(2,1)) / &
                                              (real(wm(2,2)) - real(wm(2,3)*wm(1,2))/real(wm(1,3)) )
                                        ix(2)=1

                                        ts(3)= tm(1)/real(wm(1,3)) - (real(wm(1,2))/real(wm(1,3)))*ts(2)
                                        x(3) =- (real(wm(1,2))/real(wm(1,3)))*x(2)
                                        ix(3)=1
                                     end if

                                  case (2)
                                     wc(1)=wm(2,1)-wm(2,3)*wm(1,1)/wm(1,3)
                                     wc(2)=wm(2,2)
                                     tc=tm(2)-tm(1)*wm(2,3)/real(wm(1,3))
                                     call resolv_sist_1x2(wc,tc,ts1,x1,ix1)
                                    ts(1:2)=ts1
                                    x(1:2)=x1
                                    ix(1:2)=ix1
                                    if (ix(1) == 0) then
                                       ts(3)=tm(1)/real(wm(1,3)) - ts(1)*real(wm(1,1))/real(wm(1,3))
                                       x(3)=0.0
                                    else
                                       ix(1)=1

                                       ts(2)=(tm(2) - tm(1)*wm(2,3)/real(wm(1,3)))/real(wm(2,2))
                                       x(2) =(real(wm(1,1)*wm(2,3))/real(wm(1,3)) - real(wm(2,1)))/real(wm(2,2))
                                       ix(2)=1

                                       ts(3)=tm(1)/real(wm(1,3))
                                       x(3) =-real(wm(1,1))/real(wm(1,3))
                                       ix(3)=1
                                    end if

                                 case (3)
                                    wc(1)=wm(2,1)-wm(1,1)*wm(2,2)/wm(1,2)
                                    wc(2)=wm(2,3)
                                    tc=tm(2)-tm(1)*wm(2,2)/real(wm(1,2))
                                    call resolv_sist_1x2(wc,tc,ts1,x1,ix1)
                                    ts(1)=ts1(1)
                                    ts(3)=ts1(2)
                                    x(1)=x1(1)
                                    x(3)=x1(2)
                                    if (ix1(1) == 1) ix(1)=1
                                    if (ix1(1) == 2) ix(1)=3
                                    if (ix1(2) == 1) ix(3)=1
                                    if (ix1(2) == 2) ix(3)=3
                                    if (ix(1) == 0) then
                                       ts(2)=tm(1)/real(wm(1,2)) - ts(1)*real(wm(1,1))/real(wm(1,2))
                                       x(2)=0.0
                                    else
                                       ix(1) =1

                                       ts(2)=tm(1)/real(wm(1,2))
                                       x(2) =-real(wm(1,1))/real(wm(1,2))
                                       ix(2)=1

                                       ts(3)=(tm(2) - tm(1)*wm(2,2)/real(wm(1,2)))/real(wm(2,3))
                                       x(3) =(real(wm(1,1)*wm(2,2))/real(wm(1,2)) - real(wm(2,1)))/real(wm(2,3))
                                       ix(3)=1
                                    end if
                               end select
                         end select

                      case (0)
                         call resolv_sist_1x3(wm(1,:),tm(1),ts,x,ix)
                   end select

             end select
       end select

       return
    End Subroutine Resolv_Sist_2x3

    !!----
    !!---- Subroutine Resolv_Sist_3X3(W,T,Ts,X,Ix)
    !!----    integer, dimension(3,3),          intent(in) :: w      !  In -> Input vector
    !!----    real(kind=cp), dimension(3),      intent(in) :: t      !  In -> Input value
    !!----    real(kind=cp), dimension(3),      intent(out):: ts     ! Out -> Fixed value of solution
    !!----    real(kind=cp), dimension(3),      intent(out):: x      ! Out -> Fixed value for x,y
    !!----    integer, dimension(3),            intent(out):: ix     ! Out -> determine if solution
    !!----                                                                     1: x, 2: y, 3: z
    !!--<<
    !!----              w11 x1 + w12 x2 + w13 x3 = t1
    !!----              w21 x1 + w22 x2 + w23 x3 = t2
    !!----              w31 x1 + w32 x2 + w33 x3 = t3
    !!----              x_sol(i)= ts(i) + x(i) ix(i)
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Resolv_Sist_3x3(w,t,ts,x,ix)
       !---- Arguments ----!
       integer, dimension(3,3),          intent(in) :: w
       real(kind=cp), dimension(3),      intent(in) :: t
       real(kind=cp), dimension(3),      intent(out):: ts
       real(kind=cp), dimension(3),      intent(out):: x
       integer, dimension(3),            intent(out):: ix

       !---- Local variables ----!
       integer                 :: i,j,deter
       integer, dimension(3)   :: fila
       integer, dimension(3,3) :: w1
       integer, dimension(2,3) :: wm
       real(kind=cp)                    :: rnum, rden
       real(kind=cp), dimension(3)      :: t1
       real(kind=cp), dimension(2)      :: tm
       real(kind=cp),dimension(3,3)     :: rw

       !---- Initialize ----!
       ts  = 0.0
       x   = 1.0
       ix  = 0
       call init_err_math3d()

       deter=determ_a(w)
       rden=real(deter)

       if (deter /= 0) then
          !---- X(1) ----!
          rw=real(w)
          rw(:,1)=t
          rnum=determ_a(rw)
          ts(1)=rnum/rden

          !---- X(2) ----!
          rw=real(w)
          rw(:,2)=t
          rnum=determ_a(rw)
          ts(2)=rnum/rden

          !---- X(3) ----!
          rw=real(w)
          rw(:,3)=t
          rnum=determ_a(rw)
          ts(3)=rnum/rden

          x=0.0

       else                     !  Singular Matrix
          !---- Are there zero rows? ----!
          fila=0
          do i=1,3
             if (all(w(i,:) == 0)) fila(i)=1
          end do
          select case (sum(fila))
             !---- All values are zeros ----!
             case (3)
                if (all(abs(t) < eps)) then
                   do i=1,3
                      ix(i)=i
                   end do
                else
                   ERR_Math3D=.true.
                   ERR_Math3D_Mess="Inconsistent system (3 x 3)"
                end if

             !---- Two rows with zeroes ----!
             case (2)
                do i=1,3
                   if (fila(i) == 0) exit
                end do
                call resolv_sist_1x3(w(i,:),t(i),ts,x,ix)

             !---- One row with zeroes ----!
             case (1)
                do i=1,3
                   if (fila(i) == 1) exit
                end do
                select case(i)
                   case (1)
                      wm(1,:)=w(2,:)
                      wm(2,:)=w(3,:)
                      tm=t(2:3)

                   case (2)
                      wm(1,:)=w(1,:)
                      wm(2,:)=w(3,:)
                      tm(1)=t(1)
                      tm(2)=t(3)

                   case (3)
                      wm(1,:)=w(1,:)
                      wm(2,:)=w(2,:)
                      tm=t(1:2)

                end select
                call resolv_sist_2x3(wm,tm,ts,x,ix)

             !---- Non zero rows ----!
             case (0)
                w1=w
                t1=t

                !---- Are there 2 rows proportional? ----!
                do i=1,3
                   if ( abs(w1(1,i)) > abs(w1(2,i)) ) then
                      if (w1(2,i) /= 0) then
                         j=w1(1,i)/w1(2,i)
                      else
                         j=0
                      end if
                      if (j /= 0) then
                         if (j*w1(2,1) == w1(1,1) .and. j*w1(2,2) == w1(1,2) .and. &
                             j*w1(2,3) == w1(1,3) ) then
                            w1(1,:)=w1(2,:)
                            t1(1)  =t1(2)
                            exit
                         end if
                      end if
                   else
                      if (w1(1,i) /= 0) then
                         j=w1(2,i)/w1(1,i)
                      else
                         j=0
                      end if
                      if (j /= 0) then
                         if (j*w1(1,1) == w1(2,1) .and. j*w1(1,2) == w1(2,2) .and. &
                             j*w1(1,3) == w1(2,3) ) then
                            w1(2,:)=w1(1,:)
                            t1(2)  =t1(1)
                            exit
                         end if
                      end if
                   end if
                end do

                do i=1,3
                   if ( abs(w1(1,i)) > abs(w1(3,i)) ) then
                      if (w1(3,i) /= 0) then
                         j=w1(1,i)/w1(3,i)
                      else
                         j=0
                      end if
                      if (j /= 0) then
                         if (j*w1(3,1) == w1(1,1) .and. j*w1(3,2) == w1(1,2) .and. &
                             j*w1(3,3) == w1(1,3) ) then
                            w1(1,:)=w1(3,:)
                            t1(1)  =t1(3)
                            exit
                         end if
                      end if
                   else
                      if (w1(1,i) /= 0) then
                         j=w1(3,i)/w1(1,i)
                      else
                         j=0
                      end if
                      if (j /= 0) then
                         if (j*w1(1,1) == w1(3,1) .and. j*w1(1,2) == w1(3,2) .and. &
                             j*w1(1,3) == w1(3,3) ) then
                            w1(3,:)=w1(1,:)
                            t1(3)  =t1(1)
                            exit
                         end if
                      end if
                   end if
                end do

                do i=1,3
                   if ( abs(w1(2,i)) > abs(w1(3,i)) ) then
                      if (w1(3,i) /= 0) then
                         j=w1(2,i)/w1(3,i)
                      else
                         j=0
                      end if
                      if (j /= 0) then
                         if (j*w1(3,1) == w1(2,1) .and. j*w1(3,2) == w1(2,2) .and. &
                             j*w1(3,3) == w1(2,3) ) then
                            w1(2,:)=w1(3,:)
                            t1(2)  =t1(3)
                            exit
                         end if
                      end if
                   else
                      if (w1(2,i) /= 0) then
                         j=w1(3,i)/w1(2,i)
                      else
                         j=0
                      end if
                      if (j /= 0) then
                         if (j*w1(2,1) == w1(3,1) .and. j*w1(2,2) == w1(3,2) .and. &
                             j*w1(2,3) == w1(3,3) ) then
                            w1(3,:)=w1(2,:)
                            t1(3)  =t1(2)
                            exit
                         end if
                      end if
                   end if
                end do

                !---- Are there 3 rows equal? ----!
                if ( (w1(1,1) == w1(2,1)) .and. (w1(1,1) == w1(3,1)) .and. &
                     (w1(1,2) == w1(2,2)) .and. (w1(1,2) == w1(3,2)) .and. &
                     (w1(1,3) == w1(2,3)) .and. (w1(1,3) == w1(3,3)) ) then

                   call resolv_sist_1x3(w1(1,:),t1(1),ts,x,ix)

                !---- Are there 2 rows equal? ----!
                elseif( (w1(1,1) == w1(2,1)) .and. (w1(1,2) == w1(2,2)) .and. &
                        (w1(1,3) == w1(2,3)) ) then

                   call resolv_sist_2x3(w1(2:3,:),t1(2:3),ts,x,ix)

                elseif( (w1(1,1) == w1(3,1)) .and. (w1(1,2) == w1(3,2)) .and. &
                        (w1(1,3) == w1(3,3)) ) then

                   call resolv_sist_2x3(w1(1:2,:),t1(1:2),ts,x,ix)

                elseif( (w1(2,1) == w1(3,1)) .and. (w1(2,2) == w1(3,2)) .and. &
                        (w1(2,3) == w1(3,3)) ) then

                   call resolv_sist_2x3(w1(1:2,:),t1(1:2),ts,x,ix)

                !---- Are linear combinations? ----!
                else
                   call resolv_sist_2x3(w1(1:2,:),t1(1:2),ts,x,ix)

                end if

          end select
       end if

       return
    End Subroutine Resolv_Sist_3x3

 End Module CFML_Math_3D
