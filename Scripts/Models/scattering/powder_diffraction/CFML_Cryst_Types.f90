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
!!---- MODULE: CFML_Crystal_Metrics
!!----   INFO: Module to define crystallographic types and to provide
!!----         automatic crystallographic operations.
!!----
!!---- HISTORY
!!----    Update: 05/03/2011
!!----
!!--.. INFORMATION
!!--..
!!--..    List Of Matrix Relationships For Crystallographic Applications
!!--..
!!--..    Small "t" is for transpose, inv(F) is the inverse of matrix F
!!--..
!!--..    Basis vectors as symbolic matrices
!!--..       At = (a,b,c)  At'=(a',b',c') ;  At* = (a*,b*,c*)  At*'=(a*',b*',c*')
!!--..
!!--..    Direct and reciprocal metric tensors: G, G*=inv(G)
!!--..    X  column vector in     direct space, referred to basis A
!!--..    X* column vector in reciprocal space, referred to basis A*
!!--..
!!--..       A'  = M  A           X'  = inv(Mt) X
!!--..       A*  = G* A           X*  =   G     X
!!--..       A*' = inv(Mt) A*     X*' =   M     X*
!!--..
!!--..       G' = M G Mt          G*' = inv(Mt) G* inv(M)
!!--..
!!--..   Symmetry operator defined in bases: A, A', A*, A*'
!!--..       C = (R,T), C'= (R',T'), C*= (R*,T*), C*'= (R*',T*')
!!--..
!!--..       R'  = inv(Mt) R Mt  ; T' = inv(Mt) T
!!--..       R*' =  M  R* inv(M) ; T*' = M T*
!!--..       R*  = G R G*  = inv(Rt)
!!--..
!!--..   If a change of origin is performed the positions are changed
!!--..   Ot=(o1,o2,o3) origin of the new basis A' w.r.t. old basis A
!!--..
!!--..       X' = inv(Mt) (X-O)
!!--..
!!--..   Changing just the origin   Xn  = C  X  = R  X  + T
!!--..                              Xn' = C' X' = R' X' + T'
!!--..          R=R'                X'  = X -O
!!--..                              Xn' = Xn-O
!!--..                  Xn-O = R' (X-O) + T' = R X + T - O
!!--..                   R X - R O + T' = R X + T - O
!!--..                               T' = T - (O - R O) = T - (E-R)O
!!--..
!!--..   Changing the basis (A,o) -> (A',o')
!!--..                  Xn  = C  X  = R  X  + T
!!--..                  Xn' = C' X' = R' X' + T'
!!--..                  X'= inv(Mt) (X-O), Xn' = inv(Mt) (Xn-O)
!!--..
!!--..            inv(Mt) (Xn-O) = R' inv(Mt) (X-O) + T'
!!--..            inv(Mt) (R  X  + T -O) = R' inv(Mt) (X-O) + T'
!!--..            inv(Mt) R X + inv(Mt)(T-O) = R' inv(Mt) X - R' inv(Mt) O + T'
!!--..            inv(Mt) R = R' inv(Mt)  => R' = inv(Mt) R Mt
!!--..            inv(Mt) (T-O)  = - R' inv(Mt) O + T'
!!--..            T' = R' inv(Mt) O + inv(Mt) (T-O)
!!--..            T' = inv(Mt) R Mt inv(Mt) O + inv(Mt) (T-O)
!!--..            T' = inv(Mt) R  O + inv(Mt) (T-O)
!!--..            T' = inv(Mt) R  O + inv(Mt) T - inv(Mt) O
!!--..            T' = inv(Mt)( R  O + T -  O) = inv(Mt) (T -(E-R)O)
!!--..
!!--..
!!--..                       R' = inv(Mt) R Mt
!!--..
!!--..                       T' = inv(Mt) (T -(E-R)O)
!!--..
!!--..
!!--..   A symmetry operator does not change the modulus of vectors and
!!--..   the angles between vectors (dot product is invariant):
!!--..
!!--..      X' = R X ,  Y' = R Y  =>  Xt' = Xt Rt,  Yt' = Yt Rt
!!--..
!!--..      Xt' G Y' = Xt Rt G R Y = Xt G Y  =>  G = Rt G R
!!--..
!!--..
!!--..   Second rank tensor Q and Q* defined in bases A and A*.
!!--..
!!--..      Q' = M Q Mt      Q* = G* Q G*     Q*'= inv(Mt) Q* inv(M)
!!--..
!!--..   A symmetry operator R is equivalent to a transformation
!!--..   M = inv(Rt) acting on basis vectors => G' = inv(Rt) G inv(R) = G
!!--..   The anisotropic temperature factors Beta is defined in reciprocal
!!--..   space: is a tensor like Q*, the transformation of beta under
!!--..   a symmetry operator is then :
!!--..
!!--..           Beta' = Inv(Mt) Beta inv(M) = R Beta Rt
!!--..
!!----
!!---- DEPENDENCIES
!!--++    Use CFML_GlobalDeps,    only: Cp, Eps, Pi
!!--++    Use CFML_Math_General, only: Cosd, Sind, Acosd, Co_Prime, swap, Sort, atand, &
!!--++                                 Co_Linear
!!--++    Use CFML_Math_3D,      only : Matrix_Inverse, determ_A, determ_V, Cross_Product
!!----
!!---- VARIABLES
!!----    CRYSTAL_CELL_TYPE
!!----    TWOFOLD_AXES_TYPE
!!----    ZONE_AXIS_TYPE
!!----    ERR_CRYS
!!----    ERR_CRYS_MESS
!!--++    IDENTITY                       [Private]
!!--++    TPI2                           [Private]
!!----
!!---- PROCEDURES
!!----    Functions:
!!----       CART_U_VECTOR
!!----       CART_VECTOR
!!----       CELL_VOLUME_SIGMA
!!----       CONVERT_B_BETAS
!!----       CONVERT_B_U
!!----       CONVERT_BETAS_B
!!----       CONVERT_BETAS_U
!!----       CONVERT_U_B
!!----       CONVERT_U_BETAS
!!----       GET_BETAS_FROM_BISO
!!--++       METRICS                     [Private]
!!----       ROT_MATRIX
!!----       U_EQUIV
!!----
!!----    Subroutines:
!!----       CHANGE_SETTING_CELL
!!----       GET_BASIS_FROM_UVW
!!----       GET_CONVENTIONAL_CELL
!!----       GET_CRYST_FAMILY
!!--++       GET_CRYST_ORTHOG_MATRIX     [Private]
!!----       GET_DERIV_ORTH_CELL
!!----       GET_PRIMITIVE_CELL
!!----       GET_TRANSFM_MATRIX
!!----       GET_TWOFOLD_AXES
!!----       INIT_ERR_CRYS
!!----       NIGGLI_CELL                 [Overloaded]
!!--++       NIGGLI_CELL_ABC             [Private]
!!--++       NIGGLI_CELL_NIGGLIMAT       [Private]
!!--++       NIGGLI_CELL_PARAMS          [Private]
!!--++       NIGGLI_CELL_TYPE            [Private]
!!--++       NIGGLI_CELL_VECT            [Private]
!!----       READ_BIN_CRYSTAL_CELL
!!--++       RECIP                       [Private]
!!----       SET_CRYSTAL_CELL
!!----       VOLUME_SIGMA_FROM_CELL
!!----       WRITE_BIN_CRYSTAL_CELL
!!----       WRITE_CRYSTAL_CELL
!!----
!!
 Module CFML_Crystal_Metrics

    !---- Use files ----!
    Use CFML_GlobalDeps,   only : Cp, Eps, Pi, TO_RAD
    Use CFML_Math_General, only : Cosd, Sind, Acosd, Co_Prime, swap, Sort, atand, Co_Linear
    Use CFML_Math_3D,      only : Matrix_Inverse, determ_A, determ_V, Cross_Product

    implicit none

    private

    !---- List of public variables ----!

    !---- List of public functions ----!
    public :: Cart_u_vector, Cart_vector, Convert_B_Betas, Convert_B_U, &
              Convert_Betas_B, Convert_Betas_U, Convert_U_B,            &
              Convert_U_Betas, Rot_matrix, U_Equiv, Cell_Volume_Sigma,  &
              Get_Betas_From_Biso

    !---- List of public overloaded procedures: functions ----!

    !---- List of public subroutines ----!
    public :: Init_Err_Crys, Change_Setting_Cell,Set_Crystal_Cell,           &
              Get_Cryst_Family, Write_Crystal_Cell, Get_Deriv_Orth_Cell,     &
              Get_Primitive_Cell, Get_TwoFold_Axes, Get_Conventional_Cell,   &
              Get_Transfm_Matrix, Get_basis_from_uvw, Volume_Sigma_from_Cell,&
              Read_Bin_Crystal_Cell,Write_Bin_Crystal_Cell


    !---- List of public overloaded procedures: subroutines ----!

    public  :: Niggli_Cell

    !---- List of private functions ----!
    private :: metrics

    !---- List of private Subroutines ----!
    private :: Recip, Get_Cryst_Orthog_Matrix, Niggli_Cell_Vect, Niggli_Cell_Params, &
               Niggli_Cell_type, Niggli_Cell_abc,  Niggli_Cell_nigglimat

    !---- Definitions ----!

    !!----
    !!----  TYPE :: CRYSTAL_CELL_TYPE
    !!--..
    !!----  Type, public :: Crystal_Cell_Type
    !!----     real(kind=cp),dimension(3)   :: cell, ang          ! Direct cell parameters
    !!----     integer,      dimension(3)   :: lcell, lang        ! Code number for refinement in optimization procedures
    !!----     real(kind=cp),dimension(3)   :: cell_std, ang_std  ! Standar deviations cell parameters
    !!----     real(kind=cp),dimension(3)   :: rcell,rang         ! Reciprocal cell parameters
    !!----     real(kind=cp),dimension(3,3) :: GD,GR              ! Direct and reciprocal Metric Tensors
    !!----     real(kind=cp),dimension(3,3) :: Cr_Orth_cel        ! P-Matrix transforming Orthonormal
    !!----                                                        ! basis to direct Crytal cell (as I.T.)
    !!----                                                        ! (or crystallographic components to
    !!----                                                        !  Cartesian components: XC = Cr_Orth_cel X -> XC,X: column vectors)
    !!----     real(kind=cp),dimension(3,3) :: Orth_Cr_cel        ! Inv(Cr_Orth_cel) -> Cartesian to cryst. components
    !!----     real(kind=cp),dimension(3,3) :: BL_M               ! Busing-Levy B-matrix (transforms hkl to  a
    !!----                                                          Cartesian system with x//a*, y in (a*,b*) and z//c
    !!----     real(kind=cp),dimension(3,3) :: BL_Minv            ! Inverse of the Busing-Levy B-matrix
    !!----     real(kind=cp)                :: CellVol            ! Direct and Reciprocal
    !!----     real(kind=cp)                :: RCellVol           ! Cell volumes
    !!----     Character (len=1)            :: CartType           ! Cartesian Frame type: if CartType='A'
    !!----                                                        ! the Cartesian Frame has x // a.
    !!----  End Type Crystal_Cell_Type
    !!----
    !!---- Updated: November - 2013 (adding lcell and lang components)
    !!
    Type, public :: Crystal_Cell_Type
       real(kind=cp),dimension(3)   :: cell, ang
       integer,      dimension(3)   :: lcell, lang
       real(kind=cp),dimension(3)   :: cell_std, ang_std
       real(kind=cp),dimension(3)   :: rcell, rang
       real(kind=cp),dimension(3,3) :: GD,GR
       real(kind=cp),dimension(3,3) :: Cr_Orth_cel
       real(kind=cp),dimension(3,3) :: Orth_Cr_cel
       real(kind=cp),dimension(3,3) :: BL_M
       real(kind=cp),dimension(3,3) :: BL_Minv
       real(kind=cp)                :: CellVol
       real(kind=cp)                :: RCellVol
       character (len=1)            :: CartType
    End Type Crystal_Cell_Type

    !!----
    !!----  TYPE :: TWOFOLD_AXES_TYPE
    !!--..
    !!----  Type, public :: Twofold_Axes_Type
    !!----     integer                       :: ntwo        ! Number of two-fold axes
    !!----     real(kind=cp)                 :: tol         ! Angular tolerance (ca 3 degrees)
    !!----     real(kind=cp),dimension(3,12) :: caxes       ! Cartesian components of two-fold axes
    !!----     integer,dimension(3,12)       :: dtwofold    ! Direct indices of two-fold axes
    !!----     integer,dimension(3,12)       :: rtwofold    ! Reciprocal indices of two-fold axes
    !!----     integer,dimension(12)         :: dot         ! Scalar product of reciprocal and direct indices
    !!----     real(kind=cp),dimension(12)   :: cross       ! Angle between direct and reciprocal axes ( < tol)
    !!----     real(kind=cp),dimension(12)   :: maxes       ! Modulus of the zone axes (two-fold axes) vectors
    !!----     real(kind=cp),dimension(3)    :: a,b,c       ! Cartesian components of direct cell parameters
    !!----  End Type Twofold_Axes_Type
    !!----
    !!----  All components are initialised to zero in the type declaration
    !!----
    !!---- Update: October - 2008
    !!
    Type, public :: Twofold_Axes_Type
       integer                        :: ntwo=0
       real(kind=cp)                  :: tol=3.0
       real(kind=cp) ,dimension(3,12) :: caxes=0.0
       integer,dimension(3,12)        :: dtwofold=0
       integer,dimension(3,12)        :: rtwofold=0
       integer,dimension(12)          :: dot=0
       real(kind=cp), dimension(12)   :: cross=0.0
       real(kind=cp), dimension(12)   :: maxes=0.0
       real(kind=cp), dimension(3)    :: a=0.0,b=0.0,c=0.0
    End Type Twofold_Axes_Type

    !!----
    !!----  TYPE :: ZONE_AXIS_TYPE
    !!--..
    !!----  Type, public :: Zone_Axis_Type
    !!----    Integer               :: nlayer   ! number of the reciprocal layer considered normally nlayer=0
    !!----    Integer, dimension(3) :: uvw      ! Indices of the zone axis
    !!----    Integer, dimension(3) :: rx       ! Indices (reciprocal vector) of the basis vector 1
    !!----    Integer, dimension(3) :: ry       ! Indices (reciprocal vector) of the basis vector 2
    !!----  End Type Zone_Axis_Type
    !!----
    !!----
    !!----  This type comes from ResVis. It is useful to have it as a genereal type for
    !!----  many kinds of applications. Used in the subroutine Get_Basis_From_UVW.
    !!----
    !!---- Updated: February - 2012
    !!

    Type, public :: Zone_Axis_Type
      Integer               :: nlayer
      Integer, dimension(3) :: uvw
      Integer, dimension(3) :: rx
      Integer, dimension(3) :: ry
    End Type Zone_Axis_Type



    !!----
    !!---- ERR_CRYS
    !!----    logical, public :: Err_Crys
    !!----
    !!----    Logical Variable indicating an error in CFML_Crystal_Metrics module
    !!----
    !!---- Update: February - 2005
    !!
    logical, public          :: ERR_Crys

    !!----
    !!---- ERR_CRYS_MESS
    !!----    character(len=150), public :: ERR_Crys_Mess
    !!----
    !!----    String containing information about the last error
    !!----
    !!---- Update: February - 2005
    !!
    character(len=150), public :: ERR_Crys_Mess

    !!--++
    !!--++ IDENTITY
    !!--++    real(kind=cp), dimension(3,3), parameter :: identity=reshape ((/1.0,0.0,0.0,
    !!--++                                                                    0.0,1.0,0.0,
    !!--++                                                                    0.0,0.0,1.0/),(/3,3/))
    !!--++
    !!--++    (PRIVATE)
    !!--++    Identity matrix
    !!--++
    !!--++ Update: October - 2008
    !!
    real(kind=cp),dimension(3,3), parameter  :: identity=reshape ((/1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0/),(/3,3/))

    !!--++
    !!--++ TPI2
    !!--++    real(kind=cp), parameter :: tpi2=2.0*pi*pi
    !!--++
    !!--++    (PRIVATE)
    !!--++    Two times PI squared
    !!--++
    !!--++ Update: February - 2005
    !!
    real(kind=cp), parameter, private :: tpi2=2.0*pi*pi

    !---- Interfaces - Overloaded ----!

    !!--.. Three non coplanar vectors {a,b,c} generates a lattice using integer linear combinations
    !!--.. There are an infinite number of primitive unit cells generating the same lattice L.
    !!--.. N={a,b,c} is a Buerger cell if and only if |a|+|b|+|c| is a minimal value for all primitive
    !!--.. cells of L.
    !!--.. N is a Niggli cell of L if  (i) it is as Buerger cell of L and
    !!--..                            (ii) |90-alpha| + |90-beta| + |90-gamma| -> maximum
    !!--..                  / a.a  b.b  c.c \       /  s11  s22  s33 \
    !!--..   Niggli matrix  |               |   =   |                |
    !!--..                  \ b.c  a.c  a.b /       \  s23  s13  s12 /
    !!--..

    Interface  Niggli_Cell                   ! The first(s) argument(s) is(are)
      Module Procedure Niggli_Cell_abc       ! List of cell parameters passed as a 6D vector
      Module Procedure Niggli_Cell_nigglimat ! Niggli matrix passed as a 2x3 matrix (ultimately applying the algorithm)
      Module Procedure Niggli_Cell_Params    ! List of cell parameters a,b,c,alpha,beta,gamma
      Module Procedure Niggli_Cell_type      ! The object Cell is passed as argument
      Module Procedure Niggli_Cell_Vect      ! Input three vectors in Cartesian components
    End Interface  Niggli_Cell

 Contains

    !-------------------!
    !---- Functions ----!
    !-------------------!

    !!----
    !!---- Function Cart_U_Vector(Code,V,Celda) Result(Vc)
    !!----    character (len=*),             intent(in) :: code    !  In -> D: Direct, R: Reciprocal
    !!----    real(kind=cp), dimension(3),   intent(in) :: v       !  In -> Vector
    !!----    Type (Crystal_Cell_Type),      intent(in) :: Celda   !  In -> Cell Variable
    !!----    real(kind=cp),dimension(3)                :: vc      ! Out ->
    !!----
    !!----    Convert a vector in crystal space to unitary cartesian components
    !!----
    !!---- Update: February - 2005
    !!
    Function Cart_U_Vector(Code,V,Celda) Result(Vc)
       !---- Arguments ----!
       character (len=*),           intent(in) :: code
       real(kind=cp), dimension(3), intent(in) :: v
       type (Crystal_Cell_Type),    intent(in) :: Celda
       real(kind=cp), dimension(3)             :: vc

       !---- Local Variables ----!
       real(kind=cp) :: vmod

       vc=cart_vector(code,v,celda)
       vmod=sqrt(dot_product(vc,vc))
       if (vmod > 0.0) then
          vc=vc/vmod
       end if

       return
    End Function Cart_U_Vector

    !!----
    !!---- Function Cart_Vector(Code,V,Celda) Result(Vc)
    !!----    character (len=*),             intent(in) :: code     !  In -> D: Direct, R: Reciprocal
    !!----    real(kind=cp), dimension(3),   intent(in) :: v        !  In -> Vector
    !!----    Type (Crystal_Cell_Type),      intent(in) :: Celda    !  In -> Cell variable
    !!----    real(kind=cp) dimension(3)                :: vc       ! Out ->
    !!----
    !!----    Convert a vector in crystal space to cartesian components
    !!----    The value of code has been extended to use also the Busing-Levy
    !!----    Cartesian system as reference also for direct and reciprocal space.
    !!----    Codes:
    !!----    The Cartesian frame is that defined by the setting of the "Celda" object
    !!----         D: The components are given with respect to basis (a,b,c)
    !!----         R: The components are given with respect to basis (a*,b*,c*)
    !!----        BL: The components are given with respect to basis (a*,b*,c*) but
    !!----            the Cartesian frame is that defined by Busing and Levy
    !!----       BLD: The components are given with respect to basis (a,b,c) but
    !!----            the Cartesian frame is that defined by Busing and Levy
    !!----
    !!----
    !!---- Updated: June - 2012
    !!
    Function Cart_Vector(Code,V,Celda) Result(Vc)
       !---- Arguments ----!
       character (len=*),           intent(in) :: code
       real(kind=cp), dimension(3), intent(in) :: v
       type (Crystal_Cell_Type),    intent(in) :: Celda
       real(kind=cp), dimension(3)             :: vc

       select case (trim(code))
          case("d","D")
             vc = matmul(celda%Cr_Orth_cel,v)  !Direct conversion to Cartesian frame

          case ("r","R")
             vc = matmul(celda%GR,v)            !Converts to direct space
             vc = matmul(celda%Cr_Orth_cel,vc)  !Converts to Cartesian frame

          case ("bl","BL")
             vc = matmul(celda%BL_M,vc) !Direct conversion to BL Cartesian frame

          case ("bld","BLD")
             vc = matmul(celda%GD,v)   !Converts to reciprocal space
             vc = matmul(celda%BL_M,vc)!Converts to BL Cartesian frame

       end select

       return
    End Function Cart_Vector

    !!----
    !!---- Function Cell_Volume_Sigma(Cell) result(sigma)
    !!----   type(Crystal_Cell_Type), intent(in) :: Cell   !  In  ->  Cell variable
    !!----   real(kind=cp)                       :: sigma  !  Out ->  Sigma of volume
    !!----
    !!----    Calculates the standard deviation of the unit cell volume
    !!----    from the standard deviations of cell parameters. The input
    !!----    variable is of type Crytal_Cell_Type, if the standard deviations of
    !!----    of both cell axes and cell angles are zero the result is sigma=0.0,
    !!----    otherwise the calculation is performed
    !!----    It is assumed that there is no correlation (covariance terms) between
    !!----    the standard deviations of the different cell parameters.
    !!----
    !!---- Updated: January - 2013 (JRC)
    !!
    Function Cell_Volume_Sigma(Cell) result(sigma)
       !---- Arguments ----!
       type(Crystal_Cell_Type), intent(in) :: Cell
       real(kind=cp)                       :: sigma

       !--- Local variables ---!
       integer                     :: i
       real(kind=cp)               :: q,ca,cb,cc,vc,sa,sb,sc
       real(kind=cp), dimension(3) :: var_ang

       !> Check
       sigma=0.0
       if(sum(abs(Cell%cell_std)) < eps .and. sum(abs(Cell%ang_std)) < eps ) return

       vc=0.0
       do i=1,3
          q=Cell%cell_std(i)/Cell%cell(i)
          vc=vc+q*q
       end do
       if (sum(abs(Cell%ang_std)) > eps) then
          ca=cosd(Cell%ang(1)) ;  sa=sind(Cell%ang(1))
          cb=cosd(Cell%ang(2)) ;  sb=sind(Cell%ang(2))
          cc=cosd(Cell%ang(3)) ;  sc=sind(Cell%ang(3))
          q=1.0-ca*ca-cb*cb-cc*cc+2.0*ca*cb*cc
          var_ang = (Cell%ang_std * TO_RAD)**2/q
          vc=vc+ (ca-cb*cc)*(ca-cb*cc)*sa*sa * var_ang(1)
          vc=vc+ (cb-ca*cc)*(cb-ca*cc)*sb*sb * var_ang(2)
          vc=vc+ (cc-ca*cb)*(cc-ca*cb)*sc*sc * var_ang(3)
       end if

       sigma=Cell%Cellvol*sqrt(vc)

       return
    End Function Cell_Volume_Sigma

    !!--..
    !!--.. Betas are defined by the following expression of the temperature factor:
    !!--.. Taniso= exp( -(beta11 h^2 + beta22 k^2 + beta33 l^2 + 2 (beta12 h k + beta13 h l + beta23 k l)) )
    !!--.. Taniso= exp( -(bet(1) h^2 + bet(2) k^2 + bet(3) l^2 + 2 (bet(4) h k + bet(5) h l + bet(6) k l)) )
    !!--..
    !!--.. Us are defined by the following expression of the temperature factor:
    !!--.. Taniso= exp( -2pi^2 (h^2 (a*)^2 U11+ k^2 (b*)^2 U22+ l^2 (c*)^2 U33+
    !!--..                2 (h k (a*) (b*) U12+ h l (a*) (c*) U13+  k l (b*) (c*) U23)) )
    !!--..

    !!----
    !!---- Function Convert_B_Betas(B,Cell) Result(Beta)
    !!----    real(kind=cp),dimension(6), intent(in)  :: B
    !!----    type (Crystal_cell_Type),   intent(in)  :: Cell
    !!----    real(kind=cp),dimension(6)              :: Beta
    !!----
    !!----    Convert Thermal factors from B to Betas
    !!----
    !!---- Update: February - 2003
    !!
    Function Convert_B_Betas(B,Cell) Result(Beta)
       !---- Arguments ----!
       real(kind=cp),dimension(6), intent(in)  :: B
       type (Crystal_cell_Type),   intent(in)  :: Cell
       real(kind=cp),dimension(6)              :: Beta

       beta(1)=0.25*b(1)*cell%gr(1,1)                ! beta11
       beta(2)=0.25*b(2)*cell%gr(2,2)                ! beta22
       beta(3)=0.25*b(3)*cell%gr(3,3)                ! beta33
       beta(4)=0.25*b(4)*cell%rcell(1)*cell%rcell(2) ! beta12
       beta(5)=0.25*b(5)*cell%rcell(1)*cell%rcell(3) ! beta13
       beta(6)=0.25*b(6)*cell%rcell(2)*cell%rcell(3) ! beta23

       return
    End Function Convert_B_Betas

    !!----
    !!---- Function Convert_B_U(B) Result(U)
    !!----    real(kind=cp),dimension(6), intent(in)  :: B
    !!----    real(kind=cp),dimension(6)              :: U
    !!----
    !!----    Convert Thermal factors from B to U
    !!----
    !!---- Update: February - 2003
    !!
    Function Convert_B_U(B) Result(U)
       !---- Arguments ----!
       real(kind=cp),dimension(6),  intent(in)  :: B
       real(kind=cp),dimension(6)               :: U

       u=b/(4.0*tpi2)

       return
    End Function Convert_B_U

    !!----
    !!---- Function Convert_Betas_B(Beta,Cell) Result(B)
    !!----    real(kind=cp),dimension(6), intent(in)  :: Beta
    !!----    type (Crystal_cell_Type),   intent(in)  :: Cell
    !!----    real(kind=cp),dimension(6)              :: B
    !!----
    !!----    Convert Thermal factors from Betas to B
    !!----
    !!---- Update: February - 2003
    !!
    Function Convert_Betas_B(Beta,Cell) Result(B)
       !---- Arguments ----!
       real(kind=cp),dimension(6), intent(in)  :: Beta
       type (Crystal_cell_Type),   intent(in)  :: Cell
       real(kind=cp),dimension(6)              :: B

       b(1)=4.0*beta(1)/cell%gr(1,1)                  ! B11
       b(2)=4.0*beta(2)/cell%gr(2,2)                  ! B22
       b(3)=4.0*beta(3)/cell%gr(3,3)                  ! B33
       b(4)=4.0*beta(4)/(cell%rcell(1)*cell%rcell(2)) ! B12
       b(5)=4.0*beta(5)/(cell%rcell(1)*cell%rcell(3)) ! B13
       b(6)=4.0*beta(6)/(cell%rcell(2)*cell%rcell(3)) ! B23

       return
    End Function Convert_Betas_B

    !!----
    !!---- Function Convert_Betas_U(Beta,Cell) Result(U)
    !!----    real(kind=cp),dimension(6), intent(in)  :: Beta
    !!----    type (Crystal_cell_Type),   intent(in)  :: Cell
    !!----    real(kind=cp),dimension(6)              :: U
    !!----
    !!----    Convert Thermal factors from Betas to U
    !!----
    !!---- Update: February - 2003
    !!
    Function Convert_Betas_U(Beta,Cell) Result(U)
       !---- Arguments ----!
       real(kind=cp),dimension(6),intent(in)  :: Beta
       type (Crystal_cell_Type),  intent(in)  :: Cell
       real(kind=cp),dimension(6)             :: U

       u(1)=beta(1)/(tpi2*cell%gr(1,1))                ! U11
       u(2)=beta(2)/(tpi2*cell%gr(2,2))                ! U22
       u(3)=beta(3)/(tpi2*cell%gr(3,3))                ! U33
       u(4)=beta(4)/(tpi2*cell%rcell(1)*cell%rcell(2)) ! U12
       u(5)=beta(5)/(tpi2*cell%rcell(1)*cell%rcell(3)) ! U13
       u(6)=beta(6)/(tpi2*cell%rcell(2)*cell%rcell(3)) ! U23

       return
    End Function Convert_Betas_U

    !!----
    !!---- Function Convert_U_B(U) Result(B)
    !!----    real(kind=cp),dimension(6), intent(in)  :: U
    !!----    real(kind=cp),dimension(6)              :: B
    !!----
    !!----    Convert Thermal factors from U to B
    !!----
    !!---- Update: February - 2003
    !!
    Function Convert_U_B(U) Result(B)
       !---- Arguments ----!
       real(kind=cp),dimension(6),        intent(in)  :: U
       real(kind=cp),dimension(6)                     :: B

       b=4.0*tpi2*u

       return
    End Function Convert_U_B

    !!----
    !!---- Function Convert_U_Betas(U,Cell) Result(Beta)
    !!----    real(kind=cp),dimension(6), intent(in)  :: U
    !!----    type (Crystal_cell_Type),   intent(in)  :: Cell
    !!----    real(kind=cp),dimension(6)              :: Beta
    !!----
    !!----    Convert Thermal factors from U to Betas
    !!----
    !!---- Update: February - 2003
    !!
    Function Convert_U_Betas(U,Cell) Result(Beta)
       !---- Arguments ----!
       real(kind=cp),dimension(6),intent(in)  :: U
       type (Crystal_cell_Type),  intent(in)  :: Cell
       real(kind=cp),dimension(6)             :: Beta

       beta(1)=tpi2*u(1)*cell%gr(1,1)                ! beta11
       beta(2)=tpi2*u(2)*cell%gr(2,2)                ! beta22
       beta(3)=tpi2*u(3)*cell%gr(3,3)                ! beta33
       beta(4)=tpi2*u(4)*cell%rcell(1)*cell%rcell(2) ! beta12
       beta(5)=tpi2*u(5)*cell%rcell(1)*cell%rcell(3) ! beta13
       beta(6)=tpi2*u(6)*cell%rcell(2)*cell%rcell(3) ! beta23

       return
    End Function Convert_U_Betas

    !!----
    !!---- Function Get_Betas_from_Biso(Biso,Cell) Result(Betas)
    !!----    real(kind=cp),             intent(in)  :: Biso
    !!----    type (Crystal_cell_Type),  intent(in)  :: Cell
    !!----    real(kind=cp),dimension(6)             :: Betas
    !!----
    !!----    Get Betas from Biso
    !!----
    !!---- Update: April - 2013
    !!
    Function Get_Betas_from_Biso(Biso,Cell) Result(Betas)
       !--- Argument ----!
       real(kind=cp),             intent(in)  :: Biso
       type (Crystal_cell_Type),  intent(in)  :: Cell
       real(kind=cp),dimension(6)             :: Betas

       !---- Local variables ----!
       real(kind=cp), dimension (3,3) :: L,LT,U,bet
       integer                        :: i

       betas=0.0

       l=Cell%Orth_Cr_cel
       lt=Transpose(l)
       u = 0.0
       do i=1,3
          u(i,i) = 0.25*biso
       end do
       bet= matmul (l,lt)
       bet= matmul (bet,u)
       do i=1,3
          betas(i) = bet(i,i)
       end do

       betas(4) = bet(1,2)
       betas(5) = bet(1,3)
       betas(6) = bet(2,3)

       return
    End Function Get_Betas_from_Biso

    !!--++
    !!--++ Function Metrics(A,B) Result(G)
    !!--++    real(kind=cp), dimension(3)  , intent(in ) :: a   !  In -> Cell Parameters
    !!--++    real(kind=cp), dimension(3)  , intent(in ) :: b   !  In -> Ang Parameters
    !!--++    real(kind=cp), dimension(3,3)              :: g   ! Out -> Metrics array
    !!--++
    !!--++    (PRIVATE)
    !!--++    Constructs the metric tensor
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Metrics(A,B) Result(G)
       !---- Arguments ----!
       real(kind=cp), dimension(3)  , intent(in ) :: a
       real(kind=cp), dimension(3)  , intent(in ) :: b
       real(kind=cp), dimension(3,3)              :: g

       !---- Local Variables ----!
       integer :: i

       G(1,2)= a(1)*a(2)*cosd(b(3))
       G(1,3)= a(1)*a(3)*cosd(b(2))
       G(2,3)= a(2)*a(3)*cosd(b(1))

       do i=1,3
          G(i,i)= a(i)*a(i)
       end do

       G(2,1)=G(1,2)
       G(3,1)=G(1,3)
       G(3,2)=G(2,3)

       return
    End Function Metrics

    !!----
    !!---- Function Rot_Matrix(U, Phi, Celda)
    !!----    real(kind=cp), dimension(3),        intent(in) :: U
    !!----    real(kind=cp),                      intent(in) :: Phi
    !!----    type (Crystal_Cell_Type), optional, intent(in) :: Celda
    !!----    real(kind=cp), dimension(3,3)                  :: Rm
    !!----
    !!----    Returns the matrix (Gibbs matrix) of the active rotation of "phi" degrees
    !!----    along the "U" direction: R v = v', the vector v is tranformed to vector v'
    !!----    keeping the reference frame unchanged.
    !!----
    !!----    If one wants to calculate the components of the vector "v" in a rotated
    !!----    reference frame it suffices to invoke the function using "-phi".
    !!----    If "Celda" is present, "U" is in "Celda" coordinates,
    !!----    if not "U" is in cartesian coordinates.
    !!----
    !!----
    !!---- Update: February - 2005
    !!
    Function Rot_Matrix(U,Phi,Celda) Result(Rm)
       !---- Argument ----!
       real(kind=cp), dimension(3), intent(in)        :: U
       real(kind=cp), intent(in)                      :: phi
       type (Crystal_Cell_Type), optional, intent(in) :: Celda
       real(kind=cp), dimension(3,3)                  :: RM

       !---- Local variables ----!
       real(kind=cp)               :: c, s, umc, umod
       real(kind=cp), dimension(3) :: UU


       if (present(celda)) then
          uu= matmul(celda%cr_orth_cel,u)
       else
          uu=u
       end if

       umod=sqrt(dot_product(uu,uu))

       if (umod < tiny(1.0)) then
          uu=(/0.0,0.0,1.0/)
       else
          uu= uu/umod
       end if

       c= cosd(phi)
       s= sind(phi)
       umc = 1.0-c
       rm(1,1)= c+ umc*uu(1)**2
       rm(1,2)= umc*uu(1)*uu(2)- s*uu(3)
       rm(1,3)= umc*uu(1)*uu(3)+ s*uu(2)

       rm(2,1)= umc*uu(2)*uu(1)+ s*uu(3)
       rm(2,2)= c+ umc*uu(2)**2
       rm(2,3)= umc*uu(2)*uu(3)- s*uu(1)

       rm(3,1)= umc*uu(3)*uu(1)- s*uu(2)
       rm(3,2)= umc*uu(3)*uu(2)+ s*uu(1)
       rm(3,3)= c+ umc*uu(3)**2

       return
    End Function Rot_Matrix

    !!----
    !!---- Function U_Equiv(Cell, Th_U) Result(Uequi)
    !!----    type(Crystal_Cell_Type),    intent(in)     :: Cell    !  In -> Cell variable
    !!----    real(kind=cp), dimension(6),intent(in)     :: Th_U    !  In -> U parameters
    !!----
    !!----    Subroutine to obtain the U equiv from U11 U22 U33 U12 U13 U23
    !!----
    !!---- Update: February - 2005
    !!
    Function U_Equiv(Cell, Th_U) Result(Uequi)
       !---- Arguments ----!
       type (Crystal_cell_Type),    intent(in)  :: Cell
       real(kind=cp), dimension(6), intent(in)  :: Th_U
       real(kind=cp)                            :: Uequi

       !---- Local variables ----!
       real(kind=cp)    :: a, b, c, as, bs, cs, cosa, cosb, cosg
       real(kind=cp)    :: u11, u22, u33, u23, u13, u12

       a  =cell%cell(1)
       b  =cell%cell(2)
       c  =cell%cell(3)
       as =cell%rcell(1)
       bs =cell%rcell(2)
       cs =cell%rcell(3)
       cosa=cosd(cell%ang(1))
       cosb=cosd(cell%ang(2))
       cosg=cosd(cell%ang(3))

       u11=Th_u(1)
       u22=Th_u(2)
       u33=Th_u(3)
       u12=Th_u(4)
       u13=Th_u(5)
       u23=Th_u(6)
       uequi= (1.0/3.0) * (u11 * a * a * as * as + &
                           u22 * b * b * bs * bs + &
                           u33 * c * c * cs * cs + &
                           2.0*u12 * a * b * as * bs * cosg + &
                           2.0*u13 * a * c * as * cs * cosb + &
                           2.0*u23 * b * c * bs * cs * cosa )

       return
    End Function U_Equiv

    !---------------------!
    !---- Subroutines ----!
    !---------------------!

    !!----
    !!---- Subroutine Change_Setting_Cell(Cell,Mat,Celln,Matkind)
    !!----    type (Crystal_Cell_Type),      intent( in)    :: Cell
    !!----    real(kind=cp), dimension (3,3),intent( in)    :: Mat
    !!----    type (Crystal_Cell_Type),      intent(out)    :: Celln
    !!----    character (len=*), optional,   intent (in)    :: matkind
    !!----
    !!---- Calculates a new cell giving the transformation matrix.
    !!---- The input matrix can be given as the S-matrix in International
    !!---- Tables or its transposed (default) that corresponds to the matrix
    !!---- relating formal column matrices containing the basis vectors.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Change_Setting_Cell(Cell,Mat,Celln,Matkind)
       !---- Arguments ----!
       type (Crystal_Cell_Type),      intent( in)    :: Cell
       real(kind=cp), dimension (3,3),intent( in)    :: Mat
       type (Crystal_Cell_Type),      intent(out)    :: Celln
       character(len=*),  optional,   intent (in)    :: Matkind

       !--- Local variables ---!
       integer                       :: i
       real(kind=cp), dimension(3)   :: cellv,angl
       real(kind=cp), dimension(3,3) :: S,Gn,ST

       if (present(matkind)) then
          if (matkind(1:2) == "it" .or. matkind(1:2) == "IT" ) then
             S=Mat
            ST=transpose(Mat)
          else
             S=transpose(Mat)
            ST=Mat
          end if
       else
          S=transpose(Mat)
         ST=Mat
       end if

       !---- Get the new metric tensor
       !---- GDN= Mat GD MatT  or GDN= ST GD S
       gn=matmul(ST,matmul(Cell%GD,S))

       !---- Calculate new cell parameters from the new metric tensor
       do i=1,3
          Cellv(i)=sqrt(gn(i,i))
       end do
       angl(1)=acosd(Gn(2,3)/(cellv(2)*cellv(3)))
       angl(2)=acosd(Gn(1,3)/(cellv(1)*cellv(3)))
       angl(3)=acosd(Gn(1,2)/(cellv(1)*cellv(2)))

       !---- Construct the new cell
       call Set_Crystal_Cell(cellv,angl,Celln)

       return
    End Subroutine Change_Setting_Cell

    !!----
    !!---- Subroutine Get_basis_from_uvw(dmin,u,cell,ZoneB,ok,mode)
    !!----    real(kind=cp)             intent(in) :: dmin  !minimum d-spacing (smax=1/dmin)
    !!----    integer, dimension(3),    intent(in) :: u     !Zone axis indices
    !!----    type (Crystal_Cell_Type), intent(in) :: cell
    !!----    type (Zone_Axis_Type),    intent(out):: ZoneB !Object containing u and basis vector in the plane
    !!----    logical,                  intent(out):: ok
    !!----    character(len=*),optional,intent(in) :: mode
    !!----
    !!----  Subroutine to construct ZA of type Zone_Axis. This subroutine picks up two reciprocal
    !!----  lattice vectors satisfying the equation
    !!----                            hu+kv+lw=0
    !!----  The two reciprocal lattice vectors have no coprime factors and
    !!----  constitute the basis of a reciprocal lattice plane. They are
    !!----  obtained as the shortest two reciprocal lattice vectors satisfying
    !!----  the above equation. If mode is provided and mode="R", we interpret
    !!----  that the input zone axis is a reciprocal lattice vector and what we
    !!----  obtain is the basis of a direct plane in terms of lattice vectors.
    !!----  If mode="R", dmin corresponds n(uvw)max
    !!----  This subroutine has been imported from resvis_proc.f90.
    !!----
    !!----  Created: February 2006 (Imported from old programs for electron diffraction, Thesis JRC)
    !!----  Updated: February 2012 (JRC)
    !!----
    Subroutine Get_basis_from_uvw(dmin,u,cell,ZoneB,ok,mode)
       !--- Arguments ---!
       real(kind=cp),            intent(in) :: dmin
       integer, dimension(3),    intent(in) :: u
       type (Crystal_Cell_Type), intent(in) :: cell
       type (Zone_Axis_Type),    intent(out):: ZoneB
       logical,                  intent(out):: ok
       character(len=*),optional,intent(in) :: mode

       !--- Local Variables ---!
       integer                :: n,ik,il,um,iv,i1,i2,i,coun01,coun02,coun1,coun2
       integer,dimension(1)   :: i0
       integer                :: kmin,kmax,lmin,lmax
       integer,dimension(3)   :: au,h,mu
       real, dimension(2)     :: rm
       real, dimension(3,3)   :: mat
       integer,dimension(3,2) :: bas
       real                   :: rv,s2max

       ZoneB%nlayer=0
       ZoneB%uvw=u
       ok=.false.

       au=abs(u)
       um=3*maxval(au)
       i0=maxloc(au)

       i=i0(1)
       iv=u(i)
       mu=u
       if (iv < 0) then
         mu=-u
         iv=-iv
       end if

       Select Case (i)
         Case(1)
           i1=2; i2=3
         Case(2)
           i1=1; i2=3
         Case(3)
           i1=1; i2=2
       End Select

       rm(1)=100000.0; rm(2)=rm(1)
       bas(:,1) = (/ 71,121, 113/)
       bas(:,2) = (/117, 91,-111/)

       if(present(mode)) then
         s2max=dmin*dmin   !here dmin is really n_max
         kmax=nint(dmin/Cell%cell(i1)+1.0)
         lmax=nint(dmin/Cell%cell(i2)+1.0)
         kmax=min(um,kmax)
         lmax=min(um,lmax)
         mat=cell%gd
       else
         s2max=1.0/(dmin*dmin)
         kmax=nint(Cell%cell(i1)/dmin+1.0)
         lmax=nint(Cell%cell(i2)/dmin+1.0)
         kmax=min(um,kmax)
         lmax=min(um,lmax)
         mat=cell%gr
       end if

       kmin=-kmax; lmin=-lmax
       coun1=0; coun2=0
       do ik=kmax,kmin,-1
          do il=lmax,lmin,-1
             if (ik == 0 .and. il == 0) cycle
             n=-ik*mu(i1)-il*mu(i2)
             if (mod(n,iv) == 0) then               !n is multiple of iv
                h(i)= n/iv ; h(i1)=ik ; h(i2) = il  !h is solution of hu+kv+lw=0
                rv=dot_product(real(h),matmul(mat,real(h)))
                if (rv > s2max  .or. rv < 1.0e-20) cycle
                if (rv < rm(1)) then
                   if (.not. co_linear(bas(:,1),h,3) ) then
                      bas(:,2)=bas(:,1)
                      rm(2) = rm(1)
                      if (coun1 >=1) coun2=coun2+1
                   end if
                   bas(:,1)=h
                   rm(1) = rv
                   coun1=coun1+1
                else if (rv < rm(2) .and. .not. co_linear(bas(:,1),h,3) ) then
                   bas(:,2)=h
                   rm(2) = rv
                   coun2=coun2+1
                end if
             end if
          end do
       end do
       ZoneB%rx=bas(:,1)
       ZoneB%ry=bas(:,2)
       if (coun1 >= 1 .and. coun2 >=1) ok=.true.
       coun01=0; coun02=0; coun1=0; coun2=0
       do i=1,3
          if (ZoneB%rx(i) < 0) coun1=coun1+1
          if (ZoneB%ry(i) < 0) coun2=coun2+1
          if (ZoneB%rx(i) == 0) coun01=coun01+1
          if (ZoneB%ry(i) == 0) coun02=coun02+1
       end do
       if (coun1 >= 2 .or. (coun1 == 1 .and. coun01 == 2)) ZoneB%rx=-ZoneB%rx
       if (coun2 >= 2 .or. (coun2 == 1 .and. coun02 == 2)) ZoneB%ry=-ZoneB%ry

       return
    End Subroutine Get_Basis_From_Uvw

    !!----
    !!---- Subroutine Get_Conventional_Cell(Twofold,Cell,Tr,Message,Ok,told)
    !!----   Type(Twofold_Axes_Type), intent(in)  :: twofold
    !!----   Type(Crystal_Cell_Type), intent(out) :: Cell
    !!----   integer, dimension(3,3), intent(out) :: tr
    !!----   character(len=*),        intent(out) :: message
    !!----   logical,                 intent(out) :: ok
    !!----   real(kind=cp), optional, intent(in)  :: told
    !!----
    !!----  This subroutine provides the "conventional" (or quasi! being still tested )
    !!----  from the supplied object "twofold" that has been obtained from a previous
    !!----  call to Get_TwoFold_Axes. The conventional unit cell can be deduced from
    !!----  the distribution of two-fold axes in the lattice. The cell produced in this
    !!----  procedure applies some rules for obtaining the conventional cell, for instance
    !!----  in monoclinic lattices (a single two-fold axis) the two-fold axis is along
    !!----  b and the final cell is right handed with a <= c and beta >= 90. It may be
    !!----  A,C or I centred. The convertion to the C-centred setting in the A and I
    !!----  centring, is not attempted. The angular tolerance for accepting a two-fold
    !!----  axis, or higher order axes, as such has been previously set into twofold%tol
    !!----  component. The output Tr-matrix is the transpose of the IT convention.
    !!----  It corresponds to the transformation between formal column matrices containing
    !!----  the basis vectors.
    !!----  The tolerance for comparing distances in angstroms told is optional.
    !!----- By default the used tolerance is 0.2 angstroms.
    !!----
    !!---- Update: November - 2008
    !!----
    Subroutine Get_Conventional_Cell(Twofold,Cell,Tr,Message,Ok,told)
       !---- Arguments ----!
       Type(Twofold_Axes_Type), intent(in)  :: Twofold
       Type(Crystal_Cell_Type), intent(out) :: Cell
       integer, dimension(3,3), intent(out) :: tr
       character(len=*),        intent(out) :: message
       logical,                 intent(out) :: ok
       real(kind=cp), optional, intent(in)  :: told

       !---- Local variables ----!
       integer, dimension(1)          :: ix
       integer, dimension(2)          :: ab
       integer, dimension(3)          :: rw,h1,h2
       integer, dimension(66)         :: inp
       integer, dimension(3,48)       :: row
       real(kind=cp), dimension(3)    :: u,v1,v2,v3,a,b,c,vec,vi,vj,vk
       real(kind=cp), dimension(48)   :: mv
       real(kind=cp), dimension(66)   :: ang
       integer                        :: iu,iv,iw,nax,i,j,k,m,namina,naminb,naminc,ia
       real(kind=cp)                  :: dot,ep,domina,dominb,dominc,aij,aik,ajk
       real(kind=cp)                  :: delt,tola
       logical                        :: hexap, hexac

       a=twofold%a; b=twofold%b; c=twofold%c
       delt=twofold%tol
       ep=cosd(90.0-delt)
       domina=9.0e+30; dominc=domina
       tr=reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/))
       ab=0; mv=0.0; ang=0.0; row=0; inp=0
       ok=.true.
       tola=0.2
       if(present(told)) tola=told

       Select Case(twofold%ntwo)
          Case (1)    !Monoclinic
             v2=twofold%caxes(:,1)
             u = v2/twofold%maxes(1)
             tr(2,:)=twofold%dtwofold(:,1)
             nax=0
             do iu=-3,3
                do iv=-3,3
                   do_iw: do iw=0,3
                      rw=(/iu,iv,iw/)
                      ! if(iu == 0 .and. iv == 0 .and. iw == 0) cycle
                      if (.not. Co_prime(rw,3)) cycle
                      vec=real(iu)*a+real(iv)*b+real(iw)*c
                      dot=sqrt(dot_product(vec,vec))
                      vec=vec/dot
                      if (abs(dot_product(u,vec)) < ep) then
                         do m=1,nax
                            if(co_linear(rw,row(:,m),3)) cycle do_iw
                         end do
                         nax=nax+1
                         row(:,nax) = rw
                         mv(nax) = dot
                         if (dot < domina) then
                            domina=dot
                            namina=nax
                            tr(1,:)=rw
                            v1=real(iu)*a+real(iv)*b+real(iw)*c
                         end if
                      end if
                   end do do_iw
                end do
             end do

             do i=1,nax
                if (i == namina) cycle
                if (mv(i) < dominc) then
                   dominc=mv(i)
                   naminc=i
                end if
             end do
             tr(3,:)=row(:,naminc)
             v3=row(1,naminc)*a+row(2,naminc)*b+row(3,naminc)*c

             !Length of the three basis vectors should be stored in mv(1),mv(2),mv(3)
             mv(1)=sqrt(dot_product(v1,v1))
             mv(2)=sqrt(dot_product(v2,v2))
             mv(3)=sqrt(dot_product(v3,v3))

             !The two shortest vectors perpendicular to the primary twofold axis have been found
             !and the transformation matrix has been constructed
             namina=determ_A(tr)
             if (namina < 0) then   !right handed system
                tr(2,:)=-tr(2,:)
                v2=-v2
                namina=-namina
             end if

             !Test if beta is lower than 90 in such a case invert c and b
             dominb=dot_product(v1/mv(1),v3/mv(3))
             if (dominb > 0.0) then  !angle beta < 90
                tr(2,:)=-tr(2,:)
                v2=-v2
                tr(3,:)=-tr(3,:)
                v3=-v3
             end if

             Select Case (namina)
                Case(1)
                   message="Monoclinic, primitive cell"
                Case(2)
                   rw=matmul((/0,1,1/),tr)
                   if (.not. co_prime(rw,3)) then
                      message="Monoclinic, A-centred cell"
                   else
                      rw=matmul((/1,1,1/),tr)
                      if (.not. co_prime(rw,3)) then
                         message="Monoclinic, I-centred cell"
                      else
                         rw=matmul((/1,1,0/),tr)
                         if(.not. co_prime(rw,3)) message="Monoclinic, C-centred cell"
                      end if
                   end if

                Case(3:)
                   message="Error in monoclinic cell"
                   ok=.false.
             End Select

          Case (3)    !Orthorhombic/Trigonal
             u(1:3)=twofold%maxes(1:3)
             ix=minloc(u)
             namina=ix(1)
             ix=maxloc(u)
             naminc=ix(1)
             if (naminc == namina) then
                namina=1; naminb=2; naminc=3
             else
                do i=1,3
                   if(i == namina) cycle
                   if(i == naminc) cycle
                   naminb=i
                   exit
                end do
             end if
             tr(1,:) = twofold%dtwofold(:,namina)
             tr(2,:) = twofold%dtwofold(:,naminb)
             tr(3,:) = twofold%dtwofold(:,naminc)
             v1 = twofold%caxes(:,namina)
             v2 = twofold%caxes(:,naminb)
             v3 = twofold%caxes(:,naminc)
             mv(1)=twofold%maxes(namina)
             mv(2)=twofold%maxes(naminb)
             mv(3)=twofold%maxes(naminc)

             !Check the system by verifying that the two-fold axes form 90 (orthorhombic)
             !or 120 degrees (Trigonal)
             domina=dot_product(v2/mv(2),v3/mv(3))
             dominb=dot_product(v1/mv(1),v3/mv(3))
             dominc=dot_product(v1/mv(1),v2/mv(2))

             if (abs(domina) < ep .and. abs(dominb) < ep .and. abs(dominc) < ep) then !orthorhombic
                namina=determ_A(tr)
                if (namina < 0) then
                   tr(3,:)=-tr(3,:)
                   v3=-v3
                   namina=-namina
                end if
                Select Case (namina)
                   Case(1)
                      message="Orthorhombic, Primitive cell"

                   Case(2)
                      rw=matmul((/0,1,1/),tr)
                      if (.not. co_prime(rw,3)) then
                         message="Orthorhombic, A-centred cell"
                      else
                         rw=matmul((/1,1,1/),tr)
                         if (.not. co_prime(rw,3)) then
                            message="Orthorhombic, I-centred cell"
                         else
                            rw=matmul((/1,1,0/),tr)
                            if (.not. co_prime(rw,3)) then
                               message="Orthorhombic, C-centred cell"
                            else
                               rw=matmul((/1,0,1/),tr)
                               if (.not. co_prime(rw,3)) message="Orthorhombic, B-centred cell"
                            end if
                         end if
                      end if

                   Case(3:)
                      message="Orthorhombic, F-centred cell"
                End Select

             else !Rhombohedral/Trigonal

                !In the Trigonal system the two-fold axes are in the plane perpendicular to
                !the three-fold axis, and valid a,b, vectors can be chosen among any two two-fold
                !axes forming an angle of 120 degrees
                !verify that 1 and 2 form 120
                ang(1)=acosd(domina)    !2-3
                ang(2)=acosd(dominb)    !1-3
                ang(3)=acosd(dominc)    !1-2
                dot=1.0
                iu=1
                j=0
                do i=1,3
                   if (abs(ang(i)-120.0) < delt) then
                      j=i
                      exit
                   end if
                end do

                if ( j == 0) then
                   do i=1,3
                      if (abs(ang(i)-60.0) < delt) then
                         j=i
                         dot=-1.0
                         iu=-1
                         exit
                      end if
                   end do
                End if

                if ( j == 0) then
                   message="Trigonal/Rhombohedral test failed! Supply only one two-fold axis"
                   ok=.false.
                   return
                else
                   Select Case (j)
                      case(1)
                         vi=v2
                         vj=dot*v3
                         h1=tr(2,:); h2=iu*tr(3,:)
                         tr(3,:)=tr(1,:)
                         tr(1,:)=h1
                         tr(2,:)=h2

                      case(2)
                         vi=v1
                         vj=dot*v3
                         h2=iu*tr(3,:)
                         tr(3,:)=tr(2,:)
                         tr(2,:)=h2

                      case(3)
                         vi=v1
                         vj=dot*v2
                         tr(2,:)=iu*tr(2,:)

                   End Select

                   v1 = vi
                   v2 = vj
                   mv(1)=sqrt(dot_product(v1,v1))
                   mv(2)=sqrt(dot_product(v2,v2))
                   vi=v1/mv(1)
                   vj=v2/mv(2)
                   ok=.false.

                   do_iu: do iu=-3,3
                      do iv=-3,3
                         do iw=0,3
                            rw=(/iu,iv,iw/)
                            if (.not. Co_prime(rw,3)) cycle
                            vec=real(iu)*a+real(iv)*b+real(iw)*c
                            dot=sqrt(dot_product(vec,vec))
                            vec=vec/dot
                            if (abs(dot_product(vi,vec)) < ep  .and. abs(dot_product(vj,vec)) < ep) then
                               tr(3,:)=rw
                               ok=.true.
                               exit do_iu
                            end if
                         end do
                      end do
                   end do do_iu

                   If (ok) then
                      namina=determ_A(tr)
                      if (namina < 0) then
                         tr(3,:)=-tr(3,:)
                         namina=-namina
                      end if
                      v3 = tr(3,1)*a+tr(3,2)*b+tr(3,3)*c
                      mv(3)=sqrt(dot_product(v3,v3))
                      Select Case (namina)
                         case(1)
                            message="Primitive hexagonal cell"
                         case(3)
                            rw=matmul((/2,1,1/),tr)
                            if (.not. co_prime(rw,3)) then
                               message="Rhombohedral, obverse setting cell"
                            else
                               message="Rhombohedral, reverse setting cell"
                            end if
                      End Select

                   Else
                      message="Trigonal/Rhombohedral test failed! Supply only one two-fold axis"
                      ok=.false.
                      return
                   End if
                End if !j==0
             End if  !orthorhombic test

          Case (5)    !Tetragonal
             m=0
             inp=0
             mv(1:5)=twofold%maxes(1:5)
             do i=1,twofold%ntwo-1
                vi=twofold%caxes(:,i)/twofold%maxes(i)
                do j=i+1,twofold%ntwo
                   vj=twofold%caxes(:,j)/twofold%maxes(j)
                   m=m+1
                   ang(m)=acosd(dot_product(vi,vj))
                   if (abs(ang(m)-45.0) < delt .or. abs(ang(m)-135.0) < delt) then
                      inp(i)=1
                      inp(j)=1
                      if (mv(i) > mv(j)) then
                         ia=j
                      else
                         ia=i
                      end if
                      if (ab(1) == 0) then
                         ab(1) = ia
                      else
                         ab(2) = ia
                      end if
                   end if
                end do
             end do

             !Determination of the c-axis (that making 90 degree with all the others)
             ix=minloc(inp)
             naminc=ix(1)

             !The two axes forming a,b are those of indices ab(1) and ab(2)
             namina=ab(1)
             naminb=ab(2)
             if (namina == 0 .or. naminb == 0) then
                ok=.false.
                message="Basis vectors a-b not found!"
                return
             end if

             tr(1,:) = twofold%dtwofold(:,namina)
             tr(2,:) = twofold%dtwofold(:,naminb)
             tr(3,:) = twofold%dtwofold(:,naminc)
             v1 = twofold%caxes(:,namina)
             v2 = twofold%caxes(:,naminb)
             v3 = twofold%caxes(:,naminc)
             mv(1)=twofold%maxes(namina)
             mv(2)=twofold%maxes(naminb)
             mv(3)=twofold%maxes(naminc)
             namina=determ_A(tr)
             if (namina < 0) then
                tr(3,:)=-tr(3,:)
                v3=-v3
                namina=-namina
             end if

             Select Case (namina)
                Case(1)
                   message="Tetragonal, Primitive cell"
                Case(2)
                   message="Tetragonal, I-centred cell"
                Case(3:)
                   message="Error in tetragonal cell"
                   ok=.false.
                   return
             End Select

          Case (7)    !Hexagonal

             m=0
             inp=0
             mv(1:7)=twofold%maxes(1:7)
             hexap=.false.;  hexac=.false.

             !Search tha a-b plane
             do_ii:do i=1,twofold%ntwo-1
                vi=twofold%caxes(:,i)/twofold%maxes(i)
                do j=i+1,twofold%ntwo
                   vj=twofold%caxes(:,j)/twofold%maxes(j)
                   aij=acosd(dot_product(vi,vj))
                   if (abs(aij-120.0) < delt) then
                      if (abs(mv(i)-mv(j)) < tola .and. .not. hexap ) then
                         rw(1)=i; rw(2)=j
                         u(1)=mv(i); u(2)=mv(j)
                         hexap=.true.
                         exit do_ii
                      end if
                   end if
                end do
             end do do_ii

             if (hexap) then ! Search the c-axis, it should be also a two-fold axis!
                             ! because Op(6).Op(6).Op(6)=Op(2)
                v1 = twofold%caxes(:,rw(1))
                v2 = twofold%caxes(:,rw(2))
                vj=v1/u(1)
                vk=v2/u(2)
                do i=1,twofold%ntwo
                   vi=twofold%caxes(:,i)/twofold%maxes(i)
                   aij=acosd(dot_product(vi,vj))
                   aik=acosd(dot_product(vi,vk))
                   if (abs(aij-90.0) < delt .and. abs(aik-90.0) < delt ) then
                      rw(3)=i
                      u(3)= mv(i)
                      hexac=.true.
                      exit
                   end if
                end do
             else
                ok=.false.
                return
             end if

             if (hexac) then
                do i=1,3
                   tr(i,:) = twofold%dtwofold(:,rw(i))
                   mv(i)=u(i)
                end do
                v3 = twofold%caxes(:,rw(3))
                namina=determ_A(tr)
                if (namina < 0) then
                   tr(3,:)=-tr(3,:)
                   v3=-v3
                   namina=-namina
                end if

                Select Case (namina)
                   Case(1)
                      message="Hexagonal, Primitive cell"
                   Case(2:)
                      message="Hexagonal, centred cell? possible mistake"
                End Select

             else
                ok=.false.
                message="The c-axis of a hexagonal cell was not found!"
                return
             end if

          Case (9)   !Cubic
             m=0
             inp=0
             mv(1:9)=twofold%maxes(1:9)
             do_i:do i=1,twofold%ntwo-2
                vi=twofold%caxes(:,i)/twofold%maxes(i)
                do j=i+1,twofold%ntwo-1
                   vj=twofold%caxes(:,j)/twofold%maxes(j)
                   do k=j+1,twofold%ntwo
                      vk=twofold%caxes(:,k)/twofold%maxes(k)
                      aij=acosd(dot_product(vi,vj))
                      aik=acosd(dot_product(vi,vk))
                      ajk=acosd(dot_product(vj,vk))
                      if (abs(aij-90.0) < delt .and. abs(aik-90.0) < delt .and. abs(ajk-90.0) < delt ) then
                         if (abs(mv(i)-mv(j)) < tola .and. abs(mv(i)-mv(k)) < tola .and. abs(mv(j)-mv(k)) < tola ) then
                            rw(1)=i; rw(2)=j; rw(3)=k
                            u(1)=mv(i); u(2)=mv(j); u(3)=mv(k)
                            exit do_i
                         end if
                      end if
                   end do
                end do
             end do do_i

             do i=1,3
                tr(i,:) = twofold%dtwofold(:,rw(i))
                mv(i)=u(i)
             end do
             v1 = twofold%caxes(:,rw(1))
             v2 = twofold%caxes(:,rw(2))
             v3 = twofold%caxes(:,rw(3))
             namina=determ_A(tr)
             if (namina < 0) then
                tr(3,:)=-tr(3,:)
                v3=-v3
                namina=-namina
             end if

             Select Case (namina)
                Case(0)
                  write(unit=message,fmt="(a)") "Pseudo-cubic but tolerance too small ... "
                  ok=.false.
                  return
                Case(1)
                   message="Cubic, Primitive cell"
                Case(2)
                   rw=matmul((/0,1,1/),tr)
                   if (.not. co_prime(rw,3)) then
                      message="Cubic, A-centred cell"
                   else
                      rw=matmul((/1,1,1/),tr)
                      if (.not. co_prime(rw,3)) then
                         message="Cubic, I-centred cell"
                      else
                         rw=matmul((/1,1,0/),tr)
                         if (.not. co_prime(rw,3)) then
                            message="Cubic, C-centred cell"
                         else
                            rw=matmul((/1,0,1/),tr)
                            if (.not. co_prime(rw,3)) message="Cubic, B-centred cell"
                         end if
                      end if
                   end if

                Case(3:)
                  message="Cubic, F-centred cell"
             End Select

          case default
             write(unit=message,fmt="(a,i3)") "Wrong number of two-fold axes! ",twofold%ntwo
             ok=.false.
             return

      End Select

      !Calculation of the new cell
      ang(1)=acosd(dot_product(v2/mv(2),v3/mv(3)))
      ang(2)=acosd(dot_product(v1/mv(1),v3/mv(3)))
      ang(3)=acosd(dot_product(v1/mv(1),v2/mv(2)))
      Call Set_Crystal_Cell(mv(1:3),ang(1:3),Cell)
      ok=.true.

      return
    End Subroutine Get_Conventional_Cell

    !!----
    !!---- Subroutine Get_Cryst_Family(Cell,Car_Family,Car_Symbol,Car_System)
    !!----    type(Crystal_Cell_type),         intent(in ) :: Cell
    !!----    character(len=*),                intent(out) :: Car_Family
    !!----    character(len=*),                intent(out) :: Car_Symbol
    !!----    character(len=*),                intent(out) :: Car_System
    !!----
    !!---- Obtain the Crystal Family, Symbol and System from cell parameters
    !!----
    !!---- Update: May - 2005
    !!----
    Subroutine Get_Cryst_Family(Cell,Car_Family,Car_Symbol,Car_System)
       !---- Arguments ----!
       type(Crystal_Cell_type),   intent(in ) :: Cell
       character(len=*),          intent(out) :: Car_Family
       character(len=*),          intent(out) :: Car_Symbol
       character(len=*),          intent(out) :: Car_System

       !---- Local variables ----!
       integer, dimension(3) :: icodp, icoda
       integer               :: n1,n2

       Car_Family=" "
       Car_Symbol=" "
       Car_System=" "

       icodp=0
       icoda=0

       !---- Cell Parameters ----!

       !---- a ----!
       icodp(1)=1

       !---- b ----!
       if (abs(cell%cell(2)-cell%cell(1)) <= 0.0001) then
          icodp(2)=icodp(1)
       else
          icodp(2)=2
       end if

       !---- c ----!
       if (abs(cell%cell(3)-cell%cell(1)) <= 0.0001) then
          icodp(3)=icodp(1)
       else
          icodp(3)=3
       end if

       !---- Angles Parameters ----!

       !---- alpha ----!
       icoda(1)=1

       !---- beta ----!
       if (abs(cell%ang(2)-cell%ang(1)) <= 0.0001) then
          icoda(2)=icoda(1)
       else
          icoda(2)=2
       end if

       !---- gamma ----!
       if (abs(cell%ang(3)-cell%ang(1)) <= 0.0001) then
          icoda(3)=icoda(1)
       else
          icoda(3)=3
       end if


       n1=count(icoda==icoda(1))
       n2=count(icodp==icodp(1))
       select case (n1)
          case (1) ! all are differents
             if (n2 ==1) then
                Car_Family="Triclinic"
                Car_Symbol ="a"
                Car_System ="Triclinic"
             else
                Err_Crys=.true.
                ERR_Crys_Mess=" Error obtaining Crystal Familiy"
             end if

          case (2) ! two angles are equal
             if (icoda(1) == icoda(2)) then
                if (abs(cell%ang(3)-120.0) <= 0.0001) then
                   if (icodp(1)==icodp(2)) then
                      !---- Hexagonal ----!
                      Car_Family="Hexagonal"
                      Car_Symbol ="h"
                      Car_System ="Hexagonal"
                   else
                      Err_Crys=.true.
                      ERR_Crys_Mess=" Error obtaining Crystal Familiy"
                   end if
                else
                   !---- Monoclinic ----!
                   Car_Family="Monoclinic"
                   Car_Symbol ="m"
                   Car_System ="Monoclinic"
                end if

             else
                !---- Monoclic b-unique setting ----!
                if (abs(cell%ang(1)-90.0) <= 0.0001) then
                   Car_Family="Monoclinic"
                   Car_Symbol ="m"
                   Car_System ="Monoclinic"
                else
                   Err_Crys=.true.
                   ERR_Crys_Mess=" Error obtaining Crystal Familiy"
                end if
             end if

          case (3) ! all are the same angle
             if (abs(cell%ang(1) - 90.000) <= 0.0001) then
                select case (n2)
                   case (1)
                      !---- Orthorhombic ----!
                      Car_Family="Orthorhombic"
                      Car_Symbol ="o"
                      Car_System ="Orthorhombic"

                   case (2)
                      !---- Tetragonal ----!
                      if (icodp(1)==icodp(2)) then
                         Car_Family="Tetragonal"
                         Car_Symbol ="t"
                         Car_System ="Tetragonal"
                      else
                         err_crys=.true.
                         ERR_Crys_Mess=" Error obtaining Crystal Familiy"
                      end if

                   case (3)
                      !---- Cubic ----!
                      Car_Family="Cubic"
                      Car_Symbol ="c"
                      Car_System ="Cubic"
                end select

             else
                if (n2 == 3) then
                   !---- Hexagonal with rhombohedral axes ----!
                   Car_Family="Hexagonal"
                   Car_Symbol ="h"
                   Car_System ="Trigonal"
                else
                   Err_Crys=.true.
                   ERR_Crys_Mess=" Error obtaining Crystal Familiy"
                end if
             end if

       end select ! n

       return
    End Subroutine Get_Cryst_Family

    !!--++
    !!--++ Subroutine Get_Cryst_Orthog_Matrix(Cellv,Ang, Crystort,Cartype)
    !!--++    real(kind=cp), dimension(3  ), intent (in ) :: cellv           !  In ->  a,b,c parameters
    !!--++    real(kind=cp), dimension(3  ), intent (in ) :: ang             !  In ->  angles parameters of cell
    !!--++    real(kind=cp), dimension(3,3), intent (out) :: CrystOrt        ! Out ->  Conversion matrix (a) = (e) CrystOrt
    !!--++    character (len=1), optional,   intent (in)  :: CarType         !  In ->  Type of Cartesian axes
    !!--++
    !!--++    (PRIVATE)
    !!--++    Obtains the matrix giving the crystallographic basis in
    !!--++    direct space in terms of a Cartesian basis. The output matrix
    !!--++    can be directly used for transforming crystallographic components
    !!--++    to Cartesian components of the components of a vector considered
    !!--++    as a column vector:   XC = CrystOrt X.
    !!--++
    !!--++    If CartType is not present, or if it is not equal to 'A',
    !!--++    the cartesian system is defined as:
    !!--++          z // c; y is in the bc-plane; x is y ^ z
    !!--++    a = (a sinbeta singamma*, -a sinbeta cosgamma*, a cosbeta )
    !!--++    b = (         0         ,     b sinalpha      , b cosalpha)
    !!--++    c = (         0         ,         0           , c         )
    !!--++
    !!--++    If CartType = 'A', the Cartesian system is defined as:
    !!--++         x // a; y is in the ab-plane; z is x ^ z
    !!--++    a = (       a   ,         0           ,       0             )
    !!--++    b = ( b cosgamma,    b singamma       ,       0             )
    !!--++    c = (  c cosbeta, -c sinbeta cosalpha*, c sinbeta sinalpha* )
    !!--++
    !!--++    The output matrix is the tranposed of the above one(s) so that the
    !!--++    matrix can directly be used for transforming "components" given
    !!--++    in a crystallographic basis to "components" in cartesian basis
    !!--++    when the components are used as "column" vectors.
    !!--++
    !!--++      [a] = C [e] , In [a],[e] basis vectors are in column form
    !!--++      (a) = (e) CT, In (a),(e) basis vectors are in row form
    !!--++      CrystOrt = CT  => (a) = (e) CystOrt, in ITC: (a) = (e) P
    !!--++
    !!--++    Remember that  C.CT = GD (direct cell metrics)
    !!--++
    !!--++
    !!--++      Xc = CrystOrt X (Xc Cartesian components, X crystallographic components)
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Get_Cryst_Orthog_Matrix(Cellv,Ang, Crystort,CarType)
       !---- Arguments ----!
       real(kind=cp), dimension(3  ), intent (in ) :: cellv,ang
       real(kind=cp), dimension(3,3), intent (out) :: CrystOrt
       character (len=1), optional,   intent (in ) :: CarType

       !---- Local Variables ----!
       real(kind=cp) :: cosgas, singas

       if (present(CarType)) then
          if (CarType == "A" .or. CarType == "a" ) then  ! x//a
             !  Transponse of the following matrix:
             !    a = (       a   ,         0           ,       0             )
             !    b = ( b cosgamma,    b singamma       ,       0             )
             !    c = (  c cosbeta, -c sinbeta cosalpha*, c sinbeta sinalpha* )
             cosgas =(cosd(ang(3))*cosd(ang(2))-cosd(ang(1)))/(sind(ang(3))*sind(ang(2)))
             singas = sqrt(1.0-cosgas**2)
             CrystOrt(1,1) = cellv(1)
             CrystOrt(1,2) = cellv(2)*cosd(ang(3))
             CrystOrt(1,3) = cellv(3)*cosd(ang(2))
             CrystOrt(2,1) = 0.0
             CrystOrt(2,2) = cellv(2)*sind(ang(3))
             CrystOrt(2,3) =-cellv(3)*sind(ang(2))*cosgas
             CrystOrt(3,1) = 0.0
             CrystOrt(3,2) = 0.0
             CrystOrt(3,3) = cellv(3)*sind(ang(2))*singas
             return
          end if
       end if

       !
       !  By default, the cartesian frame is such as z//c
       !  Transponse of the following matrix:
       !    a = (a sinbeta singamma*, -a sinbeta cosgamma*, a cosbeta )
       !    b = (         0         ,     b sinalpha      , b cosalpha)
       !    c = (         0         ,         0           , c         )
       cosgas =(cosd(ang(1))*cosd(ang(2))-cosd(ang(3)))/(sind(ang(1))*sind(ang(2)))
       singas = sqrt(1.0-cosgas**2)
       CrystOrt(1,1) = cellv(1)*sind(ang(2))*singas
       CrystOrt(1,2) = 0.0
       CrystOrt(1,3) = 0.0
       CrystOrt(2,1) =-cellv(1)*sind(ang(2))*cosgas
       CrystOrt(2,2) = cellv(2)*sind(ang(1))
       CrystOrt(2,3) = 0.0
       CrystOrt(3,1) = cellv(1)*cosd(ang(2))
       CrystOrt(3,2) = cellv(2)*cosd(ang(1))
       CrystOrt(3,3) = cellv(3)

       return
    End Subroutine Get_Cryst_Orthog_Matrix

    !!----
    !!---- Subroutine Get_Deriv_Orth_Cell(Cellp,De_Orthcell,Cartype)
    !!----    type(Crystal_Cell_type),         intent(in ) :: cellp
    !!----    real(kind=cp), dimension(3,3,6), intent(out) :: de_Orthcell
    !!----    character (len=1), optional,     intent(in ) :: CarType
    !!----
    !!----    Subroutine to get derivative matrix of the transformation matrix
    !!----    to orthogonal frame. Useful for calculations of standard deviations
    !!----    of distances and angles. The specialised subroutine calculating
    !!----    sigmas of distances "distance_and_sigma" is in Atom_mod.
    !!----    The output matrices "de_Orthcell" are the derivatives of, with
    !!----    respect to a(1),b(2),c(3),alpha(4),beta(5) and gamma(6) of the
    !!----    matrix   "Cellp%Cr_Orth_cel".
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Get_Deriv_Orth_Cell(Cellp,De_Orthcell,Cartype)
       !---- Arguments ----!
       type(Crystal_Cell_type),         intent(in ) :: cellp
       real(kind=cp), dimension(3,3,6), intent(out) :: de_Orthcell
       character (len=1), optional,     intent(in ) :: CarType

       !---- Local Variables ----!
       real(kind=cp) ::  ca,cb,cg,sa,sb,sg,f,g, fa,fb,fc,ga,gb,gc

       de_Orthcell=0.0
       ca=cosd(cellp%ang(1))
       cb=cosd(cellp%ang(2))
       cg=cosd(cellp%ang(3))
       sa=sind(cellp%ang(1))
       sb=sind(cellp%ang(2))
       sg=sind(cellp%ang(3))

       if (present(CarType)) then
          if (CarType == "A" .or. CarType == "a" ) then  ! x//a

             f=(ca-cb*cg)/sg    !-cosgas*sinbeta
             g=SQRT(sb*sb-f*f)  ! singas*sinbeta
             fa=-sa/sg          ! df/dalpha
             fb=sb*cg/sg        ! df/dbeta
             fc=cb/sg**2        ! df/dgamma
             ga=-f*fa/g         ! dg/dalpha
             gb=(sb*cb-f*fb)/g  ! dg/dbeta
             gc=f/g*fc          ! dg/dgamma

             ! M: Transponse of the following matrix:
             !    a = (       a   ,         0           ,       0             )
             !    b = ( b cosgamma,    b singamma       ,       0             )
             !    c = (  c cosbeta, -c sinbeta cosalpha*, c sinbeta sinalpha* )

             !
             !        (   a         b*cg        c*cb )
             !    M = (   0         b*sg        c*f  )
             !        (   0          0          c*g  )
             !
             !           (   1      0      0 )
             !  dM_da =  (   0      0      0 )
             !           (   0      0      0 )
             de_Orthcell(1,1,1) = 1.0

             !           (   0      cg     0 )
             !  dM_db =  (   0      sg     0 )
             !           (   0      0      0 )
             de_Orthcell(1,2,2) = cg
             de_Orthcell(2,2,2) = sg

             !
             !            (   0          0          cb )
             !  dM_dc =   (   0          0          f  )
             !            (   0          0          g  )
             de_Orthcell(1,3,3) = cb
             de_Orthcell(2,3,3) = f
             de_Orthcell(3,3,3) = g

             !
             !             (   0          0           0   )
             ! dM_dalpha=  (   0          0          c*fa )
             !             (   0          0          c*ga )
             !
             de_Orthcell(2,3,4) = cellp%cell(3)*fa
             de_Orthcell(3,3,4) = cellp%cell(3)*ga

             !
             !             (   0          0         -c*sb )
             ! dM_dbeta =  (   0          0          c*fb )
             !             (   0          0          c*gb )
             !
             de_Orthcell(1,3,5) = -cellp%cell(3)*sb
             de_Orthcell(2,3,5) =  cellp%cell(3)*fb
             de_Orthcell(3,3,5) =  cellp%cell(3)*gb

             !
             !              (   0        -b*sg         0   )
             ! dM_dgamma =  (   0         b*cg        c*fc )
             !              (   0          0          c*gc )
             !
             de_Orthcell(1,2,6) = -cellp%cell(2)*sg
             de_Orthcell(2,2,6) =  cellp%cell(2)*cg
             de_Orthcell(2,3,6) =  cellp%cell(3)*fc
             de_Orthcell(3,3,6) =  cellp%cell(3)*gc

             return
          end if
       end if

       !
       !  By default, the cartesian frame is such as z//c
       !  Transponse of the following matrix:
       !    a = (a sinbeta singamma*, -a sinbeta cosgamma*, a cosbeta )
       !    b = (         0         ,     b sinalpha      , b cosalpha)
       !    c = (         0         ,         0           , c         )

       !         ( a sinbeta singamma*          0             0 )
       !    M =  (-a sinbeta cosgamma*      b sinalpha        0 )
       !         ( a cosbeta                b cosalpha        c )

       f=(cg-ca*cb)/sa    !-sinbeta . cosgamma*
       g=SQRT(sb*sb-f*f)  ! sinbeta . singamma*
       fa= cb/sa**2       ! df/dalpha
       fb=sb*ca/sa        ! df/dbeta
       fc=-sb/sa          ! df/dgamma
       ga=-f*fa/g         ! dg/dalpha
       gb=(sb*cb-f*fb)/g  ! dg/dbeta
       gc=f/g*fc          ! dg/dgamma

       !         ( a*g        0         0 )
       !    M =  ( a*f      b*sa        0 )
       !         ( a*cb     b*ca        c )

       !
       !           (   g       0      0 )
       !  dM_da =  (   f       0      0 )
       !           (   cb      0      0 )
       de_Orthcell(1,1,1) = g
       de_Orthcell(1,2,1) = f
       de_Orthcell(1,3,1) = cb

       !           (   0      0      0 )
       !  dM_db =  (   0      sa     0 )
       !           (   0      ca     0 )
       de_Orthcell(1,2,2) = sa
       de_Orthcell(3,2,2) = ca

       !
       !            (   0      0      0  )
       !  dM_dc =   (   0      0      0  )
       !            (   0      0      1  )
       de_Orthcell(3,3,3) = 1

       !
       !             ( a*ga         0          0 )
       ! dM_dalpha=  ( a*fa       -b*ca        0 )
       !             (   0         b*sa        0 )
       !
       de_Orthcell(1,1,4) = cellp%cell(1)*ga
       de_Orthcell(2,1,4) = cellp%cell(1)*fa
       de_Orthcell(2,2,4) =-cellp%cell(2)*ca
       de_Orthcell(3,2,4) = cellp%cell(2)*sa

       !
       !             (  a*gb        0         0 )
       ! dM_dbeta =  (  a*fb        0         0 )
       !             ( -a*sb        0         0 )
       !
       de_Orthcell(1,1,5) = cellp%cell(1)*gb
       de_Orthcell(2,1,5) = cellp%cell(1)*fb
       de_Orthcell(3,1,5) =-cellp%cell(1)*sb

       !
       !              (  a*gc     0      0   )
       ! dM_dgamma =  (  a*fc     0      0   )
       !              (   0       0      0   )
       !
       de_Orthcell(1,1,6) = cellp%cell(1)*gc
       de_Orthcell(2,1,6) = cellp%cell(1)*fc

       return
    End Subroutine Get_Deriv_Orth_Cell

    !!----
    !!---- Subroutine Get_Primitive_Cell(Lat_Type,Centred_Cell,Primitive_Cell,Transfm)
    !!----    character(len=*),               intent(in)  :: lat_type
    !!----    type(Crystal_Cell_Type),        intent(in)  :: centred_cell
    !!----    type(Crystal_Cell_Type),        intent(out) :: primitive_cell
    !!----    real(kind=cp), dimension(3,3),  intent(out) :: transfm
    !!----
    !!----    Subroutine for getting the primitive cell from a centred cell
    !!----    On input Lat_type is the lattice type: P,A,B,C,I,R or F
    !!----    Centred_cell is the Crystal_Cell_Type of the input lattice
    !!----    The subroutine calculates the transformation matric "transfm"
    !!----    and provides the complete description of the primitive cell
    !!----    in the object, of type Crystal_Cell_Type, primitive_cell.
    !!----
    !!---- Update: April - 2008
    !!
    Subroutine Get_Primitive_Cell(Lat_Type,Centred_Cell,Primitive_Cell,Transfm)
       !---- Arguments ----!
       character(len=*),              intent(in)  :: lat_type
       type(Crystal_Cell_Type),       intent(in)  :: centred_cell
       type(Crystal_Cell_Type),       intent(out) :: primitive_cell
       real(kind=cp), dimension(3,3), intent(out) :: transfm

       !---- Local variables ----!
       integer                       :: i
       real(kind=cp), dimension(3)   :: celp,celang
       real(kind=cp), dimension(3,3) :: cart,metric
       character(len=1)              :: lat

       lat=adjustl(lat_type)
       Select Case(lat)
          case("a","A")
             transfm= reshape((/1.0,0.0,0.0,  0.0,0.5,0.5,  0.0,-0.5,0.5/),(/3,3/))
          case("b","B")
             transfm= reshape((/0.5,0.0,0.5,  0.0,1.0,0.0, -0.5, 0.0,0.5/),(/3,3/))
          case("c","C")
             transfm= reshape((/0.5,0.5,0.0, -0.5,0.5,0.0,  0.0, 0.0,1.0/),(/3,3/))
          case("i","I")
             transfm= reshape((/1.0,0.0,0.0,  0.0,1.0,0.0,  0.5, 0.5,0.5/),(/3,3/))
          case("r","R")
             transfm= reshape((/2.0/3.0, 1.0/3.0, 1.0/3.0,  &
                               -1.0/3.0, 1.0/3.0, 1.0/3.0,  &
                               -1.0/3.0,-2.0/3.0, 1.0/3.0/),(/3,3/))
          case("f","F")
             transfm= reshape((/0.5,0.0,0.5,  0.5,0.5,0.0,  0.0, 0.5,0.5/),(/3,3/))
          case default  !assumed primitive
             primitive_cell=centred_cell
             transfm= reshape((/1.0,0.0,0.0,  0.0,1.0,0.0,  0.0,0.0,1.0/),(/3,3/))
             return
       End Select
       transfm=transpose(transfm)
       cart=matmul(transfm,transpose(Centred_Cell%Cr_Orth_cel))
       metric=matmul(cart,transpose(cart))

       !---- Calculate new cell parameters from the new metric tensor
       do i=1,3
          Celp(i)=sqrt(metric(i,i))
       end do

       celang(1)=acosd(metric(2,3)/(celp(2)*celp(3)))
       celang(2)=acosd(metric(1,3)/(celp(1)*celp(3)))
       celang(3)=acosd(metric(1,2)/(celp(1)*celp(2)))
       call Set_Crystal_Cell(celp,celang,primitive_cell)

       return
    End Subroutine Get_Primitive_Cell

    !!----
    !!---- Subroutine Get_Transfm_Matrix(cella,cellb,trm,ok,tol)
    !!----    type(Crystal_Cell_Type),     intent(in) :: cella,cellb
    !!----    real(kind=cp),dimension(3,3),intent(out):: trm
    !!----    Logical,                     intent(out):: ok
    !!----    real(kind=cp),optional,      intent(in) :: tol
    !!----
    !!----    Subroutine for getting the transformation matrix between two
    !!----    primitive unit cells (the range of indices is fixed to -2 to 2)
    !!----
    !!---- Update: January - 2011
    !!
    Subroutine Get_Transfm_Matrix(cella,cellb,trm,ok,tol)
       !---- Arguments ----!
       type(Crystal_Cell_Type),     intent(in) :: cella,cellb
       real(kind=cp),dimension(3,3),intent(out):: trm
       Logical,                     intent(out):: ok
       real(kind=cp),optional,      intent(in) :: tol

       !---- Local variables ----!
       type(Crystal_Cell_Type) :: Cellt
       integer,dimension(3,3)  :: Nu
       integer                 :: j,i1,i2,i3,i4,i5,i6,i7,i8,i9
       real(kind=cp)           :: tolt

       tolt=0.3
       if(present(tol)) tolt=tol
       ok=.false.
       dox: do i1=-2,2                     !         |i1  i4  i7|
          do i2=-2,2                       !    Nu = |i2  i5  i8|
             do i3=-2,2                    !         |i3  i6  i9|
                do i4=-2,2
                   do i5=-2,2
                      do i6=-2,2
                         do i7=-2,2
                            do i8=-2,2
                               do i9=-2,2
                                  j=i1*i5*i9+i4*i8*i3+i2*i6*i7-i3*i5*i7-i8*i6*i1-i2*i4*i9     !determinant (much faster than calling determ_A)
                                  if ( j /= 1) cycle
                                  Nu=reshape((/i1,i2,i3,i4,i5,i6,i7,i8,i9/),(/3,3/))
                                  Trm=real(Nu)
                                  call Change_Setting_Cell(Cella,Trm,Cellt)
                                  if (Sum(abs(Cellt%cell(:)-Cellb%cell(:)))+Sum(abs(Cellt%ang(:)-Cellb%ang(:))) < tolt  ) then
                                     ok=.true.
                                     exit dox
                                  end if
                               end do    !i9
                            end do     !i8
                         end do      !i7
                      end do       !i6
                   end do        !i5
                end do         !i4
             end do          !i3
          end do           !i2
       end do  dox       !i1

       return
    End Subroutine Get_Transfm_Matrix

    !!----
    !!---- Subroutine Get_TwoFold_Axes(Celln,Tol,Twofold)
    !!----    type(Crystal_Cell_Type), intent (in) :: Celln
    !!----    real(kind=cp),           intent (in) :: tol !angular tolerance in degrees
    !!----    Type(Twofold_Axes_Type), intent(out) :: twofold
    !!----
    !!----    Subroutine for getting the possible two-fold axes (within an
    !!----    angular tolerance tol) existing in the lattice generated by the
    !!----    unit cell "Celln". Strictly independent two-fold axes are stored
    !!----    in the variable "twofold" that is of type Twofold_Axes_Type
    !!----    The output order of the two-fold axes is ascending in their
    !!----    modulus. Shorter vectors appears before longer ones.
    !!----    The conditions for a reciprocal or direct row to be a two-fold
    !!----    axis are discussed by Y. Le Page in J.Appl.Cryst. 15, 255 (1982).
    !!----
    !!----
    !!---- Update: November - 2008
    !!
    Subroutine Get_TwoFold_Axes(Celln,Tol,Twofold)
       !---- Arguments ----!
       type(Crystal_Cell_Type), intent (in) :: Celln
       real(kind=cp),           intent (in) :: Tol !angular tolerance in degrees
       Type(twofold_axes_type), intent(out) :: Twofold

       !---- Local variables ----!
       integer                        :: i,j,n,m, ih,ik,il,iu,iv,iw,imax,ntwo
       real(kind=cp), dimension(3)    :: dv, rv, a, b, c, as, bs, cs, cross
       real(kind=cp), dimension(  12) :: maxes,crossa
       integer, dimension(  12)       :: dota,ind
       real(kind=cp), dimension(3,12) :: caxes
       integer, dimension(3,12)       :: dtw,rtw
       integer, dimension(3)          :: v,h
       real(kind=cp)                  :: dot,crossm

       maxes=0.0; crossa=0.0; dota=0; caxes=0.0; dtw=0; rtw=0
       a=Celln%Cr_Orth_cel(:,1)
       b=Celln%Cr_Orth_cel(:,2)
       c=Celln%Cr_Orth_cel(:,3)
       twofold%a=a
       twofold%b=b
       twofold%c=c
       as=cross_product(b,c)/Celln%CellVol !Reciprocal lattice vectors in
       bs=cross_product(c,a)/Celln%CellVol !Cartesian components
       cs=cross_product(a,b)/Celln%CellVol
       ntwo=0
       imax=2   !Is inough if the input cell is the Buerger or Niggli cell

       do_iu: do iu=imax, 0,-1
          do iv=imax,-imax,-1
             do iw=imax,-imax,-1
                v=(/iu,iv,iw/)
                if (.not. Co_Prime(v,2)) cycle
                do ih=imax,0,-1
                   do ik=imax,-imax,-1
                      do_il:do il=imax,-imax,-1
                         h=(/ih,ik,il/)
                         if (.not. Co_Prime(h,2)) cycle
                         n=abs(ih*iu+ik*iv+il*iw)
                         if ( n == 2 .or. n == 1) then
                            dv=real(iu)*a+real(iv)*b+real(iw)*c
                            rv=real(ih)*as+real(ik)*bs+real(il)*cs
                            cross=cross_product(dv,rv)
                            dot=sqrt(dot_product(cross,cross))
                            crossm=atand(dot/real(n))
                            if (abs(crossm) <= tol) then
                               do m=1,ntwo
                                  if (determ_V((/17,41,71/),v,dtw(:,m) ) == 0) cycle do_il
                               end do
                               ntwo=ntwo+1
                               dtw(:,ntwo)= v
                               dv=v(1)*a+v(2)*b+v(3)*c
                               caxes(:,ntwo)=dv
                               maxes(ntwo)=sqrt(dot_product(dv,dv))
                               rtw(:,ntwo)= h
                               dota(ntwo)=n
                               crossa(ntwo)=crossm
                            end if
                            if (ntwo == 12) exit do_iu
                         end if
                      end do do_il
                   end do
                end do
             end do
          end do
       end do do_iu
       call sort(maxes,ntwo,ind)
       do i=1,ntwo
          j=ind(i)
          twofold%dtwofold(:,i)= dtw(:,j)
          twofold%caxes(:,i)= caxes(:,j)
          twofold%maxes(i)= maxes(j)
          twofold%rtwofold(:,i)= rtw(:,j)
          twofold%dot(i)= dota(j)
          twofold%cross(i)= crossa(j)
       End do
       twofold%ntwo=ntwo
       twofold%tol=tol

       return
    End Subroutine Get_TwoFold_Axes

    !!----
    !!---- SUBROUTINE INIT_ERR_CRYS()
    !!----
    !!----    Initialize Flags of Errors in this module
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Init_Err_Crys()

       Err_Crys=.false.
       ERR_Crys_Mess=" "

       return
    End Subroutine Init_Err_Crys

    !!----
    !!---- Subroutine Niggli_Cell(XXX,Niggli_Point,Celln,Trans)
    !!----   XXX is one of:
    !!----   real(kind=cp),dimension(6),              intent(in out) :: Ad             ! Cell Parameters
    !!----   or
    !!----   real(kind=cp),dimension(2,3),            intent(in out) :: N_Mat          ! Niggli Matrix
    !!----   or
    !!----   real(kind=cp)                            intent(in out) :: A, B, C, Alfa, Beta, Gamma
    !!----   or
    !!----   type(Crystal_Cell_Type),                 intent(in out ):: cell
    !!----   or
    !!----   real(kind=cp),dimension(3),              intent(in)     :: A,B,C         ! 3 vectors
    !!----   real(kind=cp),dimension(5), optional,    intent(out)    :: Niggli_Point
    !!----   type(Crystal_Cell_Type),optional,        intent(out)    :: Celln
    !!----   real(kind=cp), dimension(3,3), optional, intent(out)    :: Trans
    !!----
    !!----    Calculates the Niggli cell
    !!----
    !!---- Update: October - 2008
    !!

    !!--++
    !!--++ Subroutine Niggli_Cell_ABC(Ad,Niggli_Point,Celln,Trans)
    !!--++    real(kind=cp),dimension(6),              intent(in out) :: Ad
    !!--++    real(kind=cp),dimension(5), optional,    intent(out)    :: Niggli_Point
    !!--++    type(Crystal_Cell_Type),optional,        intent(out)    :: celln
    !!--++    real(kind=cp), dimension(3,3), optional, intent(out)    :: trans
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the Niggli cell when the input is the list of cell parameters
    !!--++    provided as a 6D vector. Calls the subroutine Niggli_Cell_Nigglimat for
    !!--++    the effective calculations
    !!--++
    !!--++ Update: October - 2008
    !!
    Subroutine Niggli_Cell_ABC(Ad,Niggli_Point,Celln,Trans)    !Scalar algorithm
       !---- Arguments ----!
       real(kind=cp),dimension(6),              intent(in out) :: ad
       real(kind=cp),dimension(5), optional,    intent(out)    :: Niggli_Point
       type(Crystal_Cell_Type),optional,        intent(out)    :: celln
       real(kind=cp), dimension(3,3), optional, intent(out)    :: trans

       !---- Local variables ----!
       real(kind=cp), dimension(2,3)    :: n_mat
       type(Crystal_Cell_Type)          :: celda

       n_mat(1,1)=ad(1)*ad(1)
       n_mat(1,2)=ad(2)*ad(2)
       n_mat(1,3)=ad(3)*ad(3)
       n_mat(2,1)=ad(2)*ad(3)*cosd(ad(4))
       n_mat(2,2)=ad(1)*ad(3)*cosd(ad(5))
       n_mat(2,3)=ad(1)*ad(2)*cosd(ad(6))

       if (present(Niggli_Point)) then
          if (present(trans)) then
             call Niggli_Cell_nigglimat(n_mat,Niggli_Point,celda,trans)
          else
             call Niggli_Cell_nigglimat(n_mat,Niggli_Point,celda)
          end if
       else if(present(trans)) then
          call Niggli_Cell_nigglimat(n_mat,celln=celda,trans=trans)
       else
          call Niggli_Cell_nigglimat(n_mat,celln=celda)
       end if

       if (Err_Crys) return
       if (present(celln)) celln=celda

       !Reconstruct the new cell (Niggli Cell)
       ad(1) = sqrt(n_mat(1,1))
       ad(2) = sqrt(n_mat(1,2))
       ad(3) = sqrt(n_mat(1,3))
       ad(4) = acosd(n_mat(2,1)/(ad(2)*ad(3)))
       ad(5) = acosd(n_mat(2,2)/(ad(1)*ad(3)))
       ad(6) = acosd(n_mat(2,3)/(ad(1)*ad(2)))

       return
    End Subroutine Niggli_Cell_abc

    !!--++
    !!--++ Subroutine Niggli_Cell_Nigglimat(N_Mat,Niggli_Point,Celln,Trans)    !Scalar algorithm
    !!--++    real(kind=cp),dimension(2,3),              intent(in out) :: n_mat
    !!--++    real(kind=cp),dimension(5),      optional, intent(out)    :: Niggli_Point
    !!--++    type(Crystal_Cell_Type), optional,         intent(out)    :: celln
    !!--++    real(kind=cp), dimension(3,3),   optional, intent(out)    :: trans
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the Niggli cell when the input is the Niggli Matrix (part of the metrics)
    !!--++    of a primitive cell. Applies the scalar algorithm of
    !!--++    I. Krivy and B. Gruber, Acta Cryst A32, 297 (1976)
    !!--++    If Trans is present, Celln should also be present.
    !!--++
    !!--++ Update: January - 2011
    !!
    Subroutine Niggli_Cell_Nigglimat(N_Mat,Niggli_Point,Celln,Trans)    !Scalar algorithm
       !---- Arguments ----!
       real(kind=cp),dimension(2,3),              intent(in out) :: n_mat
       real(kind=cp),dimension(5),      optional, intent(out)    :: Niggli_Point
       type(Crystal_Cell_Type),         optional, intent(out)    :: celln
       real(kind=cp), dimension(3,3),   optional, intent(out)    :: trans

       !--- Local variables ---!
       type(Crystal_Cell_Type)       :: Cellp
       real(kind=cp)                 :: A,B,C,u,v,w,eps
       real(kind=cp), dimension(3,3) :: trm
       real(kind=cp), dimension(3)   :: cel,ang
       integer                       :: iu,iv,iw, ncount ! ncount is the counter no more that Numiter=100
                                                         ! iterations are permitted. In case of exhausting
                                                         ! the iteration Err_Crys=.true. but the current
                                                         ! cell is output anyway
       real(kind=cp),parameter        :: epr=0.0001      !Relative epsilon
       integer, parameter             :: numiter=100
       logical                        :: ok

       ! N is a Niggli cell of L if  (i) it is as Buerger cell of L and
       !                            (ii) |90-alpha| + |90-beta| + |90-gamma| -> maximum
       !                  / a.a  b.b  c.c \       /  s11  s22  s33 \
       !   Niggli matrix  |               |   =   |                |
       !                  \ b.c  a.c  a.b /       \  s23  s13  s12 /
       !
       ! I. Krivy and B. Gruber, Acta Cryst A32, 297 (1976)
       ! Krivy-Gruber algorithms safely implemented (suggestion of Ralf Grosse-Kunsleve)
       ! R.W. Grosse-Kunstleve, N. K. Sauter and P. D. Adams, Acta Cryst A60, 1-6 (2004)
       ! Epsilon: e
       !    x < y -> x < y-e;    x > y -> y < x-e
       !   x <= y -> .not. y < x-e;   x >= y -> .not. x < y-e
       !   x == y -> .not. (x < y-e .or. y < x-e)
       !
       A=n_mat(1,1)
       B=n_mat(1,2)
       C=n_mat(1,3)
       u=2.0*n_mat(2,1)
       v=2.0*n_mat(2,2)
       w=2.0*n_mat(2,3)
       eps=epr*(A*B*C)**(1.0/6.0)
       ncount=0
       ok=.true.
       if (present(trans)) then
          !Construct the input cell Cellp from its Niggli parameters
          cel(1) = sqrt(A)
          cel(2) = sqrt(B)
          cel(3) = sqrt(C)
          ang(1) = acosd(u/(cel(2)*cel(3)*2.0))
          ang(2) = acosd(v/(cel(1)*cel(3)*2.0))
          ang(3) = acosd(w/(cel(1)*cel(2)*2.0))
          call Set_Crystal_Cell(cel,ang, Cellp)
       end if

       do
          ncount=ncount+1
          if (ncount > numiter) then
             ok=.false.
             exit
          end if

          !---- if(A > B .or. ( A == B  .and. abs(u) > abs(v)) ) then  ! A1
          if (B < A-eps .or. ( .not.( A < B-eps .or. B < A-eps)  .and. abs(v) < abs(u)-eps ) ) then  ! A1
             call swap(A,B)
             call swap(u,v)
          end if

          !---- if(B > C .or. ( B == C .and. abs(v) > abs(w)) ) then  ! A2
          if (C < B-eps .or. ( .not.( C < B-eps .or. B < C-eps) .and. abs(w) < abs(v)-eps) ) then  ! A2
             call swap(B,C)
             call swap(v,w)
             cycle
          end if

          !---- if (u*v*w > 0.0) then                                 ! A3
          iu=1; iv=1; iw=1
          if ( u < -eps) iu=-1
          if ( v < -eps) iv=-1
          if ( w < -eps) iw=-1
          if (abs(u) < eps) iu=0
          if (abs(v) < eps) iv=0
          if (abs(w) < eps) iw=0
          if (iu*iv*iw > 0) then                                      ! A3
             u=abs(u)
             v=abs(v)
             w=abs(w)
          else                                                        ! A4
             u=-abs(u)
             v=-abs(v)
             w=-abs(w)
          end if

          !---- if( abs(u) > B .or. ( u == B .and. 2.0*v < w) .or. ( u == -B .and. w < 0.0)) then  ! A5
          if ( B < abs(u)-eps  .or. ( .not.(u < B-eps .or. B < u-eps) .and. 2.0*v < w-eps) .or. &
             ( .not.(u < -B-eps .or. -B < u-eps) .and. w < -eps)) then  ! A5
             iu=1; if( u < -eps) iu=-1
             C = B+C - u * iu
             v =  v  - w * iu
             u = u - 2.0*B*iu
             cycle
          end if

          !---- if( abs(v) > A .or. ( v == A .and. 2.0*u < w) .or. ( v == -A .and. w < 0.0)) then  ! A6
          if ( A < abs(v)-eps .or. (.not. (v < A-eps .or. A < v-eps) .and. 2.0*u < w-eps) .or. &
             ( .not.( v < -A-eps .or. -A < v-eps) .and. w < -eps)) then  ! A6
             iv=1; if( v < -eps) iv=-1
             C = A+C - v * iv
             u =  u  - w * iv
             v = v - 2.0*A*iv
             cycle
          end if

          !---- if( abs(w) > A .or. ( w == A .and. 2.0*u < v) .or. ( w == -A .and. v < 0.0)) then  ! A7
          if ( A < abs(w)-eps .or. ( .not. (w < A-eps .or. A < w-eps) .and. 2.0*u < v-eps) .or. &
             ( .not. (w < -A-eps .or. -A < w-eps) .and. v < -eps)) then  ! A7
             iw=1; if( w < -eps) iw=-1
             B = A+B - w * iw
             u =  u  - v * iw
             w = w - 2.0*A*iw
             cycle
          end if

          !---- if(u+v+w+A+B < 0.0 .or. (u+v+w+A+B == 0.0 .and. 2.0*(A+v)+w > 0.0 )) then  ! A8
          if (u+v+w+A+B < -eps .or. ( abs(u+v+w+A+B) < eps .and. 2.0*(A+v)+w > eps )) then  ! A8
             C=A+B+C+u+v+w
             u=2.0*B+u+w
             v=2.0*A+v+w
             cycle
          end if
          exit
       end do

       !---- Reconstruct the new Niggli matrix
       n_mat(1,1)=A; n_mat(1,2)=B; n_mat(1,3)=C
       n_mat(2,1)=0.5*u; n_mat(2,2)=0.5*v; n_mat(2,3)=0.5*w

       if (.not. ok) Then
          Err_Crys=.true.
          ERR_Crys_Mess=" The limit of iterations in Niggli_Cell_NiggliMat has been reached!"
          return
       end if

       if (present(Niggli_point)) then
          Niggli_point(1)= A/C
          Niggli_point(2)= B/C
          Niggli_point(3)= u/C
          Niggli_point(4)= v/C
          Niggli_point(5)= w/C
       end if

       if (present(celln)) then
          !Reconstruct the new cell (Niggli Cell)
          cel(1) = sqrt(A)
          cel(2) = sqrt(B)
          cel(3) = sqrt(C)
          ang(1) = acosd(u/(cel(2)*cel(3)*2.0))
          ang(2) = acosd(v/(cel(1)*cel(3)*2.0))
          ang(3) = acosd(w/(cel(1)*cel(2)*2.0))
          call Set_Crystal_Cell(cel,ang, Celln)
          if (present(trans)) then
            Call Get_Transfm_Matrix(cellp,celln,trm,ok)
            if(ok) then
              trans=trm
            else
              trans=identity
            end if
          end if
       end if

       return
    End Subroutine Niggli_Cell_nigglimat

    !!--++
    !!--++ Subroutine Niggli_Cell_Params(A,B,C,Al,Be,Ga,Niggli_Point,Celln,Trans)
    !!--++    real(kind=cp),                           intent (in out)  :: a,b,c,al,be,ga
    !!--++    real(kind=cp),dimension(5), optional,    intent(out)      :: Niggli_Point
    !!--++    type(Crystal_Cell_Type),optional,        intent(out)      :: celln
    !!--++    real(kind=cp), dimension(3,3), optional, intent(out)      :: trans
    !!--++
    !!--++    (OVERLOAD)
    !!--++     Calculates the Niggli cell when the input is the list of cell parameters
    !!--++     provided as six scalars.
    !!--++     Calls the subroutine Niggli_Cell_Nigglimat for the effective calculations
    !!--++
    !!--++ Update: October - 2008
    !!
    Subroutine Niggli_Cell_Params(A,B,C,Al,Be,Ga,Niggli_Point,Celln,Trans)
       !---- Arguments ----!
       real(kind=cp),                           intent (in out)  :: a,b,c,al,be,ga
       real(kind=cp),dimension(5), optional,    intent(out)      :: Niggli_Point
       type(Crystal_Cell_Type), optional,       intent(out)      :: celln
       real(kind=cp), dimension(3,3), optional, intent(out)      :: trans

       !--- Local variables ---!
       type(Crystal_Cell_Type)          :: celda
       real(kind=cp), dimension(2,3)    :: n_mat


       call Init_Err_Crys()
       if ( al+be < ga+1.0  .or. al+ga < be+1.0 .or. be+ga < al+1.0) then
          Err_Crys=.true.
          ERR_Crys_Mess=" The provided angles cannot set a unit cell!"
          return
       end if

       call Set_Crystal_Cell((/a,b,c/),(/al,be,ga/), Celda)
       if (Err_Crys) return

       n_mat(1,1)=Celda%GD(1,1); n_mat(1,2)=Celda%GD(2,2); n_mat(1,3)=Celda%GD(3,3)
       n_mat(2,1)=Celda%GD(2,3); n_mat(2,2)=Celda%GD(1,3); n_mat(2,3)=Celda%GD(1,2)

       if (present(Niggli_Point)) then
          if (present(trans)) then
             call Niggli_Cell_nigglimat(n_mat,Niggli_Point,celda,trans)
          else
             call Niggli_Cell_nigglimat(n_mat,Niggli_Point,celda)
          end if
       else if(present(trans)) then
          call Niggli_Cell_nigglimat(n_mat,celln=celda,trans=trans)
       else
          call Niggli_Cell_nigglimat(n_mat,celln=celda)
       end if
       if (Err_Crys) return

       if (present(celln)) then
          celln=celda
       else
           a=celda%cell(1); b=celda%cell(2); c=celda%cell(3)
          al=celda%ang(1); be=celda%ang(2); ga=celda%ang(3)
       end if

       return
    End Subroutine Niggli_Cell_Params

    !!--++
    !!--++ Subroutine Niggli_Cell_Type(Cell,Niggli_Point,Celln,Trans)
    !!--++    type(Crystal_Cell_Type),                 intent(in out ) :: cell
    !!--++    real(kind=cp),dimension(5),    optional, intent(out)     :: Niggli_Point
    !!--++    type(Crystal_Cell_Type),       optional, intent(out)     :: celln
    !!--++    real(kind=cp), dimension(3,3), optional, intent(out)     :: trans
    !!--++
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the Niggli cell when the input is an object of type Crystal_Cell_Type
    !!--++    Calls the subroutine Niggli_Cell_Nigglimat for the effective calculations
    !!--++
    !!--++ Update: October - 2008
    !!
    Subroutine Niggli_Cell_Type(Cell,Niggli_Point,Celln,Trans)
       !---- Arguments ----!
       type(Crystal_Cell_Type),                 intent(in out ) :: cell
       real(kind=cp),dimension(5),    optional, intent(out)     :: Niggli_Point
       type(Crystal_Cell_Type),       optional, intent(out)     :: celln
       real(kind=cp), dimension(3,3), optional, intent(out)     :: trans

       !--- Local variables ---!
       type(Crystal_Cell_Type)         :: celda
       real(kind=cp), dimension(2,3)   :: n_mat

       call Init_Err_Crys()
       celda=cell
       n_mat(1,1)=Celda%GD(1,1); n_mat(1,2)=Celda%GD(2,2); n_mat(1,3)=Celda%GD(3,3)
       n_mat(2,1)=Celda%GD(2,3); n_mat(2,2)=Celda%GD(1,3); n_mat(2,3)=Celda%GD(1,2)

       if (present(Niggli_Point)) then
          if (present(trans)) then
             call Niggli_Cell_nigglimat(n_mat,Niggli_Point,celda,trans)
          else
             call Niggli_Cell_nigglimat(n_mat,Niggli_Point,celda)
          end if
       else if(present(trans)) then
          call Niggli_Cell_nigglimat(n_mat,celln=celda,trans=trans)
       else
          call Niggli_Cell_nigglimat(n_mat,celln=celda)
       end if
       if (Err_Crys) return

       if (present(celln)) then
          celln=celda
       else
          cell=celda
       end if

       return
    End Subroutine Niggli_Cell_Type

    !!--++
    !!--++ Subroutine Niggli_Cell_Vect(A,B,C,Niggli_Point,Celln,Trans)
    !!--++    real(kind=cp),dimension(3),                intent(in)     :: a,b,c
    !!--++    real(kind=cp),dimension(5),      optional, intent(out)    :: Niggli_Point
    !!--++    type(Crystal_Cell_Type),         optional, intent(out)    :: celln
    !!--++    real(kind=cp), dimension(3,3),   optional, intent(out)    :: trans
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the Niggli cell when the input is given as three vectors
    !!--++    in Cartesian components. A test of linear indenpendency is performed.
    !!--++    Calls the subroutine Niggli_Cell_Nigglimat for the effective calculations
    !!--++
    !!--++ Update: October - 2008
    !!
    Subroutine Niggli_Cell_Vect(A,B,C,Niggli_Point,Celln,Trans)
       !---- Arguments ----!
       real(kind=cp),dimension(3),                intent(in)     :: a,b,c
       real(kind=cp),dimension(5),      optional, intent(out)    :: Niggli_Point
       type(Crystal_Cell_Type),         optional, intent(out)    :: celln
       real(kind=cp), dimension(3,3),   optional, intent(out)    :: trans

       !--- Local variables ---!
       type(Crystal_Cell_Type)       :: celda
       real(kind=cp), dimension(2,3) :: n_mat
       real(kind=cp)                 :: det

       det=determ_V(a,b,c)
       if (abs(det) < 0.0001) then
          Err_Crys=.true.
          ERR_Crys_Mess=" The three input vectors are nor linearly independent!"
          return
       end if
       n_mat(1,1)=dot_product(a,a); n_mat(1,2)=dot_product(b,b); n_mat(1,3)=dot_product(c,c)
       n_mat(2,1)=dot_product(b,c); n_mat(2,2)=dot_product(a,c); n_mat(2,3)=dot_product(a,b)

       if (present(Niggli_Point)) then
          if (present(trans)) then
             call Niggli_Cell_nigglimat(n_mat,Niggli_Point,celda,trans)
          else
             call Niggli_Cell_nigglimat(n_mat,Niggli_Point,celda)
          end if
       else if(present(trans)) then
          call Niggli_Cell_nigglimat(n_mat,celln=celda,trans=trans)
       else
          call Niggli_Cell_nigglimat(n_mat,celln=celda)
       end if
       if (Err_Crys) return
       if (present(celln)) celln=celda

       return
    End Subroutine Niggli_Cell_Vect

    !!----
    !!---- Subroutine Read_Bin_Crystal_Cell(Celda,Lun,ok)
    !!----    Type (Crystal_Cell_Type),  intent(out) :: Celda   ! Out -> Cell variable
    !!----    Integer,                   intent(in)  :: lun     !  In -> Unit to write
    !!----    logical,                   intent(out) :: ok
    !!----
    !!----    Reads the cell characteristics in a binary file associated to the
    !!----    logical unit lun. The file is supposed to be opened with form="unformatted",
    !!----    access="stream" or equivalent
    !!----
    !!----    Updated: February - 2013
    !!
    Subroutine Read_Bin_Crystal_Cell(Celda,Lun,ok)
       !---- Arguments ----!
       Type (Crystal_Cell_Type),  intent(out) :: Celda
       Integer,                   intent(in)  :: Lun
       logical,                   intent(out) :: ok
       integer :: ier
       ok=.true.
       read(unit=lun,iostat=ier)             &
                       Celda%cell,           &
                       Celda%ang,            &
                       Celda%cell_std,       &
                       Celda%ang_std,        &
                       Celda%rcell,          &
                       Celda%rang,           &
                       Celda%GD,Celda%GR,    &
                       Celda%Cr_Orth_cel,    &
                       Celda%Orth_Cr_cel,    &
                       Celda%BL_M,           &
                       Celda%BL_Minv,        &
                       Celda%CellVol,        &
                       Celda%RCellVol,       &
                       Celda%CartType
       if( ier /= 0) ok=.false.
       return
    End Subroutine Read_Bin_Crystal_Cell


    !!--++
    !!--++ Subroutine Recip(A,Ang,Ar,Angr,Vol,Volr)
    !!--++    real(kind=cp), dimension(3), intent(in ) :: a        !  In -> a,b,c
    !!--++    real(kind=cp), dimension(3), intent(in ) :: ang      !  In -> alpha,beta,gamma
    !!--++    real(kind=cp), dimension(3), intent(out) :: ar       !  In -> a*,b*,c*
    !!--++    real(kind=cp), dimension(3), intent(out) :: angr     !  In -> alpha*,beta*,gamma*
    !!--++    real(kind=cp),               intent(out) :: vol      ! Out -> Vol
    !!--++    real(kind=cp),               intent(out) :: volr     ! Out -> Vol*
    !!--++
    !!--++    (PRIVATE)
    !!--++    Calculates the reciprocal lattice vectors and cell volume
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Recip(A,Ang,Ar,Angr,Vol,Volr)
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent(in ) :: a,ang
       real(kind=cp), dimension(3), intent(out) :: ar,angr
       real(kind=cp),               intent(out) :: vol,volr

       !---- Local Variables ----!
       integer        :: i
       real(kind=cp)  :: s,p,cose

       p=1.0
       s=1.0
       do i=1,3
          cose=cosd(ang(i))
          p=p*cose
          s=s-cose*cose
       end do
       vol=sqrt(abs(s+2.0*p))

       do i=1,3
          vol=vol*a(i)
       end do
       volr=1.0/vol

       ar(1)=a(2)*a(3)*sind(ang(1))/vol
       ar(2)=a(3)*a(1)*sind(ang(2))/vol
       ar(3)=a(1)*a(2)*sind(ang(3))/vol
       angr(1)=(cosd(ang(2))*cosd(ang(3))-cosd(ang(1)))/(sind(ang(2))*sind(ang(3)))
       angr(2)=(cosd(ang(1))*cosd(ang(3))-cosd(ang(2)))/(sind(ang(1))*sind(ang(3)))
       angr(3)=(cosd(ang(2))*cosd(ang(1))-cosd(ang(3)))/(sind(ang(2))*sind(ang(1)))
       do i=1,3
          angr(i)=acosd(angr(i))
       end do

       return
    End Subroutine Recip

    !!----
    !!---- Subroutine Set_Crystal_Cell(Cellv,Angl,Celda,Cartype,Scell,Sangl)
    !!----    real(kind=cp), dimension (3),        intent(in ) :: cellv   !  In -> a,b,c
    !!----    real(kind=cp), dimension (3),        intent(in ) :: angl    !  In -> angles of cell parameters
    !!----    Type (Crystal_Cell_Type),            intent(out) :: Celda   !  Out-> Celda components
    !!----    character (len=1),          optional,intent(in ) :: CarType !  In -> Type of Cartesian Frame
    !!----    real(kind=cp), dimension(3),optional,intent(in ) :: scell,sangl
    !!----
    !!----    Constructs the object "Celda" of type Crystal_Cell. Control for error
    !!----    is present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Set_Crystal_Cell(Cellv,Angl,Celda,Cartype,Scell,Sangl)
       !---- Arguments ----!
       real(kind=cp), dimension (3),        intent(in ) :: cellv, angl
       Type (Crystal_Cell_Type),            intent(out) :: Celda
       character (len=1),          optional,intent(in ) :: CarType
       real(kind=cp), dimension(3),optional,intent(in ) :: scell,sangl

       !---- Local Variables ----!
       integer :: ifail

       call Init_Err_Crys()

       if (present(scell) .and. present(sangl)) then
          Celda%cell_std=scell
          Celda%ang_std=sangl
       else
          Celda%cell_std=0.0
          Celda%ang_std=0.0
          Celda%lcell=0    !These codes are attributed in refinement programs
          Celda%lang=0     !In order to preserve the values given by these programs the
       end if              !procedure should be invoked with standard deviations

       Celda%cell=cellv
       Celda%ang=angl
       where(Celda%ang < eps) Celda%ang =90.0
       call recip(cellv,angl,Celda%rcell,Celda%rang,Celda%CellVol,Celda%RCellVol)
       if (present(CarType)) then
          call Get_Cryst_Orthog_matrix(cellv,angl,Celda%Cr_Orth_cel,CarType)
          Celda%CartType=CarType
       else
          call Get_Cryst_Orthog_matrix(cellv,angl,Celda%Cr_Orth_cel)
          Celda%CartType="C"
       end if
       call matrix_inverse(Celda%Cr_Orth_cel,Celda%Orth_Cr_cel,ifail)

       if (ifail /= 0) then
          err_crys=.true.
          ERR_Crys_Mess=" Bad cell parameters "
          return
       end if

       Celda%GD=Metrics(cellv,angl)
       Celda%GR=Metrics(Celda%rcell,Celda%rang)

       ! Busing-Levy matrix component
       !(it corresponds to the transpose of Orth_Cr_cel when Celda%CartType="C")
       If (Celda%CartType == "C") then
          Celda%bl_m=Transpose(Celda%Orth_Cr_cel)
          Celda%bl_minv=Transpose(Celda%Cr_Orth_cel)
       else
          Celda%bl_m(1,1)=celda%rcell(1)
          Celda%bl_m(1,2)=celda%rcell(2)*cosd(celda%rang(3))
          Celda%bl_m(1,3)=celda%rcell(3)*cosd(celda%rang(2))
          Celda%bl_m(2,2)=celda%rcell(2)*sind(celda%rang(3))
          Celda%bl_m(2,3)=-(celda%rcell(3)*sind(celda%rang(2))*cosd(celda%ang(1)))
          Celda%bl_m(3,3)=1.0/celda%cell(3)
          Celda%bl_m(2,1)=0.0
          Celda%bl_m(3,1)=0.0
          Celda%bl_m(3,2)=0.0
          call matrix_inverse(Celda%bl_m,Celda%bl_minv,ifail)

          if (ifail /= 0) then
             err_crys=.true.
             ERR_Crys_Mess=" Bad cell parameters "
             return
          end if
       end if

       return
    End Subroutine Set_Crystal_Cell

    !!----
    !!---- Subroutine Volume_Sigma_from_Cell(cell,ang,sigc,siga,volume,sigv)
    !!----    real(kind=cp), dimension(3),  intent(in) :: Cell   !  In  ->  a,b,c parameters
    !!----    real(kind=cp), dimension(3),  intent(in) :: Ang    !  In  -> alpha, beta, gamma
    !!----    real(kind=cp), dimension(3),  intent(in) :: SigC   !  In  -> sigmas for a ,b and c
    !!----    real(kind=cp), dimension(3),  intent(in) :: SigA   !  In  -> sigmas for angles
    !!----    real(kind=cp),                intent(out):: Volume ! Out  -> Volume from cell parameters
    !!----    real(kind=cp),                intent(out):: SigV   ! Out  -> Sigma for Volume
    !!----
    !!----    Calculates the volume and their standard deviation from unit cell
    !!----    parameters. If the standard deviations of cell parameters are zero
    !!----    the result is sigma=0.0, otherwise the calculation is performed.
    !!----    It is assumed that there is no correlation (covariance terms) between
    !!----    the standard deviations of the different cell parameters.
    !!----
    !!---- Updated: January - 2013 (JGP)
    !!
    Subroutine Volume_Sigma_from_Cell(cell,ang,sigc,siga,volume,sigv)
       !---- Arguments ----!
       real(kind=cp), dimension(3),  intent(in) :: Cell   !  In  ->  a,b,c parameters
       real(kind=cp), dimension(3),  intent(in) :: Ang    !  In  -> alpha, beta, gamma
       real(kind=cp), dimension(3),  intent(in) :: SigC   !  In  -> sigmas for a ,b and c
       real(kind=cp), dimension(3),  intent(in) :: SigA   !  In  -> sigmas for angles
       real(kind=cp),                intent(out):: Volume ! Out  -> Volume from cell parameters
       real(kind=cp),                intent(out):: SigV   ! Out  -> Sigma for Volume

       !---- Local Variables ----!
       real(kind=cp) :: a,b,c,ca,cb,cg,sa,sb,sg
       real(kind=cp) :: t, dvda, dvdb, dvdc, dvdalpha, dvdbeta, dvdgamma

       !> Init
       volume=0.0
       sigv=0.0

       a=cell(1)
       b=cell(2)
       c=cell(3)
       ca=cosd(ang(1))
       cb=cosd(ang(2))
       cg=cosd(ang(3))
       sa=sind(ang(1))
       sb=sind(ang(2))
       sg=sind(ang(3))

       t=sqrt(1.0 - ca**2 - cb**2 - cg**2 + 2.0*ca*cb*cg)

       volume=a*b*c*t

       if(sum(abs(sigc)) < eps .and. sum(abs(siga)) < eps ) return

       dvda=b*c*t
       dvdb=a*c*t
       dvdc=a*b*t

       dvdalpha=(a*b*c)*( (sa/t)*(ca-cb*cg) )
       dvdbeta= (a*b*c)*( (sb/t)*(cb-ca*cg) )
       dvdgamma=(a*b*c)*( (sg/t)*(cg-ca*cb) )

       sigv= (dvda*sigc(1))**2 + (dvdb*sigc(2))**2 + (dvdc*sigc(3))**2 +  &
             (dvdalpha*siga(1)*to_rad)**2 + (dvdbeta*siga(2)*to_rad)**2 + &
             (dvdgamma*siga(3)*to_rad)**2

       sigv=sqrt(sigv)

       return
    End Subroutine Volume_Sigma_from_Cell

    !!----
    !!---- Subroutine Write_Crystal_Cell(Celda,Lun)
    !!----    Type (Crystal_Cell_Type),  intent(in)  :: Celda   !  In -> Cell variable
    !!----    Integer,                   intent(in)  :: lun     !  In -> Unit to write
    !!----
    !!----    Writes the cell characteristics in a binary file associated to the
    !!----    logical unit lun. The file is supposed to be opened with form="unformatted",
    !!----    access="stream" or equivalent
    !!----
    !!---- Update: February - 2013
    !!
    Subroutine Write_Bin_Crystal_Cell(Celda,Lun)
       !---- Arguments ----!
       Type (Crystal_Cell_Type),  intent(in) :: Celda
       Integer,                   intent(in) :: Lun
       write(unit=lun) Celda%cell,           &
                       Celda%ang,            &
                       Celda%cell_std,       &
                       Celda%ang_std,        &
                       Celda%rcell,          &
                       Celda%rang,           &
                       Celda%GD,Celda%GR,    &
                       Celda%Cr_Orth_cel,    &
                       Celda%Orth_Cr_cel,    &
                       Celda%BL_M,           &
                       Celda%BL_Minv,        &
                       Celda%CellVol,        &
                       Celda%RCellVol,       &
                       Celda%CartType
       return
    End Subroutine Write_Bin_Crystal_Cell

    !!----
    !!---- Subroutine Write_Crystal_Cell(Celda,Lun)
    !!----    Type (Crystal_Cell_Type),  intent(in)  :: Celda   !  In -> Cell variable
    !!----    Integer,optional           intent(in)  :: lun     !  In -> Unit to write
    !!----
    !!----    Writes the cell characteristics in a file associated to the
    !!----    logical unit lun
    !!----
    !!---- Update: January - 2011
    !!
    Subroutine Write_Crystal_Cell(Celda,Lun)
       !---- Arguments ----!
       Type (Crystal_Cell_Type),  intent(in) :: Celda
       Integer,optional,          intent(in) :: Lun

       !---- Local variables ----!
       integer            :: iunit
       integer            :: i,j

       iunit=6
       if (present(lun)) iunit=lun

       Write(unit=iunit,fmt="(/,a)")    "        Metric information:"
       Write(unit=iunit,fmt="(a,/)")    "        -------------------"
       Write(unit=iunit,fmt="(a,/)")    " => Direct cell parameters:"
       Write(unit=iunit,fmt="(3(a,f12.4))")"         a = ", Celda%cell(1),"      b = ", Celda%cell(2), "      c = ", Celda%cell(3)
       Write(unit=iunit,fmt="(3(a,f12.3))")"     alpha = ", Celda%ang(1) ,"   beta = ", Celda%ang(2) , "  gamma = ", Celda%ang(3)
       Write(unit=iunit,fmt="(a,f12.4)")   "                        Direct Cell Volume = ",Celda%CellVol
       Write(unit=iunit,fmt="(/,a,/)")     " => Reciprocal cell parameters:"
       Write(unit=iunit,fmt="(3(a,f12.6))")"         a*= ", Celda%rcell(1),"      b*= ",Celda%rcell(2),"      c*= ", Celda%rcell(3)
       Write(unit=iunit,fmt="(3(a,f12.3))")"     alpha*= ", Celda%rang(1) ,"   beta*= ",Celda%rang(2) ,"  gamma*= ", Celda%rang(3)
       Write(unit=iunit,fmt="(a,f12.8)")   "                    Reciprocal Cell Volume = ",Celda%RCellVol
       Write(unit=iunit,fmt="(/,a,/)")     " => Direct and Reciprocal Metric Tensors:"
       Write(unit=iunit,fmt="(a)")         "                   GD                                       GR"

       do i=1,3
          Write(unit=iunit,fmt="(3f12.4,a,3f12.6)") (Celda%GD(i,j),j=1,3),"      ", (Celda%GR(i,j),j=1,3)
       end do

       if (Celda%CartType == "A") then
          Write(unit=iunit,fmt="(/,a,/)") " =>  Cartesian frame: x // a; y is in the ab-plane; z is x ^ y   "
       else
          Write(unit=iunit,fmt="(/,a,/)") " =>  Cartesian frame: z // c; y is in the bc-plane; x is y ^ z   "
       end if

       Write(unit=iunit,fmt="(a)")       "     Crystal_to_Orthonormal_Matrix              Orthonormal_to_Crystal Matrix"
       Write(unit=iunit,fmt="(a)")       "              Cr_Orth_cel                               Orth_Cr_cel  "
       do i=1,3
          Write(unit=iunit,fmt="(3f12.4,a,3f12.6)") (Celda%Cr_Orth_cel(i,j),j=1,3),"      ", (Celda%Orth_Cr_cel(i,j),j=1,3)
       end do

       Write(unit=iunit,fmt="(/,a)")     "     Busing-Levy B-matrix: Hc=B.H            Inverse of the Busing-Levy B-matrix"
       Write(unit=iunit,fmt="(a)")       "                BL_M                                      BL_Minv  "
       do i=1,3
          Write(unit=iunit,fmt="(3f12.6,a,3f12.4)") (Celda%BL_M(i,j),j=1,3),"      ", (Celda%BL_Minv(i,j),j=1,3)
       end do

       return
    End Subroutine Write_Crystal_Cell

 End Module CFML_Crystal_Metrics
