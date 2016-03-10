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
!!---- MODULE: CFML_Geometry_Calc
!!----   INFO: Routines for Geometry Calculations
!!----
!!---- HISTORY
!!----    Update: 06/03/2011
!!----
!!----
!!---- DEPENDENCIES
!!--++    CFML_Math_3D:  Cross_Product
!!--++    CFML_GlobalDeps: Eps, Pi, Cp, Sp, To_Rad, To_Deg
!!--++    CFML_Math_General: Acosd, Cosd, Sind
!!--++    CFML_Crystal_Metrics: Crystal_Cell_Type
!!----
!!---- VARIABLES
!!----    COORDINATION_TYPE
!!----    COORD_INFO
!!--++    EPSI
!!----    ERR_GEOM
!!----    ERR_GEOM_MESS
!!----    POINT_LIST_TYPE
!!----
!!---- PROCEDURES
!!----    Functions:
!!----       ANGLE_DIHEDRAL
!!--++       ANGLE_DIHEDRAL_IJKN            [Overloaded]
!!--++       ANGLE_DIHEDRAL_UVW             [Overloaded]
!!----       ANGLE_MOD
!!--++       ANGLE_MODN                     [Overloaded]
!!--++       ANGLE_MODV                     [Overloaded]
!!----       ANGLE_UV
!!--++       ANGLE_UVI                      [Overloaded]
!!--++       ANGLE_UVR                      [Overloaded]
!!----       COORD_MOD
!!--++       COORD_MODN                     [Overloaded]
!!--++       COORD_MODV                     [Overloaded]
!!----       DISTANCE
!!--++       DISTANCE_FR                    [Overloaded]
!!--++       DISTANCE_FR_DP                 [Overloaded]
!!--++       DISTANCE_SC                    [Overloaded]
!!----       MATRIX_PHITHECHI
!!----       MATRIX_RX
!!----       MATRIX_RY
!!----       MATRIX_RZ
!!----
!!----    Subroutines:
!!----       ALLOCATE_COORDINATION_TYPE
!!----       ALLOCATE_POINT_LIST
!!----       CALC_DIST_ANGLE
!!----       CALC_DIST_ANGLE_SIGMA
!!----       DEALLOCATE_COORDINATION_TYPE
!!----       DEALLOCATE_POINT_LIST
!!----       DISTANCE_AND_SIGMA
!!----       GET_ANGLEN_AXIS_FROM_ROTMAT
!!----       GET_EULER_FROM_FRACT
!!----       GET_MATRIX_MOVING_V_TO_U
!!----       GET_OMEGACHIPHI
!!----       GET_PHITHECHI
!!----       GET_TRANSF_LIST
!!----       INIT_ERR_GEOM
!!----       P1_DIST
!!----       PRINT_DISTANCES
!!----       SET_NEW_ASYMUNIT
!!----       SET_ORBITS_INLIST
!!----       SET_ROTATION_MATRIX
!!----       SET_TDIST_COORDINATION
!!----       SET_TDIST_PARTIAL_COORDINATION
!!----
!!
 Module CFML_Geometry_Calc

    !---- Use Modules ----!
    use CFML_GlobalDeps,                 only: Sp, Cp, dp, eps, pi, to_rad, to_deg
    use CFML_Math_General,               only: acosd, cosd, sind, Modulo_Lat
    use CFML_Math_3D,                    only: Cross_Product, Matrix_Inverse, determ_A
    use CFML_String_Utilities,           only: Frac_Trans_1Dig, L_Case,U_Case,pack_string,setnum_std, get_logunit
    use CFML_Crystal_Metrics,            only: Crystal_Cell_Type, Get_Deriv_Orth_Cell,Rot_Matrix
    use CFML_Atom_TypeDef,               only: atom_list_type,Atoms_Cell_Type,Equiv_Atm, Wrt_Lab, Atom_Equiv_List_Type, &
                                               allocate_atom_list
    use CFML_Crystallographic_Symmetry,  only: Space_Group_Type, ApplySo, Lattice_Trans, Get_Multip_Pos, &
                                               searchop, Read_SymTrans_Code, Write_SymTrans_Code, Get_Orbit

    implicit none

    private

    !---- List of public functions ----!

    !---- List of public overloaded procedures: functions ----!
    public :: Angle_Dihedral, Angle_Mod, Angle_Uv, Coord_Mod, Distance, Matrix_PhiTheChi, Matrix_Rx, &
              Matrix_Ry, Matrix_Rz

    !---- List of public subroutines ----!
    public :: Allocate_Coordination_Type, Allocate_Point_List, Calc_Dist_Angle, Calc_Dist_Angle_Sigma, &
              Deallocate_Coordination_Type, Deallocate_Point_List, Distance_and_Sigma, Get_Euler_From_Fract, &
              Get_PhiTheChi, init_err_geom, P1_Dist, Print_Distances, Set_Orbits_InList, Set_TDist_Coordination, &
              Get_Transf_List, Set_TDist_Partial_Coordination, Get_Anglen_Axis_From_RotMat, Get_Matrix_moving_v_to_u, &
              Get_OmegaChiPhi, Set_Rotation_Matrix, Set_New_AsymUnit

    !---- List of public overloaded procedures: subroutines ----!

    !---- List of private functions ----!
    private :: Angle_Dihedral_Uvw,  Angle_Dihedral_Ijkn, Angle_Uvi, Angle_Uvr, Angle_Modn, Angle_Modv, &
               Coord_Modn, Coord_Modv, Distance_fr, Distance_fr_dp, Distance_sc

    !---- Definitions ----!

    !!----
    !!---- TYPE :: COORDINATION_TYPE
    !!--..
    !!---- Type, public :: Coordination_Type
    !!----    integer                                      :: Natoms    ! number of atoms
    !!----    integer                                      :: Max_Coor  ! Maximum number of connected atoms to a given one
    !!----    integer,       dimension(:),     allocatable :: Coord_Num ! Counter of distances connected to the current atom
    !!----    integer,       dimension(:,:),   allocatable :: N_Cooatm  ! Pointer to the ordinal number in the list of the attached
    !!----                                                              ! atom to the atom given by the first index
    !!----    integer,       dimension(:,:),   allocatable :: N_Sym     !
    !!----    real(kind=cp), dimension(:,:),   allocatable :: Dist      ! List of distances related to an atom
    !!----    real(kind=cp), dimension(:,:),   allocatable :: S_Dist    ! List of Sigma(distances)
    !!----    real(kind=cp), dimension(:,:,:), allocatable :: Tr_coo    !
    !!---- End type Coordination_Type
    !!----
    !!---- Update: February - 2005
    !!
    Type, public :: Coordination_Type
       integer                                      :: Natoms    ! number of atoms
       integer                                      :: Max_Coor  ! Maximum number of connected atoms to a given one
       integer,       dimension(:),     allocatable :: Coord_Num ! Counter of distances connected to the current atom
       integer,       dimension(:,:),   allocatable :: N_Cooatm  ! Pointer to the ordinal number in the list of the attached
                                                                 ! atom to the atom given by the first index
       integer,       dimension(:,:),   allocatable :: N_Sym     ! Number of symmetry operator to apply to N_Cooatm
       real(kind=cp), dimension(:,:),   allocatable :: Dist      ! List of distances related to an atom
       real(kind=cp), dimension(:,:),   allocatable :: S_Dist    ! List of Sigma(distances)
       real(kind=cp), dimension(:,:,:), allocatable :: Tr_coo
    End type Coordination_Type

    !!----
    !!---- COORD_INFO
    !!----    type(Coordination_Type), public :: coord_info
    !!----
    !!----    Coordination Information
    !!----
    !!---- Update: March - 2005
    !!
    type(Coordination_Type), public :: coord_info

    !!--++
    !!--++ EPSI
    !!--++    real(kind=cp), parameter :: epsi=0.001
    !!--++
    !!--++    (PRIVATE)
    !!--++    Epsilon for roughly comparing distances
    !!--++
    !!--++ Update: February - 2005
    !!
    real(kind=cp), parameter, private :: epsi=0.001

    !!----
    !!---- ERR_GEOM
    !!----    logical, public  :: err_geom
    !!----
    !!----    Logical Variable indicating an error in CFML_Geometry_Calc module
    !!----
    !!---- Update: February - 2005
    !!
    logical, public  :: err_geom

    !!----
    !!---- ERR_Geom_Mess
    !!----    character(len=150), public :: ERR_Geom_Mess
    !!----
    !!----    String containing information about the last error
    !!----
    !!---- Update: February - 2005
    !!
    character(len=150), public :: ERR_Geom_Mess


    !!----
    !!---- TYPE :: POINT_LIST_TYPE
    !!--..
    !!---- Type, public :: Point_List_Type
    !!----    integer                                       :: np   !number of points in list
    !!----    character(len=12), dimension(:),  allocatable :: nam  !name/label associated to each point
    !!----    integer,           dimension(:),  allocatable :: p    !integer pointer for various purposes
    !!----    real(kind=cp)      dimension(:,:),allocatable :: x    !fractional coordinates of points
    !!---- End type Point_List_Type
    !!----
    !!---- Update: February - 2005
    !!
    Type, public :: point_list_type
       integer                                       :: np   !number of points in list
       character(len=12), dimension(:),  allocatable :: nam  !name/label associated to each point
       integer,           dimension(:),  allocatable :: p    !integer pointer for various purposes
       real(kind=cp),     dimension(:,:),allocatable :: x    !fractional coordinates of points
    End type point_list_type


    !---- Interfaces - Overlapp ----!
    Interface  Angle_Dihedral
       Module Procedure Angle_Dihedral_Ijkn
       Module Procedure Angle_Dihedral_Uvw
    End Interface

    Interface  Angle_Uv
       Module Procedure Angle_UvI
       Module Procedure Angle_UvR
    End Interface

    Interface  Angle_Mod
       Module Procedure Angle_ModN
       Module Procedure Angle_ModV
    End Interface

    Interface  Coord_Mod
       Module Procedure Coord_ModN
       Module Procedure Coord_ModV
    End Interface

    Interface  Distance
       Module Procedure Distance_FR_DP
       Module Procedure Distance_FR
       Module Procedure Distance_SC
    End Interface

 Contains

    !---- Functions ----!

    !!----
    !!---- Function Angle_Dihedral(U,V,W) Or (Ri,Rj,Rk,Rn)   Result(Angle)
    !!----    real(kind=cp), dimension(3), intent( in) :: u       !  In -> Vector 1
    !!----    real(kind=cp), dimension(3), intent( in) :: v       !  In -> Vector 2
    !!----    real(kind=cp), dimension(3), intent( in) :: w       !  In -> Vector 3
    !!----    or
    !!----    real(kind=cp), dimension(3), intent( in) :: ri      !  In -> Vector position ri
    !!----    real(kind=cp), dimension(3), intent( in) :: rj      !  In -> Vector position rj
    !!----    real(kind=cp), dimension(3), intent( in) :: rk      !  In -> Vector position rk
    !!----    real(kind=cp), dimension(3), intent( in) :: rl      !  In -> Vector position rn
    !!----    real(kind=cp)                            :: angle   ! Out -> Dihedral angle
    !!----
    !!----    Calculates the dihedral angle between planes "u-v" and "v-w", where vectors U,V,W
    !!----    are given in cartesian components.
    !!----    Calculates the dihedral angle corresponding to the four points (ri,rj,rk,rn)
    !!----    given in cartesian components. The definition used for the dihedral angle
    !!----    is the following:
    !!--<<
    !!----    Phi(i,j,k,n) = acos { (rij x rjk) (rjk x rkn) / |rij x rjk| / |rjk x rkn| }
    !!----
    !!----    with this definition the sign of Phi is positive if the vector product
    !!----    (rij x rjk) x (rjk x rkn) is in the same direction as rjk, and negative if
    !!----    the direction is opposite.
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Function Angle_Dihedral_Ijkn(Ri,Rj,Rk,Rn) Result(Angle)
    !!--++    real(kind=cp), dimension(3), intent( in) :: ri       !  In -> Vector position ri
    !!--++    real(kind=cp), dimension(3), intent( in) :: rj       !  In -> Vector position rj
    !!--++    real(kind=cp), dimension(3), intent( in) :: rk       !  In -> Vector position rk
    !!--++    real(kind=cp), dimension(3), intent( in) :: rl       !  In -> Vector position rn
    !!--++    real(kind=cp)                            :: angle    ! Out -> Dihedral angle
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the dihedral angle corresponding to the four points (ri,rj,rk,rn)
    !!--++    given in cartesian components. The definition used for the dihedral angle
    !!--++    is the following:
    !!--++
    !!--++    Phi(i,j,k,n) = acos { (rij x rjk) (rjk x rkn) / |rij x rjk| / |rjk x rkn| }
    !!--++
    !!--++    with this definition the sign of Phi is positive if the vector product
    !!--++    (rij x rjk) x (rjk x rkn) is in the same direction as rjk, and negative if
    !!--++    the direction is opposite.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Angle_Dihedral_Ijkn(ri,rj,rk,rn) result(angle)
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent( in) :: ri,rj,rk,rn
       real(kind=cp)                            :: angle

       angle=Angle_Dihedral_Uvw(rj-ri ,rk-rj, rn-rk )

       return
    End Function Angle_Dihedral_Ijkn

    !!--++
    !!--++ Function Angle_Dihedral_Uvw(U,V,W) Result(Angle)
    !!--++    real(kind=cp), dimension(3), intent( in) :: u       !  In -> Vector 1
    !!--++    real(kind=cp), dimension(3), intent( in) :: v       !  In -> Vector 2
    !!--++    real(kind=cp), dimension(3), intent( in) :: w       !  In -> Vector 3
    !!--++    real(kind=cp)                            :: angle   ! Out -> Dihedral angle
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the dihedral angle between planes u-v and v-w
    !!--++    Vectors u,v,w are given in cartesian components.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Angle_Dihedral_Uvw(u,v,w) result(angle)
       !---- Argument ----!
       real(kind=cp), dimension(3), intent( in) :: u,v,w
       real(kind=cp)                            :: angle

       !---- Local variables ----!
       real(kind=cp)               :: uvmod, vwmod, sig
       real(kind=cp), dimension(3) :: uv,vw

       angle=0.0

       uv=cross_product(u,v)
       vw=cross_product(v,w)
       sig = -sign(1.0_cp, dot_product(cross_product(uv,vw),v))
       uvmod=sqrt(dot_product(uv,uv))
       vwmod=sqrt(dot_product(vw,vw))
       if (uvmod < eps .or. vwmod < eps) return
       angle=acosd(dot_product(uv,vw)/uvmod/vwmod)*sig

       return
    End Function Angle_Dihedral_Uvw

    !!----
    !!---- Function Angle_Mod(X) Result (Y)
    !!----     real(kind=cp),               intent(in) :: x
    !!----                  or
    !!----     real(kind=cp), dimension(:), intent(in) :: x
    !!----
    !!----     Calculates the angle [-pi,pi)
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Function Angle_Modn(Angle) Result(AngMod)
    !!--++    real(kind=cp), intent(in) :: Angle    !  In/Out -> Angle
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Transforms angle in radians between -pi and +pi
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Angle_ModN(Angle) Result(AngMod)
       !---- Arguments ----!
       real(kind=cp), intent(in) :: Angle
       real(kind=cp)             :: AngMod

       AngMod=mod(angle+6.0*pi,2.0*pi)
       if (angmod > pi) angmod=angmod-2.0*pi

       return
    End Function Angle_ModN

    !!--++
    !!--++ Function Angle_Modv(V_Angle) Result(VAngMod)
    !!--++    real(kind=cp), dimension(:), intent(in) :: V_Angle
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Transforms angles in radians between -pi and +pi
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Angle_ModV(V_Angle) Result(VAngMod)
       !---- Arguments ----!
       real(kind=cp), dimension(:),intent(in) :: V_Angle
       real(kind=cp), dimension(size(V_Angle)):: VAngMod

       !---- Local Variables ----!
       integer :: i

       VAngMod=mod(V_Angle+6.0*pi,2.0*pi)
       do i=1,size(V_Angle)
          if (VAngMod(i) > pi) VAngMod(i)=VAngMod(i)-2.0*pi
       end do

       return
    End Function Angle_ModV

    !!----
    !!---- Function Angle_Uv(U,V,G) Result(Angle)
    !!----    integer/real(kind=cp), dimension(:), intent( in)     :: u      !  In -> Vector 1
    !!----    integer/real(kind=cp), dimension(:), intent( in)     :: v      !  In -> Vector 2
    !!----    real(kind=cp), dimension(:,:), intent( in), optional :: g      !  In -> Metric tensor
    !!----    real(kind=cp)                                        :: angle  ! Out -> Angle
    !!----
    !!----    Calculates the angle between vectors u and v given in cartesian
    !!----    components. If g is not given cartesian components are assumed.
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Function Angle_UvI(Ui,Vi,G) Result(Angle)
    !!--++    integer, dimension(:),                   intent(in) :: ui      !  In -> Vector 1
    !!--++    integer, dimension(:),                   intent(in) :: vi      !  In -> Vector 2
    !!--++    real(kind=cp), dimension(:,:), optional, intent(in) :: g       !  In -> Metric tensor
    !!--++    real(kind=cp)                                       :: angle   ! Out -> Angle
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the angle between vectors u and v given in cartesians
    !!--++    or fractional components. If g is not given cartesian components
    !!--++    are assumed.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Angle_UvI(Ui,Vi,G) Result(Angle)
       !---- Argument ----!
       integer, dimension(:),   intent( in)                 :: ui
       integer, dimension(:),   intent( in)                 :: vi
       real(kind=cp), dimension(:,:), intent( in), optional :: g   !metric tensor
       real(kind=cp)                                        :: angle

       !---- Local variables ----!
       real(kind=cp)                      :: umod, vmod
       real(kind=cp), dimension(size(ui)) :: u
       real(kind=cp), dimension(size(vi)) :: v

       angle=0.0

       u=real(ui)
       v=real(vi)

       if (present(g)) then
          umod = sqrt(dot_product(u,matmul(g,u)))
          vmod = sqrt(dot_product(v,matmul(g,v)))
          if (umod < eps .or. vmod < eps) return
          angle=acosd(dot_product(u,matmul(g,v))/umod/vmod)
       else
          umod=sqrt(dot_product(u,u))
          vmod=sqrt(dot_product(v,v))
          if (umod < eps .or. vmod < eps) return
          angle=acosd(dot_product(u,v)/umod/vmod)
       end if

       return
    End Function Angle_uvi

    !!--++
    !!--++ Function Angle_Uvr(U,V,G) Result(Angle)
    !!--++    real(kind=cp), dimension(:), intent( in)             :: u      !  In -> Vector 1
    !!--++    real(kind=cp), dimension(:), intent( in)             :: v      !  In -> Vector 2
    !!--++    real(kind=cp), dimension(:,:), intent( in), optional :: g      !  In -> Metric tensor
    !!--++    real(kind=cp)                                        :: angle  ! Out -> Angle
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the angle between vectors u and v given in cartesian
    !!--++    or fractional components. If g is not given cartesian components
    !!--++    are assumed.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Angle_UvR(u,v,g) result(angle)
       !---- Argument ----!
       real(kind=cp), dimension(:),   intent( in)           :: u
       real(kind=cp), dimension(:),   intent( in)           :: v
       real(kind=cp), dimension(:,:), intent( in), optional :: g   !metric tensor
       real(kind=cp)                                        :: angle

       !---- Local variables ----!
       real(kind=cp)   :: umod, vmod

       angle=0.0

       if (present(g)) then
          umod = sqrt(dot_product(u,matmul(g,u)))
          vmod = sqrt(dot_product(v,matmul(g,v)))
          if (umod < eps .or. vmod < eps) return
          angle=acosd(dot_product(u,matmul(g,v))/umod/vmod)
       else
          umod=sqrt(dot_product(u,u))
          vmod=sqrt(dot_product(v,v))
          if (umod < eps .or. vmod < eps) return
          angle=acosd(dot_product(u,v)/umod/vmod)
       end if

       return
    End Function Angle_uvr

    !!----
    !!---- Function Coord_Mod(X) Result (Y)
    !!----    Real(Kind=Cp),               intent(in) :: x
    !!----                  or
    !!----    real(kind=cp), dimension(:), intent(in) :: x
    !!----
    !!----    Calculates the coordinates between [0,1)
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Function Coord_Modn(X) Result (XMod)
    !!--++    real(kind=cp), intent(in) :: x
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Transforms reduced the value between 0 and 1
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Coord_ModN(x) result(Xmod)
       !---- Arguments ----!
       real(kind=cp), intent(in) :: x
       real(kind=cp)             :: xmod

       xmod=mod(x+10.0_cp,1.0_cp)

       return
    End Function Coord_ModN

    !!--++
    !!--++ Function Coord_Modv(X) Result(XMod)
    !!--++    real(kind=cp), dimension(:), intent(in) :: x
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Transforms reduced coordinate between 0 and 1
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Coord_ModV(x) Result(Xmod)
       !---- Arguments ----!
       real(kind=cp), dimension(:), intent(in) :: x
       real(kind=cp), dimension(size(x))       :: xmod

       xmod=mod(x+10.0_cp,1.0_cp)

       return
    End Function Coord_ModV

    !!----
    !!---- Function Distance(X0,X1,Cell or Code) Result(D)
    !!----    real(kind=cp), dimension(3),        intent(in) :: x0     !  In -> Point 1
    !!----    real(kind=cp), dimension(3),        intent(in) :: x1     !  In -> Point 2
    !!----    Type (Crystal_Cell_Type),           intent(in) :: Cell   !  In -> Cell parameters
    !!----    Or
    !!----    real(kind=dp), dimension(3),        intent(in) :: x0     !  In -> Point 1
    !!----    real(kind=dp), dimension(3),        intent(in) :: x1     !  In -> Point 2
    !!----    Type (Crystal_Cell_Type),           intent(in) :: Cell   !  In -> Cell parameters
    !!----    Or
    !!----    character(len=*), optional,         intent(in) :: Code
    !!----    real(kind=cp)                                  :: d      ! Out -> Distance
    !!----
    !!----    Calculate distance between two points.
    !!----       Fractional Coordinates: Use Cell
    !!----       Cartesian Coordiantes: Code="C" or Code=" "
    !!----       Spherical Coordinates: Code="S"
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Function Distance_Fr(X0,X1,Celda) Result(D)
    !!--++    real(kind=cp), dimension(3),  intent(in) :: x0     !  In -> Point 1
    !!--++    real(kind=cp), dimension(3),  intent(in) :: x1     !  In -> Point 2
    !!--++    Type (Crystal_Cell_Type),     intent(in) :: Celda  !  In -> Cell parameters
    !!--++    real(kind=cp)                                  :: d      ! Put -> Distance
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate distance between two points in Fractional
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Distance_Fr(X0,X1,Celda) Result(Dis)
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent(in) :: x0,x1
       type (Crystal_Cell_Type),    intent(in) :: Celda
       real(kind=cp)                           :: dis

       !---- Local Variables ----!
       real(kind=cp), dimension(3) :: xr

       xr = matmul(celda%Cr_Orth_cel,x1-x0)
       dis=sqrt(dot_product(xr,xr))

       return
    End Function Distance_Fr

    !!--++
    !!--++ Function Distance_Fr_dp(X0,X1,Celda) Result(D)
    !!--++    real(kind=dp), dimension(3),  intent(in) :: x0     !  In -> Point 1
    !!--++    real(kind=dp), dimension(3),  intent(in) :: x1     !  In -> Point 2
    !!--++    Type (Crystal_Cell_Type),     intent(in) :: Celda  !  In -> Cell parameters
    !!--++    real(kind=dp)                            :: d      ! Put -> Distance
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate distance between two points in Fractional
    !!--++
    !!--++ Update: February - 2015
    !!
    Function Distance_Fr_dp(X0,X1,Celda) Result(Dis)
       !---- Arguments ----!
       real(kind=dp), dimension(3), intent(in) :: x0,x1
       type (Crystal_Cell_Type),    intent(in) :: Celda
       real(kind=dp)                           :: dis

       !---- Local Variables ----!
       real(kind=dp), dimension(3) :: xr

       xr = matmul(celda%Cr_Orth_cel,x1-x0)
       dis=sqrt(dot_product(xr,xr))

       return
    End Function Distance_Fr_dp

    !!--++
    !!--++ Function Distance_SC(X0,X1,Code) Result(D)
    !!--++    real(kind=cp), dimension(3),        intent(in) :: x0     !  In -> Point 1
    !!--++    real(kind=cp), dimension(3),        intent(in) :: x1     !  In -> Point 2
    !!--++    character(len=*), optional,         intent(in) :: Code
    !!--++    real(kind=cp)                                  :: d      ! Put -> Distance
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate distance between two points in Cartesian or Spherical
    !!--++    If Code =="C" or Blank or not present then the coordinates are Cartesian.
    !!--++    If Code =="S" then the coordinates are spherical (R, Theta, Phi).
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Distance_SC(X0,X1,Code) Result(Dis)
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent(in) :: x0,x1
       character(len=*), optional,  intent(in) :: Code
       real(kind=cp)                           :: dis

       !---- Local Variables ----!
       real(kind=cp), dimension(3) :: xr,xi,xj

       xr=0.0
       if (present(code)) then
          select case (code(1:1))
             case("S","s") ! Spherical
                xi(1)=x0(1)*cosd(x0(3))*sind(x0(2))  ! R * cos(Phi) * sin(Theta)
                xi(2)=x0(1)*sind(x0(3))*sind(x0(2))  ! R * sin(Phi) * sin(Theta)
                xi(3)=x0(1)*cosd(x0(2))              ! R * cos(Theta)

                xj(1)=x1(1)*cosd(x1(3))*sind(x1(2))  ! R * cos(Phi) * sin(Theta)
                xj(2)=x1(1)*sind(x1(3))*sind(x1(2))  ! R * sin(Phi) * sin(Theta)
                xj(3)=x1(1)*cosd(x1(2))              ! R * cos(Theta)

                xr=xi-xj
             case("C","c") ! Cartesian
                xr=x1-x0
          end select
       else
          !---- Cartesian ----!
          xr=x1-x0
       end if
       dis=sqrt(dot_product(xr,xr))

       return
    End Function Distance_SC

    !!----
    !!---- Function Matrix_Phithechi(Phi,Theta,Chi,Code) Result(M)
    !!----    real(kind=cp),                intent(in) :: Phi
    !!----    real(kind=cp),                intent(in) :: Theta
    !!----    real(kind=cp),                intent(in) :: Chi
    !!----    character(len=*), optional,   intent(in) :: Code
    !!----    real(kind=cp), dimension(3,3)            :: M    ! Put -> Active Rotation Matrix
    !!----
    !!----    Calculate the active rotation matrix corresponding to the composition
    !!----    of a positive rotation around z of angle Chi, followed by a positive rotation
    !!----    of angle Theta around the y-axis and a subsequent positive rotation of angle Phi
    !!----    around z. "Positive" means counter-clockwise.
    !!----    The matrix is M = Rz(Phi) . Ry(Theta) . Rz(Chi)
    !!----    The colums represent the components of the unitary vectors {u,v,w} that
    !!----    may be considered as an alternative orthonormal frame to the canonical {i,j,k}.
    !!----    Applying the matrix M to a point in {i,j,k} gives another point in {i,j,k} obtained
    !!----    by the successive application of the three rotations given above. The transpose
    !!----    (inverse) of the M-matrix, when applied to a point in {i,j,k}, gives the coordinates
    !!----    of the same point referred to the frame {u,v,w}. This transpose matrix corresponds
    !!----    to a passive (change or Cartesian frame) rotation leaving the points in the same
    !!----    position with respect to the  {i,j,k} frame.
    !!----    The matrix M when applied to a column vector containing the coordinates of a point
    !!----    with respect to the {u,v,w} frame provides the coordinates of the same point with
    !!----    respect to the {i,j,k} frame.
    !!----    If Code =="R" or Blank or not present then the input angles are given in radians.
    !!----    If Code =="D" then the input angles are given in degrees (Phi, Theta, Chi).
    !!----
    !!---- Update: February - 2005
    !!
    Function Matrix_Phithechi(Phi,Theta,Chi,Code) Result(Mt)
       !---- Arguments ----!
       real(kind=cp),                intent(in) :: Phi
       real(kind=cp),                intent(in) :: Theta
       real(kind=cp),                intent(in) :: Chi
       character(len=*), optional,   intent(in) :: Code
       real(kind=cp), dimension(3,3)            :: Mt

       !---- Local Variables ----!
       real(kind=cp) :: p,t,c

       if (present(code)) then
          select case (code(1:1))
             case("D","d") ! degrees
                p=Phi*to_rad
                t=Theta*to_rad
                c=Chi*to_rad
             case default ! radians
                p=Phi
                t=Theta
                c=Chi
          end select
       else
          !---- radians ----!
          p=Phi
          t=Theta
          c=Chi
       end if
       Mt(1,1)= cos(p)*cos(t)*cos(c)-sin(p)*sin(c)    !
       Mt(2,1)= sin(p)*cos(t)*cos(c)+cos(p)*sin(c)    !  u
       Mt(3,1)=-sin(t)*cos(c)                         !
       Mt(1,2)=-cos(p)*cos(t)*sin(c)-sin(p)*cos(c)    !
       Mt(2,2)=-sin(p)*cos(t)*sin(c)+cos(p)*cos(c)    !  v
       Mt(3,2)= sin(t)*sin(c)                         !
       Mt(1,3)= cos(p)*sin(t)                         !
       Mt(2,3)= sin(p)*sin(t)                         !  w
       Mt(3,3)= cos(t)                                !

       return
    End Function Matrix_Phithechi

    !!----
    !!---- Function Matrix_Rx(Ang,Code) Result(M)
    !!----    real(kind=cp),                      intent(in) :: Ang
    !!----    character(len=*), optional,         intent(in) :: Code
    !!----    real(kind=cp), dimension(3,3)                  :: M    ! Put -> Active Rotation Matrix
    !!----
    !!----    Calculate the active rotation matrix corresponding to the positive rotation
    !!----    of an angle Phi around the x-axis. The transpose matrix corresponds to a
    !!----    passive rotation that changes the orthogonal system to {u,v,w} leaving the point
    !!----    at the same position w.r.t. the canonical {i,j,k} frame.
    !!----    If Code =="R" or Blank or not present then the input angle is given in radians.
    !!----    If Code =="D" then the input angle is given in degrees.
    !!----
    !!---- Update: February - 2005
    !!
    Function Matrix_Rx(Ang,Code) Result(Mt)
       !---- Arguments ----!
       real(kind=cp),               intent(in) :: Ang
       character(len=*), optional,  intent(in) :: Code
       real(kind=cp), dimension(3,3)           :: Mt

       !---- Local Variables ----!
       real(kind=cp) :: p

       if (present(code)) then
          select case (code(1:1))
             case("D","d") ! degrees
                p=Ang*to_rad
             case default ! radians
                p=Ang
          end select
       else
          !---- radians ----!
          p=Ang
       end if
       Mt(1,1)= 1.0        !              1  0  0
       Mt(2,1)= 0.0        !  u           0  c -s     Rx
       Mt(3,1)= 0.0        !              0  s  c
       Mt(1,2)= 0.0        !
       Mt(2,2)= cos(p)     !  v
       Mt(3,2)= sin(p)     !
       Mt(1,3)= 0.0        !
       Mt(2,3)=-sin(p)     !  w
       Mt(3,3)= cos(p)     !

       return
    End Function Matrix_Rx

    !!----
    !!---- Function Matrix_Ry(Ang,Code) Result(M)
    !!----    real(kind=cp),                      intent(in) :: Ang
    !!----    character(len=*), optional,         intent(in) :: Code
    !!----    real(kind=cp), dimension(3,3)                  :: M    ! Put -> Active Rotation Matrix
    !!----
    !!----    Calculate the active rotation matrix corresponding to the positive rotation
    !!----    of an angle Phi around the y-axis. The transpose matrix corresponds to a
    !!----    passive rotation that changes the orthogonal system to {u,v,w} leaving the point
    !!----    at the same position w.r.t. the canonical {i,j,k} frame.
    !!----    If Code =="R" or Blank or not present then the input angle is given in radians.
    !!----    If Code =="D" then the input angle is given in degrees.
    !!----
    !!---- Update: February - 2005
    !!
    Function Matrix_Ry(Ang,Code) Result(Mt)
       !---- Arguments ----!
       real(kind=cp),               intent(in) :: Ang
       character(len=*), optional,  intent(in) :: Code
       real(kind=cp), dimension(3,3)           :: Mt

       !---- Local Variables ----!
       real(kind=cp) :: p

       if (present(code)) then
          select case (code(1:1))
             case("D","d") ! degrees
                p=Ang*to_rad
             case default ! radians
                p=Ang
          end select
       else
          !---- radians ----!
          p=Ang
       end if
       Mt(1,1)= cos(p)  !             c  0  s
       Mt(2,1)= 0.0     !  u          0  1  0      Ry
       Mt(3,1)=-sin(p)  !            -s  0  c
       Mt(1,2)= 0.0     !
       Mt(2,2)= 1.0     !  v
       Mt(3,2)= 0.0     !
       Mt(1,3)= sin(p)  !
       Mt(2,3)= 0.0     !  w
       Mt(3,3)= cos(p)  !

       return
    End Function Matrix_Ry

    !!----
    !!---- Function Matrix_Rz(Ang,Code) Result(M)
    !!----    real(kind=cp),                      intent(in) :: Ang
    !!----    character(len=*), optional,         intent(in) :: Code
    !!----    real(kind=cp), dimension(3,3)                  :: M    ! Put -> Active Rotation Matrix
    !!----
    !!----    Calculate the active rotation matrix corresponding to the positive rotation
    !!----    of an angle Phi around the z-axis. The transpose matrix corresponds to a
    !!----    passive rotation that changes the orthogonal system to {u,v,w} leaving the point
    !!----    at the same position w.r.t. the canonical {i,j,k} frame.
    !!----    If Code =="R" or Blank or not present then the input angle is given in radians.
    !!----    If Code =="D" then the input angle is given in degrees.
    !!----
    !!---- Update: February - 2005
    !!
    Function Matrix_Rz(Ang,Code) Result(Mt)
       !---- Arguments ----!
       real(kind=cp),               intent(in) :: Ang
       character(len=*), optional,  intent(in) :: Code
       real(kind=cp), dimension(3,3)           :: Mt

       !---- Local Variables ----!
       real(kind=cp) :: p

       if (present(code)) then
          select case (code(1:1))
             case("D","d") ! degrees
                p=Ang*to_rad
             case default ! radians
                p=Ang
          end select
       else
          !---- radians ----!
          p=Ang
       end if
       Mt(1,1)= cos(p)  !                 c  -s  0
       Mt(2,1)= sin(p)  !  u              s   c  0    Rz
       Mt(3,1)= 0.0     !                 0   0  1
       Mt(1,2)=-sin(p)  !
       Mt(2,2)= cos(p)  !  v
       Mt(3,2)= 0.0     !
       Mt(1,3)= 0.0     !
       Mt(2,3)= 0.0     !  w
       Mt(3,3)= 1.0     !

       return
    End Function Matrix_Rz

    !---------------------!
    !---- Subroutines ----!
    !---------------------!

    !!----
    !!---- Subroutine Allocate_Coordination_Type(nasu,numops,dmax,Max_Coor)
    !!----    integer,       intent(in) :: nasu      !  In -> Number of atoms in asymmetric unit
    !!----    integer,       intent(in) :: numops    !  In -> Number of S.O. excluding lattice centerings
    !!----    real(kind=cp), intent(in) :: dmax      !  In -> Maximun distance to be calculated
    !!----    integer,      intent(out) :: Max_Coor  !  Maximum coordination allowed
    !!----
    !!----    Allocation of Coordination_Type.
    !!----    Should be called before using this module.
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Allocate_Coordination_Type(nasu,numops,dmax,Max_Coor)
       !---- Arguments ----!
       integer,       intent(in) :: nasu
       integer,       intent(in) :: numops
       real(kind=cp), intent(in) :: dmax
       integer,      intent(out) :: Max_Coor

       !---- local variables ----!
       real(kind=cp), parameter :: r_atom=0.4_cp !Radius of a typical atom

       if (allocated(Coord_Info%Coord_Num)) deallocate(Coord_Info%Coord_Num)
       if (allocated(Coord_Info%N_Cooatm))  deallocate(Coord_Info%N_Cooatm)
       if (allocated(Coord_Info%N_Sym))     deallocate(Coord_Info%N_Sym)
       if (allocated(Coord_Info%Dist))      deallocate(Coord_Info%Dist)
       if (allocated(Coord_Info%S_Dist))    deallocate(Coord_Info%S_Dist)
       if (allocated(Coord_Info%Tr_Coo))    deallocate(Coord_Info%Tr_Coo)


       max_coor= (dmax/r_atom)**3
       max_coor=max(max_coor,nasu*numops)

       !---- Assigninmg the new values ----!
       Coord_Info%Natoms=nasu
       Coord_Info%Max_Coor= max_coor

       allocate (Coord_Info%Coord_Num(nasu))
       allocate (Coord_Info%N_Cooatm(max_coor,nasu))
       allocate (Coord_Info%N_Sym(max_coor,nasu))
       allocate (Coord_Info%Dist(max_coor,nasu))
       allocate (Coord_Info%S_Dist(max_coor,nasu))
       allocate (Coord_Info%Tr_Coo(3,max_coor,nasu))

       Coord_Info%Coord_Num=0
       Coord_Info%N_Cooatm =0
       Coord_Info%N_Sym    =0
       Coord_Info%Dist     =0.0
       Coord_Info%S_Dist   =0.0
       Coord_Info%Tr_Coo   =0.0

       return
    End Subroutine Allocate_Coordination_Type

    !!----
    !!---- Subroutine Allocate_Point_List(N,Pl,Ier)
    !!----    integer,               intent(in)     :: n      !  In -> Dimension for allocating components of the type
    !!----    type(point_list_type), intent(in out) :: pl     !  In Out-> Type with allocatable components
    !!----    integer,               intent(out)    :: ier    !  Out -> if ier /= 0 an error occurred.
    !!----
    !!----    Allocation of an objet of type Point_List_Type
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Allocate_Point_List(n,Pl,Ier)
       !---- Arguments ----!
       integer,               intent(in)     :: n
       type(point_list_type), intent(in out) :: pl
       integer,               intent(out)    :: ier

       ier=0
       if (n <= 0) then
          ier=1
          return
       end if

       if ( .not. allocated(pl%nam) ) allocate(pl%nam(n),stat=ier)
       if ( .not. allocated(pl%p) )   allocate(pl%p(n),stat=ier)
       if ( .not. allocated(pl%x) )   allocate(pl%x(3,n),stat=ier)

       pl%nam= " "
       pl%np=0
       pl%p=0
       pl%x=0.0

       return
    End subroutine Allocate_Point_List

    !!----
    !!---- Subroutine Calc_Dist_Angle(Dmax, Dangl, Cell, Spg, A, Lun)
    !!----    real(kind=cp),            intent(in)             :: dmax   !  In -> Max. Distance to calculate
    !!----    real(kind=cp),            intent(in)             :: dangl  !  In -> Max. distance for angle calculations
    !!----    type (Crystal_cell_type), intent(in)             :: Cell   !  In -> Object of Crytal_Cell_Type
    !!----    type (Space_Group_type),  intent(in)             :: SpG    !  In -> Object of Space_Group_Type
    !!----    type (atom_list_type),   intent(in)             :: A      !  In -> Object of atom_list_type
    !!----    integer,                  optional, intent(in)   :: lun    !  In -> Logical Unit for writing
    !!----
    !!----    Subroutine to calculate distances and angles, below the prescribed distances
    !!----    "dmax" and "dangl" (angles of triplets at distance below "dangl" to an atom),
    !!----    without standard deviations. If dangl=0.0, no angle calculations are done.
    !!----    Needs as input the objects Cell (of type Crystal_cell), SpG (of type Space_Group)
    !!----    and A (of type atom_list, that should be allocated in the calling program).
    !!----    Writes results in file (unit=lun) if lun is present
    !!----    Control for error is present.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Calc_Dist_Angle(Dmax, Dangl, Cell, Spg, A, Lun)
       !---- Arguments ----!
       real(kind=cp),            intent(in)   :: Dmax, Dangl
       type (Crystal_cell_Type), intent(in)   :: Cell
       type (Space_Group_Type),  intent(in)   :: SpG
       type (atom_list_type),    intent(in)   :: A
       integer, optional,        intent(in)   :: lun

       !---- Local Variables ----!
       logical                            :: iprin
       integer                            :: i,j,k,lk,i1,i2,i3,jl,npeq,nn,L,nlines, max_coor,ico
       character(len= 80), dimension(12)  :: texto = " "
       character(len=  5)                 :: nam,nam1,nam2
       character(len= 40)                 :: transla
       character(len=160)                 :: form3
       character(len= 90)                 :: form2= "(a,3I4,a,a,a,a,a,f9.4,a,3F8.4,a,t85,a)"  !  JRC feb 2014 &   ! TR 4 fev. 2013
                                           !  "("" "",3I4,""  ("",a,"")-("",a,""):"",f9.4,""   "",3F8.4,""  "",a,""  "",a)"
       integer, dimension(3)              :: ic1,ic2
       real(kind=cp),    dimension(3)     :: xx,x1,xo,Tn,xr, QD
       real(kind=cp)                      :: T,dd, da1,da2,da12,cang12,ang12,cang1,ang2,ang1

       real(kind=cp), allocatable,dimension(:,:) :: uu
       real(kind=cp), allocatable,dimension(:,:) :: bcoo

       iprin=.false.
       if (present(lun)) then
          if (lun > 0) iprin=.true.
       end if

       call init_err_geom()

       call allocate_coordination_type(A%natoms,Spg%multip,Dmax,Max_coor)
       if(allocated(uu)) deallocate(uu)
       allocate(uu(3,Max_coor))
       if(allocated(bcoo)) deallocate(bcoo)
       allocate(bcoo(3,Max_coor))

       qd(:)=1.0/cell%rcell(:)
       ic2(:)= nint(dmax/cell%cell(:)+2.5_cp)
       ic1(:)=-ic2(:)
       npeq=spg%numops
       if (dangl > epsi .and. iprin ) then
          form3="(""    ("",a,"")-("",a,"")-("",a,""):"",f8.3/"
          form3=trim(form3)//"""    ("",a,"")-("",a,"")-("",a,""):"",f8.3/"
          form3=trim(form3)//"""    ("",a,"")-("",a,"")-("",a,""):"",f8.3/"
          form3=trim(form3)//"""         ("",a,"") :"",3f8.4,""  ("",a,"") :"",3f8.4)"
       end if

       if (spg%centred == 2) then
          npeq=2*npeq
          if (iprin) then
             write(unit=lun,fmt="(/,a)")" => Symmetry operators combined with inversion centre:"
             nlines=1
             do i=SpG%NumOps+1,npeq
                if (mod(i,2) == 0) then
                   write(unit=texto(nlines)(36:70),fmt="(a,i2,a,a)") &
                               " => SYMM(",i,"): ",trim(SpG%SymopSymb(i))
                   nlines=nlines+1
                else
                   write(unit=texto(nlines)( 1:34),fmt="(a,i2,a,a)")  &
                               " => SYMM(",i,"): ",trim(SpG%SymopSymb(i))
                end if
             end do
             do i=1,min(nlines,12)
                write(unit=lun,fmt="(a)") texto(i)
             end do
          end if
       end if

       do i=1,a%natoms
          xo(:)=a%atom(i)%x(:)
          nam=a%atom(i)%lab
          if (iprin) then
             write(unit=lun,fmt="(/,/,a)")"    -------------------------------------------------------------------"
             write(unit=lun,fmt="(a,f8.4,a,a,3f8.4)")   &
                       "    Distances less than",dmax,"  to atom: ",nam, xo
             write(unit=lun,fmt="(a,/,/)")"    -------------------------------------------------------------------"
             write(unit=lun,fmt="(/,/,a,/,/)") & ! TR 4 fev. 2013
             " Orig. extr. p.equiv.           Distance      x_ext   y_ext   z_ext  (tx,ty,tz)     Sym. op."
          end if
          Coord_Info%Coord_Num(i)=0
          ico=0
          do k=1,a%natoms
             lk=1
             uu(:,lk)=xo(:)
             nam1=a%atom(k)%lab
             do j=1,npeq
                xx=ApplySO(Spg%SymOp(j),a%atom(k)%x)
                do i1=ic1(1),ic2(1)
                   do i2=ic1(2),ic2(2)
                      do i3=ic1(3),ic2(3)
                         do_jl:do jl=1,Spg%NumLat
                            Tn(:)=real((/i1,i2,i3/))+Spg%Latt_trans(:,jl)
                            x1(:)=xx(:)+tn(:)
                            do l=1,3
                               t=abs(x1(l)-xo(l))*qd(l)
                               if (t > dmax) cycle do_jl
                            end do
                            do nn=1,lk
                               if (sum(abs(uu(:,nn)-x1(:)))  <= epsi) cycle do_jl
                            end do
                            xr = matmul(cell%cr_orth_cel,x1-xo)
                            dd=sqrt(dot_product(xr,xr))
                            if (dd > dmax .or. dd < 0.001) cycle
                            ico=ico+1

                            if (Coord_Info%Coord_Num(i) > Coord_Info%Max_Coor) then
                               err_geom=.true.
                               ERR_Geom_Mess=" => Too many distances around atom: "//nam
                               return
                            end if

                            lk=lk+1
                            uu(:,lk)=x1(:)
                            Coord_Info%Dist(ico,i)=dd
                            Coord_Info%N_Cooatm(ico,i)=k
                            bcoo(:,ico)=x1(:)
                            Coord_Info%Tr_Coo(:,ico,i)=tn
                            if (iprin) then
                               call Frac_Trans_1Dig(tn,transla)
                               write(unit=lun,fmt=form2) " ",i,k,j,"  (",nam,")-(",nam1,"):",dd,"   ",x1(:), "  "//transla, &
                                                         trim(Spg%SymOpSymb(j)) !JRC Feb2014
                            end if
                         end do do_jl
                      end do !i3
                   end do !i2
                end do !i1
             end do !j
          end do !k

          Coord_Info%Coord_Num(i)=ico
          if (dangl <= epsi) cycle     !loop on "i" still running

          !---- Angle calculations for bonded atoms at distance lower than DANGL

          if (iprin) then
                write(unit=lun,fmt="(/,/,a)")       "   -------------------------------------------------------"
                write(unit=lun,fmt="(a,a,3f8.4)")   "   -  Angles around atom: ",nam, xo
                write(unit=lun,fmt="(a,/)")         "   -------------------------------------------------------"
          end if
          do j=1,Coord_Info%Coord_Num(i)
             if (Coord_Info%dist(j,i) < epsi .or. Coord_Info%dist(j,i) > dangl) cycle
             da1=Coord_Info%dist(j,i)
             i1=Coord_Info%N_Cooatm(j,i)
             nam1=a%atom(i1)%lab
             do k=j+1,Coord_Info%Coord_Num(i)
                if (Coord_Info%dist(k,i) < epsi .OR. Coord_Info%dist(k,i) > dangl) cycle
                da2=Coord_Info%dist(k,i)
                i2=Coord_Info%N_Cooatm(k,i)
                nam2=a%atom(i2)%lab
                xx(:)=bcoo(:,k)-bcoo(:,j)
                xr = matmul(Cell%Cr_Orth_cel,xx)
                da12=sqrt(dot_product(xr,xr))
                cang12=0.5_cp*(da1/da2+da2/da1-da12*da12/da1/da2)
                ang12=acosd(cang12)
                cang1=0.5_cp*(da12/da2+da2/da12-da1*da1/da12/da2)
                ang1=acosd(cang1)
                ang2=180.0_cp-ang1-ang12

                if (iprin) then
                    write(unit=lun,fmt="(/,3(a,f8.4))")  &
                         "     Atm-1   Atm-2   Atm-3            d12 =",da1,"  d23 =",da2,"   d13 =",da12
                    write(unit=lun,fmt=form3)  nam1,nam,nam2,ang12,   &
                         nam,nam2,nam1,ang1, nam,nam1,nam2,ang2,  &
                         nam1,bcoo(:,j),nam2, bcoo(:,k)
                end if
             end do !k
          end do !j
       end do !i

       return
    End Subroutine Calc_Dist_Angle

    !!----
    !!---- Subroutine Calc_Dist_Angle_Sigma(Dmax, Dangl, Cell, Spg, A, Lun, Lun_cons, Lun_cif,filen,rdmax,ramin)
    !!----    real(kind=cp),             intent(in)   :: dmax     !  In -> Max. Distance to calculate
    !!----    real(kind=cp),             intent(in)   :: dangl    !  In -> Max. distance for angle calculations
    !!----    type (Crystal_cell_type),  intent(in)   :: Cell     !  In -> Object of Crytal_Cell_Type
    !!----    type (Space_Group_type),   intent(in)   :: SpG      !  In -> Object of Space_Group_Type
    !!----    type (atom_list_type),     intent(in)   :: A        !  In -> Object of atom_list_type
    !!----    integer, optional,         intent(in)   :: lun      !  In -> Logical Unit for writing
    !!----    integer, optional,         intent(in)   :: lun_cons !  In -> Logical unit for writing restraints
    !!----    integer, optional,         intent(in)   :: lun_cif  !  In -> Logical unit for writing CIF file with distances and angles
    !!----    character(len=*), optional,intent(in)   :: filrest  !  In -> Name of file for writing restraints
    !!----    real(kind=cp),    optional,intent(in)   :: rdmax,ramin  !  Maximum distan and minimum angle for output in restraints file
    !!----
    !!----    Subroutine to calculate distances and angles, below the prescribed distances
    !!----    "dmax" and "dangl" (angles of triplets at distance below "dangl" to an atom),
    !!----    with standard deviations. If dangl=0.0, no angle calculations are done.
    !!----    Needs as input the objects Cell (of type Crystal_cell), SpG (of type Space_Group)
    !!----    and A (or type atom_list, that should be allocated in the calling program).
    !!----    Writes results in file (unit=lun) if the argument lun is present. In case
    !!----    lun_cif is provided, the program writes in the already opened CIF file (in
    !!----    the calling program) the items related to distances. If lun_cons is provided
    !!----    the program writes items containing restraints to the file CFML_restraints.tpcr
    !!----    or to file "filrest" if provided as argument.
    !!----    Control for error is present.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Calc_Dist_Angle_Sigma(Dmax, Dangl, Cell, Spg, A, Lun, Lun_cons, Lun_cif,filrest,rdmax,ramin)
       !---- Arguments ----!
       real(kind=cp),             intent(in)   :: dmax, dangl
       type (Crystal_cell_Type),  intent(in)   :: Cell
       type (Space_Group_Type),   intent(in)   :: SpG
       type (Atom_list_type),     intent(in)   :: A
       integer, optional,         intent(in)   :: lun
       integer, optional,         intent(in)   :: lun_cons
       integer, optional,         intent(in)   :: lun_cif
       character(len=*), optional,intent(in)   :: filrest
       real(kind=cp),    optional,intent(in)   :: rdmax, ramin

       !---- Local Variables ----!
       logical                            :: iprin
       integer,parameter                  :: nconst=500
       integer                            :: i,j,k,lk,i1,i2,i3,jl,nn,L,&
                                             itnum1,itnum2,num_const, max_coor,num_angc,ico
       character(len=  6)                 :: nam,nam1,nam2
       character(len= 40)                 :: transla
       character(len= 20)                 :: text,tex,texton
       character(len=132)                 :: line
       character(len=160)                 :: form3
       character(len= 90)                 :: form2= "(a,3i4,a,a,a,a,a,a12,3F8.4,a,t85,a)"  !  JRC feb 2014 form2= &   ! TR 4 fev. 2013
                                             !"("" "",3I4,""  ("",a,"")-("",a,""):"",f9.4,""   "",3F8.4,""  "",a,""  "",a)"
       integer, dimension(3)              :: ic1,ic2
       integer, dimension(192)            :: itnum
       real(kind=cp),dimension(3,3,6)     :: DerM
       real(kind=cp),    dimension(3)     :: xx,x1,xo,Tn, QD,so,ss,s1,s2,x2,tr1,tr2
       real(kind=cp)                      :: T,dd, da1,da2,da12,cang12,ang12,cang1,ang2,ang1,rest_d,rest_a
       real(kind=cp)                      :: sdd,sda1,sda2,sda12,sang12,sang2,sang1,srel1,srel2,srel12

       real(kind=cp), allocatable, dimension(:,:) :: uu
       real(kind=cp), allocatable, dimension(:,:) :: bcoo
       real(kind=cp), allocatable, dimension(:,:) :: sbcoo
       real(kind=cp), allocatable, dimension(:,:) :: trcoo

       character(len=132), dimension(:), allocatable  :: const_text
       character(len=132), dimension(:), allocatable  :: dist_text
       character(len=132), dimension(:), allocatable  :: angl_text

       character(len=8) :: codesym
       logical :: esta

       !--- write CIF ---------------------------------------------------------------------
       integer, parameter                             :: max_cif_dist_text = 1500
       integer, parameter                             :: max_cif_angl_text = 6000
       integer                                        :: n_cif_dist_text
       integer                                        :: n_cif_angl_text
       character (len=12)                             :: CIF_bond_site_symm_2
       character (len=12)                             :: CIF_angle_site_symm_1
       character (len=12)                             :: CIF_angle_site_symm_3
       character (len=132), dimension(:), allocatable :: cif_dist_text
       character (len=132), dimension(:), allocatable :: cif_angl_text
       !-----------------------------------------------------------------------------------


       iprin=.false.
       if (present(lun)) then
          if (lun > 0) iprin=.true.
       end if
       rest_d=dmax
       rest_a=45.0
       if(present(rdmax)) rest_d=rdmax
       if(present(ramin)) rest_a=ramin
       call init_err_geom()
       call Allocate_Coordination_Type(A%natoms,Spg%Multip,Dmax,max_coor)

       if(allocated(uu)) deallocate(uu)
       allocate(uu(3,max_coor))
       if(allocated(bcoo)) deallocate(bcoo)
       allocate(bcoo(3,max_coor))
       if(allocated(sbcoo)) deallocate(sbcoo)
       allocate(sbcoo(3,max_coor))
       if(allocated(trcoo)) deallocate(trcoo)
       allocate(trcoo(3,max_coor))


       call get_deriv_Orth_cell(cell,DerM,"A")

       if (present(lun_cons)) then
          num_angc=0
          num_const=0
          if(present(filrest)) then
            open (unit=lun_cons, file=trim(filrest), status="replace", action="write")
          else
            open (unit=lun_cons, file="CFML_Restraints.tpcr", status="replace", action="write")
          end if
          write(unit=lun_cons,fmt="(a)") " FILE with lines for soft distance and angle constraints (restraints)."
          write(unit=lun_cons,fmt="(a)") " It is intended to help editing PCR files with restraints by pasting, "
          write(unit=lun_cons,fmt="(a)") " after correcting the values as wished, to the appropriate lines.  "
          write(unit=lun_cons,fmt="(a)") " Lines with repeated identical distances have been excluded because symmetry "
          write(unit=lun_cons,fmt="(a)") " already force a hard constraint."
          write(unit=lun_cons,fmt="(a)") " Accidental coincidences have also been excluded, check that in list of distances! "
          write(unit=lun_cons,fmt="(/,a)")   " Warning! "
          write(unit=lun_cons,fmt="(a,/,a/)") " Symmetry constrained angles have not been eliminated,",&
                                              " this has to be performed by hand!"

          !---- Set ITnum ----!
          i=0
          i1=1
          i2=24
          if (spg%hexa) then
             i1=25
             i2=36
          end if
          do j=1,Spg%multip
             call searchop(SpG%Symop(j)%Rot(:,:),i1,i2,i)
             Itnum(j)=i
          end do
          if (allocated(const_text)) deallocate(const_text)
          allocate(const_text(nconst)) !Maximum number of restraints
          const_text(:)(1:132)=" "
          if (allocated(dist_text)) deallocate(dist_text)
          allocate(dist_text(nconst)) !Maximum number of restraints
          dist_text(:)(1:132)=" "
          if (allocated(angl_text)) deallocate(angl_text)
          allocate(angl_text(nconst)) !Maximum number of restraints
          angl_text(:)(1:132)=" "
       end if

       if (present(lun_cif)) then
          write(unit=lun_cif, fmt='(a)') " "
          write(unit=lun_cif, fmt='(a)') "#=============================================================================#"
          write(unit=lun_cif, fmt='(a)') "#                      UNIT CELL INFORMATION                                  #"
          write(unit=lun_cif, fmt='(a)') "#=============================================================================#"
          write(unit=lun_cif, fmt='(a)') "_symmetry_cell_setting                "//trim(SPG%CrystalSys)
          write(unit=lun_cif, fmt='(a)') "_symmetry_space_group_name_H-M       '"//trim(SPG%SPG_symb)//"'"
          write(unit=lun_cif, fmt='(a)') "_symmetry_space_group_name_Hall      '"//trim(SPG%Hall)//"'"
          write(unit=lun_cif, fmt='(a)') " "
          write(unit=lun_cif, fmt='(a)') "loop_"
          write(unit=lun_cif, fmt='(a)') "    _symmetry_equiv_pos_as_xyz   #<--must include 'x,y,z'"

          do i=1,SPG%multip
             write(unit=lun_cif, fmt='(a)') "'"//trim(SPG%SymopSymb(i))//"'"
          end do
          write(unit=lun_cif, fmt='(a)') " "

          write(unit=lun_cif, fmt='(a)') "#=============================================================================#"
          write(unit=lun_cif, fmt='(a)') "#                       MOLECULAR GEOMETRY                                    #"
          write(unit=lun_cif, fmt='(a)') "#=============================================================================#"

          if (allocated(CIF_dist_text)) deallocate(CIF_dist_text)
          allocate(CIF_dist_text(max_cif_dist_text)) !Maximum number of distances
          CIF_dist_text(:)(1:132)=" "
          if (allocated(CIF_angl_text)) deallocate(CIF_angl_text)
          allocate(CIF_angl_text(max_cif_angl_text)) !Maximum number of angles
          CIF_angl_text(:)(1:132)=" "
          n_cif_dist_text = 0
          n_cif_angl_text = 0
       end if

       qd(:)=1.0/cell%rcell(:)
       ic2(:)= nint(dmax/cell%cell(:)+2.5_cp)
       ic1(:)=-ic2(:)
       if (dangl > epsi .and. iprin ) then
          form3=            "(""    ("",a,"")-("",a,"")-("",a,""):"",a12/"
          form3=trim(form3)//"""    ("",a,"")-("",a,"")-("",a,""):"",a12/"
          form3=trim(form3)//"""    ("",a,"")-("",a,"")-("",a,""):"",a12/"
          form3=trim(form3)//"""         ("",a,"") :"",3f9.5,""  ("",a,"") :"",3f9.5)"
       end if
       do i=1,a%natoms
          xo(:)=a%atom(i)%x(:)
          so(:)=a%atom(i)%x_std(:)
          nam=a%atom(i)%lab
          Select Case (len_trim(nam))
             case(1)
                nam="  "//trim(nam)
             case(2:5)
                nam=" "//trim(nam)
          End Select
          if (iprin) then
             write(unit=lun,fmt="(/,/,a)")"    -------------------------------------------------------------------"
             write(unit=lun,fmt="(a,f8.4,a,a,3f8.4)")   &
                       "    Distances less than",dmax,"  to atom: ",nam, xo
             write(unit=lun,fmt="(a,/,/)")"    -------------------------------------------------------------------"
             write(unit=lun,fmt="(/,/,a,/,/)") &		! TR 4 fev. 2013
                  " Orig. extr. p.equiv.           Distance      x_ext   y_ext   z_ext  (tx,ty,tz)     Sym. op."
          end if

          ico=0
          do k=1,a%natoms
             lk=1
             uu(:,lk)=xo(:)
             nam1=a%atom(k)%lab
             Select Case (len_trim(nam1))
               case(1)
                  nam1="  "//trim(nam1)
               case(2:5)
                  nam1=" "//trim(nam1)
             End Select
             ss(:)=A%atom(k)%x_std(:)
             do j=1,Spg%Multip
                xx=ApplySO(Spg%SymOp(j),a%atom(k)%x)

                do i1=ic1(1),ic2(1)
                   do i2=ic1(2),ic2(2)
                      do_i3:do i3=ic1(3),ic2(3)

                            Tn(:)=real((/i1,i2,i3/))
                            x1(:)=xx(:)+tn(:)
                            do l=1,3
                               t=abs(x1(l)-xo(l))*qd(l)
                               if (t > dmax) cycle  do_i3
                            end do
                            do nn=1,lk
                               if (sum(abs(uu(:,nn)-x1(:)))  <= epsi) cycle  do_i3
                            end do
                            call distance_and_sigma(Cell,DerM,xo,x1,so,ss,dd,sdd)
                            if (dd > dmax .or. dd < 0.001) cycle
                            ico=ico+1
                            if (Coord_Info%Coord_Num(i) > Coord_Info%Max_Coor) then
                               err_geom=.true.
                               ERR_Geom_Mess=" => Too many distances around atom: "//nam
                               return
                            end if
                            lk=lk+1
                            uu(:,lk)=x1(:)

                            Coord_Info%Dist(ico,i)=dd
                            Coord_Info%S_Dist(ico,i)=sdd
                            Coord_Info%N_Cooatm(ico,i)=k
                            Coord_Info%N_sym(ico,i)=j
                            Coord_Info%Tr_Coo(:,ico,i)=tn

                            bcoo(:,ico)=x1(:)
                            sbcoo(:,ico)=ss(:)
                            trcoo(:,ico)=Tn(:)
                            if (iprin) then
                               call Frac_Trans_1Dig(tn,transla)
                               call setnum_std(dd,sdd,text)
                               !write(unit=lun,fmt=form2) i,k,j,nam,nam1,dd,x1(:), transla, trim(Spg%SymOpSymb(j))! TR 4 fev. 2013
                               write(unit=lun,fmt=form2) " ",i,k,j,"  (",nam,")-(",nam1,"):",text,x1(:), "  "//transla, &
                                                         trim(Spg%SymOpSymb(j)) !JRC Feb2014
                            end if

                            if(present(lun_cons) .and. dd <= rest_d) then
                              esta=.false.
                              write(unit=line,fmt="(a4,tr2,a4,i5,3f10.5,tr5,2f7.4)") A%atom(i)%lab ,A%atom(k)%lab ,&
                                     Itnum(j), tn(:)+SpG%Symop(j)%tr(:) ,dd, sdd
                              if(num_const == 0) then
                                const_text(1)=line(1:132)
                                num_const=1
                                write(unit=dist_text(1),fmt="(a,2f9.5,a)") "DFIX ",dd,sdd, &
                                                                           "  "//trim(A%atom(i)%lab)//"  "//trim(A%atom(k)%lab)
                                call Write_SymTrans_Code(j,tn,codesym)
                                dist_text(1)=trim(dist_text(1))//codesym
                              else
                                do l=num_const,1,-1
                                 if( (line(1:4) == const_text(l)(1:4) .and. line(7:10) == const_text(l)(7:10)) .or. &
                                     (line(1:4) == const_text(l)(7:10) .and. line(7:10) == const_text(l)(1:4)) ) then
                                   if(line(51:132) == const_text(l)(51:132)) then
                                        esta=.true.
                                        exit
                                   end if
                                 end if
                                end do
                                if(.not. esta) then
                                  num_const=num_const+1
                                  if(num_const > NCONST) then
                                     num_const=num_const-1
                                  end if
                                  const_text(num_const)=line(1:132)
                                  write(unit=dist_text(num_const),fmt="(a,2f9.5,a)") "DFIX ",dd,sdd,&
                                        "  "//trim(A%atom(i)%lab)//"  "//trim(A%atom(k)%lab)
                                  call Write_SymTrans_Code(j,tn,codesym)
                                  dist_text(num_const)=trim(dist_text(num_const))//trim(codesym)
                                end if
                              end if
                            end if

                            if(present(lun_cif) .and. n_cif_dist_text < max_cif_dist_text) then
                               call setnum_std(dd,sdd,text)
                               n_cif_dist_text = n_cif_dist_text + 1

                               !if(i1==0 .and. i2==0 .and. i3==0 .and. j==1) then
                               ! write(unit=CIF_bond_site_symm_2, fmt='(a)') "       . ?"
                               !else
                                write(unit=CIF_bond_site_symm_2, fmt='(a,i3, a, 3i1,a)') " ", j, "_", nint(tn+5.0), " ?"
                               !end if

                               write(unit=CIF_dist_text(n_cif_dist_text), fmt='(6a)') &
                                     A%atom(i)%lab(1:4), "  ", A%atom(k)%lab(1:4), " ", text(1:12), CIF_bond_site_symm_2
                            end if
                      end do do_i3 !i3
                   end do !i2
                end do !i1
             end do !j
          end do !k

          Coord_Info%Coord_Num(i)=ico
          if (dangl <= epsi) cycle     !loop on "i" still running

          !---- Angle calculations for bonded atoms at distance lower than DANGL
          if (present(lun_cons)) write(unit=lun_cons,fmt="(a,a)")"=> Help for possible angle restraints around atom ",A%atom(i)%lab

          if (iprin) then
             write(unit=lun,fmt="(/,/,a)")       "   -------------------------------------------------------"
             write(unit=lun,fmt="(a,a,3f8.4)")   "   -  Angles around atom: ",nam, xo
             write(unit=lun,fmt="(a,/)")         "   -------------------------------------------------------"
          end if
          do j=1,Coord_Info%Coord_Num(i)
             if (Coord_Info%Dist(j,i) < epsi .or. Coord_Info%Dist(j,i) > dangl) cycle
             da1=Coord_Info%Dist(j,i)
             sda1=Coord_Info%S_Dist(j,i)
             i1=Coord_Info%N_Cooatm(j,i)
             nam1=a%atom(i1)%lab
             Select Case (len_trim(nam1))
               case(1)
                  nam1="  "//trim(nam1)
               case(2:5)
                  nam1=" "//trim(nam1)
             End Select
             if (present(lun_cons)) then
               itnum1=itnum(Coord_Info%N_sym(j,i))
               tr1(:)=trcoo(:,j)+SpG%Symop(Coord_Info%N_sym(j,i))%tr(:)
             end if
             do k=j+1,Coord_Info%Coord_Num(i)
                if (Coord_Info%Dist(k,i) < epsi .OR. Coord_Info%Dist(k,i) > dangl) cycle
                da2=Coord_Info%Dist(k,i)
                sda2=Coord_Info%S_Dist(k,i)
                i2=Coord_Info%N_Cooatm(k,i)
                nam2=a%atom(i2)%lab
                Select Case (len_trim(nam2))
                  case(1)
                     nam2="  "//trim(nam2)
                  case(2:5)
                     nam2=" "//trim(nam2)
                End Select
                if (present(lun_cons)) then
                  itnum2=itnum(Coord_Info%N_sym(k,i))
                  tr2(:)=trcoo(:,k)+SpG%Symop(Coord_Info%N_sym(k,i))%tr(:)
                end if
                x1(:)=bcoo(:,k)
                x2(:)=bcoo(:,j)
                s1(:)=sbcoo(:,k)
                s2(:)=sbcoo(:,j)
                call distance_and_sigma(Cell,derM,x1,x2,s1,s2,da12,sda12)
                if( da12 < 0.0001) cycle

                cang12=0.5_cp*(da1/da2+da2/da1-da12*da12/da1/da2)
                ang12=ACOSd(cang12)
                cang1=0.5_cp*(da12/da2+da2/da12-da1*da1/da12/da2)
                ang1=ACOSd(cang1)
                ang2=180.0_cp-ang1-ang12

               ! if(abs(abs(cang12)-1.0) < 0.0001) then
               !   sang12=0.0
               ! else
               !  dcang121=(1.0/da2-cang12/da1)**2
               !  dcang122=(1.0/da1-cang12/da2)**2
               !  dcang1212=(da12/da2/da1)**2
               !  sang12=sqrt((dcang121*sda1**2+dcang122*sda2**2+dcang1212*sda12**2)/(1.0-cang12**2))*to_deg
               ! end if
               ! if(abs(abs(cang1)-1.0) < 0.0001) then
               !   sang1=0.0
               ! else
               !  dcang112=(1.0/da2-cang1/da12)**2
               !  dcang12=(1.0/da12-cang1/da2)**2
               !  dcang11=(da1/da2/da12)**2
               !  sang1=sqrt((dcang11*sda1**2+dcang12*sda2**2+dcang112*sda12**2)/(1.0-cang1**2))*to_deg
               ! end if
               ! sang2=sqrt(sang1**2+sang12**2)


                !---- Alternative calculation of angles' sigmas ----!
                srel1=(sda1/da1)**2
                srel12=(sda12/da12)**2
                srel2=(sda2/da2)**2
                sang12=SQRT(srel1+srel2+(sda12*da12/da1/da2)**2)*to_deg
                sang1=SQRT(srel12+srel2+(sda1*da1/da2/da12)**2)*to_deg
                sang2=SQRT(srel12+srel1+(sda2*da2/da1/da12)**2)*to_deg

                if (iprin) then
                   call setnum_std(da1,sda1,tex)
                   call setnum_std(da2,sda2,text)
                   call setnum_std(da12,sda12,texton)
                   write(unit=lun,fmt="(/,a,3a21)")  &
                        "     Atm-1   Atm-2   Atm-3           "," d12 ="//tex,"  d23 ="//text,"   d13 ="//texton
                   call setnum_std(ang12,sang12,tex)
                   call setnum_std(ang1,sang1,text)
                   call setnum_std(ang2,sang2,texton)
                   write(unit=lun,fmt=form3)  nam1,nam,nam2,tex,    &
                                              nam,nam2,nam1,text,   &
                                              nam,nam1,nam2,texton, &
                                              nam1,bcoo(:,j),  nam2, bcoo(:,k)
                end if

                if (present(lun_cons)) then

                  if(ang2 >= rest_a) &
                  write(unit=lun_cons,fmt="(3(a6,tr1),i3,i4,tr1,3f8.4,tr1,3f8.4,2f7.2)") &
                  A%atom(i)%lab ,nam1 ,nam2 ,itnum1,itnum2,tr1(:),tr2(:),ang2,sang2

                  if(ang1 >= rest_a) &
                  write(unit=lun_cons,fmt="(3(a6,tr1),i3,i4,tr1,3f8.4,tr1,3f8.4,2f7.2)") &  !Another angle of the same triangle
                  A%atom(i)%lab ,nam2 ,nam1 ,itnum2,itnum1,tr2(:),tr1(:),ang1,sang1

                  if(ang12 >= rest_a .and. itnum1==1 .and. sum(abs(tr1)) < 0.001) & !Good constraint
                  write(unit=lun_cons,fmt="(3(a6,tr1),i3,i4,tr1,3f8.4,tr1,3f8.4,2f7.2)") &
                  adjustl(nam1),A%atom(i)%lab ,nam2 ,itnum1,itnum2,tr1(:),tr2(:),ang12,sang12

                  if(ang12 >= rest_a .and. itnum2==1 .and. sum(abs(tr2)) < 0.001) & !Good constraint
                  write(unit=lun_cons,fmt="(3(a6,tr1),i3,i4,tr1,3f8.4,tr1,3f8.4,2f7.2)") &  !Another angle of the same triangle
                  adjustl(nam2)," "//A%atom(i)%lab ,nam1 ,itnum2,itnum1,tr2(:),tr1(:),ang12,sang12

                  if(num_angc == 0) then

                    if(ang2 >= rest_a) then
                      num_angc=num_angc+1
                      line=" "
                      write(unit=line,fmt="(a,2f9.3,a)") "AFIX ",ang2,sang2,&
                                                         "  "//trim(A%atom(i)%lab)//" "//trim(nam1)
                      call Write_SymTrans_Code(Coord_Info%N_sym(j,i),trcoo(:,j),codesym)
                      line=trim(line)//trim(codesym)//" "//trim(nam2)
                      call Write_SymTrans_Code(Coord_Info%N_sym(k,i),trcoo(:,k),codesym)
                      line=trim(line)//trim(codesym)
                      angl_text(1)=line(1:132)
                    end if
                    !Repeating with another angle of the same triangle
                    if(ang1 >= rest_a) then
                      num_angc=num_angc+1
                      line=" "
                      write(unit=line,fmt="(a,2f9.3,a)") "AFIX ",ang1,sang1,&
                                                         "  "//trim(A%atom(i)%lab)//" "//trim(nam2)
                      call Write_SymTrans_Code(Coord_Info%N_sym(k,i),trcoo(:,k),codesym)
                      line=trim(line)//trim(codesym)//" "//trim(nam1)
                      call Write_SymTrans_Code(Coord_Info%N_sym(j,i),trcoo(:,j),codesym)
                      line=trim(line)//trim(codesym)
                      angl_text(num_angc)=line(1:132)
                    end if

                    if(ang12 >= rest_a .and. itnum1==1 .and. sum(abs(tr1)) < 0.001) then
                      num_angc=num_angc+1
                      line=" "
                      write(unit=line,fmt="(a,2f9.3,a)") "AFIX ",ang12,sang12,&
                                                         " "//trim(nam1)//"  "//trim(A%atom(i)%lab)
                      call Write_SymTrans_Code(1,(/0.0,0.0,0.0/),codesym)
                      line=trim(line)//trim(codesym)//" "//trim(nam2)
                      call Write_SymTrans_Code(Coord_Info%N_sym(k,i),trcoo(:,k),codesym)
                      line=trim(line)//trim(codesym)
                      angl_text(num_angc)=line(1:132)
                    end if

                    if(ang12 >= rest_a .and. itnum2==1 .and. sum(abs(tr2)) < 0.001) then
                      num_angc=num_angc+1
                      line=" "
                      write(unit=line,fmt="(a,2f9.3,a)") "AFIX ",ang12,sang12,&
                                                         " "//trim(nam2)//"  "//trim(A%atom(i)%lab)
                      call Write_SymTrans_Code(1,(/0.0,0.0,0.0/),codesym)
                      line=trim(line)//trim(codesym)//" "//trim(nam1)
                      call Write_SymTrans_Code(Coord_Info%N_sym(j,i),trcoo(:,j),codesym)
                      line=trim(line)//trim(codesym)
                      angl_text(num_angc)=line(1:132)
                    end if

                  else

                    if(ang2 >= rest_a) then
                      line=" "
                      write(unit=line,fmt="(a,2f9.3,a)") "AFIX ",ang2,sang2,&
                                                         "  "//trim(A%atom(i)%lab)//" "//trim(nam1)
                      call Write_SymTrans_Code(Coord_Info%N_sym(j,i),trcoo(:,j),codesym)
                      line=trim(line)//trim(codesym)//" "//trim(nam2)
                      call Write_SymTrans_Code(Coord_Info%N_sym(k,i),trcoo(:,k),codesym)
                      line=trim(line)//trim(codesym)

                      esta=.false.
                      jl=index(line,"_")
                      if(jl == 0) jl=len_trim(line)
                      do l=num_angc,1,-1
                       if( line(1:jl) == angl_text(l)(1:jl)) then
                           esta=.true.
                           exit
                       end if
                      end do
                      if(.not. esta) then
                        num_angc=num_angc+1
                        if(num_angc > NCONST) num_angc=NCONST
                        angl_text(num_angc)=line(1:132)
                      end if
                    end if

                    if(ang1 >= rest_a) then
                      !Repeating with another angle of the same triangle
                      line=" "
                      write(unit=line,fmt="(a,2f9.3,a)") "AFIX ",ang1,sang1,&
                                                         "  "//trim(A%atom(i)%lab)//" "//trim(nam2)
                      call Write_SymTrans_Code(Coord_Info%N_sym(k,i),trcoo(:,k),codesym)
                      line=trim(line)//trim(codesym)//" "//trim(nam1)
                      call Write_SymTrans_Code(Coord_Info%N_sym(j,i),trcoo(:,j),codesym)
                      line=trim(line)//trim(codesym)

                      esta=.false.
                      jl=index(line,"_")
                      if(jl == 0) jl=len_trim(line)
                      do l=num_angc,1,-1
                       if( line(1:jl) == angl_text(l)(1:jl)) then
                           esta=.true.
                           exit
                       end if
                      end do
                      if(.not. esta) then
                        num_angc=num_angc+1
                        if(num_angc > NCONST) num_angc=NCONST
                        angl_text(num_angc)=line(1:132)
                      end if
                    end if

                    if(ang12 >= rest_a .and. itnum1==1 .and. sum(abs(tr1)) < 0.001) then
                      !Repeating with another angle of the same triangle
                      line=" "
                      write(unit=line,fmt="(a,2f9.3,a)") "AFIX ",ang12,sang12,&
                                                          " "//trim(nam1)//"  "//trim(A%atom(i)%lab)
                      call Write_SymTrans_Code(1,(/0.0,0.0,0.0/),codesym)
                      line=trim(line)//trim(codesym)//" "//trim(nam2)
                      call Write_SymTrans_Code(Coord_Info%N_sym(k,i),trcoo(:,k),codesym)
                      line=trim(line)//trim(codesym)

                      esta=.false.
                      jl=index(line,"_")
                      if(jl == 0) jl=len_trim(line)
                      do l=num_angc,1,-1
                       if( line(1:jl) == angl_text(l)(1:jl)) then
                           esta=.true.
                           exit
                       end if
                      end do
                      if(.not. esta) then
                        num_angc=num_angc+1
                        if(num_angc > NCONST) num_angc=NCONST
                        angl_text(num_angc)=line(1:132)
                      end if
                    end if

                    if(ang12 >= rest_a .and. itnum1==2 .and. sum(abs(tr2)) < 0.001) then
                      !Repeating with another angle of the same triangle
                      line=" "
                      write(unit=line,fmt="(a,2f9.3,a)") "AFIX ",ang12,sang12,&
                                                          " "//trim(nam2)//"  "//trim(A%atom(i)%lab)
                      call Write_SymTrans_Code(1,(/0.0,0.0,0.0/),codesym)
                      line=trim(line)//trim(codesym)//" "//trim(nam1)
                      call Write_SymTrans_Code(Coord_Info%N_sym(j,i),trcoo(:,j),codesym)
                      line=trim(line)//trim(codesym)

                      esta=.false.
                      jl=index(line,"_")
                      if(jl == 0) jl=len_trim(line)
                      do l=num_angc,1,-1
                       if( line(1:jl) == angl_text(l)(1:jl)) then
                           esta=.true.
                           exit
                       end if
                      end do
                      if(.not. esta) then
                        num_angc=num_angc+1
                        if(num_angc > NCONST) num_angc=NCONST
                        angl_text(num_angc)=line(1:132)
                      end if
                    end if

                  end if

                end if !present(lun_cons)

                if (present(lun_cif) .and. n_cif_angl_text < max_cif_angl_text .and. ang12 > 45.0) then
                   !--- Change: I have included the condition ang12 > 45 for selecting the angle to write in
                   !            the CIF file, normally the angles below that value are irrelevant from the
                   !            chemical point of view.
                   ! j: indice de l'operateur de symetrie pour atome 1 ----! No!, Now j and k correspond to indices
                   ! k: indice de l'operateur de symetrie pour atome 2 ----!      running on coordination around atom i!
                   ! tr1: translation associee a op_j ----No!      =>  trcoo(:,j)!
                   ! tr2: translation associee a op_k ?  ----No!   =>  trcoo(:,k)!
                   n_cif_angl_text = n_cif_angl_text + 1

                   !  The commented lines correspond to wrong selections of translations!!!!!!!
                   !  Moreover the indices j and k were taken as the ordinal numbers of the
                   !  symmetry operators and that's not true!
                   !if (j==1 .and. nint(tr1(1))==0 .and. nint(tr1(2))==0 .and. nint(tr1(3))==0) then
                   !   write(unit=CIF_angle_site_symm_1, fmt='(a)') "       ."
                   !else
                   !   write(unit=CIF_angle_site_symm_1, fmt='(a,i3, a, 3I1)') " ", j, "_",  &
                   !         nint(tr1(1)+5.0), nint(tr1(2)+5.0), nint(tr1(3)+5.0)
                   !end if
                   !if (k==1 .and. nint(tr2(1))==0 .and. nint(tr2(2))==0 .and. nint(tr2(3))==0) then
                   !   write(unit=CIF_angle_site_symm_3, fmt='(a)') "  .  ?"
                   !else
                   !   write(unit=CIF_angle_site_symm_3, fmt='(a,i3, a, 3I1,a)') " ", k, "_", &
                   !         nint(tr2(1)+5.0), nint(tr2(2)+5.0), nint(tr2(3)+5.0), " ?"
                   !end if

                   write(unit=CIF_angle_site_symm_1, fmt='(a,i3, a, 3I1)') " ", &
                         Coord_Info%N_sym(j,i), "_", nint(trcoo(:,j)+5.0)
                   write(unit=CIF_angle_site_symm_3, fmt='(a,i3, a, 3I1,a)') " ", &
                         Coord_Info%N_sym(k,i), "_", nint(trcoo(:,k)+5.0), " ?"

                   write(unit=CIF_angl_text(n_cif_angl_text), fmt='(10a)')        &
                         nam1(1:4)," ", nam(1:4), " ",nam2, tex(1:12), " ",       &
                         trim(CIF_angle_site_symm_1), " ", trim(CIF_angle_site_symm_3)

                end if
             end do !k
          end do !j
       end do !i

       if (present(lun_cons)) then
          write(unit=lun_cons,fmt="(/,a,i5)")"=> Total number of independent distances: ",num_const
          write(unit=lun_cons,fmt="(a,/)")   "   List of possible restraints: "
          write(unit=lun_cons,fmt="(a)")" At1   At2  ITnum     T1        T2        T3          DIST   SIGMA"
          do i=1,num_const
             write(unit=lun_cons,fmt="(2x,a)") trim(const_text(i))
          end do

          write(unit=lun_cons,fmt="(/,a)")   "   ========================================= "
          write(unit=lun_cons,fmt="(a  )")   "   List of possible restraints in CFL format "
          write(unit=lun_cons,fmt="(a,/)")   "   ========================================= "


          write(unit=lun_cons,fmt="(/a,i5)")"=> Total number of independent distance restraints: ",num_const
          do i=1,num_const
             write(unit=lun_cons,fmt="(a)") trim(dist_text(i))
          end do
          write(unit=lun_cons,fmt="(/a,i5)")"=> Total number of possible angle restraints: ",num_angc
          do i=1,num_angc
             write(unit=lun_cons,fmt="(a)") trim(angl_text(i))
          end do
          close(unit=lun_cons)
       end if

       if (present(lun_cif)) then
          if (n_CIF_dist_text /=0) then
             write(unit=lun_cif, fmt='(a)') "loop_"
             write(unit=lun_cif, fmt='(a)') "   _geom_bond_atom_site_label_1"
             write(unit=lun_cif, fmt='(a)') "   _geom_bond_atom_site_label_2"
             write(unit=lun_cif, fmt='(a)') "   _geom_bond_distance"
             write(unit=lun_cif, fmt='(a)') "   _geom_bond_site_symmetry_2"
             write(unit=lun_cif, fmt='(a)') "   _geom_bond_publ_flag"

             do i=1, n_CIF_dist_text
                write(unit=lun_CIF, fmt='(a)') trim(CIF_dist_text(i))
             end do
          end if

          if (n_CIF_angl_text /=0) then
             write(unit=lun_cif, fmt='(a)') ""
             write(unit=lun_cif, fmt='(a)') "loop_"
             write(unit=lun_cif, fmt='(a)') "   _geom_angle_atom_site_label_1"
             write(unit=lun_cif, fmt='(a)') "   _geom_angle_atom_site_label_2"
             write(unit=lun_cif, fmt='(a)') "   _geom_angle_atom_site_label_3"
             write(unit=lun_cif, fmt='(a)') "   _geom_angle"
             write(unit=lun_cif, fmt='(a)') "   _geom_angle_site_symmetry_1"
             write(unit=lun_cif, fmt='(a)') "   _geom_angle_site_symmetry_3"
             write(unit=lun_cif, fmt='(a)') "   _geom_angle_publ_flag"

             do i=1, n_CIF_angl_text
              write(unit=lun_CIF, fmt='(a)') trim(CIF_angl_text(i))
             end do
          end if
       end if

       return
    End Subroutine Calc_Dist_Angle_Sigma

    !!----
    !!---- Subroutine Deallocate_Coordination_Type()
    !!----
    !!----    Deallocation of Coordination_Type.
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Deallocate_Coordination_Type()

       if (allocated(Coord_Info%Coord_Num)) deallocate(Coord_Info%Coord_Num)
       if (allocated(Coord_Info%N_Cooatm))  deallocate(Coord_Info%N_Cooatm)
       if (allocated(Coord_Info%N_Sym))     deallocate(Coord_Info%N_Sym)
       if (allocated(Coord_Info%Dist))      deallocate(Coord_Info%Dist)
       if (allocated(Coord_Info%S_Dist))    deallocate(Coord_Info%S_Dist)
       if (allocated(Coord_Info%Tr_Coo))    deallocate(Coord_Info%Tr_Coo)

       !---- Assigninmg the new values ----!
       Coord_Info%Natoms=0
       Coord_Info%Max_Coor= 0

       return
    End Subroutine Deallocate_Coordination_Type

    !!----
    !!---- Subroutine Deallocate_Point_List(Pl)
    !!----    type(point_list_type), intent(in out) :: pl  !  In Out-> Type with allocatable components
    !!----
    !!----     De-allocation of an objet of type point_list_type
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Deallocate_Point_List(Pl)
       !---- Arguments ----!
       type(point_list_type), intent(in out) :: pl

       if (allocated(pl%nam) ) deallocate(pl%nam)
       if (allocated(pl%p) )   deallocate(pl%p)
       if (allocated(pl%x) )   deallocate(pl%x)

       return
    End Subroutine Deallocate_Point_List

    !!----
    !!---- Subroutine Distance_and_Sigma(Cellp,DerM,x0,x1,s0,s1,dis,s)
    !!----    Type(Crystal_Cell_Type),         intent(in)  :: Cellp         ! Cell object
    !!----    real(kind=cp), dimension(3,3,6), intent(in)  :: DerM          ! Matrix of derivatives of Cellp%Cr_Orth_cel
    !!----    real(kind=cp), dimension(3),     intent(in)  :: x0,x1,s0,s1   ! Two points in fractional coordinates and sigmas
    !!----    real(kind=cp),                   intent(out) :: dis,s         ! Distance and sigma
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Distance_and_Sigma(Cellp,DerM,x0,x1,s0,s1,dis,s)
       !---- Arguments ----!
       Type(Crystal_Cell_Type),         intent(in)  :: Cellp         ! Cell object
       real(kind=cp), dimension(3,3,6), intent(in)  :: DerM          ! Matrix of derivatives of Cellp%Cr_Orth_cel
       real(kind=cp), dimension(3),     intent(in)  :: x0,x1,s0,s1   ! Two points in fractional coordinates and sigmas
       real(kind=cp),                   intent(out) :: dis,s         ! Distance and sigma

       !---- Local variables ----!
       integer                     :: i
       real(kind=cp), dimension(3) :: xc,xf
       real(kind=cp), dimension(6) :: dc,df

       xf=x1-x0
       xc = matmul(cellp%Cr_Orth_cel,xf)
       dis=sqrt(dot_product(xc,xc))
       do i=1,6
          dc(i) = dot_product(xc,matmul(DerM(:,:,i),xf))
       end do
       do i=1,3
          df(i) = dot_product(xc,Cellp%Cr_Orth_cel(:,i))
       end do
       df(4:6) =-df(1:3)
       s=0.0
       do i=1,3
          s = s + (dc(i)*Cellp%cell_std(i))**2
          s = s + (dc(i+3)*Cellp%ang_std(i)*to_rad)**2
          s = s + (df(i)*s1(i))**2 + (df(i+3)*s0(i))**2
       end do
       s=sqrt(s)/dis

       return
    End Subroutine Distance_and_Sigma

    !!----
    !!----  Subroutine Get_Anglen_Axis_From_RotMat(R,axis,angle)
    !!----    real(kind=cp), dimension(3,3), intent(in) :: R             !Input orthogonal matrix
    !!----    real(kind=cp), dimension(3),   intent(out):: axis          !Non normalized rotation axis
    !!----    real(kind=cp),                 intent(out):: angle         !Angle of rotation
    !!----
    !!----  Subroutine to obtain the axis and angle of rotation corresponding to
    !!----  an input orthogonal matrix. A Cartesian frame is assumed
    !!----
    !!---- Update: January - 2011
    !!----
    Subroutine Get_Anglen_Axis_From_RotMat(R,axis,angle)
      Real(kind=cp), dimension(3,3), intent(in) :: R
      Real(kind=cp), dimension(3),   intent(out):: axis
      Real(kind=cp),                 intent(out):: angle
      !--- Local variables ---!
      Real(kind=cp) :: va

      va=(R(1,1)+R(2,2)+R(3,3)-1.0_cp)*0.5_cp
      if(va < -1.0_cp) va=-1.0_cp
      if(va >  1.0_cp) va= 1.0_cp
      angle= acosd(va)
      if(abs(abs(angle)-180.0_cp) < epsi) then
         axis= (/                sqrt(R(1,1)+1.0_cp), &
                sign(1.0_cp,R(1,2))*sqrt(R(2,2)+1.0_cp), &
                sign(1.0_cp,R(1,3))*sqrt(R(3,3)+1.0_cp) /)
      else
         axis= (/  R(2,3)-R(3,2), &
                   R(3,1)-R(1,3), &
                   R(1,2)-R(2,1) /)
      end if
      return
    End Subroutine Get_Anglen_Axis_From_RotMat

    !!----
    !!----  Subroutine Get_Euler_From_Fract(X1,X2,X3,Mt,Phi,Theta,Chi,Eum,Code)
    !!----    real(kind=cp),           dimension(3),   intent (in) :: x1,x2,x3
    !!----    real(kind=cp),           dimension(3,3), intent (in) :: M !Matrix transforming to Cartesian coordinates
    !!----    real(kind=cp),                           intent(out) :: theta,phi,chi
    !!----    real(kind=cp), optional, dimension(3,3), intent(out) :: EuM
    !!----    character(len=*), optional,              intent (in) :: Code
    !!----
    !!----  Subroutine to obtain the Euler angles (2nd setting) of a Cartesian frame having
    !!----  as origin the point x3, the z-axis along x1-x3 and the "xz" plane coincident with
    !!----  the plane generated by the two vectors (x2-x3,x1-x3). The
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Get_Euler_From_Fract(X1,X2,X3,Mt,Phi,Theta,Chi,Eum,Code)
       !---- Arguments ----!
       real(kind=cp),           dimension(3),   intent (in) :: x1,x2,x3
       real(kind=cp),           dimension(3,3), intent (in) :: Mt
       real(kind=cp),                           intent(out) :: theta,phi,chi
       real(kind=cp), optional, dimension(3,3), intent(out) :: EuM
       character(len=*), optional,              intent (in) :: Code

       !---- Local variables ----!
       real(kind=cp), dimension(3)   :: u,v,w
       real(kind=cp), dimension(3,3) :: rot

!  U = ( cosPhi cosTheta cosChi - sinPhi sinChi,   sinPhi cosTheta cosChi+cosPhi sinChi,  -sinTheta cosChi)
!  V = (-sinPhi cosChi   - cosPhi cosTheta sinChi, cosPhi cosChi -sinPhi cosTheta sinChi,  sinTheta sinChi)
!  W = ( cosPhi sinTheta, sinPhi sinTheta,  cosTheta)
!
!     This corresponds to Euler angles defined in the following way:
!
!     In the starting position the cartesian frame (u,v,w) coincides with the crystallographic
!     cartesian frame (e1//a, e2 in the a-b plane and e3= e1 x e2). First a rotation Chi around
!     the e3 axis is applied, then a rotation Theta around the e2 axis and finally a rotation Phi
!     around e3. The total rotation matrix is
!
!          R(Phi,Theta,Chi) = R(e3,Phi) R(e2,Theta) R(e3,Chi) = [[ u, v, w]]
!
!     The columns of the active rotation matrix are the components of the unitary vectors u,v,w.

       w=matmul(Mt,x1-x3)
       w=w/sqrt(dot_product(w,w))
       u=matmul(Mt,x2-x3)
       u=u/sqrt(dot_product(u,u))
       v=cross_product(w,u)
       v=v/sqrt(dot_product(v,v))
       u=cross_product(v,w) !already normalized
       rot(:,1)=u; rot(:,2)=v;  rot(:,3)=w  !Matrix Rot ([u,v,w] columns)
       if (present(EuM)) EuM=rot
       if (present(Code)) then
          call get_PhiTheChi(rot,phi,theta,chi,Code)
       else
          call get_PhiTheChi(rot,phi,theta,chi)
       end if

       return
    End Subroutine Get_Euler_From_Fract

    !!---- Subroutine Get_Matrix_moving_v_to_u(v,u,R,w,ang)
    !!----   real(kind=cp), dimension(3),           intent(in)  :: v,u   !Starting and final vectors
    !!----   real(kind=cp), dimension(3,3),         intent(out) :: R     !Rotation matrix moving v to u:  u=Rv
    !!----   real(kind=cp), optional,               intent(out) :: ang   !angle between the two vectors
    !!----   real(kind=cp), optional,dimension(3),  intent(out) :: w     !axis normal to plane of the two vectors
    !!----
    !!----   Subroutine to get the orthogonal matrix that rotates a vector v
    !!----   to orient it along the vector u. Makes use of Cross_Product and
    !!----   Rot_matrix (Gibbs matrix)
    !!----
    !!----    Created: February 2010 (JRC)
    !!----    Updated: March 2013 (JRC)
    !!----
    !!
    Subroutine Get_Matrix_moving_v_to_u(v,u,R,w,ang)
      real(kind=cp), dimension(3),           intent(in)  :: v,u
      real(kind=cp), dimension(3,3),         intent(out) :: R
      real(kind=cp), optional,               intent(out) :: ang
      real(kind=cp), optional,dimension(3),  intent(out) :: w
      !--- Local variables ---!
      integer                        :: i,iu,iv
      real(kind=cp), parameter       :: ep=1.0e-5_cp
      real(kind=cp)                  :: mv,mu,mvu,phi,c
      logical                        :: co_linear
      real(kind=cp), dimension(3)    :: vu
      integer, dimension(1)          :: im
      real(kind=cp), parameter, dimension(3,3):: ident=reshape((/1.0_cp,0.0_cp,0.0_cp, &
                                                                 0.0_cp,1.0_cp,0.0_cp, &
                                                                 0.0_cp,0.0_cp,1.0_cp/),(/3,3/))

      if(present(ang)) ang=0.0
      if(present(w))   w=0.0
      !First determine if the two input vectors are co-linear
      im=maxloc(abs(v))
      iv=im(1)
      im=maxloc(abs(u))
      iu=im(1)
      co_linear=.true.
      if(iu == iv) then ! may be co-linear
        if(abs(u(iu)) > ep) then
          c=v(iv)/u(iu)
          do i=1,3
            if(abs( v(i)-c*u(i) ) > ep ) then
               co_linear=.false.
               exit
            end if
          end do
        end if
      else
        co_linear=.false.
      end if
      if(co_linear) then
        mvu=v(iv)*u(iu)
        if(mvu < 0.0) then   !opposed vectors
          R=-ident
        else                 !parallel vectors
          R=ident
        end if
      else
        ! non co-linear
        vu=Cross_Product(v,u)      !Rotation axis
        mv=sqrt(dot_product(v,v))
        mu=sqrt(dot_product(u,u))
        phi=dot_product(u,v)/mv/mu
        phi=acosd(phi)        !Angle between the two input vectors
        R=Rot_matrix(vu,phi)  !Gibbs matrix
        if(present(ang)) ang=phi
        if(present(w)) w=vu
      end if
      return
    End Subroutine Get_Matrix_moving_v_to_u

    !!----
    !!---- Subroutine Get_OmegaChiPhi(Mt,Omega,Chi,Phi,Code)
    !!----    real(kind=cp), dimension(3,3),intent(in)  :: Mt
    !!----    real(kind=cp),                intent(out) :: Omega
    !!----    real(kind=cp),                intent(out) :: Chi
    !!----    real(kind=cp),                intent(out) :: Phi
    !!----    character(len=*), optional,   intent(in)  :: Code
    !!----
    !!----    Calculate the Euler Angles corresponding to an orthogonal matrix
    !!----    The definition of the Euler angles in this case correspond to the
    !!----    rotation matrix of Busing and Levy for diffractometry obtained from
    !!----    the composition of a rotation around z of angle Phi, followed by a
    !!----    rotation of angle Chi around the y-axis and a subsequent rotation of angle
    !!----    Omega around z.
    !!----    The matrix is supposed to be of the form: M = Rz(Omega).Ry(Chi).Rz(Phi)
    !!----    If Code =="R" or not present then the output angles are provided in radians.
    !!----    If Code =="D" then the output angles are provided in degrees.
    !!----    A checking of the input matrix is given before calculating the angles.
    !!----    The user must check the logical variable "ERR_RotMat" after calling this
    !!----    subroutine. If ERR_RotMat=.true. it means that the input matrix is not orthogonal.
    !!----    The obtained rotations should be interpreted as changes of reference systems, the
    !!----    angles correspond to the motor settings to put a reciprocal vector in Cartesian
    !!----    coordinates w.r.t. the L-system (all angles equal to zero) in the position given
    !!----    by the active rotation matrix Mt:  z4= Mt z1.
    !!----
    !!---- Updated: March - 2013
    !!
    Subroutine Get_OmegaChiPhi(Mt,Omega,Chi,Phi,Code)  !Conventional Euler angles of diffractometry
       !---- Arguments ----!
       real(kind=cp), dimension(3,3),intent(in)  :: Mt
       real(kind=cp),                intent(out) :: Omega
       real(kind=cp),                intent(out) :: Chi
       real(kind=cp),                intent(out) :: Phi
       character(len=*), optional,   intent(in)  :: Code

       !---- Local Variables ----!
       real(kind=cp), dimension(3,3):: MTT
       real(kind=cp), parameter, dimension(3,3) :: &
                      identity = reshape ( (/1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0/),(/3,3/))

       MTT=transpose(Mt)
       MTT=matmul(MTT,Mt)-identity
       if (sum(abs(MTT)) > 5.0*eps) then
          ERR_Geom=.true.
          ERR_Geom_Mess=" Error in Get_OmegaChiPhi ... the input matrix is not orthogonal! "
          return
       end if
       if (abs(Mt(3,3)-1.0) < eps) then  !M(3,3)=cos(Chi)=1
          Chi=0.0
          Omega=0.0                       ! Omega and Phi have the same axis, we select Omega=0
          !Phi=acos(Mt(1,1))              ! M(1,1)=cos(Omega)cos(Chi)cos(Phi)-sin(Omega)sin(Phi)
          Phi=atan2(Mt(1,2),Mt(1,1))      ! M(1,2)=cos(Omega)cos(Chi)sin(Phi)+sin(Omega)cos(Phi)
       else if(abs(Mt(3,3)+1.0) < eps) then  !M(3,3)=cos(Chi)=-1
          Chi=pi
          Omega=0.0                       ! Omega and Phi have the same axis, we select Omega=0
          !Phi=acos(-Mt(1,1))             ! We use also the elements (11) and (12)
          Phi=atan2(-Mt(1,2),Mt(1,1))
       else
          !Chi=acos(Mt(3,3))  !Better use the relation below (In BL there is an error in eqn 48 for omega)
          Omega=atan2(-Mt(2,3),Mt(1,3))       !M(1,3)=  cos(Omega)sin(Chi)   M(2,3)= -sin(Omega)sin(Chi)
          Phi=atan2(-Mt(3,2),-Mt(3,1))        !M(3,1)= -sin(Chi)cos(Phi)     M(3,2)= -sin(Chi)sin(Phi)
          Chi=atan2( Sqrt(Mt(3,1)*Mt(3,1)+Mt(3,2)*Mt(3,2)), Mt(3,3) )
       end if
       if (present(Code)) then
          if (code(1:1)=="D" .or. code(1:1)=="d") then
             Phi=Phi*to_deg
             Omega=Omega*to_deg
             Chi=Chi*to_deg
          end if
       end if

       return
    End Subroutine Get_OmegaChiPhi

    !!----
    !!---- Subroutine Get_PhiTheChi(Mt,Phi,Theta,Chi,Code)
    !!----    real(kind=cp), dimension(3,3),intent(in)  :: Mt
    !!----    real(kind=cp),                intent(out) :: Phi
    !!----    real(kind=cp),                intent(out) :: Theta
    !!----    real(kind=cp),                intent(out) :: Chi
    !!----    character(len=*), optional,   intent(in)  :: Code
    !!----
    !!----    Calculate the Euler Angles corresponding to an orthogonal matrix
    !!----    The definition of the Euler angles in this case correspond to the
    !!----    active rotation matrix obtained from the composition of a rotation
    !!----    around z of angle Chi, followed by a rotation of angle Theta
    !!----    around the y-axis and a subsequent rotation of angle Phi around z.
    !!----    The matrix is supposed to be of the form: M = Rz(Phi).Ry(Theta).Rz(Chi)
    !!----    If Code =="R" or not present then the output angles are provided in radians.
    !!----    If Code =="D" then the output angles are provided in degrees.
    !!----    A checking of the input matrix is given before calculating the angles.
    !!----    The user must check the logical variable "err_geom" after calling this
    !!----    subroutine. If err_geom=.true. it means that the input matrix is not orthogonal.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Get_PhiTheChi(Mt,Phi,Theta,Chi,Code)
       !---- Arguments ----!
       real(kind=cp), dimension(3,3),intent(in)  :: Mt
       real(kind=cp),                intent(out) :: Phi
       real(kind=cp),                intent(out) :: Theta
       real(kind=cp),                intent(out) :: Chi
       character(len=*), optional,   intent(in)  :: Code

       !---- Local Variables ----!
       real(kind=cp), dimension(3,3):: MTT
       real(kind=cp), parameter, dimension(3,3) :: &
                      identity = reshape ( (/1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0/),(/3,3/))

       MTT=transpose(Mt)
       MTT=matmul(MTT,Mt)-identity
       if (sum(abs(MTT)) > 5.0*eps) then
          err_geom=.true.
          ERR_Geom_Mess=" Error in Get_PhiTheChi ... the input matrix is not orthogonal! "
          return
       end if
       if (abs(Mt(3,3)-1.0) < eps) then  !M(3,3)=cos(Theta)
          Theta=0.0
          Phi=0.0
          Chi=acos(Mt(1,1))               !M(1,1)=cos(Phi)cos(Theta)cos(Chi)-sin(Phi)sin(Chi)
       else if(abs(Mt(3,3)+1.0) < eps) then
          Theta=pi
          Phi=0.0
          Chi=acos(-Mt(1,1))
       else
          Theta=acos(Mt(3,3))
          Phi=atan2(Mt(2,3),Mt(1,3))     !M(1,3)=cos(Phi)sin(Theta)  M(2,3)=sin(phi)sin(Theta)
          Chi=atan2(Mt(3,2),-Mt(3,1))    !M(3,1)= -sin(Theta)cos(Chi)   M(3,2)= sin(Theta)sin(Chi)
       end if
       if (present(Code)) then
          if (code(1:1)=="D" .or. code(1:1)=="d") then
             Phi=Phi*to_deg
             Theta=Theta*to_deg
             Chi=Chi*to_deg
          end if
       end if

       return
    End Subroutine Get_PhiTheChi

    !!----
    !!---- Subroutine Get_Transf_List(Trans,Ox,Pl,Npl,Ifail)
    !!----   real(kind=cp), dimension(3,3), intent(in)     :: trans   !Matrix transforming the basis
    !!----   real(kind=cp), dimension(3  ), intent(in)     :: ox      !Coordinates of origin of the new basis
    !!----   type(point_list_type),         intent(in)     :: pl      !Input List of points
    !!----   type(point_list_type),         intent(in out) :: npl     !Output list of transformed points
    !!----   integer,                       intent(out)    :: ifail   !If ifail/=0 matrix inversion failed
    !!----
    !!----  Subroutine to get the fractional coordinates of the points of the input list "pl" in the
    !!----  new transformed cell ( a'= trans a) displaced to the new origing "ox". The coordinates
    !!----  are generated using only lattice translations. All coordinates are reduced to be
    !!----  between 0.0 and 1.0, so that  0.0 <= x,y,z < 1.0
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Get_Transf_List(trans,ox,pl,npl,ifail)
       !---- Arguments ----!
       real(kind=cp),         dimension(3,3), intent(in)     :: trans
       real(kind=cp),         dimension(3  ), intent(in)     :: ox
       type(point_list_type),                 intent(in)     :: pl
       type(point_list_type),                 intent(in out) :: npl
       integer,                               intent(out)    :: ifail

       !---- local variables ----!
       integer                       :: i,j,ia,ib,ic,nat,mm
       integer, dimension(3)         :: mini,maxi
       real(kind=cp), dimension(7,3) :: vecpar
       real(kind=cp), dimension(3,3) :: si
       real(kind=cp), dimension(3  ) :: xx, xxn,v

       ifail=0
       call matrix_inverse(trans,si,ifail)
       if (ifail == 1) return

       !----  Construction of the 7 vertices of the new cell
       !----  1:a, 2:b, 3:c, 4:a+b, 5:a+c, 6:b+c 7:a+b+c
       do j=1,3
          do i=1,3
             vecpar(i,j)=trans(i,j)
          end do
          vecpar(4,j)=trans(1,j)+trans(2,j)
          vecpar(5,j)=trans(1,j)+trans(3,j)
          vecpar(6,j)=trans(2,j)+trans(3,j)
          vecpar(7,j)=trans(1,j)+trans(2,j)+trans(3,j)
       end do

       !---- Exploration of the vertex matrix
       mini(:)=1000
       maxi(:)=-1000
       do j=1,3
          do i=1,7
             if (vecpar(i,j) < mini(j)) mini(j)=nint(min(vecpar(i,j),0.0_cp))
             if (vecpar(i,j) > maxi(j)) maxi(j)=nint(max(vecpar(i,j),1.0_cp))
          end do
       end do

       !
       !   Explore the region  a-> min(1)---max(1)  where atoms will be generated
       !                       b-> min(2)---max(2)
       !                       c-> min(3)---max(3)
       !   and select those belonging to the interior of the new cell before
       !   translation to the new origin.
       !   set the translation to the new origin, put the atoms inside the new
       !   unit cell and, finally, print atoms coordinates
       !
       nat=0
       do mm=1,pl%np
          do ia=mini(1),maxi(1)
             xx(1)=pl%x(1,mm)+real(ia)
             do ib=mini(2),maxi(2)
                xx(2)=pl%x(2,mm)+real(ib)
                do_ic: do ic=mini(3),maxi(3)
                   xx(3)=pl%x(3,mm)+real(ic)
                   xxn=matmul(xx-ox,si)
                   xxn=Modulo_Lat(xxn)
                   do i=nat,1,-1
                      v=npl%x(:,i)-xxn(:)
                      if (Lattice_trans(v,"P") ) cycle do_ic
                   end do
                   nat=nat+1
                   npl%x(:,nat)= xxn
                   if ( nat < 10) then
                      write(unit=npl%nam(nat),fmt="(a,i1)") trim(pl%nam(mm))//"_",nat
                   else if( nat < 100) then
                      write(unit=npl%nam(nat),fmt="(a,i2)") trim(pl%nam(mm))//"_",nat
                   else
                      write(unit=npl%nam(nat),fmt="(a,i3)") trim(pl%nam(mm))//"_",nat
                   end if
                end do do_ic
             end do
          end do
       end do
       npl%np=nat

       return
    End Subroutine Get_Transf_List

    !!----
    !!---- Subroutine Init_Err_Geom()
    !!----
    !!----    Initialize the errors flags in CFML_Geometry_Calc
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Init_Err_Geom()

       err_geom=.false.
       ERR_Geom_Mess=" "

       return
    End Subroutine Init_Err_Geom

    !!----
    !!---- Subroutine P1_Dist(Dmax, Cell, Spg, Ac, Lun)
    !!----    real(kind=cp),            intent(in)    :: dmax      !  In -> Max. Distance to be calculated
    !!----    type (Crystal_cell_Type), intent(in)    :: Cell      !  In -> Object of Crystal_cell_Type
    !!----    type (Space_Group_Type),  intent(in)    :: SpG       !  In -> Object of Space_Group_Type
    !!----    type (Atoms_Cell_Type),   intent(in out):: Ac        !  In -> Object of Atoms_Cell_Type
    !!----                                                           Out -> Updated Object of Atoms_Cell_Type
    !!----    integer,optional,         intent(in)    :: lun       !  In -> Logical Unit for writing
    !!----
    !!----    Subroutine calculate distances, below the prescribed distances "dmax",
    !!----    without standard deviations. No symmetry is applied: only lattice translations.
    !!----    Need as input the objects "Cell" (of type Crystal_cell_type), "SpG" (of type Space_Group_Type)
    !!----    and "Ac" (or type Atoms_Cell). Complete the construction of Ac.
    !!----    Control for error is present.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine P1_Dist(Dmax, Cell, Spg, Ac, Lun)
       !---- Arguments ----!
       real(kind=cp),            intent(in)       :: dmax
       type (Crystal_cell_Type), intent(in)       :: Cell
       type (Space_Group_Type),  intent(in)       :: SpG
       type (Atoms_Cell_Type),   intent(in out)   :: Ac
       integer, optional,        intent(in)       :: lun

       !---- Local Variables ----!
       logical                                :: iprint
       character(len=6 )                      :: nam,nam1
       character(len=40)                      :: transla
       character(len=90)                      :: form1,form2="(a,2i4,a,a,a,a,a,f10.4,a,t62,3F8.4)"
       integer                                :: i,k,lk,i1,i2,i3,jl,nn,L,inew,ne,id
       integer, dimension(3)                  :: ic1,ic2
       integer, dimension(Ac%nat,Ac%nat)      :: mn  !neighbouring matrix
       real(kind=cp)                          :: T,dd
       real(kind=cp), dimension(3)            :: xx,x1,xo,Tn,xr, QD
       real(kind=cp), dimension(3,Ac%nat*Ac%nat*spg%multip) :: u

       iprint=.false.
       if (present(lun)) then
          if (lun > 0) iprint=.true.
       end if
       call init_err_geom()
       id=3*nint(0.74048*(dmax/1.1)**3)

       qd(:)=1.0/cell%rcell(:)
       ic2(:)= nint(dmax/cell%cell(:)+3.0)
       ic1(:)=-ic2(:)
       mn(:,:) = 0
       inew=0
       do i=1,ac%nat
          xo(:)=Ac%xyz(:,i)
          nam= Ac%noms(i)
          if (iprint) then
             write(unit=lun,fmt="(/,/,a)")"    -------------------------------------------------------------------"
             write(unit=lun,fmt="(a,f8.4,a,a,3f8.4)")   &
                       "    Distances less than",dmax,"  to atom: ",nam, xo(:)
             write(unit=lun,fmt="(a,/,/)")"    -------------------------------------------------------------------"
             write(unit=lun,fmt="(/,/,a,/,/)") &
                       " Orig. extr.                    Distance     tx   ty   tz       x_ext   y_ext   z_ext"
          end if
          ne=0
          do k=1,Ac%nat
             lk=1
             u(:,lk)=xo(:)
             xx(:)=Ac%xyz(:,k)
             nam1= Ac%noms(k)
             do i1=ic1(1),ic2(1)
                do i2=ic1(2),ic2(2)
                   do i3=ic1(3),ic2(3)
                      do_jl:do jl=1,Spg%NumLat
                         Tn(:)=(/real(i1),real(i2),real(i3)/)+Spg%Latt_trans(:,jl)
                         x1(:)=xx(:)+tn(:)
                         do l=1,3
                            t=abs(x1(l)-xo(l))*qd(l)
                            if (t > dmax) cycle do_jl
                         end do
                         do nn=1,lk
                            if (sum(abs(u(:,nn)-x1(:)))  <= epsi) cycle do_jl
                         end do
                         xr = matmul(cell%cr_orth_cel,x1-xo)
                         dd=sqrt(dot_product(xr,xr))
                         if (dd > dmax .or. dd < 0.001) cycle
                         lk=lk+1
                         u(:,lk)=x1(:)
                         call Frac_Trans_1Dig(tn,transla)
                         if (iprint) write(unit=lun,fmt=form2)" ",i,k,"   (",nam,")-(",nam1,"):",dd,"  "//transla,x1(:)
                         mn(i,k)=mn(i,k)+1
                         ne=ne+1
                         IF (ne > id) THEN
                            err_geom=.true.
                            ERR_Geom_Mess="Too many connected atoms! in sub. P1_dist"
                            return
                         END IF
                         Ac%neighb_atom(ne,i)=k    !Pointer to the number of atom connected to i
                         Ac%distance   (ne,i)=dd   !Corresponding distance
                         Ac%trans(:,ne,i)=tn(:)    !corresponding lattice translation
                         do nn=1,inew
                            if (abs(dd-Ac%ddist(nn)) <= epsi) then
                               if (equiv_atm(nam,nam1,Ac%ddlab(nn)))  cycle do_jl
                            end if
                         end do
                         inew=inew+1
                         Ac%ddist(inew)=dd
                         Ac%ddlab(inew)=wrt_lab(nam,nam1)
                      end do do_jl
                   end do !i3
                end do !i2
             end do !i1
          end do !k
          Ac%neighb(i)=ne
       end do !i
       Ac%ndist=inew
       if (iprint) then
          write(unit=lun,fmt="(/,/,a)") " -------------------"
          write(unit=lun,fmt="(a)"  )   " Neighbouring matrix"
          write(unit=lun,fmt="(a)")     " -------------------"
          write(unit=lun,fmt="(a)")
          write(unit=form1,fmt="(a,i4,a)") "(a,",Ac%nat,"i3)"
          write(unit=lun,fmt=form1)"     ",(i,i=1,Ac%nat)
          write(unit=lun,fmt="(a)")
          write(unit=form1,fmt="(a,i4,a)") "(i3,a,",Ac%nat,"i3)"
          do i=1,ac%nat
             write(unit=lun,fmt=form1) i,"  ",(mn(i,k),k=1,Ac%nat)
          end do
          write(unit=lun,fmt="(a,/,/,/)")
       end if

       return
    End Subroutine P1_Dist

    !!----
    !!---- Subroutine Print_Distances(Lun, Dmax, Cell, Spg, A)
    !!----    integer,                  intent(in)   :: lun    !  In -> Logical Unit for writing
    !!----    real(kind=cp),            intent(in)   :: dmax   !  In -> Max. Distance to be calculated
    !!----    type (Crystal_cell_Type), intent(in)   :: Cell   !  In -> Object of Crystal_cell_Type
    !!----    type (Space_Group_Type),  intent(in)   :: SpG    !  In -> Object of Space_Group_Type
    !!----    type (atom_list_type),   intent(in)   :: A      !  In -> Object of atom_list_type
    !!----
    !!----    Subroutine to print distances, below the prescribed distances
    !!----    "dmax", without standard deviations.
    !!----    Need as input the objects "Cell" (of type Crystal_cell_type), "SpG"
    !!----    (of type Space_Group_type) and "A" (or type atom_list_type, that should be
    !!----    allocated in the calling program).
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Print_Distances(Lun, Dmax, Cell, Spg, A)
       !-- Arguments --!
       integer,                  intent(in)   :: lun
       real(kind=cp),            intent(in)   :: dmax
       type (Crystal_cell_Type), intent(in)   :: Cell
       type (Space_Group_Type),  intent(in)   :: SpG
       type (atom_list_type),    intent(in)   :: A

       !---- Local Variables ----!
       integer                           :: i,j,k,lk,i1,i2,i3,jl,npeq,nn,L,nlines
       character(len=80), dimension(12)  :: texto=" "
       character(len=5 )                 :: nam,nam1
       character(len=40)                 :: transla
       character(len=54)                 :: form2="(a,3i4,a,a,a,a,a,f10.4,a,t66,3F8.4)" !&
                                            !"("" "",3I4,""  ("",a,"")-("",a,""):"",f9.4,""   "",a,""  "",3F8.4)"
       integer,          dimension(3)    :: ic1,ic2
       real(kind=cp),    dimension(3)    :: xx,x1,xo,Tn,xr, QD
       real(kind=cp)                     :: T,dd
       real(kind=cp), dimension(3,A%Natoms*Spg%multip) :: uu

       qd(:)=1.0/cell%rcell(:)
       ic2(:)= nint(dmax/cell%cell(:)+1.0)
       ic1(:)=-ic2(:)
       npeq=spg%numops

       if (Spg%Centred == 2) then
          npeq=2*npeq
          write(unit=lun,fmt="(a)")" => Symmetry operators combined with inversion centre:"
          nlines=1
          do i=SpG%NumOps+1,npeq
             if (mod(i,2) == 0) then
                write(unit=texto(nlines)(36:70),fmt="(a,i2,a,a)") &
                                           " => SYMM(",i,"): ",trim(SpG%SymopSymb(i))
                nlines=nlines+1
             else
                write(unit=texto(nlines)( 1:34),fmt="(a,i2,a,a)")  &
                                           " => SYMM(",i,"): ",trim(SpG%SymopSymb(i))
             end if
          end do
          do i=1,min(nlines,12)
             write(unit=lun,fmt="(a)") texto(i)
          end do
       end if

       do i=1,A%natoms
          nam=a%atom(i)%lab
          xo(:)=a%atom(i)%x(:)
          write(unit=lun,fmt="(/,/,a)")"    -------------------------------------------------------------------"
          write(unit=lun,fmt="(a,f8.4,a,a,3f8.4)")   &
                    "    Distances less than",dmax,"  to atom: ",nam, xo(:)
          write(unit=lun,fmt="(a,/,/)")"    -------------------------------------------------------------------"
          write(unit=lun,fmt="(/,/,a,/,/)") &
                    " Orig. extr. p.equiv.           Distance     tx   ty   tz       x_ext   y_ext   z_ext"
          do k=1,a%natoms
             lk=1
             uu(:,lk)=xo(:)
             nam1=a%atom(k)%lab
             do j=1,npeq
                xx=ApplySO(Spg%SymOp(j),a%atom(k)%x)
                do i1=ic1(1),ic2(1)
                   do i2=ic1(2),ic2(2)
                      do i3=ic1(3),ic2(3)
                         do_jl:do jl=1,Spg%NumLat
                            Tn(:)=real((/i1,i2,i3/))+Spg%Latt_trans(:,jl)
                            x1(:)=xx(:)+tn(:)
                            do l=1,3
                               t=abs(x1(l)-xo(l))*qd(l)
                               if (t > dmax) cycle do_jl
                            end do
                            do nn=1,lk
                               if (sum(abs(uu(:,nn)-x1(:)))  <= epsi) cycle do_jl
                            end do
                            xr = matmul(cell%cr_orth_cel,x1-xo)
                            dd=sqrt(dot_product(xr,xr))
                            if (dd > dmax .or. dd < 0.001) cycle
                            lk=lk+1
                            uu(:,lk)=x1(:)
                            call Frac_Trans_1Dig(tn,transla)
                            write(unit=lun,fmt=form2)" ",i,k,j,"   (",nam,")-(",nam1,"):",dd,"  "//transla,x1(:)
                            !write(unit=lun,fmt=form2)i,k,j,nam ,nam1,dd,transla,x1(:)
                                        !    "("" "",3I4,""  ("",a,"")-("",a,""):"",f9.4,""   "",a,""  "",3F8.4)"
                         end do do_jl
                      end do !i3
                   end do !i2
                end do !i1
             end do !j
          end do !k
       end do !i

       return
    End Subroutine Print_Distances

    !!---- Subroutine Set_New_AsymUnit(SpGn,Ate,Mat,orig,A_n,matkind,debug)
    !!----    type (Space_Group_Type) ,      intent(in    ) :: SpGn    !New space group that has been previously set
    !!----    type (Atom_Equiv_List_Type),   intent(in    ) :: Ate     !In old group
    !!----    real(kind=cp), dimension (3,3),intent(in    ) :: Mat     !Transformation matrix from the old to the new setting
    !!----    real(kind=cp), dimension (  3),intent(in    ) :: orig    !Displacement of the origin in the old setting
    !!----    type (Atom_list_Type),         intent(in out) :: A_n     !New atom list
    !!----    character (len=*), optional,   intent(in    ) :: matkind !Kind of transformation matrix
    !!----    character (len=*), optional,   intent(in    ) :: debug
    !!----
    !!----    Updated: January 2014 (JRC)
    !!----
    !!----
    Subroutine Set_New_AsymUnit(SpGn,Ate,Mat,orig,A_n,matkind,debug)
       type (Space_Group_Type) ,      intent(in    ) :: SpGn
       type (Atom_Equiv_List_Type),   intent(in    ) :: Ate !In old group
       real(kind=cp), dimension (3,3),intent(in    ) :: Mat
       real(kind=cp), dimension (  3),intent(in    ) :: orig
       type (Atom_list_Type),         intent(out   ) :: A_n
       character (len=*), optional,   intent(in    ) :: matkind
       character (len=*), optional,   intent(in    ) :: debug
       ! Local variables
       integer                           :: i,j,k,m,ifail,L,n,Ls,ip,L1
       integer                           :: i1,i2,i3,maxa,maxp,maxm,mult
       real(kind=cp), dimension (3,3)    :: S,Sinv
       real(kind=cp)                     :: determ
       logical                           :: newp,fail
       real(kind=cp), dimension (  3)    :: pos
       real(kind=cp), dimension (3,192)  :: orb
       type(point_list_type)             :: pl
       type (Atom_list_Type)             :: A
       character(len=*),parameter,dimension(26) :: let=(/"a","b","c","d","e","f","g","h", &
          "i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z"/)
       real(kind=cp), allocatable, dimension (:,:) :: vec
       integer,parameter         :: lu=93
       real(kind=cp), parameter  :: epsi = 0.002

       if(present(matkind)) then
        if(matkind(1:2) == "it" .or. matkind(1:2) == "IT" ) then
          S=Mat             !Atoms positions X'=inv(Mat) (X-O)
        else
          S=transpose(Mat)  !Atoms positions X'=inv(MatT)(X-O)
        end if
       else
          S=transpose(Mat)
       end if
       call matrix_inverse(S,Sinv,ifail) !Atoms positions X'= Sinv X
       if (ifail /= 0) then
         err_geom=.true.
         err_geom_Mess= "Inversion Matrix Failed on: Change_Setting_SG"
         return
       end if
       if(present(debug)) then
         open(unit=lu,file="similar_debug.lis",status="replace",action="write")
         write(unit=lu,fmt="(a)")  "  Debugging SIMILAR calculations  "
         write(unit=lu,fmt="(a/)") "  ============================== "
       end if

       determ=determ_a(S)
       determ=abs(determ)

       m=0
       if(determ > 1.0001) then  !Generate indices for lattice translations to be applied
         i1=max(nint(maxval(abs(S(:,1))))-1,1)
         i2=max(nint(maxval(abs(S(:,2))))-1,1)
         i3=max(nint(maxval(abs(S(:,3))))-1,1)
         allocate(vec(3,(i1+1)*(i2+1)*(i3+1)))
         do i=0,i1
          do j=0,i2
           do k=0,i3
            if(i==0 .and.j==0 .and. k==0) cycle
            m=m+1
            vec(:,m) = real((/i,j,k/))
            !write(*,*) "  vect: ",m," :",vec(:,m)
           end do
          end do
         end do
       end if

       maxm=m    !maximum number of translations to be applied before changing the basis
       maxa=maxval(Ate%Atm(:)%mult)  !highest multiplicity of the atom sequence list
       !Factor 2
       maxp=2*maxa*determ    !maximum multiplicity in the new cell of a particular atom type
       !write(*,*) " Allocating atoms_list and Point list for ", maxp*Ate%nauas, " and ", maxp, " values"
       call Allocate_Atom_List(maxp*Ate%nauas,A)  !Atom list in the new cell, we must use "maxp"*Ate%nauas and not "maxa"
       Call Allocate_Point_List(maxp*Ate%nauas,Pl,Ifail)
       if(ifail /= 0) then
         !write(*,*) " Error allocating PL for ",maxp," values"
          err_geom=.true.
          write(unit=err_geom_Mess,fmt="(a,i8,a)")  " Error allocating PL for ",maxp," values"
          return
       end if
       Ls=0

       !write(*,*) " Allocating PL and A successful "

       do i=1,Ate%nauas
         Ls=Ls+1
         !Setting pl object
           ip=index(Ate%Atm(i)%Lab(1),"_")
           if(ip /= 0) then
              pl%nam(i)=Ate%Atm(i)%Lab(1)(1:ip-1)
           else
              pl%nam(i)=Ate%Atm(i)%Lab(1)
           end if
           pl%x=0.0
           pl%p=0.0
         !
         n=0
         do j=1,Ate%Atm(i)%mult
           n=n+1
           pos(:)    = Ate%Atm(i)%x(:,j)-orig(:)
           !write(*,*) " n=",n
           pl%x(:,n) = matmul(Sinv,pos)  !Complete list in new coordinate system of atoms of type i
           if(present(debug)) then
             write(unit=lu,fmt="(i4,2(a,3f8.4),a)") n,"  Atom: "//pl%nam(i)//" at (", &
                                         Ate%Atm(i)%x(:,j),") trasform to (",pl%x(:,n),")"
           end if
           pl%x(:,n) = Modulo_Lat(pl%x(:,n))  !Complete list in new coordinate system of atoms of type i
         end do
         if(determ > 1.0) then
          doj:do j=1,Ate%Atm(i)%mult
               do m=1,maxm
                 pos(:) = Ate%Atm(i)%x(:,j)-orig(:)+ vec(:,m)
                 pos(:) = Modulo_Lat(matmul(Sinv,pos))
                 newp=.true.
                 do k=1,n
                    if (sum(abs(pos(:)-pl%x(:,k))) < epsi) then
                       newp=.false.
                       exit
                    end if
                 end do
                 if (newp) then ! new position
                    n=n+1
                    !write(*,*) "  n=",n
                    pl%x(:,n) = pos(:)
                    if(present(debug)) then
                      write(unit=lu,fmt="(i4,2(a,3f8.4),a)") n,"  Atom: "//pl%nam(i)//" at (", &
                                                  Ate%Atm(i)%x(:,j),") trasform to (",pl%x(:,n),")"
                    end if
                    if(n == maxp) exit doj
                 end if
               end do
             end do doj
         end if

         pl%np=n
         A%atom(Ls)%Lab =pl%nam(i)
         A%atom(Ls)%x(:)=pl%x(:,1)
         !write(*,"(2i5,a,i5,a)") i,Ls, "  "//Ate%Atm(i)%Lab(1), Ate%Atm(i)%mult,"   "//A%atom(Ls)%Lab

         !Determine the number of independent orbits for this point
         call Set_Orbits_Inlist(Spgn,pl)
         L=1; L1=1
         do j=2,n
           if(pl%p(j) > L) then
            Ls=Ls+1
            A%atom(Ls)%x(:)=pl%x(:,j)
            !write(unit=let,fmt="(i3.3)") L
            A%atom(Ls)%Lab =trim(pl%nam(i))//let(L1)
           ! write(*,"(2i5,a,i5,a)") i,Ls, "  "//Ate%Atm(i)%Lab(1), Ate%Atm(i)%mult,"   "//A%atom(Ls)%Lab
            L=L+1
            L1=L1+1  !using a different counter for the label
            if(L1 > 26) L1=1 !re-start the labelling with the same letter
           end if
         end do

       end do  !i=1,Ate%nauas

       !write(*,*) "  Orbits correct"
       !write(*,*) "  Allocate_Atom_List for ",Ls," atoms"
       call Allocate_Atom_List(Ls,A_n,fail)
       if(fail) then
         if(present(debug)) then
          !write(*,*) "  Error on Allocate_Atom_List for ",A_n%natoms," atoms"
          write(unit=lu,fmt="(a,i4,a)") "  Error on Allocate_Atom_List for ",A_n%natoms," atoms"
         end if
       else
         ! write(*,*) "  Success on Allocate_Atom_List for ",A_n%natoms," atoms"
       end if
       do i=1,A_n%natoms
         A_n%atom(i)%x= A%atom(i)%x
         A_n%atom(i)%Lab= A%atom(i)%Lab
         call Get_Orbit(A_n%atom(i)%x,Spgn,Mult,orb)
         A_n%atom(i)%Mult=mult
         A_n%atom(i)%occ=real(mult)/real(Spgn%Multip)
       end do
       if(allocated(A%atom)) deallocate(A%atom)
       if(present(debug)) close(unit=lu)
       return
    End Subroutine Set_New_AsymUnit


    !!----
    !!---- Subroutine Set_Orbits_Inlist(Spg,Pl)
    !!----    type(space_group_type), intent(in)     :: SpG     !  In -> Space group
    !!----    type(point_list_type),  intent(in out) :: pl      !  In -> list of points
    !!----
    !!----    Set up of the integer pointer "pl%p" in the object "pl" of type point_list_type.
    !!----    Each point is associated with the number of an orbit. This pointer is useful
    !!----    to get the asymmetric unit with respect to the input space group of an arbitrary
    !!----    list of points (atom coordinates).
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Set_Orbits_Inlist(Spg,Pl)
       !---- Arguments ----!
       type(space_group_type), intent(in)     :: SpG
       type(point_list_type),  intent(in out) :: pl

       !--- Local variables ---!
       integer                     :: i,j,norb,nt
       real(kind=cp), dimension(3) :: x,xx,v

       norb=0
       pl%p=0
       do i=1,pl%np
          if (pl%p(i) == 0) then
             norb=norb+1
             pl%p(i)=norb
             x=pl%x(:,i)
             do j=1,Spg%multip
                xx=ApplySO(Spg%SymOp(j),x)
                xx=modulo_lat(xx)
                do nt=1,pl%np
                   if (pl%p(nt) /= 0) cycle
                   v=pl%x(:,nt)-xx(:)
                   if (Lattice_trans(v,Spg%spg_lat)) pl%p(nt)=norb
                end do
             end do
          end if
       end do

       return
    End Subroutine Set_Orbits_Inlist

    !!---- Subroutine Set_Rotation_Matrix(ang,Rot)
    !!----   real(kind=cp), dimension(3),   intent( in) :: ang
    !!----   real(kind=cp), dimension(3,3), intent(out) :: Rot
    !!----
    !!----  Subroutine calculating the rotation matrix Rot corresponding to
    !!----  the application (active rotations) of the following succesive rotations:
    !!----
    !!----  Rot = Rx(ang(3)) . Ry(ang(2)) . Rz(ang(1))
    !!----
    !!----    Created: October 2009  (JRC)
    !!----    Updated: March 2013 (JRC)
    !!----

    Subroutine Set_Rotation_Matrix(ang,Rot)
      real(kind=cp), dimension(3),   intent( in) :: ang
      real(kind=cp), dimension(3,3), intent(out) :: Rot
      !Local variables
      real(kind=cp), dimension(3,3) :: Rx,Ry,Rz
      Rx=Matrix_Rx(ang(1),"D")
      Ry=Matrix_Ry(ang(2),"D")
      Rz=Matrix_Rz(ang(3),"D")
      Rot=Matmul(Rx,matmul(Ry,Rz))
      return
    End Subroutine Set_Rotation_Matrix

    !!----
    !!---- Subroutine Set_TDist_Coordination(Max_coor,Dmax, Cell, Spg, A)
    !!----    integer,                  intent(in)   :: max_coor !  Maximum expected coordination
    !!----    real(kind=cp),            intent(in)   :: dmax     !  In -> Max. Distance to calculate
    !!----    real(kind=cp),            intent(in)   :: dangl    !  In -> Max. distance for angle calculations
    !!----    type (Crystal_cell_type), intent(in)   :: Cell     !  In -> Object of Crytal_Cell_Type
    !!----    type (Space_Group_type),  intent(in)   :: SpG      !  In -> Object of Space_Group_Type
    !!----    type (atom_list_type),   intent(in)    :: A        !  In -> Object of atom_list_type
    !!----
    !!----    Subroutine to calculate distances, below the prescribed distance "dmax"
    !!----    Sets up the coordination type: Coord_Info for each atom in the asymmetric unit
    !!----    Needs as input the objects Cell (of type Crystal_cell), SpG (of type Space_Group)
    !!----    and A (of type atom_list, that should be allocated in the calling program).
    !!----    The input argument Max_Coor is obtained, before calling the present procedure,
    !!----    by a call to Allocate_Coordination_Type with arguments:(A%natoms,Spg%Multip,Dmax,max_coor)
    !!----    Further calls to this routine do not need a previous call to Allocate_Coordination_Type.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Set_TDist_Coordination(max_coor,Dmax, Cell, Spg, A)
       !---- Arguments ----!
       integer,                  intent(in)   :: max_coor
       real(kind=cp),            intent(in)   :: dmax
       type (Crystal_cell_Type), intent(in)   :: Cell
       type (Space_Group_Type),  intent(in)   :: SpG
       type (atom_list_type),    intent(in)   :: A

       !---- Local Variables ----!
       integer                              :: i,j,k,lk,i1,i2,i3,nn,L,ico
       integer,       dimension(3)          :: ic1,ic2
       real(kind=cp), dimension(3)          :: xx,x1,xo,Tn,xr, QD
       real(kind=cp)                        :: T,dd
       real(kind=cp), dimension(3,max_coor) :: uu

      ! call init_err_geom()  !Control of error

       qd(:)=1.0/cell%rcell(:)
       !ic2(:)= nint(dmax/cell%cell(:)+1.5)
       ic2(:)= int(dmax/cell%cell(:))+1
       ic1(:)=-ic2(:)
       do i=1,a%natoms
          xo(:)=a%atom(i)%x(:)

          ico=0
          do k=1,a%natoms
             lk=1
             uu(:,lk)=xo(:)
             do j=1,Spg%Multip
                xx=ApplySO(Spg%SymOp(j),a%atom(k)%x)
                do i1=ic1(1),ic2(1)
                   do i2=ic1(2),ic2(2)
                      do_i3:do i3=ic1(3),ic2(3)
                            Tn(1)=real(i1); Tn(2)=real(i2); Tn(3)=real(i3)
                            x1(:)=xx(:)+tn(:)
                            do l=1,3
                               t=abs(x1(l)-xo(l))*qd(l)
                               if (t > dmax) cycle  do_i3
                            end do
                            do nn=1,lk
                               if (sum(abs(uu(:,nn)-x1(:)))  <= epsi) cycle  do_i3
                            end do
                            xr = matmul(cell%cr_orth_cel,x1-xo)
                            dd=sqrt(dot_product(xr,xr))
                            if (dd > dmax .or. dd < 0.001) cycle
                            ico=ico+1
                           ! Control not performed ... it is supposed that max_coor is large enough
                           !if (Coord_Info%Coord_Num(i) > Coord_Info%Max_Coor) then
                           !   err_geom=.true.
                           !   ERR_Geom_Mess=" => Too many distances around an atom"
                           !   return
                           !end if
                            lk=lk+1
                            uu(:,lk)=x1(:)
                            Coord_Info%Dist(ico,i)=dd
                            Coord_Info%N_Cooatm(ico,i)=k
                            Coord_Info%N_sym(ico,i)=j

                            ! Added by JGP
                            Coord_Info%Tr_Coo(:,ico,i)=tn
                      end do do_i3 !i3
                   end do !i2
                end do !i1
             end do !j
          end do !k
          Coord_Info%Coord_Num(i)=ico
       end do !i

       return
    End Subroutine Set_TDist_Coordination

    !!----
    !!---- Subroutine Set_TDist_Partial_Coordination(List,Max_coor,Dmax, Cell, Spg, A)
    !!----    integer,                  intent(in)   :: List     !  Modified atom
    !!----    integer,                  intent(in)   :: max_coor !  Maximum expected coordination
    !!----    real(kind=cp),            intent(in)   :: dmax     !  In -> Max. Distance to calculate
    !!----    real(kind=cp),            intent(in)   :: dangl    !  In -> Max. distance for angle calculations
    !!----    type (Crystal_cell_type), intent(in)   :: Cell     !  In -> Object of Crytal_Cell_Type
    !!----    type (Space_Group_type),  intent(in)   :: SpG      !  In -> Object of Space_Group_Type
    !!----    type (atom_list_type),   intent(in)    :: A        !  In -> Object of atom_list_type
    !!----
    !!----    Modify the coordination type: Coord_Info for the atoms affected by the change of atom "List"
    !!----    Needs as input the objects Cell (of type Crystal_cell), SpG (of type Space_Group)
    !!----    and A (or type atom_list, that should be allocated in the calling program).
    !!----    This routine is a modification of Set_TDist_Coordination to avoid superfluous calculations
    !!----    in global optimization methods. It assumes that Set_TDist_Coordination has previously been
    !!----    called and the object "Coord_Info" has already been set.
    !!----
    !!---- Update: May - 2009
    !!
    Subroutine Set_TDist_Partial_Coordination(List,max_coor,Dmax, Cell, Spg, A)
       !---- Arguments ----!
       integer,                  intent(in)   :: List
       integer,                  intent(in)   :: max_coor
       real(kind=cp),            intent(in)   :: dmax
       type (Crystal_cell_Type), intent(in)   :: Cell
       type (Space_Group_Type),  intent(in)   :: SpG
       type (atom_list_type),    intent(in)   :: A

       !---- Local Variables ----!
       integer                              :: i,j,k,lk,i1,i2,i3,nn,L,ic,ico
       integer,       dimension(3)          :: ic1,ic2
       integer,       dimension(A%natoms)   :: po,pn
       real(kind=cp), dimension(3)          :: xx,x1,xo,Tn,xr, QD
       real(kind=cp)                        :: T,dd
       real(kind=cp), dimension(3,max_coor) :: uu

      ! call init_err_geom()  !Control of error

       po=0; pn=0
       po(List)=1 !This atom has a modified coordination sphere
       ic=Coord_Info%Coord_Num(List)
       do i=1,ic
         po(Coord_Info%N_Cooatm(i,List))=1  !This atom has a modified coordination sphere
       end do

       qd(:)=1.0/cell%rcell(:)
       ic2(:)= int(dmax/cell%cell(:))+1
       ic1(:)=-ic2(:)
       !Determine the new coordination sphere of the changed atom
       i=List
       xo(:)=a%atom(i)%x(:)

       ico=0
       do k=1,a%natoms
          lk=1
          uu(:,lk)=xo(:)
          do j=1,Spg%Multip
             xx=ApplySO(Spg%SymOp(j),a%atom(k)%x)
             do i1=ic1(1),ic2(1)
                do i2=ic1(2),ic2(2)
                   do_i3:do i3=ic1(3),ic2(3)
                         Tn(:)=real((/i1,i2,i3/))
                         x1(:)=xx(:)+tn(:)
                         do l=1,3
                            t=abs(x1(l)-xo(l))*qd(l)
                            if (t > dmax) cycle  do_i3
                         end do
                         do nn=1,lk
                            if (sum(abs(uu(:,nn)-x1(:)))  <= epsi) cycle  do_i3
                         end do
                         xr = matmul(cell%cr_orth_cel,x1-xo)
                         dd=sqrt(dot_product(xr,xr))
                         if (dd > dmax .or. dd < 0.001) cycle
                         ico=ico+1
                         lk=lk+1
                         uu(:,lk)=x1(:)
                         Coord_Info%Dist(ico,i)=dd
                         Coord_Info%N_Cooatm(ico,i)=k
                         Coord_Info%N_sym(ico,i)=j
                   end do do_i3 !i3
                end do !i2
             end do !i1
          end do !j
       end do !k
       Coord_Info%Coord_Num(i)=ico
       pn(list)=0
       po(list)=0

       ic=Coord_Info%Coord_Num(List)    !New coordination number of atom List
       do i=1,ic
         pn(Coord_Info%N_Cooatm(i,List))=1  !This atom has now a newly modified coordination sphere
       end do
       !Look now the changed coordinaion spheres
       do i=1,a%natoms
         if(pn(i) == 0 .and. po(i) == 0) cycle
         !if(po(i) == 1 .and. pn(i) == 1) then !the atom remains in the coordination sphere, only recalculation of distance is needed
         !  ic=Coord_Info%Coord_Num(i)
         !  do k=1,ic
         !   if(List == Coord_Info%N_Cooatm(k,i)) then
         !   end if
         !  end do
         !end if
         !DO ALL WAITING FOR A MORE EFFICIENT ALGORITHM
         xo(:)=a%atom(i)%x(:)

         ico=0
         do k=1,a%natoms
            lk=1
            uu(:,lk)=xo(:)
            do j=1,Spg%Multip
               xx=ApplySO(Spg%SymOp(j),a%atom(k)%x)
               do i1=ic1(1),ic2(1)
                  do i2=ic1(2),ic2(2)
                     do_inter:do i3=ic1(3),ic2(3)
                           Tn(:)=real((/i1,i2,i3/))
                           x1(:)=xx(:)+tn(:)
                           do l=1,3
                              t=abs(x1(l)-xo(l))*qd(l)
                              if (t > dmax) cycle  do_inter
                           end do
                           do nn=1,lk
                              if (sum(abs(uu(:,nn)-x1(:)))  <= epsi) cycle  do_inter
                           end do
                           xr = matmul(cell%cr_orth_cel,x1-xo)
                           dd=sqrt(dot_product(xr,xr))
                           if (dd > dmax .or. dd < 0.001) cycle
                           ico=ico+1
                           lk=lk+1
                           uu(:,lk)=x1(:)
                           Coord_Info%Dist(ico,i)=dd
                           Coord_Info%N_Cooatm(ico,i)=k
                           Coord_Info%N_sym(ico,i)=j
                     end do do_inter !i3
                  end do !i2
               end do !i1
            end do !j
         end do !k
         Coord_Info%Coord_Num(i)=ico
       end do !i

       return
    End Subroutine Set_TDist_Partial_Coordination

 End Module CFML_Geometry_Calc
