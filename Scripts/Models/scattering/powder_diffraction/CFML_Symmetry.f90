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
!!---- MODULE: CFML_Crystallographic_Symmetry
!!----   INFO: This module constains everything needed for handling symmetry
!!----         in Crystallography. Part of the information is obtained from
!!----         tabulated items in the module Symmetry_Tables. In particular
!!----         the correspondence of non standard settings Hermann-Mauguin
!!----         symbols and Hall symbols for space groups.
!!----         The construction of variables of the public type Space_Group_Type
!!----         is done by using a variety of algorithms and methods.
!!----         Many procedures for handling symmetry (symbolic and algebraic)
!!----         are provided in this module.
!!----
!!---- HISTORY
!!----    Update: 05/03/2011
!!----
!!---- DEPENDENCIES
!!----
!!--++    Use CFML_GlobalDeps,       only: Cp
!!--++    Use CFML_Math_General,     only: Trace, Zbelong, Modulo_Lat, equal_matrix, Equal_Vector
!!--++    Use CFML_String_Utilities, only: Equal_Sets_Text, Pack_String, Get_Fraction_2Dig, &
!!--++                                     Get_Fraction_1Dig, Frac_Trans_1Dig, L_Case,     &
!!--++                                     U_case, Ucase, Getnum, Frac_Trans_2Dig
!!--++    Use CFML_Math_3D,          only: Determ_A, matrix_inverse, Resolv_Sist_3x3
!!--++    Use CFML_Symmetry_Tables
!!----
!!----
!!---- VARIABLES
!!----    CUBIC
!!--++    EPS_SYMM                     [Private]
!!----    ERR_SYMM
!!----    ERR_SYMM_MESS
!!--++    GENER_OPER_TYPE              [Private]
!!----    HEXA
!!----    HEXAG
!!----    INLAT
!!----    Lat_Ch
!!----    LATTICE_CENTRING_TYPE
!!----    LTR
!!----    MONOC
!!----    NLAT
!!----    NUM_SPGR_INFO
!!----    ORTHOR
!!----    SPACEG
!!----    SYM_OPER_TYPE
!!----    WYCK_POS_TYPE
!!----    WYCKOFF_TYPE
!!----    SPACE_GROUP_TYPE
!!----    TETRA
!!----    TRIGO
!!----
!!---- PROCEDURES
!!----    Functions:
!!----       APPLYSO
!!----       AXES_ROTATION
!!--++       EQUAL_SYMOP               [Overloaded Operator]
!!--++       EQUIV_SYMOP               [Private]
!!----       GET_LAUE_NUM
!!----       GET_MULTIP_POS
!!----       GET_OCC_SITE
!!----       GET_POINTGROUP_NUM
!!--++       IS_AXIS                   [Private]
!!--++       IS_DIGIT                  [Private]
!!--++       IS_HEXA                   [Private]
!!----       IS_NEW_OP
!!--++       IS_PLANE                  [Private]
!!--++       IS_XYZ                    [Private]
!!----       LATTICE_TRANS
!!--++       PRODUCT_SYMOP             [Overloaded Operator]
!!----       SPGR_EQUAL
!!----       SYM_PROD
!!----
!!----    Subroutines:
!!----       ALLOCATE_LATTICE_CENTRING
!!--++       CHECK_SYMBOL_HM           [Private]
!!----       CHECK_GENERATOR
!!----       COPY_NS_SPG_TO_SPG
!!----       DECODMATMAG
!!----       GET_CENTRING_VECTORS
!!----       GET_CRYSTAL_SYSTEM
!!--++       GET_CRYSTAL_SYSTEM_R_OP   [Overloaded]
!!--++       GET_CRYSTAL_SYSTEM_R_ST   [Overloaded]
!!----       GET_GENSYMB_FROM_GENER
!!----       GET_HALLSYMB_FROM_GENER
!!----       GET_LATTICE_TYPE
!!----       GET_LAUE_PG
!!----       GET_LAUE_STR
!!----       GET_ORBIT
!!----       GET_POINTGROUP_STR
!!--++       GET_SEITZ                 [Private]
!!--++       GET_SEITZ_SYMBOL
!!--++       GET_SETTING_INFO          [Private]
!!----       GET_SHUBNIKOV_OPERATOR_SYMBOL
!!----       GET_SO_FROM_FIX
!!----       GET_SO_FROM_GENER
!!----       GET_SO_FROM_HALL
!!----       GET_SO_FROM_HMS
!!----       GET_STABILIZER
!!----       GET_STRING_RESOLV
!!----       GET_SUBORBITS
!!----       GET_SYMEL
!!----       GET_SYMKOV
!!----       GET_SYMSYMB
!!--++       GET_SYMSYMBI              [Overloaded]
!!--++       GET_SYMSYMBR              [Overloaded]
!!----       GET_T_SUBGROUPS
!!----       GET_TRASFM_SYMBOL
!!----       GET_TRANSL_SYMBOL
!!----       INIT_ERR_SYMM
!!----       INVERSE_SYMM
!!----       LATSYM
!!--++       MAX_CONV_LATTICE_TYPE     [Private]
!!--++       MOD_TRANS                 [Private]
!!----       READ_BIN_SPACEGROUP
!!----       READ_MSYMM
!!----       READ_SYMTRANS_CODE
!!----       READ_XSYM
!!----       SEARCHOP
!!----       SET_SPACEGROUP
!!----       SET_SPG_MULT_TABLE
!!----       SETTING_CHANGE
!!--++       SETTING_CHANGE_CONV       [Overloaded]
!!--++       SETTING_CHANGE_NONCONV    [Overloaded]
!!----       SIMILAR_TRANSF_SG
!!----       SYM_B_RELATIONS
!!--++       SYM_B_RELATIONS_OP        [Overloaded]
!!--++       SYM_B_RELATIONS_ST        [Overloaded]
!!----       SYM_PROD_ST
!!----       SYMMETRY_SYMBOL
!!--++       SYMMETRY_SYMBOL_OP        [Overloaded]
!!--++       SYMMETRY_SYMBOL_STR       [Overloaded]
!!--++       SYMMETRY_SYMBOL_XYZ       [Overloaded]
!!----       WRITE_BIN_SPACEGROUP
!!----       WRITE_SPACEGROUP
!!----       WRITE_SYM
!!----       WRITE_SYMTRANS_CODE
!!----       WRITE_WYCKOFF
!!----       WYCKOFF_ORBIT
!!----
!!--..    Operators:
!!--..       (*)
!!--..       (==)
!!----
!!
 Module CFML_Crystallographic_Symmetry

    !---- Used External Modules ----!
    Use CFML_GlobalDeps,       only: cp
    Use CFML_Math_General,     only: Trace, Zbelong, Modulo_Lat, equal_matrix,             &
                                     Equal_Vector,Sort,Set_Epsg,Set_Epsg_Default
    Use CFML_Math_3D,          only: Determ_A, matrix_inverse, Resolv_Sist_3x3
    Use CFML_String_Utilities, only: Equal_Sets_Text, Pack_String, Get_Fraction_2Dig,      &
                                     Get_Fraction_1Dig, Frac_Trans_1Dig, L_Case,           &
                                     U_case, Ucase, Getnum, Frac_trans_2Dig, Get_Num_String
    Use CFML_Symmetry_Tables

    implicit none

    private

    !---- List of public variables and types ----!

    !---- List of public overloaded operators ----!
    public ::  operator (*), operator (==)

    !---- List of public functions ----!
    public  :: ApplySO, Axes_Rotation, Get_Laue_Num, Get_Multip_Pos, Get_Occ_Site,     &
               Get_Pointgroup_Num, Is_New_Op, Lattice_Trans, Spgr_Equal, Sym_Prod

    !---- List of public subroutines ----!
    public  :: Decodmatmag, Get_Centring_Vectors, Get_Crystal_System, Get_Lattice_Type,              &
               Get_Laue_Pg, Get_Laue_Str, Get_orbit, Get_Pointgroup_Str, Get_So_From_Fix,            &
               Get_So_From_Gener,Get_So_From_Hall, Get_So_From_Hms, Get_HallSymb_From_Gener,         &
               Get_Stabilizer,Get_String_Resolv,Get_SubOrbits,Get_Symel, Get_Symkov, Get_SymSymb,    &
               Init_Err_Symm, Inverse_Symm, Latsym, Read_Msymm, Read_Xsym, Searchop,                 &
               Set_Spacegroup, Setting_Change, Sym_B_Relations, Sym_Prod_St, Symmetry_Symbol,        &
               Write_Spacegroup, Write_Sym, Write_Wyckoff, Wyckoff_Orbit, Get_T_SubGroups,           &
               Similar_Transf_SG, Read_SymTrans_Code, Write_SymTrans_Code, Set_SpG_Mult_Table,       &
               Get_Seitz_Symbol, Get_Trasfm_Symbol,Get_Shubnikov_Operator_Symbol,                    &
               Get_Transl_Symbol, Read_Bin_Spacegroup, Write_Bin_Spacegroup, Get_GenSymb_from_Gener, &
               Check_Generator, Copy_NS_SpG_To_SpG, Allocate_Lattice_Centring

    !---- List of private Operators ----!
    private :: Equal_Symop, Product_Symop

    !---- List of private functions ----!
    private :: Is_Axis, Is_Digit, Is_Hexa, Is_Plane, Is_Xyz, Equiv_Symop

    !---- List of private subroutines ----!
    private :: Check_Symbol_Hm, Get_Seitz, Get_SymSymbI, Get_SymSymbR, Mod_Trans, Sym_B_Relations_Op  , &
               Sym_B_Relations_St, Symmetry_Symbol_Op, Symmetry_Symbol_Xyz , Symmetry_Symbol_Str,       &
               Max_Conv_Lattice_Type,Get_Setting_Info,Get_Crystal_System_R_OP,Get_Crystal_System_R_ST, &
               Setting_Change_Conv,Setting_Change_NonConv


    !---- Global Variables ----!

    !---- Definitions ----!

    !!----
    !!---- CUBIC
    !!----    integer, parameter, public :: Cubic
    !!----
    !!----    Cubic parameter index: Cubic = 554
    !!----
    !!---- Update: February - 2005
    !!
    integer, parameter, public :: Cubic = 554

    !!--++
    !!--++ eps_symm
    !!--++    real(kind=cp), parameter, private :: eps_symm
    !!--++
    !!--++    (PRIVATE)
    !!--++    Epsilon for comparisons within this module
    !!--++
    !!--++ Update: February - 2005
    !!
    real(kind=cp), parameter, private :: eps_symm  = 0.0002_cp

    !!----
    !!---- ERR_SYMM_MESS
    !!----    character(len=150), public :: ERR_Symm_Mess
    !!----
    !!----    String containing information about the last error
    !!----
    !!---- Update: February - 2005
    !!
    character(len=150), public :: ERR_Symm_Mess

    !!----
    !!---- ERR_SYMM
    !!----    logical, public :: Err_Symm
    !!----
    !!----    Logical Variable to indicate an error on this module.
    !!----
    !!---- Update: February - 2005
    !!
    logical, public :: Err_Symm

    !!--++
    !!--++ TYPE :: GENER_OPER_TYPE
    !!--..
    !!--++ Type, private :: Gener_Oper_Type
    !!--++    integer          :: orden
    !!--++    character(len=1) :: axes
    !!--++    character(len=1) :: axes2
    !!--++    character(len=2) :: tras
    !!--++ End Type Gener_Oper_Type
    !!--++
    !!--++ Update: February - 2005
    !!
    Type, private :: Gener_Oper_Type
       integer          :: orden
       character(len=1) :: axes
       character(len=1) :: axes2
       character(len=3) :: tras
    End Type Gener_Oper_Type

    !!----
    !!---- HEXA
    !!----    logical, public :: Hexa
    !!----
    !!----    .false. Rotational part of symmetry operators  belongs to m3m
    !!----    .true.  Rotational part of symmetry operators  belongs to 6/mmm
    !!----
    !!---- Update: February - 2005
    !!
    logical, public :: Hexa

    !!----
    !!---- HEXAG
    !!----    integer, parameter, public :: Hexag
    !!----
    !!----    Index parameter for hexagonal Groups: Hexag  = 527
    !!----
    !!---- Update: February - 2005
    !!
    integer, parameter, public :: Hexag         = 527

    !!----
    !!---- INLAT
    !!----    integer, public        :: Inlat
    !!----
    !!----    Ordinal index of the lattice
    !!----
    !!---- Update: February - 2005
    !!
    integer, public        :: Inlat

    !!----
    !!---- Lat_Ch
    !!----    character(len= 1), public     :: Lat_Ch
    !!----
    !!----    First character of the space group symbol
    !!----
    !!---- Update: February - 2005
    !!
    character(len= 1), public     :: Lat_Ch

    !!----
    !!---- TYPE :: Lattice_Centring_Type
    !!--..
    !!---- Type, public :: Lattice_Centring_Type
    !!----    integer                                     :: N_lat
    !!----    logical                                     :: set
    !!----    real(kind=cp), dimension(:,:),allocatable   :: LTr
    !!---- End Type Lattice_Centring_Type
    !!----
    !!----   Lattice centring translations (including anti-translation)
    !!----   symmetry operators defined with respect to arbitrary axes.
    !!----   Normally the first translation is the identity element of the translation
    !!----   group: Ltr(:,1)=[0,0,0] or [0,0,0,1] if time inversion is considered to take
    !!----   into account also the anti-translations.
    !!----   For using this type first the program should allocate the arrays by calling
    !!----   the subroutine Allocate_Lattice_Centring and then construct totally the object
    !!----   by assigning appropriate values and putting set=.true.
    !!----
    !!---- Update: October - 2014
    !!
    Type, public :: Lattice_Centring_Type
       integer                                     :: N_lat
       logical                                     :: set
       real(kind=cp), dimension(:,:),allocatable   :: LTr
    End Type Lattice_Centring_Type


    !!----
    !!---- LTR
    !!----    real(kind=cp), dimension(3,192), public  :: Ltr
    !!----
    !!----    Centering Lattice Translations, up to 192 lattice centring
    !!----    vectors are allowed. Conventional lattice centring need only 4 vectors
    !!----
    !!---- Update: February - 2005, January-2014
    !!
    real(kind=cp), dimension(3,192), public  :: Ltr    ! Centering Lattice Translations

    !!----
    !!---- MONOC
    !!----    integer, parameter, public :: Monoc
    !!----
    !!----    Index parameter for Monoclinic Groups: Monoc  =  15
    !!----
    !!---- Update: February - 2005
    !!
    integer, parameter, public :: Monoc         =  15

    !!----
    !!---- NLAT
    !!----    integer, public      :: Nlat
    !!----
    !!----    Multiplicity of the lattice
    !!----
    !!---- Update: February - 2005
    !!
    integer, public      :: Nlat

    !!----
    !!---- TYPE :: NS_SYM_OPER_TYPE
    !!--..
    !!---- Type, public :: NS_Sym_Oper_Type
    !!----    real(kind=cp), dimension(3,3) :: Rot     !  Rotational Part of Symmetry Operator
    !!----    real(kind=cp), dimension(3)   :: Tr      !  Traslational part of Symmetry Operator
    !!---- End Type  NS_Sym_Oper_Type
    !!----
    !!----   Non-standard symmetry operator. Needed for describing non-standard space groups with
    !!----   symmetry operators defined with respect to arbitrary axes
    !!----
    !!---- Update: January - 2014
    !!
    Type, public :: NS_Sym_Oper_Type
       real(kind=cp), dimension(3,3) :: Rot
       real(kind=cp), dimension(3)   :: Tr
    End Type NS_Sym_Oper_Type

    !!----
    !!---- TYPE :: NS_SPACE_GROUP_TYPE
    !!--..
    !!---- Type, public :: NS_Space_Group_Type
    !!----    Integer                                          :: NumSpg        ! Number of the Space Group
    !!----    Character(len=20)                                :: SPG_Symb      ! Hermann-Mauguin Symbol
    !!----    Character(len=16)                                :: Hall          ! Hall symbol
    !!----    Character(len=90)                                :: gHall         ! Generalised Hall symbol
    !!----    Character(len=12)                                :: CrystalSys    ! Crystal system
    !!----    Character(len= 5)                                :: Laue          ! Laue Class
    !!----    Character(len= 5)                                :: PG            ! Point group
    !!----    Character(len= 5)                                :: Info          ! Extra information
    !!----    Character(len=90)                                :: SG_setting    ! Information about the SG setting
    !!----                                                                      ! (IT,KO,ML,ZA,Table,Standard,UnConventional)
    !!----    Character(len= 1)                                :: SPG_lat       ! Lattice type
    !!----    Character(len= 2)                                :: SPG_latsy     ! Lattice type Symbol
    !!----    Integer                                          :: NumLat        ! Number of lattice points in a cell
    !!----    real(kind=cp), allocatable, dimension(:,:)       :: Latt_trans    ! Lattice translations
    !!----    Character(len=51)                                :: Bravais       ! String with Bravais symbol + translations
    !!----    Character(len=80)                                :: Centre        ! Information about Centric or Acentric
    !!----    Integer                                          :: Centred       ! =0 Centric(-1 no at origin)
    !!----                                                                      ! =1 Acentric
    !!----                                                                      ! =2 Centric(-1 at origin)
    !!----    real(kind=cp), dimension(3)                      :: Centre_coord  ! Fractional coordinates of the inversion centre
    !!----    Integer                                          :: NumOps        ! Number of reduced set of S.O.
    !!----    Integer                                          :: Multip        ! Multiplicity of the general position
    !!----    Integer                                          :: Num_gen       ! Minimum number of operators to generate the Group
    !!----    type(NS_Sym_Oper_Type), allocatable, dimension(:):: SymOp         ! Symmetry operators (192)
    !!----    Character(len=50),      allocatable, dimension(:):: SymopSymb     ! Strings form of symmetry operators (192)
    !!---- End Type NS_Space_Group_Type
    !!----
    !!----  Definition of the type NS_Space_Group_Type: Non-standard space group. This type has been created
    !!----  in order to be able to describe symmetry operators with non-integer values when they are referred
    !!----  to arbitrary settings. They are created only as intermediate variables in some calculations.
    !!----
    !!---- Updated: February - 2005, January 2014 (JRC to make some components allocatable and change the length of some strings)
    !!
    Type, public :: NS_Space_Group_Type
       integer                                          :: NumSpg=0         ! Number of the Space Group
       character(len=20)                                :: SPG_Symb=" "     ! Hermann-Mauguin Symbol
       character(len=16)                                :: Hall=" "         ! Hall symbol
       character(len=90)                                :: gHall=" "        ! Generalised Hall symbol
       character(len=12)                                :: CrystalSys=" "   ! Crystal system
       character(len= 5)                                :: Laue=" "         ! Laue Class
       character(len= 5)                                :: PG=" "           ! Point group
       character(len= 5)                                :: Info=" "         ! Extra information
       character(len=90)                                :: SG_setting=" "   ! Information about the SG setting (IT,KO,ML,ZA,Table,Standard,UnConventional)
       character(len= 1)                                :: SPG_lat=" "      ! Lattice type
       character(len= 2)                                :: SPG_latsy=" "    ! Lattice type Symbol
       integer                                          :: NumLat=1         ! Number of lattice points in a cell
       real(kind=cp), allocatable,dimension(:,:)        :: Latt_trans       ! Lattice translations (3,12)
       character(len=51)                                :: Bravais=" "      ! String with Bravais symbol + translations
       character(len=80)                                :: Centre=" "       ! Alphanumeric information about the center of symmetry
       integer                                          :: Centred=0        ! Centric or Acentric [ =0 Centric(-1 no at origin),=1 Acentric,=2 Centric(-1 at origin)]
       real(kind=cp), dimension(3)                      :: Centre_coord=0.0 ! Fractional coordinates of the inversion centre
       integer                                          :: NumOps=0         ! Number of reduced set of S.O.
       integer                                          :: Multip=0         ! Multiplicity of the general position
       integer                                          :: Num_gen=0        ! Minimum number of operators to generate the Group
       type(NS_Sym_Oper_Type), allocatable,dimension(:) :: SymOp            ! Symmetry operators (192)
       character(len=50),      allocatable,dimension(:) :: SymopSymb        ! Strings form of symmetry operators
    End Type NS_Space_Group_Type

    !!----
    !!---- NUM_SPGR_INFO
    !!----    integer, parameter, public :: Num_Spgr_Info
    !!----
    !!----    Total dimension of SPGR_INFO: Num_Spgr_Info = 612
    !!----
    !!---- Update: February - 2005
    !!
    integer, parameter, public :: Num_Spgr_Info = 612

    !!----
    !!---- ORTHOR
    !!----    integer, parameter, public :: Orthor
    !!----
    !!----    Index parameter for Orthorhombic Groups: Orthor  = 163
    !!----
    !!---- Update: February - 2005
    !!
    integer, parameter, public :: Orthor  = 163

    !!----
    !!---- SPACEG
    !!----     character(len=20), public     :: SpaceG
    !!----
    !!----     Space group symbol
    !!----
    !!---- Update: February - 2005
    !!
    character(len=20), public   :: SpaceG

    !!----
    !!---- TYPE :: SYM_OPER_TYPE
    !!--..
    !!---- Type, public :: Sym_Oper_Type
    !!----    integer,       dimension(3,3) :: Rot     !  Rotational Part of Symmetry Operator
    !!----    real(kind=cp), dimension(3)   :: Tr      !  Traslational part of Symmetry Operator
    !!---- End Type  Sym_Oper_Type
    !!----
    !!----    Definition of Variable
    !!----
    !!---- Update: February - 2005
    !!
    Type, public :: Sym_Oper_Type
       integer,       dimension(3,3) :: Rot
       real(kind=cp), dimension(3)   :: Tr
    End Type Sym_Oper_Type

    !!----
    !!---- TYPE :: WYCK_POS_TYPE
    !!--..
    !!---- Type, public :: Wyck_Pos_Type
    !!----    integer                         :: multp     ! Multiplicity
    !!----    character(len= 6)               :: site      ! Site Symmetry
    !!----    integer                         :: norb      ! Number of elements in orbit
    !!----    character(len=40)               :: orig      ! Orig
    !!----    character(len=40),dimension(48) :: str_orbit ! Orbit
    !!----    character(len=40),dimension(192):: extra_orbit
    !!---- End Type Wyck_Pos_Type
    !!----
    !!----    Definition of Variable
    !!----
    !!---- Update: February - 2005
    !!
    Type, public :: Wyck_Pos_Type
       integer                          :: multp
       character(len= 6)                :: site
       integer                          :: norb
       character(len=40)                :: str_orig
       character(len=40),dimension(48)  :: str_orbit
    End Type Wyck_Pos_Type

    !!----
    !!---- TYPE :: WYCKOFF_TYPE
    !!--..
    !!---- Type, public :: Wyckoff_Type
    !!----    integer                            :: num_orbit      ! Number of orbits
    !!----    type(wyck_pos_type), dimension(26) :: orbit          ! Orbit type
    !!---- End Type Wyckoff_Type
    !!----
    !!----    Definition of Variable
    !!----
    !!---- Update: February - 2005
    !!
    Type, public :: Wyckoff_Type
       integer                            :: num_orbit
       type(wyck_pos_type), dimension(26) :: orbit
    End Type Wyckoff_Type

    !!----
    !!---- TYPE :: SPACE_GROUP_TYPE
    !!--..
    !!---- Type, public :: Space_Group_Type
    !!----    Integer                                         :: NumSpg        ! Number of the Space Group
    !!----    Character(len=20)                               :: SPG_Symb      ! Hermann-Mauguin Symbol
    !!----    Character(len=16)                               :: Hall          ! Hall symbol
    !!----    Character(len=90)                               :: gHall         ! Generalised Hall symbol
    !!----    Character(len=12)                               :: CrystalSys    ! Crystal system
    !!----    Character(len= 5)                               :: Laue          ! Laue Class
    !!----    Character(len= 5)                               :: PG            ! Point group
    !!----    Character(len= 5)                               :: Info          ! Extra information
    !!----    Character(len=90)                               :: SG_setting    ! Information about the SG setting
    !!----                                                                     ! (IT,KO,ML,ZA,Table,Standard,UnConventional)
    !!----    Logical                                         :: Hexa          !
    !!----    Character(len= 1)                               :: SPG_lat       ! Lattice type
    !!----    Character(len= 2)                               :: SPG_latsy     ! Lattice type Symbol
    !!----    Integer                                         :: NumLat        ! Number of lattice points in a cell
    !!----    real(kind=cp), allocatable, dimension(:,:)      :: Latt_trans    ! Lattice translations
    !!----    Character(len=51)                               :: Bravais       ! String with Bravais symbol + translations
    !!----    Character(len=80)                               :: Centre        ! Information about Centric or Acentric
    !!----    Integer                                         :: Centred       ! =0 Centric(-1 no at origin)
    !!----                                                                     ! =1 Acentric
    !!----                                                                     ! =2 Centric(-1 at origin)
    !!----    real(kind=cp), dimension(3)                     :: Centre_coord  ! Fractional coordinates of the inversion centre
    !!----    Integer                                         :: NumOps        ! Number of reduced set of S.O.
    !!----    Integer                                         :: Multip        ! Multiplicity of the general position
    !!----    Integer                                         :: Num_gen       ! Minimum number of operators to generate the Group
    !!----    type(Sym_Oper_Type), allocatable, dimension(:)  :: SymOp         ! Symmetry operators (192)
    !!----    Character(len=50),   allocatable, dimension(:)  :: SymopSymb     ! Strings form of symmetry operators (192)
    !!----    type(wyckoff_type)                              :: Wyckoff       ! Wyckoff Information
    !!----    real(kind=cp), dimension(3,2)                   :: R_Asym_Unit   ! Asymmetric unit in real space
    !!---- End Type Space_Group_Type
    !!----
    !!----     Definition of a variable type Space_Group_Type
    !!----
    !!---- Updated: February - 2005, January 2014 (JRC to make some components allocatable and change the length of some strings)
    !!
    Type, public :: Space_Group_Type
       integer                                       :: NumSpg=0         ! Number of the Space Group
       character(len=20)                             :: SPG_Symb=" "     ! Hermann-Mauguin Symbol
       character(len=16)                             :: Hall=" "         ! Hall symbol
       character(len=90)                             :: gHall=" "        ! Generalised Hall symbol
       character(len=12)                             :: CrystalSys=" "   ! Crystal system
       character(len= 5)                             :: Laue=" "         ! Laue Class
       character(len= 5)                             :: PG=" "           ! Point group
       character(len= 5)                             :: Info=" "         ! Extra information
       character(len=90)                             :: SG_setting=" "   ! Information about the SG setting (IT,KO,ML,ZA,Table,Standard,UnConventional)
       logical                                       :: Hexa=.false.     !
       character(len= 1)                             :: SPG_lat=" "      ! Lattice type
       character(len= 2)                             :: SPG_latsy=" "    ! Lattice type Symbol
       integer                                       :: NumLat=0         ! Number of lattice points in a cell
       real(kind=cp), allocatable,dimension(:,:)     :: Latt_trans       ! Lattice translations (3,12)
       character(len=51)                             :: Bravais=" "      ! String with Bravais symbol + translations
       character(len=80)                             :: Centre=" "       ! Alphanumeric information about the center of symmetry
       integer                                       :: Centred=0        ! Centric or Acentric [ =0 Centric(-1 no at origin),=1 Acentric,=2 Centric(-1 at origin)]
       real(kind=cp), dimension(3)                   :: Centre_coord=0.0 ! Fractional coordinates of the inversion centre
       integer                                       :: NumOps=0         ! Number of reduced set of S.O.
       integer                                       :: Multip=0         ! Multiplicity of the general position
       integer                                       :: Num_gen          ! Minimum number of operators to generate the Group
       type(Sym_Oper_Type), allocatable,dimension(:) :: SymOp            ! Symmetry operators (192)
       character(len=50),   allocatable,dimension(:) :: SymopSymb        ! Strings form of symmetry operators
       type(Wyckoff_Type)                            :: Wyckoff          ! Wyckoff Information
       real(kind=cp),dimension(3,2)                  :: R_Asym_Unit=0.0  ! Asymmetric unit in real space
    End Type Space_Group_Type

    !!----
    !!---- TETRA
    !!----    integer, parameter, public :: Tetra
    !!----
    !!----    Index parameter for Tetragonal Groups: Tetra = 410
    !!----
    !!---- Update: February - 2005
    !!
    integer, parameter, public :: Tetra = 410

    !!----
    !!---- TRIGO
    !!----    integer, parameter, public :: Trigo
    !!----
    !!----    Index parameter for Trigonal Groups: Trigo = 495
    !!----
    !!---- Update: February - 2005
    !!
    integer, parameter, public :: Trigo = 495


    !---- Interfaces Definitions for Overload ----!

    Interface  Get_Crystal_System
       Module Procedure Get_Crystal_System_R_OP
       Module Procedure Get_Crystal_System_R_ST
    End Interface  Get_Crystal_System

    Interface  Get_SymSymb
       Module Procedure Get_SymSymbI
       Module Procedure Get_SymSymbR
    End Interface  Get_SymSymb

    Interface  Setting_Change
       Module Procedure Setting_Change_Conv
       Module Procedure Setting_Change_NonConv
    End Interface  Setting_Change

    Interface  Sym_B_Relations
       Module Procedure Sym_B_Relations_Op
       Module Procedure Sym_B_Relations_St
    End Interface  Sym_B_Relations

    Interface  Symmetry_Symbol
       Module Procedure Symmetry_Symbol_Op
       Module Procedure Symmetry_Symbol_Str
       Module Procedure Symmetry_Symbol_Xyz
    End Interface  Symmetry_Symbol

    Interface operator (*)
       Module Procedure Product_Symop
    End Interface

    Interface Operator (==)
       Module Procedure Equal_Symop
    End Interface

 Contains

    !---- Functions ----!

    !!----
    !!---- Function Applyso(Op,V) Result(Applysop)
    !!----    Type(Sym_Oper_Type),          intent(in) :: Op        !  In -> Symmetry Operator Type
    !!----    real(kind=cp), dimension(3) , intent(in) :: v         !  In -> Vector
    !!----    real(kind=cp), dimension(3)              :: ApplySOp  ! Out -> Output vector
    !!----
    !!----    Apply a symmetry operator to a vector:  Vp = ApplySO(Op,v)
    !!----
    !!---- Update: February - 2005
    !!
    Function ApplySO(Op,V) Result(Applysop)
       !---- Arguments ----!
       Type(Sym_Oper_Type),          intent(in) :: Op
       real(kind=cp), dimension(3),  intent(in) :: v
       real(kind=cp), dimension(3)              :: ApplySOp

       ApplySOp = matmul(Op%Rot,v) + Op%tr

       return
    End Function ApplySO

    !!----
    !!---- Function Axes_Rotation(R) Result(N)
    !!----    integer, dimension(3,3), intent  (in) :: r    !  In -> Rotation part of Symmetry Operator
    !!----    integer                               :: n    ! Out -> Orden of the Rotation Part
    !!----
    !!----    Determine the orden of rotation (valid for all bases). Return a zero
    !!----    if any error occurs.
    !!----
    !!---- Update: February - 2005
    !!
    Function Axes_Rotation(R) Result(N)
       !---- Arguments ----!
       integer, dimension(3,3), intent (in) :: r
       integer                              :: n

       !---- Local Variables ----!
       integer :: det,itr

       n=0

       det=determ_A(r)
       itr=trace(r)
       select case (itr)
          case (-3)
             if (det == -1) n=-1

          case (-2)
             if (det == -1) n=-6

          case (-1)
             if (det == -1) n=-4
             if (det ==  1) n= 2

          case ( 0)
             if (det == -1) n=-3
             if (det ==  1) n= 3

          case ( 1)
             if (det == -1) n=-2
             if (det ==  1) n= 4

          case ( 2)
             if (det ==  1) n= 6

          case ( 3)
             if (det ==  1) n= 1
       end select

       return
    End Function Axes_Rotation

    !!--++
    !!--++ Function Equal_Symop(Syma,Symb) Result (Aeqb)
    !!--++    type(Sym_Oper_Type), intent (in) :: syma
    !!--++    type(Sym_Oper_Type), intent (in) :: symb
    !!--++    logical                          :: aeqb
    !!--++
    !!--++    (OVERLOADED)
    !!--++    The result is .true. if syma == symb, otherwise is .false.
    !!--++    It overloads the "==" operator for objects of type Sym_Oper_Type.
    !!--++    The calling program can use a statement like: if(syma == symb) then ...
    !!--++
    !!--++  Update: February - 2005
    !!
    Function Equal_Symop(Syma,Symb) Result (Aeqb)
       !---- Arguments ----!
       type(Sym_Oper_Type), intent (in) :: syma
       type(Sym_Oper_Type), intent (in) :: symb
       logical                          :: aeqb

       !---- Local variables ----!
       integer :: i,j

       aeqb=.false.
       do i=1,3
          if (abs(Syma%tr(i)-Symb%tr(i)) > eps_symm) return
       end do

       do i=1,3
          do j=1,3
             if (abs(Syma%Rot(i,j)-Symb%Rot(i,j)) > eps_symm) return
          end do
       end do
       aeqb=.true.

       return
    End Function Equal_Symop

    !!--++
    !!--++ Equiv_Symop(Syma,Symb,Lat) Result (Aeqb)
    !!--++    type(Sym_Oper_Type), intent (in) :: syma
    !!--++    type(Sym_Oper_Type), intent (in) :: symb
    !!--++    character (len=*),   intent (in) :: lat
    !!--++    logical                          :: aeqb
    !!--++
    !!--++    The result is .true. if Syma  differ from Symb just by a lattice
    !!--++    translation. This Function is used by the subroutine constructing
    !!--++    the multiplication table of the factor group of a space group.
    !!--++
    !!--++  Update: April - 2005
    !!
    Function Equiv_Symop(Syma,Symb,Lat) Result (Aeqb)
       !---- Arguments ----!
       type(Sym_Oper_Type), intent (in) :: syma
       type(Sym_Oper_Type), intent (in) :: symb
       character (len=*),   intent (in) :: lat
       logical                          :: aeqb

       !---- Local variables ----!
       integer                     :: i,j
       real(kind=cp), dimension(3) :: tr

       aeqb=.false.
       tr= Syma%tr-Symb%tr
       if (.not. Lattice_Trans(tr,Lat)) return
       do i=1,3
          do j=1,3
             if (abs(Syma%Rot(i,j)-Symb%Rot(i,j)) > 0) return
          end do
       end do
       aeqb=.true.

       return
    End Function Equiv_Symop


    !!----
    !!---- Function Get_Laue_Num(Laueclass) Result(Lnum)
    !!----    character(len=*), intent (in) :: laueclass    !  In -> Laue Class string
    !!----    integer                       :: lnum         ! Out -> Ordinal number according LAUE_CLASS
    !!----
    !!----    Obtain the ordinal number corresponding to the Laue-Class
    !!----    symbol according to Laue_Class array. Zero if error is present
    !!----
    !!---- Update: February - 2005
    !!
    Function Get_Laue_Num(Laueclass) Result(Lnum)
       !---- Arguments ----!
       character(len=*), intent (in) :: laueclass
       integer                       :: lnum

       !---- Local Variables ----!
       integer                       :: i
       character(len=len(laueclass)) :: laue

       lnum=0
       laue=adjustl(laueclass)

       do i=1,16
          if (laue(1:5) == laue_class(i)) then
             lnum=i
             exit
          end if
       end do
       if (lnum==15) lnum=13
       if (lnum==16) lnum=14

       return
    End Function Get_Laue_Num

    !!----
    !!----  Function Get_Multip_Pos(X,Spg) Result(Mult)
    !!----    real(kind=cp), dimension(3), intent (in) :: x        !  In -> Position vector
    !!----    type(Space_Group_type),      intent (in) :: spgr     !  In -> Space Group
    !!----    integer                                  :: mult     !  Result -> Multiplicity
    !!----
    !!----    Obtain the multiplicity of a real space point given the space group.
    !!----
    !!---- Update: February - 2005
    !!
    Function Get_Multip_Pos(x,Spg) Result(mult)
       !---- Arguments ----!
       real(kind=cp), dimension(3),  intent (in) :: x
       type(Space_Group_type),       intent (in) :: spg
       integer                                   :: mult

       !---- Local variables ----!
       integer                                :: j, nt
       real(kind=cp), dimension(3)            :: xx,v
       real(kind=cp), dimension(3,Spg%multip) :: u

       !> Init Epss
       call set_epsg(1.0e-3)

       mult=1
       u(:,1)=x(:)

       ext: do j=2,Spg%multip
          xx=ApplySO(Spg%SymOp(j),x)
          xx=modulo_lat(xx)
          do nt=1,mult
             v=u(:,nt)-xx(:)
             if (Lattice_trans(v,Spg%spg_lat)) cycle ext
          end do
          mult=mult+1
          u(:,mult)=xx(:)
       end do ext

       mult=mult*Spg%Numlat

       !> Reset value for epss
       call set_epsg_default()

       return
    End Function Get_Multip_Pos

    !!----
    !!---- Function Get_Occ_Site(Pto,Spg) Result(Occ)
    !!----    real(kind=cp),dimension(3),intent (in) :: Pto ! Point for Occupancy calculation
    !!----    Type (Space_Group_Type),   intent(in)  :: Spg ! Space Group
    !!----    real(kind=cp)                          :: Occ ! Result
    !!----
    !!----    Obtain the occupancy factor (site multiplicity/multiplicity) for Pto
    !!----
    !!---- Update: February - 2005
    !!
    Function Get_Occ_Site(Pto,Spg) Result(Occ)
       !---- Arguments ----!
       real(kind=cp), dimension(3),intent(in) :: Pto
       type (Space_Group_Type),    intent(in) :: Spg
       real(kind=cp)                          :: Occ

       !---- Local Variables ----!

       !> Init Epss
       call set_epsg(1.0e-3)

       Occ=real(Get_Multip_pos(pto,Spg))/real(Spg%multip)

       !> Reset value Epss
       call set_epsg_default()

       return
    End Function Get_Occ_Site

    !!----
    !!---- Function Get_Pointgroup_Num(Pgname) Result(Ipg)
    !!----    character(len=*), intent (in) :: pgname        !  In -> String for PointGroup
    !!----    integer                       :: ipg           ! Out -> Ordinal number as POINT_GROUP
    !!----
    !!----    Obtain the ordinal number corresponding to the Point Group
    !!----    symbol according to Point_Group array. Zero if Error is present
    !!----
    !!---- Update: July - 2014: added m3 and m3m for compatibility with Laue_class
    !!
    Function Get_Pointgroup_Num(Pgname) Result(Ipg)
       !---- Arguments ----!
       character(len=*), intent (in) :: pgname
       integer                       :: ipg

       !---- Local Variables ----!
       integer                       :: i
       character(len=len(pgname))    :: pg

       ipg=0
       pg=adjustl(pgname)

       do i=1,41                ! was 39, now 41 to accomodate m3 and m3m
          if (pg(1:5) == point_group(i)) then
             ipg=i
             exit
          end if
       end do

       !> return previous numbers for m3 and m3m
       if(ipg == 40)ipg=36      ! m3 now m-3
       if(ipg == 41)ipg=39      ! m3m now m-3m

       return
    End Function Get_PointGroup_Num

    !!--++
    !!--++ Logical Function Is_Axis(Ax) Result(Is_Axiss)
    !!--++    character(len=*), intent(in) :: Ax
    !!--++
    !!--++    (PRIVATE)
    !!--++    Detect the presence of a rotation axis
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Is_Axis(Ax) Result(Is_Axiss)
       !---- Argument ----!
       character(len=*), intent(in) :: Ax
       logical                      :: Is_axiss

       !---- Local Variables ----!
       character(len=*), dimension(6), parameter :: axis=(/"1","2","3","4","5","6"/)
       integer                                   :: i

       Is_axiss=.false.
       do i=1,6
          if (Ax == axis(i))  then
             Is_axiss=.true.
             exit
          end if
       end do

       return
    End Function Is_Axis

    !!--++
    !!--++ Logical Function Is_Digit(A) Result(Is_Digitt)
    !!--++    character(len=*), intent(in) :: A    !  In ->
    !!--++
    !!--++    (PRIVATE)
    !!--++    Determine if A is a digit
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Is_Digit(A) Result(Is_Digitt)
       !---- Argument ----!
       character(len=*), intent(in) :: A
       logical                      :: Is_digitt
       character(len=*), parameter  :: digit="0123456789"

       Is_digitt=.false.
       if (index(digit,a) /= 0 ) Is_digitt=.true.

       return
    End Function Is_Digit

    !!--++
    !!--++ Logical Function Is_Hexa(Ng,Ss)
    !!--++    integer, intent (in)                  :: ng   !  In -> Number of Symmetry Operators
    !!--++    integer, dimension(:,:,:), intent(in) :: ss   !  In -> Rotation part of Symmetry Operators  (3,3,48)
    !!--++
    !!--++    (PRIVATE)
    !!--++    Calculate if the SpaceGroup is HEXAGONAL
    !!--++    Valid only for conventional bases
    !!--++
    !!--++  Update: February - 2005
    !!
    Function Is_Hexa(Ng,Ss) Result(Is_Hexag)
       !---- Argument ----!
       integer, intent (in)                   :: ng
       integer, dimension(:,:,:), intent(in)  :: ss   !(3,3,48)
       logical                                :: is_Hexag

       !---- Local Variables ----!
       integer :: i

       Is_Hexag=.false.
       do i=2,ng
          if (sum(abs(ss(:,1,i))) > 1) then
             Is_hexag=.true.
             exit
          end if
          if (sum(abs(ss(:,2,i))) > 1) then
             Is_hexag=.true.
             exit
          end if
       end do

       return
    End Function Is_Hexa

    !!----
    !!---- Logical Function Is_New_Op(Op,N,List_Op) Result(Is_New)
    !!----    type(Sym_Oper_type), intent(in)               :: op      !  In ->  Symmetry operator
    !!----    Integer,             intent(in)               :: n       !  In ->  Integer giving the number of Op i the list
    !!----    type(Sym_Oper_type), intent(in), dimension(:) :: list_op !  In ->  List of n symmetry operators
    !!----
    !!----    Determine if a symmetry operator is or not in a given list
    !!----
    !!---- Update: February - 2005
    !!
    Function Is_New_Op(Op,N,List_Op) Result(Is_New)
       !---- Argument ----!
       type(Sym_Oper_type), intent(in)               :: op
       Integer,             intent(in)               :: n
       type(Sym_Oper_type), intent(in), dimension(:) :: list_op
       logical                                       :: is_new

       !---- Local Variables ----!
       integer :: i

       is_new=.true.
       do i=1,n
          if (op == list_op(i))  then
             is_new=.false.
             exit
          end if
       end do

       return
    End Function Is_New_Op

    !!--++
    !!--++  Logical Function Is_Plane(Ax) Result(Is_Planee)
    !!--++     character(len=*), intent(in) :: Ax
    !!--++
    !!--++     (PRIVATE)
    !!--++     Detect the presence of a mirror or glide plane
    !!--++
    !!--++  Update: February - 2005
    !!
    Function Is_Plane(Ax) Result(Is_Planee)
       !---- Argument ----!
       character(len=*), intent(in) :: Ax
       logical                      :: Is_Planee

       !---- Local Variables ----!
       character(len=*), dimension(6), parameter :: plane=(/"A","B","C","D","M","N"/)
       integer                                   :: i

       Is_planee=.false.
       do i=1,6
          if (Ax == plane(i))  then
             Is_planee=.true.
             exit
          end if
       end do

       return
    End Function Is_Plane

    !!--++
    !!--++ Logical Function Is_Xyz(A) Result(Iss_Xyz)
    !!--++    character(len=*), intent(in) :: A
    !!--++
    !!--++    (PRIVATE)
    !!--++    Determine if A is a character X, Y or Z
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Is_Xyz(A) Result(Iss_Xyz)
       !---- Argument ----!
       character(len=*), intent(in) :: A
       logical                      :: Iss_xyz

       Iss_xyz=.false.
       if (A == "x" .or. A == "X" .or.   &
           A == "y" .or. A == "Y" .or.   &
           A == "z" .or. A == "Z")  Iss_xyz=.true.

       return
    End Function Is_Xyz

    !!----
    !!---- Logical Function Lattice_Trans(V,Lat) Result(Lattice_Transl)
    !!----    real(kind=cp), dimension(3), intent( in) :: v              !  In -> Vector
    !!----    character(len=*),            intent( in) :: Lat            !  In -> Lattice Character
    !!----    logical                                  :: Lattice_Transl ! Out -> .True. or .False.
    !!----
    !!----    Determine whether a vector is a lattice vector
    !!----    depending on the Bravais lattice.
    !!----
    !!---- Update: February - 2005
    !!
    Function Lattice_Trans(V,Lat) Result(Lattice_Transl)
       !---- Argument ----!
       real(kind=cp), dimension(3), intent( in) :: v
       character(len=*),            intent( in) :: Lat
       logical                                  :: Lattice_Transl

       !---- Local variables ----!
       real(kind=cp)   , dimension(3) :: vec
       integer                        :: i

       Lattice_Transl=.false.

       if (Zbelong(v)) then                      ! if v is an integral vector =>  v is a lattice vector
          Lattice_Transl=.true.
       else                                      ! if not look for lattice type
          select case (Lat)
             case("A","a")
                vec=Ltr_a(:,2)-v
                if (Zbelong(vec)) Lattice_Transl=.true.
             case("B","b")
                vec=Ltr_b(:,2)-v
                if (Zbelong(vec)) Lattice_Transl=.true.
             case("C","c")
                vec=Ltr_c(:,2)-v
                if (Zbelong(vec)) Lattice_Transl=.true.
             case("I","i")
                vec=Ltr_i(:,2)-v
                if (Zbelong(vec)) Lattice_Transl=.true.
             case("R","r")
                vec=Ltr_r(:,2)-v
                if (Zbelong(vec)) Lattice_Transl=.true.
                vec=Ltr_r(:,3)-v
                if (Zbelong(vec)) Lattice_Transl=.true.
             case("F","f")
                vec=Ltr_f(:,2)-v
                if (Zbelong(vec)) Lattice_Transl=.true.
                vec=Ltr_f(:,3)-v
                if (Zbelong(vec)) Lattice_Transl=.true.
                vec=Ltr_f(:,4)-v
                if (Zbelong(vec)) Lattice_Transl=.true.
             case("Z")
                do i=2,nlat
                  vec=Ltr(:,i)-v
                  if (Zbelong(vec)) then
                    Lattice_Transl=.true.
                    exit
                  end if
                end do
          end select
       end if

       return
    End Function  Lattice_Trans

    !!--++
    !!--++ Function Product_Symop(Syma,Symb) Result (Symab)
    !!--++    type(Sym_Oper_Type), intent (in) :: syma
    !!--++    type(Sym_Oper_Type), intent (in) :: symb
    !!--++    type(Sym_Oper_Type)              :: symab
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Obtain the symmetry operation corresponding
    !!--++    to the product of two operators by using the * operator.
    !!--++    The calling program can use a statement like: symab=syma*symb
    !!--++
    !!--++  Update: February - 2005
    !!
    Function Product_Symop(Syma,Symb) Result (Symab)
       !---- Arguments ----!
       type(Sym_Oper_Type), intent (in) :: syma
       type(Sym_Oper_Type), intent (in) :: symb
       type(Sym_Oper_Type)              :: symab

       symab%tr  = Syma%tr + matmul(real(Syma%Rot),Symb%tr)
       Symab%Rot = matmul(Syma%Rot,Symb%Rot)

       return
    End Function Product_Symop

    !!----
    !!---- Logical Function Spgr_Equal(Spacegroup1,Spacegroup2) Result(Ispgr_Equal)
    !!----    Type (Space_Group_Type),  intent(in) :: SpaceGroup1   ! In ->
    !!----    Type (Space_Group_Type),  intent(in) :: SpaceGroup2   ! In ->
    !!----
    !!----    Determine if two SpaceGroups are equal
    !!----
    !!---- Update: February - 2005
    !!
    Function Spgr_Equal(Spacegroup1, Spacegroup2) Result(Ispgr_Equal)
       !---- Arguments ----!
       type (Space_Group_Type),  intent(in) :: SpaceGroup1, SpaceGroup2
       logical                              :: iSpGr_Equal

       !---- Trivial tests----!
       iSpGr_Equal=.false.
       if (SpaceGroup1%multip == 0 .or. SpaceGroup2%multip == 0) return
       if (SpaceGroup1%multip /= SpaceGroup2%multip) return

       iSpGr_Equal=Equal_sets_text(SpaceGroup1%SymopSymb,SpaceGroup1%multip, &
                                   SpaceGroup2%SymopSymb,SpaceGroup2%multip)

       return
    End Function Spgr_Equal

    !!----
    !!---- Function Sym_Prod(Syma,Symb,Modlat) Result (Symab)
    !!----    type(Sym_Oper_Type), intent (in) :: syma
    !!----    type(Sym_Oper_Type), intent (in) :: symb
    !!----    logical,optional,    intent (in) :: modlat
    !!----    type(Sym_Oper_Type)              :: symab
    !!----
    !!----    Obtain the symmetry operation corresponding to the product of
    !!----    two operators.
    !!----    If modlat=.true. or it is not present, the traslation
    !!----    part of the resulting operator is reduced to have components <1.0
    !!----
    !!---- Update: February - 2005
    !!
    Function Sym_Prod(Syma,Symb,Modlat) Result (Symab)
       !---- Arguments ----!
       type(Sym_Oper_Type), intent (in) :: syma
       type(Sym_Oper_Type), intent (in) :: symb
       logical,optional,    intent (in) :: modlat
       type(Sym_Oper_Type)              :: symab

       if (present(modlat)) then
          if (.not. modlat) then
             symab%tr = Syma%tr + matmul(real(Syma%Rot),Symb%tr)
          else
             symab%tr = modulo_lat(Syma%tr + matmul(real(Syma%Rot),Symb%tr))
          end if
       else
          symab%tr = modulo_lat(Syma%tr + matmul(real(Syma%Rot),Symb%tr))
       end if
       Symab%Rot = matmul(Syma%Rot,Symb%Rot)

       return
    End Function Sym_Prod

    !!---- Subroutine Allocate_Lattice_Centring(Latt,n,tinv)
    !!----   Type(Lattice_Centring_Type), intent(out)  :: Latt
    !!----   integer,                     intent(in)   :: n
    !!----   logical,  optional,          intent(in)   :: tinv
    !!----
    !!----  Allocates a Lattice_Centring_Type object. If tinv is present and tinv=.true.
    !!----  four indices are selected for the first dimension for storing the presence or
    !!----  absence of time inversion once the object is constructed.
    !!----
    !!----  Updated: October 2014
    !!----
    !!
    Subroutine Allocate_Lattice_Centring(Latt,n,tinv)
      Type(Lattice_Centring_Type), intent(out)  :: Latt
      integer,                     intent(in)   :: n
      logical,  optional,          intent(in)   :: tinv
      !--- Local variables ---!
      if(present(tinv)) then
        if(tinv) then
           allocate(Latt%Ltr(4,n))
        else
           allocate(Latt%Ltr(3,n))
        end if
      else
        allocate(Latt%Ltr(3,n))
      end if
      Latt%Ltr=0.0
      Latt%N_lat=n
      Latt%set=.false.
      return
    End Subroutine Allocate_Lattice_Centring

    !!---- Subroutine Check_Generator(gen,ok,symbol)
    !!----   Character(len=*),         intent(in)  :: gen
    !!----   logical,                  intent(out) :: ok
    !!----   character(len=*),optional,intent(out) :: symbol
    !!----
    !!----  Check that the string containing a generator, contains a legal symmetry operator
    !!----  Only integer coefficients and determinant of the rotational part equal to +1 or -1
    !!----  are allowed. In the optional argument "symbol" the nature of the operator is provided.
    !!----
    !!----  Updated: January 2014
    !!----
    !!
    Subroutine Check_Generator(gen,ok,symbol)
      Character(len=*),         intent(in)  :: gen
      logical,                  intent(out) :: ok
      character(len=*),optional,intent(out) :: symbol
      !--- Local variables ---!
      integer :: i,j,k,n,itr,idet
      character(len=len(gen)), dimension(3) :: split
      character(len=len(gen))  :: symb
      character(len=*), dimension(3), parameter :: code=(/"x","y","z"/)
      real(kind=cp)  :: det
      real(kind=cp), dimension(3,3) :: Mat,iMat
      logical :: esta

      call Init_Err_Symm()
      ok=.false.
      i=index(gen,",")
      j=index(gen,",",back=.true.)
      split(1)= l_case(pack_string(gen(1:i-1)))
      split(2)= l_case(pack_string(gen(i+1:j-1)))
      split(3)= l_case(pack_string(gen(j+1:)))
      !Remove the translation part if it exists
      !write(*,"(4a)") " => Initial split: ", (trim(split(i))//"   ",i=1,3)
      do i=1,3
        n=len_trim(split(i))
        j=index(split(i),"+",back=.true.)
        if(j /= 0) then
          symb=split(i)(j+1:)
          esta=.false.
          do k=1,len_trim(symb)
            if(symb(k:k) == code(1) .or. symb(k:k) == code(2) .or. symb(k:k) == code(3) ) then
               esta = .true.  !A translation is not provided after the matrix
               exit
            end if
          end do
          if(.not. esta) then ! a translation is given in that part of the string, so remove it!
             split(i)=split(i)(1:j-1)
          else ! we have to check starting from the left of the string
             j=index(split(i),"+") !look for the first appearance of "+"
             !Check if there are x,y,z on the left of "+"
             if(j > 1) then
                symb=split(i)(1:j-1)
                esta=.false.
                do k=1,len_trim(symb)
                  if(symb(k:k) == code(1) .or. symb(k:k) == code(2) .or. symb(k:k) == code(3) ) then
                     esta = .true.  !A translation is not provided before the matrix
                     exit
                  end if
                end do
                if(.not. esta) then   !A translation exists
                  split(i)=split(i)(j+1:)
                end if
             end if
          end if
        end if
        if(len_trim(split(i)) == n) then !Check now if instead of "+" the translation is given with "-" sign
          j=index(split(i),"-",back=.true.)
          if(j /= 0) then
            symb=split(i)(j+1:)
            esta=.false.
            do k=1,len_trim(symb)
              if(symb(k:k) == code(1) .or. symb(k:k) == code(2) .or. symb(k:k) == code(3) ) then
                 esta = .true.  !A translation is not provided after the matrix
                 exit
              end if
            end do
            if(.not. esta) then ! a translation is given in that part of the string, so remove it!
               split(i)=split(i)(1:j-1)
            else ! we have to check "-" starting from the left of the string
               j=index(split(i),"-") !look for the first appearance of "+"
               !Check if there are x,y,z on the left of "-"
               if(j > 1) then
                  symb=split(i)(1:j-1)
                  esta=.false.
                  do k=1,len_trim(symb)
                    if(symb(k:k) == code(1) .or. symb(k:k) == code(2) .or. symb(k:k) == code(3) ) then
                       esta = .true.  !A translation is not provided before the matrix
                       exit
                    end if
                  end do
                  if(.not. esta) then   !A translation exists
                    split(i)=split(i)(j+1:)
                  end if
               end if
            end if
          end if
        end if
      end do
      !write(*,"(4a)") " => Final split: ", (trim(split(i))//"   ",i=1,3)
      do i=1,3
       call Get_Num_String(trim(split(i)), Mat(i,:),code)
      end do
      !Now determine if the matrix has integer components
      iMat=real(nint(Mat))
      !now calculate the determinant ... it should be equal to +/-1!
      det=determ_A(Mat)
      idet=nint(det)
      det=abs(det)
      if(present(symbol)) then
        itr=nint(trace(Mat))
        n=0
        select case (itr)
           case (-3)
              if (idet == -1) symbol="-1"
           case (-2)
              if (idet == -1) symbol="-6"
           case (-1)
              if (idet == -1) symbol="-4"
              if (idet ==  1) symbol="2 or 21"
           case ( 0)
              if (idet == -1) symbol="-3"
              if (idet ==  1) symbol="3 or 31/32"
           case ( 1)
              if (idet == -1) symbol="m or g"
              if (idet ==  1) symbol="4 or 41,42..."
           case ( 2)
              if (idet ==  1) symbol="6 or 61,62,..."
           case ( 3)
              if (idet ==  1) symbol="1"
           case default
              n=0
        end select
        symbol=trim(symbol)//"  [undet. loc.]"
      end if
      iMat=iMat-Mat
      if(sum(abs(iMat)) > eps_symm) then
        err_symm=.true.
        err_symm_mess="The matrix corresponding to a generator has non-integer values!"
        return
      else
        if(abs(det-1.0) > eps_symm) then
          err_symm=.true.
          err_symm_mess="The matrix corresponding to a generator has a determinant with absolute value different of 1.0"
          return
        end if
      end if
      ok=.true.  !arriving here the generator is ok!
      return
    End Subroutine Check_Generator

    !---- Subroutines ----!

    !!--++
    !!--++ Subroutine Check_Symbol_Hm(Hms)
    !!--++    character (len=1), dimension(3,4), intent( in):: HMS   ! In -> Hermann-Mauguin Symbol
    !!--++
    !!--++    (PRIVATE)
    !!--++    Subroutine used by Get_SO_from_HMS.
    !!--++    Check the correctness of the Herman-Mauguin Symbol (not all!!!).
    !!--++    Logical "hexa" must be defined and control error is present.
    !!--++
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Check_Symbol_Hm(Hms)
       !---- Argument ----!
       character (len=1), dimension(3,4), intent( in):: HMS

       !---- Local Variables ----!
       logical          :: is_there,axis,plane
       character(len=1) :: Item_SP
       character(len=*), dimension(16), parameter ::                    &
                         good=(/"1","2","3","4","5","6","A","B","C","D","M","N","P","/","-"," "/)
       integer          :: ncount,five,i,j,l

       !---- Check for missprinted symbols ----!
       call init_err_symm()
       do i=1,3
          do j=1,4
             is_there=.false.
             five=0
             if (HMS(i,j) == "5") five=j
             do L=1,16
                if (HMS(i,j) == good(L)) is_there=.true.
             end do
             if (.not. is_there) then
                err_symm=.true.
                ERR_Symm_Mess=" The symbol: "//HMS(i,j)//" is not allowed"
                return
             else if (five == 1) then
                err_symm=.true.
                ERR_Symm_Mess=" The fivefold axis is not allowed"
                return
             end if
          end do
       end do

       !---- Check for repetitions and axes followed by planes (and viceversa) ----!
       do i=1,3
          do j=1,3
             Item_SP=HMS(i,j)
             if (Item_SP == " ") cycle
             is_there=.false.
             axis=Is_axis(Item_SP)
             plane=Is_plane(Item_SP)
             do L=j+1,4
                if (HMS(i,L) == Item_SP)  is_there=.true.
             end do
             if (is_there) then
                err_symm=.true.
                ERR_Symm_Mess=" The symbol: "//HMS(i,j)// &
                              " has been repeated within the same symmetry direction"
                return
             end if
             if (axis .and. Is_plane(HMS(i,j+1))) then
                err_symm=.true.
                ERR_Symm_Mess=" A rotation axis cannot be immediately followed by a plane"//char(13)//&
                              " within the same symmetry direction"
                return
             end if
             if (plane .and. Is_axis(HMS(i,j+1))) then
                err_symm=.true.
                ERR_Symm_Mess=" A mirror plane cannot be immediately followed by a rotation axis"//char(13)//&
                              " within the same symmetry direction"
                return
             end if
          end do
       end do

       !---- Check for two planes in the same symmetry direction ----!
       do i=1,3
          ncount=0
          do j=1,4
             Item_SP=HMS(i,j)
             do L=7,12
                if (good(L) == Item_SP) ncount=ncount+1
             end do
          end do
          if (ncount > 1) then
             err_symm=.true.
             ERR_Symm_Mess=" There is more than one plane within the same symmetry direction"
             return
          end if
       end do

       !---- Check for ILLEGAL screw axes ----!
       do i=1,3
          ncount=0
          do j=1,4
             Item_SP=HMS(i,j)
             if (Item_SP == " ") cycle
             do L=1,6
                if (good(L) == Item_SP) ncount=ncount+1
             end do
          end do
          if (ncount > 1) then  !there is more than one axis-symbol -> Screw
          !   if (iachar(HMS(i,1)) < iachar(HMS(i,2))) then
             if (HMS(i,1) <  HMS(i,2) ) then
                err_symm=.true.
                ERR_Symm_Mess=" Screw axis: "//HMS(i,1)//" "//HMS(i,2)//" not allowed"
                return
             end if
          end if
       end do

       return
    End Subroutine Check_Symbol_HM

    !!----
    !!---- Subroutine Decodmatmag(Sim,Xyzstring)
    !!----    integer, dimension(3,3), intent(in)  :: sim          !  In -> Rotation matrix
    !!----    character (len=*),       intent(out) :: XYZstring    ! Out -> String (Mx,My,Mz)
    !!----
    !!----    Supplies a string of the form (Mx,My,Mz) for the rotation matrix Sim.
    !!----    Logical "hexa" must be defined.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Decodmatmag(Sim,Xyzstring)
       !---- Arguments ----!
       integer,dimension (3,3), intent( in) :: sim
       character (len=*),       intent(out) :: XYZstring

       !---- Local variables ----!
       integer :: Iu,j,ihex

       call SearchOp(sim,1,36,Iu)

       if (.not. hexa) then
          j=abs(Iu)
          if (Iu < 0) j=j+24
          XYZstring= MAGmat(j)
       else
          j=abs(Iu)-24
          ihex=2
          if ( j < 0 ) then
             j=j+24
             ihex=1
          end if
          if (Iu < 0) j=j+24/ihex
          XYZstring= MAGmat(j+(ihex-1)*48)
       end if

       return
    End Subroutine DecodMatMag

    !!----
    !!---- Subroutine Get_Centring_Vectors(L,Latc,LatSymb)
    !!----    integer,                        intent (in out) :: L       ! Number of centring vectors
    !!----    real(kind=cp), dimension(:,:),  intent (in out) :: latc    ! Centering vectors
    !!----    character(len=1),               intent (   out) :: LatSymb ! Lattice symbol
    !!----
    !!----    Subroutine to complete the centring vectors of a centered lattice and to provide a lattice symbol.
    !!----    It is useful when non-conventional lattices are used to obtain all lattice
    !!----    translations with components in the range [0.0 1.0). The (0,0,0) translation
    !!----    is removed in case it comes on input.
    !!----
    !!---- Update: February - 2005, January-2014 (JRC)
    !!
    Subroutine Get_Centring_Vectors(L,Latc,LatSymb)
       !---- Arguments ----!
       integer,                       intent (in out) :: L
       real(kind=cp), dimension(:,:), intent (in out) :: latc  !(3,n)
       character(len=*),              intent (out)    :: LatSymb
       !---- Local variables ----!
       logical                                  :: isnew
       real(kind=cp), dimension(3,size(latc,2)) :: latinv,newlat
       real(kind=cp), dimension(3)              :: v,v1,v2
       integer                                  :: i,j,k1,k2,n,lat_ini,lm
       real(kind=cp), parameter                 :: teps=3.0*eps_symm

       LatSymb="P"
       if(L == 0) return
       newlat=latc
       !Purge the translations
       do i=1,L-1
         v=newlat(:,i)
         if(sum(v) < teps) cycle
         do j=i+1,L
            if(sum(abs(newlat(:,j)-v)) < teps) newlat(:,j)=0.0
         end do
       end do
       n=0
       do i=1,L
         if(sum(abs(newlat(:,i))) < teps) cycle
         n=n+1
         latc(:,n)=newlat(:,i)
       end do
       L=n  !normally n < L_initial

       latinv=0.0
       where (abs(latc)> teps)
          latinv=1.0/latc
       end where
       do
          lat_ini=L
          do i=1,L    !Even for a single centring vector this loop is executed
            v1=latc(:,i)
            do j=i,L  !start on i to ensure that for a single centring vector the internal loops are executed
              v2=latc(:,j)
              do k1=0,maxval(nint(latinv(:,i)))
                do k2=0,maxval(nint(latinv(:,j)))
                  v=modulo_lat(k1*v1+k2*v2)
                  if(sum(abs(v)) < teps) cycle
                  if( any(v > 1.0-teps) ) cycle
                  isnew=.true.
                  do lm=1,L
                    if (sum(abs(v-latc(:,lm))) < teps) then
                       isnew=.false.
                       exit
                    end if
                  end do
                  if(isnew) then
                     L=L+1
                     latc(:,L)=v
                  end if
                end do
              end do
            end do
          end do
          If(L == Lat_ini) exit !No more vectors have been generated
       end do

       !Recognize the type of Lattice
       Select Case(L)

         Case(1) !Test I, A, B, C
            if(sum(abs(latc(:,1)-Ltr_i(:,2))) < teps) then
              LatSymb="I"
              return
            end if
            if(sum(abs(latc(:,1)-Ltr_a(:,2))) < teps) then
              LatSymb="A"
              return
            end if
            if(sum(abs(latc(:,1)-Ltr_b(:,2))) < teps) then
              LatSymb="B"
              return
            end if
            if(sum(abs(latc(:,1)-Ltr_c(:,2))) < teps) then
              LatSymb="C"
              return
            end if

         Case(2)  !Test R
             isnew=.false.
             if(sum(abs(latc(:,1)-Ltr_r(:,2))) < teps .or. sum(abs(latc(:,1)-Ltr_r(:,3))) < teps) isnew=.true.
             if(isnew) then
               if(sum(abs(latc(:,2)-Ltr_r(:,2))) < teps .or. sum(abs(latc(:,2)-Ltr_r(:,3))) < teps) then
                 LatSymb="R"
                 return
               end if
             end if

         Case(3)
             isnew=.false.
             do i=2,4
                if (  sum(abs(latc(:,1)-Ltr_f(:,i))) < teps  ) then
                   isnew=.true.
                   exit
                end if
             end do
             if(isnew) then
                isnew=.false.
                do i=2,4
                   if (  sum(abs(latc(:,2)-Ltr_f(:,i))) < teps  ) then
                      isnew=.true.
                      exit
                   end if
                end do
             end if
             if(isnew) then
                isnew=.false.
                do i=2,4
                   if (  sum(abs(latc(:,3)-Ltr_f(:,i))) < teps  ) then
                       LatSymb="F"
                       return
                   end if
                end do
             end if

       End Select
       LatSymb="Z"
       return
    End Subroutine Get_Centring_Vectors

    !!----
    !!---- Subroutine Get_Crystal_System(Ng, Ss / Gen, Isystm, Crys)
    !!----    integer,                      intent(in) :: Ng     !  In -> Number of Operators (not related by
    !!----                                                                inversion and lattice traslations)
    !!----    integer, dimension(:,:,:),    intent(in) :: Ss     !  In -> Rotation Part   (3,3,48)
    !!----    or
    !!----    character(len=*),dimension(:),intent(in) :: gen    !  In -> Jones Faithful form of symmetry operators
    !!----    integer,                      intent(out):: ISystm ! Out -> Number for Crystal System
    !!----                                                                 1: Triclinic       2: Monoclinic
    !!----                                                                 3: Orthorrombic    4: Tetragonal
    !!----                                                                 5: Trigonal        6: Hexagonal
    !!----                                                                 7: Cubic
    !!----    character(len=1),             intent(out):: Crys   ! Out -> Symbol of Crystal family
    !!----
    !!----    Obtain the number and string of the Crystal System from a set of operators
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine Get_Crystal_System_R_OP(Ng, Ss, Isystm, Crys)
    !!--++    integer,                   intent(in) :: Ng       !  In -> Number of Operators (not related by
    !!--++                                                               inversion and lattice traslations)
    !!--++    integer, dimension(:,:,:), intent(in) :: Ss       !  In -> Rotation Part   (3,3,48)
    !!--++    integer,                   intent(out):: ISystm   ! Out -> Number for Crystal System
    !!--++                                                                1: Triclinic       2: Monoclinic
    !!--++                                                                3: Orthorrombic    4: Tetragonal
    !!--++                                                                5: Trigonal        6: Hexagonal
    !!--++                                                                7: Cubic
    !!--++    character(len=1),          intent(out):: Crys     ! Out -> Symbol of Crystal family
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Obtain the number and string of the Crystal System from a set of operators
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Get_Crystal_System_R_OP(Ng,Ss, Isystm, Crys)
       !---- Arguments ----!
       integer,                   intent(in) :: Ng
       integer, dimension(:,:,:), intent(in) :: Ss    !(3,3,48)
       integer,                   intent(out):: ISystm
       character(len=1),          intent(out):: Crys

       !---- Local variables ----!
       integer   :: i, ndet
       integer   :: nrot_1, nrot_2, nrot_3, nrot_4, nrot_6
       integer   :: nrot_1b, nrot_2b, nrot_3b, nrot_4b, nrot_6b

       nrot_1  = 0
       nrot_2  = 0
       nrot_3  = 0
       nrot_4  = 0
       nrot_6  = 0
       nrot_1b = 0
       nrot_2b = 0
       nrot_3b = 0
       nrot_4b = 0
       nrot_6b = 0

       do i=1,ng
          ndet= Axes_Rotation(ss(:,:,i))
          select case (ndet)
              case (-6)
                 nrot_6b=nrot_6b +1
              case (-4)
                 nrot_4b=nrot_4b +1
              case (-3)
                 nrot_3b=nrot_3b +1
              case (-2)
                 nrot_2b=nrot_2b +1
              case (-1)
                 nrot_1b=nrot_1b +1
              case ( 1)
                 nrot_1 =nrot_1  +1
              case ( 2)
                 nrot_2 =nrot_2  +1
              case ( 3)
                 nrot_3 =nrot_3  +1
              case ( 4)
                 nrot_4 =nrot_4  +1
              case ( 6)
                 nrot_6 =nrot_6  +1
              case default
                 err_symm=.true.
                 ERR_Symm_Mess= " Axes rotation wrong"
                 return
          end select
       end do

       !---- Cubic ----!
       if ( (nrot_3 + nrot_3b == 8) ) then
          isystm = 7
          crys="c"

       !---- Hexagonal ----!
       else if ( (nrot_6 + nrot_6b == 2) ) then
          isystm = 6
          crys="h"

       !---- Trigonal ----!
       else if ( (nrot_3 + nrot_3b == 2) ) then
          isystm = 5
          crys="h"

       !---- Tetragonal ----!
       else if ( (nrot_4 + nrot_4b == 2) ) then
          isystm = 4
          crys="t"

       !---- Orthorhombic ----!
       else if ( (nrot_2 + nrot_2b == 3) ) then
          isystm = 3
          crys="o"

       !---- Monoclinic  ----!
       else if ( (nrot_2 + nrot_2b == 1) ) then
          isystm = 2
          crys="m"

       !---- Triclinic  ----!
       else
          isystm = 1
          crys="a"

       end if

       return
    End Subroutine Get_Crystal_System_R_OP

    !!--++
    !!--++ Subroutine Get_Crystal_System_R_ST(Ng,Gen,Isystm, Crys)
    !!--++    integer,                      intent(in) :: Ng     !  In -> Number of Operators
    !!--++    character(len=*),dimension(:),intent(in) :: gen    !  In -> Jones Faithful form of symmetry operators
    !!--++    integer,                      intent(out):: ISystm ! Out -> Number for Crystal System
    !!--++                                                                1: Triclinic       2: Monoclinic
    !!--++                                                                3: Orthorrombic    4: Tetragonal
    !!--++                                                                5: Trigonal        6: Hexagonal
    !!--++                                                                7: Cubic
    !!--++    character(len=1),             intent(out):: Crys   ! Out -> Symbol of Crystal family
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Obtain the number and string of the Crystal System from a set of operators
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Get_Crystal_System_R_ST(Ng,gen, Isystm, Crys)
       !---- Arguments ----!
       integer,                        intent(in) :: Ng
       character(len=*), dimension(:), intent(in) :: Gen
       integer,                        intent(out):: ISystm
       character(len=1),               intent(out):: Crys

       !---- Local variables ----!
       integer, dimension(3,3,Ng) :: Ss    !(3,3,48)
       integer                    :: i

       do i=1,Ng
          call Read_Xsym(gen(i),1,Ss(:,:,i))
       end do
       call Get_Crystal_System_R_OP(Ng,Ss, Isystm, Crys)

       return
    End Subroutine Get_Crystal_System_R_ST

    !!----
    !!---- Subroutine Get_GenSymb_from_Gener(gen,ngen, SpaceH)
    !!----    character(len=*),dimension(:),  intent(in) :: gen     !  In -> list of generators is string mode
    !!----    integer,                        intent(in) :: ngen    !  In -> number of generators provided
    !!----    character(len=*),              intent(out) :: SpaceH  ! Out -> Generalised Hall Symbol
    !!----
    !!----    Determines a generalised Hall symbol for a space group formed by the symmetry symbols of
    !!----    the provided generators.
    !!----
    !!---- Updated: January - 2014 (JRC)
    !!
    Subroutine Get_GenSymb_from_Gener(gen,ngen,SpaceH)
       !---- Arguments ----!
       character(len=*),dimension(:),  intent(in) :: gen
       integer,                        intent(in) :: ngen
       character(len=*),              intent(out) :: SpaceH

       !----Local variables ----!
       character(len= 1)          :: LatSymb
       character(len=20)          :: centr
       character(len=40)          :: gen_symb
       integer                    :: ng, ini, i, orden,L,j
       integer, dimension(3,3,24) :: ss
       integer, dimension(3,3)    :: nulo,unitm

       real(kind=cp), dimension(3,24)  :: ts
       real(kind=cp), dimension(3,192) :: latc
       real(kind=cp), dimension(3)     :: ts_centre
       logical                         :: centred

       !---- Initial Values ----!
       nulo=0
       unitm=0
       unitm(1,1)=1;  unitm(2,2)=1;  unitm(3,3)=1
       latc=0.0
       centred=.false.
       centr=" "
       SpaceH=" "

       ! --- Test if lattice translations are provide with a symbol in the first generator
       if(index(gen(1),"-I") /= 0) then       !Centric with -1 at 000
         SpaceH="-I"
       else if(index(gen(1),"-A") /= 0) then
         SpaceH="-A"
       else if(index(gen(1),"-B") /= 0) then
         SpaceH="-B"
       else if(index(gen(1),"-C") /= 0) then
         SpaceH="-C"
       else if(index(gen(1),"-R") /= 0) then
         SpaceH="-R"
       else if(index(gen(1),"-F") /= 0) then
         SpaceH="-F"
       else if(index(gen(1),"-Z") /= 0) then
         SpaceH="-Z"
       else if(index(gen(1),"-P") /= 0) then
         SpaceH="-P"
       end if
       if(len_trim(SpaceH) == 0) then           !centric with -1 out of 000 or acentric
           if(index(gen(1),"I") /= 0) then
             SpaceH="I"
           else if(index(gen(1),"A") /= 0) then
             SpaceH="A"
           else if(index(gen(1),"B") /= 0) then
             SpaceH="B"
           else if(index(gen(1),"C") /= 0) then
             SpaceH="C"
           else if(index(gen(1),"R") /= 0) then
             SpaceH="R"
           else if(index(gen(1),"F") /= 0) then
             SpaceH="F"
           else if(index(gen(1),"P") /= 0) then
             SpaceH="P"
           else if(index(gen(1),"Z") /= 0) then
             SpaceH="Z"
           end if
       end if
       LatSymb="P"
       if(len_trim(SpaceH) == 0) then
         ini=1  !If there is a centring lattice it must be given in the list of the generators
       else
         ini=2
         if(len_trim(SpaceH)==1) LatSymb=trim(SpaceH)
       end if
       ng=0
       do i=ini,ngen
         ng=ng+1
         call Read_Xsym(gen(i),1,ss(:,:,ng),ts(:,ng))
       end do
       !Look for lattice translations as generators
       if(ini == 1) then
         L=0
         do i=1,ng
           if(equal_matrix(ss(:,:,i),unitm,3)) then
             L=L+1
             latc(:,L)=ts(:,i)
             ss(:,:,i)=0
           end if
         end do
         if(L > 0) then !There are lattice translations
           call Get_Centring_Vectors(L,latc,LatSymb)
         end if
       end if
       !Look for centre of symmetry as generator
       do i=1,ng
         if(equal_matrix(ss(:,:,i),-unitm,3)) then !Centre of symmetry
           ts_centre=ts(:,i)
           ss(:,:,i)=0
           centred=.true.
           exit
         end if
       end do
       if(centred) then
         if(sum(abs(ts_centre)) < eps_symm) then
             SpaceH="-"//LatSymb
          else
            ts_centre=0.5*ts_centre
            call Frac_Trans_2Dig(ts_centre,centr)
            centr="-1"//trim(centr)
          end if
       else
         if(ini ==1) SpaceH=LatSymb
       end if
       !Construct the symbol
       do i=1,ng
          if(equal_matrix(ss(:,:,i), nulo,3)) cycle
          call symmetry_symbol(ss(:,:,i),ts(:,i),gen_symb)

          if(len_trim(gen_symb) == 0) then
            orden=axes_rotation(ss(:,:,i))
            write(unit=gen_symb,fmt="(i2)") orden
            gen_symb=adjustl(gen_symb)//"[]"
          else
            j=index(gen_symb,")")
            if( j /= 0) then
               gen_symb=gen_symb(1:j)
               j=index(gen_symb,"+")
               gen_symb(j:j)=" "
               gen_symb=pack_string(gen_symb)
            else
               j=index(gen_symb," ")
               gen_symb=gen_symb(1:j)
               j=index(gen_symb,"+")
               gen_symb=gen_symb(1:j-1)
            end if
          end if
          SpaceH=trim(SpaceH)//" "//trim(gen_symb)
       end do
       SpaceH=trim(SpaceH)//" "//trim(centr)
       return
    End Subroutine Get_GenSymb_from_Gener

    !!----
    !!---- Subroutine Get_HallSymb_From_Gener(Spacegroup, Spaceh)
    !!----    type(Space_Group_Type),   intent(in out) :: SpaceGroup   !  In -> SpaceGroup type variable
    !!----                                                               Out -> SpaceGroup type variable
    !!----    character(len=*), intent(out), optional  :: SpaceH       ! Out -> Hall Symbol
    !!----
    !!----    Determines the Hall symbol. In general this routine try to obtain
    !!----    the Hall symbol from generators so you need call Get_So_from_Gener
    !!----    before and call Set_Spgr_Info.It doesn't work for arbitrary settings.
    !!----    If one wants to use arbitrary settings the subroutine Get_GenSymb_from_Gener
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Get_HallSymb_from_Gener(SpaceGroup,SpaceH)
       !---- Arguments ----!
       type(Space_Group_Type), intent(in out)        :: SpaceGroup
       character(len=*),       intent(out), optional :: SpaceH

       !----Local variables ----!
       character(len= 1)        :: axes,axes2
       character(len= 3)        :: tras
       character(len=20)        :: Hall
       character(len=*), dimension(13), parameter :: traslacion =&
                           (/"N","A","B","C","D","U","V","W","1","2","3","4","5"/)

       integer                    :: ng,ngen, ini, i, j, k, orden, nt, npos
       integer, dimension(3)      :: tt, tt1, tt2, tt3
       integer, dimension(6)      :: norden
       integer, dimension(3,3,24) :: ss
       integer, dimension(3,3)    :: ss1
       integer, dimension(3,6), parameter :: lattice=reshape((/0,6,6, 6,0,6, &
                                                     6,6,0, 6,6,6, 8,4,4, 4,8,8/),(/3,6/))
       integer, dimension(3,13), parameter :: tras_val=reshape((/6,6,6, 6,0,0, &
                                       0,6,0, 0,0,6, 3,3,3, 3,0,0, 0,3,0, 0,0,3, &
                                       1,0,0, 2,0,0, 3,0,0, 4,0,0, 5,0,0/),(/3,13/))
       integer, dimension(3,3), parameter  :: x_1   = reshape( &
                                 (/ 1, 0, 0,  0, 1, 0,  0, 0, 1/),(/3,3/))
       integer, dimension(3,3), parameter  :: z_1   = reshape( &
                                 (/ 1, 0, 0,  0, 1, 0,  0, 0, 1/),(/3,3/))
       integer, dimension(3,3), parameter  :: x_2   = reshape( &
                                 (/ 1, 0, 0,  0,-1, 0,  0, 0,-1/),(/3,3/))
       integer, dimension(3,3), parameter  :: y_2   = reshape( &
                                 (/-1, 0, 0,  0, 1, 0,  0, 0,-1/),(/3,3/))
       integer, dimension(3,3), parameter  :: z_2   = reshape( &
                                 (/-1, 0, 0,  0,-1, 0,  0, 0, 1/),(/3,3/))
       integer, dimension(3,3), parameter  :: x_3   = reshape( &
                                 (/ 1, 0, 0,  0, 0, 1,  0,-1,-1/),(/3,3/))
       integer, dimension(3,3), parameter  :: y_3   = reshape( &
                                 (/-1, 0,-1,  0, 1, 0,  1, 0, 0/),(/3,3/))
       integer, dimension(3,3), parameter  :: z_3   = reshape( &
                                 (/ 0, 1, 0, -1,-1, 0,  0, 0, 1/),(/3,3/))
       integer, dimension(3,3), parameter  :: x_4   = reshape( &
                                 (/ 1, 0, 0,  0, 0, 1,  0,-1, 0/),(/3,3/))
       integer, dimension(3,3), parameter  :: y_4   = reshape( &
                                 (/ 0, 0,-1,  0, 1, 0,  1, 0, 0/),(/3,3/))
       integer, dimension(3,3), parameter  :: z_4   = reshape( &
                                 (/ 0, 1, 0, -1, 0, 0,  0, 0, 1/),(/3,3/))
       integer, dimension(3,3), parameter  :: x_6   = reshape( &
                                 (/ 1, 0, 0,  0, 1, 1,  0,-1, 0/),(/3,3/))
       integer, dimension(3,3), parameter  :: y_6   = reshape( &
                                 (/ 0, 0,-1,  0, 1, 0,  1, 0, 1/),(/3,3/))
       integer, dimension(3,3), parameter  :: z_6   = reshape( &
                                 (/ 1, 1, 0, -1, 0, 0,  0, 0, 1/),(/3,3/))
       integer, dimension(3,3), parameter  :: x_2p  = reshape( &
                                 (/-1, 0, 0,  0, 0,-1,  0,-1, 0/),(/3,3/))
       integer, dimension(3,3), parameter  :: y_2p  = reshape( &
                                 (/ 0, 0,-1,  0,-1, 0, -1, 0, 0/),(/3,3/))
       integer, dimension(3,3), parameter  :: z_2p  = reshape( &
                                 (/ 0,-1, 0, -1, 0, 0,  0, 0,-1/),(/3,3/))
       integer, dimension(3,3), parameter  :: x_2pp = reshape( &
                                 (/-1, 0, 0,  0, 0, 1,  0, 1, 0/),(/3,3/))
       integer, dimension(3,3), parameter  :: y_2pp = reshape( &
                                 (/ 0, 0, 1,  0,-1, 0,  1, 0, 0/),(/3,3/))
       integer, dimension(3,3), parameter  :: z_2pp = reshape( &
                                 (/ 0, 1, 0,  1, 0, 0,  0, 0,-1/),(/3,3/))
       integer, dimension(3,3), parameter  :: xyz_3 = reshape( &
                                 (/ 0, 1, 0,  0, 0, 1,  1, 0, 0/),(/3,3/))
       integer, dimension(4,4), parameter :: nulo      = reshape((/0, 0, 0, 0, &
                                                                   0, 0, 0, 0, &
                                                                   0, 0, 0, 0, &
                                                                   0, 0, 0, 0/),(/4,4/))

       real(kind=cp), dimension(3,24)          :: ts
       real(kind=cp), dimension(3)             :: ts1
       type (Gener_Oper_Type),dimension(5) :: generador

       !---- Initial Values ----!
       Hall          = " "
       ngen          = 0
       generador = gener_oper_type(0," "," "," ")

       !---- Load Operators ----!
       ng=SpaceGroup%NumOps
       do i=1,ng
          ss(:,:,i) = SpaceGroup%Symop(i)%rot
          ts(:,  i) = SpaceGroup%Symop(i)%tr
       end do

       !---- Tipo de Red ----!
       select case (SpaceGroup%centred)
          case (0,1)
             hall(1:1)=SpaceGroup%SPG_lat
          case (2)
             hall(1:2)="-"//SpaceGroup%SPG_lat
       end select

       !---- Origen del Centro de inversion ----!
       if (SpaceGroup%centred == 0) then
          ngen=ngen+1
          ini=1
          tras=" "
          tt=nint(12.0*2.0*SpaceGroup%centre_coord)

          select case (SpaceGroup%Bravais)
             case ("A")
                tt1=mod(tt-lattice(:,1)+48,12)
                if (sum(tt1) < sum(tt) ) tt=tt1

             case ("B")
                tt1=mod(tt-lattice(:,2)+48,12)
                if (sum(tt1) < sum(tt) ) tt=tt1

             case ("C")
                tt1=mod(tt-lattice(:,3)+48,12)
                if (sum(tt1) < sum(tt) ) tt=tt1

             case ("I")
                tt1=mod(tt-lattice(:,4)+48,12)
                if (sum(tt1) < sum(tt) ) tt=tt1

             case ("R")
                tt1=mod(tt-lattice(:,5)+48,12)
                tt2=mod(tt-lattice(:,6)+48,12)
                if (sum(tt1) < sum(tt) ) tt=tt1
                if (sum(tt2) < sum(tt) ) tt=tt2

             case ("F")
                tt1=mod(tt-lattice(:,3)+48,12)
                tt2=mod(tt-lattice(:,2)+48,12)
                tt3=mod(tt-lattice(:,1)+48,12)
                if (sum(tt1) < sum(tt) ) tt=tt1
                if (sum(tt2) < sum(tt) ) tt=tt2
                if (sum(tt3) < sum(tt) ) tt=tt3

          end select

          do k=1,3     ! 3 pasadas
             if (tt(1) == 0 .and. tt(2) == 0 .and. tt(3) == 0) exit
             do j=1,8
                tt1=tt-tras_val(:,j)
                if ( all(tt1 >= 0) ) then
                   tras(ini:ini)=l_case(traslacion(j))
                   tt=tt1
                   ini=ini+1
                   exit
                end if
             end do
          end do

          generador(ngen)%orden= -1
          generador(ngen)%axes = "z"
          generador(ngen)%tras=tras
       end if

       !---- Suppress non needed contributions ----!
       do i=1,ng
          if (equal_matrix(ss(1:3,1:3,i), x_2,3) .or.  equal_matrix(ss(1:3,1:3,i),-x_2,3) .or. &
              equal_matrix(ss(1:3,1:3,i), y_2,3) .or.  equal_matrix(ss(1:3,1:3,i),-y_2,3) .or. &
              equal_matrix(ss(1:3,1:3,i), z_2,3) .or.  equal_matrix(ss(1:3,1:3,i),-z_2,3) .or. &
              equal_matrix(ss(1:3,1:3,i), x_3,3) .or.  equal_matrix(ss(1:3,1:3,i),-x_3,3) .or. &
              equal_matrix(ss(1:3,1:3,i), y_3,3) .or.  equal_matrix(ss(1:3,1:3,i),-y_3,3) .or. &
              equal_matrix(ss(1:3,1:3,i), z_3,3) .or.  equal_matrix(ss(1:3,1:3,i),-z_3,3) .or. &
              equal_matrix(ss(1:3,1:3,i), x_4,3) .or.  equal_matrix(ss(1:3,1:3,i),-x_4,3) .or. &
              equal_matrix(ss(1:3,1:3,i), y_4,3) .or.  equal_matrix(ss(1:3,1:3,i),-y_4,3) .or. &
              equal_matrix(ss(1:3,1:3,i), z_4,3) .or.  equal_matrix(ss(1:3,1:3,i),-z_4,3) .or. &
              equal_matrix(ss(1:3,1:3,i), x_6,3) .or.  equal_matrix(ss(1:3,1:3,i),-x_6,3) .or. &
              equal_matrix(ss(1:3,1:3,i), y_6,3) .or.  equal_matrix(ss(1:3,1:3,i),-y_6,3) .or. &
              equal_matrix(ss(1:3,1:3,i), z_6,3) .or.  equal_matrix(ss(1:3,1:3,i),-z_6,3) .or. &
              equal_matrix(ss(1:3,1:3,i), x_2p,3).or.  equal_matrix(ss(1:3,1:3,i),-x_2p,3).or. &
              equal_matrix(ss(1:3,1:3,i), y_2p,3).or.  equal_matrix(ss(1:3,1:3,i),-y_2p,3).or. &
              equal_matrix(ss(1:3,1:3,i), z_2p,3).or.  equal_matrix(ss(1:3,1:3,i),-z_2p,3).or. &
              equal_matrix(ss(1:3,1:3,i), x_2pp,3).or. equal_matrix(ss(1:3,1:3,i),-x_2pp,3).or. &
              equal_matrix(ss(1:3,1:3,i), y_2pp,3).or. equal_matrix(ss(1:3,1:3,i),-y_2pp,3).or. &
              equal_matrix(ss(1:3,1:3,i), z_2pp,3).or. equal_matrix(ss(1:3,1:3,i),-z_2pp,3).or. &
              equal_matrix(ss(1:3,1:3,i), xyz_3,3).or. equal_matrix(ss(1:3,1:3,i),-xyz_3,3) ) cycle

              ss(:,:,i)=0
       end do

       !---- Ordering following Order of rotations ----!
       norden=0
       do i=1,ng
          if (equal_matrix(ss(1:3,1:3,i),nulo(1:3,1:3),3)) cycle
           orden=axes_rotation(ss(1:3,1:3,i))
          norden(abs(orden))=norden(abs(orden))+1
       end do

       npos=0
       do j=6,1,-1
          if (norden(j) == 0) cycle
          do i=1,ng
             if (equal_matrix(ss(1:3,1:3,i),nulo(1:3,1:3),3)) cycle
              orden=axes_rotation(ss(1:3,1:3,i))
             if (abs(orden) == j) then
                ss1=ss(:,:,i)
                ts1=ts(:,i)
                npos=npos+1
                ss(:,:,i)=ss(:,:,npos)
                ts(:,i)  =ts(:,npos)
                ss(:,:,npos)=ss1
                ts(:,npos)  =ts1
             end if
          end do
       end do
       nt=npos

       if (nt == 0) then
          ngen=ngen+1
          generador(ngen)%orden= 1
          generador(ngen)%axes = "z"
          generador(ngen)%tras= " "
       end if

       !---- Ordering following axes Direction ----!
       norden=0
       do i=1,nt
          if (equal_matrix(ss(1:3,1:3,i),nulo(1:3,1:3),3)) cycle
          orden=axes_rotation(ss(1:3,1:3,i))
          norden(abs(orden))=norden(abs(orden))+1
       end do
       if (norden(6) > 0) norden(3)=0

       do i=1,nt
          orden=0
          axes =" "
          axes2=" "
          tras ="  "
          orden=axes_rotation(ss(1:3,1:3,i))
          if (norden(abs(orden)) == 0) cycle
          select case (abs(orden))
              case (1)
                 if (orden > 0) then
                    if (equal_matrix(ss(1:3,1:3,i),z_1,3)) then
                       axes="z"
                    end if
                 else
                    if (equal_matrix(ss(1:3,1:3,i),-z_1,3)) then
                       axes="z"
                    end if
                 end if

              case (2)
                 if (orden > 0) then
                    if (equal_matrix(ss(1:3,1:3,i),x_2,3)) then
                       axes="x"
                    else if (equal_matrix(ss(1:3,1:3,i),y_2,3)) then
                       axes="y"
                    else if (equal_matrix(ss(1:3,1:3,i),z_2,3)) then
                       axes="z"
                    else if (equal_matrix(ss(1:3,1:3,i),x_2p,3)) then
                       axes="'"
                       axes2="x"
                    else if (equal_matrix(ss(1:3,1:3,i),y_2p,3)) then
                       axes="'"
                       axes2="y"
                    else if (equal_matrix(ss(1:3,1:3,i),z_2p,3)) then
                       axes="'"
                       axes2="z"
                    else if (equal_matrix(ss(1:3,1:3,i),x_2pp,3)) then
                       axes=""""
                       axes2="x"
                    else if (equal_matrix(ss(1:3,1:3,i),y_2pp,3)) then
                       axes=""""
                       axes2="y"
                    else if (equal_matrix(ss(1:3,1:3,i),z_2pp,3)) then
                       axes=""""
                       axes2="z"
                    end if
                 else
                    if (equal_matrix(ss(1:3,1:3,i),-x_2,3)) then
                       axes="x"
                    else if (equal_matrix(ss(1:3,1:3,i),-y_2,3)) then
                       axes="y"
                    else if (equal_matrix(ss(1:3,1:3,i),-z_2,3)) then
                       axes="z"
                    else if (equal_matrix(ss(1:3,1:3,i),-x_2p,3)) then
                       axes="'"
                       axes2="x"
                    else if (equal_matrix(ss(1:3,1:3,i),-y_2p,3)) then
                       axes="'"
                       axes2="y"
                    else if (equal_matrix(ss(1:3,1:3,i),-z_2p,3)) then
                       axes="'"
                       axes2="z"
                    else if (equal_matrix(ss(1:3,1:3,i),-x_2pp,3)) then
                       axes=""""
                       axes2="x"
                    else if (equal_matrix(ss(1:3,1:3,i),-y_2pp,3)) then
                       axes=""""
                       axes2="y"
                    else if (equal_matrix(ss(1:3,1:3,i),-z_2pp,3)) then
                       axes=""""
                       axes2="z"
                    end if
                 end if

              case (3)
                 if (orden > 0) then
                    if (equal_matrix(ss(1:3,1:3,i),x_3,3)) then
                       axes="x"
                    else if (equal_matrix(ss(1:3,1:3,i),y_3,3)) then
                       axes="y"
                    else if (equal_matrix(ss(1:3,1:3,i),z_3,3)) then
                       axes="z"
                    else if (equal_matrix(ss(1:3,1:3,i),xyz_3,3)) then
                       axes="*"
                    end if
                 else
                    if (equal_matrix(ss(1:3,1:3,i),-x_3,3)) then
                       axes="x"
                    else if (equal_matrix(ss(1:3,1:3,i),-y_3,3)) then
                       axes="y"
                    else if (equal_matrix(ss(1:3,1:3,i),-z_3,3)) then
                       axes="z"
                    else if (equal_matrix(ss(1:3,1:3,i),-xyz_3,3)) then
                       axes="*"
                    end if
                 end if

              case (4)
                 if (orden > 0) then
                    if (equal_matrix(ss(1:3,1:3,i),x_4,3)) then
                       axes="x"
                    else if (equal_matrix(ss(1:3,1:3,i),y_4,3)) then
                       axes="y"
                    else if (equal_matrix(ss(1:3,1:3,i),z_4,3)) then
                       axes="z"
                    end if
                 else
                    if (equal_matrix(ss(1:3,1:3,i),-x_1,3)) then
                       axes="x"
                    else if (equal_matrix(ss(1:3,1:3,i),-y_4,3)) then
                       axes="y"
                    else if (equal_matrix(ss(1:3,1:3,i),-z_4,3)) then
                       axes="z"
                    end if
                 end if

              case (6)
                 if (orden > 0) then
                    if (equal_matrix(ss(1:3,1:3,i),x_6,3)) then
                       axes="x"
                    else if (equal_matrix(ss(1:3,1:3,i),y_6,3)) then
                       axes="y"
                    else if (equal_matrix(ss(1:3,1:3,i),z_6,3)) then
                       axes="z"
                    end if
                 else
                    if (equal_matrix(ss(1:3,1:3,i),-x_6,3)) then
                       axes="x"
                    else if (equal_matrix(ss(1:3,1:3,i),-y_6,3)) then
                       axes="y"
                    else if (equal_matrix(ss(1:3,1:3,i),-z_6,3)) then
                       axes="z"
                    end if
                 end if

          end select

          !---- Translations ----!
          tt=mod(nint(ts(:,i)*12.0)+48,12)

           select case (SpaceGroup%Bravais)
             case ("A")
                tt1=mod(tt-lattice(:,1)+48,12)
                if (sum(tt1) < sum(tt) ) tt=tt1

             case ("B")
                tt1=mod(tt-lattice(:,2)+48,12)
                if (sum(tt1) < sum(tt) ) tt=tt1

             case ("C")
                tt1=mod(tt-lattice(:,3)+48,12)
                if (sum(tt1) < sum(tt) ) tt=tt1

             case ("I")
                tt1=mod(tt-lattice(:,4)+48,12)
                if (sum(tt1) < sum(tt) ) tt=tt1

             case ("R")
                tt1=mod(tt-lattice(:,5)+48,12)
                tt2=mod(tt-lattice(:,6)+48,12)
                if (sum(tt1) < sum(tt) ) tt=tt1
                if (sum(tt2) < sum(tt) ) tt=tt2

             case ("F")
                tt1=mod(tt-lattice(:,3)+48,12)
                tt2=mod(tt-lattice(:,2)+48,12)
                tt3=mod(tt-lattice(:,1)+48,12)
                if (sum(tt1) < sum(tt) ) tt=tt1
                if (sum(tt2) < sum(tt) ) tt=tt2
                if (sum(tt3) < sum(tt) ) tt=tt3

          end select

          ini=1

          !---- Fractional translation ----!
          select case (abs(orden))
              case (3)
                 select case (axes)
                     case ("x")
                        if (tt(2) == 0 .and. tt(3) == 0) then
                           select case (tt(1))
                              case (4)              ! 31
                                 tras(ini:ini)="1"
                                 tt(1)=0

                              case (8)              ! 32
                                 tras(ini:ini)="2"
                                 tt(1)=0
                           end select
                        end if

                     case ("y")
                        if (tt(1) == 0 .and. tt(3) == 0) then
                           select case (tt(2))
                              case (4)              ! 31
                                 tras(ini:ini)="1"
                                 tt(2)=0

                              case (8)              ! 32
                                 tras(ini:ini)="2"
                                 tt(2)=0
                           end select
                        end if
                     case ("z")
                        if (tt(1) == 0 .and. tt(2) == 0) then
                           select case (tt(3))
                              case (4)              ! 31
                                 tras(ini:ini)="1"
                                 tt(3)=0

                              case (8)              ! 32
                                 tras(ini:ini)="2"
                                 tt(3)=0
                           end select
                        end if

                 end select

              case (6)
                 select case (axes)
                     case ("x")
                        if (tt(2) == 0 .and. tt(3) ==0) then
                           select case (tt(1))
                              case (2)              ! 61
                                 tras(ini:ini)="1"
                                 tt(1)=0

                              case (4)              ! 62
                                 tras(ini:ini)="2"
                                 tt(1)=0

                              case (8)              ! 64
                                 tras(ini:ini)="4"
                                 tt(1)=0

                              case(10)
                                 tras(ini:ini)="5"  ! 65
                                 tt(1)=0

                           end select
                        end if

                     case ("y")
                        if (tt(1) == 0 .and. tt(3) == 0) then
                           select case (tt(2))
                              case (2)              ! 61
                                 tras(ini:ini)="1"
                                 tt(2)=0

                              case (4)              ! 62
                                 tras(ini:ini)="2"
                                 tt(2)=0

                              case (8)              ! 64
                                 tras(ini:ini)="4"
                                 tt(2)=0

                              case(10)
                                 tras(ini:ini)="5"  ! 65
                                 tt(2)=0

                           end select
                        end if

                     case ("z")
                        if (tt(1) == 0 .and. tt(2) == 0) then
                           select case (tt(3))
                              case (2)              ! 61
                                 tras(ini:ini)="1"
                                 tt(3)=0

                              case (4)              ! 62
                                 tras(ini:ini)="2"
                                 tt(3)=0

                              case (8)              ! 64
                                 tras(ini:ini)="4"
                                 tt(3)=0

                              case(10)
                                 tras(ini:ini)="5"  ! 65
                                 tt(3)=0

                           end select
                        end if

                 end select
          end select

          !---- Translation vector ----!
          do k=1,3     ! 3 pasadas
             if (tt(1) == 0 .and. tt(2) == 0 .and. tt(3) == 0) exit

             do j=1,8
                tt1=tt-tras_val(:,j)
                if ( all(tt1 >= 0) ) then
                   tras(ini:ini)=l_case(traslacion(j))
                   tt(:)=tt1(:)
                   ini=ini+1
                   exit
                end if
             end do

          end do

          !---- Last check ----!
          if (nt == 1) then
             ngen=ngen+1
             generador(ngen)%orden= orden
             generador(ngen)%axes = axes
             generador(ngen)%tras = tras
          else
             if (norden(6) > 0) then
                if (abs(orden) == 6 .and. axes =="z") then
                   ngen=ngen+1

                   if (ngen > 4) then
                      err_symm=.true.
                      ERR_Symm_Mess=" Error in generators"
                      return
                   end if
                   generador(ngen)%orden= orden
                   generador(ngen)%axes = axes
                   generador(ngen)%tras = tras
                end if

                if (abs(orden) == 2 .and. (axes == "'" .or. axes =="""") .and. &
                    axes2 == "z") then
                   ngen=ngen+1

                   if (ngen > 4) then
                      err_symm=.true.
                      ERR_Symm_Mess=" Error in generators"
                      return
                   end if
                   generador(ngen)%orden= orden
                   generador(ngen)%axes = axes
                   generador(ngen)%axes2= axes2
                   generador(ngen)%tras = tras
                end if
             end if

             if (norden(4) > 0) then
                if (abs(orden) == 4 .and. axes =="z") then
                   ngen=ngen+1

                   if (ngen > 4) then
                      err_symm=.true.
                      ERR_Symm_Mess=" Error in generators"
                      return
                   end if
                   generador(ngen)%orden= orden
                   generador(ngen)%axes = axes
                   generador(ngen)%tras = tras
                end if

                if (abs(orden) == 3 .and. axes == "*") then
                   ngen=ngen+1

                   if (ngen > 4) then
                      err_symm=.true.
                      ERR_Symm_Mess=" Error in generators"
                      return
                   end if
                   generador(ngen)%orden= orden
                   generador(ngen)%axes = axes
                   generador(ngen)%tras = tras
                end if

                if (abs(orden) == 2 .and. axes == "x") then
                   ngen=ngen+1

                   if (ngen > 4) then
                      err_symm=.true.
                      ERR_Symm_Mess=" Error in generators"
                      return
                   end if
                   generador(ngen)%orden= orden
                   generador(ngen)%axes = axes
                   generador(ngen)%tras = tras
                end if
             end if

             if (norden(3) > 0 .and. norden(4) == 0) then
                if (abs(orden) == 3 .and. (axes =="z" .or. axes == "*")) then
                   ngen=ngen+1

                   if (ngen > 4) then
                      err_symm=.true.
                      ERR_Symm_Mess=" Error in generators"
                      return
                   end if
                   generador(ngen)%orden= orden
                   generador(ngen)%axes = axes
                   generador(ngen)%tras = tras
                end if

                if ( (abs(orden) == 2 .and. axes == "z")  .or. &
                     (abs(orden) == 2 .and. axes == "x")  .or. &
                     (abs(orden) == 2 .and. axes == "'")  .or. &
                     (abs(orden) == 2 .and. axes == """")) then
                   ngen=ngen+1

                   if (ngen > 4) then
                      err_symm=.true.
                      ERR_Symm_Mess=" Error in generators"
                      return
                   end if
                   generador(ngen)%orden= orden
                   generador(ngen)%axes = axes
                   generador(ngen)%axes2= axes2
                   generador(ngen)%tras = tras
                end if
             end if

             if (norden(2) > 0 .and. norden(3) == 0 .and. norden(4) == 0  &
                .and. norden(6) == 0) then
                if (abs(orden) == 2 .and. axes == "z") then
                   ngen=ngen+1

                   if (ngen > 4) then
                      err_symm=.true.
                      ERR_Symm_Mess=" Error in generators"
                      return
                   end if
                   generador(ngen)%orden= orden
                   generador(ngen)%axes = axes
                   generador(ngen)%tras = tras
                end if
                if (abs(orden) == 2 .and. axes == "x") then
                   ngen=ngen+1

                   if (ngen > 4) then
                      err_symm=.true.
                      ERR_Symm_Mess=" Error in generators"
                      return
                   end if
                   generador(ngen)%orden= orden
                   generador(ngen)%axes = axes
                   generador(ngen)%tras = tras
                end if
             end if

          end if
       end do

       !---- Purge Generators ----!
       j=0
       k=0
       if (ngen > 1) then
          do i=1,ngen
             if (generador(i)%axes =="'") j=i
             if (generador(i)%axes =="""") k=i
          end do
          if (j /= 0 .and. k /=0) then
             if (generador(j)%axes2 =="z") then
                do i=k+1,ngen
                   generador(i-1)=generador(i)
                end do
                ngen=ngen-1
             else
                do i=j+1,ngen
                   generador(i-1)=generador(i)
                end do
                ngen=ngen-1
             end if
          end if
       end if

       !---- Order Generators ----!
       select case (ngen)
          case (2)
             if (abs(generador(1)%orden) < abs(generador(2)%orden)) then
                generador(5)=generador(1)
                generador(1)=generador(2)
                generador(2)=generador(5)
             else if (abs(generador(1)%orden) == abs(generador(2)%orden)) then
                if (generador(2)%axes == "z") then
                   generador(5)=generador(1)
                   generador(1)=generador(2)
                   generador(2)=generador(5)
                end if
             end if
          case (3)
             do i=1,3
                if (abs(generador(i)%orden) == 1 .or. abs(generador(i)%orden) == 3) then
                   generador(5)=generador(i)
                   do j=i+1,3
                      generador(j-1)=generador(j)
                   end do
                   generador(3)=generador(5)
                end if
             end do
             if (abs(generador(1)%orden) < abs(generador(2)%orden)) then
                generador(5)=generador(1)
                generador(1)=generador(2)
                generador(2)=generador(5)
             else if (abs(generador(1)%orden) == abs(generador(2)%orden)) then
                if (generador(2)%axes == "z") then
                   generador(5)=generador(1)
                   generador(1)=generador(2)
                   generador(2)=generador(5)
                end if
             end if

          case (4)
             do i=1,4
                if (abs(generador(i)%orden) == 1) then
                   generador(5)=generador(i)
                   do j=i+1,4
                      generador(j-1)=generador(j)
                   end do
                   generador(4)=generador(5)
                end if
             end do
             do i=1,3
                if (abs(generador(i)%orden) == 3) then
                   generador(5)=generador(i)
                   do j=i+1,3
                      generador(j-1)=generador(j)
                   end do
                   generador(3)=generador(5)
                end if
             end do
             if (abs(generador(1)%orden) < abs(generador(2)%orden)) then
                generador(5)=generador(1)
                generador(1)=generador(2)
                generador(2)=generador(5)
             else if (abs(generador(1)%orden) == abs(generador(2)%orden)) then
                if (generador(2)%axes == "z") then
                   generador(5)=generador(1)
                   generador(1)=generador(2)
                   generador(2)=generador(5)
                end if
             end if

       end select

       !---- Hall Symbol ----!
       ini=len_trim(hall)
       ini=ini+1

       do i=1,ngen
          !---- Rotation ----!
          if (generador(i)%orden >0) then
             ini=ini+1
             write(unit=hall(ini:ini),fmt="(i1)") generador(i)%orden
          else
             ini=ini+1
             write(unit=hall(ini:ini+1),fmt="(i2)") generador(i)%orden
             ini=ini+1
          end if

          !---- Axis ----!
          select case (i)
             case (1)
                if (generador(i)%axes /= "z") then
                   ini=ini+1
                   hall(ini:ini)=generador(i)%axes
                end if

             case (2)
                if (abs(generador(i)%orden) == 2) then
                   if (abs(generador(1)%orden) == 2 .or. abs(generador(1)%orden) == 4) then
                      if (generador(i)%axes /= "x") then
                         ini=ini+1
                         hall(ini:ini)=generador(i)%axes
                      end if
                   else if (abs(generador(1)%orden) == 3 .or. abs(generador(1)%orden) == 6) then
                      if (generador(i)%axes /= "'") then
                         ini=ini+1
                         hall(ini:ini)=generador(i)%axes
                      end if
                   end if

                else
                   if (abs(generador(i)%orden) /= 1) then
                      ini=ini+1
                      hall(ini:ini)=generador(i)%axes
                   end if
                end if

             case (3)
                if (abs(generador(i)%orden) /= 3 .and. abs(generador(i)%orden) /= 1) then
                   ini=ini+1
                   hall(ini:ini)=generador(i)%axes
                end if

             case (4)

          end select

          !---- Translation ----!
          select case (len_trim(generador(i)%tras))
             case (1)
                ini=ini+1
                hall(ini:ini)=generador(i)%tras

             case (2)
                ini=ini+1
                hall(ini:ini+1)=generador(i)%tras
                ini=ini+1

             case (3)
                ini=ini+1
                hall(ini:ini+2)=generador(i)%tras
                ini=ini+2

          end select
          ini=ini+1
       end do

       !---- Check the Hall Symbol for repetitions of minus sign ----!

       i=index(hall,"-")
       if(i /= 0 ) then
         k=index(hall,"-",back=.true.)
         if(k /= i ) then
           hall=hall(1:k-1)//hall(k+1:)
         end if
       end if

       !---- Is the Hall Symbol in the table? ----!
       k=0
       do i=1,num_spgr_info
          if (hall(1:16) == spgr_info(i)%hall) then
             k=i
             exit
          end if
       end do

       if(hall(1:1) /= "-") hall=" "//hall
       Spacegroup%Hall=hall

       if (k /= 0) then
          SpaceGroup%NumSpg       = spgr_info(k)%n
          SpaceGroup%Spg_Symb     = spgr_info(k)%hm
                call get_laue_str(spgr_info(k)%laue,SpaceGroup%Laue)
                call get_PointGroup_str(spgr_info(k)%pg,SpaceGroup%PG)
          SpaceGroup%Info         = spgr_info(k)%inf_extra
          SpaceGroup%R_Asym_Unit(1,1) = real(spgr_info(k)%asu(1))/24.0
          SpaceGroup%R_Asym_Unit(2,1) = real(spgr_info(k)%asu(2))/24.0
          SpaceGroup%R_Asym_Unit(3,1) = real(spgr_info(k)%asu(3))/24.0
          SpaceGroup%R_Asym_Unit(1,2) = real(spgr_info(k)%asu(4))/24.0
          SpaceGroup%R_Asym_Unit(2,2) = real(spgr_info(k)%asu(5))/24.0
          SpaceGroup%R_Asym_Unit(3,2) = real(spgr_info(k)%asu(6))/24.0
       else
          SpaceGroup%Spg_Symb     = "Unknown"
          SpaceGroup%Info         = "User-provided generators "

       end if

       if (present(SpaceH) ) SpaceH=hall

       return
    End Subroutine Get_HallSymb_from_Gener

    !!----
    !!---- Subroutine Get_Lattice_Type(L, Latc, Lattyp)
    !!----    integer,                        intent(in)  :: L         !  number of centring vectors
    !!----    real(kind=cp), dimension(:,:),  intent(in)  :: Latc      ! (3,11) centring vectors
    !!----    character(len=*),               intent(out) :: lattyp    ! Lattice symbol
    !!----
    !!----    Subroutine to get the lattice symbol from a set of centring vectors.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Get_Lattice_Type(L, Latc, lattyp)
       !---- Arguments ----!
       integer,                        intent( in) :: L
       real(kind=cp), dimension(:,:),  intent( in) :: Latc
       character(len=*),               intent(out) :: lattyp

       !---- Local variables ----!
       logical :: latt_p, latt_a, latt_b, latt_c, latt_i, latt_r, latt_f, latt_z
       integer, dimension(6) :: latt_given
       integer, dimension(3) :: tt
       integer               :: i, j
       integer, dimension(3,6), parameter :: lattice=reshape((/0,6,6, 6,0,6, &
                                                     6,6,0, 6,6,6, 8,4,4, 4,8,8/),(/3,6/))

       if (l > 3) then  !non conventional centring
          lattyp="Z"
          return
       else if(l == 0) then !primitive lattice
          lattyp="P"
          return
       end if

       latt_p=.true.
       latt_a=.false.
       latt_b=.false.
       latt_c=.false.
       latt_i=.false.
       latt_r=.false.
       latt_f=.false.
       latt_z=.false.

       do i=1,L
          tt(1:3)=nint(12.0 * Latc(1:3,i))   ! Translations x 12

          !---- Compare the translation part of the operator with tabulated array ----!
          latt_given(:) = 0
          do j=1,6
             if (equal_vector(tt,lattice(:,j),3)) then
                latt_given(j) = 1
                select case (j)
                   case (1)
                      latt_a=.true.
                   case (2)
                      latt_b=.true.
                   case (3)
                      latt_c=.true.
                   case (4)
                      latt_i=.true.
                   case (5,6)
                      latt_r=.true.
                end select
                exit
             end if
          end do
          if (sum(latt_given) == 0) then
             latt_z = .true.
             exit
          end if
       end do

       !---- Lattice Type ----!
       if (latt_z) then
           lattyp="Z"
           return
       end if
       if ( (latt_a .and. latt_b .and. latt_c) .or. (latt_a .and. latt_b) .or. &
            (latt_a .and. latt_c) .or. (latt_b .and. latt_c) ) then
            latt_f=.true.
            latt_a=.false.
            latt_b=.false.
            latt_c=.false.
            latt_p=.false.
            latt_i=.false.
       end if
       if (latt_p) lattyp="P"
       if (latt_a) lattyp="A"
       if (latt_b) lattyp="B"
       if (latt_c) lattyp="C"
       if (latt_i) lattyp="I"
       if (latt_r) lattyp="R"
       if (latt_f) lattyp="F"

       return
    End Subroutine Get_Lattice_Type

    !!----
    !!---- Subroutine Get_Laue_Pg(Spacegroup, Laue_Car, Point_Car)
    !!----    type (Space_Group_Type),  intent( in) :: SpaceGroup   !  In -> Space Group type variable
    !!----    character(len=*),         intent(out) :: Laue_car     ! Out -> String with Laue symbol
    !!----    character(len=*),         intent(out) :: Point_car    ! Out -> String with Point Group symbol
    !!----
    !!----    Subroutine to get the information of Laue and Point Group.
    !!----    Vvalid only for conventional bases for Point Group
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Get_Laue_PG(SpaceGroup, Laue_car, Point_car)
       !---- Arguments ----!
       type (Space_Group_Type),  intent( in) :: SpaceGroup
       character (len=*),        intent(out) :: Laue_car
       character (len=*),        intent(out) :: Point_car

       !---- Local variables ----!
       integer :: nrot_1, nrot_1b
       integer :: nrot_2, nrot_2b
       integer :: nrot_3, nrot_3b
       integer :: nrot_4, nrot_4b
       integer :: nrot_6, nrot_6b
       integer :: i,n_m,ndet,ind

       !---- Initializing ----!
       point_car=" "
       laue_car=" "

       nrot_1  = 0
       nrot_2  = 0
       nrot_3  = 0
       nrot_4  = 0
       nrot_6  = 0
       nrot_1b = 0
       nrot_2b = 0
       nrot_3b = 0
       nrot_4b = 0
       nrot_6b = 0
       n_m = 0

       call init_err_symm()
       if (spacegroup%numops == 0) then
          err_symm=.true.
          ERR_Symm_Mess=" No symmetry operators are given"
          return
       end if
       do i=1,spacegroup%numops
          ndet= Axes_Rotation(SpaceGroup%Symop(i)%Rot(:,:))
          select case (ndet)
             case (-6)
                nrot_6b=nrot_6b +1
             case (-4)
                nrot_4b=nrot_4b +1
             case (-3)
                nrot_3b=nrot_3b +1
             case (-2)
                nrot_2b=nrot_2b +1
             case (-1)
                nrot_1b=nrot_1b +1
             case ( 1)
                nrot_1 =nrot_1  +1
             case ( 2)
                nrot_2 =nrot_2  +1
             case ( 3)
                nrot_3 =nrot_3  +1
             case ( 4)
                nrot_4 =nrot_4  +1
             case ( 6)
                nrot_6 =nrot_6  +1
             case default
                err_symm=.true.
                ERR_Symm_Mess=" Rotation Not Determined"
                return
          end select
       end do

       n_m = nrot_1  + nrot_2  + nrot_3  + nrot_4  + nrot_6  + &
             nrot_1b + nrot_2b + nrot_3b + nrot_4b + nrot_6b

       !---- Cubic ----!
       if ( (nrot_3 + nrot_3b == 8) ) then
          select case (n_m)
             case (12)
                if (SpaceGroup%Centred ==1) then
                   point_car="23"
                else
                   point_car="m-3"
                end if
                laue_car="m-3"

             case (24)
                if (SpaceGroup%Centred /=1) then
                   point_car="m-3m"
                else
                   if (nrot_4  == 6) point_car="432"
                   if (nrot_4b == 6) point_car="-43m"
                end if
                laue_car="m-3m"
          end select

       !---- Hexagonal ----!
       else if ( (nrot_6 + nrot_6b == 2) ) then
          select case (n_m)
             case (6)
                if (SpaceGroup%Centred /=1) then
                   point_car="6/m"
                else
                   if (nrot_6  == 2) point_car="6"
                   if (nrot_6b == 2) point_car="-6"
                end if
                laue_car="6/m"

             case (12)
                if (SpaceGroup%Centred /=1) then
                   point_car="6/mmm"
                else
                   if (nrot_6b == 2) then
                      do i=1,spacegroup%numops
                         ndet= Axes_Rotation(SpaceGroup%Symop(i)%Rot(:,:))
                         if (ndet /= 2) cycle
                         !---- This is only valid for conventional bases ---!
                         call SearchOp(SpaceGroup%Symop(i)%Rot(:,:),25,36,ind)
                         if (ind < 0) then
                            ind=-ind-12
                         end if
                         select case (ind)
                            case (31)
                               point_car="-62m"
                            case default
                               point_car="-6m2"
                         end select
                         exit
                      end do
                   end if
                   if ( (nrot_6  == 2 .and. nrot_2 == 7) ) point_car="622"
                   if ( (nrot_6  == 2 .and. nrot_2b== 6) ) point_car="6mm"
                end if
                laue_car="6/mmm"
          end select

       !---- Trigonal ----!
       else if ( (nrot_3 + nrot_3b == 2) ) then
          select case (n_m)
             case (3)
                if (SpaceGroup%Centred /=1) then
                   point_car="-3"
                else
                   point_car="3"
                end if
                laue_car="-3"

             case (6)
                if (SpaceGroup%Hexa) then
                   if (SpaceGroup%Centred /=1) then
                      do i=1,spacegroup%numops
                           ndet=Axes_Rotation(SpaceGroup%Symop(i)%Rot(:,:))
                           if (ndet /= -2) cycle
                         !---- This is only valid for conventional bases ---!
                         call SearchOp(SpaceGroup%Symop(i)%Rot(:,:),25,36,ind)
                         if (ind < 0) then
                            ind=-ind-12
                         end if
                         select case (ind)
                            case (22)
                               point_car="-31m"
                               laue_car ="-31m"
                            case default
                               point_car="-3m"
                               laue_car ="-3m"
                         end select
                         exit
                      end do
                   else
                      if (nrot_2  == 3 ) then
                         do i=1,spacegroup%numops
                            ndet=Axes_Rotation(SpaceGroup%Symop(i)%Rot(:,:))
                            if (ndet /= 2) cycle
                            !---- This is only valid for conventional bases ---!
                            call SearchOp(SpaceGroup%Symop(i)%Rot(:,:),25,36,ind)
                            if (ind < 0) then
                               ind=-ind-12
                            end if
                            select case (ind)
                               case (34)
                                  point_car="-312"
                                  laue_car ="-31m"
                               case default
                                  point_car="-32"
                                  laue_car ="-3m"
                            end select
                            exit
                         end do
                      end if

                      if (nrot_2b == 3 ) then
                         do i=1,spacegroup%numops
                            ndet=Axes_Rotation(SpaceGroup%Symop(i)%Rot(:,:))
                            if (ndet /= -2) cycle
                            !---- This is only valid for conventional bases ---!
                            call SearchOp(SpaceGroup%Symop(i)%Rot(:,:),25,36,ind)
                            if (ind < 0) then
                               ind=-ind-12
                            end if
                            select case (ind)
                               case (22)
                                  point_car="31m"
                                  laue_car ="-31m"
                               case default
                                  point_car="3m"
                                  laue_car ="-3m"
                            end select
                            exit
                         end do
                      end if
                   end if
                else
                   if (SpaceGroup%Centred /=1) then
                      point_car="-3m"
                   else
                      if (nrot_2  == 3 ) point_car="32"
                      if (nrot_2b == 3 ) point_car="3m"
                   end if
                   laue_car="-3m"
                end if

          end select

       !---- Tetragonal ----!
       else if ( (nrot_4 + nrot_4b == 2) ) then
          select case (n_m)
             case (4)
                if (SpaceGroup%Centred /=1) then
                   point_car="4/m"
                else
                   if (nrot_4  == 2 ) point_car="4"
                   if (nrot_4b == 2 ) point_car="-4"
                end if
                laue_car="4/m"

             case (8)
                if (SpaceGroup%Centred /=1) then
                   point_car="4/mmm"
                else
                   if (nrot_4b == 2 ) then
                   do i=1,spacegroup%numops
                         ndet=Axes_Rotation(SpaceGroup%Symop(i)%Rot(:,:))
                      if (ndet /= -2) cycle
                      !---- This is only valid for conventional bases ---!
                      call SearchOp(SpaceGroup%Symop(i)%Rot(:,:),1,24,ind)
                      if (ind < 0) then
                         ind=24-ind
                      end if
                      select case (ind)
                         case (28)
                            point_car="-4m2"
                         case default
                            point_car="-42m"
                      end select
                      exit
                   end do
                end if
                if ( (nrot_4  == 2 .and. nrot_2 == 5) ) point_car="422"
                if ( (nrot_4  == 2 .and. nrot_2b== 4) ) point_car="4mm"
             end if
             laue_car="4/mmm"

          end select

       !---- Orthorhombic ----!
       else if ( (nrot_2 + nrot_2b == 3) ) then
          if (SpaceGroup%Centred /=1) then
             point_car="mmm"
          else
             if (nrot_2  == 3 ) point_car="222"
             if (nrot_2b == 2 ) then
                do i=1,spacegroup%numops
                       ndet=Axes_Rotation(SpaceGroup%Symop(i)%Rot(:,:))
                   if (ndet /= 2) cycle
                   !---- This is only valid for conventional bases ---!
                   call SearchOp(SpaceGroup%Symop(i)%Rot(:,:),1,24,ind)
                   select case (ind)
                      case (4)
                         point_car="2mm"
                      case (3)
                         point_car="m2m"
                      case default
                         point_car="mm2"
                   end select
                   exit
                end do
             end if
          end if
          laue_car="mmm"

       !---- Monoclinic ----!
       else if ( (nrot_2 + nrot_2b == 1)  ) then
          if (SpaceGroup%Centred /=1) then
             point_car="2/m"
          else
             if (nrot_2  == 1 ) point_car="2"
             if (nrot_2b == 1 ) point_car="m"
          end if
          laue_car="2/m"

       !---- Triclinic ----!
       else if (n_m == 1) then
          if (SpaceGroup%Centred /=1) then
             point_car="-1"
          else
             point_car="1"
          end if
          laue_car="-1"

       end if

       return
    End Subroutine Get_Laue_PG

    !!----
    !!---- Subroutine Get_Laue_Str(Ilaue,Laue_Str)
    !!----    integer,          intent( in) :: ilaue         !  In -> Ordinal number in LAUE_CLASS
    !!----    character(len=*), intent(out) :: Laue_Str      ! Out -> String with the Laue class
    !!----
    !!----    Obtain the string for the Laue-Class. Control of error is
    !!----    present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Get_Laue_Str(Ilaue,Str)
       !---- Arguments ----!
       integer,          intent( in) :: ilaue
       character(len=*), intent(out) :: str

       call init_err_symm()
       if (ilaue < 1 .or. ilaue > 16) then
          err_symm=.true.
          ERR_Symm_Mess=" Laue Number Incorrect"
       else
          str=laue_class(ilaue)
       end if

       return
    End Subroutine Get_Laue_Str

    !!----
    !!----  Subroutine Get_Orbit(X,Spg,Mult,orb,ptr,prim,symm)
    !!----    real(kind=cp), dimension(3),  intent (in) :: x     !  In -> Position vector
    !!----    type(Space_Group_type),       intent (in) :: spgr  !  In -> Space Group
    !!----    integer,                      intent(out) :: mult  !  Out -> Multiplicity
    !!----    real(kind=cp), dimension(:,:),intent(out) :: orb   !  Out -> List of equivalent positions
    !!----    integer,dimension(:),optional,intent(out) :: ptr   !  Out -> Pointer to the symmetry elements
    !!----    character(len=*),    optional,intent( in) :: prim  !  In  -> If given, only the primitive cell is considered
    !!----    character(len=*),    optional,intent( in) :: symm  !  In  -> If given, the coordinates are normalized as to be -1/2 <= x <1/2
    !!----
    !!----    Obtain the multiplicity and list of equivalent positions
    !!----    (including centring!) modulo integer lattice translations or within the range [-1/2,1/2) if symm is given.
    !!----
    !!---- Update: June - 2011 (JRC - removing pointer to stabilizer)
    !!
    Subroutine Get_Orbit(x,Spg,Mult,orb,ptr,prim,symm)
       !---- Arguments ----!
       real(kind=cp), dimension(3),  intent (in) :: x
       type(Space_Group_type),       intent (in) :: spg
       integer,                      intent(out) :: mult
       real(kind=cp),dimension(:,:), intent(out) :: orb
       integer,dimension(:),optional,intent(out) :: ptr
       character(len=*),    optional,intent( in) :: prim
       character(len=*),    optional,intent( in) :: symm

       !---- Local variables ----!
       integer                                :: j, nt
       real(kind=cp), dimension(3)            :: xx,v
       character(len=1)                       :: laty

       laty="P"
       if(present(prim)) laty=Spg%spg_lat
       mult=1
       orb(:,1)=x(:)
       if(present(ptr)) ptr(mult) = 1
       ext: do j=2,Spg%multip
          xx=ApplySO(Spg%SymOp(j),x)
          xx=modulo_lat(xx)
          do nt=1,mult
             v=orb(:,nt)-xx(:)
             if (Lattice_trans(v,Spg%spg_lat)) then
               if (.not. Lattice_trans(v,laty)) cycle  !Count in orbit the centred related atoms
               cycle ext
             end if
          end do
          mult=mult+1
          orb(:,mult)=xx(:)
          if(present(ptr)) ptr(mult) = j   !Effective symop
       end do ext

       if(present(symm)) then
         !Normalize the coordinates to be -1/2 <= x < 1/2
         do j=1,Mult
           do nt=1,3
              if(Orb(nt,j) >= 0.5) Orb(nt,j)= Orb(nt,j) - 1.0
           end do
         end do
       end if

       return
    End Subroutine Get_Orbit

    !!----
    !!---- Subroutine Get_Pointgroup_Str(Ipg,Str)
    !!----    integer,          intent( in) :: ipg        !  In -> Ordinal number for POINT_GROUP
    !!----    character(len=*), intent(out) :: Str        ! Out -> String for Point Group
    !!----
    !!----    Obtain the string for the Point Group. Error control is present
    !!----
    !!---- Update: Update: July - 2014: added m3 and m3m for compatibility with Laue_class
    !!
    Subroutine Get_Pointgroup_Str(Ipg,Str)
       !---- Arguments ----!
       integer,          intent( in) :: ipg
       character(len=*), intent(out) :: str

       call init_err_symm()
       if (ipg < 1 .or. ipg > 41) then
          err_symm=.true.
          ERR_Symm_Mess=" Point Group Number Incorrect"
       else
          str=point_group(ipg)
       end if

       return
    End Subroutine Get_PointGroup_Str

    !!--++
    !!--++ Subroutine Get_Seitz(N_Op,Tt,Seitz_Symb)
    !!--++    integer,                     intent( in) :: n_op          !  In -> Number of the rotational matrix
    !!--++    real(kind=cp), dimension(3), intent( in) :: tt            !  In -> Translation part
    !!--++    character (len=*),           intent(out) :: Seitz_symb    ! Out -> Seitz Symbol
    !!--++
    !!--++    (PRIVATE)
    !!--++    Provide the Seitz symbol of a symmetry operator.
    !!--++    This is mainly for internal use in the module.
    !!--++    Run before SearchOp.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Get_Seitz(n_op,tt,Seitz_symb)
       !---- Arguments ----!
       integer,                     intent( in) :: n_op
       real(kind=cp), dimension(3), intent( in) :: tt
       character (len=*),           intent(out) :: Seitz_symb

       !---- Local variables ----!
       character (len=*), dimension(16), parameter  :: fracc =(/" 0 ","1/2","1/3","2/3",    &
                        "1/4","3/4","1/6","5/6","1/8","3/8","5/8","7/8","1  ","2  ","3  ","4  "/)
       integer :: i,j,ini
       real(kind=cp), dimension(16), parameter :: frac= (/0.0, 0.5,1.0/3.0,2.0/3.0,0.25,0.75, &
                                                  1.0/6.0,5.0/6.0,0.125,0.375,0.625,0.875,1.0,2.0,3.0,4.0/)

       if (hexa) then
          Seitz_symb(1:14) ="{"//X_d6h(n_op)(2:13)//"|"
          ini=15
       else
          Seitz_symb(1:10) ="{"//X_Oh(n_op)(2:9)//"|"
          ini=11
       end if
       xyz:do i=1,3
          do j=1,16
             if (abs(frac(j)-abs(tt(i))) < eps_symm) then
                if (tt(i) < 0.0) then
                   Seitz_symb(ini:ini+3)="-"//fracc(j)
                else
                   Seitz_symb(ini:ini+3)=" "//fracc(j)
                end if
                ini=ini+4
                cycle xyz
             end if
          end do
       end do xyz

       Seitz_symb(ini:ini)="}"

       return
    End Subroutine Get_Seitz

    !!----
    !!---- Subroutine Get_Seitz_Symbol(iop,itim,tr,Seitz_symb)
    !!----    integer,                   intent(in) :: iop,itim      !  In -> Number of the rotational matrix, time inversion
    !!----    real(kind=cp),dimension(3),intent(in) :: tr            !  In -> Translation part
    !!----    character(len=*),          intent(out):: Seitz_symb    ! Out -> Seitz Symbol
    !!----
    !!----    Provide the Seitz symbol of a symmetry operator. It uses the Litvin notation and
    !!----    the ordering is that of Table given by Harold T. Stokes and Branton J. Campbell.
    !!----    Hexa should be defined before using this subroutine. This subroutine is intended
    !!----    to be used with the reading of Magnetic Space Groups (see CFML_Magnetic_Symmetry)
    !!----
    !!---- Update: November 2012
    !!
    Subroutine Get_Seitz_Symbol(iop,itim,tr,Seitz_symb)
      integer,                     intent(in) :: iop,itim
      real(kind=cp), dimension(3), intent(in) :: tr
      character(len=*),            intent(out):: Seitz_symb
      !---- Local variables ----!
      integer :: i
      character(len=25) :: transl
      character(len=8)  :: operator_symb
      character(len=6)  :: Fracc

      if(hexa) then
        Operator_symb=Litvin_point_op_hex_label(iop)
      else
        Operator_symb=Litvin_point_op_label(iop)
      end if
      transl=" "
      do i=1,3
        call Get_Fraction_2Dig(tr(i),Fracc)
        transl=trim(transl)//trim(Fracc)//","
      end do
      i=len_trim(transl)
      transl(i:i)=" "
      do i=1,len_trim(transl)
        if(transl(i:i) == "+") transl(i:i)=" "
      end do
      Seitz_symb="("//trim(operator_symb)//" | "//trim(transl)//")"
      Seitz_symb=Pack_String(Seitz_symb)
      if(itim == -1)  Seitz_symb=trim(Seitz_symb)//"'"
      return
    End Subroutine Get_Seitz_Symbol


    !!--++
    !!--++ Subroutine Get_Setting_Info(Mat,orig,setting,matkind)
    !!--++    real(kind=cp), dimension (3,3),intent( in)    :: Mat     ! Matrix transforming the basis
    !!--++    real(kind=cp), dimension (  3),intent( in)    :: orig    ! Coordinates of the new origin
    !!--++    character (len=*),             intent(out)    :: setting ! String with the new setting
    !!--++    character (len=*), optional,   intent( in)    :: matkind ! Type of the input matrix
    !!--++
    !!--++    (PRIVATE)
    !!--++    Subroutine to construct a string with the transformation of the basis
    !!--++    corresponding to the matrix "Mat" and new origin "orig"
    !!--++    If matkind is given and matkind="it"/"IT", the input matrix is given
    !!--++    as in International Tables: (a' b' c') = (a b c) Mat
    !!--++    If matkind is not given or if it is not equal to "it"/"IT" the input matrix
    !!--++    is the transpose of the International convention (column matrices for basis vectors)
    !!--++    An example of the output is: a'=a+c, b'=2b, c'=-a+c  -> Origin: (0,1/4,0)
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Get_Setting_Info(Mat,orig,setting,matkind)
       !---- Arguments ----!
       real(kind=cp), dimension (3,3),intent( in)    :: Mat
       real(kind=cp), dimension (  3),intent( in)    :: orig
       character (len=*),             intent(out)    :: setting
       character (len=*), optional,   intent( in)    :: matkind

       !---- local variables ----!
       real(kind=cp), dimension (  3), parameter  :: nul = (/ 0.0, 0.0, 0.0/)
       real(kind=cp), dimension (3,3)  :: S
       character (len=22)     :: tro
       integer                :: i

       if (present(matkind)) then
          if (matkind(1:2) == "it" .or. matkind(1:2) == "IT" ) then
             S=transpose(Mat)
          else
             S=Mat
          end if
       else
          S=Mat
       end if

       call Get_SymSymb(S,nul,setting)
       i=index(setting,",")
       setting="a'="//setting(1:i)//" b'="//setting(i+1:)
       i=index(setting,",",back=.true.)
       setting=setting(1:i)//" c'="//setting(i+1:)
       do i=1,len_trim(setting)
          if (setting(i:i) == "x")  setting(i:i) = "a"
          if (setting(i:i) == "y")  setting(i:i) = "b"
          if (setting(i:i) == "z")  setting(i:i) = "c"
       end do

       call Frac_Trans_2Dig(Orig,tro)
       i=len_trim(setting)
       setting(i+2:)=" -> Origin: "//trim(tro)

       return
    End Subroutine Get_Setting_Info

    !!----
    !!---- Subroutine Get_Shubnikov_Operator_Symbol(Mat,Rot,tr,ShOp_symb,mcif)
    !!----   integer, dimension(3,3), intent(in) :: Mat,Rot     ! Symmetry operators for positions and magnetic moments
    !!----   real,    dimension(3),   intent(in) :: tr          ! Translation associated to the symmetry operator
    !!----   character(len=*),        intent(out):: ShOp_symb   ! String with the Shubnikov operator symbol
    !!----   logical,  optional,      intent(in) :: mcif        ! if present the Shubnikov operator is like in mcif: -x,y+1/2,z  mx,-my,-mz +1
    !!----
    !!---- Subroutine to construct a string with the Shubnikov operator
    !!---- in the following form: (-x,y+1/2,-z;u,-v,w)
    !!---- It also working for Wyckoff positions, when the matrices Mat and Rot
    !!---- are not symmetry operators (det=0). It is extensively used when reading
    !!---- the database containing the Magnetic Space Groups provided by
    !!---- < Harold T. Stokes and Branton J. Campbell
    !!----   Brigham Young University, Provo, Utah, USA
    !!----   June 2010 >
    !!----
    !!---- Updated: November 2012, January 2014
    !!----
    Subroutine Get_Shubnikov_Operator_Symbol(Mat,Rot,tr,ShOp_symb,mcif)
      integer,       dimension(3,3), intent(in) :: Mat,Rot
      real(kind=cp), dimension(3),   intent(in) :: tr
      character(len=*),              intent(out):: ShOp_symb
      logical, optional,             intent(in) :: mcif
      !---- Local variables ----!
      integer                 :: i,i1,i2,idet
      integer, dimension(3,3) :: sMat
      character(len=25)       :: xyz_op, uvw_op, mxmymz_op
      character(len=2)        :: time_inv

      call Get_SymSymb(Mat,tr,xyz_op)
      call Get_SymSymb(Rot,(/0.0_cp,0.0_cp,0.0_cp/),uvw_op)

      do i=1,len_trim(uvw_op)
        if(uvw_op(i:i) == "x")  uvw_op(i:i)="u"
        if(uvw_op(i:i) == "y")  uvw_op(i:i)="v"
        if(uvw_op(i:i) == "z")  uvw_op(i:i)="w"
      end do
      i1=index(xyz_op,",")
      if(i1 == 1) xyz_op="0"//trim(xyz_op)
      i2=index(xyz_op,",",back=.true.)
      if(i2 == len_trim(xyz_op)) xyz_op=trim(xyz_op)//"0"
      i1=index(xyz_op,",,")
      if(i1 /= 0) xyz_op=xyz_op(1:i1)//"0"//xyz_op(i1+1:)

      i1=index(uvw_op,",")
      if(i1 == 1) uvw_op="0"//trim(uvw_op)
      i2=index(uvw_op,",",back=.true.)
      if(i2 == len_trim(uvw_op)) uvw_op=trim(uvw_op)//"0"
      i1=index(uvw_op,",,")
      if(i1 /= 0) uvw_op=uvw_op(1:i1)//"0"//uvw_op(i1+1:)
      xyz_op=Pack_string(xyz_op)
      uvw_op=Pack_string(uvw_op)
      if(present(mcif)) then
        idet=determ_A(Mat)
        sMat=(idet*Mat-Rot)
        if(sum(sMat) == 0) then
          time_inv="+1"
        else
          time_inv="-1"
        end if
        !Expand the operator uvw_op to convert it to mx,my,mz like
        mxmymz_op=" "
        do i=1,len_trim(uvw_op)
          Select Case(uvw_op(i:i))
            case("u")
               mxmymz_op=trim(mxmymz_op)//"mx"
            case("v")
               mxmymz_op=trim(mxmymz_op)//"my"
            case("w")
               mxmymz_op=trim(mxmymz_op)//"mz"
            case default
               mxmymz_op=trim(mxmymz_op)//uvw_op(i:i)
          End Select
        end do
        ShOp_symb=trim(xyz_op)//" "//trim(mxmymz_op)//" "//time_inv
      else
        ShOp_symb="("//trim(xyz_op)//";"//trim(uvw_op)//")"
      end if
      return
    End Subroutine Get_Shubnikov_Operator_Symbol

    !!----
    !!---- Subroutine Get_So_From_Fix(Isystm,Isymce,Ibravl,Ng,Ss,Ts,Latsy,Co,Spacegen,lsym)
    !!----    integer,                     intent(out) :: ISYSTM    ! Out -> Number of the crystalline system
    !!----                                                          ! Out    (1:T, 2:M, 3:O, 4:T, 5:R-Trg, 6:H, 7:C)
    !!----    integer,                     intent(out) :: ISYMCE    ! Out -> 0 Centric (-1 not at origin)
    !!----                                                                   1 Acentric
    !!----                                                                   2 Centric (-1 at origin)
    !!----    integer,                     intent(out) :: IBRAVL    ! Out -> Index of the Bravais Lattice type
    !!----                                                                   1   2   3   4   5   6   7   8
    !!----                                                                   "P","A","B","C","I","R","F","Z"
    !!----    integer,                     intent(in ) :: NG        !  In -> Number of symmetry operators
    !!----    real(kind=cp),dimension(:,:),intent(in ) :: TS        !  In -> Translation parts of the symmetry operators(3,48)
    !!----    integer, dimension(:,:,:),   intent(in ) :: SS        !  In -> Rotation parts of the symmetry operators (3,3,48)
    !!----    character (len=2),           intent(out) :: latsy     ! Out -> Bravais Lattice symbol
    !!----    real(kind=cp),dimension(3)  ,intent(out) :: Co        ! Out -> Coordinates of origin
    !!----    character (len=1),           intent(out) :: SpaceGen  ! Out -> Type of Cell
    !!----    character (len=1),           intent(in)  :: lsym      ! In  -> Type of Cell forced
    !!----
    !!----    Determines some of items of the object Space_Group_Type from FIXed
    !!----    symmetry operators given by user.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Get_SO_from_FIX(Isystm,Isymce,Ibravl,Ng,Ss,Ts,Latsy,Co,SpaceGen,lsym)
       !---- Arguments ----!
       integer,                      intent(out) :: Isystm
       integer,                      intent(out) :: Isymce
       integer,                      intent(out) :: Ibravl
       integer,                      intent(in ) :: Ng
       integer, dimension(:,:,:),    intent(in ) :: Ss  !(3,3,48)
       real(kind=cp),dimension(:,:), intent(in ) :: Ts  !(3  ,48)
       character (len= 2),           intent(out) :: Latsy
       real(kind=cp),dimension(3),   intent(out) :: Co
       character (len= 1),           intent(out) :: SpaceGen
       character (len= 1),optional,  intent(in ) :: lsym

       !---- Local Variables ----!
       logical :: latt_p, latt_a, latt_b, latt_c, latt_i, latt_r, latt_f, latt_z

       character(len=*), dimension(8),  parameter :: red = &
                         (/"P","A","B","C","I","R","F","Z"/)

       integer, dimension(3,6), parameter :: lattice=reshape((/0,6,6, 6,0,6, &
                                                   6,6,0, 6,6,6, 8,4,4, 4,8,8/),(/3,6/))

       integer, dimension(4,4), parameter :: identidad = reshape((/1, 0, 0, 0, &
                                                                   0, 1, 0, 0, &
                                                                   0, 0, 1, 0, &
                                                                   0, 0, 0, 1/),(/4,4/))
       integer, dimension(6)          :: latt_given
       real(kind=cp), dimension(3,11) :: latc
       integer, dimension(3)          :: tt
       integer                        :: i,j,l

       !---- Initializing ----!

       isystm  = 0
       isymce  = 1
       ibravl  = 0
       co      = 0.0
       latsy   = " "
       SpaceGen= " "
       if(present(lsym)) then
         if(.not. (lsym == "P" .or. lsym=="p")) SpaceGen=lsym
       end if

       if(len_trim(SpaceGen) == 0) then  !Test lattice translation only if lsym has not been provided
          latt_p=.true.                  !or if lsym="P"
          latt_a=.false.
          latt_b=.false.
          latt_c=.false.
          latt_i=.false.
          latt_r=.false.
          latt_f=.false.
          latt_z=.false.


          !---- Determine the type of lattice ----!
          !---- This is only in case an explicit translation generator is given.
          l=0
          do i=1,ng
             if (equal_matrix(ss(:,:,i),identidad(1:3,1:3),3)) then
                tt(1)=nint(12.0 * ts(1,i))   ! Translations x 12
                tt(2)=nint(12.0 * ts(2,i))
                tt(3)=nint(12.0 * ts(3,i))

                !---- Identity (I,0) is being processed ----!
                if (tt(1) == 0 .and. tt(2) == 0 .and. tt(3) == 0) cycle

                !---- Compare the translation part of the operator with tabulated array ----!
                l=l+1  !counts the number of non trivial centring vectors
                latt_given(:) = 0
                do j=1,6
                   if (equal_vector(tt,lattice(:,j),3)) then
                      latt_given(j)=1
                      select case (j)
                         case (1)
                            latt_a=.true.
                         case (2)
                            latt_b=.true.
                         case (3)
                            latt_c=.true.
                         case (4)
                            latt_i=.true.
                         case (5,6)
                            latt_r=.true.
                      end select
                      exit
                   end if
                end do
                latc(:,L) = ts(:,i)
                if(sum(latt_given) == 0) latt_z=.true.
             end if
          end do

          !---- Lattice Type ----!
          if ( (latt_a .and. latt_b .and. latt_c) .or. (latt_a .and. latt_b) .or. &
               (latt_a .and. latt_c) .or. (latt_b .and. latt_c) ) then
             latt_f=.true.
             latt_a=.false.
             latt_b=.false.
             latt_c=.false.
             latt_p=.false.
             latt_i=.false.
          end if

          if (latt_p) then
             SpaceGen="P"
             Ibravl  = 1
          end if
          if (latt_a) then
             SpaceGen="A"
             Ibravl  = 2
          end if
          if (latt_b) then
             SpaceGen="B"
             Ibravl  = 3
          end if
          if (latt_c) then
             SpaceGen="C"
             Ibravl  = 4
          end if
          if (latt_i) then
             SpaceGen="I"
             Ibravl  = 5
          end if
          if (latt_r) then
             SpaceGen="R"
             Ibravl  = 6
          end if
          if (latt_f) then
             SpaceGen="F"
             Ibravl  = 7
          end if
          if (latt_z) then
             SpaceGen="Z"
             Ibravl  = 8
          end if
       end if

       if (len_trim(SpaceGen) /= 0) then
          call latsym(SpaceGen,L,latc)
       else
          err_symm=.true.
          ERR_Symm_Mess=" Lattice Type couldn't be determined"
          return
       end if

       !---- Centre of symmetry ? ----!
       do i=1,ng
          if (equal_matrix(ss(:,:,i),-identidad(1:3,1:3),3)) then
             tt(1)=nint(12.0 * ts(1,i))
             tt(2)=nint(12.0 * ts(2,i))
             tt(3)=nint(12.0 * ts(3,i))

             if (tt(1) == 0 .and. tt(2) == 0 .and. tt(3) == 0) then
                isymce=2       ! Centric with -1 at origin
             else
                isymce=0       ! Centric without -1 at origin
                co=0.5*ts(:,i)
             end if
          end if
       end do

       !---- Determination of the crystalline system and Bravais lattice ----!
       call get_crystal_System(ng,ss,isystm,latsy(1:1))
       latsy(2:)=red(ibravl)

       return
    End Subroutine Get_So_From_Fix

    !!----
    !!---- Subroutine Get_So_From_Gener(Isystm,Isymce,Ibravl,Ng,Ss,Ts,Latsy,Co,Num_G,Spacegen)
    !!----    integer,                      intent(out)   :: ISYSTM    ! Out -> Number of the crystalline system
    !!----                                                                      (1:T, 2:M, 3:O, 4:T, 5:R-Trg, 6:H, 7:C)
    !!----    integer,                      intent(out)   :: ISYMCE    ! Out -> 0 Centric (-1 not at origin)
    !!----                                                                      1 Acentric
    !!----                                                                      2 Centric (-1 at origin)
    !!----    integer,                      intent(out)   :: IBRAVL    ! Out -> Index of the Bravais Lattice type
    !!----                                                                      1   2   3   4   5   6   7   8
    !!----                                                                      "P","A","B","C","I","R","F","Z"
    !!----    integer,                      intent(in out):: NG        !  In -> Number of defined generators
    !!----                                                             ! Out -> Number of symmetry operators
    !!----    integer, dimension(:,:,:),    intent(in out):: SS        !  In -> Rotation parts of the given generators  (3,3,48)
    !!----                                                             ! Out -> Rotation parts of the symmetry operators
    !!----    real(kind=cp),dimension(:,:), intent(in out):: TS        !  In -> Translation parts of the given generators  (3,48)
    !!----                                                             ! Out -> Translation parts of the symmetry operators
    !!----    character (len=2),            intent(out)   :: latsy     ! Out -> Bravais Lattice symbol
    !!----    real(kind=cp),dimension(3),   intent(out)   :: Co        ! Out -> Coordinates of origin
    !!----    integer,                      intent(out)   :: Num_g     ! Out -> Minimum number of generators
    !!----    character (len=1),            intent(out)   :: SpaceGen  ! Out -> Type of Cell
    !!----
    !!----    Calculates the whole set of symmetry operators from a set of given generators.
    !!----
    !!---- Update: February - 2005, February-2014 (JRC)
    !!
    Subroutine Get_SO_from_Gener(Isystm,Isymce,Ibravl,Ng,Ss,Ts,Latsy,Co,Num_g,SpaceGen,num_lat,lat_cent)
       !---- Arguments ----!
       integer,                                              intent(   out) :: Isystm
       integer,                                              intent(   out) :: Isymce
       integer,                                              intent(   out) :: Ibravl
       integer,                                              intent(in out) :: Ng
       integer, dimension(:,:,:),                            intent(in out) :: Ss ! (3,3,48)
       real(kind=cp),dimension(:,:),                         intent(in out) :: Ts ! (3  ,48)
       character (len=*),                                    intent(   out) :: Latsy
       real(kind=cp),dimension(3),                           intent(   out) :: Co
       integer,                                              intent(   out) :: Num_g
       character (len=*),                                    intent(   out) :: SpaceGen
       integer, optional,                                    intent(   out) :: num_lat
       real(kind=cp), dimension(:,:), allocatable, optional, intent(   out) :: lat_cent

       !---- Local Variables ----!
       real(kind=cp),dimension(3,192)  :: latc
       integer                         :: nlat_t
       logical :: latt_p, latt_a, latt_b, latt_c, latt_i, latt_r, latt_f, latt_z
       integer, dimension(6) :: latt_given
       character(len=*), dimension(8),  parameter :: red = &
                           (/"P","A","B","C","I","R","F","Z"/)
       integer, dimension(3,192)          :: lat_trans
       integer, dimension(3)              :: txyz
       integer, dimension(3,6), parameter :: lattice=reshape((/0,6,6, 6,0,6, &
                                                     6,6,0, 6,6,6, 8,4,4, 4,8,8/),(/3,6/))
       integer, dimension(4,4), parameter :: identidad = reshape((/1, 0, 0, 0, &
                                                                   0, 1, 0, 0, &
                                                                   0, 0, 1, 0, &
                                                                   0, 0, 0, 1/),(/4,4/))
       integer, dimension(4,4), parameter :: nulo      = reshape((/0, 0, 0, 0, &
                                                                   0, 0, 0, 0, &
                                                                   0, 0, 0, 0, &
                                                                   0, 0, 0, 0/),(/4,4/))
       integer, parameter             :: num_tab=24
       integer, dimension(4,4,num_tab):: tabla
       integer, dimension(4,4)        :: m1,m2
       integer, dimension(3)          :: tt,tt1
       real(kind=cp), parameter       :: lat_norm=12.0  !a multiple of 12 is compulsory
       integer,       parameter       :: ilat_norm=12, ilat_fact=ilat_norm*4
       integer :: i,j,k,n,nop,nopp,ipos,nt,ntp,L,jcen
       integer :: tx, ty, tz
       logical :: cen_found

       !---- Initializing ----!
       nop     = ng
       isystm  = 0
       isymce  = 1
       ibravl  = 0
       co      = 0.0
       latsy   = " "
       SpaceGen=" "

       latt_p=.true.
       latt_a=.false.
       latt_b=.false.
       latt_c=.false.
       latt_i=.false.
       latt_r=.false.
       latt_f=.false.
       latt_z=.false.
       !Set to zero all previous centring vectors
       nlat=0
       Ltr=0.0_cp

       !---- Redundances ----!
       do i=1,nop
          if (equal_matrix(ss(:,:,i),      nulo(1:3,1:3),3)) cycle  !ignore zero matrices
          if (equal_matrix(ss(:,:,i), identidad(1:3,1:3),3)) cycle  !Do not consider lattice translations
          if (equal_matrix(ss(:,:,i),-identidad(1:3,1:3),3)) cycle  !Do not consider inversion centre
          do j=i+1,nop
             if (equal_matrix(ss(:,:,j),      nulo(1:3,1:3),3)) cycle
             if (equal_matrix(ss(:,:,j), identidad(1:3,1:3),3)) cycle
             if (equal_matrix(ss(:,:,j),-identidad(1:3,1:3),3)) cycle

             !---- Traslation part ----!
             if (equal_matrix(ss(:,:,i),ss(:,:,j),3)) then
                tt =nint(lat_norm * ts(:,i))
                tt1=nint(lat_norm * ts(:,j))

                tx=mod(tt(1)-tt1(1)+ilat_fact,ilat_norm)        !?????????
                ty=mod(tt(2)-tt1(2)+ilat_fact,ilat_norm)
                tz=mod(tt(3)-tt1(3)+ilat_fact,ilat_norm)

                if (tx == 0 .and. ty == 0 .and. tz == 0) then
                   ss(:,:,j)=0
                   ts(:,j)=0.0
                   cycle
                else
                   tx=mod(tt(1)+tt1(1)+ilat_fact,ilat_norm)
                   ty=mod(tt(2)+tt1(2)+ilat_fact,ilat_norm)
                   tz=mod(tt(3)+tt1(3)+ilat_fact,ilat_norm)

                   if (tx == 0 .and. ty == 0 .and. tz == 0) then
                      ss(:,:,j)=0
                      ts(:,j)=0.0
                      cycle
                   else
                      ss(:,:,j)=identidad(1:3,1:3)
                      ts(1,j)=real(tx)/lat_norm
                      ts(2,j)=real(ty)/lat_norm
                      ts(3,j)=real(tz)/lat_norm
                      cycle
                   end if
                end if
             end if

             !---- Inversion part ----!
             if (equal_matrix(ss(:,:,i),-ss(:,:,j),3)) then
                tt =nint(lat_norm * ts(:,i))
                tt1=nint(lat_norm * ts(:,j))

                tx=mod(tt(1)-tt1(1)+ilat_fact,ilat_norm)        !?????????
                ty=mod(tt(2)-tt1(2)+ilat_fact,ilat_norm)
                tz=mod(tt(3)-tt1(3)+ilat_fact,ilat_norm)

                if (tx == 0 .and. ty == 0 .and. tz == 0) then
                   ss(:,:,j)=-identidad(1:3,1:3)
                   ts(:,j)=0.0
                   cycle
                else
                   tx=mod(tt(1)+tt1(1)+ilat_fact,ilat_norm)
                   ty=mod(tt(2)+tt1(2)+ilat_fact,ilat_norm)
                   tz=mod(tt(3)+tt1(3)+ilat_fact,ilat_norm)

                   if (tx == 0 .and. ty == 0 .and. tz == 0) then
                      ss(:,:,j)=-identidad(1:3,1:3)
                      ts(:,j)=0.0
                      cycle
                   else
                      ss(:,:,j)=-identidad(1:3,1:3)
                      ts(1,j)=real(tx)/lat_norm
                      ts(2,j)=real(ty)/lat_norm
                      ts(3,j)=real(tz)/lat_norm
                      cycle
                   end if
                end if
             end if

          end do
       end do

       !---- Determine the type of lattice before starting generation ----!
       !---- This is only in case an explicit translation generator is given.
       !---- This block construct the lattice vectors and remove the corresponding operators
       L=0
       do i=1,nop
          if (equal_matrix(ss(:,:,i),identidad(1:3,1:3),3)) then
             tt=nint(lat_norm * ts(:,i))   ! Translations x 12

             !---- Identity (I,0) is being processed ----!
             if (tt(1) == 0 .and. tt(2) == 0 .and. tt(3) == 0) then
                ss(:,:,i)=0
                ts(:,i)=0.0
                cycle
             end if
             !---- Compare the translation part of the operator with tabulated array ----!
             L=L+1  !counts the number of non trivial centring vectors
             latt_given(:) = 0
             do j=1,6
                if (equal_vector(tt,(ilat_norm/12)*lattice(:,j),3)) then
                   latt_given(j) = 1
                   select case (j)
                      case (1)
                         latt_a=.true.
                      case (2)
                         latt_b=.true.
                      case (3)
                         latt_c=.true.
                      case (4)
                         latt_i=.true.
                      case (5,6)
                         latt_r=.true.
                   end select
                   exit
                end if
             end do
             latc(:,L) = ts(:,i)
             ss(:,:,i)=0   !Removing the operators with pure translations
             ts(:,i)=0.0
             if(sum(latt_given) == 0) latt_z = .true.
          end if
       end do

       !---- Lattice Type ----!
       if ( (latt_a .and. latt_b .and. latt_c) .or. (latt_a .and. latt_b) .or. &
            (latt_a .and. latt_c) .or. (latt_b .and. latt_c) ) then
          latt_f=.true.
          latt_a=.false.
          latt_b=.false.
          latt_c=.false.
          latt_p=.false.
          latt_i=.false.
       end if
       if(latt_f .and. latt_r) latt_z=.true.
       if(latt_i .and. latt_r) latt_z=.true.
       if(latt_a .and. latt_r) latt_z=.true.
       if(latt_b .and. latt_r) latt_z=.true.
       if(latt_c .and. latt_r) latt_z=.true.

       if(latt_z) then
          latt_f=.false.
          latt_a=.false.
          latt_b=.false.
          latt_c=.false.
          latt_p=.false.
          latt_i=.false.
          latt_r=.false.
       end if

       if (latt_p) then
          SpaceGen="P"
          Ibravl  = 1
          nlat_t=0
       end if
       if (latt_a) then
          SpaceGen="A"
          Ibravl  = 2
          nlat_t=1
          latc(:,1)=(/0.0,0.5,0.5/)
       end if
       if (latt_b) then
          SpaceGen="B"
          Ibravl  = 3
          nlat_t=1
          latc(:,1)=(/0.5,0.0,0.5/)
       end if
       if (latt_c) then
          SpaceGen="C"
          Ibravl  = 4
          nlat_t=1
          latc(:,1)=(/0.5,0.5,0.0/)
       end if
       if (latt_i) then
          SpaceGen="I"
          Ibravl  = 5
          nlat_t=1
          latc(:,1)=(/0.5,0.5,0.5/)
       end if
       if (latt_r) then
          SpaceGen="R"
          Ibravl  = 6
          nlat_t=2
          latc(:,1)=(/ 2.0/3.0, 1.0/3.0, 1.0/3.0 /)
          latc(:,2)=(/ 1.0/3.0, 2.0/3.0, 2.0/3.0 /)
       end if
       if (latt_f) then
          SpaceGen="F"
          Ibravl  = 7
          nlat_t=3
          latc(:,1)=(/0.5,0.5,0.0/)
          latc(:,2)=(/0.5,0.0,0.5/)
          latc(:,3)=(/0.0,0.5,0.5/)
       end if
       if (latt_z) then
          Ibravl  = 8
          !Determine here the total number of non-trivial centring vectors
          call get_centring_vectors(L,latc,SpaceGen)
          nlat_t=L
          do i=1,nlat_t
             lat_trans(:,i)=maxval(nint(lat_norm*latc(:,i)))
          end do
       end if

       if(present(num_lat) .and. present(lat_cent)) then
         num_lat=nlat_t
         if(allocated(lat_cent)) deallocate(lat_cent)
         allocate(lat_cent(3,num_lat))
         lat_cent(:,1:num_lat)=latc(:,1:num_lat)
       end if

       if (len_trim(SpaceGen) /= 0) then
          !write(*,"(a)") " => Lattice centrings at Get_SO_from_Gener: "
          !do i=1,nlat_t
          !   write(*,"(i8,3f14.5)")i,latc(:,i)
          !end do
          call latsym(SpaceGen,nlat_t,latc)
       else
          err_symm=.true.
          ERR_Symm_Mess=" Lattice Type couldn't be determined"
          return
       end if

       !---- Removing Centre of symmetry if found ----!
       do i=1,nop
          if (equal_matrix(ss(:,:,i),-identidad(1:3,1:3),3)) then
             tt=nint(lat_norm * ts(:,i))

             if (tt(1) == 0 .and. tt(2) == 0 .and. tt(3) == 0) then
                isymce=2       ! Centric with -1 at origin
             else
                if(isymce /= 2) then !do that only if a centre has not been found at the origin
                  isymce=0       ! Centric without -1 at origin
                  co=0.5*ts(:,i)
                end if
             end if
             ss(:,:,i)=0
             ts(:,i)=0.0
          end if
       end do

       !---- Purge the starting list of generators ----!
       !     by diminishing the number of generators to the strictly
       !     assymmetric block without translations and centre
       nopp=nop
       ipos=1
       do
          do i=ipos,nopp
             if (equal_matrix(ss(:,:,i),nulo(1:3,1:3),3)) then
                do j=i+1,nopp
                   ss(:,:,j-1)= ss(:,:,j)
                   ts(:,j-1)  = ts(:,j)
                end do
                ss(:,:,nopp)=0
                ts(:,nopp)=0.0
                nopp=nopp-1
                ipos=i
                exit
             end if
          end do
          if (i >= nopp) exit
       end do
       if (equal_matrix(ss(:,:,nopp),nulo(1:3,1:3),3)) nopp=nopp-1
       nop=nopp    !nop is the number of generators without centre of symmetry
                   !and lattice centrings

       !---- Now we have an eventually shorter list of generators and we
       !---- know if a centre of symmetry or lattice centrings were
       !---- among the given generators.

       !---- Creation of the symmetry operators table ----!
       !write(*,"(a,i6)") " => Number of generators (no centre/no lattice) to start the table: ",nop
       !do i=1,nop
       !  write(*,"(9i4,3f8.4)") ss(:,:,i), ts(:,i)
       !end do
       !---- Initializing ----!
       nt=1
       tabla=0
       tabla(4,4,:)=1
       do i=1,4               !---- Identity operator is the first one ----!
          tabla(i,i,1)=1
       end do

       !---- Put operators in the table ----!
       do i=1,nop
          nt=nt+1
          tabla(1:3,1:3,nt)= ss(:,:,i)
          tabla(1:3,4,nt)  = mod(nint(lat_norm*ts(1:3,i))+ilat_fact,ilat_norm)
       end do

       num_g=nop   !Minimum number of generators (except inversion)

       !---- Generate power operations from generators ----!
       do i=2,nt
             ntp=axes_rotation(tabla(1:3,1:3,i))    ! Determine the order of the generator
          if (ntp == -1 .or. ntp == -3) then
             ntp=-2*ntp
          else
             ntp=abs(ntp)
          end if

          m1=tabla(:,:,i)
          m2=identidad

          p1:do j=1,ntp-1
             m2=matmul(m2,m1)
             m2(:,4)=mod(m2(:,4)+ilat_fact,ilat_norm)

             !---- Check if the generated operation is already in the table
             do k=1,nt
                if (equal_matrix(tabla(:,:,k),m2,4)) cycle p1
             end do

             !---- Eliminating lattice contribution if necessary ----!
             select case (ibravl)
                 case (2)
                    ty=m2(2,4)
                    tz=m2(3,4)
                    if (ty >= 6 .and. tz >= 6) then
                       ty=mod(m2(2,4),6)
                       tz=mod(m2(3,4),6)

                       if (ty == 0 .and. tz == 0) then
                          cycle p1
                       else
                          m2(2,4)=ty
                          m2(3,4)=tz
                       end if
                    end if

                 case (3)
                    tx=m2(1,4)
                    tz=m2(3,4)
                    if (tx >= 6 .and. tz >=6) then
                       tx=mod(m2(1,4),6)
                       tz=mod(m2(3,4),6)

                       if (tx == 0 .and. tz == 0) then
                          cycle p1
                       else
                          m2(1,4)=tx
                          m2(3,4)=tz
                       end if
                    end if

                 case (4)
                    tx=m2(1,4)
                    ty=m2(2,4)
                    if (tx >=6 .and. ty >= 6) then
                       tx=mod(m2(1,4),6)
                       ty=mod(m2(2,4),6)

                       if (tx == 0 .and. ty == 0) then
                          cycle p1
                       else
                          m2(1,4)=tx
                          m2(2,4)=ty
                       end if
                    end if

                 case (5)
                    tx=m2(1,4)
                    ty=m2(2,4)
                    tz=m2(3,4)
                    if (tx >=6 .and. ty >=6 .and. tz >=6) then
                       tx=mod(m2(1,4),6)
                       ty=mod(m2(2,4),6)
                       tz=mod(m2(3,4),6)

                       if (tx == 0 .and. ty == 0 .and. tz == 0) then
                          cycle p1
                       else
                          m2(1,4)=tx
                          m2(2,4)=ty
                          m2(3,4)=tz
                       end if
                    end if

                 case (6)
                    tx=m2(1,4)
                    ty=m2(2,4)
                    tz=m2(3,4)
                    if (tx >=8 .and. ty >=4 .and. tz >=4) then
                       tx=mod(m2(1,4),8)
                       ty=mod(m2(2,4),4)
                       tz=mod(m2(3,4),4)
                       if (tx == 0 .and. ty == 0 .and. tz == 0) then
                          cycle p1
                       else
                          m2(1,4)=tx
                          m2(2,4)=ty
                          m2(3,4)=tz
                       end if

                    else if (tx >=4 .and. ty >=8 .and. tz >=8) then
                       tx=mod(m2(1,4),4)
                       ty=mod(m2(2,4),8)
                       tz=mod(m2(3,4),8)
                       if (tx == 0 .and. ty == 0 .and. tz == 0) then
                          cycle p1
                       else
                          m2(1,4)=tx
                          m2(2,4)=ty
                          m2(3,4)=tz
                       end if
                    end if

                 case (7)
                    tx=m2(1,4)
                    ty=m2(2,4)
                    tz=m2(3,4)
                    if (ty >= 6 .and. tz >=6) then
                       ty=mod(m2(2,4),6)
                       tz=mod(m2(3,4),6)
                       if (tx == 0 .and. ty == 0 .and. tz == 0) then
                          cycle p1
                       else
                          m2(2,4)=ty
                          m2(3,4)=tz
                       end if
                    else if (tx >=6 .and. tz >=6) then
                       tx=mod(m2(1,4),6)
                       tz=mod(m2(3,4),6)
                       if (tx == 0 .and. ty == 0 .and. tz == 0) then
                          cycle p1
                       else
                          m2(1,4)=tx
                          m2(3,4)=tz
                       end if
                    else if (tx >=6 .and. ty >=6) then
                       tx=mod(m2(1,4),6)
                       ty=mod(m2(2,4),6)
                       if (tx == 0 .and. ty == 0 .and. tz == 0) then
                          cycle p1
                       else
                          m2(1,4)=tx
                          m2(2,4)=ty
                       end if
                    end if

                 Case(8)
                    !eliminate translations
                    txyz=m2(1:3,4)
                    do k=1,nlat_t
                       if(txyz(1) >= lat_trans(1,k) .and. txyz(2) >= lat_trans(2,k) .and. txyz(3) >= lat_trans(3,k)) then
                          txyz(:)=mod(m2(1:3,4), lat_trans(:,k))
                          if(txyz(1) == 0 .and. txyz(2) == 0 .and.txyz(3) == 0 ) then
                             cycle p1
                          else
                              m2(1:3,4)=txyz(:)
                          end if
                       end if
                    end do

             end select

             nt=nt+1
             if (nt > num_tab) then
                err_symm=.true.
                ERR_Symm_Mess=" Dimension of Table exceeded (I)"
                return
             end if
             tabla(:,:,nt)=m2
          end do p1

       end do

       !write(*,"(a,i6)") " => Number of terms in the table at stage I: ",nt
       !do i=1,nt
       !  write(*,"(9i3,3i6)") tabla(1:3,1:3,i),tabla(1:3,4,i)
       !end do
       !write(*,"(//a//)") " => MULTIPLICATION OF GENERATORS: "

       !---- Multiplications between generators ----!
       do
          if (nt == 1) exit
          n=nt

          do i=1,n
             p2:do j=i,n

                m2=matmul(tabla(:,:,i),tabla(:,:,j))
                m2(:,4)=mod(m2(:,4)+ilat_fact,ilat_norm)

                !---- Eliminating lattice contribution if necessary ----!
                select case (ibravl)
                   case (2)
                      ty=m2(2,4)
                      tz=m2(3,4)
                      if (ty >= 6 .and. tz >= 6) then
                         m2(2,4)=mod(m2(2,4),6)
                         m2(3,4)=mod(m2(3,4),6)
                      end if

                   case (3)
                      tx=m2(1,4)
                      tz=m2(3,4)
                      if (tx >= 6 .and. tz >= 6) then
                         m2(1,4)=mod(m2(1,4),6)
                         m2(3,4)=mod(m2(3,4),6)
                      end if

                   case (4)
                      tx=m2(1,4)
                      ty=m2(2,4)
                      if (tx >= 6 .and. ty >= 6) then
                         m2(1,4)=mod(m2(1,4),6)
                         m2(2,4)=mod(m2(2,4),6)
                      end if

                   case (5)
                      tx=m2(1,4)
                      ty=m2(2,4)
                      tz=m2(3,4)
                      if (tx >= 6 .and. ty >= 6 .and. tz >= 6) then
                         m2(1,4)=mod(m2(1,4),6)
                         m2(2,4)=mod(m2(2,4),6)
                         m2(3,4)=mod(m2(3,4),6)
                      end if

                   case (6)
                      tx=m2(1,4)
                      ty=m2(2,4)
                      tz=m2(3,4)
                      if (tx >=8 .and. ty >=4 .and. tz >=4) then
                         m2(1,4)=mod(m2(1,4),8)
                         m2(2,4)=mod(m2(2,4),4)
                         m2(3,4)=mod(m2(3,4),4)
                      else if (tx >=4 .and. ty >=8 .and. tz >=8) then
                         m2(1,4)=mod(m2(1,4),4)
                         m2(2,4)=mod(m2(2,4),8)
                         m2(3,4)=mod(m2(3,4),8)
                      end if

                   case (7)
                      tx=m2(1,4)
                      ty=m2(2,4)
                      tz=m2(3,4)
                      if (ty >= 6 .and. tz >=6) then
                         m2(2,4)=mod(m2(2,4),6)
                         m2(3,4)=mod(m2(3,4),6)
                      else if (tx >=6 .and. tz >=6) then
                         m2(1,4)=mod(m2(1,4),6)
                         m2(3,4)=mod(m2(3,4),6)
                      else if (tx >=6 .and. ty >=6) then
                         m2(1,4)=mod(m2(1,4),6)
                         m2(2,4)=mod(m2(2,4),6)
                      end if

                   case(8)

                      txyz=m2(1:3,4)
                      do k=1,nlat_t
                         if(txyz(1) >= lat_trans(1,k) .and. txyz(2) >= lat_trans(2,k) .and. txyz(3) >= lat_trans(3,k)) then
                            m2(1:3,4)=mod(m2(1:3,4), lat_trans(:,k))
                         end if
                      end do

                end select

                do k=1,nt
                   if ( equal_matrix(m2(:,:),tabla(:,:,k),4) ) cycle p2
                   if ( equal_matrix(m2(:,:),tabla(:,:,k),3) ) then
                      tx=m2(1,4)+tabla(1,4,k)
                      ty=m2(2,4)+tabla(2,4,k)
                      tz=m2(3,4)+tabla(3,4,k)
                      tx=mod(tx,ilat_norm)
                      ty=mod(ty,ilat_norm)
                      tz=mod(tz,ilat_norm)
                      select case (ibravl)
                          case (2)
                             if (ty == 6 .and. tz == 6) cycle p2
                             if (ty == 0 .and. tz == 0) cycle p2

                          case (3)
                             if (tx == 6 .and. tz == 6) cycle p2
                             if (tx == 0 .and. tz == 0) cycle p2

                          case (4)
                             if (tx == 6 .and. ty == 6) cycle p2
                             if (tx == 0 .and. ty == 0) cycle p2

                          case (5)
                             if (tx == 6 .and. ty == 6 .and. tz == 6) cycle p2
                             if (tx == 0 .and. ty == 0 .and. tz == 0) cycle p2

                          case (6)
                             if (tx == 8 .and. ty == 4 .and. tz == 4) cycle p2
                             if (tx == 4 .and. ty == 8 .and. tz == 8) cycle p2
                             if (tx == 0 .and. ty == 0 .and. tz == 0) cycle p2

                          case (7)
                             if (ty == 6 .and. tz == 6) cycle p2
                             if (tx == 6 .and. tz == 6) cycle p2
                             if (tx == 6 .and. ty == 6) cycle p2

                             if (ty == 0 .and. tz == 0) cycle p2
                             if (tx == 0 .and. tz == 0) cycle p2
                             if (tx == 0 .and. ty == 0) cycle p2

                          case (8)
                             txyz=(/tx,ty,tz/)
                             if(txyz(1) == 0 .and. txyz(2) == 0 .and. txyz(3) == 0) cycle p2
                             do l=1,nlat_t
                                if(txyz(1) == lat_trans(1,l) .and. txyz(2) == lat_trans(2,l) .and. txyz(3) == lat_trans(3,l)) cycle p2
                             end do

                      end select

                      tx=m2(1,4)-tabla(1,4,k)
                      ty=m2(2,4)-tabla(2,4,k)
                      tz=m2(3,4)-tabla(3,4,k)
                      tx=mod(tx+ilat_fact,ilat_norm)
                      ty=mod(ty+ilat_fact,ilat_norm)
                      tz=mod(tz+ilat_fact,ilat_norm)

                      select case (ibravl)
                          case (2)
                             if (ty == 6 .and. tz == 6) cycle p2
                             if (ty == 0 .and. tz == 0) cycle p2

                          case (3)
                             if (tx == 6 .and. tz == 6) cycle p2
                             if (tx == 0 .and. tz == 0) cycle p2

                          case (4)
                             if (tx == 6 .and. ty == 6) cycle p2
                             if (tx == 0 .and. ty == 0) cycle p2

                          case (5)
                             if (tx == 6 .and. ty == 6 .and. tz == 6) cycle p2
                             if (tx == 0 .and. ty == 0 .and. tz == 0) cycle p2

                          case (6)
                             if (tx == 8 .and. ty == 4 .and. tz == 4) cycle p2
                             if (tx == 4 .and. ty == 8 .and. tz == 8) cycle p2
                             if (tx == 0 .and. ty == 0 .and. tz == 0) cycle p2

                          case (7)
                             if (ty == 6 .and. tz == 6) cycle p2
                             if (tx == 6 .and. tz == 6) cycle p2
                             if (tx == 6 .and. ty == 6) cycle p2

                             if (ty == 0 .and. tz == 0) cycle p2
                             if (tx == 0 .and. tz == 0) cycle p2
                             if (tx == 0 .and. ty == 0) cycle p2

                          case (8)
                             txyz=(/tx,ty,tz/)
                             if(txyz(1) == 0 .and. txyz(2) == 0 .and. txyz(3) == 0) cycle p2
                             do l=1,nlat_t
                                if(txyz(1) == lat_trans(1,l) .and. txyz(2) == lat_trans(2,l) .and. txyz(3) == lat_trans(3,l)) cycle p2
                             end do

                      end select
                   end if
                end do

                nt=nt+1
                if (nt > num_tab) then
                   err_symm=.true.
                   write(unit=ERR_Symm_Mess,fmt="(a,i5)") " Dimension of Table exceeded (II): ",nt
                   return
                end if
                !new operator
                tabla(:,:,nt)=m2
                !write(*,"(a,i4,a,9i3,3i6)")  "  Op:",nt," -> ", tabla(1:3,1:3,nt),tabla(1:3,4,nt)
             end do p2
          end do

          if (n == nt) exit

       end do

       !---- Carga Final ----!
       ng=nt
       do i=1,nt
          ss(:,:,i)=tabla(1:3,1:3,i)
          ts(:,i)  = real(tabla(1:3,4,i))/lat_norm
       end do

       ! Check anomalous cases where the order of generators has produced
       ! a centre of symmetry at the end.
       ! In some cases the number of operators were wrong => the list of
       ! operators contained centrosymmetric related items.

       cen_found = .false.
       do i=1,nt
          if ( equal_matrix(ss(:,:,i),-identidad,3) ) then
             cen_found=.true.
             jcen=i
             exit
          end if
       end do

       if(cen_found) then
           ng=nt/2     !in all cases only half operators are needed
           if (isymce == 1) then
              tt=tabla(1:3,4,jcen)
              if (tt(1) == 0 .and. tt(2) == 0 .and. tt(3) == 0) then
                 isymce=2       ! Centric with -1 at origin
              else
                 isymce=0       ! Centric without -1 at origin
                 co=0.5*ts(:,jcen)
              end if
           end if
       end if

       !---- Determination of the crystalline system and Bravais lattice ----!
       call get_crystal_System(ng,ss,isystm,latsy(1:1))
       latsy(2:)=red(ibravl)

       return
    End Subroutine Get_SO_from_Gener

    !!----
    !!---- Subroutine Get_So_From_Hall(Isystm,Isymce,Ibravl,Ng,Ss,Ts,Latsy,Co,Num_G,Hall)
    !!----    integer,                   intent(out)  :: ISYSTM    ! Out -> Number of the crystalline system
    !!----                                                                  (1:T, 2:M, 3:O, 4:T, 5:R-Trg, 6:H, 7:C)
    !!----    integer,                   intent(out)  :: ISYMCE    ! Out -> 0 Centric (-1 not at origin)
    !!----                                                                  1 Acentric
    !!----                                                                  2 Centric (-1 at origin)
    !!----    integer,                   intent(out)  :: IBRAVL    ! Out -> Index of the Bravais Lattice type
    !!----                                                                  1   2   3   4   5   6   7
    !!----                                                                 "P","A","B","C","I","R","F"
    !!----    integer,                   intent(out)  :: NG        ! Out -> Number of symmetry operators
    !!----    real(kind=cp),    dimension(:,:),   intent(out)  :: TS        ! Out -> Translation parts of the symmetry operators  (3,48)
    !!----    integer, dimension(:,:,:), intent(out)  :: SS        ! Out -> Rotation parts of the symmetry operators     (3,3,48)
    !!----    character (len=2),         intent(out)  :: latsy     ! Out -> Bravais lattice symbol
    !!----    real(kind=cp), dimension(3),        intent(out)  :: Co        ! Out -> Coordinates of symmetry center
    !!----    integer,                   intent(out)  :: num_g     ! Out -> Number of generators
    !!----    character (len=20),        intent( in)  :: Hall      !  In -> Hall Spacegroup symbol
    !!----
    !!----    Subroutine to get all the information contained in the Hall symbol. This
    !!----    routine to interpret the Hall symbol for a space group.
    !!--..    (Author:Javier Gonzalez-Platas)
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Get_SO_from_Hall(Isystm,Isymce,Ibravl,Ng,Ss,Ts,Latsy,Co,Num_g,Hall)
       !---- Arguments ----!
       integer,                   intent(out) :: Isystm
       integer,                   intent(out) :: Isymce
       integer,                   intent(out) :: Ibravl
       integer,                   intent(out) :: Ng
       integer, dimension(:,:,:), intent(out) :: Ss  !(3,3,48)
       real(kind=cp),    dimension(:,:),   intent(out) :: Ts  !(3,48)
       character (len= 2),        intent(out) :: Latsy
       real(kind=cp),    dimension(3),     intent(out) :: Co
       integer,                   intent(out) :: Num_g
       character (len=*),         intent( in) :: Hall

       !----Local variables ----!
       character (len=16)                         :: group
       character(len=*), dimension(7),  parameter :: red = &
                         (/"P","A","B","C","I","R","F"/)
       character(len=*), dimension(13), parameter :: traslacion =&
                         (/"A","B","C","N","U","V","W","D","1","2","3","4","5"/)
       character(len=*), dimension(6),  parameter :: ejes_rotacion = &
                         (/"X","Y","Z","'","""","*"/)
       character(len=*), dimension(5),  parameter :: rotacion=(/"1","2","3","4","6"/)
       character(len=*), dimension(5),  parameter :: shift=(/"1","2","3","4","5"/)
       integer, dimension(3,13), parameter :: tras_val=reshape((/6,0,0, 0,6,0, &
                                       0,0,6, 6,6,6, 3,0,0, 0,3,0, 0,0,3, 3,3,3, &
                                       1,0,0, 2,0,0, 3,0,0, 4,0,0, 5,0,0/),(/3,13/))
       integer, dimension(3,3), parameter  :: x_1   = reshape( &
                                 (/ 1, 0, 0,  0, 1, 0,  0, 0, 1/),(/3,3/))
       integer, dimension(3,3), parameter  :: y_1   = reshape( &
                                 (/ 1, 0, 0,  0, 1, 0,  0, 0, 1/),(/3,3/))
       integer, dimension(3,3), parameter  :: z_1   = reshape( &
                                 (/ 1, 0, 0,  0, 1, 0,  0, 0, 1/),(/3,3/))
       integer, dimension(3,3), parameter  :: x_2   = reshape( &
                                 (/ 1, 0, 0,  0,-1, 0,  0, 0,-1/),(/3,3/))
       integer, dimension(3,3), parameter  :: y_2   = reshape( &
                                 (/-1, 0, 0,  0, 1, 0,  0, 0,-1/),(/3,3/))
       integer, dimension(3,3), parameter  :: z_2   = reshape( &
                                 (/-1, 0, 0,  0,-1, 0,  0, 0, 1/),(/3,3/))
       integer, dimension(3,3), parameter  :: x_3   = reshape( &
                                 (/ 1, 0, 0,  0, 0, 1,  0,-1,-1/),(/3,3/))
       integer, dimension(3,3), parameter  :: y_3   = reshape( &
                                 (/-1, 0,-1,  0, 1, 0,  1, 0, 0/),(/3,3/))
       integer, dimension(3,3), parameter  :: z_3   = reshape( &
                                 (/ 0, 1, 0, -1,-1, 0,  0, 0, 1/),(/3,3/))
       integer, dimension(3,3), parameter  :: x_4   = reshape( &
                                 (/ 1, 0, 0,  0, 0, 1,  0,-1, 0/),(/3,3/))
       integer, dimension(3,3), parameter  :: y_4   = reshape( &
                                 (/ 0, 0,-1,  0, 1, 0,  1, 0, 0/),(/3,3/))
       integer, dimension(3,3), parameter  :: z_4   = reshape( &
                                 (/ 0, 1, 0, -1, 0, 0,  0, 0, 1/),(/3,3/))
       integer, dimension(3,3), parameter  :: x_6   = reshape( &
                                 (/ 1, 0, 0,  0, 1, 1,  0,-1, 0/),(/3,3/))
       integer, dimension(3,3), parameter  :: y_6   = reshape( &
                                 (/ 0, 0,-1,  0, 1, 0,  1, 0, 1/),(/3,3/))
       integer, dimension(3,3), parameter  :: z_6   = reshape( &
                                 (/ 1, 1, 0, -1, 0, 0,  0, 0, 1/),(/3,3/))
       integer, dimension(3,3), parameter  :: x_2p  = reshape( &
                                 (/-1, 0, 0,  0, 0,-1,  0,-1, 0/),(/3,3/))
       integer, dimension(3,3), parameter  :: y_2p  = reshape( &
                                 (/ 0, 0,-1,  0,-1, 0, -1, 0, 0/),(/3,3/))
       integer, dimension(3,3), parameter  :: z_2p  = reshape( &
                                 (/ 0,-1, 0, -1, 0, 0,  0, 0,-1/),(/3,3/))
       integer, dimension(3,3), parameter  :: x_2pp = reshape( &
                                 (/-1, 0, 0,  0, 0, 1,  0, 1, 0/),(/3,3/))
       integer, dimension(3,3), parameter  :: y_2pp = reshape( &
                                 (/ 0, 0, 1,  0,-1, 0,  1, 0, 0/),(/3,3/))
       integer, dimension(3,3), parameter  :: z_2pp = reshape( &
                                 (/ 0, 1, 0,  1, 0, 0,  0, 0,-1/),(/3,3/))
       integer, dimension(3,3), parameter  :: xyz_3 = reshape( &
                                 (/ 0, 1, 0,  0, 0, 1,  1, 0, 0/),(/3,3/))
       integer, dimension(4,4), parameter :: identidad = reshape((/1, 0, 0, 0, &
                                                                     0, 1, 0, 0, &
                                                                     0, 0, 1, 0, &
                                                                     0, 0, 0, 1/),(/4,4/))
       integer, parameter             :: num_tab=24
       integer, dimension(4,4,num_tab):: tabla
       integer, dimension(4,4,4) :: gener
       integer, dimension(4,4)   :: sn, snp
       integer, dimension(4,4)   :: m1,m2
       integer, dimension(4)     :: num_rot
       integer, dimension(4)     :: num_eje
       integer, dimension(4)     :: num_tras
       integer, dimension(3)     :: vtras

       logical                   :: only_rot
       integer                   :: i,j,k,n,nt,ntp,npos
       integer                   :: pos_ini,pos_act,pos_fin,ini,fin
       integer                   :: ngen, neje, nrot, signo
       integer                   :: tx,ty,tz

       !---- Inicio ----!
       isystm=0
       isymce=1
       ibravl=0
       ng=0
       ss=0
       ts=0.0
       co=0.0
       latsy=" "
       num_g=1
       call init_err_symm()

       !---- Convert to Upper case ----!
       group = hall
       call ucase(group)
       group=adjustl(group)

       pos_ini=1
       pos_act=1
       pos_fin=len_trim(group)

       !---- Centric / Acentric ----!
       if (group(pos_ini:pos_ini) == "-") then
          isymce=2
          pos_ini=pos_ini+1
       else
       !
       ! Determine first if there are parenthesis
       !
          i=index(group(1:pos_fin),"(")
          if(i == 0) then
             npos=index(group(1:pos_fin),"-1",back=.true.)
          else
             npos=index(group(1:i),"-1",back=.true.)
          end if
          if (npos /= 0) then
             vtras=0
             do i=npos+2,pos_fin
                do j=1,13
                   if (group(i:i) == traslacion(j)) then
                      if (j < 9) then       ! a b c n u v w d
                         vtras=vtras+tras_val(:,j)
                      end if
                   end if
                end do
             end do

             vtras=mod(vtras,12)
             co=real(vtras)/12.0
             co=0.5*co
             if (vtras(1) == 0 .and. vtras(2) == 0 .and. vtras(3) == 0) then
                isymce=2
             else
                isymce=0
                pos_fin=npos-2
             end if
          end if
       end if

       !---- Tipo de Celda ----!
       ibravl=0
       do i=1,7
          if (group(pos_ini:pos_ini) /= red(i)) cycle
          ibravl=i
          exit
       end do
       if (ibravl == 0) then
          err_symm=.true.
          ERR_Symm_Mess=" IBRAVL is Equal Zero"
          return
       end if
       pos_ini=pos_ini+2

       !---- Determinacion de Generadores ----!
       gener=0
       gener(4,4,:)=1

       num_rot=0
       num_eje=0
       num_tras=1
       ngen=0

       do
          if (pos_ini > pos_fin) exit        ! Fin de caracterizacion

          pos_act=index(group(pos_ini:pos_fin)," ")

          only_rot=.false.
          nrot=0
          vtras=0
          signo=1
          neje=0
          if (ngen==0) neje=3                ! Eje C

          if (pos_act /=0) then
             ini=pos_ini
             if (pos_act /= 1) then
                fin=pos_ini+pos_act-2
             else
                fin=pos_ini+pos_act-1
             end if
          else
             ini=pos_ini
             fin=pos_fin
          end if

          !---- Desplazamiento del origen ----!
          if (group(ini:ini)=="(") then
             npos=0
             do i=ini+1,pos_fin-1            ! Eliminamos parentesis
                if (group(i:i)==" ") cycle
                if (group(i:i)=="-") then
                   signo=-1
                   cycle
                end if

                npos=npos+1
                do j=1,5
                   if (group(i:i)==shift(j)) then
                      vtras(npos)=j*signo
                      signo=1
                      exit
                   end if
                end do

             end do

             sn=0
             snp=0
             do i=1,4
                sn(i,i)=1
                snp(i,i)=1
             end do
             do i=1,3
                sn(i,4) =  vtras(i)
                snp(i,4)= -vtras(i)
             end do

             gener(:,:,ngen)=matmul(sn,gener(:,:,ngen))
             gener(:,:,ngen)=matmul(gener(:,:,ngen),snp)

             exit         ! Fin de busqueda
          end if

          !---- Eje de rotacion Propio/Impropio ----!
          if (group(ini:ini)=="-") then
             signo=-1
             ini=ini+1
          end if

          !---- Eje de rotacion ----!
          do j=1,5
             if (group(ini:ini) /= rotacion(j)) cycle
             nrot=j
             exit
          end do
          if (nrot==0) then
             err_symm=.true.
             return
          end if

          if (ini ==fin) only_rot=.true.
          ini=ini+1
          ini=min(ini,fin)

          !---- Direccion de Rotacion y Traslaciones ----!
          do i=ini,fin
             do j=1,6
                if (group(i:i) /= ejes_rotacion(j)) cycle
                neje=j
                exit
             end do

             if (neje == 0) then
                select case (ngen)
                    case (0)
                       neje=3
                    case (1)
                       neje=1
                       if (nrot == 2) then
                          if (num_rot(1)==2 .or. num_rot(1)==4) neje=1
                          if (num_rot(1)==3 .or. num_rot(1)==5) neje=4
                       end if
                    case (2)
                       neje=1
                       if (nrot == 3) neje=6
                    case (3)
                       neje=1
                end select
             end if


             if (only_rot) cycle    ! Solo eje de rotacion

             do j=1,13
                if (group(i:i) == traslacion(j)) then
                   if (j < 9) then       ! a b c n u v w d
                      vtras=vtras+tras_val(:,j)
                      select case (j)
                          case (5:8)
                          num_tras(ngen+1)=num_tras(ngen+1)*3
                      end select

                   else                  ! 1 2 3 4 6
                      if (nrot ==3 .or. nrot==4 .or. nrot==5) then
                         n=j-8
                         if (nrot==5) then
                            n=n*2
                         else
                            n=n*12/nrot
                         end if
                         vtras=0
                         vtras(neje)=n

                         num_tras(ngen+1)=num_tras(ngen+1)*(nrot-1)
                      else
                         err_symm=.true.
                         return
                      end if
                   end if
                end if
             end do

          end do

          !---- Cargando informacion ----!
          ngen=ngen+1
          num_rot(ngen)=nrot
          num_eje(ngen)=neje

          select case (nrot)
              case (1)
                 select case (neje)
                     case (1)
                        gener(1:3,1:3,ngen)=x_1*signo
                     case (2)
                        gener(1:3,1:3,ngen)=y_1*signo
                     case (3)
                        gener(1:3,1:3,ngen)=z_1*signo
                     case (4:6)
                        err_symm=.true.
                        return
                 end select
                 gener(1:3,4,ngen)=vtras

              case (2)
                 select case (neje)
                     case (1)
                        gener(1:3,1:3,ngen)=x_2*signo
                     case (2)
                        gener(1:3,1:3,ngen)=y_2*signo
                     case (3)
                        gener(1:3,1:3,ngen)=z_2*signo
                     case (4)
                        select case (num_eje(1))
                            case (1)
                               gener(1:3,1:3,ngen)=x_2p*signo
                            case (2)
                               gener(1:3,1:3,ngen)=y_2p*signo
                            case (3)
                               gener(1:3,1:3,ngen)=z_2p*signo
                            case (6)
                               gener(1:3,1:3,ngen)=z_2p*signo
                            case default
                               err_symm=.true.
                               return
                        end select
                     case (5)
                        select case (num_eje(1))
                            case (1)
                               gener(1:3,1:3,ngen)=x_2pp*signo
                            case (2)
                               gener(1:3,1:3,ngen)=y_2pp*signo
                            case (3)
                               gener(1:3,1:3,ngen)=z_2pp*signo
                            case default
                               err_symm=.true.
                               return
                        end select
                     case (6)
                        err_symm=.true.
                        return
                 end select
                 gener(1:3,4,ngen)=vtras

              case (3)
                 select case (neje)
                     case (1)
                        gener(1:3,1:3,ngen)=x_3*signo
                     case (2)
                        gener(1:3,1:3,ngen)=y_3*signo
                     case (3)
                        gener(1:3,1:3,ngen)=z_3*signo
                     case (4:5)
                        err_symm=.true.
                        return
                     case (6)
                        gener(1:3,1:3,ngen)=xyz_3*signo
                 end select
                 gener(1:3,4,ngen)=vtras

              case (4)
                 select case (neje)
                     case (1)
                        gener(1:3,1:3,ngen)=x_4*signo
                     case (2)
                        gener(1:3,1:3,ngen)=y_4*signo
                     case (3)
                        gener(1:3,1:3,ngen)=z_4*signo
                     case (4:6)
                        err_symm=.true.
                        return
                 end select
                 gener(1:3,4,ngen)=vtras

              case (5)
                 select case (neje)
                     case (1)
                        gener(1:3,1:3,ngen)=x_6*signo
                     case (2)
                        gener(1:3,1:3,ngen)=y_6*signo
                     case (3)
                        gener(1:3,1:3,ngen)=z_6*signo
                     case (4:6)
                        err_symm=.true.
                        return
                 end select
                 gener(1:3,4,ngen)=vtras

          end select

          pos_ini=fin+2
       end do

       !---- Tabla de caracteres ----!
       tabla=0
       tabla(:,:,1)=identidad

       !---- Put Generators on the table ----!
       nt=1
       do i=1,ngen
          if (ngen == 1 .and. num_rot(1)==1) exit    ! only triclinic
          nt=nt+1
          tabla(:,:,nt)=gener(:,:,i)
          tabla(:,4,nt)=mod(tabla(:,4,nt)+48,12)
       end do

       !num_g=nt-1 !Minimum number of generators
        num_g=ngen

            !---- Generate power operations from generators ----!
       do i=2,nt
           ntp=axes_rotation(tabla(:,:,i))    ! Determine the order of the generator
          if (ntp == -1 .or. ntp == -3) then
             ntp=-2*ntp
          else
             ntp=abs(ntp)
          end if

          m1=tabla(:,:,i)
          m2=identidad

          p1:do j=1,ntp-1
             m2=matmul(m2,m1)
             m2(:,4)=mod(m2(:,4)+48,12)

             !---- Check if the generated operation is already in the table
             do k=1,nt
                if (equal_matrix(tabla(:,:,k),m2,4)) cycle p1
             end do

             !---- Eliminating lattice contribution if necessary ----!
             select case (ibravl)
                 case (2)
                    ty=m2(2,4)
                    tz=m2(3,4)
                    if (ty >= 6 .and. tz >= 6) then
                       ty=mod(m2(2,4),6)
                       tz=mod(m2(3,4),6)

                       if (ty == 0 .and. tz == 0) then
                          cycle p1
                       else
                          m2(2,4)=ty
                          m2(3,4)=tz
                       end if
                    end if

                 case (3)
                    tx=m2(1,4)
                    tz=m2(3,4)
                    if (tx >= 6 .and. tz >=6) then
                       tx=mod(m2(1,4),6)
                       tz=mod(m2(3,4),6)

                       if (tx == 0 .and. tz == 0) then
                          cycle p1
                       else
                          m2(1,4)=tx
                          m2(3,4)=tz
                       end if
                    end if

                 case (4)
                    tx=m2(1,4)
                    ty=m2(2,4)
                    if (tx >=6 .and. ty >= 6) then
                       tx=mod(m2(1,4),6)
                       ty=mod(m2(2,4),6)

                       if (tx == 0 .and. ty == 0) then
                          cycle p1
                       else
                          m2(1,4)=tx
                          m2(2,4)=ty
                       end if
                    end if

                 case (5)
                    tx=m2(1,4)
                    ty=m2(2,4)
                    tz=m2(3,4)
                    if (tx >=6 .and. ty >=6 .and. tz >=6) then
                       tx=mod(m2(1,4),6)
                       ty=mod(m2(2,4),6)
                       tz=mod(m2(3,4),6)

                       if (tx == 0 .and. ty == 0 .and. tz == 0) then
                          cycle p1
                       else
                          m2(1,4)=tx
                          m2(2,4)=ty
                          m2(3,4)=tz
                       end if
                    end if

                 case (6)
                    tx=m2(1,4)
                    ty=m2(2,4)
                    tz=m2(3,4)
                    if (tx >=8 .and. ty >=4 .and. tz >=4) then
                       tx=mod(m2(1,4),8)
                       ty=mod(m2(2,4),4)
                       tz=mod(m2(3,4),4)
                       if (tx == 0 .and. ty == 0 .and. tz == 0) then
                          cycle p1
                       else
                          m2(1,4)=tx
                          m2(2,4)=ty
                          m2(3,4)=tz
                       end if

                    else if (tx >=4 .and. ty >=8 .and. tz >=8) then
                       tx=mod(m2(1,4),4)
                       ty=mod(m2(2,4),8)
                       tz=mod(m2(3,4),8)
                       if (tx == 0 .and. ty == 0 .and. tz == 0) then
                          cycle p1
                       else
                          m2(1,4)=tx
                          m2(2,4)=ty
                          m2(3,4)=tz
                       end if
                    end if

                 case (7)
                    tx=m2(1,4)
                    ty=m2(2,4)
                    tz=m2(3,4)
                    if (ty >= 6 .and. tz >=6) then
                       ty=mod(m2(2,4),6)
                       tz=mod(m2(3,4),6)
                       if (tx == 0 .and. ty == 0 .and. tz == 0) then
                          cycle p1
                       else
                          m2(2,4)=ty
                          m2(3,4)=tz
                       end if
                    else if (tx >=6 .and. tz >=6) then
                       tx=mod(m2(1,4),6)
                       tz=mod(m2(3,4),6)
                       if (tx == 0 .and. ty == 0 .and. tz == 0) then
                          cycle p1
                       else
                          m2(1,4)=tx
                          m2(3,4)=tz
                       end if
                    else if (tx >=6 .and. ty >=6) then
                       tx=mod(m2(1,4),6)
                       ty=mod(m2(2,4),6)
                       if (tx == 0 .and. ty == 0 .and. tz == 0) then
                          cycle p1
                       else
                          m2(1,4)=tx
                          m2(2,4)=ty
                       end if
                    end if
             end select

             nt=nt+1
             if (nt > num_tab) then
                err_symm=.true.
                ERR_Symm_Mess=" Dimension of Table exceeded (I)"
                return
             end if
             tabla(:,:,nt)=m2
          end do p1
       end do

       !---- Multiplications between generators ----!
       do
          if (nt == 1) exit
          n=nt

          do i=1,n
             p2:do j=i,n

                m2=matmul(tabla(:,:,i),tabla(:,:,j))
                m2(:,4)=mod(m2(:,4)+48,12)

                !---- Eliminating lattice contribution if necessary ----!
                select case (ibravl)
                   case (2)
                      ty=m2(2,4)
                      tz=m2(3,4)
                      if (ty >= 6 .and. tz >= 6) then
                         ty=mod(m2(2,4),6)
                         tz=mod(m2(3,4),6)
                         m2(2,4)=ty
                         m2(3,4)=tz
                      end if

                   case (3)
                      tx=m2(1,4)
                      tz=m2(3,4)
                      if (tx >= 6 .and. tz >= 6) then
                         tx=mod(m2(1,4),6)
                         tz=mod(m2(3,4),6)
                         m2(1,4)=tx
                         m2(3,4)=tz
                      end if

                   case (4)
                      tx=m2(1,4)
                      ty=m2(2,4)
                      if (tx >= 6 .and. ty >= 6) then
                         tx=mod(m2(1,4),6)
                         ty=mod(m2(2,4),6)
                         m2(1,4)=tx
                         m2(2,4)=ty
                      end if

                   case (5)
                      tx=m2(1,4)
                      ty=m2(2,4)
                      tz=m2(3,4)
                      if (tx >= 6 .and. ty >= 6 .and. tz >= 6) then
                         tx=mod(m2(1,4),6)
                         ty=mod(m2(2,4),6)
                         tz=mod(m2(3,4),6)
                         m2(1,4)=tx
                         m2(2,4)=ty
                         m2(3,4)=ty
                      end if

                   case (6)
                      tx=m2(1,4)
                      ty=m2(2,4)
                      tz=m2(3,4)
                      if (tx >=8 .and. ty >=4 .and. tz >=4) then
                         tx=mod(m2(1,4),8)
                         ty=mod(m2(2,4),4)
                         tz=mod(m2(3,4),4)
                         m2(1,4)=tx
                         m2(2,4)=ty
                         m2(3,4)=tz
                      else if (tx >=4 .and. ty >=8 .and. tz >=8) then
                         tx=mod(m2(1,4),4)
                         ty=mod(m2(2,4),8)
                         tz=mod(m2(3,4),8)
                         m2(1,4)=tx
                         m2(2,4)=ty
                         m2(3,4)=tz
                      end if

                   case (7)
                      tx=m2(1,4)
                      ty=m2(2,4)
                      tz=m2(3,4)
                      if (ty >= 6 .and. tz >=6) then
                         ty=mod(m2(2,4),6)
                         tz=mod(m2(3,4),6)
                         m2(2,4)=ty
                         m2(3,4)=tz
                      else if (tx >=6 .and. tz >=6) then
                         tx=mod(m2(1,4),6)
                         tz=mod(m2(3,4),6)
                         m2(1,4)=tx
                         m2(3,4)=tz
                      else if (tx >=6 .and. ty >=6) then
                         tx=mod(m2(1,4),6)
                         ty=mod(m2(2,4),6)
                         m2(1,4)=tx
                         m2(2,4)=ty
                      end if
                end select

                do k=1,nt
                   if ( equal_matrix(m2(:,:),tabla(:,:,k),4) ) cycle p2
                   if ( equal_matrix(m2(:,:),tabla(:,:,k),3) ) then
                      tx=m2(1,4)+tabla(1,4,k)
                      ty=m2(2,4)+tabla(2,4,k)
                      tz=m2(3,4)+tabla(3,4,k)
                      tx=mod(tx,12)
                      ty=mod(ty,12)
                      tz=mod(tz,12)
                      select case (ibravl)
                          case (2)
                             if (ty == 6 .and. tz == 6) cycle p2
                             if (ty == 0 .and. tz == 0) cycle p2

                          case (3)
                             if (tx == 6 .and. tz == 6) cycle p2
                             if (tx == 0 .and. tz == 0) cycle p2

                          case (4)
                             if (tx == 6 .and. ty == 6) cycle p2
                             if (tx == 0 .and. ty == 0) cycle p2

                          case (5)
                             if (tx == 6 .and. ty == 6 .and. tz == 6) cycle p2
                             if (tx == 0 .and. ty == 0 .and. tz == 0) cycle p2

                          case (6)
                             if (tx == 8 .and. ty == 4 .and. tz == 4) cycle p2
                             if (tx == 4 .and. ty == 8 .and. tz == 8) cycle p2
                             if (tx == 0 .and. ty == 0 .and. tz == 0) cycle p2

                          case (7)
                             if (ty == 6 .and. tz == 6) cycle p2
                             if (tx == 6 .and. tz == 6) cycle p2
                             if (tx == 6 .and. ty == 6) cycle p2

                             if (ty == 0 .and. tz == 0) cycle p2
                             if (tx == 0 .and. tz == 0) cycle p2
                             if (tx == 0 .and. ty == 0) cycle p2

                      end select

                      tx=m2(1,4)-tabla(1,4,k)
                      ty=m2(2,4)-tabla(2,4,k)
                      tz=m2(3,4)-tabla(3,4,k)
                      tx=mod(tx+48,12)
                      ty=mod(ty+48,12)
                      tz=mod(tz+48,12)

                      select case (ibravl)
                          case (2)
                             if (ty == 6 .and. tz == 6) cycle p2
                             if (ty == 0 .and. tz == 0) cycle p2

                          case (3)
                             if (tx == 6 .and. tz == 6) cycle p2
                             if (tx == 0 .and. tz == 0) cycle p2

                          case (4)
                             if (tx == 6 .and. ty == 6) cycle p2
                             if (tx == 0 .and. ty == 0) cycle p2

                          case (5)
                             if (tx == 6 .and. ty == 6 .and. tz == 6) cycle p2
                             if (tx == 0 .and. ty == 0 .and. tz == 0) cycle p2

                          case (6)
                             if (tx == 8 .and. ty == 4 .and. tz == 4) cycle p2
                             if (tx == 4 .and. ty == 8 .and. tz == 8) cycle p2
                             if (tx == 0 .and. ty == 0 .and. tz == 0) cycle p2

                          case (7)
                             if (ty == 6 .and. tz == 6) cycle p2
                             if (tx == 6 .and. tz == 6) cycle p2
                             if (tx == 6 .and. ty == 6) cycle p2

                             if (ty == 0 .and. tz == 0) cycle p2
                             if (tx == 0 .and. tz == 0) cycle p2
                             if (tx == 0 .and. ty == 0) cycle p2
                      end select
                   end if
                end do

                nt=nt+1
                if (nt > num_tab) then
                   err_symm=.true.
                   ERR_Symm_Mess=" Dimension of Table exceeded (II)"
                   return
                end if
                tabla(:,:,nt)=m2

             end do p2
          end do

          if (n == nt) exit

       end do

       !---- Carga Final ----!
       ng=nt
       do i=1,nt
          ss(:,:,i)=tabla(1:3,1:3,i)
          ts(:,i)  = real(tabla(1:3,4,i))/12.0
       end do

       !---- Determination of the crystalline system and Bravais lattice ----!
       call get_crystal_System(ng,ss,isystm,latsy(1:1))
       latsy(2:)=red(ibravl)
       call latsym(red(ibravl))

       return
    End Subroutine Get_SO_from_Hall

    !!----
    !!---- Subroutine Get_So_From_Hms(Isystm,Isymce,Ibravl,Ng,Ss,Ts,Latsy,Spaceh)
    !!----    integer,                        intent(out)  :: ISYSTM    ! Out -> Number of the crystalline system
    !!----                                                                       (1:T, 2:M, 3:O, 4:T, 5:R-Trg, 6:H, 7:C)
    !!----    integer,                        intent(out)  :: ISYMCE    ! Out -> 0 Centric (-1 not at origin)
    !!----                                                                       1 Acentric
    !!----                                                                       2 Centric (-1 at origin)
    !!----    integer,                        intent(out)  :: IBRAVL    ! Out -> Index of the Bravais Lattice type
    !!----                                                                       1   2   3   4   5   6   7
    !!----                                                                      "P","A","B","C","F","I","R"
    !!----    integer,                        intent(out)  :: NG        ! Out -> Number of symmetry operators
    !!----    real(kind=cp),dimension(:,:),   intent(out)  :: TS        ! Out -> Translation parts of the symmetry operators
    !!----    integer, dimension(:,:,:),      intent(out)  :: SS        ! Out -> Rotation parts of the symmetry operators
    !!----    character (len=2),              intent(out)  :: latsy     ! Out -> Bravais lattice symbol
    !!----    character (len=20),             intent( in)  :: SpaceH    !  In -> H-M Spacegroup symbol
    !!----
    !!----    Subroutine to get all the information contained in the H-M symbol.
    !!----    Routine to interpret Hermann-Mauguin symbol for space group.
    !!--..    This routine has been adapted from a program supplied by prof. Burzlaff,
    !!--..    University of Erlangen, Germany.
    !!--..    (Author:Juan Rodriguez-Carvajal)
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Get_SO_from_HMS(Isystm,Isymce,Ibravl,Ng,Ss,Ts,Latsy,SpaceH)
       !---- Arguments ----!
       integer,                   intent(out) :: ISYSTM
       integer,                   intent(out) :: ISYMCE
       integer,                   intent(out) :: IBRAVL
       integer,                   intent(out) :: NG
       integer, dimension(:,:,:), intent(out) :: Ss  !(3,3,48)
       real(kind=cp),    dimension(:,:),   intent(out) :: Ts  !(3,48)
       character (len= 2),        intent(out) :: Latsy
       character (len=*),         intent( in) :: SpaceH

       !---- Local variables ----!
       character (len=20):: GROUP
       character (len=1) :: M,N, Item_SP
       character (len=1), dimension(3,4) :: HMS
       character (len=*), dimension(7), parameter :: IBRA=(/"P","A","B","C","F","I","R"/)
       integer :: SYS,i,j,k,l,NBR,NE,MM,ID,IC,NS,IND,NBL,NLQ,MA
       integer, dimension(3) ::  NMA
       integer, dimension(3,3) :: E
       real(kind=cp) :: TC
       real(kind=cp), dimension(   3):: TE, SH

       NMA = 0
       SH  = 0.25

       !---- Convert to upper case  SpaceG -> GROUP ----!
       group=adjustl(SpaceH)
       call ucase(group)

       NBL=-1
       NG=1
       NS=1
       IBRAVL=0
       TE(1:3) = 0.0
       TS(1:3,1:24) = 0.0
       E(1:3 ,1:3) = 0
       SS(1:3,1:3,1:24) = 0
       SS(1,1,1) = 1
       SS(2,2,1) = 1
       SS(3,3,1) = 1
       HMS(1:3,1:4) = " "
       call init_err_symm()

       do i=1,len_trim(group)         !scanning the Upper-ed case space group symbol
          item_sp=group(i:i)
          if (item_sp == " ") then     !If blank cycle after initializing the indices
             nlq=0                     !for HMS
             ic=0
             cycle
          end if

          if (nbl < 0) then
             do j=1,7
                nbr=j
                if (item_sp == ibra(j)) then
                   nbl=0
                   ibravl=nbr        !Bravais Lattice symbol (first non-blank item of GROUP)
                   exit
                end if
             end do

             if (ibravl == 0) then
                err_symm=.true.
                ERR_Symm_Mess=" Wrong space-group symbol: "//SpaceH
                return
             else if (IBRAVL == 5) then     !These changes are to conform with
                IBRAVL=7                    !the definition of LAT in SYMMETRY
             else if (IBRAVL == 6) then     !modules
                IBRAVL=5
             else if (IBRAVL == 7) THEN
                IBRAVL=6
             end if                         !                1   2   3   4   5   6   7
             cycle                          !               "P","A","B","C","F","I","R"
          end if                            !                p   a   b   c   i   r   f

          if (nlq == 0) nbl=nbl+1      !New blank separating symmetry directions
          ic=ic+1                      !Maximum =4  eg. 63/m
          if (ic > 4) ic=4             !Protection against bad typing of
          if (nbl > 3) nbl=3           !space group symbol
          hms(nbl,ic)=item_sp
          nlq=1
          if (item_sp == "/") ns=0
       end do                          !End loop of Scanning the Upper-ed case space group symbol

       !---- Determination of the crystal system ----!
       SYS=0
       do i=1,4
          if (hms(2,i) == "3") then    !cubic
             sys=6
             isystm=7
             latsy="c"//ibra(NBR)
             exit
          else if (hms(1,i) == "3") then
             sys=5                        !trigonal (rhombohedral)
             isystm=5
             latsy="h"//ibra(NBR)
             exit
          else if (hms(1,i) == "6") then
             sys=5                        !hexagonal (same block as trigonal)
             isystm=6
             latsy="hP"
             exit
          else if(hms(1,i) == "4") then
             sys=4                        !tetragonal
             isystm=4
             latsy="t"//ibra(NBR)
             exit
          end if
       end do

       if (nbl <= 1 .and. sys==0) then
          if (hms(1,1) == "1" .or. hms(1,1) == "-") then
             sys=1           !triclinic
             isystm=1
             latsy="a"//ibra(nbr)
          else
             sys=2           !monoclinic
             isystm=2
             latsy="m"//ibra(nbr)
             do i=1,4
                hms(2,i)=hms(1,i)  ! put the symbol in the form l 1 2/m 1
                hms(1,i)=" "
             end do
             hms(1,1)="1"   !complete the symbol with 1 along a and c
             hms(3,1)="1"
          end if
       end if

       if (nbl > 1 .and. sys==0) then
          sys=3     !orthorhombic
          isystm=3
          latsy="o"//ibra(nbr)
          if (hms(1,1) == "1".or.hms(2,1) == "1") then
             sys=2
             isystm=2
             latsy="m"//ibra(nbr)
          end if
       end if

       call check_symbol_hm(HMS)
       if (err_symm) return
       call latsym(ibra(nbr))


       SELECT CASE (SYS)     !SYS is the crystal family
          CASE (1)      !  TRICLINIC
             IF (HMS(1,1) == "-") NS=0

          CASE (2)      !  MONOCLINIC
             NG=2
             DO I=1,3
                IF (HMS(I,1)/="1") IND=I
             END DO
             ID=1
             IF (HMS(IND,1) == "2") ID=-1
             DO I=1,3
                SS(I,I,2)=SS(I,I,1)*ID
             END DO
             SS(IND,IND,2)=-SS(IND,IND,2)
             DO I=1,3
                IF (HMS(I,1) == "2".AND.HMS(I,2) == "1") TS(I,2)=0.5
                DO J=1,4
                   IF (HMS(I,J) == "A") TS(1,2)=0.5
                   IF (HMS(I,J) == "B") TS(2,2)=0.5
                   IF (HMS(I,J) == "C") TS(3,2)=0.5
                   IF (HMS(I,J) == "N") THEN
                      K=I+1
                      IF (K > 3) K=K-3
                      TS(K    ,2)=0.5
                      TS(6-K-I,2)=0.5
                   END IF
                END DO
             END DO

          CASE (3)   !  ORTHORHOMBIC
             NG=4
             IC=0
             IND=1
             IF (HMS(1,1)/="2".AND.HMS(2,1)/="2".AND.HMS(3,1) /= "2") IND=-1
             IF (IND == -1) NS=0
             DO I=1,3
                ID=1
                IF (HMS(I,1) == "2") ID=-1
                DO J=1,3
                   SS(J,J,1+I)=SS(J,J,1)*ID*IND
                END DO
                SS(I,I,1+I)=-SS(I,I,1+I)
             END DO
             DO I=1,3
                IF (HMS(I,1) == "2" .AND. HMS(I,2) == "1") TS(I,1+I)=0.5
                DO J=1,4
                   IF (HMS(I,J) == "A") TS(1,1+I)=0.5
                   IF (HMS(I,J) == "B") TS(2,1+I)=0.5
                   IF (HMS(I,J) == "C") TS(3,1+I)=0.5
                   IF (HMS(I,J) == "N" .OR. HMS(I,J) == "D") THEN
                      K=I+1
                      IF (K > 3) K=K-3
                      IF (HMS(I,J) == "D") THEN
                         IC=1
                         IF (NS == 1) THEN
                            TS(1,1+I)=0.25
                            TS(2,1+I)=0.25
                            TS(3,1+I)=0.25
                         ELSE                !was missing
                            TS(    K,1+I)=0.25
                            TS(6-K-I,1+I)=0.25
                         END IF
                      ELSE                   !was missing
                        TS(K    ,1+I)=0.5
                        TS(6-K-I,1+I)=0.5
                      END IF
                   END IF
                END DO
             END DO

             if (ic == 1) then
                call mod_trans(ng,ns,ts,isymce)
                return
             end if

             if (ns == 1) then
                ic=0
                do i=1,3
                   if (ss(1,1,1+i)*ss(2,2,1+i)*ss(3,3,1+i) == -1) ic=1  !there are planes
                end do

                if (ic == 1) then
                   do i=1,3
                      if (ss(1,1,1+i)*ss(2,2,1+i)*ss(3,3,1+i) == 1) id=i
                   end do
                   do i=1,3
                      tc=ts(i,2)+ts(i,3)+ts(i,4)
                      if (abs(tc) < eps_symm .or. abs(tc-1.0) < eps_symm) cycle
                      IF (HMS(1,1) == "M" .AND. HMS(2,1) == "N".OR.         &
                          HMS(2,1) == "M" .AND. HMS(3,1) == "N".OR.         &
                          HMS(3,1) == "M" .AND. HMS(1,1) == "N")    THEN
                         k=i-1
                         if (k == 0) k=k+3
                         ts(i,1+k)=0.5
                         cycle !was missing
                      end if
                      do j=1,3
                         if (id/=j) then
                            if (abs(ts(i,1+j)-0.5) > eps_symm) ts(i,1+j)=0.5
                         end if
                      end do
                   end do
                      call mod_trans(ng,ns,ts,isymce)
                      return
                end if   ! it was else

                   tc=ts(1,2)+ts(2,3)+ts(3,4)
                   if (abs(tc) < eps_symm) then
                      call mod_trans(ng,ns,ts,isymce)
                      return
                   end if
                   do i=1,3
                      k=i+1
                      if (k > 3) k=k-3
                      if (tc > 0.5) then
                         if (tc > 1.0) then
                            ts(k,1+i)=0.5
                            cycle
                         end if
                         if (abs(ts(i,1+i)) >=0.000001) cycle
                         l=k+1
                         if (l > 3) l=l-3
                         ts(l,1+k)=0.5
                         ts(k,1+l)=0.5
                         cycle
                      end if
                      if (abs(ts(i,1+i)) < eps_symm) cycle
                      mm=i-1
                      if (mm == 0) mm=mm+3
                      ts(i,1+mm)=0.5
                   end do
                   call mod_trans(ng,ns,ts,isymce)
                   return
             end if

                do i=1,3
                   k=1+i
                   if (k > 3) k=k-3
                   tc=ts(i,1+k)+ts(i,1+6-i-k)
                   if (abs(tc-1.0) < eps_symm) tc=0.0
                   ts(i,1+i)=tc
                end do

                !---- special treatment of c m m a, c m c a, i m m a ---- !
                if (nbr == 1 .or. nbr == 5) then
                   call mod_trans(ng,ns,ts,isymce)
                   return
                end if
                ma=0
                do i=1,3
                   do j=1,4
                      IF (HMS(I,J) == "M") NMA(I)=1
                   end do
                   ma=ma+nma(i)
                end do

                if (.not. (nbr == 6 .and. ma == 2) ) then

                   if (ma == 0 .or. ma == 3 .or. nbr == 6) then
                      call mod_trans(ng,ns,ts,isymce)
                      return
                   end if
                   do i=1,3
                      if (nma(nbr-1) == 1) then
                         call mod_trans(ng,ns,ts,isymce)
                         return
                      end if
                      sh(nbr-1)=0.0
                   end do

                end if

                   !---- origin shift ----!
                   do i=1,ng
                      do j=1,3
                         do k=1,3
                            id=1
                            if (j/=k) id=0
                            ts(j,i)=ts(j,i)+(id-ss(j,k,i))*sh(k)
                         end do
                         if (ts(j,i) > 1.0) ts(j,i)=ts(j,i)-1.0
                         if (ts(j,i) < 0.0) ts(j,i)=ts(j,i)+1.0
                      end do
                   end do
                   call mod_trans(ng,ns,ts,isymce)
                   return

          CASE (4)   !  TETRAGONAL
             NG=4
             IF (NBL == 3) NG=8
             SS(1,2,2)=-1
             SS(2,1,2)=1
             SS(3,3,2)=1
             M=HMS(1,1)
             N=HMS(1,2)
             DO I=1,3
                DO J=1,3
                   IF (M == "-") SS(I,J,2)=-SS(I,J,2)
                END DO
             END DO
             IF (M/="-") THEN
                IF (N == "1") TS(3,2)=0.25
                IF (N == "2") TS(3,2)=0.5
                IF (N == "3") TS(3,2)=0.75
                IF (HMS(1,3) == "N" .OR. (HMS(1,4) =="N" .AND. NBL == 3)) TS(1,2)=0.5
                IF ((HMS(1,4) == "N" .AND. NBL == 1).OR.(N == "1" .AND.     &
                     NS == 1 .AND. NBR == 6)) TS(2,2)=0.5
                IF (N == "1" .AND. NS == 0 .AND. NBR == 6) THEN
                   TS(1,2)=0.25
                   TS(2,2)=0.75
                   IF (NBL == 1) TS(1,2)=0.75
                   IF (NBL == 1) TS(2,2)=0.25
                ELSE IF (HMS(2,2) == "1".OR.(HMS(1,4)/="N" .AND. HMS(2,1)    &
                                 ==   "N" .AND. HMS(3,1) == "M")) THEN
                   TS(1,2)=0.5
                   TS(2,2)=0.5
                END IF
             END IF
             ss(1,1,3)=-1
             ss(2,2,3)=-1
             ss(3,3,3)=1
             ts(1,3)=ss(1,2,2)*ts(2,2)+ts(1,2)
             ts(2,3)=ss(2,1,2)*ts(1,2)+ts(2,2)
             ts(3,3)=ss(3,3,2)*ts(3,2)+ts(3,2)
             do i=1,3
                if (nbr == 6                                 &
                            .and. abs(ts(1,3)-0.5) < eps_symm     &
                            .and. abs(ts(2,3)-0.5) < eps_symm     &
                            .and. abs(ts(3,3)-0.5) < eps_symm) ts(i,3)=0.0
             end do
             do i=1,3
                ts(i,4)=ts(i,2)
                do j=1,3
                   ts(i,4)=ts(i,4)+ss(i,j,2)*ts(j,3)
                   do k=1,3
                      ss(i,j,4)=ss(i,j,4)+ss(i,k,2)*ss(k,j,3)
                   end do
                end do
             end do
             if (nbl == 1) then
                call mod_trans(ng,ns,ts,isymce)
                return
             end if
             m=hms(2,1)
             n=hms(3,1)
             ne=4
             IF (NS == 0) THEN
                E(1,1)=-1
                E(2,2)=1
                E(3,3)=1
                IF (M == "C".OR.M == "N") TE(3)=0.5
                IF (M == "B".OR.M == "N") TE(2)=0.5
                IF (M == "B".OR.M == "N") TE(1)=0.5
                IF (HMS(1,3) == "N".OR.HMS(1,4) == "N") TE(1)=TE(1)+0.5
             ELSE IF (M/="2".AND.N/="2")  THEN
                M=HMS(2,1)
                E(1,1)=-1
                E(2,2)=1
                E(3,3)=1
                IF (M == "C".OR.M == "N") TE(3)=0.5
                IF (M == "N".OR.M == "B") TE(1)=0.5
                IF (M == "N".OR.M == "B") TE(2)=0.5
             ELSE IF (M == "2".AND.N == "2") THEN
                E(1,2)=1
                E(2,1)=1
                E(3,3)=-1
    !            IF (.NOT.(HMS(2,2)/="0" .OR. NBR == 6 .OR. HMS(1,2) == "0")) THEN
                IF (.NOT.(HMS(2,2)/=" " .OR. NBR == 6 .OR. HMS(1,2) == " ")) THEN
                   IF (HMS(1,2) == "1") TE(3)=0.75
                   IF (HMS(1,2) == "2") TE(3)=0.5
                   IF (HMS(1,2) == "3") TE(3)=0.25
                END IF
             ELSE IF (M == "2") THEN
                E(1,1)=1
                E(2,2)=-1
                E(3,3)=-1
                IF (N == "C") TE(3)=0.5
                IF (N == "D") TE(3)=0.25
                IF (N == "D") TE(2)=0.5
                IF (.NOT. (HMS(2,2) /= "1")) THEN
                   TE(1)=0.5
                   TE(2)=0.5
                END IF
             ELSE
                IF (M == "C" .OR. M == "N") TE(3)=0.5
                E(1,1)=-1
                E(2,2)=1
                E(3,3)=1
                IF (.NOT.(M /= "N" .AND. M /= "B") ) THEN
                   TE(1)=0.5
                   TE(2)=0.5
                END IF
             END IF

          CASE(5)   !  HEXAGONAL  and TRIGONAL (RHOMBOHEDRAL)
             NG=3
             NE=6
             IF (HMS(1,1) == "3".OR.(HMS(1,2) == "3".AND.HMS(1,1) == "-")) NE=3
             M=HMS(1,1)
             N=HMS(1,2)
             IF (M == "-".AND.N == "3") NS=0
             IF (M == "6") THEN
                NG=NG+NG
                SS(1,1,2)=1
                SS(1,2,2)=-1
                SS(2,1,2)=1
                SS(3,3,2)=1
                IF (N == "1") TS(3,2)=1.0/6.0
                IF (N == "2") TS(3,2)=2.0/6.0
                IF (N == "3") TS(3,2)=3.0/6.0
                IF (N == "4") TS(3,2)=4.0/6.0
                IF (N == "5") TS(3,2)=5.0/6.0
                DO I=1,4
                   DO J=1,3
                      TS(J,2+I)=TS(J,2)
                      DO K=1,3
                         TS(J,2+I)=TS(J,2+I)+SS(J,K,2)*TS(K,1+I)
                         IF (TS(J,2+I) > 1.0) TS(J,2+I)=TS(J,2+I)-1.0
                         DO L=1,3
                            SS(J,K,2+I)=SS(J,K,2+I)+SS(J,L,2)*SS(L,K,1+I)
                         END DO
                      END DO
                   END DO
                END DO
                IF (NBL == 1) THEN
                   CALL MOD_TRANS(NG,NS,TS,ISYMCE)
                   RETURN
                END IF

             ELSE

                SS(1,2,2)=-1
                SS(2,1,2)=1
                SS(2,2,2)=-1
                SS(3,3,2)=1
                IF (N == "1") TS(3,2)=1.0/3.0
                IF (N == "2") TS(3,2)=2.0/3.0
                SS(1,1,3)=-1
                SS(2,1,3)=-1
                SS(1,2,3)=1
                SS(3,3,3)=1
                TS(3,3)=2.0*TS(3,2)
                IF (TS(3,3) >= 1.0) TS(3,3)=TS(3,3)-1.0
                IF (NBL == 1 .AND. N /= "6") THEN
                   CALL MOD_TRANS(NG,NS,TS,ISYMCE)
                   RETURN
                END IF
                IF (N == "6") THEN
                   NG=NG+NG
                   DO I=1,3
                      DO J=1,3
                         DO K=1,3
                            SS(J,K,3+I)=SS(J,K,I)
                            SS(3,3,3+I)=-1
                         END DO
                      END DO
                   END DO
                END IF
                IF (NBL == 1) THEN
                   CALL MOD_TRANS(NG,NS,TS,ISYMCE)
                   RETURN
                END IF
                IF (.NOT.(HMS(2,1)/="C".AND.HMS(3,1)/="C"))  THEN
                   TS(3,4)=0.5
                   TS(3,5)=0.5
                   TS(3,6)=0.5
                END IF
             END IF

             NG=NG+NG
             M=HMS(2,1)
             N=HMS(3,1)
             IF (M == "1") THEN
                IF (N == "2") THEN
                   E(1,2)=-1
                   E(2,1)=-1
                   E(3,3)=-1
                   TE(3)=2.0*TS(3,2)
                   IF (TE(3) > 1.0) TE(3)=TE(3)-1.0
                ELSE
                   E(1,2)=1
                   E(2,1)=1
                   E(3,3)=1
                   IF (N == "C") TE(3)=0.5
                END IF
             ELSE IF (M == "2") THEN
                E(1,2)=1
                E(2,1)=1
                E(3,3)=-1
                TE(3)=2.0*TS(3,2)
                !---- GROUP P 31 I 2 AND P 32 I 2 ----!
                IF(HMS(1,1)=="3".AND.(HMS(1,2)=="1".OR.HMS(1,2)=="2")) TE(3)=0.0
             ELSE
                E(1,2)=-1
                E(2,1)=-1
                E(3,3)=1
                IF (M == "C") TE(3)=0.5
             END IF

          CASE (6)    !  CUBIC
             NG=12
             IF (NBL == 3) NG=24
             IF (HMS(1,1)/="2".AND.HMS(1,1)/="4".AND.HMS(1,1)/="-") NS=0
             DO I=1,3
                DO J=1,3
                   SS(J,J,1+I)=1
                   IF (I/=J) THEN
                      SS(J,J,1+I)=-1
                      IF (HMS(1,1) == "N") TS(J,1+I)=0.5
                      IF (HMS(1,1) == "D") TS(J,1+I)=0.25
                   END IF
                END DO
             END DO
             IF (.NOT.((HMS(1,1) /="A".AND.HMS(3,1)/="D".AND.HMS(1,2)    &
                                 /="3" .AND.HMS(1,2)/="1".OR.NBR == 5))) THEN
                DO I=1,3
                   TS(I,1+I)=0.5
                   K=I+1
                   IF (K == 4) K=1
                   TS(K,1+I)=0.5
                END DO
             END IF
             DO I=1,4
                DO J=1,3
                   DO K=1,3
                      L=J+1
                      IF (L == 4) L=1
                      MM=J-1
                      IF (MM == 0) MM=3
                      SS(J,K,4+I)=SS(L,K ,I)
                      SS(J,K,8+I)=SS(MM,K,I)
                      TS(J,4+I)=TS(L ,I)
                      TS(J,8+I)=TS(MM,I)
                   END DO
                END DO
             END DO
             IF (NG == 12) THEN
                CALL MOD_TRANS(NG,NS,TS,ISYMCE)
                RETURN
             END IF
             NE=12
             E(1,2)=1
             E(2,1)=1
             E(3,3)=1
             IF (HMS(3,1) == "2") E(3,3)=-1
             IF (HMS(3,1) == "C") TE(3)=0.5
             DO I=1,3
                IF (HMS(3,1)=="N".OR.HMS(1,2)=="2") TE(I)=0.5
                IF (HMS(3,1)=="D".OR.HMS(1,2)=="1".OR.HMS(1,2)=="3") TE(I)=0.25
             END DO
             IF (HMS(1,2) == "1".AND.NBR == 1) TE(1)=0.75
             IF (.NOT.((HMS(1,2) /="1".OR.NBR/=6).AND.(HMS(1,2)   &
                                 /="3" .OR.NBR/=1))) THEN
                TE(2)=0.75
                TE(3)=0.75
             END IF
       END SELECT     ! On crystal system

       if (sys == 4 .or. sys == 5 .or. sys == 6) then
          do i=1,ne
             do j=1,3
                ts(j,ne+i)=te(j)
                do k=1,3
                   ts(j,ne+i)=ts(j,ne+i)+e(j,k)*ts(k,i)
                   do l=1,3
                      ss(j,k,ne+i)=ss(j,k,ne+i)+e(j,l)*ss(l,k,i)
                   end do
                end do
             end do
          end do
       end if

       call mod_trans(ng,ns,ts,isymce)

       return
    End Subroutine Get_SO_from_HMS

    !!----
    !!---- Subroutine Get_Stabilizer(X,Spg,Order,Ptr,Atr)
    !!----    real(kind=cp), dimension(3),  intent (in)  :: x     ! real(kind=cp) space position (fractional coordinates)
    !!----    type(Space_Group_type),       intent (in)  :: Spg   ! Space group
    !!----    integer,                      intent(out)  :: order ! Number of sym.op. keeping invariant the position x
    !!----    integer, dimension(:),        intent(out)  :: ptr   ! Array pointing to the symmetry operators numbers
    !!----                                                        ! of the stabilizer (point group) of x
    !!----    real(kind=cp), dimension(:,:),intent(out)  :: atr   ! Associated additional translation to the symmetry operator
    !!----
    !!----    Subroutine to obtain the list of symmetry operator of a space group that leaves
    !!----    invariant an atomic position. This subroutine provides a pointer to the symmetry
    !!----    operators of the site point group and the additional translation with respect to
    !!----    the canonical representant.
    !!----
    !!---- Update: June - 2011 (JRC)
    !!
    Subroutine Get_Stabilizer(X,Spg,Order,Ptr,Atr)
       !---- Arguments ----!
       real(kind=cp), dimension(3),  intent (in)  :: x     ! real space position (fractional coordinates)
       type(Space_Group_type),       intent (in)  :: Spg   ! Space group
       integer,                      intent(out)  :: order ! Number of sym.op. keeping invariant the position x
       integer, dimension(:),        intent(out)  :: ptr   ! Array pointing to the symmetry operators numbers
                                                           ! of the stabilizer of x
       real(kind=cp), dimension(:,:),intent(out)  :: atr   ! Associated additional translation to the symmetry operator
       !---- Local variables ----!
       real(kind=cp), dimension(3)    :: xx, tr

       integer                        :: j,n1,n2,n3

       order    = 1    !Identity belongs always to the stabilizer
       ptr(:)   = 0
       atr(:,:) = 0.0
       ptr(1)   = 1

       do n1=-1,1
        do n2=-1,1
          do n3=-1,1
            tr=real((/n1,n2,n3/))
             do j=2,Spg%multip
                xx=ApplySO(Spg%SymOp(j),x)+tr-x
                if (sum(abs(xx)) > 2.0 * eps_symm) cycle
                order=order+1
                ptr(order)=j
                atr(:,order)=tr
             end do
          end do
        end do
       end do

       return
    End Subroutine Get_Stabilizer

    !!----
    !!---- Subroutine Get_String_Resolv(T,X,Ix,Symb)
    !!----    real(kind=cp), dimension(3), intent( in) :: t      !  In -> Traslation part
    !!----    real(kind=cp), dimension(3), intent( in) :: x      !  In -> real(kind=cp) part of variable
    !!----    integer, dimension(3),       intent( in) :: ix     !  In -> Frags: 1:x, 2:y, 3:z
    !!----    character (len=*),           intent(out) :: symb   ! Out -> String
    !!----
    !!----    Returning a string for point, axes or plane give as
    !!----    written in fractional form from Resolv_sist procedures.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Get_String_Resolv(t,x,ix,symb)
       !---- Arguments ----!
       real(kind=cp), dimension(3),      intent( in) :: t
       real(kind=cp), dimension(3),      intent( in) :: x
       integer, dimension(3),   intent( in) :: ix
       character (len=*),       intent(out) :: symb

       !---- Local Variables ----!
       character(len=60) :: car
       integer           :: i, np, npos
       real(kind=cp),dimension(3) :: xx

       !---- Main ----!
       xx=x
       do i=1,3
          call get_fraction_2dig(x(i),car)
          np=index(car,"1/2")
          if (np > 0) then
             xx=2.0*x
             exit
          end if
       end do

       symb=" "
       npos=1
       do i=1,3
          !---- Only t value ----!
          if (abs(xx(i)) <= eps_symm) then
             call get_fraction_2dig(t(i),car)
             car=adjustl(car)
             if (car(1:1) == "+") car=car(2:)
             np=len_trim(car)
             if (i < 3) then
                symb(npos:)=car(1:np)//", "
                npos=npos+np+2
             else
                symb(npos:)=car(1:np)
             end if
             cycle
          end if

          call get_fraction_2dig(xx(i),car)
          car=adjustl(car)
          if (abs(abs(xx(i)) - 1.0) <= eps_symm) then
             if (car(1:2) == "+1") car=car(3:)
             if (car(1:2) == "-1") car(2:)=car(3:)
          else
             if (car(1:1) == "+") car=car(2:)
          end if
          np=len_trim(car)
          symb(npos:)=car(1:np)
          npos=npos+np
          select case (ix(i))
             case (1)
                symb(npos:)="x"
             case (2)
                symb(npos:)="y"
             case (3)
                symb(npos:)="z"
          end select
          npos=npos+1
          if (abs(t(i)) > 0.0 ) then
             call get_fraction_2dig(t(i),car)
             car=adjustl(car)
             np=len_trim(car)
             symb(npos:)=car(1:np)
             npos=npos+np
          end if
          if (i < 3) then
             symb(npos:)=", "
             npos=npos+2
          end if

       end do
       symb=pack_string(symb)

       return
    End Subroutine Get_String_Resolv

    !!----
    !!----  Subroutine Get_SubOrbits(X,Spg,ptr,Mult,orb,ind,conv)
    !!----    real(kind=cp), dimension(3),  intent (in) :: x     !  In -> Position vector
    !!----    type(Space_Group_type),       intent (in) :: spgr  !  In -> Space Group
    !!----    integer,dimension(:),         intent( in) :: ptr   !  In -> Pointer to symops of a subgroup
    !!----    integer,                      intent(out) :: mult  !  Out -> Multiplicity
    !!----    real, dimension(:,:),         intent(out) :: orb   !  Out -> List of equivalent positions
    !!----    integer,dimension(:),         intent(out) :: ind   !  Out -> Integer giving the number of the suborbits
    !!----    character(len=*), optional,   intent( in) :: conv  !  In  -> If present centring transl. are considered
    !!----
    !!----    Obtain the multiplicity and list of equivalent positions
    !!----    modulo lattice translations (including centring!) of a
    !!----    position. When symmetry operators of a subgroup of Spg is given
    !!----    an index vector (ind) gives the division in subOrbits.
    !!----    The pointer ptr indicates the symmetry operators of Spg belonging
    !!----    to the subgroup. The first zero value of ptr terminates the search.
    !!----    If the optional argument "conv" is given the centring translations
    !!----    are considered. The orbits are formed by all atoms within a
    !!----    conventional unit cell. Otherwise the orbit is formed only with
    !!----    the content of a primitive cell.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Get_SubOrbits(x,Spg,ptr,mult,orb,ind,conv)
       !---- Arguments ----!
       real(kind=cp), dimension(3),    intent (in) :: x
       type(Space_Group_type),         intent (in) :: spg
       integer,dimension(:),           intent( in) :: ptr
       integer,                        intent(out) :: mult
       real(kind=cp),dimension(:,:),   intent(out) :: orb
       integer,dimension(:),           intent(out) :: ind
       character(len=*), optional,     intent( in) :: conv

       !---- Local variables ----!
       integer                                 :: i,j, nt,is, numorb
       real(kind=cp), dimension(3)             :: xx,v,xi
       character(len=1)                        :: laty

       laty=Spg%spg_lat
       if(present(conv)) laty="P"
       ! First obtain the equivalent positions in the full group
       mult=1
       orb(:,1)=x(:)
       ext: do j=2,Spg%multip
          xx=ApplySO(Spg%SymOp(j),x)
          xx=modulo_lat(xx)
          do nt=1,mult
             v=orb(:,nt)-xx(:)
             if (Lattice_trans(v,laty)) cycle ext
          end do
          mult=mult+1
          orb(:,mult)=xx(:)
       end do ext

       numorb=1
       ind=0
       do i=1,mult
        if(ind(i) /= 0) cycle
        xi=orb(:,i)
        do j=1,Spg%multip
           is= ptr(j)
           if(is == 0) exit
           xx=ApplySO(Spg%SymOp(is),xi)
           xx=modulo_lat(xx)
           do nt=1,mult
              if(ind(nt) /= 0) cycle
              v=orb(:,nt)-xx(:)
              if (Lattice_trans(v,laty)) then
                ind(nt)=numorb
                exit
              end if
           end do
        end do !j
        numorb=numorb+1
       end do !i

       return
    End Subroutine Get_SubOrbits

    !!----
    !!---- Subroutine Get_Symel(Sim,Xyzstring)
    !!----    integer, dimension(3,3), intent( in) :: sim         !  In -> Rotational part
    !!----    character (len=*),       intent(out) :: XYZstring   ! Out -> String
    !!----
    !!----    Supplies a string with the "symmetry element" (I.T.) for the
    !!----    rotation matrix Sim. They correspond to the symbols given in
    !!----    I.T. for space groups Pm3m and P6/mmm.
    !!----    Logical "hexa" must be defined
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Get_SymEl(Sim,Xyzstring)
       !---- Arguments ----!
       integer,dimension (3,3), intent( in) :: sim
       character (len=*),       intent(out) :: XYZstring

       !---- Local Variables ----!
       integer :: Iu,i1,i2,j

       if (.not. hexa) then
          i1=1
          i2=24
       else
          i1=25
          i2=36
       end if
       call SearchOp(sim,i1,i2,Iu)

       if (.not. hexa) then
          j=abs(Iu)
          if (Iu < 0) j=j+24
          XYZstring=IntSymOh(j)
       else
          j=abs(Iu)-24
          if (Iu < 0) j=j+12
          XYZstring=IntSymD6h(j)
       end if

       return
    End Subroutine Get_SymEl

    !!----
    !!---- Subroutine Get_Symkov(Sim,Xyzstring)
    !!----    integer, dimension(3,3), intent( in) :: sim        !  In -> Rotational part
    !!----    character (len=*),       intent(out) :: XYZstring
    !!----
    !!----    Supplies a string with the "symmetry element" (I.T.) for the rotation
    !!----    matrix Sim. They correspond to the symbols Kovalev.
    !!----    Logical "hexa" must be defined
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Get_SymKov(Sim,Xyzstring)
       !---- Arguments ----!
       integer,dimension (3,3), intent( in) :: sim
       character (len=*),       intent(out) :: XYZstring

       !---- Local variables ----!
       integer :: Iu,i1,i2,j

       if (.not. hexa) then
          i1=1
          i2=24
       else
          i1=25
          i2=36
       end if
       call SearchOp(sim,i1,i2,Iu)

       if (.not. hexa) then
          j=abs(Iu)
          if (Iu < 0) j=j+24
          XYZstring=IntSymOh(j)//" -> "//Kov_Oh(j)
       else
          j=abs(Iu)-24
          if (Iu < 0) j=j+12
          XYZstring=IntSymD6h(j)//" -> "//Kov_D6h(j)
       end if

       return
    End Subroutine Get_SymKov

    !!----
    !!---- Subroutine Get_SymSymb(Sim,Tt,Strsym)
    !!----    real(kind=cp)/integer, dimension(3,3), intent( in)    :: sim      !  In -> Rotational part of the S.O.
    !!----    real(kind=cp), dimension( 3),          intent( in)    :: tt       !  In -> Translational part of the S.O.
    !!----    character (len=*),                     intent(out)    :: Strsym   ! Out -> String in th form X,Y,-Z, ...
    !!----
    !!----    Obtain the Jones Faithful representation of a symmetry operator
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine Get_SymsymbI(Sim,Tt,Strsym)
    !!--++    integer, dimension(3,3),      intent( in)    :: sim      !  In -> Rotational part of the S.O.
    !!--++    real(kind=cp), dimension( 3), intent( in)    :: tt       !  In -> Translational part of the S.O.
    !!--++    character (len=*),            intent(out)    :: Strsym   ! Out -> String in th form X,Y,-Z, ...
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Obtain the Jones Faithful representation of a symmetry operator
    !!--++
    !!--++ Update: February - 2005, January-2014 (changed for a more robust algorithm,JRC)
    !!
    Subroutine Get_SymSymbI(X,T,Symb)
       !---- Arguments ----!
       integer,       dimension(3,3), intent( in) :: x
       real(kind=cp), dimension(3),   intent( in) :: t
       character (len=*),          intent(out) :: symb

       !---- Local Variables ----!
       character(len=*),dimension(3),parameter :: xyz=(/"x","y","z"/)
       character(len= 25)              :: car
       character(len= 25),dimension(3) :: sym
       integer           :: i,j

       !---- Main ----!
       symb=" "
       do i=1,3
          sym(i)=" "
          do j=1,3
             if(x(i,j) == 1) then
                sym(i) = trim(sym(i))//"+"//xyz(j)
             else if(x(i,j) == -1) then
                sym(i) =  trim(sym(i))//"-"//xyz(j)
             else if(x(i,j) /= 0) then
               car=" "
               write(unit=car,fmt="(i3,a)") x(i,j),xyz(j)
               if(x(i,j) > 0) car="+"//trim(car)
               sym(i)=trim(sym(i))//pack_string(car)
             end if
          end do
          if (abs(t(i)) > eps_symm ) then
             call get_fraction_2dig(t(i),car)
             sym(i)=trim(sym(i))//trim(car)
          end if
          sym(i)=adjustl(sym(i))
          if(sym(i)(1:1) == "+")  then
            sym(i)(1:1) = " "
            sym(i)=adjustl(sym(i))
          end if
          sym(i)=pack_string(sym(i))
       end do
       symb=trim(sym(1))//","//trim(sym(2))//","//trim(sym(3))
       return
    End Subroutine Get_SymSymbI

    !!--++
    !!--++  Subroutine Get_SymSymbR(X,T,Symb)
    !!--++     real(kind=cp),    dimension(3,3),    intent( in) :: x
    !!--++     real(kind=cp),    dimension(3),      intent( in) :: t
    !!--++     character (len=*),                   intent(out) :: symb
    !!--++
    !!--++     (OVERLOADED)
    !!--++     Returning a string for symmetry operators or for points, axes or plane give as
    !!--++     written in fractional form
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Get_SymSymbR(X,T,Symb)
       !---- Arguments ----!
       real(kind=cp),    dimension(3,3), intent( in) :: x
       real(kind=cp),    dimension(3),   intent( in) :: t
       character (len=*),                intent(out) :: symb

       !---- Local Variables ----!
       character(len= 25):: car
       integer           :: i,j,k, np,npp,npos
       real(kind=cp)     :: suma

       !---- Main ----!
       symb=" "
       npos=1
       do i=1,3
          npp=0
          do j=1,3
             if (abs(x(i,j)) > 0.0 ) then
                call get_fraction_2dig(x(i,j),car)
                car=adjustl(car)
                if (abs(abs(x(i,j))-1.0) <= eps_symm) then
                     if (npp == 0) then
                        select case (car(1:2))
                           case ("-1")
                              car(2:)=car(3:)//"  "
                           case ("+1")
                              car=car(3:)//"  "
                        end select
                     else
                        car(2:)=car(3:)//"  "
                     end if
                else
                   if (npp == 0) then
                      if (car(1:1) =="+") then
                         car=car(2:)//"  "
                      end if
                   end if
                end if

                np=len_trim(car)
                select case (j)
                   case (1)
                      k=index(car(1:np),"/")
                      if( k /= 0) then
                        if(car(k-1:k-1) == "1") then
                          car(k-1:k-1) = "x"
                          symb(npos:)=car(1:np)
                        else
                          symb(npos:)=car(1:k-1)//"x"//car(k:np)
                        end if
                      else
                        symb(npos:)=car(1:np)//"x"
                      end if
                   case (2)
                      k=index(car(1:np),"/")
                      if( k /= 0) then
                        if(car(k-1:k-1) == "1") then
                          car(k-1:k-1) = "y"
                          symb(npos:)=car(1:np)
                        else
                          symb(npos:)=car(1:k-1)//"y"//car(k:np)
                        end if
                      else
                        symb(npos:)=car(1:np)//"y"
                      end if
                   case (3)
                      k=index(car(1:np),"/")
                      if( k /= 0) then
                        if(car(k-1:k-1) == "1") then
                          car(k-1:k-1) = "z"
                          symb(npos:)=car(1:np)
                        else
                          symb(npos:)=car(1:k-1)//"z"//car(k:np)
                        end if
                      else
                        symb(npos:)=car(1:np)//"z"
                      end if
                end select
                npos=len_trim(symb)+1
                npp=npos
             end if
          end do

          if (abs(t(i)) <= eps_symm .and. npp /= 0) then
             if (i < 3) then
                symb(npos:)=", "
                npos=len_trim(symb)+2
             end if
             cycle
          end if

          call get_fraction_2dig(t(i),car)
          car=adjustl(car)
          suma=0.0
          do j=1,3
             suma=suma+abs(x(i,j))
          end do
          np=len_trim(car)
          if (suma <= 3.0*eps_symm) then
             if (car(1:1) == "+") car=car(2:np)//" "
          end if

          if (i < 3) then
             symb(npos:)=car(1:np)//", "
             npos=len_trim(symb)+2
          else
             symb(npos:)=car(1:np)
          end if
       end do

       symb=pack_string(symb)

       return
    End Subroutine Get_SymSymbR

    !!----
    !!---- Subroutine Get_T_SubGroups(SpG,SubG,nsg,point)
    !!----    type (Space_Group_Type) ,             intent( in) :: SpG
    !!----    type (Space_Group_Type) ,dimension(:),intent(out) :: SubG
    !!----    integer,                              intent(out) :: nsg
    !!----    logical, dimension(:,:), optional,    intent(out) :: point
    !!----
    !!----    Subroutine to obtain the list of all non-trivial translationengleiche
    !!----    subgroups (t-subgroups) of a given space group. The unit cell setting
    !!----    is supposed to be the same as that of the input space group "SpG"
    !!----    The search of space sub-groups is performed using a systematic combination
    !!----    of the symmetry operators of the group.
    !!----    The optional argument point has dimensions (SpG%multip,nsg) and contains
    !!----    true point(i,j)=.true. if the operator i of the space group SpG belongs
    !!----    to the subgroup SubG(j).
    !!----
    !!---- Update: February - 2005, April 2015
    !!
    Subroutine Get_T_SubGroups(SpG,SubG,nsg,point)
       !---- Arguments ----!
       type (Space_Group_Type) ,             intent( in) :: SpG
       type (Space_Group_Type) ,dimension(:),intent(out) :: SubG
       integer,                              intent(out) :: nsg
       logical, dimension(:,:), optional,    intent(out) :: point
       !--- Local variables ---!
       integer                            :: i,L,j,k, nc, maxg,ng , nla, i1,i2,nop
       character (len=30), dimension(192) :: gen
       logical                            :: newg, cen_added

       maxg=size(SubG)
       !---- Construct first the generators of centring translations ----!
       ng=0
       nop=SpG%numops !number of symmetry operators excluding lattice centrings
       if (SpG%centred /= 1) nop=nop*2
       do i=2,SpG%numlat
          ng=ng+1
          gen(ng)= SpG%SymopSymb(1+nop*(i-1))
       end do

       nla=ng
       nc=SpG%Numops+1  !Position of the centre of symmetry if it exist
       L=0
       !---- Determine first the triclinic subgroups
       cen_added=.false.
       do
           L=L+1
           newg=.true.
           call set_spacegroup(" ",SubG(L),gen,ng,"gen")
           do j=1,L-1
              if (SpGr_Equal(SubG(L), SubG(j))) then
                 newg=.false.
                 exit
              end if
           end do
           if (newg) then
              call get_HallSymb_from_gener(SubG(L))
           else
              L=L-1
           end if
           if (SpG%centred /= 1 .and. newg .and. .not. cen_added) then !add the centre of symmetry if needed
              ng=ng+1
              gen(ng)=SpG%SymopSymb(nc)
              cen_added=.true.
           else
              exit
           end if
       end do

       !---- Determine first the groups with only one rotational generator
       do i=2,nop
          ng=nla+1
          gen(ng) = SpG%SymopSymb(i)
          cen_added=.false.
          do
             L=L+1
             if (L > maxg) then
                nsg=maxg
                return
             end if
             newg=.true.
             call set_spacegroup(" ",SubG(L),gen,ng,"gen")
             do j=1,L-1
                if (SpGr_Equal(SubG(L), SubG(j))) then
                   newg=.false.
                   exit
                end if
             end do
             if (newg) then
                call get_HallSymb_from_gener(SubG(L))
             else
                L=L-1
             end if
             if (SpG%centred /= 1 .and. newg .and. .not. cen_added) then !add the centre of symmetry if needed
                ng=ng+1
                gen(ng)=SpG%SymopSymb(nc)
                cen_added=.true.
             else
                exit
             end if
          end do
       end do

       !---- Determine now the groups with two rotational generator ----!

       do i1=2,nop-1
          gen(nla+1) = SpG%SymopSymb(i1)
          do i2 = i1+1,nop
             gen(nla+2) = SpG%SymopSymb(i2)
             ng=nla+2
             cen_added=.false.
             do
                L=L+1
                if (L > maxg) then
                   nsg=maxg
                   return
                end if
                newg=.true.
                call set_spacegroup(" ",SubG(L),gen,ng,"gen")
                do j=1,L-1
                   if (SpGr_Equal(SubG(L), SubG(j))) then
                      newg=.false.
                      exit
                   end if
                end do
                if (newg) then
                   call get_HallSymb_from_gener(SubG(L))
                else
                   L=L-1
                end if
                if (SpG%centred /= 1 .and. newg .and. .not. cen_added) then !add the centre of symmetry if needed
                   ng=ng+1
                   gen(ng)=SpG%SymopSymb(nc)
                   cen_added=.true.
                else
                   exit
                end if
             end do
          end do
       end do
       nsg=L
       if(present(point)) then
         point=.false.
         do j=1,nsg
           L=1
           do i=1,SpG%multip
              do k=L,SubG(j)%multip
               if(SubG(j)%SymopSymb(k) == SpG%SymopSymb(i)) then
                  point(i,j) = .true.
                  L=k+1
                  exit
               end if
              end do
           end do
         end do
       end if

       return
    End Subroutine Get_T_SubGroups

    !!----
    !!---- Subroutine Get_Trasfm_Symbol(Mat,tr,abc_symb,oposite)
    !!----    integer, dimension(3,3), intent(in) :: Mat
    !!----    real,    dimension(3),   intent(in) :: tr
    !!----    character(len=*),        intent(out):: abc_symb
    !!----    logical,optional,        intent(in) :: oposite
    !!----
    !!----    Provides the short symbol for a setting change defined by
    !!----    the transfomation matrix Mat and origin given by the translation
    !!----    vector tr. For instance given the matrix:
    !!----
    !!----     1  0 -1                      a'=a-c
    !!----     0  2  0   corresponding to   b'=2b
    !!----     1  0  1                      c'=a+c
    !!----     And the change of origin given by (0.5,0.0,0.5)
    !!----     The subroutine provide the symbol: (1/2,0,1/2; a-c,2b,a+c)
    !!----     If "oposite" is provided theh the symbol: (a-c,2b,a+c; 1/2,0,1/2)
    !!----
    !!----
    !!---- Update: November - 2012, February 2016 (optional argument)
    !!
    Subroutine Get_Trasfm_Symbol(Mat,tr,abc_symb,oposite)
      integer,       dimension(3,3), intent(in) :: Mat
      real(kind=cp), dimension(3),   intent(in) :: tr
      character(len=*),              intent(out):: abc_symb
      logical,optional,              intent(in) :: oposite
      !---- Local variables ----!
      integer :: i
      character(len=25) :: xyz_op, transl
      character(len=6)  :: Fracc
      call Get_SymSymb(Mat,(/0.0_cp,0.0_cp,0.0_cp/),xyz_op)
      do i=1,len_trim(xyz_op)
        if(xyz_op(i:i) == "x")  xyz_op(i:i)="a"
        if(xyz_op(i:i) == "y")  xyz_op(i:i)="b"
        if(xyz_op(i:i) == "z")  xyz_op(i:i)="c"
      end do
      transl=" "
      do i=1,3
        call Get_Fraction_2Dig(tr(i),Fracc)
        transl=trim(transl)//trim(Fracc)//","
      end do
      i=len_trim(transl)
      transl(i:i)=";"
      do i=1,len_trim(transl)-2
        if(transl(i:i) == "+") transl(i:i)=" "
      end do
      transl=Pack_string(transl)
      abc_symb="("//trim(transl)//" "//trim(xyz_op)//")"
      if(present(oposite)) then
        i=len_trim(transl)
        abc_symb="("//trim(xyz_op)//"; "//transl(1:i-1)//")"
      end if
      return
    End Subroutine Get_Trasfm_Symbol

    !!----
    !!---- Subroutine Get_Transl_Symbol(tr,Transl_symb)
    !!----   real,    dimension(3),   intent(in) :: tr
    !!----   character(len=*),        intent(out):: Transl_symb
    !!----
    !!----    Provides the short symbol for a translation vector
    !!----    for which the coordinates are given as fractional symbols
    !!----
    !!---- Update: November - 2012
    !!
    Subroutine Get_Transl_Symbol(tr,Transl_symb)
      real(kind=cp), dimension(3),   intent(in) :: tr
      character(len=*),              intent(out):: Transl_symb
      !---- Local variables ----!
      integer :: i
      character(len=25) :: transl
      character(len=6)  :: Fracc

      transl=" "
      do i=1,3
        call Get_Fraction_2Dig(tr(i),Fracc)
        transl=trim(transl)//trim(Fracc)//","
      end do
      i=len_trim(transl)
      transl(i:i)=" "
      do i=1,len_trim(transl)
        if(transl(i:i) == "+") transl(i:i)=" "
      end do
      Transl_symb="("//trim(transl)//")"
      return
    End Subroutine Get_Transl_Symbol

    !!----
    !!---- Subroutine Init_Err_Symm()
    !!----
    !!----    Initialize the errors flags in this Module
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Init_Err_Symm()

       err_symm=.false.
       ERR_Symm_Mess=" "

       return
    End Subroutine Init_Err_Symm

    !!----
    !!---- Subroutine Inverse_Symm(R,T,S,U)
    !!----    integer, dimension(3,3),     intent(in)  :: R     !  In -> Rotational Part
    !!----    real(kind=cp), dimension(3), intent(in)  :: t     !  In -> Traslational part
    !!----    integer, dimension(3,3),     intent(out) :: S     ! Out -> New Rotational part
    !!----    real(kind=cp), dimension(3), intent(out) :: u     ! Out -> new traslational part
    !!----
    !!----    Calculates the inverse of the symmetry operator (R,t)
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Inverse_Symm(R,t,S,u)
       !---- Arguments ----!
       integer, dimension(3,3),     intent(in)  :: R
       real(kind=cp), dimension(3), intent(in)  :: t
       integer, dimension(3,3),     intent(out) :: S
       real(kind=cp), dimension(3), intent(out) :: u

       !---- Local variables ----!
       integer                        :: ifail
       real(kind=cp), dimension(3,3)  :: a,b

       call init_err_symm()
       a=real(r)
       s=0
       u=0.0

       call matrix_inverse(a,b,ifail)
       if (ifail /= 0) then
          err_symm=.true.
          ERR_Symm_Mess= "Inversion Matrix Failed"
          return
       end if
       s=nint(b)
       u=matmul(-b,t)

       return
    End Subroutine Inverse_Symm

    !!----
    !!---- Subroutine Latsym(Symb,Numl,Latc)
    !!----    character (len=*),                       intent(in)  :: SYMB  !  In -> Space Group H-M/Hall symbol
    !!----    integer, optional,                       intent(in)  :: numL  !  Number of centring vectors
    !!----    real(kind=cp),optional, dimension(:,:),  intent(in)  :: latc  !  Centering vectors
    !!----
    !!--<<        Inlat  Lattice type & associated translations
    !!----          1     P: { 000 }
    !!----          2     A: { 000;  0  1/2 1/2 }+
    !!----          3     B: { 000; 1/2  0  1/2 }+
    !!----          4     C: { 000; 1/2 1/2  0  }+
    !!----          5     I: { 000; 1/2 1/2 1/2 }+
    !!----          6     R: { 000; 2/3 1/3 1/3; 1/3 2/3 2/3   } +
    !!----          7     F: { 000;  0  1/2 1/2; 1/2  0  1/2; 1/2 1/2  0 } +
    !!----          8     Z: { 000;  user-given centring vectors } +
    !!-->>
    !!----    Provides the Lattice type of the S.G. SYMB. Also gives the index (Inlat)
    !!----    of the lattice, the multiplicity (Nlat) and the fractionnal lattice translations
    !!----    ((Ltr(in,j)j=1,3),in=1,Nlat) and Lat_Ch.
    !!----
    !!---- Update: February - 2005, January 2014 (JRC)
    !!
    Subroutine LatSym(SYMB,numL,latc)
       !---- Argument ----!
       character(len=*),                        intent(in)  :: SYMB
       integer, optional,                       intent(in)  :: numL
       real(kind=cp),optional, dimension(:,:),  intent(in)  :: latc  !general vector (JRC, Jan2014)

       !---- Local variables ----!
       character(len=1)                        :: LAT
       character(len=len(symb))                :: SYMBB
       integer                                 :: i

       call init_err_symm()
       symbb=adjustl(symb)
       do i=1,len_trim(symbb)
          if (symbb(i:i) == "-" .or. symbb(i:i) == " ") cycle
          lat=symbb(i:i)
          exit
       end do

       nlat=1
       ltr(:,1)=0.0
       select case (lat)
          case ("P","p")
             lat="P"
             nlat=1
             inlat=1

          case ("A","a")
             lat="A"
             nlat=2
             inlat=2
             ltr(1,2)=0.0
             ltr(2,2)=0.5
             ltr(3,2)=0.5

          case ("B","b")
             lat="B"
             nlat=2
             inlat=3
             ltr(1,2)=0.5
             ltr(2,2)=0.0
             ltr(3,2)=0.5

          case ("C","c")
             lat="C"
             nlat=2
             inlat=4
             ltr(1,2)=0.5
             ltr(2,2)=0.5
             ltr(3,2)=0.0

          case ("I","i")
             lat="I"
             nlat=2
             inlat=5
             ltr(:,2)=0.5

          case ("R","r")
             lat="R"
             nlat=3
             inlat=6
             ltr(1,2)=2.0/3.0
             ltr(2,2)=1.0/3.0
             ltr(3,2)=1.0/3.0
             ltr(1,3)=1.0/3.0
             ltr(2,3)=2.0/3.0
             ltr(3,3)=2.0/3.0

          case ("F","f")
             lat="F"
             nlat=4
             inlat=7
             ltr(1,2)=0.5
             ltr(2,2)=0.5
             ltr(3,2)=0.0
             ltr(1,3)=0.5
             ltr(2,3)=0.0
             ltr(3,3)=0.5
             ltr(1,4)=0.0
             ltr(2,4)=0.5
             ltr(3,4)=0.5

          case ("Z","z","X","x")
             if(present(numL) .and. present(latc)) then
              lat="Z"
              nlat=numL+1
              !nlat=min(nlat,12) !restriction removed in January 2014
              inlat=8
              do i=2,nlat
                ltr(:,i)=latc(:,i-1)
              end do
             else
               err_symm=.true.
               ERR_Symm_Mess="Unconventional Lattice Symbol Z needs centring vectors"
             end if
          case default
             err_symm=.true.
             ERR_Symm_Mess="Wrong Lattice Symbol "//LAT
       end select

       Lat_Ch=LAT

       return
    End Subroutine Latsym

    !!--++
    !!--++ Subroutine Max_Conv_Lattice_Type(L, Latc, Lattyp)
    !!--++    integer,                        intent(in)  :: L         !  number of centring vectors
    !!--++    real(kind=cp), dimension(:,:),  intent(in)  :: Latc      ! (3,11) centring vectors
    !!--++    character(len=*),               intent(out) :: lattyp    ! Lattice symbol
    !!--++
    !!--++    (PRIVATE)
    !!--++    Subroutine to get the maximum conventional lattice symbol from
    !!--++    a set of possible centring vectors.
    !!--++    Used by subroutine: Similar_Transf_SG
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Max_Conv_Lattice_Type(L, Latc, lattyp)
       !---- Arguments ----!
       integer,                        intent( in) :: L
       real(kind=cp), dimension(:,:),  intent( in) :: Latc
       character(len=*),               intent(out) :: lattyp

       !---- Local variables ----!
       logical                            :: latt_p, latt_a, latt_b, latt_c, latt_i, latt_r, latt_f
       integer, dimension(3)              :: tt
       integer                            :: i, j
       integer, dimension(3,6), parameter :: lattice=reshape((/0,6,6, 6,0,6, &
                                                     6,6,0, 6,6,6, 8,4,4, 4,8,8/),(/3,6/))

       if (l == 0) then !primitive lattice
          lattyp="P"
          return
       end if

       latt_p=.true.
       latt_a=.false.
       latt_b=.false.
       latt_c=.false.
       latt_i=.false.
       latt_r=.false.
       latt_f=.false.

       do i=1,L
          tt(1:3)=nint(12.0 * Latc(1:3,i))   ! Translations x 12

          !---- Compare the translation part of the operator with tabulated array ----!
          do j=1,6
             if (equal_vector(tt,lattice(:,j),3)) then
                select case (j)
                   case (1)
                      latt_a=.true.
                   case (2)
                      latt_b=.true.
                   case (3)
                      latt_c=.true.
                   case (4)
                      latt_i=.true.
                   case (5,6)
                      latt_r=.true.
                end select
                exit
             end if
          end do
       end do

       !---- Lattice Type ----!
       if ( (latt_a .and. latt_b .and. latt_c) .or. (latt_a .and. latt_b) .or. &
            (latt_a .and. latt_c) .or. (latt_b .and. latt_c) ) then
            latt_f=.true.
            latt_a=.false.
            latt_b=.false.
            latt_c=.false.
            latt_p=.false.
            latt_i=.false.
       end if
       if (latt_p) lattyp="P"
       if (latt_a) lattyp="A"
       if (latt_b) lattyp="B"
       if (latt_c) lattyp="C"
       if (latt_i) lattyp="I"
       if (latt_r) lattyp="R"
       if (latt_f) lattyp="F"

       return
    End Subroutine Max_Conv_Lattice_Type

    !!--++
    !!--++ Subroutine Mod_Trans(Ng, Ns, Ts, Isymce)
    !!--++    integer, intent( in)                           :: ng      ! In -> Number of operators
    !!--++    integer, intent( in)                           :: ns      ! In ->
    !!--++    real(kind=cp), dimension(3,24), intent(in out) :: ts      ! In -> Traslation part
    !!--++                                                                Out ->
    !!--++    integer, intent(out),optional                  :: isymce  ! Out -> Origin information
    !!--++                                                                0= Ccenter of Inversion in the Origin
    !!--++                                                                1= Non centrosymmetric
    !!--++                                                                2= Center of Inversion out of origin
    !!--++
    !!--++    (PRIVATE)
    !!--++    Subroutine used by Get_SO_from_HMS.
    !!--++    Put all tranlations in conventional form (positive and less than 1)
    !!--++    Provides Isymce
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Mod_Trans(Ng,Ns,Ts,Isymce)
       !---- Arguments ----!
       integer, intent(          in)                  :: ng,ns
       real(kind=cp), dimension(3,24), intent(in out) :: ts
       integer, intent(out),optional                  :: isymce

       !---- Local Variables ----!
       integer :: i

       do i=1,ng
          ts(:,i)=modulo_lat(ts(:,i))
       end do
       if (present(isymce)) isymce=2-ns

       return
    End Subroutine Mod_Trans

    !!----
    !!---- Subroutine Read_Bin_Spacegroup(SpG,Lun,ok)
    !!----    type (Space_Group),  intent(out) :: SpG   !  Out -> SpaceGroup Variable
    !!----    integer,             intent(in)  :: Lun   !  In -> Logical unit of the file
    !!----    logical,             intent(out) :: ok    !  .true. if everything is OK
    !!----
    !!----    Reading in file of logical unit "lun" the full structure of Space_Group_Type, SpG
    !!----    The file should have been opened with the access="stream" attribute. The procedure
    !!----    reads in the given order a series of bytes corresponding to the components of the
    !!----    type SpG.
    !!----
    !!---- Update: February - 2013
    !!
    Subroutine Read_Bin_SpaceGroup(SpG,lun,ok)
       !---- Arguments ----!
       type (Space_Group_Type),intent(out) :: SpG
       integer,                intent(in)  :: lun
       logical,                intent(out) :: ok

       !---- Local variables ----!
       integer                           :: i,j,ier

       ok=.true.
       read(unit=Lun,iostat=ier) SpG%NumSpg,        &   ! Number of the Space Group
                                 SpG%SPG_Symb,      &   ! Hermann-Mauguin Symbol
                                 SpG%Hall,          &   ! Hall symbol
                                 SpG%CrystalSys,    &   ! Crystal system
                                 SpG%Laue,          &   ! Laue Class
                                 SpG%PG,            &   ! Point group
                                 SpG%Info,          &   ! Extra information
                                 SpG%SG_setting,    &   ! Information about the SG setting (IT,KO,ML,ZA,Table,Standard,UnConventional)
                                 SpG%Hexa,          &   !
                                 SpG%SPG_lat,       &   ! Lattice type
                                 SpG%SPG_latsy,     &   ! Lattice type Symbol
                                 SpG%NumLat             ! Number of lattice points in a cell
       if(ier /= 0) then
        ok=.false.
        return
       end if

       if(allocated(SpG%Latt_trans)) deallocate(SpG%Latt_trans)
       allocate(SpG%Latt_trans(3,SpG%NumLat))

       read(unit=Lun,iostat=ier) SpG%Latt_trans,    &   ! Lattice translations
                                 SpG%Bravais,       &   ! String with Bravais symbol + translations
                                 SpG%Centre,        &   ! Alphanumeric information about the center of symmetry
                                 SpG%Centred,       &   ! Centric or Acentric [ =0 Centric(-1 no at origin),=1 Acentric,=2 Centric(-1 at origin)]
                                 SpG%Centre_coord,  &   ! Fractional coordinates of the inversion centre
                                 SpG%NumOps,        &   ! Number of reduced set of S.O.
                                 SpG%Multip,        &   ! Multiplicity of the general position
                                 SpG%Num_gen            ! Minimum number of operators to generate the Group
       if(ier /= 0) then
        ok=.false.
        return
       end if

       if(allocated(SpG%SymOp)) deallocate(SpG%SymOp)
       allocate(SpG%SymOp(SpG%Multip))
       if(allocated(SpG%SymOpSymb)) deallocate(SpG%SymOpSymb)
       allocate(SpG%SymOpSymb(SpG%Multip))

       do i=1,SpG%Multip
         read(unit=Lun,iostat=ier) SpG%SymOp(i)%Rot,SpG%SymOp(i)%tr, & ! Symmetry operators
                                   SpG%SymopSymb(i)                    ! Strings form of symmetry operators
         if(ier /= 0) then
          ok=.false.
          return
         end if
       end do
       read(unit=Lun,iostat=ier) SpG%R_Asym_Unit      ! Asymmetric unit in real space
       if(ier /= 0) then
         ok=.false.
         return
       end if
       read(unit=Lun,iostat=ier) SpG%Wyckoff%num_orbit              ! Wyckoff Information
       if(ier /= 0) then
        ok=.false.
        return
       end if
       if (SpG%Wyckoff%num_orbit == 0) return
       do i=1,SpG%Wyckoff%num_orbit
         read(unit=Lun,iostat=ier) SpG%Wyckoff%orbit(i)%norb
         read(unit=Lun,iostat=ier) SpG%Wyckoff%orbit(i)%str_Orig
         do j=1,SpG%Wyckoff%orbit(i)%norb
           read(unit=Lun,iostat=ier) SpG%Wyckoff%orbit(i)%str_orbit(j)
         end do
         if(ier /= 0) then
          ok=.false.
          return
         end if
       end do
       return
    End Subroutine Read_Bin_SpaceGroup

    !!----
    !!---- Subroutine Read_Msymm(Info,Sim,P_Mag,ctrl)
    !!----    character (len=*),       intent( in) :: Info   !  In -> Input string with S.Op.
    !!----                                                            in the form: MSYM  u,w,w,p_mag
    !!----    integer, dimension(3,3), intent(out) :: sim    ! Out -> Rotation matrix
    !!----    real(kind=cp),           intent(out) :: p_mag  ! Out -> magnetic phase
    !!----    logical, optional,       intent(in)  :: ctrl   ! in  -> If provided and .true. an error condition
    !!----                                                            is raised if the det(Sim)=0
    !!----    Read magnetic symmetry operators in the form U,V,W, etc...
    !!----    Provides the magnetic rotational matrix and phase associated to a MSYM symbol
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Msymm(Info,Sim,P_Mag,ctrl)
       !---- Arguments ----!
       character (len=*),       intent( in) :: Info
       integer, dimension(3,3), intent(out) :: sim
       real(kind=cp),           intent(out) :: p_mag
       logical, optional,       intent(in)  :: ctrl

       !---- Local variables ----!
       integer ::  i,imax,nop,s,ifound,j,ioerr,istart,mod_istart
       character(len=len(info)) :: aux
       logical :: control

       control=.false.
       if(present(ctrl)) control=ctrl
       call init_err_symm()
       do j=len(Info),1,-1
          if (info(j:j) == ",") exit
       end do
       p_mag=0.0
       imax=j-1
       read(unit=info(j+1:),fmt=*,iostat=ioerr) p_mag
       if (ioerr /= 0) then
          p_mag=0.0
       end if
       sim = 0
       aux=adjustl(l_case(Info))
       if(aux(1:4) == "msym" .or. aux(1:4) == "dsym") then
         istart=6
       else
         istart=1
       end if

       do nop=1,3
          s=1
          mod_istart=0
          ifound=0
          do i=istart,imax
             if (aux(i:i) == " ") cycle
             if (aux(i:i) == "," .or. info(i:i) == "*") then
                mod_istart=1
                exit
             end if
             ifound=1
             if (aux(i:i) == "u" ) then
                sim(nop,1)=s
                s=1
             else if (aux(i:i) == "v") then
                sim(nop,2)=s
                s=1
             else if(aux(i:i) == "w") then
                sim(nop,3)=s
                s=1
             else if(aux(i:i) == "+") then
                s=1
             else if(aux(i:i) == "-") then
                s=-1
             else
                err_symm=.true.
                ERR_Symm_Mess=" Invalid character... "//aux(I:I)//" in Sym. Op."
                return
             end if
          end do    !End loop through the string

          if (mod_istart == 1) then
            istart=i+1
          end if

          if (ifound == 0) then
             err_symm=.true.
             ERR_Symm_Mess=" Blank operator field "//info
             return
          end if
       end do    !End external loop over the three expected items

       if (determ_A(sim) == 0 .and. control) then      !Verify it is a suitable s.o.
          err_symm=.true.
          ERR_Symm_Mess=" The above operator is wrong "//info
          return
       end if

       if (ifound == 1) return

       err_symm=.true.
       ERR_Symm_Mess=" The above operator is wrong "//info

       return
    End Subroutine Read_Msymm

    !!----
    !!---- Subroutine Read_SymTrans_Code(Code,N,Tr)
    !!----    character (len=*),          intent( in) :: Code
    !!----    integer,                    intent(out) :: N
    !!----    real(kind=cp),dimension(3), intent(out) :: Tr
    !!----
    !!----    Read a Code string for reference the symmetry operator and the
    !!----    Traslation applied.
    !!--<<        _2.555     : N_Op = 2, Tr=( 0.0, 0.0, 0.0)
    !!----        _3.456     : N_Op = 3, Tr=(-1.0, 0.0, 1.0)
    !!-->>
    !!----
    !!---- Update: April - 2005
    !!
    Subroutine Read_SymTrans_Code(Code,N,Tr)
       !---- Arguments ----!
       character (len=*),          intent( in) :: Code
       integer,                    intent(out) :: N
       real(kind=cp),dimension(3), intent(out) :: Tr

       !---- Local variables ----!
       character(len=20) :: car
       integer          :: i,j,k,n_ini,n_end,nt

       N=1
       Tr=0.0
       if (len_trim(code) <= 0) return

       car=adjustl(code)
       n_ini=index(car,"_")
       n_ini=n_ini+1

       !---- Found Number of Symmetry Operator ----!
       n_end=index(car,".")
       if (n_end ==0) n_end=len_trim(car)+1
       read (unit=car(n_ini:n_end-1),fmt=*) n

       !---- Found the Traslation ----!
       n_ini=index(car,".")
       if (n_ini /= 0) then
          n_ini=n_ini+1
          n_end=len_trim(car)
          read (unit=car(n_ini:n_end),fmt=*) nt
          i=nt/100
          j=mod(nt,100)/10
          k=nt-(i*100+j*10)
          i=i-5
          j=j-5
          k=k-5
          tr(1)=real(i)
          tr(2)=real(j)
          tr(3)=real(k)
       end if

       return
    End Subroutine Read_SymTrans_Code

    !!----
    !!---- Subroutine Read_Xsym(Info,Istart,Sim,Tt,ctrl)
    !!----    character (len=*),                     intent( in)    :: Info   !  In -> String with the symmetry symbol
    !!----                                                                             in the form: SYMM  x,-y+1/2,z
    !!----    integer,                               intent(in)     :: istart !  In -> Starting index of info to read in.
    !!----    integer, dimension(3,3),               intent(out)    :: sim    ! Out -> Rotational part of S.O.
    !!----    real(kind=cp), optional, dimension(3), intent(out)    :: tt     ! Out -> Traslational part of S.O.
    !!----
    !!----
    !!----    Read symmetry or transformation operators in the form X,Y,Z, etc...
    !!----    Provides the rotational matrix and translation associated a to SYMM symbol
    !!----    in the Jones Faithful representation.
    !!----
    !!---- Update: June - 2011 (JRC, adding ctrl for controlling if a real symmetry operator is needed)
    !!
    Subroutine Read_Xsym(Info,Istart,Sim,Tt,ctrl)
       !---- Arguments ----!
       character (len=*),                     intent(in)     :: Info
       integer,                               intent(in)     :: istart
       integer, dimension(3,3),               intent(out)    :: sim
       real(kind=cp), optional, dimension(3), intent(out)    :: tt
       logical,       optional,               intent(in)     :: ctrl

       !---- Local variables ----!
       character (len=*), dimension(10), parameter :: ANUM=(/"1","2","3","4","5","6","7","8","9","0"/)
       integer, dimension(10), parameter           :: NUM =(/1,2,3,4,5,6,7,8,9,0/)
       integer :: i,imax,nop,s,np,isl,ifound,ip,k,mod_istart,ST=0,I_P,ist
       real(kind=cp) :: t,a
       logical       :: control

       control=.true.
       if(present(ctrl)) control=ctrl
       call init_err_symm()
       imax=len_trim(info)
       if (present(tt)) tt=0.0
       sim = 0
       ist=istart
       do nop=1,3
          s=1
          t=0.0
          ip=0
          i_p=1
          np=0
          isl=0
          ifound=0
          mod_istart=0
          loop_string: do i=ist,imax
             if (info(i:i) == " ") cycle
             if (info(i:i) == "," .or. info(i:i) == "*") then
                mod_istart=1
                exit
             end if
             ifound=1
             if (info(i:i) == "X" .or. info(i:i) == "x") then
                sim(nop,1)=s*i_p
                i_p=1
                s=1
             else if (info(i:i) == "Y" .or. info(i:i) == "y") then
                sim(nop,2)=s*i_p
                i_p=1
                s=1
             else if(info(i:i) == "Z" .or. info(i:i) == "z") then
                sim(nop,3)=s*i_p
                i_p=1
                s=1
             else if(info(i:i) == "+") then
                s=1
             else if(info(i:i) == "-") then
                s=-1
             else if(info(i:i) == "/") then
                isl=1
             else if(info(i:i) == ".") then
                ip=1
             else
                st=s
                do k=1,10
                   if (info(i:i) == anum(k))  then
                      if (is_xyz(info(i+1:i+1))) then
                         i_p=num(k)
                         cycle loop_string
                      else
                         a=num(k)
                         if (isl == 1) then
                            t=t/a
                         else if(ip == 1) then
                            np=np+1
                            t=t+a/10**np
                         else
                            t=10.0*t+a
                         end if
                         cycle loop_string
                      end if
                   end if
                end do
                err_symm=.true.
                ERR_Symm_Mess=" Invalid character... "//INFO(I:I)//" in operator string"
                return
             end if
          end do  loop_string   !end loop through the string (index:i= ist,imax)

          if (mod_istart == 1) then
             ist=i+1
          end if

          t=t*st
          if (present(tt)) tt(nop)=t

          if (ifound == 0) then
             err_symm=.true.
             ERR_Symm_Mess=" Blank operator field"
             return
          end if

       end do    !End external loop over the three expected items (index:NOP)

       if (determ_A(sim) == 0 .and. control) then      !Verify it is a suitable s.o.
          err_symm=.true.
          ERR_Symm_Mess=" The above operator is wrong: "//info
          return
       end if

       if (ifound == 1) return

       err_symm=.true.
       ERR_Symm_Mess=" The above operator is wrong: "//info

       return
    End Subroutine Read_Xsym

    !!----
    !!---- Subroutine Searchop(Sim,I1,I2,Isl)
    !!----    integer , dimension(3,3), Intent(in)  :: sim      !  In -> Rotational part of a symmetry operator
    !!----    integer ,                 Intent(in)  :: i1       !  In -> i1=1,  i2=24  if not hexagonal  (matrices of m3m )
    !!----    integer ,                 Intent(in)  :: i2       !  In -> i1=25, i2=36  if     hexagonal  (matrices of 6/mmm)
    !!----    integer ,                 Intent(out) :: Isl      ! Out -> Index of the matrix Mod6(Isl,:,:)=sim.
    !!----                                                               This index allow to get the corresponding tabulated symmetry symbol.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Searchop(Sim,I1,I2,Isl)
       !---- Arguments ----!
       integer , dimension(3,3), Intent(in) :: sim
       integer , Intent(in)                 :: i1,i2
       integer , Intent(out)                :: Isl

       !---- Local variables ----!
       integer               :: iss,ipass,j,k,im

       iss=1
       ipass=0
       call init_err_symm()
       do
          ipass=ipass+1
          imdo:  do im=i1,i2
             Isl=0
             do j=1,3
                do k=1,3
                   if (sim(j,k) /= iss*Mod6(im,j,k)) cycle imdo
                end do
             end do
             Isl=iss*im
             exit
          end do imdo

          if (Isl /= 0) return

          if (ipass >=2 ) then
             ERR_Symm_Mess=" Try to re-write your S.O. using a rotational part"
             if (i1 == 1 .and.  i2 == 24) then
                ERR_Symm_Mess=trim(ERR_Symm_Mess)//" identical to a S.O. of the space group P m -3 m"
             else if(i1 == 25 .and.  i2 == 36) then
                ERR_Symm_Mess=trim(ERR_Symm_Mess)//" identical to a S.O. of the space group P 6/m m m"
             else
                ERR_Symm_Mess=trim(ERR_Symm_Mess)//" identical to a S.O. of the space group P m -3 m or P 6/m m m"
             end if
             err_symm=.true.
             return
          end if
          iss=-1
       end do

       return
    End Subroutine Searchop

    !!----
    !!---- Subroutine Set_Spacegroup(Spacegen, Spacegroup, Gen, Ngen, Mode, Force_Hall)
    !!----    character (len=*),                       intent(in)     :: SpaceGen     !  In -> String with Number, Hall or Hermman-Mauguin
    !!----    Type (Space_Group),                         intent(out) :: SpaceGroup   ! Out -> SpaceGroup variable
    !!----    character (len=*), dimension(:),  intent(in ), optional :: gen          !  In -> String Generators
    !!----    Integer,                          intent(in ), optional :: ngen         !  In -> Number of Generators
    !!----    character (len=*),                intent(in ), optional :: Mode         !  In -> HMS, ITC, Hall, Gen, Fix
    !!----    character (len=*),                intent(in ), optional :: force_hall   !  In -> f_hall (if present force generation from Hall)
    !!----
    !!----    Subroutine that construct the object SpaceGroup from the H-M or Hall symbol.
    !!----    Expand the set of operators including centre of symmetry and non integer
    !!----    translations for centred cells.
    !!----    If the optional argument Gen is given, then Ngen and Mode="GEN" should be given.
    !!----    If the optional argument mode="ITC", the space group will be generated using the
    !!----    the generators given in the International Tables for the standard setting. In this
    !!----    case the string in SpaceGen should correspond to the Hermann-Mauguin symbol.
    !!----    If the optional argument mode="HMS","HALL" is given the string in SpaceGen
    !!----    should correspond to the desired symbol.
    !!----    If Gen,NGen and Mode are not given but force_hall="f_hall" is given, the generation
    !!----    of the symmetry operators from the symbol of the space group is according to the Hall
    !!----    symbol even if the provided symbol is of Hermann-Maugin type.
    !!----    The use of the different options give rise to different ordering of the symmetry
    !!----    operators or different origins and settings for the same space group.
    !!----
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Set_SpaceGroup(Spacegen,Spacegroup,Gen,Ngen,Mode,Force_Hall)
       !----Arguments ----!
       character (len=*),                intent(in )           :: SpaceGen
       Type (Space_Group_Type),          intent(out)           :: SpaceGroup
       character (len=*), dimension(:),  intent(in ), optional :: gen
       Integer,                          intent(in ), optional :: ngen
       character (len=*),                intent(in ), optional :: Mode
       character (len=*),                intent(in ), optional :: force_hall

       !---- Local variables ----!
       character (len=*),dimension(0:2), parameter  :: Centro = &
                                          (/"Centric (-1 not at origin)", &
                                            "Acentric                  ", &
                                            "Centric (-1 at origin)    "/)
       character (len=20)               :: Spgm
       character (len=20)               :: ssymb
       character (len=130)              :: gener
       character (len=3)                :: opcion
       character (len=2)                :: Latsy
       integer                          :: num, i, j, iv, istart
       integer,      dimension(1)       :: ivet
       integer,      dimension(5)       :: poscol
       integer                          :: isystm,isymce,ibravl,Num_g
       integer                          :: m,l,ngm,k,ier
       integer                          :: ng
       integer,      dimension(3,3,384) :: ss
       real(kind=cp),dimension(3,384)   :: ts
       real(kind=cp),dimension(3)       :: co
       real(kind=cp),dimension(1)       :: vet
       real(kind=cp),dimension(3)       :: vec
       logical                          :: ok

       !---- Inicializing Space Group ----!
       call init_err_symm()

       !Constructor eliminated because some components are nowadays allocatable
       !There's no risk of undefined fields because in this procedure everything is set.
       !
       !SpaceGroup=Space_Group_Type(0,"unknown","unknown","unknown","?","?","?","?",.false.,"?","?", &
       !                       0,0.0,"?","?", -1, 0.0,  0, 0,  0, Sym_Oper_Type(0, 0.0),"?",         &
       !                       wyckoff_type(0,wyck_pos_type(0," ",0," "," ")),0.0)
       SpaceGroup%R_Asym_Unit(:,2)=1.0
       SpaceGroup%gHall=" "

       !---- Loading Tables ----!
       call Set_Spgr_Info()
       call Set_Wyckoff_Info()

       !---- Mode Option ----!
       opcion=" "
       spgm=adjustl(SpaceGen)
       spgm=u_case(spgm)
       num=-1

       if (present(mode)) then
          opcion=adjustl(mode)
          call ucase(opcion)

          Select case (opcion(1:3))

            case("HMS")
                do i=1,num_spgr_info
                   if (spgm(1:12) == spgr_info(i)%hm) then
                      num=i
                      exit
                   end if
                end do

            case("ITC")
                read(unit=spgm(1:12),fmt=*,iostat=ier) ivet(1)

                if( ier == 0) then
                   do i=1,num_spgr_info
                     if (ivet(1) == spgr_info(i)%n) then
                        num=i
                      exit
                     end if
                   end do
                else
                   do i=1,num_spgr_info
                     if (spgm(1:12) == spgr_info(i)%hm) then
                       num=i
                       exit
                     end if
                   end do
                end if

            case("HAL")
                do i=1,num_spgr_info
                   if (spgm(1:16) == u_case(spgr_info(i)%hall)) then
                      num=i
                      exit
                   end if
                end do
          End select

       else       ! detect automatically the symbol of the group

          call getnum(spgm,vet,ivet,iv)
          if (iv /= 1) then
             !---- Is HM Symbol ? ----!
             do i=1,num_spgr_info
                if (spgm(1:12) == spgr_info(i)%hm) then
                   num=i
                   opcion="HMS"
                   if(present(force_hall)) then
                     opcion="HAL"
                     spgm=spgr_info(i)%hall
                   end if
                   exit
                end if
             end do
              !Special treatment of groups P N M 21 (211), P 21 N M (213),
              !P M 21 N (215), B A M B (375), C 4 2 21 (429), C -4 B 2 (457)
              ! and F -4 D 2 (463)   => Force HALL in all these cases
              if(opcion(1:3) == "HMS" .and.  &
              (num==211 .or. num==213 .or. num==215  &
              .or. num==375 .or. num==429 .or. num==457 .or. num==463)) then
                  opcion(1:3) = "HAL"
                  spgm=spgr_info(num)%hall
              end if
             !---- Is a standard Hall Symbol ? ----!
             if (num < 0) then
                do i=1,num_spgr_info
                   if (spgm(1:16) == u_case(spgr_info(i)%hall)) then
                      num=i
                      opcion="HAL"
                      exit
                   end if
                end do
             end if

             !---- Using Generators ----!
             if (num <=0) then
                 if(present(gen)) then
                   opcion="GEN"
                 else  !The last option is a non-standard Hall symbol
                    opcion="HAL"
                 end if
             end if
          else
             if (ivet(1) > 0 .and. ivet(1) < 231) then
               do i=1,num_spgr_info
                   if (ivet(1) == spgr_info(i)%n) then
                      num=i
                      spgm=spgr_info(i)%hall
                      opcion="HAL"
                      exit
                   end if
                end do
             else
                err_symm=.true.
                ERR_Symm_Mess=" Number of Space Group out of limits"
                return
             end if
          end if
       end if  ! present(mode)

       select case (opcion(1:3))
          case ("FIX")
             if (present(gen) .and. present(ngen))  then
                ng=ngen
                istart=1
                num_g=ng-1
                do i=1,ngen
                   call Read_Xsym(gen(i),istart,ss(:,:,i),ts(:,i))
                end do
             else
                err_symm=.true.
                ERR_Symm_Mess=" Generators should be provided if FIX option is Used"
                return
             end if
             call Get_SO_from_FIX(isystm,isymce,ibravl,ng,ss,ts,latsy,co,Spgm,Spacegen(1:1))

             SpaceGroup%Spg_Symb     = "unknown "
             SpaceGroup%Hall         = "unknown "
             SpaceGroup%Laue         = " "
             SpaceGroup%Info         = "Fixed symmetry operators (no info)"
             SpaceGroup%SPG_lat      = Lat_Ch
             SpaceGroup%NumLat       = nlat

             if(allocated(SpaceGroup%Latt_trans)) deallocate(SpaceGroup%Latt_trans)
             allocate(SpaceGroup%Latt_trans(3,nlat))

             SpaceGroup%Latt_trans   = Ltr(:,1:nlat)
             SpaceGroup%Num_gen      = max(0,num_g)
             SpaceGroup%Centre_coord = co
             SpaceGroup%SG_setting   = "Non-Conventional (user-given operators)"
             SpaceGroup%CrystalSys   = " "
             SpaceGroup%Bravais      = Latt(ibravl)
             SpaceGroup%SPG_latsy    = latsy
             SpaceGroup%centred      = isymce
             SpaceGroup%centre       = Centro(isymce)
             SpaceGroup%Numops       = NG

          case ("GEN")
             if (present(gen) .and. present(ngen))  then
                do i=1,ngen
                   call Check_Generator(gen(i),ok)
                   !write(*,"(a,i3,a,tr2,L)") " => Generator # ",i,"  "//trim(gen(i)), ok
                   if(.not. ok) return
                end do
                ng=ngen
                istart=1
                num_g=ng
                call Get_GenSymb_from_Gener(gen,ng,SpaceGroup%gHall)
                do i=1,ngen
                   call Read_Xsym(gen(i),istart,ss(:,:,i),ts(:,i))
                end do
             else
                err_symm=.true.
                ERR_Symm_Mess=" Generators should be provided in GEN calling Set_SpaceGroup"
                return
             end if
             call Get_SO_from_Gener(Isystm,Isymce,Ibravl,Ng,Ss,Ts,Latsy, &
                                    Co,Num_g,Spgm)

             SpaceGroup%CrystalSys   = sys_cry(isystm)
             SpaceGroup%SG_setting   = "Non-Conventional (user-given operators)"
             SpaceGroup%SPG_lat      = Lat_Ch
             SpaceGroup%SPG_latsy    = latsy
             SpaceGroup%NumLat       = nlat
             if(allocated(SpaceGroup%Latt_trans)) deallocate(SpaceGroup%Latt_trans)
             allocate(SpaceGroup%Latt_trans(3,nlat))
             SpaceGroup%Latt_trans   = Ltr(:,1:nlat)
             SpaceGroup%Bravais      = Latt(ibravl)
             SpaceGroup%centre       = Centro(isymce)
             SpaceGroup%centred      = isymce
             SpaceGroup%Centre_coord = co
             SpaceGroup%Numops       = NG
             SpaceGroup%Num_gen      = max(0,num_g)

          case ("HAL")
             call Get_SO_from_Hall (Isystm,Isymce,Ibravl,Ng,Ss,Ts,Latsy, &
                                    Co,Num_g,Spgm)

             if (num > 0) then
                SpaceGroup%NumSpg       = spgr_info(num)%n
                SpaceGroup%Spg_Symb     = spgr_info(num)%hm
                SpaceGroup%Hall         = spgr_info(num)%hall
                call get_laue_str(spgr_info(num)%laue,SpaceGroup%Laue)
                call get_PointGroup_str(spgr_info(num)%pg,SpaceGroup%PG)
                SpaceGroup%Info         = spgr_info(num)%inf_extra
                SpaceGroup%R_Asym_Unit(1,1) = real(spgr_info(num)%asu(1))/24.0
                SpaceGroup%R_Asym_Unit(2,1) = real(spgr_info(num)%asu(2))/24.0
                SpaceGroup%R_Asym_Unit(3,1) = real(spgr_info(num)%asu(3))/24.0
                SpaceGroup%R_Asym_Unit(1,2) = real(spgr_info(num)%asu(4))/24.0
                SpaceGroup%R_Asym_Unit(2,2) = real(spgr_info(num)%asu(5))/24.0
                SpaceGroup%R_Asym_Unit(3,2) = real(spgr_info(num)%asu(6))/24.0
             else
                SpaceGroup%Hall         = Spgm
             end if
             SpaceGroup%CrystalSys   = sys_cry(isystm)
             SpaceGroup%SG_setting   = "Generated from Hall symbol"
             SpaceGroup%SPG_lat      = Lat_Ch
             SpaceGroup%SPG_latsy    = latsy
             SpaceGroup%NumLat       = nlat
             if(allocated(SpaceGroup%Latt_trans)) deallocate(SpaceGroup%Latt_trans)
             allocate(SpaceGroup%Latt_trans(3,nlat))
             SpaceGroup%Latt_trans   = Ltr(:,1:nlat)
             SpaceGroup%Bravais      = Latt(ibravl)
             SpaceGroup%centre       = Centro(isymce)
             SpaceGroup%centred      = isymce
             SpaceGroup%Centre_coord = co
             SpaceGroup%Numops       = NG
             SpaceGroup%Num_gen      = max(0,num_g)

          case ("HMS")
             i=index(SpaceGen,":")
             co=0.0
             if (i /=0 .and. num > 0) then
                spgm=spgr_info(num)%hall
                call Get_SO_from_Hall (Isystm,Isymce,Ibravl,Ng,Ss,Ts,Latsy, &
                                       Co,Num_g,Spgm)
             else
                if (i /= 0) then
                   Spgm=SpaceGen(1:i-1)
                else
                   Spgm=SpaceGen
                end if
                call Get_SO_from_HMS  (Isystm,Isymce,Ibravl,Ng,Ss,Ts,Latsy, &
                                       Spgm)
             end if
             if (num > 0) then
                SpaceGroup%NumSpg       = spgr_info(num)%n
                SpaceGroup%Spg_Symb     = spgr_info(num)%hm
                SpaceGroup%Hall         = spgr_info(num)%hall
                call get_laue_str(spgr_info(num)%laue,SpaceGroup%Laue)
                call get_PointGroup_str(spgr_info(num)%pg,SpaceGroup%PG)
                SpaceGroup%Info         = spgr_info(num)%inf_extra
                SpaceGroup%R_Asym_Unit(1,1) = real(spgr_info(num)%asu(1))/24.0
                SpaceGroup%R_Asym_Unit(2,1) = real(spgr_info(num)%asu(2))/24.0
                SpaceGroup%R_Asym_Unit(3,1) = real(spgr_info(num)%asu(3))/24.0
                SpaceGroup%R_Asym_Unit(1,2) = real(spgr_info(num)%asu(4))/24.0
                SpaceGroup%R_Asym_Unit(2,2) = real(spgr_info(num)%asu(5))/24.0
                SpaceGroup%R_Asym_Unit(3,2) = real(spgr_info(num)%asu(6))/24.0
             else
                SpaceGroup%Spg_Symb     = SpaceGen
                SpaceGroup%Num_gen= 0    !unknown
             end if
             SpaceGroup%CrystalSys   = sys_cry(isystm)
             if (i /=0 .and. num > 0) then
                SpaceGroup%SG_setting   = "Generated from Hall symbol"
                SpaceGroup%Num_gen=max(0,num_g)
             else
                SpaceGroup%SG_setting   ="IT (Generated from Hermann-Mauguin symbol)"
                if(num > 0) then
                   Select Case (spgr_info(num)%n)
                     case(1:2)
                        SpaceGroup%Num_gen= 0    !triclinic
                     case(3:15)
                        SpaceGroup%Num_gen= 1    !monoclinic
                     case(16:74)
                        SpaceGroup%Num_gen= 2    !orthorhombic
                     case(75:88)
                        SpaceGroup%Num_gen= 1    !tetragonal
                     case(89:142)
                        SpaceGroup%Num_gen= 2    !tetragonal
                     case(143:148)
                        SpaceGroup%Num_gen= 1    !trigonal
                     case(149:167)
                        SpaceGroup%Num_gen= 2    !trigonal
                     case(168:176)
                        SpaceGroup%Num_gen= 1    !hexagonal
                     case(177:194)
                        SpaceGroup%Num_gen= 2    !hexagonal
                     case(195:230)
                        SpaceGroup%Num_gen= 3    !cubic
                     case default
                        SpaceGroup%Num_gen= 0    !unknown
                   End Select
                end if
             end if
             SpaceGroup%SPG_lat      = Lat_Ch
             SpaceGroup%SPG_latsy    = latsy
             SpaceGroup%NumLat       = nlat
             if(allocated(SpaceGroup%Latt_trans)) deallocate(SpaceGroup%Latt_trans)
             allocate(SpaceGroup%Latt_trans(3,nlat))
             SpaceGroup%Latt_trans   = Ltr(:,1:nlat)
             SpaceGroup%Bravais      = Latt(ibravl)
             SpaceGroup%centre       = Centro(isymce)
             SpaceGroup%centred      = isymce
             SpaceGroup%Centre_coord = co
             SpaceGroup%Numops       = NG

          case ("ITC")

             call get_generators(spgm,gener)
             if (err_symtab) then
                err_symm=.true.
                ERR_Symm_Mess=" Problems in SpaceGroup: "//trim(spgm)//" => the HM symbol or the number is incorrect "
                return
             else  !Decode gener in generators to construct the space group
                k=0
                do i=1,len_trim(gener)
                   if (gener(i:i) == ";") then
                      k=k+1
                      poscol(k)=i
                   end if
                end do
                if (k /= 0) then
                   ssymb=" "
                   ssymb= adjustl(gener(1:poscol(1)-1))
                   call Read_Xsym(ssymb,1,ss(:,:,1),ts(:,1))
                   do i=2,k
                      ssymb=" "
                      ssymb=adjustl(gener(poscol(i-1)+1:poscol(i)-1))
                      call Read_Xsym(ssymb,1,ss(:,:,i),ts(:,i))
                   end do
                   ssymb=" "
                   ssymb=gener(poscol(k)+1:)
                   call Read_Xsym(ssymb,1,ss(:,:,k+1),ts(:,k+1))
                else
                   ssymb=gener
                   call Read_Xsym(ssymb,1,ss(:,:,k+1),ts(:,k+1))
                end if
             end if

             ng=k+1     !k+1 is the number of generators
             num_g=ng
             call Get_SO_from_Gener(Isystm,Isymce,Ibravl,Ng,Ss,Ts,Latsy, &
                                    Co,Num_g,Spgm)
             if (num > 0) then
                SpaceGroup%NumSpg       = spgr_info(num)%n
                SpaceGroup%Spg_Symb     = spgr_info(num)%hm
                SpaceGroup%Hall         = spgr_info(num)%hall
                call get_laue_str(spgr_info(num)%laue,SpaceGroup%Laue)
                call get_PointGroup_str(spgr_info(num)%pg,SpaceGroup%PG)
                SpaceGroup%Info         = spgr_info(num)%inf_extra
                SpaceGroup%R_Asym_Unit(1,1) = real(spgr_info(num)%asu(1))/24.0
                SpaceGroup%R_Asym_Unit(2,1) = real(spgr_info(num)%asu(2))/24.0
                SpaceGroup%R_Asym_Unit(3,1) = real(spgr_info(num)%asu(3))/24.0
                SpaceGroup%R_Asym_Unit(1,2) = real(spgr_info(num)%asu(4))/24.0
                SpaceGroup%R_Asym_Unit(2,2) = real(spgr_info(num)%asu(5))/24.0
                SpaceGroup%R_Asym_Unit(3,2) = real(spgr_info(num)%asu(6))/24.0
                SpaceGroup%SG_setting   = "Generated from explicit IT generators"
             else
                SpaceGroup%Spg_Symb     = SpaceGen
                SpaceGroup%Num_gen= 0    !unknown
             end if

             SpaceGroup%CrystalSys   = sys_cry(isystm)
             SpaceGroup%SPG_lat      = Lat_Ch
             SpaceGroup%SPG_latsy    = latsy
             SpaceGroup%NumLat       = nlat
             if(allocated(SpaceGroup%Latt_trans)) deallocate(SpaceGroup%Latt_trans)
             allocate(SpaceGroup%Latt_trans(3,nlat))
             SpaceGroup%Latt_trans   = Ltr(:,1:nlat)
             SpaceGroup%Bravais      = Latt(ibravl)
             SpaceGroup%centre       = Centro(isymce)
             SpaceGroup%centred      = isymce
             SpaceGroup%Centre_coord = co
             SpaceGroup%Numops       = NG
             SpaceGroup%Num_gen      = max(0,num_g)

          case default
             err_symm=.true.
             ERR_Symm_Mess=" Problems in SpaceGroup"
             return
       end select

       if (err_symm) return
       if (Is_Hexa(ng,ss)) SpaceGroup%Hexa=.true.

       hexa=SpaceGroup%Hexa  !added 24/05/2007

       if (opcion(1:3) /= "FIX") then              !This has been changed of place for allocating
           select case (SpaceGroup%centred)        !the allocatable components properly
              case (0)
                 SpaceGroup%Multip = 2*NG*nlat
              case (1)
                 SpaceGroup%Multip =   NG*nlat
              case (2)
                 SpaceGroup%Multip = 2*NG*nlat
           end select
       else
           SpaceGroup%Multip =   NG
       end if

       !Allocate here the total number of symmetry operators (JRC, Jan2014)

       if(allocated(SpaceGroup%Symop)) deallocate(SpaceGroup%Symop)
       allocate(SpaceGroup%Symop(SpaceGroup%Multip))
       if(allocated(SpaceGroup%SymopSymb)) deallocate(SpaceGroup%SymopSymb)
       allocate(SpaceGroup%SymopSymb(SpaceGroup%Multip))

       do i=1,SpaceGroup%Numops
          SpaceGroup%Symop(i)%Rot(:,:) = ss(:,:,i)
          SpaceGroup%Symop(i)%tr(:)    = ts(:,i)
       end do

       if (opcion(1:3) /= "FIX") then
          m=SpaceGroup%Numops
          if (SpaceGroup%centred == 0) then
             do i=1,SpaceGroup%Numops
                m=m+1
                vec=-ts(:,i)+2.0*SpaceGroup%Centre_coord(:)
                SpaceGroup%Symop(m)%Rot(:,:) = -ss(:,:,i)
                SpaceGroup%Symop(m)%tr(:)    =  modulo_lat(vec)
             end do
          end if
          if (SpaceGroup%centred == 2) then
             do i=1,SpaceGroup%Numops
                m=m+1
                SpaceGroup%Symop(m)%Rot(:,:) = -ss(:,:,i)
                SpaceGroup%Symop(m)%tr(:)    =  modulo_lat(-ts(:,i))
             end do
          end if
          ngm=m
          if (SpaceGroup%NumLat > 1) then

             do L=2,SpaceGroup%NumLat  ! min(SpaceGroup%NumLat,4)  restriction removed Jan2014 (JRC)
                do i=1,ngm
                   m=m+1
                   vec=SpaceGroup%Symop(i)%tr(:) + SpaceGroup%Latt_trans(:,L)
                   SpaceGroup%Symop(m)%Rot(:,:) = SpaceGroup%Symop(i)%Rot(:,:)
                   SpaceGroup%Symop(m)%tr(:)    = modulo_lat(vec)
                end do
             end do
          end if

       end if
       !write(*,"(a)") " => Generating the symetry operators symbols"
       do i=1,SpaceGroup%multip  ! min(SpaceGroup%multip,192) restriction removed Jan2014 (JRC)
          call Get_SymSymb(SpaceGroup%Symop(i)%Rot(:,:), &
                           SpaceGroup%Symop(i)%tr(:)   , &
                           SpaceGroup%SymopSymb(i))
       end do
       !write(*,"(a)") " => done"

       if (num <= 0) then
          call Get_Laue_PG(SpaceGroup, SpaceGroup%Laue, SpaceGroup%PG)
       end if
       !write(*,"(a)") " => Point group done"

       if(isymce == 0) then
          SpaceGroup%centre = trim(SpaceGroup%centre)//"  Gen(-1):"//SpaceGroup%SymopSymb(NG+1)
       end if

       if(opcion(1:3)=="GEN") call Get_HallSymb_from_Gener(SpaceGroup)

       !write(*,"(a)") " => Wyckoff information"

       !---- Wyckoff information ----!
       if (len_trim(SpaceGroup%Spg_Symb) /= 0) then
          do i=1,273
             if (SpaceGroup%Spg_Symb(1:12) /= wyckoff_info(i)%hm) cycle
             SpaceGroup%Wyckoff%num_orbit=wyckoff_info(i)%norbit
             do j=1,wyckoff_info(i)%norbit

                call wyckoff_orbit(SpaceGroup,wyckoff_info(i)%corbit(j), &
                                   SpaceGroup%Wyckoff%Orbit(j)%norb,     &
                                   SpaceGroup%Wyckoff%Orbit(j)%Str_Orbit)
                SpaceGroup%Wyckoff%Orbit(j)%multp=SpaceGroup%Wyckoff%Orbit(j)%norb*spacegroup%numlat
             end do
             exit
          end do
          SpaceGroup%Spg_Symb(2:)=l_case(SpaceGroup%Spg_Symb(2:))  !Make lowercase the HM generators of the group
       end if
       !write(*,"(a)") " => Wyckoff done"

       return
    End Subroutine Set_SpaceGroup

    !!----
    !!---- Subroutine Set_SpG_Mult_Table(SpG,tab,complete)
    !!----   Type(Space_Group_Type),    intent (in)    :: SpG
    !!----   integer, dimension(:,:),   intent (out)   :: tab
    !!----   logical, optional,         intent (in)    :: complete
    !!----
    !!----   Subroutine to construct the multiplication table of the factor group of
    !!----   a space group. Two operators are equal if they differ only in a lattice
    !!----   translation. The multiplication table is a square matrix with integer
    !!----   numbers corresponding to the ordering of operators in the space group
    !!----   If "complete" is not present, or if complete=.false., we consider only
    !!----   the symmetry operators corresponding to the "primitive" content of the
    !!----   unit cell, so a maximun 48x48 matrix is needed to hold the table in this
    !!----   case. If complete is present and .true., the full table is constructed.
    !!----
    !!----
    !!----  Update: April 2005
    !!----

    Subroutine Set_SpG_Mult_Table(SpG,tab,complete)
      Type(Space_Group_Type),    intent (in)    :: SpG
      integer, dimension(:,:),   intent (out)   :: tab
      logical, optional,         intent (in)    :: complete

       !---- Local Variables ----!
       Type(Sym_Oper_Type) :: Opi,Opj,Opk
       integer :: i,j, ng, k
       logical :: eqvo
       character(len=1) :: lat

       tab=0
       lat=SpG%SPG_lat
       ng=SpG%Numops
       if(SpG%Centred /= 1) ng=2*ng
       if(present(complete)) then
         if(complete) then
           lat="P"
           ng=SpG%Multip
         end if
       end if

       do i=1,ng
         Opi=SpG%SymOp(i)
         do j=1,ng
           Opj=SpG%SymOp(j)
           Opk=Opi*Opj
           do k=1,ng
             eqvo= Equiv_Symop(Opk,SpG%SymOp(k),lat)
             if(eqvo) then
               tab(i,j)=k
               exit
             end if
           end do
           if(tab(i,j) == 0) then
             err_symm=.true.
             ERR_Symm_Mess=" Problems constructing the multiplication Table of the space group: "//trim(spg%spg_symb)
             return
           end if
         end do
       end do

      return
    End Subroutine Set_SpG_Mult_Table

    !!--++
    !!--++ Subroutine Setting_Change_Conv(From_Syst,To_Syst,Spacegroup, Car_Sym, Icar_Sym)
    !!--++    character(len=2),    intent(in)     :: From_Syst   !  In -> IT : International Tables
    !!--++                                                                ML : Miller & Love
    !!--++                                                                KO : Kovalev
    !!--++                                                                BC : Bradley & Cracknell
    !!--++                                                                ZA : Zack
    !!--++    character(len=2),    intent(in)     :: To_Syst     !  In -> (Idem to From_Syst)
    !!--++    type (Space_Group),  intent(in out) :: SpaceGroup  !  In ->
    !!--++                                                         Out ->
    !!--++    character(len=35),    intent(out)   :: car_sym     ! Out ->
    !!--++    character(len=35),    intent(out)   :: icar_sym    ! Out ->
    !!--++
    !!--++    Traslate From From_Syst to To_syst the set of symmetry operators
    !!--++
    !!--++   Update: February - 2005 (Name changed and overloaded by JRC in Jan2014)
    !!
    Subroutine Setting_Change_Conv(From_Syst, To_Syst, SpaceGroup, car_sym, icar_sym)
       !---- Arguments ----!
       character(len=2),          intent(in)     :: From_Syst, To_Syst
       type (Space_Group_Type),   intent(in out) :: SpaceGroup
       character(len=35),         intent(out)    :: car_sym, icar_sym

       !---- Local Variables ----!
       character(len=2) :: car1, car2
       integer                 :: i,j,num
       integer, dimension(4,4) :: s, si, st, sti, w
       integer, dimension(3,3) :: r, r_inv, rt, rt_inv
       real(kind=cp), dimension(3)      :: t, t_inv, tt, tt_inv

       !---- Initializing variables ----!
       call init_err_symm()
       call Set_System_Equiv()

       car1=From_Syst
       car2=To_Syst
       car_sym=" "
       icar_sym=" "
       call ucase(car1)
       call ucase(car1)

       !---- Checking data ----!
       if (len_trim (car1) == 0) then
          err_symm=.true.
          ERR_Symm_Mess=" Blank Option"
          return
       end if
       if (len_trim (car2) == 0) then
          err_symm=.true.
          ERR_Symm_Mess=" Blank Option"
          return
       end if
       if (SpaceGroup%NumSpg <= 0 .or. SpaceGroup%NumSpg > 230 ) then
          err_symm=.true.
          ERR_Symm_Mess=" Space Group Not Defined..."
          return
       end if
       SpaceGroup%SG_setting="Changed from "//car1//" to "//car2
       num=SpaceGroup%NumSpg
       r     = 0
       r_inv = 0
       rt    = 0
       rt_inv= 0
       t     = 0.0
       t_inv = 0.0
       tt    = 0.0
       tt_inv= 0.0
       do i=1,3
          r(i,i)      = 1
          r_inv(i,i)  = 1
          rt(i,i)     = 1
          rt_inv(i,i) = 1
       end do
       s     = 0
       si    = 0
       st    = 0
       sti   = 0
       w     = 0
       do i=1,4
          s(i,i)   = 1
          si(i,i)  = 1
          st(i,i)  = 1
          sti(i,i) = 1
          w(i,i)   = 1
       end do

       select case (car1)
          case ("IT")    !---- International Tables ----!
             select case (car2)
                case ("IT")
                   return
                case ("ML")
                   car_sym=system_equiv(num)%ml
                case ("KO")
                   car_sym=system_equiv(num)%ko
                case ("BC")
                   car_sym=system_equiv(num)%bc
                case ("ZA")
                   car_sym=system_equiv(num)%za
             end select
             j=1
             call read_Xsym(car_sym,j,r,t)
             call inverse_symm(r,t,r_inv,t_inv)
             if (err_symm) return

          case ("ML")    !---- Miller & Love ----!
             select case (car2)
                case ("IT")
                   car_sym=system_equiv(num)%ml
                   j=1
                   call read_Xsym(car_sym,j,r_inv,t_inv)
                   call inverse_symm(r_inv,t_inv,r,t)
                   if (err_symm) return

                case ("ML")
                   return

                case ("KO")
                   car_sym=system_equiv(num)%ml
                   j=1
                   call read_Xsym(car_sym,j,r_inv,t_inv)
                   call inverse_symm(r_inv,t_inv,r,t)
                   if (err_symm) return
                   car_sym=system_equiv(num)%ko
                   j=1
                   call read_Xsym(car_sym,j,rt,tt)
                   call inverse_symm(rt,tt,rt_inv,t_inv)
                   if (err_symm) return

                case ("BC")
                   car_sym=system_equiv(num)%ml
                   j=1
                   call read_Xsym(car_sym,j,r_inv,t_inv)
                   call inverse_symm(r_inv,t_inv,r,t)
                   if (err_symm) return
                   car_sym=system_equiv(num)%bc
                   j=1
                   call read_Xsym(car_sym,j,rt,tt)
                   call inverse_symm(rt,tt,rt_inv,t_inv)
                   if (err_symm) return

                case ("ZA")
                   car_sym=system_equiv(num)%ml
                   j=1
                   call read_Xsym(car_sym,j,r_inv,t_inv)
                   call inverse_symm(r_inv,t_inv,r,t)
                   if (err_symm) return
                   car_sym=system_equiv(num)%za
                   j=1
                   call read_Xsym(car_sym,j,rt,tt)
                   call inverse_symm(rt,tt,rt_inv,t_inv)
                   if (err_symm) return
             end select

          case ("KO")    !---- Kovalev ----!
             select case (car2)
                case ("IT")
                   car_sym=system_equiv(num)%ko
                   j=1
                   call read_Xsym(car_sym,j,r_inv,t_inv)
                   call inverse_symm(r_inv,t_inv,r,t)
                   if (err_symm) return

                case ("ML")
                   car_sym=system_equiv(num)%ko
                   j=1
                   call read_Xsym(car_sym,j,r_inv,t_inv)
                   call inverse_symm(r_inv,t_inv,r,t)
                   if (err_symm) return
                   car_sym=system_equiv(num)%ml
                   j=1
                   call read_Xsym(car_sym,j,rt,tt)
                   call inverse_symm(rt,tt,rt_inv,t_inv)
                   if (err_symm) return

                case ("KO")
                   return

                case ("BC")
                   car_sym=system_equiv(num)%ko
                   j=1
                   call read_Xsym(car_sym,j,r_inv,t_inv)
                   call inverse_symm(r_inv,t_inv,r,t)
                   if (err_symm) return
                   car_sym=system_equiv(num)%bc
                   j=1
                   call read_Xsym(car_sym,j,rt,tt)
                   call inverse_symm(rt,tt,rt_inv,t_inv)
                   if (err_symm) return

                case ("ZA")
                   car_sym=system_equiv(num)%ko
                   j=1
                   call read_Xsym(car_sym,j,r_inv,t_inv)
                   call inverse_symm(r_inv,t_inv,r,t)
                   if (err_symm) return
                   car_sym=system_equiv(num)%za
                   j=1
                   call read_Xsym(car_sym,j,rt,tt)
                   call inverse_symm(rt,tt,rt_inv,t_inv)
                   if (err_symm) return
             end select

          case ("BC")    !---- Bradley & Cracknell ----!
             select case (car2)
                case ("IT")
                   car_sym=system_equiv(num)%bc
                   j=1
                   call read_Xsym(car_sym,j,r_inv,t_inv)
                   call inverse_symm(r_inv,t_inv,r,t)
                   if (err_symm) return

                case ("ML")
                   car_sym=system_equiv(num)%bc
                   j=1
                   call read_Xsym(car_sym,j,r_inv,t_inv)
                   call inverse_symm(r_inv,t_inv,r,t)
                   if (err_symm) return
                   car_sym=system_equiv(num)%ml
                   j=1
                   call read_Xsym(car_sym,j,rt,tt)
                   call inverse_symm(rt,tt,rt_inv,t_inv)
                   if (err_symm) return

                case ("KO")
                   car_sym=system_equiv(num)%bc
                   j=1
                   call read_Xsym(car_sym,j,r_inv,t_inv)
                   call inverse_symm(r_inv,t_inv,r,t)
                   if (err_symm) return
                   car_sym=system_equiv(num)%ko
                   j=1
                   call read_Xsym(car_sym,j,rt,tt)
                   call inverse_symm(rt,tt,rt_inv,t_inv)
                   if (err_symm) return

                case ("BC")
                   return

                case ("ZA")
                   car_sym=system_equiv(num)%bc
                   j=1
                   call read_Xsym(car_sym,j,r_inv,t_inv)
                   call inverse_symm(r_inv,t_inv,r,t)
                   if (err_symm) return
                   car_sym=system_equiv(num)%za
                   j=1
                   call read_Xsym(car_sym,j,rt,tt)
                   call inverse_symm(rt,tt,rt_inv,t_inv)
                   if (err_symm) return
             end select

          case ("ZA")    !---- Zak ----!
             select case (car2)
                case ("IT")
                   car_sym=system_equiv(num)%za
                   j=1
                   call read_Xsym(car_sym,j,r_inv,t_inv)
                   call inverse_symm(r_inv,t_inv,r,t)
                   if (err_symm) return

                case ("ML")
                   car_sym=system_equiv(num)%za
                   j=1
                   call read_Xsym(car_sym,j,r_inv,t_inv)
                   call inverse_symm(r_inv,t_inv,r,t)
                   if (err_symm) return
                   car_sym=system_equiv(num)%ml
                   j=1
                   call read_Xsym(car_sym,j,rt,tt)
                   call inverse_symm(rt,tt,rt_inv,t_inv)
                   if (err_symm) return

                case ("KO")
                   car_sym=system_equiv(num)%za
                   j=1
                   call read_Xsym(car_sym,j,r_inv,t_inv)
                   call inverse_symm(r_inv,t_inv,r,t)
                   if (err_symm) return
                   car_sym=system_equiv(num)%ko
                   j=1
                   call read_Xsym(car_sym,j,rt,tt)
                   call inverse_symm(rt,tt,rt_inv,t_inv)
                   if (err_symm) return

                case ("BC")
                   car_sym=system_equiv(num)%za
                   j=1
                   call read_Xsym(car_sym,j,r_inv,t_inv)
                   call inverse_symm(r_inv,t_inv,r,t)
                   if (err_symm) return
                   car_sym=system_equiv(num)%bc
                   j=1
                   call read_Xsym(car_sym,j,rt,tt)
                   call inverse_symm(rt,tt,rt_inv,t_inv)
                   if (err_symm) return

                case ("ZA")
                   return
             end select

       end select

       call Get_SymSymb(rt_inv,t_inv,icar_sym)
       s(1:3,1:3)  = r
       s(1:3,4)    = mod(nint(t*24.0)+48,24)
       si(1:3,1:3) = r_inv
       si(1:3,4)   = mod(nint(t_inv*24.0)+48,24)

       st(1:3,1:3) = rt
       st(1:3,4)   = mod(nint(tt*24.0)+48,24)
       sti(1:3,1:3)= rt_inv
       sti(1:3,4)  = mod(nint(tt_inv*24.0)+48,24)

       s=matmul(st,s)
       si=matmul(si,sti)

       do i=1,SpaceGroup%multip
          w(1:3,1:3) = SpaceGroup%Symop(i)%Rot
          w(1:3,4)   = mod(nint(SpaceGroup%Symop(i)%Tr*24.0)+48,24)
          w=matmul(s,w)
          w=matmul(w,si)
          SpaceGroup%Symop(i)%Rot = w(1:3,1:3)
          SpaceGroup%Symop(i)%Tr  = mod(real(w(1:3,4)/24.0)+10.0_cp,1.0_cp)
       end do
       do i=1,SpaceGroup%numops
          call Get_SymSymb(SpaceGroup%Symop(i)%Rot,  &
                           SpaceGroup%Symop(i)%tr,SpaceGroup%SymopSymb(i))
       end do

       return
    End Subroutine Setting_Change_Conv

    !!--++
    !!--++ Subroutine Setting_Change_NonConv(Mat,Orig,Spg,Spgn,Matkind)
    !!--++   real(kind=cp), dimension(3,3),intent(in) :: mat      ! Basis transformation matrix
    !!--++   real(kind=cp), dimension(3),  intent(in) :: orig     ! New origing in the old basis
    !!--++   type (Space_Group_Type),      intent(in) :: SpG      ! Input space group
    !!--++   type (NS_Space_Group_Type),  intent(out) :: SpGn     ! New space group in the new setting.
    !!--++   character (len=*), optional,  intent(in) :: matkind  ! Kind of transformation matrix
    !!--++
    !!--++    Transform the symmetry operators of the space group to a new basis given by
    !!--++    the matrix "mat" and vector "orig"
    !!--++    If matkind is given and matkind="it"/"IT", the input matrix is given
    !!--++    as in International Tables: Mat=P =>  (a,b,c)'=(a,b,c)P
    !!--++    Otherwise it is the trasposed matrix Mat=Pt
    !!--++
    !!--++ Created: January - 2014 (JRC)
    !!
    Subroutine Setting_Change_NonConv(Mat,Orig,Spg,Spgn,Matkind)
       !---- Arguments ----!
       real(kind=cp), dimension(3,3),intent(in) :: Mat
       real(kind=cp), dimension(3),  intent(in) :: Orig
       type (Space_Group_Type),      intent(in) :: SpG
       type (NS_Space_Group_Type),  intent(out) :: SpGn
       character (len=*), optional,  intent(in) :: Matkind

       !--- Local variables ---!
       integer                 :: ifail, i, j, k, L, im, nc, m, ngm,n,ngen
       real(kind=cp)           :: det
       character(len=40)       :: transla
       character(len=1)        :: LatSymb
       real(kind=cp), dimension (3,3), parameter :: e = reshape ((/1.0,0.0,0.0,  &
                                                                   0.0,1.0,0.0,  &
                                                                   0.0,0.0,1.0/),(/3,3/))
       real(kind=cp), dimension (3,192)    :: newlat = 0.0 !big enough number of centring tranlations
       real(kind=cp), dimension (3,3)      :: S, Sinv, rot, rotn  !S is the ITC matrix P.
       integer,       dimension (3,3)      :: nulo
       real(kind=cp), dimension (  3)      :: tr, trn, v
       logical                             :: lattl,change_only_origin
       character(len=80)                   :: symbsg
       character(len=60),dimension(15)     :: gen
       character(len=180)                  :: setting
       real(kind=cp),  dimension(3,3,Spg%Multip) :: sm
       real(kind=cp),  dimension(3,Spg%Multip)   :: tm

       call Init_Err_Symm()
       change_only_origin=.false.
       nulo=0
       call get_setting_info(Mat,orig,setting,matkind)
       symbsg=Pack_String(SpG%spg_symb)
       if (present(matkind)) then
          if (matkind(1:2) == "it" .or. matkind(1:2) == "IT" ) then
             S=Mat
          else
             S=transpose(Mat)
          end if
       else
          S=transpose(Mat)
       end if
       setting = trim(setting)//" det:"
       if(equal_matrix(S,e,3)) change_only_origin=.true.
       det=determ_a(Mat)
       i=len_trim(setting)
       write(unit=setting(i+2:),fmt="(f6.2)") det
       !write(unit=*,fmt="(a)") " => Setting Symbol: "//trim(setting)
       call matrix_inverse(S,Sinv,ifail)
       if (ifail /= 0) then
          err_symm=.true.
          ERR_Symm_Mess= "Inversion Matrix Failed on: Setting_Change_NonConv"
          return
       end if

       L=0
       if (SpG%NumLat > 1) then  !Original lattice is centered
          do i=2,SpG%NumLat      !Transform the centring vectors to the new lattice
             v=Modulo_Lat(matmul(Sinv,SpG%Latt_trans(:,i)))
             if (sum(v) < eps_symm) cycle
             L=L+1
             newlat(:,L)=v
          end do
       end if
       do i=1,3  !Test the basis vectors of the original setting
         rot(:,i)=Modulo_Lat(Sinv(:,i))
         if (sum(rot(:,i)) < eps_symm) cycle
         L=L+1
         newlat(:,L)=rot(:,i)
       end do

       if (det > 1 ) then  !The new lattice is centred
          im=nint(det)-1         !Determine the new lattice translations
          ngm=L+im
          doi: do i=0,im
             v(1) = i
             do j=0,im
                v(2) = j
                do k=0,im
                   v(3) = k
                   if (nint(sum(v)) == 0) cycle
                   tr=Modulo_Lat(matmul(Sinv,v))
                   if (sum(tr) < eps_symm) cycle
                   lattl =.true.
                   do m=1,L
                      if (sum(abs(tr-newlat(:,m))) < eps_symm) then
                         lattl =.false.
                         exit
                      end if
                   end do
                   if (lattl) then ! new lattice translation
                      L=L+1
                      newlat(:,L) = tr(:)
                      if (L == ngm) exit doi
                   end if
                end do !k
             end do !j
          end do doi !i
       end if

       call get_centring_vectors(L,newlat,LatSymb)  !Complete the centring vectors
       !Now we have L centring translations
       call LatSym(LatSymb,L,newlat)  !provides the value of the global variable inlat: index of the type of lattice
       SpGn%SPG_lat      = LatSymb
       SpGn%SPG_latsy    = SpG%SPG_latsy(1:1)//LatSymb

       !---- Change of symmetry operator under a change of basis and origin
       !----  A'= M A,  origin O =>  X'=inv(Mt)(X-O)
       !----  Symmetry operator C = (R,T)  -> C' = (R',T')
       !----   R' = inv(Mt) R Mt                 ITC:    R'= inv(P) R P
       !----   T' = inv(Mt) (T -(E-R)O)                  T'= inv(P) (T-(E-R)O)
       sm=0.0
       tm=0.0
       sm(:,:,1)=SpG%SymOp(1)%Rot
       tm(:,1)=SpG%SymOp(1)%tr
       n=1
       do_i:do i=2,SpG%NumOps
          Rot=SpG%SymOp(i)%rot
          Rotn=matmul(matmul(Sinv,Rot),S)
          !irot=nint(Rotn)
          do k=n,1,-1
            if(equal_matrix(Rotn,sm(:,:,k),3))  cycle do_i
          end do
          n=n+1
          sm(:,:,N)=Rotn
          tr=SpG%SymOp(i)%tr
          trn=matmul(Sinv,tr-matmul(e-Rot,orig))
          tm(:,n)=Modulo_Lat(trn)
       end do do_i

       SpGn%Centred=SpG%Centred
       SpGn%Centre_coord=SpG%Centre_coord
       if (SpG%Centred /= 1) then !the space group is centro-symmetric
          nc=SpG%NumOps+1
          Rot=SpG%SymOp(nc)%rot
          tr=SpG%SymOp(nc)%tr
          trn=matmul(Sinv,tr-matmul(e-Rot,orig)) ! matmul(Sinv,tr-2*orig)
          trn= Modulo_Lat(trn)
          if(sum(abs(trn)) > 3.0*eps_symm) then
            SpGn%Centred=0
            SpGn%Centre_coord=0.5*trn
          end if
       end if

       !Do another thing we conserve the transformations and generate ourself the new group
       !The new multiplicity is
       i=1
       if(SpGn%Centred /= 1) i=2
       SpGn%multip= n * i * nlat  !nlat=L+1

       allocate(SpGn%SymOp(SpGn%multip), SpGn%SymOpSymb(SpGn%multip))
       SpGn%NumOps=n
       do i=1,SpGn%NumOps
         SpGn%SymOp(i)%Rot=sm(:,:,i)
         SpGn%SymOp(i)%tr=tm(:,i)
       end do

       allocate(SpGn%Latt_trans(3,nlat))
       SpGn%NumLat    = nlat
       SpGn%Latt_trans= Ltr(:,1:nlat)
       SpGn%CrystalSys   = SpG%CrystalSys
       SpGn%SG_setting   = setting
       SpGn%Bravais      = Latt(inlat)
       Select Case (SpGn%Centred)
           Case(0,2)
             call Frac_Trans_2Dig(SpGn%Centre_coord,transla)
             SpGn%centre="Centric, -1 at "//trim(transla)
           Case Default
             SpGn%centre="Acentric"
       End Select
       SpGn%Num_gen      = SpG%Num_gen
       SpGn%PG           = SpG%PG
       SpGn%Laue         = SpG%laue
       SpGn%NumSpg=SpG%NumSpg
       m=SpGn%Numops
       if (SpGn%centred /= 1) then
          do i=1,SpGn%Numops
             m=m+1
             SpGn%Symop(m)%Rot(:,:) = -SpGn%Symop(i)%Rot(:,:)
             SpGn%Symop(m)%tr(:)    =  modulo_lat(-SpGn%Symop(i)%tr(:)+2.0*SpGn%Centre_coord)
          end do
       end if
       ngm=m
       if (SpGn%NumLat > 1) then
          do L=2,SpGn%NumLat
             do i=1,ngm
                m=m+1
                trn=SpGn%Symop(i)%tr(:) + SpGn%Latt_trans(:,L)
                SpGn%Symop(m)%Rot(:,:) = SpGn%Symop(i)%Rot(:,:)
                SpGn%Symop(m)%tr(:)    = modulo_lat(trn)
             end do
          end do
       end if
       do i=1,SpGn%multip
          call Get_SymSymb(SpGn%Symop(i)%Rot(:,:), &
                           SpGn%Symop(i)%tr(:)   , &
                           SpGn%SymopSymb(i))
       end do
       !Try to assign a Hall symbol to the space group in the new setting
       !If the hall symbol has been found and the symbol exists in the table the H-M symbol is also set.
       SpGn%hall="From:"//trim(SpG%hall)
       SpGn%spg_symb="From:"//trim(SpG%spg_symb)
       if(change_only_origin) then
         SpGn%spg_symb=trim(symbsg)
       else
         if(SpGn%NumSpg == 0) then
            SpGn%spg_symb="From:"//trim(symbsg)
         end if
       end if
       !Generate a general symbol, first select generators as a function of SpGn$Numops
       n=SpGn%NumOps
       Select Case(SpGn%centred)
         Case(0)
           gen(1)=SpGn%SPG_lat
           gen(2)=SpGn%SymopSymb(n+1)
           ngen=2
         Case(1)
           gen(1)=SpGn%SPG_lat
           ngen=1
         Case(2)
           gen(1)="-"//SpGn%SPG_lat
           ngen=1
       End Select

       Select Case(n)
         case(1:3)
           ngen=ngen+1
           gen(ngen)=SpGn%SymopSymb(2)
         case(4:)
           ngen=ngen+1
           gen(ngen)=SpGn%SymopSymb(2)
           ngen=ngen+1
           gen(ngen)=SpGn%SymopSymb(3)
       End Select
       !
       call Get_GenSymb_from_Gener(gen,ngen,SpGn%ghall)

       return
    End Subroutine Setting_Change_NonConv

    !!----
    !!---- Subroutine Similar_Transf_Sg(Mat,Orig,Spg,Spgn,Matkind)
    !!----    real(kind=cp), dimension (3,3),   intent( in)    :: Mat     ! Matrix transforming the basis
    !!----    real(kind=cp), dimension (  3),   intent( in)    :: orig    ! Coordinates of the new origin
    !!----    type (Space_Group_Type) ,         intent( in)    :: SpG     ! Initial space group
    !!----    type (Space_Group_Type) ,         intent(out)    :: SpGn    ! Maximum subgroup of SpG
    !!----    character (len=*), optional,      intent( in)    :: matkind ! Type of the input matrix
    !!----    character (len=*), optional,      intent( in)    :: Fix_lat ! Fixing Lattice type
    !!----
    !!----    Subroutine to construct a space group "SpGn" that is a maximal subgroup
    !!----    of the input space group "SpG" compatible with the transformation
    !!----    of the basis corresponding to the matrix "Mat" and the new origin "orig".
    !!----    The transformed SpGn will have (if it is the case) conventional centring vectors.
    !!----    If matkind is given and matkind="it"/"IT", the input matrix is given
    !!----    as in International Tables:
    !!--<<
    !!----                      (a' b' c') = (a b c) Mat
    !!-->>
    !!----    If matkind is not given or if it is not equal to "it"/"IT" the input matrix
    !!----    is the transpose of the International convention (column matrices for basis vectors)
    !!----    The new space group is obtained using the properties of conventional Bravais
    !!----    lattices and symmetry operators. Only the symmetry operators of the conventionnal
    !!----    form are retained to construct the new space group. If the Hermann-Mauguin symbol
    !!----    is not given, that means it correspond to a special setting. The Hall symbol is
    !!----    always given.
    !!----    The coordinates of the origin is always given with respect to the (a b c) basis.
    !!----    If Fix_lat is given a conventional lattice centring, this is fixed irrespective
    !!----    of the centring obtained by applying the similarity transformation. For instance
    !!----    is Fix_lat="P" and the transformation implies new centring vectors or the input
    !!----    group is centred, the generators with fraccional translations are removed from
    !!----    the group. If Fix_lat="A" (or whatever) the program will add the corresponding
    !!----    generators irrespective that the generator is in the original/transformed group.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Similar_Transf_SG(Mat,orig,SpG,SpGn,matkind,Fix_lat)
       !---- Arguments ----!
       real(kind=cp), dimension (3,3),   intent( in)    :: Mat
       real(kind=cp), dimension (  3),   intent( in)    :: orig
       type (Space_Group_Type) ,         intent( in)    :: SpG
       type (Space_Group_Type) ,         intent(out)    :: SpGn
       character (len=*), optional,      intent( in)    :: matkind
       character (len=*), optional,      intent( in)    :: Fix_lat

       !--- Local variables ---!
       integer                 :: ifail, i, j, k, det, L, im, nc, m, ngm, ngen, Isystm
       real(kind=cp), dimension (3,3), parameter :: e = reshape ((/1.0,0.0,0.0,  &
                                                                   0.0,1.0,0.0,  &
                                                                   0.0,0.0,1.0/),(/3,3/))
       real(kind=cp), dimension (3,192):: newlat = 0.0 !big enough number of centring tranlations
       real(kind=cp), dimension (3,3)  :: S, Sinv, rot, rotn
       integer,       dimension (3,3)  :: irot
       real(kind=cp), dimension (  3)  :: tr, trn, v
       real(kind=cp)                   :: rmin,rmax
       logical                         :: latt
       character(len=40),dimension(60) :: gen
       integer,       dimension (60)   :: pt
       character(len=40)               :: string
       character(len=80)               :: setting, symbsg
       character(len=12)               :: csys
       character(len=1)                :: lattsymb, crys

       err_symm=.false.
       call get_setting_info(Mat,orig,setting,matkind)
       symbsg=Pack_String(SpG%spg_symb)
       csys=SpG%CrystalSys
       if (present(matkind)) then
          if (matkind(1:2) == "it" .or. matkind(1:2) == "IT" ) then
             S=Mat
          else
             S=transpose(Mat)
          end if
       else
          S=transpose(Mat)
       end if

       setting = trim(setting)
       det=determ_a(Mat)
       call matrix_inverse(S,Sinv,ifail)
       if (ifail /= 0) then
          err_symm=.true.
          ERR_Symm_Mess= "Inversion Matrix Failed on: similar_SG"
          return
       end if

       if(present(Fix_lat)) then
          lattsymb=Fix_lat
       else
          L=0
          if (SpG%NumLat > 1) then  !Original lattice is centered
             do i=2,SpG%NumLat      !Transform the centring vectors to the new lattice
                v=Modulo_Lat(matmul(Sinv,SpG%Latt_trans(:,i)))
                if (sum(v) < eps_symm) cycle
                L=L+1
                newlat(:,L)=v
             end do
          end if

          do i=1,3  !Test also the basis vectors of the original setting
            rot(:,i)=Modulo_Lat(Sinv(:,i))
            if (sum(rot(:,i)) < eps_symm) cycle
            L=L+1
            newlat(:,L)=rot(:,i)
          end do

          if (det > 1 ) then  !The new lattice is centred
             im=det-1         !Determine the new lattice translations
             ngm=L+im
             doi: do i=0,im
                v(1) = i
                do j=0,im
                   v(2) = j
                   do k=0,im
                      v(3) = k
                      if (nint(sum(v)) == 0) cycle
                      tr=Modulo_Lat(matmul(Sinv,v))
                      if (sum(tr) < eps_symm) cycle
                      latt =.true.
                      do m=1,L
                         if (sum(abs(tr-newlat(:,m))) < eps_symm) then
                            latt =.false.
                            exit
                         end if
                      end do
                      if (latt) then ! new lattice translation
                         L=L+1
                         newlat(:,L) = tr(:)
                         if (L == ngm) exit doi
                      end if
                   end do !k
                end do !j
             end do doi !i
          end if

          call get_centring_vectors(L,newlat,lattsymb)
          if(lattsymb == "Z") call Max_Conv_Lattice_Type(L, newlat, lattsymb)
          !newlat is not used anymore
       end if
       !---- Select the generators of the maximum conventional lattice
       ngen=0
       Select Case(lattsymb)
          Case("A")
             ngen=1
             gen(1)="x,y+1/2,z+1/2"
          Case("B")
             ngen=1
             gen(1)="x+1/2,y,z+1/2"
          Case("C")
             ngen=1
             gen(1)="x+1/2,y+1/2,z"
          Case("I")
             ngen=1
             gen(1)="x+1/2,y+1/2,z+1/2"
          Case("F")
             ngen=2
             gen(1)="x+1/2,y+1/2,z"
             gen(2)="x+1/2,y,z+1/2"
          Case("R")
             ngen=1
             gen(1)="x+2/3,y+1/3,z+1/3"
       End Select

       !---- Up to here all "conventionnal" translational generators have been obtained
       !---- Set the minimum and maximum admissible component of translations
       select case (csys)
          Case("Triclinic")
             rmin=0.0
             rmax=1.0
          Case("Monoclinic")
             rmin=0.5
             rmax=0.5
          Case("Orthorhombic")
             rmin=0.5
             rmax=0.5
             if (lattsymb == "F") then
                rmin=0.25
                rmax=0.75
             end if
          Case("Tetragonal")
             rmin=0.25
             rmax=0.75
          Case("Rhombohedral","Hexagonal","Trigonal")
             rmin=1.0/6.0
             rmax=5.0/6.0
          Case("Cubic")
             rmin=0.25
             rmax=0.75
          Case default
             rmin=0.5
             rmax=0.5
       end select

       !---- Change of symmetry operator under a change of basis and origin
       !----  A'= M A,  origin O =>  X'=inv(Mt)(X-O)
       !----  Symmetry operator C = (R,T)  -> C' = (R',T')
       !----   R' = inv(Mt) R Mt                 ITC:    R'= inv(P) R P
       !----   T' = inv(Mt) (T -(E-R)O)                  T'= inv(P) (T-(E-R)O)
       do i=2,SpG%NumOps
          Rot=SpG%SymOp(i)%rot
          tr=SpG%SymOp(i)%tr
          Rotn=matmul(matmul(Sinv,Rot),S)
          irot=abs(nint(rotn))
          if ( any(irot > 1) ) cycle    !Conserve only the conventional forms  |aij|=1,0
          if (.not. Zbelong(Rotn)) cycle
          ! Verify is the associated translation is admissible in the crystal system of
          ! the parent space group.
          trn=matmul(Sinv,tr-matmul(e-Rot,orig))
          trn=Modulo_Lat(trn)
          if ( any((trn < rmin .and. trn > 0.0) .or. trn > rmax) ) cycle  !internal compiler error in gfortran
          call Get_SymSymb(nint(Rotn),trn,string)
          ngen=ngen+1
          gen(ngen)=string
       end do

       !----Obtain the maximum expected crystal system after going to the new setting
       call Get_Crystal_System(Ngen,Gen, Isystm, Crys)

       select case (Isystm)
          Case(1)
             rmin=0.0
             rmax=1.0
          Case(2)
             rmin=0.5
             rmax=0.5
          Case(3)
             rmin=0.5
             rmax=0.5
             if (lattsymb == "F") then
                rmin=0.25
                rmax=0.75
             end if
          Case(4)
             rmin=0.25
             rmax=0.75
          Case(5,6)
             rmin=1.0/6.0
             rmax=5.0/6.0
          Case(7)
             rmin=0.25
             rmax=0.75
          Case default
             rmin=0.0
             rmax=1.0
       end select

       pt(1:ngen) = 1
       do i=1,ngen
          string=gen(i)

          !---- Test if the generator is still compatible with the crystal system
          call Read_Xsym(string,1,iRot,tr)
          if ( any((tr < rmin .and. tr > 0.0) .or. tr > rmax) ) then
             pt(i)=0
             cycle
          end if
          j=index(string,",")
          k=index(string,",",back=.true.)
          select case (Isystm)
             Case(1,2,3)  ! "Triclinic","Monoclinic","Orthorhombic"
                if (index(string(1:j),"y") /= 0 .or. index(string(1:j),"z") /= 0) pt(i)=0
                if (index(string(j:k),"x") /= 0 .or. index(string(1:j),"z") /= 0) pt(i)=0
                if (index(string(k: ),"x") /= 0 .or. index(string(1:j),"y") /= 0) pt(i)=0
             Case(4,5,6)  ! "Tetragonal","Rhombohedral","Hexagonal","Trigonal"
                if (index(string(1:k),"z") /= 0 ) pt(i)=0
                if (index(string(k: ),"x") /= 0 .or. index(string(k: ),"y") /= 0) pt(i)=0
          end select
       end do

       m=0
       do i=1,ngen
          string=gen(i)
          if (pt(i) == 1) then
             m=m+1
             gen(m)=string
          end if
       end do
       ngen=m
       if (SpG%Centred /= 1) then !the space group is centro-symmetric
          nc=SpG%NumOps+1
          Rot=SpG%SymOp(nc)%rot
          tr=SpG%SymOp(nc)%tr
          trn=matmul(Sinv,tr-matmul(e-Rot,orig)) ! matmul(Sinv,tr-2*orig)
          trn= Modulo_Lat(trn)
          if(Lattice_Trans(trn,lattsymb)) trn=(/0.0,0.0,0.0/) !Check Lattice centring

          if (.not. any((trn < rmin .and. trn > 0.0) .or. trn > rmax) ) then
             ngen=ngen+1
             call Get_SymSymb(SpG%SymOp(nc)%rot,trn,gen(ngen))
          end if

       end if

       !---- Check if non conventionnal centring vectors have been generated from
       !---- the given generators. In such a case reduce by one unit the number of
       !---- generators and restart the generation
       dob:do
          call set_spacegroup("  ",SpGn,gen,ngen,"GEN")
          do i=1,SpGn%multip
             call symmetry_symbol(SpGn%SymOp(i),string)
             string=adjustl(string)
             if (string(1:1) == "t") then
                if (lattice_trans(SpGn%SymOp(i)%tr,SpGn%SPG_lat)) cycle
                ngen=ngen-1
                cycle dob
             end if
          end do
          exit
       end do dob

       If(present(Fix_Lat)) then
          SpGn%spg_symb="From("//trim(symbsg)//") Lat:"//Fix_lat
       else
          SpGn%spg_symb="From("//trim(symbsg)//")"
       end if
       call get_HallSymb_from_gener(SpGn)
       SpGn%SG_setting=setting

       return
    End Subroutine Similar_Transf_SG


    !!----
    !!---- Subroutine Sym_B_Relations(Op/Symb,B_Ind,B_Fac)
    !!----    integer, dimension(3,3),     intent (in) :: Op      !  In  -> Rotation Matrix
    !!----    character(len=*),            intent (in) :: Symb    !  In  -> Symmetry string
    !!----
    !!----    integer, dimension(6),       intent(out) :: B_Ind   !  Out -> B Index
    !!----    real(kind=cp), dimension(6), intent(out) :: B_Fac   !  Out -> B Factor
    !!----
    !!----    Symmetry relations among coefficients of the anisotropic temperature
    !!----    factor.
    !!----
    !!----    Order for B is: B11 B22 B33 B12 B13 B23
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine Sym_B_Relations_Op(R,B_Ind,B_Fac)
    !!--++    integer,dimension(3,3),      intent (in) :: R
    !!--++    integer, dimension(6),       intent(out) :: B_Ind
    !!--++    real(kind=cp), dimension(6), intent(out) :: B_Fac
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Symmetry relations among coefficients of the anisotropic temperature
    !!--++    factor.
    !!--++
    !!--++    Order for B is: B11 B22 B33 B12 B13 B23
    !!--++
    !!--++    B is considered as a 6-D vector and a single 6x6 matrix RB is constructed
    !!--++    in such a way as the matrix relation  B'ij = Sum{kh}[Rik Bkh Rjh] = Bij
    !!--++    is writen as B'= RB B = B  => (RB-I) B = 0
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Sym_B_Relations_OP(R,B_Ind,B_Fac)
       !---- Arguments ----!
       integer,dimension(3,3),      intent (in) :: R
       integer, dimension(6),       intent(out) :: B_Ind
       real(kind=cp), dimension(6), intent(out) :: B_Fac

       !---- Local variables ----!
       integer, dimension(6,6) :: rb
       integer                 :: i,j,k,nvar
       integer                 :: i1,i2

       !---- Init variables ----!
       err_symm=.false.
       ERR_Symm_Mess=" "

        rb(1,1)=r(1,1)*r(1,1)
        rb(1,2)=r(2,1)*r(2,1)
        rb(1,3)=r(3,1)*r(3,1)
        rb(1,4)=2*r(1,1)*r(2,1)
        rb(1,5)=2*r(1,1)*r(3,1)
        rb(1,6)=2*r(2,1)*r(3,1)

        rb(2,1)=r(1,2)*r(1,2)
        rb(2,2)=r(2,2)*r(2,2)
        rb(2,3)=r(3,2)*r(3,2)
        rb(2,4)=2*r(1,2)*r(2,2)
        rb(2,5)=2*r(1,2)*r(3,2)
        rb(2,6)=2*r(2,2)*r(3,2)

        rb(3,1)=r(1,3)*r(1,3)
        rb(3,2)=r(2,3)*r(2,3)
        rb(3,3)=r(3,3)*r(3,3)
        rb(3,4)=2*r(1,3)*r(2,3)
        rb(3,5)=2*r(1,3)*r(3,3)
        rb(3,6)=2*r(2,3)*r(3,3)

        rb(4,1)=r(1,1)*r(1,2)
        rb(4,2)=r(2,1)*r(2,2)
        rb(4,3)=r(3,1)*r(3,2)
        rb(4,4)=r(1,1)*r(2,2)+r(1,2)*r(2,1)
        rb(4,5)=r(1,1)*r(3,2)+r(3,1)*r(1,2)
        rb(4,6)=r(2,1)*r(3,2)+r(3,1)*r(2,2)

        rb(5,1)=r(1,1)*r(1,3)
        rb(5,2)=r(2,1)*r(2,3)
        rb(5,3)=r(3,1)*r(3,3)
        rb(5,4)=r(1,1)*r(2,3)+r(2,1)*r(1,3)
        rb(5,5)=r(1,1)*r(3,3)+r(1,3)*r(3,1)
        rb(5,6)=r(2,1)*r(3,3)+r(3,1)*r(2,3)

        rb(6,1)=r(1,2)*r(1,3)
        rb(6,2)=r(2,2)*r(2,3)
        rb(6,3)=r(3,2)*r(3,3)
        rb(6,4)=r(1,2)*r(2,3)+r(2,2)*r(1,3)
        rb(6,5)=r(1,2)*r(3,3)+r(3,2)*r(1,3)
        rb(6,6)=r(2,2)*r(3,3)+r(3,2)*r(2,3)

      !---- (Rb-1) Array ----!

       do i=1,6
          rb(i,i)=rb(i,i)-1
       end do

       !---- Init Output variables ----!
       b_ind=-1
       b_fac= 0.0
       nvar = 0

       !---- Free B parameters ----!
       do i=1,6
          if (all(rb(i,:)==0)) then
             b_ind(i)=i
             b_fac(i)=1.0
             nvar=nvar+1
          end if
       end do

       do j=1,6
          if (all(rb(:,j)==0)) then
             if (b_ind(j) < 0 ) then
                b_ind(j)=j
                b_fac(j)=1.0
                nvar=nvar+1
             end if
          end if
       end do

       !---- Zero B parameters ----!
       if (nvar /= 6) then
          do i=1,6
             j=count(rb(i,:)/=0)
             if (j /= 1) cycle
             do k=1,6
                if (rb(i,k)/=0 .and. b_ind(k) < 0) then
                   b_ind(k)=k
                   nvar=nvar+1
                   exit
                end if
             end do
          end do
       end if

       !---- Other relations ----!
       if (nvar /=6) then
          do i=1,6
             j=count(rb(i,:)/=0)
             if (j /= 2) cycle
             do j=1,6
                if (rb(i,j)/=0) then
                   i1=j
                   exit
                end if
             end do
             do k=i1+1,6
                if (rb(i,k)/=0) then
                   i2=k
                   exit
                end if
             end do

             if (b_ind(i1) < 0 .and. b_ind(i2) < 0) then
                b_ind(i1)=i1
                b_ind(i2)=i1
                b_fac(i1)=1.0
                b_fac(i2)=-real(rb(i,i1))/real(rb(i,i2))
                nvar=nvar+2
             else
                if (b_ind(i1) < 0) then
                   b_fac(i1)=-real(rb(i,i2))/real(rb(i,i1))
                   b_ind(i1)=i2
                else
                   b_fac(i2)=-real(rb(i,i1))/real(rb(i,i2))
                   b_ind(i2)=i1
                end if
                nvar=nvar+1
             end if
          end do
       end if

       if (any(b_ind==-1)) then
          err_symm=.true.
          ERR_Symm_Mess="Symmetry relations in B Factors are wrong! "
       end if

       return
    End Subroutine Sym_B_Relations_OP

    !!--++
    !!--++ Subroutine Sym_B_Relations_St(Symmcar,B_Ind,B_Fac)
    !!--++    character(len=*),            intent (in) :: Symmcar
    !!--++    integer, dimension(6),       intent(out) :: B_Ind
    !!--++    real(kind=cp), dimension(6), intent(out) :: B_Fac
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Symmetry relations among coefficients of the anisotropic temperature
    !!--++    factor.
    !!--++
    !!--++    Order for B is: B11 B22 B33 B12 B13 B23
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Sym_B_Relations_ST(Symmcar,B_Ind,B_Fac)
       !---- Arguments ----!
       character(len=*),            intent (in) :: Symmcar
       integer, dimension(6),       intent(out) :: B_Ind
       real(kind=cp), dimension(6), intent(out) :: B_Fac

       !---- Local variables ----!
       integer, dimension(3,3) :: a
       real(kind=cp), dimension(3)      :: t

       call read_xsym(symmcar,1,a,t)
       call sym_b_relations_op(a,b_ind,b_fac)

       return
    End Subroutine Sym_B_Relations_ST

    !!----
    !!---- Subroutine Sym_Prod_St(Syma,Symb,Symab,Modlat)
    !!----    character(len=*),         intent (in)  :: syma
    !!----    character(len=*),         intent (in)  :: symb
    !!----    character(len=len(syma)), intent (out) :: symab
    !!----    logical, optional,        intent (in)  :: modlat
    !!----
    !!----    Obtain the symbol/Op/Matrix+trans of the  symmetry operation corresponding
    !!----    to the product of two operators given in the Jone's Faithful(symbol)
    !!----    representation or in Symmetry Operator type.
    !!--<<
    !!----     Op_a =  (Sa,ta) ;  Op_b =  (Sb,tb)
    !!----
    !!----     Op_ab =  (Sa,ta) (Sb,tb)  = (Sa Sb,  Sa tb + ta)
    !!-->>
    !!----    If modlat=.true. or it is not present, the traslation
    !!----    part of the resulting operator is reduced to have components < 1.0
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Sym_Prod_St(Syma,Symb,Symab,Modlat)
       !---- Arguments ----!
       character(len=*),         intent (in) :: syma
       character(len=*),         intent (in) :: symb
       character(len=len(syma)), intent(out) :: symab
       logical,optional,         intent (in) :: modlat

       !--- Local variables ---!
       integer, dimension (3,3)      :: Sa,Sb
       real(kind=cp),dimension (3)   :: ta,tb

       call Read_Xsym(syma,1,Sa,ta)
       call Read_Xsym(symb,1,Sb,tb)

       if(present(modlat)) then
         if(.not. modlat) then
           ta = ta + matmul(real(Sa),tb)
         else
           ta = modulo_lat(ta + matmul(real(Sa),tb))
         end if
       else
         ta = modulo_lat(ta + matmul(real(Sa),tb))
       end if
       Sa = matmul(Sa,Sb)
       call Get_symsymb(Sa,ta,symab)

       return
    End Subroutine Sym_Prod_St

    !!----
    !!---- Subroutine Symmetry_Symbol(Op,Symb), (S,T,Symb), (Symm,Symb)
    !!----    type(Sym_Oper_type),         intent (in) :: Op
    !!----
    !!----    integer, dimension(3,3),     intent (in) :: S
    !!----    real(kind=cp), dimension(3), intent (in) :: t
    !!----
    !!----    character(len=*),            intent (in) :: Symm
    !!----
    !!----    character(len=*),            intent (out):: symb
    !!----
    !!----    Obtain the symbol of the symmetry element of the operator Op
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine Symmetry_Symbol_Op(Op,Symb)
    !!--++    type(Sym_Oper_type), intent (in)  :: Op
    !!--++    character(len=*),    intent (out) :: symb
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Obtain the symbol of the symmetry element of the operator Op
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Symmetry_Symbol_OP(Op,symb)
       !---- Arguments ----!
       type(Sym_Oper_Type),   intent (in)  :: Op
       character(len=*),     intent (out)  :: symb

       call symmetry_symbol_str(Op%Rot,Op%tr,symb)

       return
    End Subroutine Symmetry_Symbol_OP

    !!--++
    !!--++ Subroutine Symmetry_Symbol_Str(S,T,Symb)
    !!--++    integer, dimension(3,3),     intent( in) :: s
    !!--++    real(kind=cp), dimension(3), intent( in) :: t
    !!--++    character (len=*),           intent(out) :: symb
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Obtain the symbol of the symmetry element corresponding to operator (S,T)
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Symmetry_Symbol_Str(S,T,Symb)
       !---- Arguments ----!
       integer,       dimension(3,3),    intent( in) :: s
       real(kind=cp), dimension(3),      intent( in) :: t
       character (len=*),                intent(out) :: symb

       !---- Local variables ----!
       character (len=80)      :: carsym
       character (len=1)       :: signo
       integer                 :: i, n, npos
       integer, dimension(3)   :: ix1, ix2, ix3
       integer, dimension(3,3) :: w
       integer, dimension(3,3), parameter :: identidad = reshape((/1, 0, 0, &
                                                                   0, 1, 0, &
                                                                   0, 0, 1/),(/3,3/))
       real(kind=cp)                    :: rnum
       real(kind=cp), dimension(3)      :: t0,t1,t2,t3
       real(kind=cp), dimension(3)      :: x1,x2,x3
       real(kind=cp), dimension(3)      :: p0,p1,p2,p3
       real(kind=cp), dimension(3,3)    :: ww

       !---- Initialize ----!
       symb=" "
       n=axes_rotation(s)
       !t0=mod(t+10.0_cp,1.0_cp)  !Attempt to use the given translation
       t0=t                       !of the symmetry operator
       x1 =0.0
       ix1=0
       call Init_Err_Symm()

       select case (n)
          case (1) ! Traslation or identity
             if (sum(abs(t)) <= 3.0*eps_symm) then
                symb(1:1) ="1"
             else
                symb(1:3)="t ("
                npos=4
                call get_string_resolv(t0,x1,ix1,carsym)
                symb(npos:)=carsym(1:len_trim(carsym))//")"
             end if

          case (:-3) ! Rotoinversion
             !---- Inversion point ----!
             w=s-identidad
             call resolv_sist_3x3(w,-t0,t3,x3,ix3)

             !---- Axes rotation ----!
             w=matmul(s,s)-identidad
             t1=matmul(real(s),t0)+t0
             call resolv_sist_3x3(w,-t1,t2,x2,ix2)

             !---- Sense of rotation ----!
             !---- P0, P1 ----!
             p0=0.0
             p1=1.0
             do i=1,3
                if (ix2(i) == 0) then
                   p0(i)=t2(i)
                   p1(i)=t2(i)
                else
                   p0(i)=t2(i)+x2(i)*p0(ix2(i))
                   p1(i)=t2(i)+x2(i)*p1(ix2(i))
                end if
             end do

             !---- P2 ----!
             do i=1,3
                if (p1(i) > 0.0 ) exit
             end do
             select case (i)
                case (1)
                   p2(3)=0.5*p1(3)
                   p2(2)=0.7*p1(2)
                   p2(1)=-(p2(2)*p1(2) + p2(3)*p1(3))/p1(1)

                case (2)
                   p2(1)=0.5*p1(1)
                   p2(3)=0.7*p1(3)
                   p2(2)=-(p2(1)*p1(1) + p2(3)*p1(3))/p1(2)

                case (3)
                   p2(1)=0.5*p1(1)
                   p2(2)=0.7*p1(2)
                   p2(3)=-(p2(1)*p1(1) + p2(2)*p1(2))/p1(3)
             end select
             do i=1,3
                if (abs(p2(i) - p0(i)) <= eps_symm) p2(i)=p2(i)*p2(i)+0.5*real(i)
             end do

             !---- P3 ----!
             p3=matmul(real(s),p2)+t0
             ww(1,:)=p1-p0
             ww(2,:)=p2-p0
             ww(3,:)=p3-p0
             rnum=determ_a(ww)
             if (rnum > 0.0) then
                signo="-"
             else
                signo="+"
             end if

             !---- Determine the final symbol ----!
             write(unit=symb,fmt="(i2)") n
             symb=adjustl(symb)
             npos=len_trim(symb)
             npos=npos+1
             symb(npos:npos)=signo
             npos=npos+2
             call get_string_resolv(t2,x2,ix2,carsym)
             symb(npos:)=carsym(1:len_trim(carsym))//";"
             npos=len_trim(symb)+2
             call get_string_resolv(t3,x3,ix3,carsym)
             symb(npos:)=carsym(1:len_trim(carsym))

          case (-2)  ! Reflection or glide reflection
             t1=matmul(s,t0)+t0
             if (t1(1) <= eps_symm .and. t1(2) <= eps_symm .and. &
                 t1(3) <= eps_symm) then        ! Pure Reflection

                !----Mirror Plane ----!
                w=s-identidad
                call resolv_sist_3x3(w,-t0,t3,x3,ix3)
                symb(1:2)="m "
                npos=3
                call get_string_resolv(t3,x3,ix3,carsym)
                symb(npos:)=carsym(1:len_trim(carsym))
             else                          ! Glide Reflection
                t3=0.5*t1
                w=s-identidad
                t1=t0-t3
                call resolv_sist_3x3(w,-t1,t2,x2,ix2)

                !---- Determine the final symbol ----!
                symb(1:2)="g "

                !---- a: (1/2, 0, 0) ----!
                if ( (abs(t3(1) - 0.5) <= eps_symm) .and. (abs(t3(2)) <= eps_symm) .and. &
                     (abs(t3(3)) <= eps_symm) ) then
                   symb(1:2)="a "
                end if

                !---- b: (0, 1/2, 0) ----!
                if ( (abs(t3(2) - 0.5) <= eps_symm) .and. (abs(t3(1)) <= eps_symm) .and. &
                     (abs(t3(3)) <= eps_symm) ) then
                   symb(1:2)="b "
                end if

                !---- c: (0, 0, 1/2) ----!
                if ( (abs(t3(3) - 0.5) <= eps_symm) .and. (abs(t3(2)) <= eps_symm) .and. &
                     (abs(t3(1)) <= eps_symm) ) then
                   symb(1:2)="c "
                end if

                !---- n: ( 1/2, 1/2, 0); (0, 1/2, 1/2); (1/2, 0, 1/2) ----!
                !---- n: ( 1/2, 1/2, 1/2) ----!
                !---- n: (-1/2, 1/2, 1/2); (1/2, -1/2, 1/2); (1/2, 1/2, -1/2) ----!
                if ( (abs(t3(1) - 0.5) <= eps_symm) .and. (abs(t3(2) - 0.5) <= eps_symm) .and. &
                     (abs(t3(3)) <= eps_symm) ) then
                   symb(1:2)="n "
                end if
                if ( (abs(t3(2) - 0.5) <= eps_symm) .and. (abs(t3(3) - 0.5) <= eps_symm) .and. &
                     (abs(t3(1)) <= eps_symm) ) then
                   symb(1:2)="n "
                end if
                if ( (abs(t3(1) - 0.5) <= eps_symm) .and. (abs(t3(3) - 0.5) <= eps_symm) .and. &
                     (abs(t3(2)) <= eps_symm) ) then
                   symb(1:2)="n "
                end if
                if ( (abs(t3(1) - 0.5) <= eps_symm) .and. (abs(t3(2) - 0.5) <= eps_symm) .and. &
                     (abs(t3(3) - 0.5) <= eps_symm) ) then
                   symb(1:2)="n "
                end if
                if ( (abs(t3(1) + 0.5) <= eps_symm) .and. (abs(t3(2) - 0.5) <= eps_symm) .and. &
                     (abs(t3(3) - 0.5) <= eps_symm) ) then
                   symb(1:2)="n "
                end if
                if ( (abs(t3(1) - 0.5) <= eps_symm) .and. (abs(t3(2) + 0.5) <= eps_symm) .and. &
                     (abs(t3(3) - 0.5) <= eps_symm) ) then
                   symb(1:2)="n "
                end if
                if ( (abs(t3(1) - 0.5) <= eps_symm) .and. (abs(t3(2) - 0.5) <= eps_symm) .and. &
                     (abs(t3(3) + 0.5) <= eps_symm) ) then
                   symb(1:2)="n "
                end if

                !---- d: ( 1/4,+-1/4, 0); (0, 1/4,+-1/4); (+-1/4, 0, 1/4) ----!
                !---- d: ( 1/4, 1/4,+-1/4); (+-1/4, 1/4, 1/4); (1/4,+-1/4, 1/4) ----!
                !---- d: (-1/4, 1/4,+-1/4); (+-1/4,-1/4, 1/4); (1/4,+-1/4,-1/4) ----!
                p3=t3
                p3=mod(p3+10.0_cp,1.0_cp)
                do i=1,3
                   if (p3(i) > 0.5) p3(i)=p3(i) -1.0
                end do
                if ( (abs(p3(1) - 0.25) <= eps_symm) .and. (abs(abs(p3(2)) - 0.25) <= eps_symm) .and. &
                     (abs(p3(3)) <= eps_symm) ) then
                   symb(1:2)="d "
                end if
                if ( (abs(p3(2) - 0.25) <= eps_symm) .and. (abs(abs(p3(3)) - 0.25) <= eps_symm) .and. &
                     (abs(p3(1)) <= eps_symm) ) then
                   symb(1:2)="d "
                end if
                if ( (abs(p3(3) - 0.25) <= eps_symm) .and. (abs(abs(p3(1)) - 0.25) <= eps_symm) .and. &
                     (abs(p3(2)) <= eps_symm) ) then
                   symb(1:2)="d "
                end if
                if ( (abs(p3(1) - 0.25) <= eps_symm) .and. (abs(abs(p3(3)) - 0.25) <= eps_symm) .and. &
                     (abs(p3(2) - 0.25) <= eps_symm) ) then
                   symb(1:2)="d "
                end if
                if ( (abs(p3(2) - 0.25) <= eps_symm) .and. (abs(abs(p3(1)) - 0.25) <= eps_symm) .and. &
                     (abs(p3(3) - 0.25) <= eps_symm) ) then
                   symb(1:2)="d "
                end if
                if ( (abs(p3(1) - 0.25) <= eps_symm) .and. (abs(abs(p3(2)) - 0.25) <= eps_symm) .and. &
                     (abs(p3(3) - 0.25) <= eps_symm) ) then
                   symb(1:2)="d "
                end if
                if ( (abs(p3(1) + 0.25) <= eps_symm) .and. (abs(abs(p3(3)) - 0.25) <= eps_symm) .and. &
                     (abs(p3(2) - 0.25) <= eps_symm) ) then
                   symb(1:2)="d "
                end if
                if ( (abs(p3(2) + 0.25) <= eps_symm) .and. (abs(abs(p3(1)) - 0.25) <= eps_symm) .and. &
                     (abs(p3(3) - 0.25) <= eps_symm) ) then
                   symb(1:2)="d "
                end if
                if ( (abs(p3(3) + 0.25) <= eps_symm) .and. (abs(abs(p3(2)) - 0.25) <= eps_symm) .and. &
                     (abs(p3(1) - 0.25) <= eps_symm) ) then
                   symb(1:2)="d "
                end if
                npos=3

                !---- Glide Part ----!
                if ( symb(1:1) == "n" .or. symb(1:1) == "d" .or. &
                     symb(1:1) == "g" ) then
                   symb(npos:)="("
                   npos=npos+1
                   x1 =0.0
                   ix1=0
                   call get_string_resolv(t3,x1,ix1,carsym)
                   symb(npos:)=carsym(1:len_trim(carsym))//")"
                   npos=len_trim(symb)+2
                end if

                !---- Location of Glide Plane ----!
                call get_string_resolv(t2,x2,ix2,carsym)
                symb(npos:)=carsym(1:len_trim(carsym))
             end if

          case (-1)  ! Inversion
             t1=0.5*t0
             symb(1:3)="-1 "
             npos=4
             x1 =0.0
             ix1=0
             call get_string_resolv(t1,x1,ix1,carsym)
             symb(npos:)=carsym(1:len_trim(carsym))

          case (2:)  ! Rotation / Screw Rotation
             w=identidad
             t1=t0
             do i=1,n-1
                w=matmul(w,s)
                t1=t1+matmul(w,t0)
             end do
             if (abs(t1(1)) <= eps_symm .and. abs(t1(2)) <= eps_symm &
                 .and. abs(t1(3)) <= eps_symm) then              ! Pure rotation

                !---- Rotations axes ----!
                w=s-identidad
                call resolv_sist_3x3(w,-t0,t2,x2,ix2)

                !---- Sense of rotation ----!
                !---- P0, P1 ----!
                p0=0.0
                p1=1.0
                do i=1,3
                   if (ix2(i) == 0) then
                      p0(i)=t2(i)
                      p1(i)=t2(i)
                   else
                      p0(i)=t2(i)+x2(i)*p0(ix2(i))
                      p1(i)=t2(i)+x2(i)*p1(ix2(i))
                   end if
                end do

                !---- P2 ----!
                do i=1,3
                   if (p1(i) > 0.0 ) exit
                end do
                select case (i)
                   case (1)
                      p2(3)=0.5*p1(3)
                      p2(2)=0.7*p1(2)
                      p2(1)=-(p2(2)*p1(2) + p2(3)*p1(3))/p1(1)

                   case (2)
                      p2(1)=0.5*p1(1)
                      p2(3)=0.7*p1(3)
                      p2(2)=-(p2(1)*p1(1) + p2(3)*p1(3))/p1(2)

                   case (3)
                      p2(1)=0.5*p1(1)
                      p2(2)=0.7*p1(2)
                      p2(3)=-(p2(1)*p1(1) + p2(2)*p1(2))/p1(3)
                end select
                do i=1,3
                   if (abs(p2(i) - p0(i)) <= eps_symm) p2(i)=p2(i)*p2(i)+0.5*real(i)
                end do

                !---- P3 ----!
                p3=matmul(real(s),p2)+t0
                ww(1,:)=p1-p0
                ww(2,:)=p2-p0
                ww(3,:)=p3-p0

                rnum=determ_a(ww)
                if (rnum > 0.0) then
                   signo="+"
                else
                   signo="-"
                end if

                !---- Determine the final symbol ----!
                write(unit=symb,fmt="(i2)") n
                symb=adjustl(symb)
                npos=len_trim(symb)
                if ( n /= 2) then
                   npos=npos+1
                   symb(npos:)=signo
                end if
                npos=npos+2
                call get_string_resolv(t2,x2,ix2,carsym)
                symb(npos:)=carsym(1:len_trim(carsym))
             else                     ! Screw Rotation
                t3=(1.0/real(n))*t1
                w=s-identidad
                t1=t0-t3
                call resolv_sist_3x3(w,-t1,t2,x2,ix2)

                !---- Sense of rotation ----!
                !---- P0, P1 ----!
                p0=0.0
                p1=1.0
                do i=1,3
                   if (ix2(i) == 0) then
                      p0(i)=t2(i)
                      p1(i)=t2(i)
                   else
                      p0(i)=t2(i)+x2(i)*p0(ix2(i))
                      p1(i)=t2(i)+x2(i)*p1(ix2(i))
                   end if
                end do

                !---- P2 ----!
                do i=1,3
                   if (p1(i) > 0.0 ) exit
                end do
                select case (i)
                   case (1)
                      p2(3)=0.5*p1(3)
                      p2(2)=0.7*p1(2)
                      p2(1)=-(p2(2)*p1(2) + p2(3)*p1(3))/p1(1)

                   case (2)
                      p2(1)=0.5*p1(1)
                      p2(3)=0.7*p1(3)
                      p2(2)=-(p2(1)*p1(1) + p2(3)*p1(3))/p1(2)

                   case (3)
                      p2(1)=0.5*p1(1)
                      p2(2)=0.7*p1(2)
                      p2(3)=-(p2(1)*p1(1) + p2(2)*p1(2))/p1(3)
                end select
                do i=1,3
                   if (abs(p2(i) - p0(i)) <= eps_symm) p2(i)=p2(i)*p2(i)+0.5*real(i)
                end do

                !---- P3 ----!
                p3=matmul(real(s),p2)+t0
                ww(1,:)=p1-p0
                ww(2,:)=p2-p0
                ww(3,:)=p3-p0
                rnum=determ_a(ww)
                if (rnum > 0.0) then
                   signo="+"
                else
                   signo="-"
                end if

                !---- Determine the final symbol ----!
                write(unit=symb,fmt="(i2)") n
                symb=adjustl(symb)
                npos=len_trim(symb)
                if ( n /= 2) then
                   npos=npos+1
                   symb(npos:npos)=signo
                end if
                npos=npos+2

                !---- Screw Part ----!
                symb(npos:)="("
                npos=npos+1
                x1 =0.0
                ix1=0
                call get_string_resolv(t3,x1,ix1,carsym)
                symb(npos:)=carsym(1:len_trim(carsym))//")"
                npos=len_trim(symb)+2
                call get_string_resolv(t2,x2,ix2,carsym)
                symb(npos:)=carsym(1:len_trim(carsym))
             end if

       end select

       return
    End Subroutine Symmetry_Symbol_Str

    !!--++
    !!--++ Subroutine Symmetry_Symbol_Xyz(Symm,Symb)
    !!--++    character(len=*), intent (in)  :: symm
    !!--++    character(len=*), intent (out) :: symb
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Obtain the symbol of the  symmetry element corresponding
    !!--++    to an operator given in the Jone's Faithful representation
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Symmetry_Symbol_Xyz(Symm,Symb)
       !---- Arguments ----!
       character(len=*), intent (in)  :: symm
       character(len=*), intent (out) :: symb

       !--- Local variables ---!
       integer, dimension (3,3)      :: s
       real(kind=cp),    dimension (3)        :: t

       call Read_Xsym(symm,1,s,t)
       call symmetry_symbol_str(s,t,symb)

       return
    End Subroutine Symmetry_Symbol_Xyz

    !!----
    !!---- Subroutine Write_Bin_Spacegroup(SpG,Lun)
    !!----    type (Space_Group),  intent(in)  :: SpG   !  In -> SpaceGroup Variable
    !!----    integer,             intent(in)  :: Lun   !  In -> Logical unit of the file
    !!----
    !!----    Writing in file of logical unit "lun" the full structure of Space_Group_Type, SpG
    !!----    The file should have been opened with the access="stream" attribute. The procedure
    !!----    writes in the given order a series of bytes corresponding to the components of the
    !!----    type SpG. For reading back a Space Group structure from a binary file the subroutine
    !!----    Read_Bin_Spacegroup has to be used.
    !!----
    !!---- Update: February - 2013
    !!
    Subroutine Write_Bin_SpaceGroup(SpG,lun)
       !---- Arguments ----!
       type (Space_Group_Type),intent(in) :: SpG
       integer,                intent(in) :: lun

       !---- Local variables ----!
       integer                           :: i,j

       !---- Writing variables ----!
       write(unit=Lun) SpG%NumSpg,        &   ! Number of the Space Group
                       SpG%SPG_Symb,      &   ! Hermann-Mauguin Symbol
                       SpG%Hall,          &   ! Hall symbol
                       SpG%CrystalSys,    &   ! Crystal system
                       SpG%Laue,          &   ! Laue Class
                       SpG%PG,            &   ! Point group
                       SpG%Info,          &   ! Extra information
                       SpG%SG_setting,    &   ! Information about the SG setting (IT,KO,ML,ZA,Table,Standard,UnConventional)
                       SpG%Hexa,          &   !
                       SpG%SPG_lat,       &   ! Lattice type
                       SpG%SPG_latsy,     &   ! Lattice type Symbol
                       SpG%NumLat,        &   ! Number of lattice points in a cell
                       SpG%Latt_trans,    &   ! Lattice translations
                       SpG%Bravais,       &   ! String with Bravais symbol + translations
                       SpG%Centre,        &   ! Alphanumeric information about the center of symmetry
                       SpG%Centred,       &   ! Centric or Acentric [ =0 Centric(-1 no at origin),=1 Acentric,=2 Centric(-1 at origin)]
                       SpG%Centre_coord,  &   ! Fractional coordinates of the inversion centre
                       SpG%NumOps,        &   ! Number of reduced set of S.O.
                       SpG%Multip,        &   ! Multiplicity of the general position
                       SpG%Num_gen            ! Minimum number of operators to generate the Group
       do i=1,SpG%Multip
         write(unit=Lun) SpG%SymOp(i)%Rot,SpG%SymOp(i)%tr ! Symmetry operators
         write(unit=Lun) SpG%SymopSymb(i)                 ! Strings form of symmetry operators
       end do
       write(unit=Lun) SpG%R_Asym_Unit                    ! Asymmetric unit in real space
       write(unit=Lun) SpG%Wyckoff%num_orbit              ! Wyckoff Information
       do i=1,SpG%Wyckoff%num_orbit
         write(unit=Lun) SpG%Wyckoff%orbit(i)%norb
         write(unit=Lun) SpG%Wyckoff%orbit(i)%str_Orig
         do j=1,SpG%Wyckoff%orbit(i)%norb
           write(unit=Lun) SpG%Wyckoff%orbit(i)%str_orbit(j)
         end do
       end do
       return
    End Subroutine Write_Bin_SpaceGroup

    !!----
    !!---- Subroutine Write_Spacegroup(Spacegroup,Iunit,Full)
    !!----    type (Space_Group),  intent(in)  :: SpaceGroup !  In -> SpaceGroup Variable
    !!----    integer,  optional,  intent(in)  :: iunit      !  In -> Write information on Iunit
    !!----    logical,  optional,  intent(in)  :: full       !  In -> Full operator or not
    !!----
    !!----    Writing in file of logical unit "lun" the characteristics of
    !!----    the space group "SpaceG". Part of the information contained
    !!----    in  SpaceGroup may be undefined, depending on the tabulated
    !!----    nature of the item. If full=.true. is present the whole group
    !!----    is output including the symmetry symbol associated to each
    !!----    operator.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Write_SpaceGroup(SpaceGroup,Iunit,Full)
       !---- Arguments ----!
       type (Space_Group_Type),intent(in) :: SpaceGroup
       integer,   optional,    intent(in) :: iunit
       logical,   optional,    intent(in) :: full

       !---- Local variables ----!
       integer,  parameter                      :: max_lines=192
       character (len=100), dimension(max_lines):: texto
       character (len=40)                       :: aux
       integer                                  :: lun
       integer                                  :: i, nlines
       logical                                  :: print_latt

       !---- Initializing variables ----!
       lun=6
       print_latt=.false.
       if (present(iunit)) lun=iunit
       if (present(full))  print_latt=.true.

       !---- Printing ----!
       write(unit=lun,fmt="(/,/,a)")          "        Information on Space Group: "
       write(unit=lun,fmt="(a,/ )")           "        --------------------------- "
       write(unit=lun,fmt="(a,i3)")          " =>   Number of Space group: ", SpaceGroup%NumSpg
       write(unit=lun,fmt="(a,a)")           " =>  Hermann-Mauguin Symbol: ", trim(SpaceGroup%SPG_Symb)
       write(unit=lun,fmt="(a,a)")           " =>             Hall Symbol: ", trim(SpaceGroup%Hall)
       if(len_trim(SpaceGroup%gHall) > 1) &
       write(unit=lun,fmt="(a,a)")           " => Generalized Hall Symbol: ", trim(SpaceGroup%gHall)
       if(len_trim(SpaceGroup%info) > 1) &
       write(unit=lun,fmt="(a,a)")           " =>    Table Setting Choice: ", trim(SpaceGroup%info)
       write(unit=lun,fmt="(a,a)")           " =>            Setting Type: ", trim(SpaceGroup%SG_setting)

       write(unit=lun,fmt="(a,a)")           " =>          Crystal System: ", trim(SpaceGroup%CrystalSys)
       write(unit=lun,fmt="(a,a)")           " =>              Laue Class: ", trim(SpaceGroup%Laue)
       write(unit=lun,fmt="(a,a)")           " =>             Point Group: ", trim(SpaceGroup%Pg)

       write(unit=lun,fmt="(a,a)")           " =>         Bravais Lattice: ", trim(SpaceGroup%SPG_Lat)
       write(unit=lun,fmt="(a,a)")           " =>          Lattice Symbol: ", trim(SpaceGroup%SPG_Latsy)

       write(unit=lun,fmt="(a,i3)")          " =>  Reduced Number of S.O.: ", SpaceGroup%NumOps
       write(unit=lun,fmt="(a,i3)")          " =>    General multiplicity: ", SpaceGroup%Multip
       write(unit=lun,fmt="(a,a)")           " =>          Centrosymmetry: ", trim(SpaceGroup%Centre)
       write(unit=lun,fmt="(a,i3)")          " =>  Generators (exc. -1&L): ", SpaceGroup%num_gen
       write(unit=lun,fmt="(a,f6.3,a,f6.3)") " =>         Asymmetric unit: ", SpaceGroup%R_Asym_Unit(1,1), &
                                                                  " <= x <= ",SpaceGroup%R_Asym_Unit(1,2)
       write(unit=lun,fmt="(a,f6.3,a,f6.3)") "                             ", SpaceGroup%R_Asym_Unit(2,1), &
                                                                  " <= y <= ",SpaceGroup%R_Asym_Unit(2,2)
       write(unit=lun,fmt="(a,f6.3,a,f6.3)") "                             ", SpaceGroup%R_Asym_Unit(3,1), &
                                                                  " <= z <= ",SpaceGroup%R_Asym_Unit(3,2)

       if (SpaceGroup%centred == 0) then
          call Frac_Trans_1Dig(SpaceGroup%Centre_coord,texto(1))
          write(unit=lun,fmt="(a,a)")        " =>               Centre at: ", trim(texto(1))
       end if
       if (SpaceGroup%SPG_Lat == "Z" .or. print_latt) then
          texto(:) (1:100) = " "
          if (SpaceGroup%SPG_Lat == "Z") then
            write(unit=lun,fmt="(a,i3)")          " => Non-conventional Centring vectors:",SpaceGroup%Numlat
          else
            write(unit=lun,fmt="(a,i3)")          " => Centring vectors:",SpaceGroup%Numlat-1
          end if
          nlines=1
          do i=2,SpaceGroup%Numlat
             call Frac_Trans_1Dig(SpaceGroup%Latt_trans(:,i),aux)
             if (mod(i-1,2) == 0) then
                write(unit=texto(nlines)(51:100),fmt="(a,i2,a,a)") &
                                           " => Latt(",i-1,"): ",trim(aux)
                nlines=nlines+1
             else
                write(unit=texto(nlines)( 1:50),fmt="(a,i2,a,a)")  &
                                           " => Latt(",i-1,"): ",trim(aux)
             end if
          end do
          do i=1,nlines
             write(unit=lun,fmt="(a)") texto(i)
          end do
       end if

       !---- Symmetry Operators ----!
       if (present(full)) then
          write(unit=lun,fmt="(/,a,/)")        " => List of all Symmetry Operators and Symmetry Symbols"

          do i=1,SpaceGroup%Multip
             texto(1)=" "
             call Symmetry_Symbol(SpaceGroup%SymopSymb(i),texto(1))
             write(unit=lun,fmt="(a,i3,2a,t50,2a)") " => SYMM(",i,"): ",trim(SpaceGroup%SymopSymb(i)), &
                                                     "Symbol: ",trim(texto(1))
          end do

          !---- Wyckoff Information ----!
          call Write_Wyckoff(SpaceGroup%Wyckoff, SpaceGroup%SPG_Symb,lun)

       else
          write(unit=lun,fmt="(/,a)") " => List of S.O. without inversion and lattice centring translations"

          texto(:) (1:100) = " "
          nlines=1
          do i=1,SpaceGroup%NumOps
             if (mod(i,2) == 0) then
                write(unit=texto(nlines)(51:100),fmt="(a,i3,a,a)") &
                                           " => SYMM(",i,"): ",trim(SpaceGroup%SymopSymb(i))
                nlines=nlines+1
             else
                write(unit=texto(nlines)( 1:50),fmt="(a,i3,a,a)")  &
                                           " => SYMM(",i,"): ",trim(SpaceGroup%SymopSymb(i))
             end if
             if(nlines == max_lines) then
                texto(nlines)=trim(texto(nlines))//"   <= Maximum number of lines exhausted!"
                exit
             end if
          end do
          do i=1,nlines
             write(unit=lun,fmt="(a)") trim(texto(i))
          end do

       end if


       return
    End Subroutine Write_SpaceGroup

    !!----
    !!---- Subroutine Write_Sym(Lun,Indx,Sim,Tt,P_Mag,Mag)
    !!----    integer,                     intent(in) :: lun       !  In -> Logical unit of the file to write
    !!----    integer,dimension(3,3),      intent(in) :: sim       !  In -> Rotational part of the S.O.
    !!----    integer,                     intent(in) :: indx      !  In -> Ordinal of the current Symm.Operator
    !!----    real(kind=cp), dimension(3), intent(in) :: tt        !  In -> Translation part of the S.O.
    !!----    real(kind=cp),               intent(in) :: p_mag     !  In -> Magnetic phase of the magnetic S.O.
    !!----    logical,                     intent(in) :: mag       !  In -> .true. if it is a magnetic S.O.
    !!----
    !!----    Writing the reduced set of symmetry operators
    !!----    Logical hexa must be defined (valid for conventional bases)
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Write_Sym(Lun,Indx,Sim,Tt,P_Mag,Mag)
       !---- Arguments ----!
       integer,                     intent(in) :: lun,indx
       integer, dimension(3,3),     intent(in) :: sim
       real(kind=cp), dimension(3), intent(in) :: tt
       real(kind=cp),               intent(in) :: p_mag
       logical,                     intent(in) :: mag

       !---- Local variables ----!
       character (len=35)             :: symcod
       character (len=40)             :: Seitz_symb
       integer                        :: j,ihex,i1,i2,isl

       if (.not. hexa) then
          i1=1
          i2=24
       else
          i1=25
          i2=36
       end if
       call SearchOp(sim,i1,i2,Isl)
       call Get_SymSymb(sim,tt,Symcod)

       if (hexa) then
          j=abs(isl)-24
          if(Isl < 0) j=j+12
          call  Get_Seitz(j,tt,Seitz_symb)
          write(unit=lun,fmt="(i4,4(a,a))") indx," :: ",trim(IntSymD6h(j))," :: ", &
                                      trim(Kov_D6h(j))," :: ",trim(SymCod)," :: ",trim(Seitz_symb)

          if (mag) then
             j=abs(isl)-24
             ihex=2
             if (j < 0) then
                j=j+24
                ihex=1
             end if
             if (isl < 0) j=j+24/ihex
             write(unit=lun,fmt="(a,i2,a,a19,a,f12.4)") "      (",indx,"): ",  &
                                                  MAGmat(J+(ihex-1)*48)," MPhas: ",P_MAG
          end if

       else              ! No hexa
          j=abs(isl)
          if (isl < 0) j=j+24
          call  Get_Seitz(j,tt,Seitz_symb)
          write(unit=lun,fmt="(i4,4(a,a))") indx," :: ",trim(IntSymOh(j))," :: ", &
                                      trim(Kov_Oh(j))," :: ",trim(SymCod)," :: ",trim(Seitz_symb)
          if (mag) then
             j=abs(isl)
             if (isl < 0) j=j+24
             write(unit=lun,fmt="(a,i2,a,a13,a,f12.4)") "      (",indx,"): ",   &
                                                  MAGmat(J)," MPhas: ",P_MAG
          end if
       end if            ! End if(Hexa)

       return
    End Subroutine Write_Sym

    !!----
    !!---- Subroutine Write_SymTrans_Code(N,Tr,Code)
    !!----    integer,                    intent(in)  :: N
    !!----    real(kind=cp),dimension(3), intent(in)  :: Tr
    !!----    character (len=*),          intent(out) :: Code
    !!----
    !!----    Write the code string for reference the symmetry operator and the
    !!----    Traslation applied.
    !!--<<        _2.555     : N_Op = 2, Tr=( 0.0, 0.0, 0.0)
    !!----        _3.456     : N_Op = 3, Tr=(-1.0, 0.0, 1.0)
    !!-->>
    !!----
    !!---- Update: April - 2005
    !!
    Subroutine Write_SymTrans_Code(N,Tr,Code)
       !---- Arguments ----!
       integer,                    intent(in)  :: N
       real(kind=cp),dimension(3), intent(in)  :: Tr
       character (len=*),          intent(out) :: Code

       !---- Local Variables ----!
       character(len=3)      :: car
       integer, dimension(3) :: i

       Code=" "
       if (N <=0) return
       car="   "
       !---- Number of the Symmetry Operator ----!
       write(unit=car,fmt="(i3)") n
       car=adjustl(car)
       Code="_"//trim(car)
       car="   "
       !---- Traslation Part ----!
       i=5+nint(tr)
       if (any(i /= 5)) then
          write(unit=car(1:1),fmt="(i1)") i(1)
          write(unit=car(2:2),fmt="(i1)") i(2)
          write(unit=car(3:3),fmt="(i1)") i(3)
          code=trim(code)//"."//trim(car)
       else
          if(len_trim(code)==2 .and. code(2:2) == "1") Code=" "
       end if


       return
    End Subroutine Write_SymTrans_Code

    !!----
    !!---- Subroutine Write_Wyckoff(Wyckoff,Spg_Name,Lun, Sorting)
    !!----    type(wyckoff_type), intent(in) :: Wyckoff     !  In -> Wyckoff Type variable
    !!----    character(len=*),   intent(in) :: Spg_Name    !  In -> SpaceGroup Name
    !!----    integer,optional,   intent(in) :: Lun         !  In -> Unit to write the information
    !!----    logical, optional,  intent(in) :: Sorting     !  In -> .true. for sorting list
    !!----
    !!----    Print/Write the Wyckoff positions in Lun unit
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Write_Wyckoff(Wyckoff,Spg, Lun, Sorting)
       !---- Arguments ----!
       type(wyckoff_type), intent(in) :: wyckoff
       character(len=*),   intent(in) :: Spg
       integer, optional,  intent(in) :: Lun
       logical, optional,  intent(in) :: Sorting

       !---- Local variables ----!
       character(len=3)      :: carm
       character(len=12)     :: site
       integer               :: i,j,iunit
       integer,dimension(26) :: list,order
       character(len=*), dimension(26),parameter :: alphabet = (/  &
       "a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z"/)

       if (wyckoff%num_orbit == 0) return
       iunit=6
       if (present(lun)) iunit=lun

       !---- Sorting the final Wyckoff List ----!
       do i=1, wyckoff%num_orbit
          list(i)=wyckoff%orbit(i)%norb
          order(i)=i
       end do

       if (present(sorting)) then
          if (sorting) call sort(list,wyckoff%num_orbit,order)
       end if

       !---- Info ----!
       write(unit=iunit,fmt="(/,a)") " => Special Wyckoff Positions for "//trim(spg)
       write(unit=iunit,fmt="(a)") " "
       write(unit=iunit,fmt="(a)") "    Multp     Site        Representative Coordinates (centring translations excluded)"
       do i=wyckoff%num_orbit,1,-1
          write(unit=carm,fmt="(i3)") wyckoff%orbit(order(i))%multp
          site=alphabet(i)
          do j=1,wyckoff%orbit(order(i))%norb,3
             write(unit=iunit,fmt="(a,a,t15,a,t30,a,t50,a,t70,a)") "    ",&
                   carm,site,wyckoff%orbit(order(i))%str_orbit(j:j+2)
             carm=" "
             site=" "
          end do
          write(unit=iunit,fmt="(a)") " "
       end do

       return
     End Subroutine Write_Wyckoff

    !!----
    !!---- Subroutine Wyckoff_Orbit(Spacegroup,Wyckoffstr,N_Orbit,Orbitstr)
    !!----    type (Space_Group_Type),       intent( in) :: SpaceGroup !  In -> SpaceGroup Variable
    !!----    character(len=*),              intent( in) :: WyckoffStr !  In -> Representative of the Orbit
    !!----    integer,                       intent(out) :: N_Orbit    ! Out -> Number of Components in the Orbit
    !!----    character(len=*),dimension(:), intent(out) :: OrbitStr   ! Out -> Wyckoff Positions Strings
    !!----
    !!----    Calculation of the Wyckoff positions from the representative element
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Wyckoff_Orbit(SGrp,Wyckoff_Car, N, Wyckoff_Orb)
       !---- Arguments ----!
       type (Space_Group_Type),        intent( in) :: SGrp           !  In -> Space Group Information
       character(len=*),               intent( in) :: Wyckoff_Car    !  In -> Representative of the Orbit to calculate
       integer,                        intent(out) :: N              ! Out -> Number of components in the orbit
       character(len=*),dimension(:),  intent(out) :: Wyckoff_Orb    ! Out -> Wyckoff positions for this Orbit

       !---- Local Variables ----!
       logical                         :: delete
       character(len=40)               :: symb,symb2
       integer                         :: i,j,k,num
       integer,dimension(3,3)          :: w
       real(kind=cp),   dimension(3,3) :: w1
       real(kind=cp), dimension(3)     :: t,t1,t2

       Wyckoff_Orb=" "
       n=0
       if (len_trim(wyckoff_car) <= 0) return

       n=1
       wyckoff_orb(n)=adjustl(wyckoff_car)
       call Read_Xsym(wyckoff_car,1,w,t)
       err_symm=.false.

       num=sgrp%multip/sgrp%numlat

       do i=2,num
          w1=real(sgrp%symop(i)%rot)
          t1=sgrp%symop(i)%tr
          t1=applyso(sgrp%symop(i),t)
          t1=mod(t1+10.0_cp,1.0_cp)
          w1=matmul(w1,real(w))
          call Get_SymSymb(w1,t1,symb)
          delete=.false.
          do j=1,n
             if (symb == wyckoff_orb(j)) then
                delete=.true.
                exit
             end if
          end do
          if (delete) cycle

          !---- Lattice Contribution ----!
          do j=2,sgrp%numlat
             t2=t1+sgrp%latt_trans(:,j)
             t2=mod(t2+10.0_cp,1.0_cp)
             call Get_SymSymb(w1,t2,symb2)
             delete=.false.
             do k=1,n
                if (symb2 == wyckoff_orb(k)) then
                   delete=.true.
                   exit
                end if
             end do
             if (delete) exit
          end do
          if (delete) cycle

          n=n+1
          wyckoff_orb(n)=adjustl(symb)
       end do

       return
    End Subroutine Wyckoff_Orbit

    Subroutine Copy_NS_SpG_To_SpG(SpGN,SpG)
       !---- Arguments ----!
       type(NS_Space_Group_type), intent(in)    :: SpGN
       type(Space_Group_type),    intent(out)   :: SpG

       !---- Local Variables ----!
       logical              :: change
       integer              :: i,j,k
       real, dimension(3,3) :: w

       !> Init
       call init_err_symm()
       change=.true.

       !> Check if the copy is possible
       loop_1: do k=1,SpG%Multip
          w=SpGn%Symop(k)%Rot
          w=abs(w)*100.0
          do i=1,3
             do j=1,3
                if (w(i,j) > 0.5 .and. w(i,j) < 99.5) then
                   change=.false.
                   exit loop_1
                end if
             end do
          end do
       end do loop_1

       if (.not. change) then
          err_symm=.true.
          ERR_Symm_Mess="No copy was possible for SpgN to Spg "
          return
       end if

       SpG%NumSpg      = SpGn%NumSpg
       SpG%SPG_Symb    = SpGn%SPG_Symb
       SpG%Hall        = SpGn%Hall
       SpG%gHall       = SpGn%gHall
       SpG%CrystalSys  = SpGn%CrystalSys
       SpG%Laue        = SpGn%Laue
       SpG%PG          = SpGn%PG
       SpG%Info        = SpGn%Info
       SpG%SG_setting  = SpGn%SG_setting
       SpG%SPG_lat     = SpGn%SPG_lat
       SpG%SPG_latsy   = SpGn%SPG_latsy
       SpG%NumLat      = SpGn%NumLat
       if(allocated(SpG%Latt_Trans)) deallocate(SpG%Latt_Trans)
       allocate(SpG%Latt_Trans(3,SpG%NumLat))
       SpG%Latt_Trans  = SpGn%Latt_Trans
       SpG%Bravais     = SpGn%Bravais
       SpG%Centre      = SpGn%Centre
       SpG%Centred     = SpGn%Centred
       SpG%Centre_coord= SpGn%Centre_coord
       SpG%NumOps      = SpGn%NumOps
       SpG%Multip      = SpGn%Multip
       SpG%Num_gen     = SpGn%Num_gen
       if(allocated(SpG%SymopSymb)) deallocate(SpG%SymopSymb)
       allocate(SpG%SymopSymb(SpG%Multip))
       SpG%SymopSymb=SpGn%SymopSymb
       if(allocated(SpG%Symop)) deallocate(SpG%Symop)
       allocate(SpG%Symop(SpG%Multip))
       do i=1,SpG%Multip
         SpG%Symop(i)%Rot(:,:) = nint(SpGn%Symop(i)%Rot(:,:))
         SpG%Symop(i)%tr(:) =  SpGn%Symop(i)%tr(:)
       end do

       return
    End Subroutine Copy_NS_SpG_To_SpG

 End Module CFML_Crystallographic_Symmetry
