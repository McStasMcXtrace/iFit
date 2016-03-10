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
!!---- MODULE: CFML_Symmetry_Tables
!!----   INFO: Tabulated information on Crystallographic Symmetry
!!----
!!---- HISTORY
!!----    Update: 04/03/2011
!!----
!!----
!!---- DEPENDENCIES
!!--++    Use CFML_GlobalDeps,       only: cp
!!--++    Use CFML_String_Utilities, only: U_case
!!----
!!---- VARIABLES
!!----    BC_D6H
!!----    BC_OH
!!----    DEPMAT
!!----    ERR_SYMTAB
!!----    ERR_SYMTAB_MESS
!!--++    IT_SET                       [Private]
!!----    INTSYMD6H
!!----    INTSYMOH
!!----    KOV_D6H
!!----    KOV_OH
!!----    LATT
!!----    LAUE_CLASS
!!----    LTR_A
!!----    LTR_B
!!----    LTR_C
!!----    LTR_F
!!----    LTR_I
!!----    LTR_R
!!----    MAGMAT
!!----    ML_D6H
!!----    ML_OH
!!----    MOD6
!!----    POINT_GROUP
!!--++    SPG_GEN                      [Private]
!!----    SPGR_INFO_TYPE
!!----    SPGR_INFO
!!----    SYS_CRY
!!----    TABLE_EQUIV_TYPE
!!----    SYSTEM_EQUIV
!!----    WYCK_INFO_TYPE
!!----    WYCKOFF_INFO
!!----    X_D6H
!!----    X_OH
!!----    ZAK_D6H
!!----    ZAK_OH
!!----
!!---- PROCEDURES
!!----    Functions:
!!----
!!----    Subroutines:
!!----       GET_GENERATORS
!!----       REMOVE_SPGR_INFO
!!----       REMOVE_SYSTEM_EQUIV
!!----       REMOVE_WYCKOFF_INFO
!!--++       SET_IT_GEN                [Private]
!!----       SET_SPGR_INFO
!!----       SET_SYSTEM_EQUIV
!!----       SET_WYCKOFF_INFO
!!----
!!
 Module CFML_Symmetry_Tables
    !---- Use modules ----!
    Use CFML_GlobalDeps,        only: cp
    Use CFML_String_Utilities, only: U_Case

    !---- Variables ----!
    implicit none

    private

    !---- List of public subroutines ----!
    public :: get_generators
    public :: set_spgr_info, set_system_equiv, set_wyckoff_info
    public :: remove_spgr_info, remove_system_equiv, remove_wyckoff_info

    !---- List of private subroutines ----!
    private :: set_IT_gen

    !---- Definitions ----!

    !!----
    !!---- BC_D6H
    !!----    character (len=*), dimension(24), parameter, public :: BC_D6h
    !!----
    !!----    Bradley & Cracknell Notation
    !!----
    !!---- Update: February - 2005
    !!
    character (len=*), dimension(24), parameter, public  :: BC_D6h =(/                  &
       "  E  "," C+_3"," C-_3"," C_2 "," C-_6"," C+_6","C'_23","C'_21","C'_22", &
       "C`_23","C`_21","C`_22","  I  "," S-_6"," S+_6"," s_h "," S+_3"," S-_3", &
       " s_v3"," s_v1"," s_v2"," s_d3"," s_d1"," s_d2" /)

    !!----
    !!---- BC_OH
    !!----    character(len=*), dimension(48), parameter, public :: BC_Oh
    !!----
    !!----    Bradley & Cracknell Notation
    !!----
    !!---- Update: February - 2005
    !!
    character(len=*), dimension(48), parameter, public :: BC_Oh =(/             &
       "  E  "," C_2z"," C_2y"," C_2x","C+_31","C+_34","C+_33","C+_32","C-_31", &
       "C-_33","C-_32","C-_34"," C_2a"," C_2b","C-_4z","C+_4z","C-_4x"," C_2d", &
       " C_2f","C+_4x","C+_4y"," C_2c","C-_4y"," C_2e","  I  "," s_z "," s_y ", &
       " s_x ","S-_61","S-_64","S-_63","S-_62","S+_61","S+_63","S+_62","S+_64", &
       " s_da"," s_db","S+_4z","S-_4z","S+_4x"," s_dd"," s_df","S-_4x","S-_4y", &
       " s_dc","S+_4y"," s_de"  /)

    !!----
    !!---- DEPMAT
    !!----    character(len=*), dimension(72), parameter, public :: Depmat
    !!----
    !!----    Magnetic array
    !!----
    !!---- Update: February - 2005
    !!
    character(len=*), dimension(72), parameter, public :: Depmat = (/       &
       "( Dx, Dy, Dz)      ","(-Dx,-Dy, Dz)      ","(-Dx, Dy,-Dz)      ",   &
       "( Dx,-Dy,-Dz)      ","( Dz, Dx, Dy)      ","( Dz,-Dx,-Dy)      ",   &
       "(-Dz,-Dx, Dy)      ","(-Dz, Dx,-Dy)      ","( Dy, Dz, Dx)      ",   &
       "(-Dy, Dz,-Dx)      ","( Dy,-Dz,-Dx)      ","(-Dy,-Dz, Dx)      ",   &
       "( Dy, Dx,-Dz)      ","(-Dy,-Dx,-Dz)      ","( Dy,-Dx, Dz)      ",   &
       "(-Dy, Dx, Dz)      ","( Dx, Dz,-Dy)      ","(-Dx, Dz, Dy)      ",   &
       "(-Dx,-Dz,-Dy)      ","( Dx,-Dz, Dy)      ","( Dz, Dy,-Dx)      ",   &
       "( Dz,-Dy, Dx)      ","(-Dz, Dy, Dx)      ","(-Dz,-Dy,-Dx)      ",   &
       "(-Dx,-Dy,-Dz)      ","( Dx, Dy,-Dz)      ","( Dx,-Dy, Dz)      ",   &
       "(-Dx, Dy, Dz)      ","(-Dz,-Dx,-Dy)      ","(-Dz, Dx, Dy)      ",   &
       "( Dz, Dx,-Dy)      ","( Dz,-Dx, Dy)      ","(-Dy,-Dz,-Dx)      ",   &
       "( Dy,-Dz, Dx)      ","(-Dy, Dz, Dx)      ","( Dy, Dz,-Dx)      ",   &
       "(-Dy,-Dx, Dz)      ","( Dy, Dx, Dz)      ","(-Dy, Dx,-Dz)      ",   &
       "( Dy,-Dx,-Dz)      ","(-Dx,-Dz, Dy)      ","( Dx,-Dz,-Dy)      ",   &
       "( Dx, Dz, Dy)      ","(-Dx, Dz,-Dy)      ","(-Dz,-Dy, Dx)      ",   &
       "(-Dz, Dy,-Dx)      ","( Dz,-Dy,-Dx)      ","( Dz, Dy, Dx)      ",   &
       "( Dx   ,    Dy, Dz)","(   -Dy, Dx-Dy, Dz)","(-Dx+Dy,-Dx   , Dz)",   &
       "(-Dx   ,   -Dy, Dz)","(    Dy,-Dx+Dy, Dz)","( Dx-Dy, Dx   , Dz)",   &
       "(    Dy, Dx   ,-Dz)","( Dx-Dy,   -Dy,-Dz)","(-Dx   ,-Dx+Dy,-Dz)",   &
       "(   -Dy,-Dx   ,-Dz)","(-Dx+Dy,    Dy,-Dz)","( Dx   , Dx-Dy,-Dz)",   &
       "(-Dx   ,   -Dy,-Dz)","(    Dy,-Dx+Dy,-Dz)","( Dx-Dy, Dx   ,-Dz)",   &
       "( Dx   ,    Dy,-Dz)","(   -Dy, Dx-Dy,-Dz)","(-Dx+Dy,-Dx   ,-Dz)",   &
       "(   -Dy,-Dx   , Dz)","(-Dx+Dy,    Dy, Dz)","( Dx   , Dx-Dy, Dz)",   &
       "(    Dy, Dx   , Dz)","( Dx-Dy,   -Dy, Dz)","(-Dx   ,-Dx+Dy, Dz)"   /)

    !!----
    !!---- ERR_SYMTAB
    !!----    logical, public :: Err_Symtab
    !!----
    !!----    Logical Variable to indicate an error on this module.
    !!----
    !!---- Update: January - 2005
    !!
    logical, public :: ERR_Symtab=.false.

    !!----
    !!---- ERR_SYMTAB_MESS
    !!----    character(len=150), public :: ERR_SymTab_Mess
    !!----
    !!----    String containing information about the last error
    !!----
    !!---- Update: February - 2005
    !!
    character(len=150), public :: ERR_SymTab_Mess=" "

    !!--++
    !!--++ IT_SET
    !!--++    logical, private :: it_set=.false.
    !!--++
    !!--++    (PRIVATE)
    !!--++    Variable to test if generators have been set
    !!--++
    !!--++ Update: February - 2005
    !!
    logical, private :: it_set=.false.

    !!----
    !!---- INTSYMD6H
    !!----    character(len=* ), dimension(24), parameter, public:: IntSymD6h
    !!----
    !!----    International Symbols For Point Group Elements Of 6/mmm (D6h)
    !!----
    !!---- Update: February - 2005
    !!
    character(len=* ), dimension(24), parameter, public :: IntSymD6h =(/     &
       "  1           "," 3+ ( 0, 0, z)"," 3- ( 0, 0, z)","  2 ( 0, 0, z)",  &
       " 6- ( 0, 0, z)"," 6+ ( 0, 0, z)","  2 ( x, x, 0)","  2 ( x, 0, 0)",  &
       "  2 ( 0, y, 0)","  2 ( x,-x, 0)","  2 ( x,2x, 0)","  2 (2x, x, 0)",  &
       " -1           ","-3+ ( 0, 0, z)","-3- ( 0, 0, z)","  m ( x, y, 0)",  &
       "-6- ( 0, 0, z)","-6+ ( 0, 0, z)","  m ( x,-x, z)","  m ( x,2x, z)",  &
       "  m (2x, x, z)","  m ( x, x, z)","  m ( x, 0, z)","  m ( 0, y, z)"   /)

    !!----
    !!---- INTSYMOH
    !!----    character(len=* ), dimension(48), parameter, public :: IntSymOh
    !!----
    !!----    International Symbols For Point Group Elements Of M3M (Oh)
    !!----
    !!---- Update: February - 2005
    !!
    character(len=* ), dimension(48), parameter, public :: IntSymOh = (/     &
       "  1           ","  2 ( 0, 0, z)","  2 ( 0, y, 0)","  2 ( x, 0, 0)",  &
       " 3+ ( x, x, x)"," 3+ (-x, x,-x)"," 3+ ( x,-x,-x)"," 3+ (-x,-x, x)",  &
       " 3- ( x, x, x)"," 3- ( x,-x,-x)"," 3- (-x,-x, x)"," 3- (-x, x,-x)",  &
       "  2 ( x, x, 0)","  2 ( x,-x, 0)"," 4- ( 0, 0, z)"," 4+ ( 0, 0, z)",  &
       " 4- ( x, 0, 0)","  2 ( 0, y, y)","  2 ( 0, y,-y)"," 4+ ( x, 0, 0)",  &
       " 4+ ( 0, y, 0)","  2 ( x, 0, x)"," 4- ( 0, y, 0)","  2 (-x, 0, x)",  &
       " -1           ","  m ( x, y, 0)","  m ( x, 0, z)","  m ( 0, y, z)",  &
       "-3+ ( x, x, x)","-3+ (-x, x,-x)","-3+ ( x,-x,-x)","-3+ (-x,-x, x)",  &
       "-3- ( x, x, x)","-3- ( x,-x,-x)","-3- (-x,-x, x)","-3- (-x, x,-x)",  &
       "  m ( x,-x, z)","  m ( x, x, z)","-4- ( 0, 0, z)","-4+ ( 0, 0, z)",  &
       "-4- ( x, 0, 0)","  m ( x, y,-y)","  m ( x, y, y)","-4+ ( x, 0, 0)",  &
       "-4+ ( 0, y, 0)","  m (-x, y, x)","-4- ( 0, y, 0)","  m ( x, y, x)"   /)

    !!----
    !!---- KOV_D6H
    !!----    character(len=*), dimension(24), parameter, public :: Kov_D6h
    !!----
    !!----    Kovalev Notation
    !!----
    !!---- Update: February - 2005
    !!
    character(len=*), dimension(24), parameter, public :: Kov_d6h=(/       &
       " h1"," h3"," h5"," h4"," h6"," h2","h11"," h9"," h7"," h8","h12",  &
       "h10","h13","h15","h17","h16","h18","h14","h23",                    &
       "h21","h19","h20","h24","h22"/)

    !!----
    !!---- KOV_OH
    !!----    character(len=*), dimension(48), parameter, public :: Kov_Oh
    !!----
    !!----    Kovalev Notation
    !!----
    !!---- Update: February - 2005
    !!
    character(len=*), dimension(48), parameter, public :: Kov_Oh=(/               &
       " h1"," h4"," h3"," h2"," h9","h10","h12","h11"," h5"," h7"," h6"," h8",   &
       "h16","h13","h15","h14","h20","h18","h17","h19","h24","h23",               &
       "h22","h21","h25","h28","h27","h26","h33","h34","h36","h35",               &
       "h29","h31","h30","h32","h40","h37","h39","h38","h44","h42",               &
       "h41","h43","h48","h47","h46","h45"/)

    !!----
    !!---- LATT
    !!----    character(len=* ), dimension( 8) , parameter, public :: Latt
    !!----
    !!----    Lattice Traslations
    !!----
    !!---- Update: February - 2005
    !!
    character(len=* ), dimension( 8) , parameter, public  :: Latt =(/  &
       "  P: { 000 }                                       ",          &
       "  A: { 000;  0  1/2 1/2 }+                         ",          &
       "  B: { 000; 1/2  0  1/2 }+                         ",          &
       "  C: { 000; 1/2 1/2  0  }+                         ",          &
       "  I: { 000; 1/2 1/2 1/2 }+                         ",          &
       "  R: { 000; 2/3 1/3 1/3; 1/3 2/3 2/3   }+          ",          &
       "  F: { 000;  0  1/2 1/2; 1/2  0  1/2; 1/2 1/2  0 }+",          &
       "  Z: { 000;  Unconventional Z-centering vectors  }+"   /)

    !!----
    !!---- LAUE_CLASS
    !!----    character(len=*), dimension(16), parameter, public :: Laue_class
    !!----
    !!----    Laue symbols
    !!----
    !!---- Update: February - 2005
    !!
    character(len=*), dimension(16), parameter, public :: laue_class=(/ &
       "-1   ","2/m  ","mmm  ","4/m  ","4/mmm","-3 R ","-3m R","-3   ", &
       "-3m1 ","-31m ","6/m  ","6/mmm","m-3  ","m-3m ","m3   ","m3m  "/)

    !!----
    !!---- Litvin_point_op_label
    !!----    character(len=*), dimension(48), parameter, public :: Litvin_point_op_label
    !!----
    !!----    Symbols of point operators as given by Litvin (Non-hexagonal)
    !!----    The order corresponds to the Table given by Harold T. Stokes and Branton J. Campbell
    !!----
    !!---- Update: November - 2012, reordered according to the last tables 15/2/2016
    !!
    character(len=*), dimension(48), parameter, public :: Litvin_point_op_label=(/ &
       "1       ","2x      ","2y      ","2z      ","3xyz-1  ","3xy-z   ","3-xyz   ","3x-yz   ", &
       "3xyz    ","3x-yz-1 ","3xy-z-1 ","3-xyz-1 ","2-xy    ","4z      ","4z-1    ","2xy     ", &
       "2-yz    ","2yz     ","4x      ","4x-1    ","2-xz    ","4y-1    ","2xz     ","4y      ", &
       "-1      ","mx      ","my      ","mz      ","-3xyz-1 ","-3xy-z  ","-3-xyz  ","-3x-yz  ", &
       "-3xyz   ","-3x-yz-1","-3xy-z-1","-3-xyz-1","m-xy    ","-4z     ","-4z-1   ","mxy     ", &
       "m-yz    ","myz     ","-4x     ","-4x-1   ","m-xz    ","-4y-1   ","mxz     ","-4y     "/)

    !!----
    !!---- Litvin_point_op
    !!----    character(len=*), dimension(48), parameter, public :: Litvin_point_op
    !!----
    !!----    Jones Faithful symbols of point operators as given by Litvin (Non-hexagonal)
    !!----    The order corresponds to the Table given by Harold T. Stokes and Branton J. Campbell
    !!----
    !!---- Update: November - 2012, reordered according to the last tables 15/2/2016
    !!

    character(len=*), dimension(48), parameter, public :: Litvin_point_op=(/ &
       "x,y,z   ", "x,-y,-z ", "-x,y,-z ", "-x,-y,z ", "y,z,x   ",           &
       "y,-z,-x ", "-y,z,-x ", "-y,-z,x ", "z,x,y   ", "z,-x,-y ",           &
       "-z,x,-y ", "-z,-x,y ", "-y,-x,-z", "-y,x,z  ", "y,-x,z  ",           &
       "y,x,-z  ", "-x,-z,-y", "-x,z,y  ", "x,-z,y  ", "x,z,-y  ",           &
       "-z,-y,-x", "-z,y,x  ", "z,-y,x  ", "z,y,-x  ", "-x,-y,-z",           &
       "-x,y,z  ", "x,-y,z  ", "x,y,-z  ", "-y,-z,-x", "-y,z,x  ",           &
       "y,-z,x  ", "y,z,-x  ", "-z,-x,-y", "-z,x,y  ", "z,-x,y  ",           &
       "z,x,-y  ", "y,x,z   ", "y,-x,-z ", "-y,x,-z ", "-y,-x,z ",           &
       "x,z,y   ", "x,-z,-y ", "-x,z,-y ", "-x,-z,y ", "z,y,x   ",           &
       "z,-y,-x ", "-z,y,-x ", "-z,-y,x "/)


    !!----
    !!---- Litvin_point_op_hex_label
    !!----    character(len=*), dimension(24), parameter, public :: Litvin_point_op_hex_label
    !!----
    !!----    Symbols of point operators as given by Litvin (Hexagonal)
    !!----    The order corresponds to the Table given by Harold T. Stokes and Branton J. Campbell
    !!----
    !!---- Update: November - 2012, reordered according to the last tables 15/2/2016
    !!
    character(len=*), dimension(24), parameter, public :: Litvin_point_op_hex_label=(/ &
       "1    ","6z   ","3z   ","2z   ","3z-1 ","6z-1 ","2x   ","21   ",                &
       "2xy  ","22   ","2y   ","23   ","-1   ","-6z  ","-3z  ","mz   ",                &
       "-3z-1","-6z-1","mx   ","m1   ","mxy  ","m2   ","my   ","m3   "/)


    !!----
    !!---- Litvin_point_op_hex
    !!----    character(len=*), dimension(24), parameter, public :: Litvin_point_op_hex
    !!----
    !!----    Jones Faithful symbols of point operators as given by Litvin (Hexagonal)
    !!----    The order corresponds to the Table given by Harold T. Stokes and Branton J. Campbell
    !!----
    !!---- Update: November - 2012, reordered according to the last tables 15/2/2016
    !!

    character(len=*), dimension(24), parameter, public :: Litvin_point_op_hex=(/      &
       "x,y,z     ","x-y,x,z   ","-y,x-y,z  ","-x,-y,z   ","-x+y,-x,z ","y,-x+y,z  ", &
       "x-y,-y,-z ","x,x-y,-z  ","y,x,-z    ","-x+y,y,-z ","-x,-x+y,-z","-y,-x,-z  ", &
       "-x,-y,-z  ","-x+y,-x,-z","y,-x+y,-z ","x,y,-z    ","x-y,x,-z  ","-y,x-y,-z ", &
       "-x+y,y,z  ","-x,-x+y,z ","-y,-x,z   ","x-y,-y,z  ","x,x-y,z   ","y,x,z     "/)


    !!----
    !!---- LTR_A
    !!----    real(kind=cp), dimension(3,2), parameter, public :: Ltr_A
    !!----
    !!----    Lattice Traslations of type A
    !!----
    !!---- Update: February - 2005
    !!
    real(kind=cp), dimension(3,2), parameter, public :: Ltr_a =reshape ( (/0.0,0.0,0.0, 0.0,0.5,0.5/), (/3,2/) )

    !!----
    !!---- LTR_B
    !!----    real(kind=cp), dimension(3,2), parameter, public :: Ltr_B
    !!----
    !!----    Lattice Traslations of type B
    !!----
    !!---- Update: February - 2005
    !!
    real(kind=cp), dimension(3,2), parameter, public :: Ltr_b =reshape ( (/0.0,0.0,0.0, 0.5,0.0,0.5/), (/3,2/) )

    !!----
    !!---- LTR_C
    !!----    real(kind=cp), dimension(3,2), parameter, public :: Ltr_C
    !!----
    !!----    Lattice Traslations of type C
    !!----
    !!---- Update: February - 2005
    !!
    real(kind=cp), dimension(3,2), parameter, public :: Ltr_c =reshape ( (/0.0,0.0,0.0, 0.5,0.5,0.0/), (/3,2/) )

    !!----
    !!---- LTR_F
    !!----    real(kind=cp), dimension(3,4), parameter, public
    !!----
    !!----    Lattice Traslations of type F
    !!----
    !!---- Update: February - 2005
    !!
    real(kind=cp), dimension(3,4), parameter, public :: &
                   Ltr_f =reshape( (/0.0,0.0,0.0, 0.0,0.5,0.5, 0.5,0.0,0.5, 0.5,0.5,0.0 /),(/3,4/) )

    !!----
    !!---- LTR_I
    !!----    real(kind=cp), dimension(3,2), parameter, public :: Ltr_I
    !!----
    !!----    Lattice Traslations of type I
    !!----
    !!---- Update: February - 2005
    !!
    real(kind=cp), dimension(3,2), parameter, public :: Ltr_i =reshape ( (/0.0,0.0,0.0, 0.5,0.5,0.5/), (/3,2/) )

    !!----
    !!---- LTR_R
    !!----    real(kind=cp), dimension(3,3), parameter, public :: Ltr_R
    !!----
    !!----    Lattice Traslations of type R
    !!----
    !!---- Update: February - 2005
    !!
    real(kind=cp), dimension(3,3), parameter, public :: &
                   Ltr_r =reshape( (/0.0,0.0,0.0, 2.0/3.0,1.0/3.0,1.0/3.0,  1.0/3.0,2.0/3.0,2.0/3.0/),(/3,3/) )

    !!----
    !!---- MAGMAT
    !!----    character(len=* ), dimension(72), parameter, public :: Magmat
    !!----
    !!----    Magnetic array
    !!----
    !!---- Update: February - 2005
    !!
    character(len=* ), dimension(72), parameter, public :: Magmat = (/      &
       "( Mx, My, Mz)      ","(-Mx,-My, Mz)      ","(-Mx, My,-Mz)      ",   &
       "( Mx,-My,-Mz)      ","( Mz, Mx, My)      ","( Mz,-Mx,-My)      ",   &
       "(-Mz,-Mx, My)      ","(-Mz, Mx,-My)      ","( My, Mz, Mx)      ",   &
       "(-My, Mz,-Mx)      ","( My,-Mz,-Mx)      ","(-My,-Mz, Mx)      ",   &
       "( My, Mx,-Mz)      ","(-My,-Mx,-Mz)      ","( My,-Mx, Mz)      ",   &
       "(-My, Mx, Mz)      ","( Mx, Mz,-My)      ","(-Mx, Mz, My)      ",   &
       "(-Mx,-Mz,-My)      ","( Mx,-Mz, My)      ","( Mz, My,-Mx)      ",   &
       "( Mz,-My, Mx)      ","(-Mz, My, Mx)      ","(-Mz,-My,-Mx)      ",   &
       "(-Mx,-My,-Mz)      ","( Mx, My,-Mz)      ","( Mx,-My, Mz)      ",   &
       "(-Mx, My, Mz)      ","(-Mz,-Mx,-My)      ","(-Mz, Mx, My)      ",   &
       "( Mz, Mx,-My)      ","( Mz,-Mx, My)      ","(-My,-Mz,-Mx)      ",   &
       "( My,-Mz, Mx)      ","(-My, Mz, Mx)      ","( My, Mz,-Mx)      ",   &
       "(-My,-Mx, Mz)      ","( My, Mx, Mz)      ","(-My, Mx,-Mz)      ",   &
       "( My,-Mx,-Mz)      ","(-Mx,-Mz, My)      ","( Mx,-Mz,-My)      ",   &
       "( Mx, Mz, My)      ","(-Mx, Mz,-My)      ","(-Mz,-My, Mx)      ",   &
       "(-Mz, My,-Mx)      ","( Mz,-My,-Mx)      ","( Mz, My, Mx)      ",   &
       "( Mx   ,    My, Mz)","(   -My, Mx-My, Mz)","(-Mx+My,-Mx   , Mz)",   &
       "(-Mx   ,   -My, Mz)","(    My,-Mx+My, Mz)","( Mx-My, Mx   , Mz)",   &
       "(    My, Mx   ,-Mz)","( Mx-My,   -My,-Mz)","(-Mx   ,-Mx+My,-Mz)",   &
       "(   -My,-Mx   ,-Mz)","(-Mx+My,    My,-Mz)","( Mx   , Mx-My,-Mz)",   &
       "(-Mx   ,   -My,-Mz)","(    My,-Mx+My,-Mz)","( Mx-My, Mx   ,-Mz)",   &
       "( Mx   ,    My,-Mz)","(   -My, Mx-My,-Mz)","(-Mx+My,-Mx   ,-Mz)",   &
       "(   -My,-Mx   , Mz)","(-Mx+My,    My, Mz)","( Mx   , Mx-My, Mz)",   &
       "(    My, Mx   , Mz)","( Mx-My,   -My, Mz)","(-Mx   ,-Mx+My, Mz)"   /)

    !!----
    !!---- ML_D6H
    !!----    character(len=*), dimension(24), parameter, public:: ML_D6h
    !!----
    !!----    Miller & Love Notation
    !!----
    !!---- Update: February - 2005
    !!
    character(len=*), dimension(24), parameter, public :: ML_d6h=(/               &
       " 1"," 3"," 5"," 4"," 6"," 2"," 9"," 7","11","12","10"," 8","13","15","17",&
       "16","18","14","21","19","23","24","22","20"/)

    !!----
    !!---- ML_OH
    !!----     character(len=*), dimension(48), parameter, public :: ML_Oh
    !!----
    !!----     Miller & Love Notation
    !!----
    !!---- Update: February - 2005
    !!
    character(len=*), dimension(48), parameter, public :: ML_Oh=(/                &
       " 1"," 4"," 3"," 2"," 9","10","12","11"," 5"," 7"," 6"," 8","16","13","15",&
       "14","20","18","17","19","24","23","22","21","25","28","27","26","33","34",&
       "36","35","29","31","30","32","40","37","39","38","44","42","41","43","48",&
       "47","46","45"/)

    !!----
    !!---- MOD6
    !!----    Integer,  dimension(36,3,3), parameter, public :: Mod6
    !!----
    !!----    Matrix Types For Rotational Operators In Conventional Basis
    !!----    1->24 Oh, 25->36 D6h
    !!----
    !!---- Update: February - 2005
    !!
    Integer,  dimension(36,3,3), parameter, public :: Mod6 = reshape (  (/     &
       1,-1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,                   &
      -1, 1, 0, 0, 0, 0, 1, 0,-1,-1, 0, 1, 0, 1,-1, 0,-1, 1,                   &
       0, 0, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0, 1,-1,-1, 1, 0, 0,                   &
       0, 0, 0, 0, 0, 0, 0, 1,-1, 0,-1, 1, 1, 0,-1,-1, 0, 1,                   &
       0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0, 0, 0,                   &
       0, 0,-1, 1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                   &
       0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0,                   &
       0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 1,-1, 1,-1, 0,-1, 1, 0,                   &
       1,-1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                   &
       0, 0, 1,-1, 1,-1, 1,-1, 0,-1, 1, 0, 0,-1, 1, 0, 1,-1,                   &
       0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1,                   &
      -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                   &
       0, 0, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                   &
       0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                   &
       0, 0, 0, 0, 0, 0, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 1, 1,                   &
      -1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                   &
       1, 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1, 1, 1, 0, 0,                   &
       0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1 /), (/36,3,3/) )

    !!----
    !!---- POINT_GROUP
    !!----    character(len=*), dimension(39), parameter, public :: Point_group
    !!----
    !!----    Point Group Symbols
    !!----
    !!---- Update: July - 2014: added m3 and m3m for compatibility with Laue_class
    !!
    character(len=*), dimension(41), parameter, public :: point_group=(/  &
       "1    ","-1   ","2    ","m    ","2/m  ","222  ","mm2  ","m2m  ",   &
       "2mm  ","mmm  ","4    ","-4   ","4/m  ","422  ","4mm  ","-42m ",   &
       "-4m2 ","4/mmm","3    ","-3   ","32   ","3m   ","-3m  ","312  ",   &
       "31m  ","-31m ","6    ","-6   ","6/m  ","622  ","6mm  ","-62m ",   &
       "-6m2 ","6/mmm","23   ","m-3  ","432  ","-43m ","m-3m ","m3   ",   &
       "m3m  "/)

    !!--++
    !!--++ SPG_GEN
    !!--++    character(len=120), private, dimension(230) :: spg_gen
    !!--++
    !!--++    (PRIVATE)
    !!--++    Variable to hold the generators of all space groups in the standard setting
    !!--++
    !!--++ Update: February - 2005
    !!
    character(len=120), private, dimension(230) :: spg_gen

    !!----
    !!---- TYPE :: SPGR_INFO_TYPE
    !!--..
    !!---- Type, public :: Spgr_Info_Type
    !!----    integer                 :: N           ! Number of the Spacegroup
    !!----    character (len=12)      :: HM          ! Hermann-Mauguin
    !!----    character (len=16)      :: Hall        ! Hall
    !!----    integer                 :: Laue        ! Laue Group
    !!----    integer                 :: Pg          ! Point group
    !!----    integer, dimension(6)   :: Asu         ! Asymmetric unit * 24
    !!----    character (len= 5)      :: Inf_extra   ! Extra information
    !!---- End Type Spgr_Info_Type
    !!----
    !!----    Definition for General Info about Space Groups
    !!----
    !!---- Update: February - 2005
    !!
    Type, public :: Spgr_Info_Type
       integer                 :: N
       character (len=12)      :: HM
       character (len=16)      :: Hall
       integer                 :: Laue
       integer                 :: Pg
       integer, dimension(6)   :: Asu
       character (len= 5)      :: Inf_Extra
    End Type Spgr_Info_Type

    !!----
    !!---- SPGR_INFO
    !!----    Type(Spgr_Info_Type), allocatable, dimension(:), public :: Spgr_info
    !!----
    !!----    General Info about Space Groups
    !!----    Present dimension: 612
    !!----
    !!---- Update: February - 2005
    !!
    Type(Spgr_Info_Type), allocatable, dimension(:), public :: Spgr_Info

    !!----
    !!---- SYS_CRY
    !!----    character(len=* ), dimension(7) , parameter, public :: Sys_cry
    !!----
    !!----    System Type
    !!----
    !!---- Update: February - 2005
    !!
    character(len=* ), dimension(7) , parameter, public:: sys_cry =(/  &
       "Triclinic   ","Monoclinic  ","Orthorhombic","Tetragonal  ",    &
       "Trigonal    ","Hexagonal   ","Cubic       " /)

    !!----
    !!---- TYPE :: TABLE_EQUIV_TYPE
    !!--..
    !!---- Type, public :: Table_Equiv_Type
    !!----    character(len= 6)      :: SC     ! Schoenflies
    !!----    character(len=17)      :: ML     ! Miller & Love
    !!----    character(len=18)      :: KO     ! Kovalev
    !!----    character(len=32)      :: BC     ! Bradley & Cracknell
    !!----    character(len=18)      :: ZA     ! Zak
    !!---- End Type Table_Equiv_Type
    !!----
    !!----    Definition for Equivalences on a Table
    !!----
    !!---- Update: February - 2005
    !!
    Type, public :: Table_Equiv_Type
       character(len= 6)      :: SC                ! Schoenflies
       character(len=17)      :: ML                ! Miller & Love
       character(len=18)      :: KO                ! Kovalev
       character(len=32)      :: BC                ! Bradley & Cracknell
       character(len=18)      :: ZA                ! Zak
    End Type Table_Equiv_Type

    !!----
    !!---- SYSTEM_EQUIV
    !!----    Type(Table_Equiv_Type), allocatable, dimension(:), public :: System_Equiv
    !!----
    !!----    General Info about Space Groups
    !!----
    !!---- Update: February - 2005
    !!
    Type(Table_Equiv_Type), allocatable, dimension(:), public :: System_Equiv

    !!----
    !!---- TYPE :: WYCK_INFO_TYPE
    !!--..
    !!---- Type, public :: Wyck_Info_Type
    !!----    character (len=12)                :: HM          ! Hermann-Mauguin
    !!----    integer                           :: Norbit      ! Number of orbites
    !!----    character (len= 15),dimension(24) :: Corbit      ! Generator of the orbit
    !!---- End Type Wyck_Info_Type
    !!----
    !!----    Definition for Wyckoff Positions acording to IT
    !!----
    !!---- Update: February - 2005
    !!
    Type, public :: Wyck_Info_Type
       character (len=12)               :: HM
       integer                          :: Norbit
       character (len=15),dimension(26) :: Corbit
    End Type Wyck_Info_Type

    !!----
    !!---- WYCKOFF_INFO
    !!----    Type(Wyck_Info_Type), allocatable, dimension(:), public :: Wyckoff_info
    !!----
    !!----    General Info about Wyckoff Positions on IT
    !!----    Present dimension:
    !!----
    !!---- Update: February - 2005
    !!
    Type(Wyck_Info_Type), allocatable, dimension(:), public :: Wyckoff_Info

    !!----
    !!---- X_D6H
    !!----    character(len=* ), dimension(24), parameter, public:: X_D6h
    !!----
    !!---- Update: February - 2005
    !!
    character(len=* ), dimension(24), parameter, public   :: X_d6h = (/      &
       "( x  ,   y, z)","(  -y, x-y, z)","(-x+y,-x  , z)","(-x  ,  -y, z)",  &
       "(   y,-x+y, z)","( x-y, x  , z)","(   y, x  ,-z)","( x-y,  -y,-z)",  &
       "(-x  ,-x+y,-z)","(  -y,-x  ,-z)","(-x+y,   y,-z)","( x  , x-y,-z)",  &
       "(-x  ,  -y,-z)","(   y,-x+y,-z)","( x-y, x  ,-z)","( x  ,   y,-z)",  &
       "(  -y, x-y,-z)","(-x+y,-x  ,-z)","(  -y,-x  , z)","(-x+y,   y, z)",  &
       "( x  , x-y, z)","(   y, x  , z)","( x-y,  -y, z)","(-x  ,-x+y, z)"   /)

    !!----
    !!---- X_OH
    !!----    character(len=* ), dimension(48), parameter, public :: X_oh
    !!----
    !!---- Update: February - 2005
    !!
    character(len=* ), dimension(48), parameter, public  :: X_oh = (/                 &
       "( x, y, z)","(-x,-y, z)","(-x, y,-z)","( x,-y,-z)","( z, x, y)","( z,-x,-y)", &
       "(-z,-x, y)","(-z, x,-y)","( y, z, x)","(-y, z,-x)","( y,-z,-x)","(-y,-z, x)", &
       "( y, x,-z)","(-y,-x,-z)","( y,-x, z)","(-y, x, z)","( x, z,-y)","(-x, z, y)", &
       "(-x,-z,-y)","( x,-z, y)","( z, y,-x)","( z,-y, x)","(-z, y, x)","(-z,-y,-x)", &
       "(-x,-y,-z)","( x, y,-z)","( x,-y, z)","(-x, y, z)","(-z,-x,-y)","(-z, x, y)", &
       "( z, x,-y)","( z,-x, y)","(-y,-z,-x)","( y,-z, x)","(-y, z, x)","( y, z,-x)", &
       "(-y,-x, z)","( y, x, z)","(-y, x,-z)","( y,-x,-z)","(-x,-z, y)","( x,-z,-y)", &
       "( x, z, y)","(-x, z,-y)","(-z,-y, x)","(-z, y,-x)","( z,-y,-x)","( z, y, x)"  /)

    !!----
    !!---- ZAK_D6H
    !!----    character (len=*), dimension(24), parameter, public :: Zak_D6h
    !!----
    !!----    Zak Notation
    !!----
    !!---- Update: February - 2005
    !!
    character (len=*), dimension(24), parameter, public :: Zak_D6h =(/          &
       "   E   "," C(z)_3","C(2z)_3","  C_2  ","C(5z)_6"," C(z)_6","  U(xy)",   &
       "  U(x) ","  U(y) ","  U(3) ","  U(2) ","  U(1) ","   I   ","S(5z)_6",   &
       " S(z)_6","  s(z) "," S(z)_3","S(2z)_3"," s(xy) ","  s(x) ","  s(y) ",   &
       "  s(3) ","  s(2) ","  s(1) " /)

    !!----
    !!---- ZAK_OH
    !!----    character(len=* ), dimension(48), parameter, public :: Zak_Oh
    !!----
    !!----    Zak Notation
    !!----
    !!---- Update: February - 2005
    !!
    character(len=* ), dimension(48), parameter, public :: Zak_Oh =(/           &
       "     E     ","    U(z)   ","    U(y)   ","    U(x)   ","  C(xyz)_3 ",   &
       " C(-xy-z)_3"," C(x-y-z)_3"," C(-x-yz)_3"," C(2xyz)_3 ","C(2x-y-z)_3",   &
       " C(2x-yz)_3","C(-2xy-z)_3","    U(xy)  ","   U(-xy)  ","   C(3z)_4 ",   &
       "   C(z)_4  ","   C(3x)_4 ","    U(yz)  ","   U(y-z)  ","   C(x)_4  ",   &
       "   C(y)_4  ","    U(xz)  ","   C(3y)_4 ","   U(x-z)  ","      I    ",   &
       "    s(z)   ","    s(y)   ","    s(x)   "," S(5xyz)_6 ","S(-5xy-z)_6",   &
       "S(5x-y-z)_6","S(-5x-yz)_6","  S(xyz)_6 "," S(x-y-z)_6"," S(-x-yz)_6",   &
       " S(-xy-z)_6","    s(xy)  ","   s(-xy)  ","   S(z)_4  ","  S(3z)_4  ",   &
       "   S(x)_4  ","    s(yz)  ","   s(y-z)  ","  S(3x)_4  ","  S(3y)_4  ",   &
       "    s(xz)  ","   S(y)_4  ","   s(x-z)  " /)

 Contains

    !---------------------!
    !---- Subroutines ----!
    !---------------------!

    !!----
    !!---- Subroutine Get_Generators(Spg,Gener)
    !!----    character (len=*), intent(in)  :: spg     !  In -> Hermann_Mauguin symbol or number of S.Group
    !!----    character (len=*), intent(out) :: gener   ! Out -> String with all generators
    !!----
    !!----    Provides the string "gener" containing the list of the generators
    !!----    (as given in the IT Crystallography) corresponding to the space group
    !!----    of symbol "spg". In "spg" the Hermann-Mauguin symbol or the number of the
    !!----    space group should be given. The calling program is responsible of decoding
    !!----    the string "gener". Generator are given in the Jone's Faithful notation and
    !!----    the separator is the symbol ";". An example, corresponding to the space
    !!----    group "R 3 c" is  gener = " x+1/3,y+2/3,z+2/3; -y,x-y,z; -y,-x,z+1/2"
    !!----    The variable is the string contained between the quotes.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Get_Generators(Spg,Gener)
       !---- Arguments ----!
       character (len=*), intent(in)  :: spg
       character (len=*), intent(out) :: gener

       !----  Local variables ----!
       logical                 :: ok
       integer                 :: i, ier, numg
       character(len=len(spg)) :: symb,sp

       err_symtab=.false.
       if (.not. it_set) call set_IT_gen()
       ok=.false.

       read(unit=spg,fmt=*,iostat=ier) numg
       if (ier == 0) then
          if (numg > 0 .and. numg <= 230) then
             gener=spg_gen(numg)(12:)
             ok=.true.
          else
             gener=spg_gen(1)(12:)
          end if
       else
          symb=u_case(spg)
          do i=1,230
             sp=u_case(spg_gen(i)(1:10))
             if (symb == sp) then
                gener=spg_gen(i)(12:)
                ok=.true.
                exit
             end if
          end do
       end if

       if (.not. ok) then
          err_symtab=.true.
          ERR_SymTab_Mess=" Error in the symbol or number of the space group"
       end if

       return
    End Subroutine Get_Generators

    !!----
    !!---- Subroutine Remove_Spgr_Info()
    !!----
    !!----    Deallocating SPGR_INFO Data
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Remove_Spgr_Info()

       if (allocated(spgr_info)) deallocate(spgr_info)

       return
    End Subroutine Remove_Spgr_Info

    !!----
    !!---- Subroutine Remove_System_Equiv()
    !!----
    !!----    Deallocating SPGR_INFO Data
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Remove_System_Equiv()

       if (allocated(System_Equiv)) deallocate(System_Equiv)

       return
    End Subroutine Remove_System_Equiv

    !!----
    !!---- Subroutine Remove_Wyckoff_Info()
    !!----
    !!----    Deallocating WYCKOFF_INFO Data
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Remove_Wyckoff_Info()

       if (allocated(wyckoff_info)) deallocate(wyckoff_info)

       return
    End Subroutine Remove_Wyckoff_Info

    !!--++
    !!--++ Subroutine Set_It_Gen()
    !!--++
    !!--++    (PRIVATE)
    !!--++    Fills the components of the Spg_Gen character variable
    !!--++    Called once by the public subroutine Get_Generators
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Set_It_Gen()

       spg_gen(  1) =  "P 1       : x,y,z "
       spg_gen(  2) =  "P -1      : -x,-y,-z "
       spg_gen(  3) =  "P 2       : -x,y,-z "
       spg_gen(  4) =  "P 21      : -x,y+1/2,-z "
       spg_gen(  5) =  "C 2       : x+1/2,y+1/2,z; -x,y,-z "
       spg_gen(  6) =  "P m       : x,-y,z "
       spg_gen(  7) =  "P c       : x,-y,z+1/2 "
       spg_gen(  8) =  "C m       : x+1/2,y+1/2,z; x,-y,z "
       spg_gen(  9) =  "C c       : x+1/2,y+1/2,z; x,-y,z+1/2 "
       spg_gen( 10) =  "P 2/m     : -x,y,-z; -x,-y,-z "
       spg_gen( 11) =  "P 21/m    : -x,y+1/2,-z; -x,-y,-z "
       spg_gen( 12) =  "C 2/m     : x+1/2,y+1/2,z; -x,y,-z; -x,-y,-z "
       spg_gen( 13) =  "P 2/c     : -x,y,-z+1/2; -x,-y,-z "
       spg_gen( 14) =  "P 21/c    : -x,y+1/2,-z+1/2; -x,-y,-z "
       spg_gen( 15) =  "C 2/c     : x+1/2,y+1/2,z; -x,y,-z+1/2; -x,-y,-z "
       spg_gen( 16) =  "P 2 2 2   : -x,-y,z; -x,y,-z "
       spg_gen( 17) =  "P 2 2 21  : -x,-y,z+1/2; -x,y,-z+1/2 "
       spg_gen( 18) =  "P 21 21 2 : -x,-y,z; -x+1/2,y+1/2,-z "
       spg_gen( 19) =  "P 21 21 21: -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2 "
       spg_gen( 20) =  "C 2 2 21  : x+1/2,y+1/2,z; -x,-y,z+1/2; -x,y,-z+1/2 "
       spg_gen( 21) =  "C 2 2 2   : x+1/2,y+1/2,z; -x,-y,z; -x,y,-z "
       spg_gen( 22) =  "F 2 2 2   : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; -x,y,-z "
       spg_gen( 23) =  "I 2 2 2   : x+1/2,y+1/2,z+1/2; -x,-y,z; -x,y,-z "
       spg_gen( 24) =  "I 21 21 21: x+1/2,y+1/2,z+1/2; -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2 "
       spg_gen( 25) =  "P m m 2   : -x,-y,z; x,-y,z "
       spg_gen( 26) =  "P m c 21  : -x,-y,z+1/2; x,-y,z+1/2 "
       spg_gen( 27) =  "P c c 2   : -x,-y,z; x,-y,z+1/2 "
       spg_gen( 28) =  "P m a 2   : -x,-y,z; x+1/2,-y,z "
       spg_gen( 29) =  "P c a 21  : -x,-y,z+1/2; x+1/2,-y,z "
       spg_gen( 30) =  "P n c 2   : -x,-y,z; x,-y+1/2,z+1/2 "
       spg_gen( 31) =  "P m n 21  : -x+1/2,-y,z+1/2; x+1/2,-y,z+1/2 "
       spg_gen( 32) =  "P b a 2   : -x,-y,z; x+1/2,-y+1/2,z "
       spg_gen( 33) =  "P n a 21  : -x,-y,z+1/2; x+1/2,-y+1/2,z "
       spg_gen( 34) =  "P n n 2   : -x,-y,z; x+1/2,-y+1/2,z+1/2 "
       spg_gen( 35) =  "C m m 2   : x+1/2,y+1/2,z; -x,-y,z; x,-y,z "
       spg_gen( 36) =  "C m c 21  : x+1/2,y+1/2,z; -x,-y,z+1/2; x,-y,z+1/2 "
       spg_gen( 37) =  "C c c 2   : x+1/2,y+1/2,z; -x,-y,z; x,-y,z+1/2 "
       spg_gen( 38) =  "A m m 2   : x,y+1/2,z+1/2; -x,-y,z; x,-y,z "
       spg_gen( 39) =  "A b m 2   : x,y+1/2,z+1/2; -x,-y,z; x,-y+1/2,z "
       spg_gen( 40) =  "A m a 2   : x,y+1/2,z+1/2; -x,-y,z; x+1/2,-y,z "
       spg_gen( 41) =  "A b a 2   : x,y+1/2,z+1/2; -x,-y,z; x+1/2,-y+1/2,z "
       spg_gen( 42) =  "F m m 2   : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; x,-y,z "
       spg_gen( 43) =  "F d d 2   : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; x+1/4,-y+1/4,z+1/4 "
       spg_gen( 44) =  "I m m 2   : x+1/2,y+1/2,z+1/2; -x,-y,z; x,-y,z "
       spg_gen( 45) =  "I b a 2   : x+1/2,y+1/2,z+1/2; -x,-y,z; x+1/2,-y+1/2,z "
       spg_gen( 46) =  "I m a 2   : x+1/2,y+1/2,z+1/2; -x,-y,z; x+1/2,-y,z "
       spg_gen( 47) =  "P m m m   : -x,-y,z; -x,y,-z; -x,-y,-z "
       spg_gen( 48) =  "P n n n   : -x+1/2,-y+1/2,z; -x+1/2,y,-z+1/2; -x,-y,-z "
       spg_gen( 49) =  "P c c m   : -x,-y,z; -x,y,-z+1/2; -x,-y,-z "
       spg_gen( 50) =  "P b a n   : -x+1/2,-y+1/2,z; -x+1/2,y,-z; -x,-y,-z "
       spg_gen( 51) =  "P m m a   : -x+1/2,-y,z; -x,y,-z; -x,-y,-z "
       spg_gen( 52) =  "P n n a   : -x+1/2,-y,z; -x+1/2,y+1/2,-z+1/2; -x,-y,-z "
       spg_gen( 53) =  "P m n a   : -x+1/2,-y,z+1/2; -x+1/2,y,-z+1/2; -x,-y,-z "
       spg_gen( 54) =  "P c c a   : -x+1/2,-y,z; -x,y,-z+1/2; -x,-y,-z "
       spg_gen( 55) =  "P b a m   : -x,-y,z; -x+1/2,y+1/2,-z; -x,-y,-z "
       spg_gen( 56) =  "P c c n   : -x+1/2,-y+1/2,z; -x,y+1/2,-z+1/2; -x,-y,-z "
       spg_gen( 57) =  "P b c m   : -x,-y,z+1/2; -x,y+1/2,-z+1/2; -x,-y,-z "
       spg_gen( 58) =  "P n n m   : -x,-y,z; -x+1/2,y+1/2,-z+1/2; -x,-y,-z "
       spg_gen( 59) =  "P m m n   : -x+1/2,-y+1/2,z; -x,y+1/2,-z; -x,-y,-z "
       spg_gen( 60) =  "P b c n   : -x+1/2,-y+1/2,z+1/2; -x,y,-z+1/2; -x,-y,-z "
       spg_gen( 61) =  "P b c a   : -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; -x,-y,-z "
       spg_gen( 62) =  "P n m a   : -x+1/2,-y,z+1/2; -x,y+1/2,-z; -x,-y,-z "
       spg_gen( 63) =  "C m c m   : x+1/2,y+1/2,z; -x,-y,z+1/2; -x,y,-z+1/2; -x,-y,-z "
       spg_gen( 64) =  "C m c a   : x+1/2,y+1/2,z; -x,-y+1/2,z+1/2; -x,y+1/2,-z+1/2; -x,-y,-z "
       spg_gen( 65) =  "C m m m   : x+1/2,y+1/2,z; -x,-y,z; -x,y,-z; -x,-y,-z "
       spg_gen( 66) =  "C c c m   : x+1/2,y+1/2,z; -x,-y,z; -x,y,-z+1/2; -x,-y,-z "
       spg_gen( 67) =  "C m m a   : x+1/2,y+1/2,z; -x,-y+1/2,z; -x,y+1/2,-z; -x,-y,-z "
       spg_gen( 68) =  "C c c a   : x+1/2,y+1/2,z; -x+1/2,-y,z; -x,y,-z+1/2; -x,-y,-z "
       spg_gen( 69) =  "F m m m   : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; -x,y,-z; -x,-y,-z "
       spg_gen( 70) =  "F d d d   : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x+3/4,-y+3/4,z; -x+3/4,y,-z+3/4; -x,-y,-z "
       spg_gen( 71) =  "I m m m   : x+1/2,y+1/2,z+1/2; -x,-y,z; -x,y,-z; -x,-y,-z "
       spg_gen( 72) =  "I b a m   : x+1/2,y+1/2,z+1/2; -x,-y,z; -x+1/2,y+1/2,-z; -x,-y,-z "
       spg_gen( 73) =  "I b c a   : x+1/2,y+1/2,z+1/2; -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; -x,-y,-z "
       spg_gen( 74) =  "I m m a   : x+1/2,y+1/2,z+1/2; -x,-y+1/2,z; -x,y+1/2,-z; -x,-y,-z "
       spg_gen( 75) =  "P 4       : -x,-y,z; -y,x,z "
       spg_gen( 76) =  "P 41      : -x,-y,z+1/2; -y,x,z+1/4 "
       spg_gen( 77) =  "P 42      : -x,-y,z; -y,x,z+1/2 "
       spg_gen( 78) =  "P 43      : -x,-y,z+1/2; -y,x,z+3/4 "
       spg_gen( 79) =  "I 4       : x+1/2,y+1/2,z+1/2; -x,-y,z; -y,x,z "
       spg_gen( 80) =  "I 41      : x+1/2,y+1/2,z+1/2; -x+1/2,-y+1/2,z+1/2; -y,x+1/2,z+1/4 "
       spg_gen( 81) =  "P -4      : -x,-y,z; y,-x,-z "
       spg_gen( 82) =  "I -4      : x+1/2,y+1/2,z+1/2; -x,-y,z; y,-x,-z "
       spg_gen( 83) =  "P 4/m     : -x,-y,z; -y,x,z; -x,-y,-z "
       spg_gen( 84) =  "P 42/m    : -x,-y,z; -y,x,z+1/2; -x,-y,-z "
       spg_gen( 85) =  "P 4/n     : -x+1/2,-y+1/2,z; -y+1/2,x,z; -x,-y,-z "
       spg_gen( 86) =  "P 42/n    : -x+1/2,-y+1/2,z; -y,x+1/2,z+1/2; -x,-y,-z "
       spg_gen( 87) =  "I 4/m     : x+1/2,y+1/2,z+1/2; -x,-y,z; -y,x,z; -x,-y,-z "
       spg_gen( 88) =  "I 41/a    : x+1/2,y+1/2,z+1/2; -x+1/2,-y,z+1/2; -y+3/4,x+1/4,z+1/4; -x,-y,-z "
       spg_gen( 89) =  "P 4 2 2   : -x,-y,z; -y,x,z; -x,y,-z "
       spg_gen( 90) =  "P 4 21 2  : -x,-y,z; -y+1/2,x+1/2,z; -x+1/2,y+1/2,-z "
       spg_gen( 91) =  "P 41 2 2  : -x,-y,z+1/2; -y,x,z+1/4; -x,y,-z "
       spg_gen( 92) =  "P 41 21 2 : -x,-y,z+1/2; -y+1/2,x+1/2,z+1/4; -x+1/2,y+1/2,-z+1/4 "
       spg_gen( 93) =  "P 42 2 2  : -x,-y,z; -y,x,z+1/2; -x,y,-z "
       spg_gen( 94) =  "P 42 21 2 : -x,-y,z; -y+1/2,x+1/2,z+1/2; -x+1/2,y+1/2,-z+1/2 "
       spg_gen( 95) =  "P 43 2 2  : -x,-y,z+1/2; -y,x,z+3/4; -x,y,-z "
       spg_gen( 96) =  "P 43 21 2 : -x,-y,z+1/2; -y+1/2,x+1/2,z+3/4; -x+1/2,y+1/2,-z+3/4 "
       spg_gen( 97) =  "I 4 2 2   : x+1/2,y+1/2,z+1/2; -x,-y,z; -y,x,z; -x,y,-z "
       spg_gen( 98) =  "I 41 2 2  : x+1/2,y+1/2,z+1/2; -x+1/2,-y+1/2,z+1/2; -y,x+1/2,z+1/4; -x+1/2,y,-z+3/4 "
       spg_gen( 99) =  "P 4 m m   : -x,-y,z; -y,x,z; x,-y,z "
       spg_gen(100) =  "P 4 b m   : -x,-y,z; -y,x,z; x+1/2,-y+1/2,z "
       spg_gen(101) =  "P 42 c m  : -x,-y,z; -y,x,z+1/2; x,-y,z+1/2 "
       spg_gen(102) =  "P 42 n m  : -x,-y,z; -y+1/2,x+1/2,z+1/2; x+1/2,-y+1/2,z+1/2 "
       spg_gen(103) =  "P 4 c c   : -x,-y,z; -y,x,z; x,-y,z+1/2 "
       spg_gen(104) =  "P 4 n c   : -x,-y,z; -y,x,z; x+1/2,-y+1/2,z+1/2 "
       spg_gen(105) =  "P 42 m c  : -x,-y,z; -y,x,z+1/2; x,-y,z "
       spg_gen(106) =  "P 42 b c  : -x,-y,z; -y,x,z+1/2; x+1/2,-y+1/2,z "
       spg_gen(107) =  "I 4 m m   : x+1/2,y+1/2,z+1/2; -x,-y,z; -y,x,z; x,-y,z "
       spg_gen(108) =  "I 4 c m   : x+1/2,y+1/2,z+1/2; -x,-y,z; -y,x,z; x,-y,z+1/2 "
       spg_gen(109) =  "I 41 m d  : x+1/2,y+1/2,z+1/2; -x+1/2,-y+1/2,z+1/2; -y,x+1/2,z+1/4; x,-y,z "
       spg_gen(110) =  "I 41 c d  : x+1/2,y+1/2,z+1/2; -x+1/2,-y+1/2,z+1/2; -y,x+1/2,z+1/4; x,-y,z+1/2 "
       spg_gen(111) =  "P -4 2 m  : -x,-y,z; y,-x,-z; -x,y,-z "
       spg_gen(112) =  "P -4 2 c  : -x,-y,z; y,-x,-z; -x,y,-z+1/2 "
       spg_gen(113) =  "P -4 21 m : -x,-y,z; y,-x,-z; -x+1/2,y+1/2,-z "
       spg_gen(114) =  "P -4 21 c : -x,-y,z; y,-x,-z; -x+1/2,y+1/2,-z+1/2 "
       spg_gen(115) =  "P -4 m 2  : -x,-y,z; y,-x,-z; x,-y,z "
       spg_gen(116) =  "P -4 c 2  : -x,-y,z; y,-x,-z; x,-y,z+1/2 "
       spg_gen(117) =  "P -4 b 2  : -x,-y,z; y,-x,-z; x+1/2,-y+1/2,z "
       spg_gen(118) =  "P -4 n 2  : -x,-y,z; y,-x,-z; x+1/2,-y+1/2,z+1/2 "
       spg_gen(119) =  "I -4 m 2  : x+1/2,y+1/2,z+1/2; -x,-y,z; y,-x,-z; x,-y,z "
       spg_gen(120) =  "I -4 c 2  : x+1/2,y+1/2,z+1/2; -x,-y,z; y,-x,-z; x,-y,z+1/2 "
       spg_gen(121) =  "I -4 2 m  : x+1/2,y+1/2,z+1/2; -x,-y,z; y,-x,-z; -x,y,-z "
       spg_gen(122) =  "I -4 2 d  : x+1/2,y+1/2,z+1/2; -x,-y,z; y,-x,-z; -x+1/2,y,-z+3/4 "
       spg_gen(123) =  "P 4/m m m : -x,-y,z; -y,x,z; -x,y,-z; -x,-y,-z "
       spg_gen(124) =  "P 4/m c c : -x,-y,z; -y,x,z; -x,y,-z+1/2; -x,-y,-z "
       spg_gen(125) =  "P 4/n b m : -x+1/2,-y+1/2,z; -y+1/2,x,z; -x+1/2,y,-z; -x,-y,-z "
       spg_gen(126) =  "P 4/n n c : -x+1/2,-y+1/2,z; -y+1/2,x,z; -x+1/2,y,-z+1/2; -x,-y,-z "
       spg_gen(127) =  "P 4/m b m : -x,-y,z; -y,x,z; -x+1/2,y+1/2,-z; -x,-y,-z "
       spg_gen(128) =  "P 4/m n c : -x,-y,z; -y,x,z; -x+1/2,y+1/2,-z+1/2; -x,-y,-z "
       spg_gen(129) =  "P 4/n m m : -x+1/2,-y+1/2,z; -y+1/2,x,z; -x,y+1/2,-z; -x,-y,-z "
       spg_gen(130) =  "P 4/n c c : -x+1/2,-y+1/2,z; -y+1/2,x,z; -x,y+1/2,-z+1/2; -x,-y,-z "
       spg_gen(131) =  "P 42/m m c: -x,-y,z; -y,x,z+1/2; -x,y,-z; -x,-y,-z "
       spg_gen(132) =  "P 42/m c m: -x,-y,z; -y,x,z+1/2; -x,y,-z+1/2; -x,-y,-z "
       spg_gen(133) =  "P 42/n b c: -x+1/2,-y+1/2,z; -y+1/2,x,z+1/2; -x+1/2,y,-z; -x,-y,-z "
       spg_gen(134) =  "P 42/n n m: -x+1/2,-y+1/2,z; -y+1/2,x,z+1/2; -x+1/2,y,-z+1/2; -x,-y,-z "
       spg_gen(135) =  "P 42/m b c: -x,-y,z; -y,x,z+1/2; -x+1/2,y+1/2,-z; -x,-y,-z "
       spg_gen(136) =  "P 42/m n m: -x,-y,z; -y+1/2,x+1/2,z+1/2; -x+1/2,y+1/2,-z+1/2; -x,-y,-z "
       spg_gen(137) =  "P 42/n m c: -x+1/2,-y+1/2,z; -y+1/2,x,z+1/2; -x,y+1/2,-z; -x,-y,-z "
       spg_gen(138) =  "P 42/n c m: -x+1/2,-y+1/2,z; -y+1/2,x,z+1/2; -x,y+1/2,-z+1/2; -x,-y,-z "
       spg_gen(139) =  "I 4/m m m : x+1/2,y+1/2,z+1/2; -x,-y,z; -y,x,z; -x,y,-z; -x,-y,-z "
       spg_gen(140) =  "I 4/m c m : x+1/2,y+1/2,z+1/2; -x,-y,z; -y,x,z; -x,y,-z+1/2; -x,-y,-z "
       spg_gen(141) =  "I 41/a m d: x+1/2,y+1/2,z+1/2; -x+1/2,-y,z+1/2; -y+1/4,x+3/4,z+1/4; -x+1/2,y,-z+1/2; -x,-y,-z "
       spg_gen(142) =  "I 41/a c d: x+1/2,y+1/2,z+1/2; -x+1/2,-y,z+1/2; -y+1/4,x+3/4,z+1/4; -x+1/2,y,-z; -x,-y,-z "
       spg_gen(143) =  "P 3       : -y,x-y,z "
       spg_gen(144) =  "P 31      : -y,x-y,z+1/3 "
       spg_gen(145) =  "P 32      : -y,x-y,z+2/3 "
       spg_gen(146) =  "R 3       : x+1/3,y+2/3,z+2/3; -y,x-y,z "
       spg_gen(147) =  "P -3      : -y,x-y,z; -x,-y,-z "
       spg_gen(148) =  "R -3      : x+1/3,y+2/3,z+2/3; -y,x-y,z; -x,-y,-z "
       spg_gen(149) =  "P 3 1 2   : -y,x-y,z; -y,-x,-z "
       spg_gen(150) =  "P 3 2 1   : -y,x-y,z; y,x,-z "
       spg_gen(151) =  "P 31 1 2  : -y,x-y,z+1/3; -y,-x,-z+2/3 "
       spg_gen(152) =  "P 31 2 1  : -y,x-y,z+1/3; y,x,-z "
       spg_gen(153) =  "P 32 1 2  : -y,x-y,z+2/3; -y,-x,-z+1/3 "
       spg_gen(154) =  "P 32 2 1  : -y,x-y,z+2/3; y,x,-z "
       spg_gen(155) =  "R 3 2     : x+1/3,y+2/3,z+2/3; -y,x-y,z; y,x,-z "
       spg_gen(156) =  "P 3 m 1   : -y,x-y,z; -y,-x,z "
       spg_gen(157) =  "P 3 1 m   : -y,x-y,z; y,x,z "
       spg_gen(158) =  "P 3 c 1   : -y,x-y,z; -y,-x,z+1/2 "
       spg_gen(159) =  "P 3 1 c   : -y,x-y,z; y,x,z+1/2 "
       spg_gen(160) =  "R 3 m     : x+1/3,y+2/3,z+2/3; -y,x-y,z; -y,-x,z "
       spg_gen(161) =  "R 3 c     : x+1/3,y+2/3,z+2/3; -y,x-y,z; -y,-x,z+1/2 "
       spg_gen(162) =  "P -3 1 m  : -y,x-y,z; -y,-x,-z; -x,-y,-z "
       spg_gen(163) =  "P -3 1 c  : -y,x-y,z; -y,-x,-z+1/2; -x,-y,-z "
       spg_gen(164) =  "P -3 m 1  : -y,x-y,z; y,x,-z; -x,-y,-z "
       spg_gen(165) =  "P -3 c 1  : -y,x-y,z; y,x,-z+1/2; -x,-y,-z "
       spg_gen(166) =  "R -3 m    : x+1/3,y+2/3,z+2/3; -y,x-y,z; y,x,-z; -x,-y,-z "
       spg_gen(167) =  "R -3 c    : x+1/3,y+2/3,z+2/3; -y,x-y,z; y,x,-z+1/2; -x,-y,-z "
       spg_gen(168) =  "P 6       : -y,x-y,z; -x,-y,z "
       spg_gen(169) =  "P 61      : -y,x-y,z+1/3; -x,-y,z+1/2 "
       spg_gen(170) =  "P 65      : -y,x-y,z+2/3; -x,-y,z+1/2 "
       spg_gen(171) =  "P 62      : -y,x-y,z+2/3; -x,-y,z "
       spg_gen(172) =  "P 64      : -y,x-y,z+1/3; -x,-y,z "
       spg_gen(173) =  "P 63      : -y,x-y,z; -x,-y,z+1/2 "
       spg_gen(174) =  "P -6      : -y,x-y,z; x,y,-z "
       spg_gen(175) =  "P 6/m     : -y,x-y,z; -x,-y,z; -x,-y,-z "
       spg_gen(176) =  "P 63/m    : -y,x-y,z; -x,-y,z+1/2; -x,-y,-z "
       spg_gen(177) =  "P 6 2 2   : -y,x-y,z; -x,-y,z; y,x,-z "
       spg_gen(178) =  "P 61 2 2  : -y,x-y,z+1/3; -x,-y,z+1/2; y,x,-z+1/3 "
       spg_gen(179) =  "P 65 2 2  : -y,x-y,z+2/3; -x,-y,z+1/2; y,x,-z+2/3 "
       spg_gen(180) =  "P 62 2 2  : -y,x-y,z+2/3; -x,-y,z; y,x,-z+2/3 "
       spg_gen(181) =  "P 64 2 2  : -y,x-y,z+1/3; -x,-y,z; y,x,-z+1/3 "
       spg_gen(182) =  "P 63 2 2  : -y,x-y,z; -x,-y,z+1/2; y,x,-z "
       spg_gen(183) =  "P 6 m m   : -y,x-y,z; -x,-y,z; -y,-x,z "
       spg_gen(184) =  "P 6 c c   : -y,x-y,z; -x,-y,z; -y,-x,z+1/2 "
       spg_gen(185) =  "P 63 c m  : -y,x-y,z; -x,-y,z+1/2; -y,-x,z+1/2 "
       spg_gen(186) =  "P 63 m c  : -y,x-y,z; -x,-y,z+1/2; -y,-x,z "
       spg_gen(187) =  "P -6 m 2  : -y,x-y,z; x,y,-z; -y,-x,z "
       spg_gen(188) =  "P -6 c 2  : -y,x-y,z; x,y,-z+1/2; -y,-x,z+1/2 "
       spg_gen(189) =  "P -6 2 m  : -y,x-y,z; x,y,-z; y,x,-z "
       spg_gen(190) =  "P -6 2 c  : -y,x-y,z; x,y,-z+1/2; y,x,-z "
       spg_gen(191) =  "P 6/m m m : -y,x-y,z; -x,-y,z; y,x,-z; -x,-y,-z "
       spg_gen(192) =  "P 6/m c c : -y,x-y,z; -x,-y,z; y,x,-z+1/2; -x,-y,-z "
       spg_gen(193) =  "P 63/m c m: -y,x-y,z; -x,-y,z+1/2; y,x,-z+1/2; -x,-y,-z "
       spg_gen(194) =  "P 63/m m c: -y,x-y,z; -x,-y,z+1/2; y,x,-z; -x,-y,-z "
       spg_gen(195) =  "P 2 3     : -x,-y,z; -x,y,-z; z,x,y "
       spg_gen(196) =  "F 2 3     : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; -x,y,-z; z,x,y "
       spg_gen(197) =  "I 2 3     : x+1/2,y+1/2,z+1/2; -x,-y,z; -x,y,-z; z,x,y "
       spg_gen(198) =  "P 21 3    : -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y "
       spg_gen(199) =  "I 21 3    : x+1/2,y+1/2,z+1/2; -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y "
       spg_gen(200) =  "P m -3    : -x,-y,z; -x,y,-z; z,x,y; -x,-y,-z "
       spg_gen(201) =  "P n -3    : -x+1/2,-y+1/2,z; -x+1/2,y,-z+1/2; z,x,y; -x,-y,-z "
       spg_gen(202) =  "F m -3    : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; -x,y,-z; z,x,y; -x,-y,-z "
       spg_gen(203) =  "F d -3    : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x+1/4,-y+1/4,z; -x+1/4,y,-z+1/4; z,x,y; -x,-y,-z "
       spg_gen(204) =  "I m -3    : x+1/2,y+1/2,z+1/2; -x,-y,z; -x,y,-z; z,x,y; -x,-y,-z "
       spg_gen(205) =  "P a -3    : -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y; -x,-y,-z "
       spg_gen(206) =  "I a -3    : x+1/2,y+1/2,z+1/2; -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y; -x,-y,-z "
       spg_gen(207) =  "P 4 3 2   : -x,-y,z; -x,y,-z; z,x,y; y,x,-z "
       spg_gen(208) =  "P 42 3 2  : -x,-y,z; -x,y,-z; z,x,y; y+1/2,x+1/2,-z+1/2 "
       spg_gen(209) =  "F 4 3 2   : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; -x,y,-z; z,x,y; y,x,-z "
       spg_gen(210) =  "F 41 3 2  : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y+1/2,z+1/2; -x+1/2,y+1/2,-z; z,x,y; y+3/4,x+1/4,-z+3/4 "
       spg_gen(211) =  "I 4 3 2   : x+1/2,y+1/2,z+1/2; -x,-y,z; -x,y,-z; z,x,y; y,x,-z "
       spg_gen(212) =  "P 43 3 2  : -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y; y+1/4,x+3/4,-z+3/4 "
       spg_gen(213) =  "P 41 3 2  : -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y; y+3/4,x+1/4,-z+1/4 "
       spg_gen(214) =  "I 41 3 2  : x+1/2,y+1/2,z+1/2; -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y; y+3/4,x+1/4,-z+1/4 "
       spg_gen(215) =  "P -4 3 m  : -x,-y,z; -x,y,-z; z,x,y; y,x,z "
       spg_gen(216) =  "F -4 3 m  : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; -x,y,-z; z,x,y; y,x,z "
       spg_gen(217) =  "I -4 3 m  : x+1/2,y+1/2,z+1/2; -x,-y,z; -x,y,-z; z,x,y; y,x,z "
       spg_gen(218) =  "P -4 3 n  : -x,-y,z; -x,y,-z; z,x,y; y+1/2,x+1/2,z+1/2 "
       spg_gen(219) =  "F -4 3 c  : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; -x,y,-z; z,x,y; y+1/2,x+1/2,z+1/2 "
       spg_gen(220) =  "I -4 3 d  : x+1/2,y+1/2,z+1/2; -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y; y+1/4,x+1/4,z+1/4 "
       spg_gen(221) =  "P m -3 m  : -x,-y,z; -x,y,-z; z,x,y; y,x,-z; -x,-y,-z "
       spg_gen(222) =  "P n -3 n  : -x+1/2,-y+1/2,z; -x+1/2,y,-z+1/2; z,x,y; y,x,-z+1/2; -x,-y,-z "
       spg_gen(223) =  "P m -3 n  : -x,-y,z; -x,y,-z; z,x,y; y+1/2,x+1/2,-z+1/2; -x,-y,-z "
       spg_gen(224) =  "P n -3 m  : -x+1/2,-y+1/2,z; -x+1/2,y,-z+1/2; z,x,y; y+1/2,x+1/2,-z; -x,-y,-z "
       spg_gen(225) =  "F m -3 m  : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; -x,y,-z; z,x,y; y,x,-z; -x,-y,-z "
       spg_gen(226) =  "F m -3 c  : x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x,-y,z; -x,y,-z; z,x,y; y+1/2,x+1/2,-z+1/2; -x,-y,-z "

       spg_gen(227) =  &
       "F d -3 m  : x+1/2,y+1/2,z; x+1/2,y,z+1/2;-x+3/4,-y+1/4,z+1/2;-x+1/4,y+1/2,-z+3/4;z,x,y;y+3/4,x+1/4,-z+1/2;-x,-y,-z"

       spg_gen(228) =  &
       "F d -3 c  : x+1/2,y+1/2,z; x+1/2,y,z+1/2;-x+1/4,-y+3/4,z+1/2;-x+3/4,y+1/2,-z+1/4;z,x,y;y+3/4,x+1/4,-z;-x,-y,-z"

       spg_gen(229) =  "I m -3 m  : x+1/2,y+1/2,z+1/2; -x,-y,z; -x,y,-z; z,x,y; y,x,-z; -x,-y,-z "
       spg_gen(230) =  "I a -3 d  : x+1/2,y+1/2,z+1/2; -x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y; y+3/4,x+1/4,-z+1/4; -x,-y,-z "

       return
    End Subroutine Set_It_Gen


    !!----
    !!---- Subroutine Set_Spgr_Info()
    !!----    Number of the Space Group
    !!----    Hermann-Mauguin Symbol
    !!----    Hall symbol
    !!----    Laue Group                                                                                                 ----
    !!----    Point Group
    !!----    Asymmetric unit in direct space.
    !!----    Miscellaneous Information depending on crystal system:
    !!----        Monoclinic         b           c           a
    !!----                        abc  c-ba   abc  ba-c   abc -acb
    !!----                        ---------   ---------   --------
    !!----        cell choice 1    b1   -b1    c1   -c1    a1  -a1
    !!----        cell choice 2    b2   -b2    c2   -c2    a2  -a2
    !!----        cell choice 3    b3   -b3    c3   -c3    a3  -a3
    !!----        Orthorhombic     ba-c   change of basis abc -> ba-c
    !!----                         1      origin choice 1
    !!----                         2ba-c  origin choice 2, change basis
    !!----                                abc -> ba-c
    !!----        Tetragonal       1      origin choice 1
    !!----        Cubic            2      origin choice 2
    !!----        Trigonal         H      hexagonal axes
    !!----                         R      rhombohedral axes
    !!----
    !!----    Set Information on Spgr_info array
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Set_Spgr_Info()

       if (.not. allocated(spgr_info) ) allocate(spgr_info(612) )

       !---- Triclinic ----!
       spgr_info(1:14)= (/                                           &
            spgr_info_type(  1,"P 1         ","P 1             ", 1, 1, (/ 0, 0, 0, 24, 24, 24/),"     ") , &
            spgr_info_type(  1,"A 1         ","A 1             ", 1, 1, (/ 0, 0, 0, 24, 24, 24/),"     ") , &
            spgr_info_type(  1,"B 1         ","B 1             ", 1, 1, (/ 0, 0, 0, 24, 24, 24/),"     ") , &
            spgr_info_type(  1,"C 1         ","C 1             ", 1, 1, (/ 0, 0, 0, 24, 24, 24/),"     ") , &
            spgr_info_type(  1,"I 1         ","I 1             ", 1, 1, (/ 0, 0, 0, 24, 24, 24/),"     ") , &
            spgr_info_type(  1,"R 1         ","R 1             ", 1, 1, (/ 0, 0, 0, 24, 24, 24/),"     ") , &
            spgr_info_type(  1,"F 1         ","F 1             ", 1, 1, (/ 0, 0, 0, 24, 24, 24/),"     ") , &
            spgr_info_type(  2,"P -1        ","-P 1            ", 1, 2, (/ 0, 0, 0, 12, 24, 24/),"     ") , &
            spgr_info_type(  2,"A -1        ","-A 1            ", 1, 2, (/ 0, 0, 0, 12, 24, 24/),"     ") , &
            spgr_info_type(  2,"B -1        ","-B 1            ", 1, 2, (/ 0, 0, 0, 12, 24, 24/),"     ") , &
            spgr_info_type(  2,"C -1        ","-C 1            ", 1, 2, (/ 0, 0, 0, 12, 24, 24/),"     ") , &
            spgr_info_type(  2,"I -1        ","-I 1            ", 1, 2, (/ 0, 0, 0, 12, 24, 24/),"     ") , &
            spgr_info_type(  2,"R -1        ","-R 1            ", 1, 2, (/ 0, 0, 0, 12, 24, 24/),"     ") , &
            spgr_info_type(  2,"F -1        ","-F 1            ", 1, 2, (/ 0, 0, 0, 12, 24, 24/),"     ") /)

       !---- Monoclinic ----!
       spgr_info(15:44)= (/                                           &
            spgr_info_type(  3,"P 1 2 1     ","P 2y            ", 2, 3, (/ 0, 0, 0, 24, 24, 12/),"b    ") , &
            spgr_info_type(  3,"P 2         ","P 2y            ", 2, 3, (/ 0, 0, 0, 24, 24, 12/),"b    ") , &
            spgr_info_type(  3,"P 1 1 2     ","P 2             ", 2, 3, (/ 0, 0, 0, 12, 24, 24/),"c    ") , &
            spgr_info_type(  3,"P 2 1 1     ","P 2x            ", 2, 3, (/ 0, 0, 0, 24, 12, 24/),"a    ") , &
            spgr_info_type(  4,"P 1 21 1    ","P 2yb           ", 2, 3, (/ 0, 0, 0, 24, 24, 12/),"b    ") , &
            spgr_info_type(  4,"P 21        ","P 2yb           ", 2, 3, (/ 0, 0, 0, 24, 24, 12/),"b    ") , &
            spgr_info_type(  4,"P 1 1 21    ","P 2c            ", 2, 3, (/ 0, 0, 0, 12, 24, 24/),"c    ") , &
            spgr_info_type(  4,"P 21 1 1    ","P 2xa           ", 2, 3, (/ 0, 0, 0, 24, 12, 24/),"a    ") , &
            spgr_info_type(  5,"C 1 2 1     ","C 2y            ", 2, 3, (/ 0, 0, 0, 12, 12, 24/),"b1   ") , &
            spgr_info_type(  5,"C 2         ","C 2y            ", 2, 3, (/ 0, 0, 0, 12, 12, 24/),"b1   ") , &
            spgr_info_type(  5,"A 1 2 1     ","A 2y            ", 2, 3, (/ 0, 0, 0, 12, 12, 24/),"b2   ") , &
            spgr_info_type(  5,"A 2         ","A 2y            ", 2, 3, (/ 0, 0, 0, 12, 12, 24/),"b2   ") , &
            spgr_info_type(  5,"I 1 2 1     ","I 2y            ", 2, 3, (/ 0, 0, 0, 12, 12, 24/),"b2   ") , &
            spgr_info_type(  5,"I 2         ","I 2y            ", 2, 3, (/ 0, 0, 0, 12, 12, 24/),"b2   ") , &
            spgr_info_type(  5,"A 1 1 2     ","A 2             ", 2, 3, (/ 0, 0, 0, 24, 12, 12/),"c1   ") , &
            spgr_info_type(  5,"B 1 1 2     ","B 2             ", 2, 3, (/ 0, 0, 0, 24, 12, 12/),"c2   ") , &
            spgr_info_type(  5,"I 1 1 2     ","I 2             ", 2, 3, (/ 0, 0, 0, 24, 12, 12/),"c3   ") , &
            spgr_info_type(  5,"B 2 1 1     ","B 2x            ", 2, 3, (/ 0, 0, 0, 12, 24, 12/),"a1   ") , &
            spgr_info_type(  5,"C 2 1 1     ","C 2x            ", 2, 3, (/ 0, 0, 0, 12, 24, 12/),"a2   ") , &
            spgr_info_type(  5,"I 2 1 1     ","I 2x            ", 2, 3, (/ 0, 0, 0, 12, 24, 12/),"a3   ") , &
            spgr_info_type(  6,"P 1 M 1     ","P -2y           ", 2, 4, (/ 0, 0, 0, 24, 12, 24/),"b    ") , &
            spgr_info_type(  6,"P M         ","P -2y           ", 2, 4, (/ 0, 0, 0, 24, 12, 24/),"b    ") , &
            spgr_info_type(  6,"P 1 1 M     ","P -2            ", 2, 4, (/ 0, 0, 0, 24, 24, 12/),"c    ") , &
            spgr_info_type(  6,"P M 1 1     ","P -2x           ", 2, 4, (/ 0, 0, 0, 12, 24, 24/),"a    ") , &
            spgr_info_type(  7,"P 1 C 1     ","P -2yc          ", 2, 4, (/ 0, 0, 0, 24, 12, 24/),"b1   ") , &
            spgr_info_type(  7,"P C         ","P -2yc          ", 2, 4, (/ 0, 0, 0, 24, 12, 24/),"b1   ") , &
            spgr_info_type(  7,"P 1 N 1     ","P -2yac         ", 2, 4, (/ 0, 0, 0, 24, 12, 24/),"b2   ") , &
            spgr_info_type(  7,"P N         ","P -2yac         ", 2, 4, (/ 0, 0, 0, 24, 12, 24/),"b2   ") , &
            spgr_info_type(  7,"P 1 A 1     ","P -2ya          ", 2, 4, (/ 0, 0, 0, 24, 12, 24/),"b3   ") , &
            spgr_info_type(  7,"P A         ","P -2ya          ", 2, 4, (/ 0, 0, 0, 24, 12, 24/),"b3   ") /)

       spgr_info(45:74)= (/                                           &
            spgr_info_type(  7,"P 1 1 A     ","P -2a           ", 2, 4, (/ 0, 0, 0, 24, 24, 12/),"c1   ") , &
            spgr_info_type(  7,"P 1 1 N     ","P -2ab          ", 2, 4, (/ 0, 0, 0, 24, 24, 12/),"c2   ") , &
            spgr_info_type(  7,"P 1 1 B     ","P -2b           ", 2, 4, (/ 0, 0, 0, 24, 24, 12/),"c3   ") , &
            spgr_info_type(  7,"P B 1 1     ","P -2xb          ", 2, 4, (/ 0, 0, 0, 12, 24, 24/),"a1   ") , &
            spgr_info_type(  7,"P N 1 1     ","P -2xbc         ", 2, 4, (/ 0, 0, 0, 12, 24, 24/),"a2   ") , &
            spgr_info_type(  7,"P C 1 1     ","P -2xc          ", 2, 4, (/ 0, 0, 0, 12, 24, 24/),"a3   ") , &
            spgr_info_type(  8,"C 1 M 1     ","C -2y           ", 2, 4, (/ 0, 0, 0, 24,  6, 24/),"b1   ") , &
            spgr_info_type(  8,"C M         ","C -2y           ", 2, 4, (/ 0, 0, 0, 24,  6, 24/),"b1   ") , &
            spgr_info_type(  8,"A 1 M 1     ","A -2y           ", 2, 4, (/ 0, 0, 0, 24,  6, 24/),"b2   ") , &
            spgr_info_type(  8,"A M         ","A -2y           ", 2, 4, (/ 0, 0, 0, 24,  6, 24/),"b2   ") , &
            spgr_info_type(  8,"I 1 M 1     ","I -2y           ", 2, 4, (/ 0, 0, 0, 24,  6, 24/),"b3   ") , &
            spgr_info_type(  8,"I M         ","I -2y           ", 2, 4, (/ 0, 0, 0, 24,  6, 24/),"b3   ") , &
            spgr_info_type(  8,"A 1 1 M     ","A -2            ", 2, 4, (/ 0, 0, 0, 24, 24,  6/),"c1   ") , &
            spgr_info_type(  8,"B 1 1 M     ","B -2            ", 2, 4, (/ 0, 0, 0, 24, 24,  6/),"c2   ") , &
            spgr_info_type(  8,"I 1 1 M     ","I -2            ", 2, 4, (/ 0, 0, 0, 24, 24,  6/),"c3   ") , &
            spgr_info_type(  8,"B M 1 1     ","B -2x           ", 2, 4, (/ 0, 0, 0,  6, 24, 24/),"a1   ") , &
            spgr_info_type(  8,"C M 1 1     ","C -2x           ", 2, 4, (/ 0, 0, 0,  6, 24, 24/),"a2   ") , &
            spgr_info_type(  8,"I M 1 1     ","I -2x           ", 2, 4, (/ 0, 0, 0,  6, 24, 24/),"a3   ") , &
            spgr_info_type(  9,"C 1 C 1     ","C -2yc          ", 2, 4, (/ 0, 0, 0, 24,  6, 24/),"b1   ") , &
            spgr_info_type(  9,"C C         ","C -2yc          ", 2, 4, (/ 0, 0, 0, 24,  6, 24/),"b1   ") , &
            spgr_info_type(  9,"A 1 N 1     ","A -2yac         ", 2, 4, (/ 0, 0, 0, 24,  6, 24/),"b2   ") , &
            spgr_info_type(  9,"A N         ","A -2yac         ", 2, 4, (/ 0, 0, 0, 24,  6, 24/),"b2   ") , &
            spgr_info_type(  9,"I 1 A 1     ","I -2ya          ", 2, 4, (/ 0, 0, 0, 24,  6, 24/),"b3   ") , &
            spgr_info_type(  9,"I A         ","I -2ya          ", 2, 4, (/ 0, 0, 0, 24,  6, 24/),"b3   ") , &
            spgr_info_type(  9,"A 1 A 1     ","A -2ya          ", 2, 4, (/ 0, 0, 0, 24,  6, 24/),"-b1  ") , &
            spgr_info_type(  9,"A A         ","A -2ya          ", 2, 4, (/ 0, 0, 0, 24,  6, 24/),"-b1  ") , &
            spgr_info_type(  9,"C 1 N 1     ","C -2ybc         ", 2, 4, (/ 0, 0, 0, 24,  6, 24/),"-b2  ") , &
            spgr_info_type(  9,"C N         ","C -2ybc         ", 2, 4, (/ 0, 0, 0, 24,  6, 24/),"-b2  ") , &
            spgr_info_type(  9,"I 1 C 1     ","I -2yc          ", 2, 4, (/ 0, 0, 0, 24,  6, 24/),"-b3  ") , &
            spgr_info_type(  9,"I C         ","I -2yc          ", 2, 4, (/ 0, 0, 0, 24,  6, 24/),"-b3  ") /)

       spgr_info(75:104)= (/                                           &
            spgr_info_type(  9,"A 1 1 A     ","A -2a           ", 2, 4, (/ 0, 0, 0, 24, 24,  6/),"c1   ") , &
            spgr_info_type(  9,"B 1 1 N     ","B -2bc          ", 2, 4, (/ 0, 0, 0, 24, 24,  6/),"c2   ") , &
            spgr_info_type(  9,"I 1 1 B     ","I -2b           ", 2, 4, (/ 0, 0, 0, 24, 24,  6/),"c3   ") , &
            spgr_info_type(  9,"B 1 1 B     ","B -2b           ", 2, 4, (/ 0, 0, 0, 24, 24,  6/),"-c1  ") , &
            spgr_info_type(  9,"A 1 1 N     ","A -2ac          ", 2, 4, (/ 0, 0, 0, 24, 24,  6/),"-c2  ") , &
            spgr_info_type(  9,"I 1 1 A     ","I -2a           ", 2, 4, (/ 0, 0, 0, 24, 24,  6/),"-c3  ") , &
            spgr_info_type(  9,"B B 1 1     ","B -2xb          ", 2, 4, (/ 0, 0, 0,  6, 24, 24/),"a1   ") , &
            spgr_info_type(  9,"C N 1 1     ","C -2xbc         ", 2, 4, (/ 0, 0, 0,  6, 24, 24/),"a2   ") , &
            spgr_info_type(  9,"I C 1 1     ","I -2xc          ", 2, 4, (/ 0, 0, 0,  6, 24, 24/),"a3   ") , &
            spgr_info_type(  9,"C C 1 1     ","C -2xc          ", 2, 4, (/ 0, 0, 0,  6, 24, 24/),"-a1  ") , &
            spgr_info_type(  9,"B N 1 1     ","B -2xbc         ", 2, 4, (/ 0, 0, 0,  6, 24, 24/),"-a2  ") , &
            spgr_info_type(  9,"I B 1 1     ","I -2xb          ", 2, 4, (/ 0, 0, 0,  6, 24, 24/),"-a3  ") , &
            spgr_info_type( 10,"P 1 2/M 1   ","-P 2y           ", 2, 5, (/ 0, 0, 0, 12, 12, 24/),"b    ") , &
            spgr_info_type( 10,"P 2/M       ","-P 2y           ", 2, 5, (/ 0, 0, 0, 12, 12, 24/),"b    ") , &
            spgr_info_type( 10,"P 1 1 2/M   ","-P 2            ", 2, 5, (/ 0, 0, 0, 24, 12, 12/),"c    ") , &
            spgr_info_type( 10,"P 2/M 1 1   ","-P 2x           ", 2, 5, (/ 0, 0, 0, 12, 24, 12/),"a    ") , &
            spgr_info_type( 11,"P 1 21/M 1  ","-P 2yb          ", 2, 5, (/ 0, 0, 0, 24,  6, 24/),"b    ") , &
            spgr_info_type( 11,"P 21/M      ","-P 2yb          ", 2, 5, (/ 0, 0, 0, 24,  6, 24/),"b    ") , &
            spgr_info_type( 11,"P 1 1 21/M  ","-P 2c           ", 2, 5, (/ 0, 0, 0, 24, 24,  6/),"c    ") , &
            spgr_info_type( 11,"P 21/M 1 1  ","-P 2xa          ", 2, 5, (/ 0, 0, 0,  6, 24, 24/),"a    ") , &
            spgr_info_type( 11,"B 1 21/M 1  ","-B 2yb          ", 2, 5, (/ 0, 0, 0, 24,  6, 24/),"b    ") , &
            spgr_info_type( 11,"B 21/M      ","-B 2yb          ", 2, 5, (/ 0, 0, 0, 24,  6, 24/),"b    ") , &
            spgr_info_type( 12,"C 1 2/M 1   ","-C 2y           ", 2, 5, (/ 0, 0, 0, 12,  6, 24/),"b1   ") , &
            spgr_info_type( 12,"C 2/M       ","-C 2y           ", 2, 5, (/ 0, 0, 0, 12,  6, 24/),"b1   ") , &
            spgr_info_type( 12,"A 1 2/M 1   ","-A 2y           ", 2, 5, (/ 0, 0, 0, 12,  6, 24/),"b2   ") , &
            spgr_info_type( 12,"A 2/M       ","-A 2y           ", 2, 5, (/ 0, 0, 0, 12,  6, 24/),"b2   ") , &
            spgr_info_type( 12,"I 1 2/M 1   ","-I 2y           ", 2, 5, (/ 0, 0, 0, 12,  6, 24/),"b3   ") , &
            spgr_info_type( 12,"I 2/M       ","-I 2y           ", 2, 5, (/ 0, 0, 0, 12,  6, 24/),"b3   ") , &
            spgr_info_type( 12,"A 1 1 2/M   ","-A 2            ", 2, 5, (/ 0, 0, 0, 24, 12,  6/),"c1   ") , &
            spgr_info_type( 12,"B 1 1 2/M   ","-B 2            ", 2, 5, (/ 0, 0, 0, 24, 12,  6/),"c2   ") /)

       spgr_info(105:134)= (/                                           &
            spgr_info_type( 12,"I 1 1 2/M   ","-I 2            ", 2, 5, (/ 0, 0, 0, 24, 12,  6/),"c3   ") , &
            spgr_info_type( 12,"B 2/M 1 1   ","-B 2x           ", 2, 5, (/ 0, 0, 0,  6, 24, 12/),"a1   ") , &
            spgr_info_type( 12,"C 2/M 1 1   ","-C 2x           ", 2, 5, (/ 0, 0, 0,  6, 24, 12/),"a2   ") , &
            spgr_info_type( 12,"I 2/M 1 1   ","-I 2x           ", 2, 5, (/ 0, 0, 0,  6, 24, 12/),"a3   ") , &
            spgr_info_type( 12,"F 1 2/M 1   ","-F 2y           ", 2, 5, (/ 0, 0, 0, 12,  6, 24/),"b1   ") , &
            spgr_info_type( 12,"F 2/M       ","-F 2y           ", 2, 5, (/ 0, 0, 0, 12,  6, 24/),"b1   ") , &
            spgr_info_type( 13,"P 1 2/C 1   ","-P 2yc          ", 2, 5, (/ 0, 0, 0, 12, 24, 12/),"b1   ") , &
            spgr_info_type( 13,"P 2/C       ","-P 2yc          ", 2, 5, (/ 0, 0, 0, 12, 24, 12/),"b1   ") , &
            spgr_info_type( 13,"P 1 2/C 1   ","-P 2yc          ", 2, 5, (/ 0, 0, 0, 12, 24, 12/),"b1   ") , &
            spgr_info_type( 13,"P 1 2/N 1   ","-P 2yac         ", 2, 5, (/ 0, 0, 0, 24, 24,  6/),"b2   ") , &
            spgr_info_type( 13,"P 2/N       ","-P 2yac         ", 2, 5, (/ 0, 0, 0, 24, 24,  6/),"b2   ") , &
            spgr_info_type( 13,"P 1 2/A 1   ","-P 2ya          ", 2, 5, (/ 0, 0, 0, 12, 24, 12/),"b3   ") , &
            spgr_info_type( 13,"P 2/A       ","-P 2ya          ", 2, 5, (/ 0, 0, 0, 12, 24, 12/),"b3   ") , &
            spgr_info_type( 13,"P 1 1 2/A   ","-P 2a           ", 2, 5, (/ 0, 0, 0, 12, 12, 24/),"c1   ") , &
            spgr_info_type( 13,"C 1 1 2/A   ","-C 2a           ", 2, 5, (/ 0, 0, 0, 12, 12, 24/),"c1   ") , &
            spgr_info_type( 13,"P 1 1 2/N   ","-P 2ab          ", 2, 5, (/ 0, 0, 0,  6, 24, 24/),"c2   ") , &
            spgr_info_type( 13,"P 1 1 2/B   ","-P 2b           ", 2, 5, (/ 0, 0, 0, 12, 12, 24/),"c3   ") , &
            spgr_info_type( 13,"P 2/B 1 1   ","-P 2xb          ", 2, 5, (/ 0, 0, 0, 24, 12, 12/),"a1   ") , &
            spgr_info_type( 13,"P 2/N 1 1   ","-P 2xbc         ", 2, 5, (/ 0, 0, 0, 24,  6, 24/),"a2   ") , &
            spgr_info_type( 13,"P 2/C 1 1   ","-P 2xc          ", 2, 5, (/ 0, 0, 0, 24, 12, 12/),"a3   ") , &
            spgr_info_type( 14,"P 1 21/C 1  ","-P 2ybc         ", 2, 5, (/ 0, 0, 0, 24,  6, 24/),"b1   ") , &
            spgr_info_type( 14,"P 21/C      ","-P 2ybc         ", 2, 5, (/ 0, 0, 0, 24,  6, 24/),"b1   ") , &
            spgr_info_type( 14,"B 1 21/C 1  ","-B 2ybc         ", 2, 5, (/ 0, 0, 0, 24,  6, 24/),"b1   ") , &
            spgr_info_type( 14,"B 21/C      ","-B 2ybc         ", 2, 5, (/ 0, 0, 0, 24,  6, 24/),"b1   ") , &
            spgr_info_type( 14,"P 1 21/N 1  ","-P 2yn          ", 2, 5, (/ 0, 0, 0, 24,  6, 24/),"b2   ") , &
            spgr_info_type( 14,"P 21/N      ","-P 2yn          ", 2, 5, (/ 0, 0, 0, 24,  6, 24/),"b2   ") , &
            spgr_info_type( 14,"P 1 21/A 1  ","-P 2yab         ", 2, 5, (/ 0, 0, 0, 24,  6, 24/),"b3   ") , &
            spgr_info_type( 14,"P 21/A      ","-P 2yab         ", 2, 5, (/ 0, 0, 0, 24,  6, 24/),"b3   ") , &
            spgr_info_type( 14,"P 1 1 21/A  ","-P 2ac          ", 2, 5, (/ 0, 0, 0, 24, 24,  6/),"c1   ") , &
            spgr_info_type( 14,"P 1 1 21/N  ","-P 2n           ", 2, 5, (/ 0, 0, 0, 24, 24,  6/),"c2   ") /)

       spgr_info(135:162)= (/                                           &
            spgr_info_type( 14,"P 1 1 21/B  ","-P 2bc          ", 2, 5, (/ 0, 0, 0, 24, 24,  6/),"c3   ") , &
            spgr_info_type( 14,"P 21/B 1 1  ","-P 2xab         ", 2, 5, (/ 0, 0, 0,  6, 24, 24/),"a1   ") , &
            spgr_info_type( 14,"P 21/N 1 1  ","-P 2xn          ", 2, 5, (/ 0, 0, 0,  6, 24, 24/),"a2   ") , &
            spgr_info_type( 14,"P 21/C 1 1  ","-P 2xac         ", 2, 5, (/ 0, 0, 0,  6, 24, 24/),"a3   ") , &
            spgr_info_type( 15,"C 1 2/C 1   ","-C 2yc          ", 2, 5, (/ 0, 0, 0, 12, 12, 12/),"b1   ") , &
            spgr_info_type( 15,"C 2/C       ","-C 2yc          ", 2, 5, (/ 0, 0, 0, 12, 12, 12/),"b1   ") , &
            spgr_info_type( 15,"A 1 2/N 1   ","-A 2yac         ", 2, 5, (/ 0, 0, 0, 12, 24,  6/),"b2   ") , &
            spgr_info_type( 15,"A 2/N       ","-A 2yac         ", 2, 5, (/ 0, 0, 0, 12, 24,  6/),"b2   ") , &
            spgr_info_type( 15,"I 1 2/A 1   ","-I 2ya          ", 2, 5, (/ 0, 0, 0, 24, 12,  6/),"b3   ") , &
            spgr_info_type( 15,"I 2/A       ","-I 2ya          ", 2, 5, (/ 0, 0, 0, 24, 12,  6/),"b3   ") , &
            spgr_info_type( 15,"A 1 2/A 1   ","-A 2ya          ", 2, 5, (/ 0, 0, 0, 12, 12, 12/),"-b1  ") , &
            spgr_info_type( 15,"A 2/A       ","-A 2ya          ", 2, 5, (/ 0, 0, 0, 12, 12, 12/),"-b1  ") , &
            spgr_info_type( 15,"C 1 2/N 1   ","-C 2ybc         ", 2, 5, (/ 0, 0, 0,  6, 24, 12/),"-b2  ") , &
            spgr_info_type( 15,"C 2/N       ","-C 2ybc         ", 2, 5, (/ 0, 0, 0,  6, 24, 12/),"-b2  ") , &
            spgr_info_type( 15,"I 1 2/C 1   ","-I 2yc          ", 2, 5, (/ 0, 0, 0,  6, 12, 24/),"-b3  ") , &
            spgr_info_type( 15,"I 2/C       ","-I 2yc          ", 2, 5, (/ 0, 0, 0,  6, 12, 24/),"-b3  ") , &
            spgr_info_type( 15,"A 1 1 2/A   ","-A 2a           ", 2, 5, (/ 0, 0, 0, 12, 12, 12/),"c1   ") , &
            spgr_info_type( 15,"B 1 1 2/N   ","-B 2bc          ", 2, 5, (/ 0, 0, 0,  6, 12, 24/),"c2   ") , &
            spgr_info_type( 15,"I 1 1 2/B   ","-I 2b           ", 2, 5, (/ 0, 0, 0,  6, 24, 12/),"c3   ") , &
            spgr_info_type( 15,"B 1 1 2/B   ","-B 2b           ", 2, 5, (/ 0, 0, 0, 12, 12, 12/),"-c1  ") , &
            spgr_info_type( 15,"A 1 1 2/N   ","-A 2ac          ", 2, 5, (/ 0, 0, 0, 12,  6, 24/),"-c2  ") , &
            spgr_info_type( 15,"I 1 1 2/A   ","-I 2a           ", 2, 5, (/ 0, 0, 0, 24,  6, 12/),"-c3  ") , &
            spgr_info_type( 15,"B 2/B 1 1   ","-B 2xb          ", 2, 5, (/ 0, 0, 0, 12, 12, 12/),"a1   ") , &
            spgr_info_type( 15,"C 2/N 1 1   ","-C 2xbc         ", 2, 5, (/ 0, 0, 0, 24,  6, 12/),"a2   ") , &
            spgr_info_type( 15,"I 2/C 1 1   ","-I 2xc          ", 2, 5, (/ 0, 0, 0, 24,  6, 12/),"a3   ") , &
            spgr_info_type( 15,"C 2/C 1 1   ","-C 2xc          ", 2, 5, (/ 0, 0, 0, 12, 12, 12/),"-a1  ") , &
            spgr_info_type( 15,"B 2/N 1 1   ","-B 2xbc         ", 2, 5, (/ 0, 0, 0, 24, 12,  6/),"-a2  ") , &
            spgr_info_type( 15,"I 2/B 1 1   ","-I 2xb          ", 2, 5, (/ 0, 0, 0, 24, 12,  6/),"-a3  ") /)

       !---- Orthorhombic ----!
       spgr_info(163:192)= (/                                           &
            spgr_info_type( 16,"P 2 2 2     ","P 2 2           ", 3, 6, (/ 0, 0, 0, 12, 12, 24/),"     ") , &
            spgr_info_type( 17,"P 2 2 21    ","P 2c 2          ", 3, 6, (/ 0, 0, 0, 12, 12, 24/),"     ") , &
            spgr_info_type( 17,"P 21 2 2    ","P 2a 2a         ", 3, 6, (/ 0, 0, 0, 24, 12, 12/),"cab  ") , &
            spgr_info_type( 17,"P 2 21 2    ","P 2 2b          ", 3, 6, (/ 0, 0, 0, 12, 24, 12/),"bca  ") , &
            spgr_info_type( 18,"P 21 21 2   ","P 2 2ab         ", 3, 6, (/ 0, 0, 0, 12, 12, 24/),"     ") , &
            spgr_info_type( 18,"P 2 21 21   ","P 2bc 2         ", 3, 6, (/ 0, 0, 0, 24, 12, 12/),"cab  ") , &
            spgr_info_type( 18,"P 21 2 21   ","P 2ac 2ac       ", 3, 6, (/ 0, 0, 0, 12, 24, 12/),"bca  ") , &
            spgr_info_type( 19,"P 21 21 21  ","P 2ac 2ab       ", 3, 6, (/ 0, 0, 0, 12, 12, 24/),"     ") , &
            spgr_info_type( 20,"C 2 2 21    ","C 2c 2          ", 3, 6, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type( 20,"A 21 2 2    ","A 2a 2a         ", 3, 6, (/ 0, 0, 0, 12, 12, 12/),"cab  ") , &
            spgr_info_type( 20,"B 2 21 2    ","B 2 2b          ", 3, 6, (/ 0, 0, 0, 12, 12, 12/),"bca  ") , &
            spgr_info_type( 21,"C 2 2 2     ","C 2 2           ", 3, 6, (/ 0, 0, 0,  6, 12, 24/),"     ") , &
            spgr_info_type( 21,"A 2 2 2     ","A 2 2           ", 3, 6, (/ 0, 0, 0, 24,  6, 12/),"cab  ") , &
            spgr_info_type( 21,"B 2 2 2     ","B 2 2           ", 3, 6, (/ 0, 0, 0, 12, 24,  6/),"bca  ") , &
            spgr_info_type( 22,"F 2 2 2     ","F 2 2           ", 3, 6, (/ 0, 0, 0,  6,  6, 24/),"     ") , &
            spgr_info_type( 23,"I 2 2 2     ","I 2 2           ", 3, 6, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type( 24,"I 21 21 21  ","I 2b 2c         ", 3, 6, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type( 25,"P M M 2     ","P 2 -2          ", 3, 7, (/ 0, 0, 0, 12, 12, 24/),"     ") , &
            spgr_info_type( 25,"P 2 M M     ","P -2 2          ", 3, 9, (/ 0, 0, 0, 24, 12, 12/),"cab  ") , &
            spgr_info_type( 25,"P M 2 M     ","P -2 -2         ", 3, 8, (/ 0, 0, 0, 12, 24, 12/),"bca  ") , &
            spgr_info_type( 26,"P M C 21    ","P 2c -2         ", 3, 7, (/ 0, 0, 0, 12, 12, 24/),"     ") , &
            spgr_info_type( 26,"P C M 21    ","P 2c -2c        ", 3, 7, (/ 0, 0, 0, 12, 12, 24/),"ba-c ") , &
            spgr_info_type( 26,"P 21 M A    ","P -2a 2a        ", 3, 9, (/ 0, 0, 0, 24, 12, 12/),"cab  ") , &
            spgr_info_type( 26,"P 21 A M    ","P -2 2a         ", 3, 9, (/ 0, 0, 0, 24, 12, 12/),"-cba ") , &
            spgr_info_type( 26,"P B 21 M    ","P -2 -2b        ", 3, 8, (/ 0, 0, 0, 12, 24, 12/),"bca  ") , &
            spgr_info_type( 26,"P M 21 B    ","P -2b -2        ", 3, 8, (/ 0, 0, 0, 12, 24, 12/),"a-cb ") , &
            spgr_info_type( 27,"P C C 2     ","P 2 -2c         ", 3, 7, (/ 0, 0, 0, 12, 12, 24/),"     ") , &
            spgr_info_type( 27,"P 2 A A     ","P -2a 2         ", 3, 9, (/ 0, 0, 0, 24, 12, 12/),"cab  ") , &
            spgr_info_type( 27,"P B 2 B     ","P -2b -2b       ", 3, 8, (/ 0, 0, 0, 12, 24, 12/),"bca  ") , &
            spgr_info_type( 28,"P M A 2     ","P 2 -2a         ", 3, 7, (/ 0, 0, 0,  6, 24, 24/),"     ") /)

       spgr_info(193:222)= (/                                           &
            spgr_info_type( 28,"P B M 2     ","P 2 -2b         ", 3, 7, (/ 0, 0, 0, 24,  6, 24/),"ba-c ") , &
            spgr_info_type( 28,"P 2 M B     ","P -2b 2         ", 3, 9, (/ 0, 0, 0, 24,  6, 24/),"cab  ") , &
            spgr_info_type( 28,"P 2 C M     ","P -2c 2         ", 3, 9, (/ 0, 0, 0, 24, 24,  6/),"-cba ") , &
            spgr_info_type( 28,"P C 2 M     ","P -2c -2c       ", 3, 8, (/ 0, 0, 0, 24, 24,  6/),"bca  ") , &
            spgr_info_type( 28,"P M 2 A     ","P -2a -2a       ", 3, 8, (/ 0, 0, 0,  6, 24, 24/),"a-cb ") , &
            spgr_info_type( 29,"P C A 21    ","P 2c -2ac       ", 3, 7, (/ 0, 0, 0,  6, 24, 24/),"     ") , &
            spgr_info_type( 29,"P B C 21    ","P 2c -2b        ", 3, 7, (/ 0, 0, 0, 24,  6, 24/),"ba-c ") , &
            spgr_info_type( 29,"P 21 A B    ","P -2b 2a        ", 3, 9, (/ 0, 0, 0, 24,  6, 24/),"cab  ") , &
            spgr_info_type( 29,"P 21 C A    ","P -2ac 2a       ", 3, 9, (/ 0, 0, 0, 24, 24,  6/),"-cba ") , &
            spgr_info_type( 29,"P C 21 B    ","P -2bc -2c      ", 3, 8, (/ 0, 0, 0, 24, 24,  6/),"bca  ") , &
            spgr_info_type( 29,"P B 21 A    ","P -2a -2ab      ", 3, 8, (/ 0, 0, 0,  6, 24, 24/),"a-cb ") , &
            spgr_info_type( 30,"P N C 2     ","P 2 -2bc        ", 3, 7, (/ 0, 0, 0, 12, 24, 12/),"     ") , &
            spgr_info_type( 30,"P C N 2     ","P 2 -2ac        ", 3, 7, (/ 0, 0, 0, 24, 12, 12/),"ba-c ") , &
            spgr_info_type( 30,"P 2 N A     ","P -2ac 2        ", 3, 9, (/ 0, 0, 0, 12, 12, 24/),"cab  ") , &
            spgr_info_type( 30,"P 2 A N     ","P -2ab 2        ", 3, 9, (/ 0, 0, 0, 12, 24, 12/),"-cba ") , &
            spgr_info_type( 30,"P B 2 N     ","P -2ab -2ab     ", 3, 8, (/ 0, 0, 0, 24, 12, 12/),"bca  ") , &
            spgr_info_type( 30,"P N 2 B     ","P -2bc -2bc     ", 3, 8, (/ 0, 0, 0, 12, 12, 24/),"a-cb ") , &
            spgr_info_type( 31,"P M N 21    ","P 2ac -2        ", 3, 7, (/ 0, 0, 0, 12, 12, 24/),"     ") , &
            spgr_info_type( 31,"P N M 21    ","P 2bc -2bc      ", 3, 7, (/ 0, 0, 0, 12, 12, 24/),"ba-c ") , &
            spgr_info_type( 31,"P 21 M N    ","P -2ab 2ab      ", 3, 9, (/ 0, 0, 0, 24, 12, 12/),"cab  ") , &
            spgr_info_type( 31,"P 21 N M    ","P -2 2ac        ", 3, 9, (/ 0, 0, 0, 24, 12, 12/),"-cba ") , &
            spgr_info_type( 31,"P N 21 M    ","P -2 -2bc       ", 3, 8, (/ 0, 0, 0, 12, 24, 12/),"bca  ") , &
            spgr_info_type( 31,"P M 21 N    ","P -2ab -2       ", 3, 8, (/ 0, 0, 0, 12, 24, 12/),"a-cb ") , &
            spgr_info_type( 32,"P B A 2     ","P 2 -2ab        ", 3, 7, (/ 0, 0, 0, 12, 12, 24/),"     ") , &
            spgr_info_type( 32,"P 2 C B     ","P -2bc 2        ", 3, 9, (/ 0, 0, 0, 24, 12, 12/),"cab  ") , &
            spgr_info_type( 32,"P C 2 A     ","P -2ac -2ac     ", 3, 8, (/ 0, 0, 0, 12, 24, 12/),"bca  ") , &
            spgr_info_type( 33,"P N A 21    ","P 2c -2n        ", 3, 7, (/ 0, 0, 0, 12, 12, 24/),"     ") , &
            spgr_info_type( 33,"P B N 21    ","P 2c -2ab       ", 3, 7, (/ 0, 0, 0, 12, 12, 24/),"ba-c ") , &
            spgr_info_type( 33,"P 21 N B    ","P -2bc 2a       ", 3, 9, (/ 0, 0, 0, 24, 12, 12/),"cab  ") , &
            spgr_info_type( 33,"P 21 C N    ","P -2n 2a        ", 3, 9, (/ 0, 0, 0, 24, 12, 12/),"-cba ") /)

       spgr_info(223:252)= (/                                           &
            spgr_info_type( 33,"P C 21 N    ","P -2n -2ac      ", 3, 8, (/ 0, 0, 0, 12, 24, 12/),"bca  ") , &
            spgr_info_type( 33,"P N 21 A    ","P -2ac -2n      ", 3, 8, (/ 0, 0, 0, 12, 24, 12/),"a-cb ") , &
            spgr_info_type( 34,"P N N 2     ","P 2 -2n         ", 3, 7, (/ 0, 0, 0, 12, 12, 24/),"     ") , &
            spgr_info_type( 34,"P 2 N N     ","P -2n 2         ", 3, 9, (/ 0, 0, 0, 24, 12, 12/),"cab  ") , &
            spgr_info_type( 34,"P N 2 N     ","P -2n -2n       ", 3, 8, (/ 0, 0, 0, 12, 24, 12/),"bca  ") , &
            spgr_info_type( 35,"C M M 2     ","C 2 -2          ", 3, 7, (/ 0, 0, 0,  6, 12, 24/),"     ") , &
            spgr_info_type( 35,"A 2 M M     ","A -2 2          ", 3, 9, (/ 0, 0, 0, 24,  6, 12/),"cab  ") , &
            spgr_info_type( 35,"B M 2 M     ","B -2 -2         ", 3, 8, (/ 0, 0, 0, 12, 24,  6/),"bca  ") , &
            spgr_info_type( 36,"C M C 21    ","C 2c -2         ", 3, 7, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type( 36,"C C M 21    ","C 2c -2c        ", 3, 7, (/ 0, 0, 0, 12, 12, 12/),"ba-c ") , &
            spgr_info_type( 36,"A 21 M A    ","A -2a 2a        ", 3, 9, (/ 0, 0, 0, 12, 12, 12/),"cab  ") , &
            spgr_info_type( 36,"A 21 A M    ","A -2 2a         ", 3, 9, (/ 0, 0, 0, 12, 12, 12/),"-cba ") , &
            spgr_info_type( 36,"B B 21 M    ","B -2 -2b        ", 3, 8, (/ 0, 0, 0, 12, 12, 12/),"bca  ") , &
            spgr_info_type( 36,"B M 21 B    ","B -2b -2        ", 3, 8, (/ 0, 0, 0, 12, 12, 12/),"a-cb ") , &
            spgr_info_type( 37,"C C C 2     ","C 2 -2c         ", 3, 7, (/ 0, 0, 0,  6, 12, 24/),"     ") , &
            spgr_info_type( 37,"A 2 A A     ","A -2a 2         ", 3, 9, (/ 0, 0, 0, 24,  6, 12/),"cab  ") , &
            spgr_info_type( 37,"B B 2 B     ","B -2b -2b       ", 3, 8, (/ 0, 0, 0, 12, 24,  6/),"bca  ") , &
            spgr_info_type( 38,"A M M 2     ","A 2 -2          ", 3, 7, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type( 38,"B M M 2     ","B 2 -2          ", 3, 7, (/ 0, 0, 0, 12, 12, 12/),"ba-c ") , &
            spgr_info_type( 38,"B 2 M M     ","B -2 2          ", 3, 9, (/ 0, 0, 0, 12, 12, 12/),"cab  ") , &
            spgr_info_type( 38,"C 2 M M     ","C -2 2          ", 3, 9, (/ 0, 0, 0, 12, 12, 12/),"-cba ") , &
            spgr_info_type( 38,"C M 2 M     ","C -2 -2         ", 3, 8, (/ 0, 0, 0, 12, 12, 12/),"bca  ") , &
            spgr_info_type( 38,"A M 2 M     ","A -2 -2         ", 3, 8, (/ 0, 0, 0, 12, 12, 12/),"a-cb ") , &
            spgr_info_type( 39,"A B M 2     ","A 2 -2c         ", 3, 7, (/ 0, 0, 0, 12,  6, 24/),"     ") , &
            spgr_info_type( 39,"B M A 2     ","B 2 -2c         ", 3, 7, (/ 0, 0, 0,  6, 12, 24/),"ba-c ") , &
            spgr_info_type( 39,"B 2 C M     ","B -2c 2         ", 3, 9, (/ 0, 0, 0, 24, 12,  6/),"cab  ") , &
            spgr_info_type( 39,"C 2 M B     ","C -2b 2         ", 3, 9, (/ 0, 0, 0, 24,  6, 12/),"-cba ") , &
            spgr_info_type( 39,"C M 2 A     ","C -2b -2b       ", 3, 8, (/ 0, 0, 0,  6, 24, 12/),"bca  ") , &
            spgr_info_type( 39,"A C 2 M     ","A -2c -2c       ", 3, 8, (/ 0, 0, 0, 12, 24,  6/),"a-cb ") , &
            spgr_info_type( 40,"A M A 2     ","A 2 -2a         ", 3, 7, (/ 0, 0, 0,  6, 12, 24/),"     ") /)

       spgr_info(253:282)= (/                                           &
            spgr_info_type( 40,"B B M 2     ","B 2 -2b         ", 3, 7, (/ 0, 0, 0, 12,  6, 24/),"ba-c ") , &
            spgr_info_type( 40,"B 2 M B     ","B -2b 2         ", 3, 9, (/ 0, 0, 0, 24,  6, 12/),"cab  ") , &
            spgr_info_type( 40,"C 2 C M     ","C -2c 2         ", 3, 9, (/ 0, 0, 0, 24, 12,  6/),"-cba ") , &
            spgr_info_type( 40,"C C 2 M     ","C -2c -2c       ", 3, 8, (/ 0, 0, 0, 12, 24,  6/),"bca  ") , &
            spgr_info_type( 40,"A M 2 A     ","A -2a -2a       ", 3, 8, (/ 0, 0, 0,  6, 24, 12/),"a-cb ") , &
            spgr_info_type( 41,"A B A 2     ","A 2 -2ac        ", 3, 7, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type( 41,"B B A 2     ","B 2 -2bc        ", 3, 7, (/ 0, 0, 0, 12, 12, 12/),"ba-c ") , &
            spgr_info_type( 41,"B 2 C B     ","B -2bc 2        ", 3, 9, (/ 0, 0, 0, 12, 12, 12/),"cab  ") , &
            spgr_info_type( 41,"C 2 C B     ","C -2bc 2        ", 3, 9, (/ 0, 0, 0, 12, 12, 12/),"-cba ") , &
            spgr_info_type( 41,"C C 2 A     ","C -2bc -2bc     ", 3, 8, (/ 0, 0, 0, 12, 12, 12/),"bca  ") , &
            spgr_info_type( 41,"A C 2 A     ","A -2ac -2ac     ", 3, 8, (/ 0, 0, 0, 12, 12, 12/),"a-cb ") , &
            spgr_info_type( 42,"F M M 2     ","F 2 -2          ", 3, 7, (/ 0, 0, 0,  6,  6, 24/),"     ") , &
            spgr_info_type( 42,"F 2 M M     ","F -2 2          ", 3, 9, (/ 0, 0, 0, 24,  6,  6/),"cab  ") , &
            spgr_info_type( 42,"F M 2 M     ","F -2 -2         ", 3, 8, (/ 0, 0, 0,  6, 24,  6/),"bca  ") , &
            spgr_info_type( 43,"F D D 2     ","F 2 -2d         ", 3, 7, (/ 0, 0, 0,  6,  6, 24/),"     ") , &
            spgr_info_type( 43,"F 2 D D     ","F -2d 2         ", 3, 9, (/ 0, 0, 0, 24,  6,  6/),"cab  ") , &
            spgr_info_type( 43,"F D 2 D     ","F -2d -2d       ", 3, 8, (/ 0, 0, 0,  6, 24,  6/),"bca  ") , &
            spgr_info_type( 44,"I M M 2     ","I 2 -2          ", 3, 7, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type( 44,"I 2 M M     ","I -2 2          ", 3, 9, (/ 0, 0, 0, 12, 12, 12/),"cab  ") , &
            spgr_info_type( 44,"I M 2 M     ","I -2 -2         ", 3, 8, (/ 0, 0, 0, 12, 12, 12/),"bca  ") , &
            spgr_info_type( 45,"I B A 2     ","I 2 -2c         ", 3, 7, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type( 45,"I 2 C B     ","I -2a 2         ", 3, 9, (/ 0, 0, 0, 12, 12, 12/),"cab  ") , &
            spgr_info_type( 45,"I C 2 A     ","I -2b -2b       ", 3, 8, (/ 0, 0, 0, 12, 12, 12/),"bca  ") , &
            spgr_info_type( 46,"I M A 2     ","I 2 -2a         ", 3, 7, (/ 0, 0, 0,  6, 24, 12/),"     ") , &
            spgr_info_type( 46,"I B M 2     ","I 2 -2b         ", 3, 7, (/ 0, 0, 0, 24,  6, 12/),"ba-c ") , &
            spgr_info_type( 46,"I 2 M B     ","I -2b 2         ", 3, 9, (/ 0, 0, 0, 12,  6, 24/),"cab  ") , &
            spgr_info_type( 46,"I 2 C M     ","I -2c 2         ", 3, 9, (/ 0, 0, 0, 12, 24,  6/),"-cba ") , &
            spgr_info_type( 46,"I C 2 M     ","I -2c -2c       ", 3, 8, (/ 0, 0, 0, 24, 12,  6/),"bca  ") , &
            spgr_info_type( 46,"I M 2 A     ","I -2a -2a       ", 3, 8, (/ 0, 0, 0,  6, 12, 12/),"a-cb ") , &
            spgr_info_type( 47,"P M M M     ","-P 2 2          ", 3,10, (/ 0, 0, 0, 12, 12, 12/),"     ") /)

       spgr_info(283:312)= (/                                           &
            spgr_info_type( 48,"P N N N:1   ","P 2 2 -1n       ", 3,10, (/ 0, 0, 0,  6, 12, 24/),"1    ") , &
            spgr_info_type( 48,"P N N N     ","-P 2ab 2bc      ", 3,10, (/ 0,-6, 0,  6,  6, 24/),"2    ") , &
            spgr_info_type( 49,"P C C M     ","-P 2 2c         ", 3,10, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type( 49,"P M A A     ","-P 2a 2         ", 3,10, (/ 0, 0, 0, 12, 12, 12/),"cab  ") , &
            spgr_info_type( 49,"P B M B     ","-P 2b 2b        ", 3,10, (/ 0, 0, 0, 12, 12, 12/),"bca  ") , &
            spgr_info_type( 50,"P B A N:1   ","P 2 2 -1ab      ", 3,10, (/ 0, 0, 0, 12, 12, 12/),"1    ") , &
            spgr_info_type( 50,"P B A N     ","-P 2ab 2b       ", 3,10, (/ 0, 0, 0,  6, 24, 12/),"2    ") , &
            spgr_info_type( 50,"P N C B:1   ","P 2 2 -1bc      ", 3,10, (/ 0, 0, 0, 12, 12, 12/),"1cab ") , &
            spgr_info_type( 50,"P N C B     ","-P 2b 2bc       ", 3,10, (/ 0, 0, 0, 12,  6, 24/),"2cab ") , &
            spgr_info_type( 50,"P C N A:1   ","P 2 2 -1ac      ", 3,10, (/ 0, 0, 0, 12, 12, 12/),"1bca ") , &
            spgr_info_type( 50,"P C N A     ","-P 2a 2c        ", 3,10, (/ 0, 0, 0, 24, 12,  6/),"2bca ") , &
            spgr_info_type( 51,"P M M A     ","-P 2a 2a        ", 3,10, (/ 0, 0, 0,  6, 12, 24/),"     ") , &
            spgr_info_type( 51,"P M M B     ","-P 2b 2         ", 3,10, (/ 0, 0, 0, 12,  6, 24/),"ba-c ") , &
            spgr_info_type( 51,"P B M M     ","-P 2 2b         ", 3,10, (/ 0, 0, 0, 24,  6, 12/),"cab  ") , &
            spgr_info_type( 51,"P C M M     ","-P 2c 2c        ", 3,10, (/ 0, 0, 0, 24, 12,  6/),"-cba ") , &
            spgr_info_type( 51,"P M C M     ","-P 2c 2         ", 3,10, (/ 0, 0, 0, 12, 24,  6/),"bca  ") , &
            spgr_info_type( 51,"P M A M     ","-P 2 2a         ", 3,10, (/ 0, 0, 0,  6, 24, 12/),"a-cb ") , &
            spgr_info_type( 52,"P N N A     ","-P 2a 2bc       ", 3,10, (/ 0, 0, 0, 24,  6, 12/),"     ") , &
            spgr_info_type( 52,"P N N B     ","-P 2b 2n        ", 3,10, (/ 0, 0, 0,  6, 24, 12/),"ba-c ") , &
            spgr_info_type( 52,"P B N N     ","-P 2n 2b        ", 3,10, (/ 0, 0, 0, 12, 24,  6/),"cab  ") , &
            spgr_info_type( 52,"P C N N     ","-P 2ab 2c       ", 3,10, (/ 0, 0, 0, 12,  6, 24/),"-cba ") , &
            spgr_info_type( 52,"P N C N     ","-P 2ab 2n       ", 3,10, (/ 0, 0, 0,  6, 12, 24/),"bca  ") , &
            spgr_info_type( 52,"P N A N     ","-P 2n 2bc       ", 3,10, (/ 0, 0, 0, 24, 12,  6/),"a-cb ") , &
            spgr_info_type( 53,"P M N A     ","-P 2ac 2        ", 3,10, (/ 0, 0, 0, 12, 24,  6/),"     ") , &
            spgr_info_type( 53,"P N M B     ","-P 2bc 2bc      ", 3,10, (/ 0, 0, 0, 24, 12,  6/),"ba-c ") , &
            spgr_info_type( 53,"P B M N     ","-P 2ab 2ab      ", 3,10, (/ 0, 0, 0,  6, 12, 24/),"cab  ") , &
            spgr_info_type( 53,"P C N M     ","-P 2 2ac        ", 3,10, (/ 0, 0, 0,  6, 24, 12/),"-cba ") , &
            spgr_info_type( 53,"P N C M     ","-P 2 2bc        ", 3,10, (/ 0, 0, 0, 24,  6, 12/),"bca  ") , &
            spgr_info_type( 53,"P M A N     ","-P 2ab 2        ", 3,10, (/ 0, 0, 0, 12,  6, 24/),"a-cb ") , &
            spgr_info_type( 54,"P C C A     ","-P 2a 2ac       ", 3,10, (/ 0, 0, 0, 12, 12, 12/),"     ") /)

       spgr_info(313:342)= (/                                           &
            spgr_info_type( 54,"P C C B     ","-P 2b 2c        ", 3,10, (/ 0, 0, 0, 12, 12, 12/),"ba-c ") , &
            spgr_info_type( 54,"P B A A     ","-P 2a 2b        ", 3,10, (/ 0, 0, 0, 12, 12, 12/),"cab  ") , &
            spgr_info_type( 54,"P C A A     ","-P 2ac 2c       ", 3,10, (/ 0, 0, 0, 12, 12, 12/),"-cba ") , &
            spgr_info_type( 54,"P B C B     ","-P 2bc 2b       ", 3,10, (/ 0, 0, 0, 12, 12, 12/),"bca  ") , &
            spgr_info_type( 54,"P B A B     ","-P 2b 2ab       ", 3,10, (/ 0, 0, 0, 12, 12, 12/),"a-cb ") , &
            spgr_info_type( 55,"P B A M     ","-P 2 2ab        ", 3,10, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type( 55,"P M C B     ","-P 2bc 2        ", 3,10, (/ 0, 0, 0, 12, 12, 12/),"cab  ") , &
            spgr_info_type( 55,"P C M A     ","-P 2ac 2ac      ", 3,10, (/ 0, 0, 0, 12, 12, 12/),"bca  ") , &
            spgr_info_type( 56,"P C C N     ","-P 2ab 2ac      ", 3,10, (/ 0, 0, 0,  6, 24, 12/),"     ") , &
            spgr_info_type( 56,"P N A A     ","-P 2ac 2bc      ", 3,10, (/ 0, 0, 0, 12,  6, 24/),"cab  ") , &
            spgr_info_type( 56,"P B N B     ","-P 2bc 2ab      ", 3,10, (/ 0, 0, 0, 24, 12,  6/),"bca  ") , &
            spgr_info_type( 57,"P B C M     ","-P 2c 2b        ", 3,10, (/ 0, 0, 0, 12, 24,  6/),"     ") , &
            spgr_info_type( 57,"P C A M     ","-P 2c 2ac       ", 3,10, (/ 0, 0, 0, 24, 12,  6/),"ba-c ") , &
            spgr_info_type( 57,"P M C A     ","-P 2ac 2a       ", 3,10, (/ 0, 0, 0,  6, 12, 24/),"cab  ") , &
            spgr_info_type( 57,"P M A B     ","-P 2b 2a        ", 3,10, (/ 0, 0, 0,  6, 24, 12/),"-cba ") , &
            spgr_info_type( 57,"P B M A     ","-P 2a 2ab       ", 3,10, (/ 0, 0, 0, 24,  6, 12/),"bca  ") , &
            spgr_info_type( 57,"P C M B     ","-P 2bc 2c       ", 3,10, (/ 0, 0, 0, 12,  6, 24/),"a-cb ") , &
            spgr_info_type( 58,"P N N M     ","-P 2 2n         ", 3,10, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type( 58,"P M N N     ","-P 2n 2         ", 3,10, (/ 0, 0, 0, 12, 12, 12/),"cab  ") , &
            spgr_info_type( 58,"P N M N     ","-P 2n 2n        ", 3,10, (/ 0, 0, 0, 12, 12, 12/),"bca  ") , &
            spgr_info_type( 59,"P M M N:1   ","P 2 2ab -1ab    ", 3,10, (/ 0, 0, 0, 12, 12, 12/),"1    ") , &
            spgr_info_type( 59,"P M M N     ","-P 2ab 2a       ", 3,10, (/ 0,-6, 0,  6,  6, 24/),"2    ") , &
            spgr_info_type( 59,"P N M M:1   ","P 2bc 2 -1bc    ", 3,10, (/ 0, 0, 0, 12, 12, 12/),"1cab ") , &
            spgr_info_type( 59,"P N M M     ","-P 2c 2bc       ", 3,10, (/ 0, 0,-6, 24,  6,  6/),"2cab ") , &
            spgr_info_type( 59,"P M N M:1   ","P 2ac 2ac -1ac  ", 3,10, (/ 0, 0, 0, 12, 12, 12/),"1bca ") , &
            spgr_info_type( 59,"P M N M     ","-P 2c 2a        ", 3,10, (/-6, 0, 0,  6, 24,  6/),"2bca ") , &
            spgr_info_type( 60,"P B C N     ","-P 2n 2ab       ", 3,10, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type( 60,"P C A N     ","-P 2n 2c        ", 3,10, (/ 0, 0, 0, 12, 12, 12/),"ba-c ") , &
            spgr_info_type( 60,"P N C A     ","-P 2a 2n        ", 3,10, (/ 0, 0, 0, 12, 12, 12/),"cab  ") , &
            spgr_info_type( 60,"P N A B     ","-P 2bc 2n       ", 3,10, (/ 0, 0, 0, 12, 12, 12/),"-cba ") /)

       spgr_info(343:372)= (/                                           &
            spgr_info_type( 60,"P B N A     ","-P 2ac 2b       ", 3,10, (/ 0, 0, 0, 12, 12, 12/),"bca  ") , &
            spgr_info_type( 60,"P C N B     ","-P 2b 2ac       ", 3,10, (/ 0, 0, 0, 12, 12, 12/),"a-cb ") , &
            spgr_info_type( 61,"P B C A     ","-P 2ac 2ab      ", 3,10, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type( 61,"P C A B     ","-P 2bc 2ac      ", 3,10, (/ 0, 0, 0, 12, 12, 12/),"ba-c ") , &
            spgr_info_type( 62,"P N M A     ","-P 2ac 2n       ", 3,10, (/ 0, 0, 0, 12,  6, 24/),"     ") , &
            spgr_info_type( 62,"P M N B     ","-P 2bc 2a       ", 3,10, (/ 0, 0, 0,  6, 12, 24/),"ba-c ") , &
            spgr_info_type( 62,"P M N B:1   ","P 2ac 2ab -1ab  ", 3,10, (/ 0, 0, 0,  6, 12, 24/),"     ") , &
            spgr_info_type( 62,"P B N M     ","-P 2c 2ab       ", 3,10, (/ 0, 0, 0, 24, 12,  6/),"cab  ") , &
            spgr_info_type( 62,"P B N M:1   ","P 2c 2n -1c     ", 3,10, (/ 0, 0, 0, 24, 12,  6/),"     ") , &
            spgr_info_type( 62,"P C M N     ","-P 2n 2ac       ", 3,10, (/ 0, 0, 0, 24,  6, 12/),"-cba ") , &
            spgr_info_type( 62,"P M C N     ","-P 2n 2a        ", 3,10, (/ 0, 0, 0,  6, 24, 12/),"bca  ") , &
            spgr_info_type( 62,"P M C N:1   ","P 2bc 2a -1a    ", 3,10, (/ 0, 0, 0,  6, 24, 12/),"     ") , &
            spgr_info_type( 62,"P N A M     ","-P 2c 2n        ", 3,10, (/ 0, 0, 0, 12, 24,  6/),"a-cb ") , &
            spgr_info_type( 63,"C M C M     ","-C 2c 2         ", 3,10, (/ 0, 0, 0, 12, 12, 24/),"     ") , &
            spgr_info_type( 63,"C C M M     ","-C 2c 2c        ", 3,10, (/ 0, 0, 0, 12, 12, 24/),"ba-c ") , &
            spgr_info_type( 63,"A M M A     ","-A 2a 2a        ", 3,10, (/ 0, 0, 0, 24, 12, 12/),"cab  ") , &
            spgr_info_type( 63,"A M A M     ","-A 2 2a         ", 3,10, (/ 0, 0, 0, 24, 12, 12/),"-cba ") , &
            spgr_info_type( 63,"B B M M     ","-B 2 2b         ", 3,10, (/ 0, 0, 0, 12, 24, 12/),"bca  ") , &
            spgr_info_type( 63,"B M M B     ","-B 2b 2         ", 3,10, (/ 0, 0, 0, 12, 24, 12/),"a-cb ") , &
            spgr_info_type( 63,"B M M B:1   ","B 2ab 2c -1ac   ", 3,10, (/ 0, 0, 0, 12, 24, 12/),"     ") , &
            spgr_info_type( 64,"C M C A     ","-C 2bc 2        ", 3,10, (/ 0, 0, 0,  6, 12, 12/),"     ") , &
            spgr_info_type( 64,"C C M B     ","-C 2bc 2bc      ", 3,10, (/ 0, 0, 0, 12,  6, 12/),"ba-c ") , &
            spgr_info_type( 64,"C C M B:1   ","C 2bc 2n -1ab   ", 3,10, (/ 0, 0, 0, 12,  6, 12/),"     ") , &
            spgr_info_type( 64,"A B M A     ","-A 2ac 2ac      ", 3,10, (/ 0, 0, 0, 12,  6, 12/),"cab  ") , &
            spgr_info_type( 64,"A C A M     ","-A 2 2ac        ", 3,10, (/ 0, 0, 0, 12, 12,  6/),"-cba ") , &
            spgr_info_type( 64,"B B C M     ","-B 2 2bc        ", 3,10, (/ 0, 0, 0, 12, 12,  6/),"bca  ") , &
            spgr_info_type( 64,"B M A B     ","-B 2bc 2        ", 3,10, (/ 0, 0, 0,  6, 12, 12/),"a-cb ") , &
            spgr_info_type( 65,"C M M M     ","-C 2 2          ", 3,10, (/ 0, 0, 0,  6, 12, 12/),"     ") , &
            spgr_info_type( 65,"A M M M     ","-A 2 2          ", 3,10, (/ 0, 0, 0, 12,  6, 12/),"cab  ") , &
            spgr_info_type( 65,"B M M M     ","-B 2 2          ", 3,10, (/ 0, 0, 0, 12, 12,  6/),"bca  ") /)

       spgr_info(373:402)= (/                                           &
            spgr_info_type( 66,"C C C M     ","-C 2 2c         ", 3,10, (/ 0, 0, 0,  6, 12, 12/),"     ") , &
            spgr_info_type( 66,"A M A A     ","-A 2a 2         ", 3,10, (/ 0, 0, 0, 12,  6, 12/),"cab  ") , &
            spgr_info_type( 66,"B A M B     ","-B 2b 2b        ", 3,10, (/ 0, 0, 0, 12, 12,  6/),"bca  ") , &
            spgr_info_type( 67,"C M M A     ","-C 2b 2         ", 3,10, (/ 0, 0, 0, 12,  6, 12/),"     ") , &
            spgr_info_type( 67,"C M M B     ","-C 2b 2b        ", 3,10, (/ 0, 0, 0,  6, 12, 12/),"ba-c ") , &
            spgr_info_type( 67,"A B M M     ","-A 2c 2c        ", 3,10, (/ 0, 0, 0, 12, 12,  6/),"cab  ") , &
            spgr_info_type( 67,"A C M M     ","-A 2 2c         ", 3,10, (/ 0, 0, 0, 12,  6, 12/),"-cba ") , &
            spgr_info_type( 67,"B M C M     ","-B 2 2c         ", 3,10, (/ 0, 0, 0,  6, 12, 12/),"bca  ") , &
            spgr_info_type( 67,"B M A M     ","-B 2c 2         ", 3,10, (/ 0, 0, 0, 12, 12,  6/),"a-cb ") , &
            spgr_info_type( 68,"C C C A:1   ","C 2 2 -1bc      ", 3,10, (/ 0, 0, 0,  6, 12, 12/),"1    ") , &
            spgr_info_type( 68,"C C C A     ","-C 2b 2bc       ", 3,10, (/ 0, 0, 0, 12,  6, 12/),"2    ") , &
            spgr_info_type( 68,"C C C B:1   ","C 2 2 -1bc      ", 3,10, (/ 0, 0, 0, 12,  6, 12/),"1ba-c") , &
            spgr_info_type( 68,"C C C B     ","-C 2b 2c        ", 3,10, (/ 0, 0, 0,  6, 12, 12/),"2ba-c") , &
            spgr_info_type( 68,"A B A A:1   ","A 2 2 -1ac      ", 3,10, (/ 0, 0, 0, 12, 12,  6/),"1cab ") , &
            spgr_info_type( 68,"A B A A     ","-A 2a 2c        ", 3,10, (/ 0, 0, 0, 12, 12,  6/),"2cab ") , &
            spgr_info_type( 68,"A C A A:1   ","A 2 2 -1ac      ", 3,10, (/ 0, 0, 0, 12,  6, 12/),"1-cba") , &
            spgr_info_type( 68,"A C A A     ","-A 2ac 2c       ", 3,10, (/ 0, 0, 0, 12,  6, 12/),"2-cba") , &
            spgr_info_type( 68,"B B C B:1   ","B 2 2 -1bc      ", 3,10, (/ 0, 0, 0, 12, 12,  6/),"1bca ") , &
            spgr_info_type( 68,"B B C B     ","-B 2bc 2b       ", 3,10, (/ 0, 0, 0,  6, 12, 12/),"2bca ") , &
            spgr_info_type( 68,"B B A B:1   ","B 2 2 -1bc      ", 3,10, (/ 0, 0, 0,  6, 12, 12/),"1a-cb") , &
            spgr_info_type( 68,"B B A B     ","-B 2b 2bc       ", 3,10, (/ 0, 0, 0, 12, 12,  6/),"2a-cb") , &
            spgr_info_type( 69,"F M M M     ","-F 2 2          ", 3,10, (/ 0, 0, 0,  6,  6, 12/),"     ") , &
            spgr_info_type( 70,"F D D D:1   ","F 2 2 -1d       ", 3,10, (/ 0, 0, 0,  3,  6, 24/),"1    ") , &
            spgr_info_type( 70,"F D D D     ","-F 2uv 2vw      ", 3,10, (/ 0,-3, 0,  3,  3, 24/),"2    ") , &
            spgr_info_type( 71,"I M M M     ","-I 2 2          ", 3,10, (/ 0, 0, 0,  6, 12, 12/),"     ") , &
            spgr_info_type( 72,"I B A M     ","-I 2 2c         ", 3,10, (/ 0, 0, 0,  6, 12, 12/),"     ") , &
            spgr_info_type( 72,"I M C B     ","-I 2a 2         ", 3,10, (/ 0, 0, 0, 12,  6, 12/),"cab  ") , &
            spgr_info_type( 72,"I C M A:1   ","I 2 2 -1b       ", 3,10, (/ 0, 0, 0, 12, 12,  6/),"     ") , &
            spgr_info_type( 72,"I C M A     ","-I 2b 2b        ", 3,10, (/ 0, 0, 0, 12, 12,  6/),"bca  ") , &
            spgr_info_type( 73,"I B C A     ","-I 2b 2c        ", 3,10, (/ 0, 0, 0,  6, 12, 12/),"     ") /)

       spgr_info(403:409)= (/                                           &
            spgr_info_type( 73,"I C A B     ","-I 2a 2b        ", 3,10, (/ 0, 0, 0, 12,  6, 12/),"ba-c ") , &
            spgr_info_type( 74,"I M M A     ","-I 2b 2         ", 3,10, (/ 0, 0, 0,  6,  6, 24/),"     ") , &
            spgr_info_type( 74,"I M M B     ","-I 2a 2a        ", 3,10, (/ 0, 0, 0,  6,  6, 24/),"ba-c ") , &
            spgr_info_type( 74,"I B M M     ","-I 2c 2c        ", 3,10, (/ 0, 0, 0, 24,  6,  6/),"cab  ") , &
            spgr_info_type( 74,"I C M M     ","-I 2 2b         ", 3,10, (/ 0, 0, 0, 24,  6,  6/),"-cba ") , &
            spgr_info_type( 74,"I M C M     ","-I 2 2a         ", 3,10, (/ 0, 0, 0,  6, 24,  6/),"bca  ") , &
            spgr_info_type( 74,"I M A M     ","-I 2c 2         ", 3,10, (/ 0, 0, 0,  6, 24,  6/),"a-cb ") /)

       !---- Tetragonal ----!
       spgr_info(410:439)= (/                                           &
            spgr_info_type( 75,"P 4         ","P 4             ", 4,11, (/ 0, 0, 0, 12, 12, 24/),"     ") , &
            spgr_info_type( 76,"P 41        ","P 4w            ", 4,11, (/ 0, 0, 0, 12, 12, 24/),"     ") , &
            spgr_info_type( 77,"P 42        ","P 4c            ", 4,11, (/ 0, 0, 0, 12, 12, 24/),"     ") , &
            spgr_info_type( 78,"P 43        ","P 4cw           ", 4,11, (/ 0, 0, 0, 12, 12, 24/),"     ") , &
            spgr_info_type( 79,"I 4         ","I 4             ", 4,11, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type( 80,"I 41        ","I 4bw           ", 4,11, (/ 0, 0, 0, 12, 24,  6/),"     ") , &
            spgr_info_type( 81,"P -4        ","P -4            ", 4,12, (/ 0, 0, 0, 12, 12, 24/),"     ") , &
            spgr_info_type( 82,"I -4        ","I -4            ", 4,12, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type( 83,"P 4/M       ","-P 4            ", 4,13, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type( 84,"P 42/M      ","-P 4c           ", 4,13, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type( 85,"P 4/N:1     ","P 4ab -1ab      ", 4,13, (/ 0, 0, 0, 12, 12, 12/),"1    ") , &
            spgr_info_type( 85,"P 4/N       ","-P 4a           ", 4,13, (/-6,-6, 0,  6,  6, 12/),"2    ") , &
            spgr_info_type( 86,"P 42/N:1    ","P 4n -1n        ", 4,13, (/ 0, 0, 0, 12, 24,  6/),"1    ") , &
            spgr_info_type( 86,"P 42/N      ","-P 4bc          ", 4,13, (/-6,-6, 0,  6,  6, 12/),"2    ") , &
            spgr_info_type( 87,"I 4/M       ","-I 4            ", 4,13, (/ 0, 0, 0, 12, 12,  6/),"     ") , &
            spgr_info_type( 88,"I 41/A:1    ","I 4bw -1bw      ", 4,13, (/ 0, 0, 0,  6,  6, 24/),"1    ") , &
            spgr_info_type( 88,"I 41/A      ","-I 4ad          ", 4,13, (/ 0, 0, 0,  6,  6, 24/),"2    ") , &
            spgr_info_type( 89,"P 4 2 2     ","P 4 2           ", 5,14, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type( 90,"P 4 21 2    ","P 4ab 2ab       ", 5,14, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type( 90,"C 4 2 21    ","C 4b 2          ", 5,14, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type( 91,"P 41 2 2    ","P 4w 2c         ", 5,14, (/ 0, 0, 0, 24, 24,  3/),"     ") , &
            spgr_info_type( 92,"P 41 21 2   ","P 4abw 2nw      ", 5,14, (/ 0, 0, 0, 24, 24,  3/),"     ") , &
            spgr_info_type( 93,"P 42 2 2    ","P 4c 2          ", 5,14, (/ 0, 0, 0, 12, 24,  6/),"     ") , &
            spgr_info_type( 94,"P 42 21 2   ","P 4n 2n         ", 5,14, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type( 95,"P 43 2 2    ","P 4cw 2c        ", 5,14, (/ 0, 0, 0, 24, 24,  3/),"     ") , &
            spgr_info_type( 96,"P 43 21 2   ","P 4nw 2abw      ", 5,14, (/ 0, 0, 0, 24, 24,  3/),"     ") , &
            spgr_info_type( 97,"I 4 2 2     ","I 4 2           ", 5,14, (/ 0, 0, 0, 12, 12,  6/),"     ") , &
            spgr_info_type( 98,"I 41 2 2    ","I 4bw 2bw       ", 5,14, (/ 0, 0, 0, 12, 24,  3/),"     ") , &
            spgr_info_type( 99,"P 4 M M     ","P 4 -2          ", 5,15, (/ 0, 0, 0, 12, 12, 24/),"     ") , &
            spgr_info_type(100,"P 4 B M     ","P 4 -2ab        ", 5,15, (/ 0, 0, 0, 12, 12, 24/),"     ") /)

       spgr_info(440:469)= (/                                           &
            spgr_info_type(101,"P 42 C M    ","P 4c -2c        ", 5,15, (/ 0, 0, 0, 12, 12, 24/),"     ") , &
            spgr_info_type(102,"P 42 N M    ","P 4n -2n        ", 5,15, (/ 0, 0, 0, 12, 12, 24/),"     ") , &
            spgr_info_type(103,"P 4 C C     ","P 4 -2c         ", 5,15, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type(104,"P 4 N C     ","P 4 -2n         ", 5,15, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type(105,"P 42 M C    ","P 4c -2         ", 5,15, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type(106,"P 42 B C    ","P 4c -2ab       ", 5,15, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type(107,"I 4 M M     ","I 4 -2          ", 5,15, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type(108,"I 4 C M     ","I 4 -2c         ", 5,15, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type(109,"I 41 M D    ","I 4bw -2        ", 5,15, (/ 0, 0, 0, 12, 12,  6/),"     ") , &
            spgr_info_type(110,"I 41 C D    ","I 4bw -2c       ", 5,15, (/ 0, 0, 0, 12, 12,  6/),"     ") , &
            spgr_info_type(111,"P -4 2 M    ","P -4 2          ", 5,16, (/ 0, 0, 0, 12, 12, 24/),"     ") , &
            spgr_info_type(112,"P -4 2 C    ","P -4 2c         ", 5,16, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type(113,"P -4 21 M   ","P -4 2ab        ", 5,16, (/ 0, 0, 0, 12, 12, 24/),"     ") , &
            spgr_info_type(114,"P -4 21 C   ","P -4 2n         ", 5,16, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type(115,"P -4 M 2    ","P -4 -2         ", 5,17, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type(116,"P -4 C 2    ","P -4 -2c        ", 5,17, (/ 0, 0, 0, 12, 24,  6/),"     ") , &
            spgr_info_type(117,"P -4 B 2    ","P -4 -2ab       ", 5,17, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type(117,"C -4 B 2    ","C -4 2b         ", 5,17, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type(118,"P -4 N 2    ","P -4 -2n        ", 5,17, (/ 0, 0, 0, 12, 12,  6/),"     ") , &
            spgr_info_type(119,"I -4 M 2    ","I -4 -2         ", 5,17, (/ 0, 0, 0, 12, 12,  6/),"     ") , &
            spgr_info_type(120,"I -4 C 2    ","I -4 -2c        ", 5,17, (/ 0, 0, 0, 12, 12,  6/),"     ") , &
            spgr_info_type(121,"I -4 2 M    ","I -4 2          ", 5,16, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type(122,"I -4 2 D    ","I -4 2bw        ", 5,16, (/ 0, 0, 0, 12, 24,  3/),"     ") , &
            spgr_info_type(122,"F -4 D 2    ","F -4 -2cd       ", 5,16, (/ 0, 0, 0, 12, 24,  3/),"     ") , &
            spgr_info_type(123,"P 4/M M M   ","-P 4 2          ", 5,18, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type(124,"P 4/M C C   ","-P 4 2c         ", 5,18, (/ 0, 0, 0, 12, 12,  6/),"     ") , &
            spgr_info_type(125,"P 4/N B M:1 ","P 4 2 -1ab      ", 5,18, (/ 0, 0, 0, 12, 12, 12/),"1    ") , &
            spgr_info_type(125,"P 4/N B M   ","-P 4a 2b        ", 5,18, (/-6,-6, 0,  6,  6, 12/),"2    ") , &
            spgr_info_type(126,"P 4/N N C:1 ","P 4 2 -1n       ", 5,18, (/ 0, 0, 0, 12, 12,  6/),"1    ") , &
            spgr_info_type(126,"P 4/N N C   ","-P 4a 2bc       ", 5,18, (/-6,-6, 0,  6,  6,  6/),"2    ") /)

       spgr_info(470:494)= (/                                           &
            spgr_info_type(127,"P 4/M B M   ","-P 4 2ab        ", 5,18, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type(128,"P 4/M N C   ","-P 4 2n         ", 5,18, (/ 0, 0, 0, 12, 12,  6/),"     ") , &
            spgr_info_type(129,"P 4/N M M:1 ","P 4ab 2ab -1ab  ", 5,18, (/ 0, 0, 0, 12, 12, 12/),"1    ") , &
            spgr_info_type(129,"P 4/N M M   ","-P 4a 2a        ", 5,18, (/-6,-6, 0,  6,  6, 12/),"2    ") , &
            spgr_info_type(130,"P 4/N C C:1 ","P 4ab 2n -1ab   ", 5,18, (/ 0, 0, 0, 12, 12,  6/),"1    ") , &
            spgr_info_type(130,"P 4/N C C   ","-P 4a 2ac       ", 5,18, (/-6,-6, 0,  6,  6,  6/),"2    ") , &
            spgr_info_type(131,"P 42/M M C  ","-P 4c 2         ", 5,18, (/ 0, 0, 0, 12, 12,  6/),"     ") , &
            spgr_info_type(132,"P 42/M C M  ","-P 4c 2c        ", 5,18, (/ 0, 0, 0, 12, 12,  6/),"     ") , &
            spgr_info_type(133,"P 42/N B C:1","P 4n 2c -1n     ", 5,18, (/ 0, 0, 0, 12, 12,  6/),"1    ") , &
            spgr_info_type(133,"P 42/N B C  ","-P 4ac 2b       ", 5,18, (/-6,-6, 0,  6,  6,  6/),"2    ") , &
            spgr_info_type(134,"P 42/N N M:1","P 4n 2 -1n      ", 5,18, (/ 0, 0, 0, 12, 24,  6/),"1    ") , &
            spgr_info_type(134,"P 42/N N M  ","-P 4ac 2bc      ", 5,18, (/-6,-6, 0,  6,  6, 12/),"2    ") , &
            spgr_info_type(135,"P 42/M B C  ","-P 4c 2ab       ", 5,18, (/ 0, 0, 0, 12, 12,  6/),"     ") , &
            spgr_info_type(136,"P 42/M N M  ","-P 4n 2n        ", 5,18, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type(137,"P 42/N M C:1","P 4n 2n -1n     ", 5,18, (/ 0, 0, 0, 12, 12,  6/),"1    ") , &
            spgr_info_type(137,"P 42/N M C  ","-P 4ac 2a       ", 5,18, (/-6,-6, 0,  6,  6,  6/),"2    ") , &
            spgr_info_type(138,"P 42/N C M:1","P 4n 2ab -1n    ", 5,18, (/ 0, 0, 0,  6, 12, 24/),"1    ") , &
            spgr_info_type(138,"P 42/N C M  ","-P 4ac 2ac      ", 5,18, (/-6,-6, 0,  6,  6, 12/),"2    ") , &
            spgr_info_type(139,"I 4/M M M   ","-I 4 2          ", 5,18, (/ 0, 0, 0, 12, 12,  6/),"     ") , &
            spgr_info_type(139,"F 4/M M M   ","-F 4 2          ", 5,18, (/ 0, 0, 0, 12, 12,  6/),"     ") , &
            spgr_info_type(140,"I 4/M C M   ","-I 4 2c         ", 5,18, (/ 0, 0, 0, 12, 12,  6/),"     ") , &
            spgr_info_type(141,"I 41/A M D:1","I 4bw 2bw -1bw  ", 5,18, (/ 0, 0, 0, 12, 12,  3/),"1    ") , &
            spgr_info_type(141,"I 41/A M D  ","-I 4bd 2        ", 5,18, (/ 0,-6, 0, 12,  6,  3/),"2    ") , &
            spgr_info_type(142,"I 41/A C D:1","I 4bw 2aw -1bw  ", 5,18, (/ 0, 0, 0, 12, 12,  3/),"1    ") , &
            spgr_info_type(142,"I 41/A C D  ","-I 4bd 2c       ", 5,18, (/ 0,-6, 0, 12,  6,  3/),"2    ") /)

       !---- Trigonal/Rhombohedral ----!
       spgr_info(495:526)= (/                                           &
            spgr_info_type(143,"P 3         ","P 3             ", 8,19, (/ 0, 0, 0, 16, 16, 24/),"     ") , &
            spgr_info_type(144,"P 31        ","P 31            ", 8,19, (/ 0, 0, 0, 24, 24,  8/),"     ") , &
            spgr_info_type(145,"P 32        ","P 32            ", 8,19, (/ 0, 0, 0, 24, 24,  8/),"     ") , &
            spgr_info_type(146,"R 3         ","R 3             ", 8,19, (/ 0, 0, 0, 16, 16,  8/),"H    ") , &
            spgr_info_type(146,"R 3:R       ","P 3*            ", 6,19, (/ 0, 0, 0, 24, 24, 24/),"R    ") , &
            spgr_info_type(147,"P -3        ","-P 3            ", 8,20, (/ 0, 0, 0, 16, 16, 12/),"     ") , &
            spgr_info_type(148,"R -3        ","-R 3            ", 8,20, (/ 0, 0, 0, 16, 16,  4/),"H    ") , &
            spgr_info_type(148,"R -3:R      ","-P 3*           ", 6,20, (/ 0, 0, 0, 24, 24, 12/),"R    ") , &
            spgr_info_type(149,"P 3 1 2     ","P 3 2           ",10,24, (/ 0, 0, 0, 16, 16, 12/),"     ") , &
            spgr_info_type(150,"P 3 2 1     ","P 3 2""         ", 9,21, (/ 0, 0, 0, 16, 16, 12/),"     ") , &
            spgr_info_type(151,"P 31 1 2    ","P 31 2c (0 0 1) ",10,24, (/ 0, 0, 0, 24, 24,  4/),"     ") , &
            spgr_info_type(152,"P 31 2 1    ","P 31 2""        ", 9,21, (/ 0, 0, 0, 24, 24,  4/),"     ") , &
            spgr_info_type(153,"P 32 1 2    ","P 32 2c (0 0 -1)",10,24, (/ 0, 0, 0, 24, 24,  4/),"     ") , &
            spgr_info_type(154,"P 32 2 1    ","P 32 2""        ", 9,21, (/ 0, 0, 0, 24, 24,  4/),"     ") , &
            spgr_info_type(155,"R 3 2       ","R 3 2""         ", 9,21, (/ 0, 0, 0, 16, 16,  4/),"H    ") , &
            spgr_info_type(155,"R 3 2:R     ","P 3* 2          ", 7,21, (/ 0, 0, 0, 24, 24, 12/),"R    ") , &
            spgr_info_type(156,"P 3 M 1     ","P 3 -2""        ", 9,22, (/ 0, 0, 0, 16, 16, 24/),"     ") , &
            spgr_info_type(157,"P 3 1 M     ","P 3 -2          ",10,25, (/ 0, 0, 0, 16, 12, 24/),"     ") , &
            spgr_info_type(158,"P 3 C 1     ","P 3 -2""c       ", 9,22, (/ 0, 0, 0, 16, 16, 12/),"     ") , &
            spgr_info_type(159,"P 3 1 C     ","P 3 -2c         ",10,25, (/ 0, 0, 0, 16, 16, 12/),"     ") , &
            spgr_info_type(160,"R 3 M       ","R 3 -2""        ", 9,22, (/ 0, 0, 0, 16, 16,  8/),"H    ") , &
            spgr_info_type(160,"R 3 M:R     ","P 3* -2         ", 7,22, (/ 0, 0, 0, 24, 24, 24/),"R    ") , &
            spgr_info_type(161,"R 3 C       ","R 3 -2""c       ", 9,22, (/ 0, 0, 0, 16, 16,  4/),"H    ") , &
            spgr_info_type(161,"R 3 C:R     ","P 3* -2n        ", 7,22, (/ 0, 0, 0, 24, 24, 24/),"R    ") , &
            spgr_info_type(162,"P -3 1 M    ","-P 3 2          ",10,26, (/ 0, 0, 0, 16, 12, 12/),"     ") , &
            spgr_info_type(163,"P -3 1 C    ","-P 3 2c         ",10,26, (/ 0, 0, 0, 16, 16,  6/),"     ") , &
            spgr_info_type(164,"P -3 M 1    ","-P 3 2""        ", 9,23, (/ 0, 0, 0, 16,  8, 24/),"     ") , &
            spgr_info_type(165,"P -3 C 1    ","-P 3 2""c       ", 9,23, (/ 0, 0, 0, 16, 16,  6/),"     ") , &
            spgr_info_type(166,"R -3 M      ","-R 3 2""        ", 9,23, (/ 0, 0, 0, 16, 16,  4/),"H    ") , &
            spgr_info_type(166,"R -3 M:R    ","-P 3* 2         ", 7,23, (/ 0, 0, 0, 24, 24, 12/),"R    ") , &
            spgr_info_type(167,"R -3 C      ","-R 3 2""c       ", 9,23, (/ 0, 0, 0, 16, 16,  2/),"H    ") , &
            spgr_info_type(167,"R -3 C:R    ","-P 3* 2n        ", 7,23, (/ 6, 6, 6, 30, 30, 18/),"R    ") /)

       !---- Hexagonal ----!
       spgr_info(527:553)= (/                                           &
            spgr_info_type(168,"P 6         ","P 6             ",11,27, (/ 0, 0, 0, 16, 12, 24/),"     ") , &
            spgr_info_type(169,"P 61        ","P 61            ",11,27, (/ 0, 0, 0, 24, 24,  4/),"     ") , &
            spgr_info_type(170,"P 65        ","P 65            ",11,27, (/ 0, 0, 0, 24, 24,  4/),"     ") , &
            spgr_info_type(171,"P 62        ","P 62            ",11,27, (/ 0, 0, 0, 24, 24,  8/),"     ") , &
            spgr_info_type(172,"P 64        ","P 64            ",11,27, (/ 0, 0, 0, 24, 24,  8/),"     ") , &
            spgr_info_type(173,"P 63        ","P 6c            ",11,27, (/ 0, 0, 0, 16, 16, 12/),"     ") , &
            spgr_info_type(174,"P -6        ","P -6            ",11,28, (/ 0, 0, 0, 16, 16, 12/),"     ") , &
            spgr_info_type(175,"P 6/M       ","-P 6            ",11,29, (/ 0, 0, 0, 16, 12, 12/),"     ") , &
            spgr_info_type(176,"P 63/M      ","-P 6c           ",11,29, (/ 0, 0, 0, 16, 16,  6/),"     ") , &
            spgr_info_type(177,"P 6 2 2     ","P 6 2           ",12,30, (/ 0, 0, 0, 16, 12, 12/),"     ") , &
            spgr_info_type(178,"P 61 2 2    ","P 61 2 (0 0 -1) ",12,30, (/ 0, 0, 0, 24, 24,  2/),"     ") , &
            spgr_info_type(179,"P 65 2 2    ","P 65 2 (0 0 1)  ",12,30, (/ 0, 0, 0, 24, 24,  2/),"     ") , &
            spgr_info_type(180,"P 62 2 2    ","P 62 2c (0 0 1) ",12,30, (/ 0, 0, 0, 24, 24,  4/),"     ") , &
            spgr_info_type(181,"P 64 2 2    ","P 64 2c (0 0 -1)",12,30, (/ 0, 0, 0, 24, 24,  4/),"     ") , &
            spgr_info_type(182,"P 63 2 2    ","P 6c 2c         ",12,30, (/ 0, 0, 0, 16, 16,  3/),"     ") , &
            spgr_info_type(183,"P 6 M M     ","P 6 -2          ",12,31, (/ 0, 0, 0, 16,  8, 24/),"     ") , &
            spgr_info_type(184,"P 6 C C     ","P 6 -2c         ",12,31, (/ 0, 0, 0, 16, 12, 12/),"     ") , &
            spgr_info_type(185,"P 63 C M    ","P 6c -2         ",12,31, (/ 0, 0, 0, 16, 12, 12/),"     ") , &
            spgr_info_type(186,"P 63 M C    ","P 6c -2c        ",12,31, (/ 0, 0, 0, 16,  8, 24/),"     ") , &
            spgr_info_type(187,"P -6 M 2    ","P -6 2          ",12,33, (/ 0, 0, 0, 16, 16, 12/),"     ") , &
            spgr_info_type(188,"P -6 C 2    ","P -6c 2         ",12,33, (/ 0, 0, 0, 16, 16,  6/),"     ") , &
            spgr_info_type(189,"P -6 2 M    ","P -6 -2         ",12,32, (/ 0, 0, 0, 16, 12, 12/),"     ") , &
            spgr_info_type(190,"P -6 2 C    ","P -6c -2c       ",12,32, (/ 0, 0, 0, 16, 16,  6/),"     ") , &
            spgr_info_type(191,"P 6/M M M   ","-P 6 2          ",12,34, (/ 0, 0, 0, 16,  8, 12/),"     ") , &
            spgr_info_type(192,"P 6/M C C   ","-P 6 2c         ",12,34, (/ 0, 0, 0, 16, 12,  6/),"     ") , &
            spgr_info_type(193,"P 63/M C M  ","-P 6c 2         ",12,34, (/ 0, 0, 0, 16, 12,  6/),"     ") , &
            spgr_info_type(194,"P 63/M M C  ","-P 6c 2c        ",12,34, (/ 0, 0, 0, 16, 16,  6/),"     ") /)

       !---- Cubic ----!
       spgr_info(554:583)= (/                                           &
            spgr_info_type(195,"P 2 3       ","P 2 2 3         ",13,35, (/ 0, 0,  0, 24, 24, 12/),"     ") , &
            spgr_info_type(196,"F 2 3       ","F 2 2 3         ",13,35, (/ 0, 0, -6, 12, 12,  6/),"     ") , &
            spgr_info_type(197,"I 2 3       ","I 2 2 3         ",13,35, (/ 0, 0,  0, 24, 12, 12/),"     ") , &
            spgr_info_type(198,"P 21 3      ","P 2ac 2ab 3     ",13,35, (/ 0, 0,-12, 12, 12, 12/),"     ") , &
            spgr_info_type(199,"I 21 3      ","I 2b 2c 3       ",13,35, (/ 0, 0,  0, 12, 12, 12/),"     ") , &
            spgr_info_type(200,"P M -3      ","-P 2 2 3        ",13,36, (/ 0, 0,  0, 12, 12, 12/),"     ") , &
            spgr_info_type(200,"P M 3       ","-P 2 2 3        ",13,36, (/ 0, 0,  0, 12, 12, 12/),"     ") , &
            spgr_info_type(201,"P N -3:1    ","P 2 2 3 -1n     ",13,36, (/ 0, 0, 0, 24, 12,  12/),"1    ") , &
            spgr_info_type(201,"P N -3      ","-P 2ab 2bc 3    ",13,36, (/-6,-6,-6, 18,  6,   6/),"2    ") , &
            spgr_info_type(201,"P N 3       ","-P 2ab 2bc 3    ",13,36, (/-6,-6,-6, 18,  6,   6/),"2    ") , &
            spgr_info_type(202,"F M -3      ","-F 2 2 3        ",13,36, (/ 0, 0, 0, 12, 12,   6/),"     ") , &
            spgr_info_type(202,"F M 3       ","-F 2 2 3        ",13,36, (/ 0, 0, 0, 12, 12,   6/),"     ") , &
            spgr_info_type(203,"F D -3:1    ","F 2 2 3 -1d     ",13,36, (/ 0, 0,-6, 12,  6,   6/),"1    ") , &
            spgr_info_type(203,"F D -3      ","-F 2uv 2vw 3    ",13,36, (/-3,-3,-9,  9,  3,   3/),"2    ") , &
            spgr_info_type(203,"F D 3       ","-F 2uv 2vw 3    ",13,36, (/-3,-3,-9,  9,  3,   3/),"2    ") , &
            spgr_info_type(204,"I M -3      ","-I 2 2 3        ",13,36, (/ 0, 0, 0, 12, 12,  12/),"     ") , &
            spgr_info_type(204,"I M 3       ","-I 2 2 3        ",13,36, (/ 0, 0, 0, 12, 12,  12/),"     ") , &
            spgr_info_type(205,"P A -3      ","-P 2ac 2ab 3    ",13,36, (/ 0, 0, 0, 12, 12,  12/),"     ") , &
            spgr_info_type(205,"P A 3       ","-P 2ac 2ab 3    ",13,36, (/ 0, 0, 0, 12, 12,  12/),"     ") , &
            spgr_info_type(206,"I A -3      ","-I 2b 2c 3      ",13,36, (/ 0, 0, 0, 12, 12,   6/),"     ") , &
            spgr_info_type(206,"I A 3       ","-I 2b 2c 3      ",13,36, (/ 0, 0, 0, 12, 12,   6/),"     ") , &
            spgr_info_type(207,"P 4 3 2     ","P 4 2 3         ",14,37, (/ 0, 0, 0, 24, 12,  12/),"     ") , &
            spgr_info_type(208,"P 42 3 2    ","P 4n 2 3        ",14,37, (/ 0, 0,-6, 12, 12,   6/),"     ") , &
            spgr_info_type(209,"F 4 3 2     ","F 4 2 3         ",14,37, (/ 0, 0,-6, 12,  6,   6/),"     ") , &
            spgr_info_type(210,"F 41 3 2    ","F 4d 2 3        ",14,37, (/ 0,-3,-3, 12,  3,   3/),"     ") , &
            spgr_info_type(211,"I 4 3 2     ","I 4 2 3         ",14,37, (/ 0, 0, 0, 12, 12,   6/),"     ") , &
            spgr_info_type(212,"P 43 3 2    ","P 4acd 2ab 3    ",14,37, (/ 0, 0,-12, 12, 18,  6/),"     ") , &
            spgr_info_type(213,"P 41 3 2    ","P 4bd 2ab 3     ",14,37, (/-6, 0, 0, 12, 18,  12/),"     ") , &
            spgr_info_type(214,"I 41 3 2    ","I 4bd 2c 3      ",14,37, (/-9,-3,-3,  3,  3,   9/),"     ") , &
            spgr_info_type(215,"P -4 3 M    ","P -4 2 3        ",14,38, (/ 0, 0, 0, 24, 12,  12/),"     ") /)

       spgr_info(584:612)= (/                                           &
            spgr_info_type(216,"F -4 3 M    ","F -4 2 3        ",14,38, (/ 0, 0,-6, 12,  6,  6/),"     ") , &
            spgr_info_type(217,"I -4 3 M    ","I -4 2 3        ",14,38, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type(218,"P -4 3 N    ","P -4n 2 3       ",14,38, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type(219,"F -4 3 C    ","F -4c 2 3       ",14,38, (/ 0, 0,-6, 12,  6,  6/),"     ") , &
            spgr_info_type(220,"I -4 3 D    ","I -4bd 2c 3     ",14,38, (/ 6, 6, 0, 12, 12, 12/),"     ") , &
            spgr_info_type(221,"P M -3 M    ","-P 4 2 3        ",14,39, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type(221,"P M 3 M     ","-P 4 2 3        ",14,39, (/ 0, 0, 0, 12, 12, 12/),"     ") , &
            spgr_info_type(222,"P N -3 N:1  ","P 4 2 3 -1n     ",14,39, (/ 0, 0, 0, 12, 12, 12/),"1    ") , &
            spgr_info_type(222,"P N -3 N    ","-P 4a 2bc 3     ",14,39, (/ 6, 6, 6, 18, 18, 18/),"2    ") , &
            spgr_info_type(222,"P N 3 N     ","-P 4a 2bc 3     ",14,39, (/ 6, 6, 6, 18, 18, 18/),"2    ") , &
            spgr_info_type(223,"P M -3 N    ","-P 4n 2 3       ",14,39, (/ 0, 0, 0, 12, 12,  6/),"     ") , &
            spgr_info_type(223,"P M 3 N     ","-P 4n 2 3       ",14,39, (/ 0, 0, 0, 12, 12,  6/),"     ") , &
            spgr_info_type(224,"P N -3 M:1  ","P 4n 2 3 -1n    ",14,39, (/ 0, 0,-6, 12, 12,  6/),"1    ") , &
            spgr_info_type(224,"P N -3 M    ","-P 4bc 2bc 3    ",14,39, (/ 6, 6, 0, 18, 18, 12/),"2    ") , &
            spgr_info_type(224,"P N 3 M     ","-P 4bc 2bc 3    ",14,39, (/ 6, 6, 0, 18, 18, 12/),"2    ") , &
            spgr_info_type(225,"F M -3 M    ","-F 4 2 3        ",14,39, (/ 0, 0, 0, 12,  6,  6/),"     ") , &
            spgr_info_type(225,"F M 3 M     ","-F 4 2 3        ",14,39, (/ 0, 0, 0, 12,  6,  6/),"     ") , &
            spgr_info_type(226,"F M -3 C    ","-F 4c 2 3       ",14,39, (/ 0, 0, 0, 12,  6,  6/),"     ") , &
            spgr_info_type(226,"F M 3 C     ","-F 4c 2 3       ",14,39, (/ 0, 0, 0, 12,  6,  6/),"     ") , &
            spgr_info_type(227,"F D -3 M:1  ","F 4d 2 3 -1d    ",14,39, (/ 0, 0,-3, 12,  3,  3/),"1    ") , &
            spgr_info_type(227,"F D -3 M    ","-F 4vw 2vw 3    ",14,39, (/-3,-3,-6,  9,  0,  0/),"2    ") , &
            spgr_info_type(227,"F D 3 M     ","-F 4vw 2vw 3    ",14,39, (/-3,-3,-6,  9,  0,  0/),"2    ") , &
            spgr_info_type(228,"F D -3 C:1  ","F 4d 2 3 -1cd   ",14,39, (/ 0, 0,-3, 12,  3,  3/),"1    ") , &
            spgr_info_type(228,"F D -3 C    ","-F 4cvw 2vw 3   ",14,39, (/-3,-3,-6,  9,  0,  0/),"2    ") , &
            spgr_info_type(228,"F D 3 C     ","-F 4cvw 2vw 3   ",14,39, (/-3,-3,-6,  9,  0,  0/),"2    ") , &
            spgr_info_type(229,"I M -3 M    ","-I 4 2 3        ",14,39, (/ 0, 0, 0, 12, 12,  6/),"     ") , &
            spgr_info_type(229,"I M 3 M     ","-I 4 2 3        ",14,39, (/ 0, 0, 0, 12, 12,  6/),"     ") , &
            spgr_info_type(230,"I A -3 D    ","-I 4bd 2c 3     ",14,39, (/-3,-3, 0,  3,  3,  6/),"     ") , &
            spgr_info_type(230,"I A 3 D     ","-I 4bd 2c 3     ",14,39, (/-3,-3, 0,  3,  3,  6/),"     ") /)

       return
    End Subroutine Set_Spgr_Info

    !!----
    !!---- Subroutine Set_System_Equiv()
    !!----
    !!----    Conversion Table    IT - ML - Kov - BC - Zak
    !!----
    !!--..   The information given in this file corresponds to that of TABLE 6 of
    !!--..   "Isotropy Subgroups of the 230 Crystallographic Space Groups", by
    !!--..   Harold T Stokes and Dorian M Hatch, World Scientific, Singapore (1988).
    !!--..
    !!--..   The transformation operators that take space group elements in the
    !!--..   International setting (International Tables of Crystallography, Hahn 1983)
    !!--..   to space-groups elements in the Miller and Love ( ML, 1967), Kovalev
    !!--..   (Kov,1986) Bradley anb Cracknell (BC, 1972) and Zak (Zak, 1969) settings.
    !!--..
    !!--..   In the international setting the basis vectors are always those of the
    !!--..   conventional unit cell. In the Trigonal system the primitive basis
    !!--..   vectors are in an obverse relationship given by (2/3 1/3 1/3),
    !!--..   (-1/3 1/3 1/3) and (-1/3, -2/3 1/3).
    !!--..   In ML the same basis vectors are chosen except that for trigonal/rhombohedral
    !!--..   system the reverse setting is adopted, so the primitive basis vectors
    !!--..   are: t1=(1/3 -1/3 1/3), t2=(1/3, 2/3 1/3) and t3=(2/3 1/3 1/3)
    !!--..   In Kovalev the a,b,c axes of the coordinate system are along the
    !!--..   conventional basis vectors of the lattice, however in the trigonal
    !!--..   system an hexagonal system is chosen so that the primitive basis vectors
    !!--..   are a1=(-1 -1 1/3), a2=(1 0 1/3) and a3=(0 1 1/3).
    !!--..   In the setting of BC the axes a,b,c of the coordinate system are chosen
    !!--..   to be the primitive basis vectors t1,t2,t3 as defined in their book.
    !!--..   The setting of Zak the basis vectors are as in the international setting,
    !!--..   but for trigonal/rhombohedral system the primitive basis vectors w.r.t. the selected
    !!--..   hexagonal coordinate system are given by: (1/3 2/3 1) (1/3 -1/3 1)
    !!--..   (-2/3 -1/3 1)
    !!--..
    !!--..   Symmetry and transformation operators of Space Groups can be given as
    !!--..   4 x 4 Seitz matrices or as a character string called Jones Faithful
    !!--..   representation. This last representation is that used in this file.
    !!--..
    !!--..   To transform a symmetry operator "gI" in the international setting into
    !!--..   a symmetry element "g" in one of the other settings, we simply perform
    !!--..   the following operation:  g = gT gI gT(-1), where gT is the transformation
    !!--..   given tabulated below.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Set_System_Equiv()

       if (.not. allocated(system_equiv) ) allocate(system_equiv(230))

       system_equiv(1:10) = (/         &
          table_equiv_type("C1_1  ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("C1_i  ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("C1_2  ","z,x,y            ","-z,x,-y           ",        &
                      "-x,z,y                          "," z,x,y            "), &
          table_equiv_type("C2_2  ","z,x,y            ","-z,x,-y           ",        &
                      "-x,z,y                          "," z,x,y            "), &
          table_equiv_type("C3_2  ","z,x,y            ","-z,x,-y           ",        &
                      " z,-x+y,-x-y                    ","-x,z,y            "), &
          table_equiv_type("C1_s  ","z,x,y            ","-z,x,-y           ",        &
                      "-x,z,y                          "," z,x,y            "), &
          table_equiv_type("C2_s  ","z,x,y            ","-z,x,-y           ",        &
                      " z,-x,-y                        "," z,x,y            "), &
          table_equiv_type("C3_s  ","z,x,y            ","-z,x,-y           ",        &
                      " z,-x+y,-x-y                    ","-x,z,y            "), &
          table_equiv_type("C4_s  ","z,x,y            ","-z,x,-y           ",        &
                      " z,-x+y,-x-y                    ","-x,z,y            "), &
          table_equiv_type("C1_2h ","z,x,y            ","-z,x,-y           ",        &
                      "-x,z,y                          "," z,x,y            ") /)

       system_equiv(11:20)= (/         &
          table_equiv_type("C2_2h ","z,x,y            ","-z,x,-y+1/4       ",        &
                      "-x,z,y+1/4                      "," z,x,y            "), &
          table_equiv_type("C3_2h ","z,x,y            ","-z,x,-y           ",        &
                      " z,-x+y,-x-y                    ","-x,z,y            "), &
          table_equiv_type("C4_2h ","z,x,y            ","-z+1/4,x,-y       ",        &
                      " z-1/4,-x,-y                    ","-x,z,y            "), &
          table_equiv_type("C5_2h ","z,x,y            ","-z+1/4,x,-y+1/4   ",        &
                      " z-1/4,-x,-y+1/4                ","-x,z,y            "), &
          table_equiv_type("C6_2h ","z,x,y            ","-z+1/4,x,-y       ",        &
                      " z-1/4,-x+y,-x-y                ","-x,z,y            "), &
          table_equiv_type("D1_2  ","x,y,z            "," x,y,z            ",        &
                      "-y,x,z                          "," x,y,z            "), &
          table_equiv_type("D2_2  ","x,y,z            "," x,y,z            ",        &
                      "-y,x,z+1/4                      "," x,y,z            "), &
          table_equiv_type("D3_2  ","x,y,z            "," x,y,z            ",        &
                      "-y,x,z                          "," x,y,z            "), &
          table_equiv_type("D4_2  ","x,y,z            "," x,y,z            ",        &
                      "-x,-y,z                         "," x,y,z            "), &
          table_equiv_type("D5_2  ","x,y,z            "," x,y,z+1/4        ",        &
                      " x-y,x+y,z                      "," x,y,z            ")/)

       system_equiv(21:30)= (/         &
          table_equiv_type("D6_2  ","x,y,z            "," x,y,z            ",        &
                      " x-y,x+y,z                      "," x,y,z            "), &
          table_equiv_type("D7_2  ","x,y,z            "," x,y,z            ",        &
                      " x+y+z,-x-y+z,x-y-z             "," x,y,z            "), &
          table_equiv_type("D8_2  ","x,y,z            "," x,y,z            ",        &
                      " x+z,-y+z,x-y                   "," x,y,z            "), &
          table_equiv_type("D9_2  ","x,y,z            "," x,y,z            ",        &
                      "-y+z,-x+z,-x-y                  "," x,y,z            "), &
          table_equiv_type("C1_2v ","x,y,z            "," x,y,z            ",        &
                      "-y,x,z                          "," x,y,z            "), &
          table_equiv_type("C2_2v ","x,y,z            "," y,x,z            ",        &
                      "-y,x,z                          "," x,y,z            "), &
          table_equiv_type("C3_2v ","x,y,z            "," x,y,z            ",        &
                      "-y,x,z                          "," x,y,z            "), &
          table_equiv_type("C4_2v ","x,y,z            "," x+1/4,y,z        ",        &
                      "-x-1/4,-y,z                     "," x,y,z            "), &
          table_equiv_type("C5_2v ","x,y,z            "," x+1/4,y,z        ",        &
                      "-x-1/4,-y,z                     "," x,y,z            "), &
          table_equiv_type("C6_2v ","x,y,z            "," y+1/4,x,z        ",        &
                      "-y-1/4,x,z                      "," x,y,z            ") /)

       system_equiv(31:40)= (/         &
          table_equiv_type("C7_2v ","x,y,z            "," x,y,z            ",        &
                      "-x,-y,z                         "," x,y,z            "), &
          table_equiv_type("C8_2v ","x,y,z            "," x+1/4,y+1/4,z    ",        &
                      "-y-1/4,x+1/4,z                  "," x,y,z            "), &
          table_equiv_type("C9_2v ","x,y,z            "," x+1/4,y+1/4,z    ",        &
                      "-x-1/4,-y+1/4,z                 "," x,y,z            "), &
          table_equiv_type("C10_2v","x,y,z            "," x+1/4,y+1/4,z    ",        &
                      "-y-1/4,x+1/4,z                  "," x,y,z            "), &
          table_equiv_type("C11_2v","x,y,z            "," x,y,z            ",        &
                      " x-y,x+y,z                      "," x,y,z            "), &
          table_equiv_type("C12_2v","x,y,z            "," x,y,z            ",        &
                      "-x-y,x-y,z                      "," x,y,z            "), &
          table_equiv_type("C13_2v","x,y,z            "," x,y,z            ",        &
                      " x-y,x+y,z                      "," x,y,z            "), &
          table_equiv_type("C14_2v","-z,y,x           "," -z,y,x           ",        &
                      "-y+z,-y-z,x                     ","-y,-z,x           "), &
          table_equiv_type("C15_2v","-z,y,x           "," -z,y,x           ",        &
                      "-y+z,-y-z,x                     ","-y,-z,x           "), &
          table_equiv_type("C16_2v","-z,y,x           "," -z,y,x           ",        &
                      "-y+z,-y-z,x                     ","-y,-z,x           ") /)

       system_equiv(41:50)= (/         &
          table_equiv_type("C17_2v","-z,y,x           "," -z,y,x           ",        &
                      "-y+z,-y-z,x                     ","-y,-z,x           "), &
          table_equiv_type("C18_2v","x,y,z            "," x,y,z            ",        &
                      " x+y+z,-x-y+z,x-y-z             "," x,y,z            "), &
          table_equiv_type("C19_2v","x,y,z            "," x-1/8,y-1/8,z    ",        &
                      " x+y+z+1/2,-x-y+z-1/2,x-y-z-1/4 "," x,y,z            "), &
          table_equiv_type("C20_2v","x,y,z            "," x,y,z            ",        &
                      " x+z,-y+z,x-y                   "," x,y,z            "), &
          table_equiv_type("C21_2v","x,y,z            "," x,y,z            ",        &
                      " x+z,-y+z,x-y                   "," x,y,z            "), &
          table_equiv_type("C22_2v","x,y,z            "," x,y,z            ",        &
                      "-y+z,-x+z,-x-y                  "," x,y,z            "), &
          table_equiv_type("D1_2h ","x,y,z            "," x,y,z            ",        &
                      "-y,x,z                          "," x,y,z            "), &
          table_equiv_type("D2_2h ","x-1/4,y-1/4,z-1/4"," x-1/4,y-1/4,z-1/4",        &
                      "-y+1/4,x-1/4,z-1/4              "," x-1/4,y-1/4,z-1/4"), &
          table_equiv_type("D3_2h ","x,y,z            "," x,y,z+1/4        ",        &
                      "-y,x,z+1/4                      "," x,y,z            "), &
          table_equiv_type("D4_2h ","x-1/4,y-1/4,z    "," x-1/4,y-1/4,z    ",        &
                      "-y+1/4,x-1/4,z                  "," x-1/4,y-1/4,z    ") /)

       system_equiv(51:60)= (/         &
          table_equiv_type("D5_2h ","x,y,z            "," y,z,x            ",        &
                      "-y,z,-x                         "," x,y,z            "), &
          table_equiv_type("D6_2h ","x,y,z            "," z+1/4,x+1/4,y    ",        &
                      " z-1/4,x+1/4,y                  "," x,y,z            "), &
          table_equiv_type("D7_2h ","x,y,z            "," x-1/4,y,z        ",        &
                      "-x-1/4,-y,z                     "," x,y,z            "), &
          table_equiv_type("D8_2h ","x,y,z            "," y,z+1/4,x        ",        &
                      "-y,z+1/4,-x                     "," x,y,z            "), &
          table_equiv_type("D9_2h ","x,y,z            "," x,y,z            ",        &
                      "-y,x,z                          "," x,y,z            "), &
          table_equiv_type("D10_2h","x,y,z            "," x+1/4,y+1/4,z+1/4",        &
                      "-y-1/4,x+1/4,z+1/4              "," x,y,z            "), &
          table_equiv_type("D11_2h","x,y,z            "," -z,-y-1/4,-x     ",        &
                      "-z,y+1/4,x                      "," x,y,z            "), &
          table_equiv_type("D12_2h","x,y,z            "," x,y,z-1/4        ",        &
                      "-y,x,z+1/4                      "," x,y,z            "), &
          table_equiv_type("D13_2h","x-1/4,y-1/4,z    "," x-1/4,y-1/4,z    ",        &
                      "-y+1/4,x-1/4,z                  "," x-1/4,y-1/4,z+1/4"), &
          table_equiv_type("D14_2h","x,y,z            "," z+1/4,x,y+1/4    ",        &
                      " z-1/4,x,y+1/4                  "," x,y,z            ") /)

       system_equiv(61:70)= (/         &
          table_equiv_type("D15_2h","x,y,z            "," x,y,z            ",        &
                      "-x,-y,z                         "," x,y,z            "), &
          table_equiv_type("D16_2h","x,y,z            "," y+1/4,x+1/4,z    ",        &
                      "-y-1/4,x+1/4,z                  "," x,y,z            "), &
          table_equiv_type("D17_2h","x,y,z            "," y,x,z            ",        &
                      " x-y,x+y,z                      "," x,y,z            "), &
          table_equiv_type("D18_2h","x,y,z            "," y,x+1/4,z        ",        &
                      " x-y+1/4,x+y+1/4,z              "," x,y,z            "), &
          table_equiv_type("D19_2h","x,y,z            "," x,y,z            ",        &
                      " x-y,x+y,z                      "," x,y,z            "), &
          table_equiv_type("D20_2h","x,y,z            "," x,y,z+1/4        ",        &
                      " x-y,x+y,z+1/4                  "," x,y,z            "), &
          table_equiv_type("D21_2h","x,y,z            "," x+1/4,y,z        ",        &
                      " x-y+1/4,x+y+1/4,z              "," x,y,z            "), &
          table_equiv_type("D22_2h","x,y-1/4,z-1/4    "," x,y-1/4,z-1/4    ",        &
                      " x-y+1/4,x+y-1/4,z-1/4          "," x,y-1/4,z-1/4    "), &
          table_equiv_type("D23_2h","x,y,z            "," x,y,z            ",        &
                      " x+y+z,-x-y+z,x-y-z             "," x,y,z            "), &
          table_equiv_type("D24_2h","x-7/8,y-7/8,z-7/8"," x-7/8,y-7/8,z-7/8",        &
                      " x+y+z-15/8,-x-y+z+5/8,x-y-z+5/8"," x-7/8,y-7/8,z-7/8") /)

       system_equiv(71:80)= (/         &
          table_equiv_type("D25_2h","x,y,z            "," x,y,z            ",        &
                      " x+z,-y+z,x-y                   "," x,y,z            "), &
          table_equiv_type("D26_2h","x,y,z            "," x,y,z-1/4        ",        &
                      " x+z+1/4,-y+z+1/4,x-y           "," x,y,z            "), &
          table_equiv_type("D27_2h","x,y,z            "," x,y,z            ",        &
                      " x+z+1/2,-y+z,x-y               "," x,y,z            "), &
          table_equiv_type("D28_2h","x,y,z            "," x,y,z+1/4        ",        &
                      " x+z+1/4,-y+z-1/4,x-y           "," x,y,z            "), &
          table_equiv_type("C1_4  ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("C2_4  ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("C3_4  ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("C4_4  ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("C5_4  ","x,y,z            "," x,y,z            ",        &
                      " y+z,x+z,x+y                    "," x,y,z            "), &
          table_equiv_type("C6_4  ","x,y,z            "," x,y,z            ",        &
                      " y+z,x+z,x+y                    "," x,y,z            ") /)

       system_equiv(81:90)= (/         &
          table_equiv_type("S1_4  ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("S2_4  ","x,y,z            "," x,y,z            ",        &
                      " y+z,x+z,x+y                    "," x,y,z            "), &
          table_equiv_type("C1_4h ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("C2_4h ","x,y,z            "," x,y,z+1/4        ",        &
                      " x,y,z+1/4                      "," x,y,z            "), &
          table_equiv_type("C3_4h ","x-3/4,y-1/4,z    "," x-3/4,y-1/4,z    ",        &
                      " x-3/4,y-1/4,z                  "," x-3/4,y-1/4,z    "), &
          table_equiv_type("C4_4h ","x-3/4,y-3/4,z-3/4"," x-3/4,y-3/4,z-3/4",        &
                      " x-3/4,y-3/4,z-3/4              "," x-3/4,y-3/4,z-3/4"), &
          table_equiv_type("C5_4h ","x,y,z            "," x,y,z            ",        &
                      " y+z,x+z,x+y                    "," x,y,z            "), &
          table_equiv_type("C6_4h ","x,y-3/4,z-7/8    "," x,y-3/4,z-7/8    ",        &
                      " y+z-13/8,x+z-7/8,x+y-3/4       "," x,y-3/4,z-7/8    "), &
          table_equiv_type("D1_4  ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("D2_4  ","x,y,z            "," x,y-1/2,z        ",        &
                      " x+1/2,y,z                      "," x,y,z            ") /)

       system_equiv(91:100) = (/         &
          table_equiv_type("D3_4  ","x,y,z            "," x,y,z+1/4        ",        &
                      " x,y,z+1/4                      "," x,y,z            "), &
          table_equiv_type("D4_4  ","x,y,z            "," x,y-1/2,z+1/8    ",        &
                      " x+1/2,y,z+1/8                  "," x,y,z            "), &
          table_equiv_type("D5_4  ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("D6_4  ","x,y,z            "," x,y+1/2,z+1/4    ",        &
                      " x+1/2,y,z+1/4                  "," x,y,z            "), &
          table_equiv_type("D7_4  ","x,y,z            "," x,y,z+1/4        ",        &
                      " x,y,z+1/4                      "," x,y,z            "), &
          table_equiv_type("D8_4  ","x,y,z            "," x,y-1/2,z-1/8    ",        &
                      " x+1/2,y,z+3/8                  "," x,y,z            "), &
          table_equiv_type("D9_4  ","x,y,z            "," x,y,z            ",        &
                      " y+z,x+z,x+y                    "," x,y,z            "), &
          table_equiv_type("D10_4 ","x,y,z            "," x,y,z            ",        &
                      " y+z+1/8,x+z+1/8,x+y            "," x,y,z            "), &
          table_equiv_type("C1_4v ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("C2_4v ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            ") /)

       system_equiv(101:110)= (/         &
          table_equiv_type("C3_4v ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("C4_4v ","x,y,z            "," x,y-1/2,z        ",        &
                      " x+1/2,y,z                      "," x,y,z            "), &
          table_equiv_type("C5_4v ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("C6_4v ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("C7_4v ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("C8_4v ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("C9_4v ","x,y,z            "," x,y,z            ",        &
                      " y+z,x+z,x+y                    "," x,y,z            "), &
          table_equiv_type("C10_4v","x,y,z            "," x,y,z            ",        &
                      " y+z,x+z,x+y                    "," x,y,z            "), &
          table_equiv_type("C11_4v","x,y,z            "," x,y,z            ",        &
                      " y+z,x+z,x+y                    "," x,y,z            "), &
          table_equiv_type("C12_4v","x,y,z            "," x,y,z            ",        &
                      " y+z,x+z,x+y                    "," x,y,z            ") /)

       system_equiv(111:120)= (/         &
          table_equiv_type("D1_2d ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("D2_2d ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("D3_2d ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("D4_2d ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("D5_2d ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("D6_2d ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("D7_2d ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("D8_2d ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("D9_2d ","x,y,z            "," x,y,z            ",        &
                      " y+z,x+z,x+y                    "," x,y,z            "), &
          table_equiv_type("D10_2d","x,y,z            "," x,y,z            ",        &
                      " y+z,x+z,x+y                    "," x,y,z            ")  /)

       system_equiv(121:130)= (/         &
          table_equiv_type("D11_2d","x,y,z            "," x,y,z            ",        &
                      " y+z,x+z,x+y                    "," x,y,z            "), &
          table_equiv_type("D12_2d","x,y,z            "," x+1/2,y,z+1/4    ",        &
                      " x+z,-y+z,x-y                   "," x,y,z            "), &
          table_equiv_type("D1_4h ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("D2_4h ","x,y,z            "," x,y,z-1/4        ",        &
                      " x,y,z+1/4                      "," x,y,z            "), &
          table_equiv_type("D3_4h ","x-3/4,y-3/4,z    "," x-3/4,y-3/4,z    ",        &
                      " x-1/4,y-3/4,z                  "," x-3/4,y-3/4,z    "), &
          table_equiv_type("D4_4h ","x-3/4,y-3/4,z-3/4"," x-3/4,y-3/4,z-3/4",        &
                      " x-1/4,y-3/4,z-3/4              "," x-3/4,y-3/4,z-3/4"), &
          table_equiv_type("D5_4h ","x,y,z            "," x,y,z            ",        &
                      " x+1/2,y,z                      "," x,y,z            "), &
          table_equiv_type("D6_4h ","x,y,z            "," x,y,z+1/4        ",        &
                      " x+1/2,y,z+1/4                  "," x,y,z            "), &
          table_equiv_type("D7_4h ","x-3/4,y-1/4,z    "," x-3/4,y+1/4,z    ",        &
                      " x-1/4,y-1/4,z                  "," x-3/4,y-1/4,z    "), &
          table_equiv_type("D8_4h ","x-3/4,y-1/4,z    "," x-3/4,y+1/4,z+1/4",        &
                      " x-1/4,y-1/4,z+1/4              "," x,y,z            ") /)

       system_equiv(131:140)= (/         &
          table_equiv_type("D9_4d ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("D10_4d","x,y,z            "," x,y,z-1/4        ",        &
                      " x,y,z+1/4                      "," x,y,z            "), &
          table_equiv_type("D11_4d","x-3/4,y-1/4,z-3/4"," x-3/4,y+1/4,z-1/2",        &
                      " x-3/4,y-1/4,z-1/2              "," x,y,z            "), &
          table_equiv_type("D12_4d","x-3/4,y-1/4,z-3/4"," x-3/4,y+1/4,z-3/4",        &
                      " x-3/4,y-1/4,z-3/4              "," x-3/4,y-1/4,z-3/4"), &
          table_equiv_type("D13_4d","x,y,z            "," x,y,z            ",        &
                      " x+1/2,y,z                      "," x,y,z            "), &
          table_equiv_type("D14_4d","x,y,z            "," x,y+1/2,z+1/4    ",        &
                      " x,y,z+1/4                      "," x+1/2,y,z        "), &
          table_equiv_type("D15_4d","x-3/4,y-1/4,z-3/4"," x-3/4,y+1/4,z-1/2",        &
                      " x-1/4,y-1/4,z-1/2              "," x,y,z            "), &
          table_equiv_type("D16_4d","x-3/4,y-1/4,z-3/4"," x-3/4,y+1/4,z-3/4",        &
                      " x-1/4,y-1/4,z-3/4              "," x,y,z            "), &
          table_equiv_type("D17_4d","x,y,z            "," x,y,z            ",        &
                      " y+z,x+z,x+y                    "," x,y,z            "), &
          table_equiv_type("D18_4d","x,y,z            "," x,y,z+1/4        ",        &
                      " y+z+1/4,x+z+3/4,x+y+1/2        "," x,y,z            ") /)

       system_equiv(141:150)= (/         &
          table_equiv_type("D19_4d","x,y-1/4,z-7/8    "," x,y-1/4,z-7/8    ",        &
                      " y+z-3/4,x+z-3/4,x+y            "," x,y-1/4,z-7/8    "), &
          table_equiv_type("D20_4d","x,y-1/4,z-7/8    "," x,y-1/4,z-9/8    ",        &
                      " y+z,x+z-1/2,x+y+1/2            "," x,y,z            "), &
          table_equiv_type("C1_3  ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("C2_3  ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("C3_3  ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("C4_3  ","y,-x+y,z         ","-2x+y,-x-y,z      ",        &
                      " x+z,-x+y+z,-y+z                "," y,-x+y,3z        "), &
          table_equiv_type("C1_3i ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("C2_3i ","y,-x+y,z         ","-2x+y,-x-y,z      ",        &
                      " x+z,-x+y+z,-y+z                "," y,-x+y,3z        "), &
          table_equiv_type("D1_3  ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("D2_3  ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            ") /)

       system_equiv(151:160)= (/         &
          table_equiv_type("D3_3  ","x,y,z-1/6        "," x,y,z+1/6        ",        &
                      " x,y,z+1/6                      "," x,y,z+1/6        "), &
          table_equiv_type("D4_3  ","x,y,z            "," x,y,z+1/3        ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("D5_3  ","x,y,z-5/6        "," x,y,z+1/3        ",        &
                      " x,y,z-1/6                      "," x,y,z-1/6        "), &
          table_equiv_type("D6_3  ","x,y,z-1/6        "," x,y,z+1/6        ",        &
                      " x,y,z+1/2                      "," x,y,z+1/2        "), &
          table_equiv_type("D7_3  ","y,-x+y,z         ","-2x+y,-x-y,z      ",        &
                      " x+z,-x+y+z,-y+z                "," y,-x+y,3z        "), &
          table_equiv_type("C1_3v ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("C2_3v ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("C3_3v ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("C4_3v ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("C5_3v ","y,-x+y,z         ","-2x+y,-x-y,z      ",        &
                      " x+z,-x+y+z,-y+z                "," y,-x+y,3z        ")/)

       system_equiv(161:170)= (/         &
          table_equiv_type("C6_3v ","y,-x+y,z         ","-2x+y,-x-y,z      ",        &
                      " x+z,-x+y+z,-y+z                "," y,-x+y,3z        "), &
          table_equiv_type("D1_3d ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("D2_3d ","x,y,z            "," x,y,z+1/4        ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("D3_3d ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("D4_3d ","x,y,z            "," x,y,z+1/4        ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("D5_3d ","y,-x+y,z         ","-2x+y,-x-y,z      ",        &
                      " x+z,-x+y+z,-y+z                "," y,-x+y,3z        "), &
          table_equiv_type("D6_3d ","y,-x+y,z         ","-2x+y,-x-y,z      ",        &
                      " x+z,-x+y+z,-y+z                "," y,-x+y,3z        "), &
          table_equiv_type("C1_6  ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("C2_6  ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("C3_6  ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            ") /)

       system_equiv(171:180)= (/         &
          table_equiv_type("C4_6  ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("C5_6  ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("C6_6  ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("C7_6  ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("C1_6h ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("C2_6h ","x,y,z            "," x,y,z            ",        &
                      " x,y,z+1/4                      "," x,y,z            "), &
          table_equiv_type("D1_6  ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("D2_6  ","x,y,z            "," x,y,z+1/6        ",        &
                      " x,y,z+1/4                      "," x,y,z            "), &
          table_equiv_type("D3_6  ","x,y,z            "," x,y,z+1/3        ",        &
                      " x,y,z+1/4                      "," x,y,z            "), &
          table_equiv_type("D4_6  ","x,y,z            "," x,y,z+1/3        ",        &
                      " x,y,z                          "," x,y,z            ") /)

       system_equiv(181:190)= (/         &
          table_equiv_type("D5_6  ","x,y,z            "," x,y,z-1/3        ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("D6_6  ","x,y,z            "," x,y,z            ",        &
                      " x,y,z+1/4                      "," x,y,z            "), &
          table_equiv_type("C1_6v ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("C2_6v ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("C3_6v ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("C4_6v ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("D1_3h ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("D2_3h ","x,y,z            "," x,y,z            ",        &
                      " x,y,z+1/4                      "," x,y,z            "), &
          table_equiv_type("D3_3h ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("D4_3h ","x,y,z            "," x,y,z            ",        &
                      " x,y,z+1/4                      "," x,y,z            ") /)

       system_equiv(191:200)= (/         &
          table_equiv_type("D1_6h ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("D2_6h ","x,y,z            "," x,y,z-1/4        ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("D3_6h ","x,y,z            "," x,y,z-1/4        ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("D4_6h ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("T1    ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("T2    ","x,y,z            "," x,y,z            ",        &
                      "-x+y+z,x-y+z,x+y-z              "," x,y,z            "), &
          table_equiv_type("T3    ","x,y,z            "," x,y,z            ",        &
                      " y+z,x+z,x+y                    "," x,y,z            "), &
          table_equiv_type("T4    ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("T5    ","x,y,z            "," x,y,z            ",        &
                      " y+z,x+z,x+y                    "," x,y,z            "), &
          table_equiv_type("T1_h  ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            ") /)

       system_equiv(201:210)= (/         &
          table_equiv_type("T2_h  ","x-3/4,y-3/4,z-3/4"," x-3/4,y-3/4,z-3/4",        &
                      " x-3/4,y-3/4,z-3/4              "," x-3/4,y-3/4,z-3/4"), &
          table_equiv_type("T3_h  ","x,y,z            "," x,y,z            ",        &
                      "-x+y+z,x-y+z,x+y-z              "," x,y,z            "), &
          table_equiv_type("T4_h  ","x-7/8,y-7/8,z-7/8"," x-7/8,y-7/8,z-7/8",        &
                      "-x+y+z-7/8,x-y+z-7/8,x+y-z-7/8  "," x-7/8,y-7/8,z-7/8"), &
          table_equiv_type("T5_h  ","x,y,z            "," x,y,z            ",        &
                      " y+z,x+z,x+y                    "," x,y,z            "), &
          table_equiv_type("T6_h  ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("T7_h  ","x,y,z            "," x,y,z            ",        &
                      " y+z,x+z,x+y                    "," x,y,z            "), &
          table_equiv_type("O1    ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("O2    ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("O3    ","x,y,z            "," x,y,z            ",        &
                      "-x+y+z,x-y+z,x+y-z              "," x,y,z            "), &
          table_equiv_type("O4    ","x,y,z            "," x,y,z            ",        &
                      "-x+y+z,x-y+z,x+y-z              "," x,y,z            ") /)

       system_equiv(211:220)= (/         &
          table_equiv_type("O5    ","x,y,z            "," x,y,z            ",        &
                      " y+z,x+z,x+y                    "," x,y,z            "), &
          table_equiv_type("O6    ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("O7    ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("O8    ","x,y,z            "," x,y,z            ",        &
                      " y+z,x+z,x+y                    "," x,y,z            "), &
          table_equiv_type("T1_d  ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("T2_d  ","x,y,z            "," x,y,z            ",        &
                      "-x+y+z,x-y+z,x+y-z              "," x,y,z            "), &
          table_equiv_type("T3_d  ","x,y,z            "," x,y,z            ",        &
                      " y+z,x+z,x+y                    "," x,y,z            "), &
          table_equiv_type("T4_d  ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("T5_d  ","x,y,z            "," x,y,z            ",        &
                      "-x+y+z,x-y+z,x+y-z              "," x,y,z            "), &
          table_equiv_type("T6_d  ","x,y,z            "," x,y,z            ",        &
                      " y+z,x+z,x+y                    "," x,y,z            ") /)

       system_equiv(221:230)= (/         &
          table_equiv_type("O1_h  ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("O2_h  ","x-3/4,y-3/4,z-3/4"," x-3/4,y-3/4,z-3/4",        &
                      " x-3/4,y-3/4,z-3/4              "," x-3/4,y-3/4,z-3/4"), &
          table_equiv_type("O3_h  ","x,y,z            "," x,y,z            ",        &
                      " x,y,z                          "," x,y,z            "), &
          table_equiv_type("O4_h  ","x-3/4,y-3/4,z-3/4"," x-3/4,y-3/4,z-3/4",        &
                      " x-3/4,y-3/4,z-3/4              "," x-3/4,y-3/4,z-3/4"), &
          table_equiv_type("O5_h  ","x,y,z            "," x,y,z            ",        &
                      "-x+y+z,x-y+z,x+y-z              "," x,y,z            "), &
          table_equiv_type("O6_h  ","x,y,z            "," x-1/4,y-1/4,z-1/4",        &
                      "-x+y+z+1/4,x-y+z+1/4,x+y-z+1/4  "," x,y,z            "), &
          table_equiv_type("O7_h  ","x,y,z            "," x+1/8,y+1/8,z+1/8",        &
                      "-x+y+z+1/8,x-y+z+1/8,x+y-z+1/8  "," x+1/8,y+1/8,z+1/8"), &
          table_equiv_type("O8_h  ","x,y,z            "," x+3/8,y+3/8,z+3/8",        &
                      "-x+y+z+3/8,x-y+z+3/8,x+y-z+3/8  "," x+3/8,y+3/8,z+3/8"), &
          table_equiv_type("O9_h  ","x,y,z            "," x,y,z            ",        &
                      " y+z,x+z,x+y                    "," x,y,z            "), &
          table_equiv_type("O10_h ","x,y,z            "," x,y,z            ",        &
                      " y+z,x+z,x+y                    "," x,y,z            ") /)

       return
    End Subroutine Set_System_Equiv

    !!----
    !!---- Subroutine Set_Wyckoff_Info()
    !!----
    !!----    Set Information on Wyckoff_info array
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Set_Wyckoff_Info()

       if (.not. allocated(wyckoff_info) ) allocate(wyckoff_info(273) )

       wyckoff_info(  1)= wyck_info_type("P 1         ", 0,     &
                    (/"               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(  2)= wyck_info_type("P -1        ", 8,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,1/2,0        ",    &
                      "1/2,0,0        ", "1/2,1/2,0      ", "1/2,0,1/2      ",    &
                      "0,1/2,1/2      ", "1/2,1/2,1/2    ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(  3)= wyck_info_type("P 2         ", 4,     &
                    (/"0,y,0          ", "0,y,1/2        ", "1/2,y,0        ",    &
                      "1/2,y,1/2      ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(  4)= wyck_info_type("P 21        ", 0,     &
                    (/"               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(  5)= wyck_info_type("C 2         ", 2,     &
                    (/"0,y,0          ", "0,y,1/2        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(  6)= wyck_info_type("A 2         ", 2,     &
                    (/"0,y,0          ", "1/2,y,1/2      ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(  7)= wyck_info_type("I 2         ", 2,     &
                    (/"0,y,0          ", "1/2,y,0        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(  8)= wyck_info_type("P M         ", 2,     &
                    (/"x,0,z          ", "x,1/2,z        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(  9)= wyck_info_type("P C         ", 0,     &
                    (/"               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 10)= wyck_info_type("C M         ", 1,     &
                    (/"x,0,z          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 11)= wyck_info_type("A M         ", 1,     &
                    (/"x,0,z          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 12)= wyck_info_type("I M         ", 1,     &
                    (/"x,0,z          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 13)= wyck_info_type("C C         ", 0,     &
                    (/"               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 14)= wyck_info_type("P 2/M       ",14,     &
                    (/"0,0,0          ", "0,1/2,0        ", "0,0,1/2        ",    &
                      "1/2,0,0        ", "1/2,1/2,0      ", "0,1/2,1/2      ",    &
                      "1/2,0,1/2      ", "1/2,1/2,1/2    ", "0,y,0          ",    &
                      "1/2,y,0        ", "0,y,1/2        ", "1/2,y,1/2      ",    &
                      "x,0,z          ", "x,1/2,z        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 15)= wyck_info_type("P 21/M      ", 5,     &
                    (/"0,0,0          ", "1/2,0,0        ", "0,0,1/2        ",    &
                      "1/2,0,1/2      ", "x,1/4,z        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 16)= wyck_info_type("C 2/M       ", 9,     &
                    (/"0,0,0          ", "0,1/2,0        ", "0,0,1/2        ",    &
                      "0,1/2,1/2      ", "1/4,1/4,0      ", "1/4,1/4,1/2    ",    &
                      "0,y,0          ", "0,y,1/2        ", "x,0,z          ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 17)= wyck_info_type("A 2/M       ", 9,     &
                    (/"0,0,0          ", "0,1/2,0        ", "1/2,0,1/2      ",    &
                      "1/2,1/2,1/2    ", "0,1/4,1/4      ", "1/2,1/4,3/4    ",    &
                      "0,y,0          ", "1/2,y,1/2      ", "x,0,z          ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 18)= wyck_info_type("I 2/M       ", 9,     &
                    (/"0,0,0          ", "0,1/2,0        ", "1/2,0,0        ",    &
                      "1/2,1/2,0      ", "3/4,1/4,3/4    ", "1/4,1/4,3/4    ",    &
                      "0,y,0          ", "1/2,y,0        ", "x,0,z          ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 19)= wyck_info_type("P 2/C       ", 6,     &
                    (/"0,0,0          ", "1/2,1/2,0      ", "0,1/2,0        ",    &
                      "1/2,0,0        ", "0,y,1/4        ", "1/2,y,1/4      ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 20)= wyck_info_type("P 2/N       ", 6,     &
                    (/"0,0,0          ", "0,1/2,1/2      ", "0,1/2,0        ",    &
                      "0,0,1/2        ", "3/4,y,3/4      ", "3/4,y,1/4      ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 21)= wyck_info_type("P 2/A       ", 6,     &
                    (/"0,0,0          ", "1/2,1/2,1/2    ", "0,1/2,0        ",    &
                      "1/2,0,1/2      ", "1/4,y,0        ", "3/4,y,1/2      ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 22)= wyck_info_type("P 21/C      ", 4,     &
                    (/"0,0,0          ", "1/2,0,0        ", "0,0,1/2        ",    &
                      "1/2,0,1/2      ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 23)= wyck_info_type("P 21/N      ", 4,     &
                    (/"0,0,0          ", "0,0,1/2        ", "1/2,0,1/2      ",    &
                      "1/2,0,0        ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 24)= wyck_info_type("P 21/A      ", 4,     &
                    (/"0,0,0          ", "1/2,0,1/2      ", "1/2,0,0        ",    &
                      "0,0,1/2        ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 25)= wyck_info_type("C 2/C       ", 5,     &
                    (/"0,0,0          ", "0,1/2,0        ", "1/4,1/4,0      ",    &
                      "1/4,1/4,1/2    ", "0,y,1/4        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 26)= wyck_info_type("A 2/N       ", 5,     &
                    (/"0,0,0          ", "0,1/2,0        ", "0,1/4,1/4      ",    &
                      "1/2,1/4,3/4    ", "3/4,y,3/4      ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 27)= wyck_info_type("I 2/A       ", 5,     &
                    (/"0,0,0          ", "0,1/2,0        ", "3/4,1/4,3/4    ",    &
                      "1/4,1/4,3/4    ", "1/4,y,0        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 28)= wyck_info_type("P 2 2 2     ",20,     &
                    (/"0,0,0          ", "1/2,0,0        ", "0,1/2,0        ",    &
                      "0,0,1/2        ", "1/2,1/2,0      ", "1/2,0,1/2      ",    &
                      "0,1/2,1/2      ", "1/2,1/2,1/2    ", "x,0,0          ",    &
                      "x,0,1/2        ", "x,1/2,0        ", "x,1/2,1/2      ",    &
                      "0,y,0          ", "0,y,1/2        ", "1/2,y,0        ",    &
                      "1/2,y,1/2      ", "0,0,z          ", "1/2,0,z        ",    &
                      "0,1/2,z        ", "1/2,1/2,z      ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 29)= wyck_info_type("P 2 2 21    ", 4,     &
                    (/"x,0,0          ", "x,1/2,0        ", "0,y,1/4        ",    &
                      "1/2,y,1/4      ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 30)= wyck_info_type("P 21 21 2   ", 2,     &
                    (/"0,0,z          ", "0,1/2,z        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 31)= wyck_info_type("P 21 21 21  ", 0,     &
                    (/"               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 32)= wyck_info_type("C 2 2 21    ", 2,     &
                    (/"x,0,0          ", "0,y,1/4        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 33)= wyck_info_type("C 2 2 2     ",11,     &
                    (/"0,0,0          ", "0,1/2,0        ", "1/2,0,1/2      ",    &
                      "0,0,1/2        ", "x,0,0          ", "x,0,1/2        ",    &
                      "0,y,0          ", "0,y,1/2        ", "0,0,z          ",    &
                      "0,1/2,z        ", "1/4,1/4,z      ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 34)= wyck_info_type("F 2 2 2     ",10,     &
                    (/"0,0,0          ", "0,0,1/2        ", "1/4,1/4,1/4    ",    &
                      "1/4,1/4,3/4    ", "x,0,0          ", "0,y,0          ",    &
                      "0,0,z          ", "1/4,1/4,z      ", "1/4,y,1/4      ",    &
                      "x,1/4,1/4      ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 35)= wyck_info_type("I 2 2 2     ",10,     &
                    (/"0,0,0          ", "1/2,0,0        ", "0,0,1/2        ",    &
                      "0,1/2,0        ", "x,0,0          ", "x,0,1/2        ",    &
                      "0,y,0          ", "1/2,y,0        ", "0,0,z          ",    &
                      "0,1/2,z        ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 36)= wyck_info_type("I 21 21 21  ", 3,     &
                    (/"x,0,1/4        ", "1/4,y,0        ", "0,1/4,z        ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 37)= wyck_info_type("P M M 2     ", 8,     &
                    (/"0,0,z          ", "0,1/2,z        ", "1/2,0,z        ",    &
                      "1/2,1/2,z      ", "x,0,z          ", "x,1/2,z        ",    &
                      "0,y,z          ", "1/2,y,z        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 38)= wyck_info_type("P M C 21    ", 2,     &
                    (/"0,y,z          ", "1/2,y,z        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 39)= wyck_info_type("P C C 2     ", 4,     &
                    (/"0,0,z          ", "0,1/2,z        ", "1/2,0,z        ",    &
                      "1/2,1/2,z      ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 40)= wyck_info_type("P M A 2     ", 3,     &
                    (/"0,0,z          ", "0,1/2,z        ", "1/4,y,z        ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 41)= wyck_info_type("P C A 21    ", 0,     &
                    (/"               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 42)= wyck_info_type("P N C 2     ", 2,     &
                    (/"0,0,z          ", "1/2,0,z        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 43)= wyck_info_type("P M N 21    ", 1,     &
                    (/"0,y,z          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 44)= wyck_info_type("P B A 2     ", 2,     &
                    (/"0,0,z          ", "0,1/2,z        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 45)= wyck_info_type("P N A 21    ", 0,     &
                    (/"               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 46)= wyck_info_type("P N N 2     ", 2,     &
                    (/"0,0,z          ", "0,1/2,z        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 47)= wyck_info_type("C M M 2     ", 5,     &
                    (/"0,0,z          ", "0,1/2,z        ", "1/4,1/4,z      ",    &
                      "x,0,z          ", "0,y,z          ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 48)= wyck_info_type("C M C 21    ", 1,     &
                    (/"0,y,z          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 49)= wyck_info_type("C C C 2     ", 3,     &
                    (/"0,0,z          ", "0,1/2,z        ", "1/4,1/4,z      ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 50)= wyck_info_type("A M M 2     ", 5,     &
                    (/"0,0,z          ", "1/2,0,z        ", "x,0,z          ",    &
                      "0,y,z          ", "1/2,y,z        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 51)= wyck_info_type("A B M 2     ", 3,     &
                    (/"0,0,z          ", "1/2,0,z        ", "x,1/4,z        ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 52)= wyck_info_type("A M A 2     ", 2,     &
                    (/"0,0,z          ", "1/4,y,z        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 53)= wyck_info_type("A B A 2     ", 1,     &
                    (/"0,0,z          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 54)= wyck_info_type("F M M 2     ", 4,     &
                    (/"0,0,z          ", "1/4,1/4,z      ", "0,y,z          ",    &
                      "x,0,z          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 55)= wyck_info_type("F D D 2     ", 1,     &
                    (/"0,0,z          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 56)= wyck_info_type("I M M 2     ", 4,     &
                    (/"0,0,z          ", "0,1/2,z        ", "x,0,z          ",    &
                      "0,y,z          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 57)= wyck_info_type("I B A 2     ", 2,     &
                    (/"0,0,z          ", "0,1/2,z        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 58)= wyck_info_type("I M A 2     ", 2,     &
                    (/"0,0,z          ", "1/4,y,z        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 59)= wyck_info_type("P M M M     ",26,     &
                    (/"0,0,0          ", "1/2,0,0        ", "0,0,1/2        ",    &
                      "1/2,0,1/2      ", "0,1/2,0        ", "1/2,1/2,0      ",    &
                      "0,1/2,1/2      ", "1/2,1/2,1/2    ", "x,0,0          ",    &
                      "x,0,1/2        ", "x,1/2,0        ", "x,1/2,1/2      ",    &
                      "0,y,0          ", "0,y,1/2        ", "1/2,y,0        ",    &
                      "1/2,y,1/2      ", "0,0,z          ", "0,1/2,z        ",    &
                      "1/2,0,z        ", "1/2,1/2,z      ", "0,y,z          ",    &
                      "1/2,y,z        ", "x,0,z          ", "x,1/2,z        ",    &
                      "x,y,0          ", "x,y,1/2        "/) )
       wyckoff_info( 60)= wyck_info_type("P N N N:1   ",12,     &
                    (/"0,0,0          ", "1/2,0,0        ", "0,0,1/2        ",    &
                      "0,1/2,0        ", "1/4,1/4,1/4    ", "3/4,3/4,3/4    ",    &
                      "x,0,0          ", "x,0,1/2        ", "0,y,0          ",    &
                      "1/2,y,0        ", "0,0,z          ", "0,1/2,z        ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 61)= wyck_info_type("P N N N     ",12,     &
                    (/"1/4,1/4,1/4    ", "3/4,1/4,1/4    ", "1/4,1/4,3/4    ",    &
                      "1/4,3/4,1/4    ", "1/2,1/2,1/2    ", "0,0,0          ",    &
                      "x,1/4,1/4      ", "x,1/4,3/4      ", "1/4,y,1/4      ",    &
                      "3/4,y,1/4      ", "1/4,1/4,z      ", "1/4,3/4,z      ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 62)= wyck_info_type("P C C M     ",17,     &
                    (/"0,0,0          ", "1/2,1/2,0      ", "0,1/2,0        ",    &
                      "1/2,0,0        ", "0,0,1/4        ", "1/2,0,1/4      ",    &
                      "0,1/2,1/4      ", "1/2,1/2,1/4    ", "x,0,1/4        ",    &
                      "x,1/2,1/4      ", "0,y,1/4        ", "1/2,y,1/4      ",    &
                      "0,0,z          ", "1/2,1/2,z      ", "0,1/2,z        ",    &
                      "1/2,0,z        ", "x,y,0          ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 63)= wyck_info_type("P B A N:1   ",12,     &
                    (/"0,0,0          ", "1/2,0,0        ", "1/2,0,1/2      ",    &
                      "0,0,1/2        ", "1/4,1/4,0      ", "1/4,1/4,1/2    ",    &
                      "x,0,0          ", "x,0,1/2        ", "0,y,0          ",    &
                      "0,y,1/2        ", "0,0,z          ", "0,1/2,z        ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 64)= wyck_info_type("P B A N     ",12,     &
                    (/"1/4,1/4,0      ", "3/4,1/4,0      ", "3/4,1/4,1/2    ",    &
                      "1/4,1/4,1/2    ", "0,0,0          ", "0,0,1/2        ",    &
                      "x,1/4,0        ", "x,1/4,1/2      ", "1/4,y,0        ",    &
                      "1/4,y,1/2      ", "1/4,1/4,z      ", "1/4,3/4,z      ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 65)= wyck_info_type("P M M A     ",11,     &
                    (/"0,0,0          ", "0,1/2,0        ", "0,0,1/2        ",    &
                      "0,1/2,1/2      ", "1/4,0,z        ", "1/4,1/2,z      ",    &
                      "0,y,0          ", "0,y,1/2        ", "x,0,z          ",    &
                      "x,1/2,z        ", "1/4,y,z        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 66)= wyck_info_type("P N N A     ", 4,     &
                    (/"0,0,0          ", "0,0,1/2        ", "1/4,0,z        ",    &
                      "x,1/4,1/4      ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 67)= wyck_info_type("P M N A     ", 8,     &
                    (/"0,0,0          ", "1/2,0,0        ", "1/2,1/2,0      ",    &
                      "0,1/2,0        ", "x,0,0          ", "x,1/2,0        ",    &
                      "1/4,y,1/4      ", "0,y,z          ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 68)= wyck_info_type("P C C A     ", 5,     &
                    (/"0,0,0          ", "0,1/2,0        ", "0,y,1/4        ",    &
                      "1/4,0,z        ", "1/4,1/2,z      ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 69)= wyck_info_type("P B A M     ", 8,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,1/2,0        ",    &
                      "0,1/2,1/2      ", "0,0,z          ", "0,1/2,z        ",    &
                      "x,y,0          ", "x,y,1/2        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 70)= wyck_info_type("P C C N     ", 4,     &
                    (/"0,0,0          ", "0,0,1/2        ", "1/4,1/4,z      ",    &
                      "1/4,3/4,z      ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 71)= wyck_info_type("P B C M     ", 4,     &
                    (/"0,0,0          ", "1/2,0,0        ", "x,1/4,0        ",    &
                      "x,y,1/4        ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 72)= wyck_info_type("P N N M     ", 7,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,1/2,0        ",    &
                      "0,1/2,1/2      ", "0,0,z          ", "0,1/2,z        ",    &
                      "x,y,0          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 73)= wyck_info_type("P M M N:1   ", 6,     &
                    (/"0,0,z          ", "0,1/2,z        ", "1/4,1/4,0      ",    &
                      "1/4,1/4,1/2    ", "0,y,z          ", "x,0,z          ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 74)= wyck_info_type("P M M N     ", 6,     &
                    (/"1/4,1/4,z      ", "1/4,3/4,z      ", "0,0,0          ",    &
                      "0,0,1/2        ", "1/4,y,z        ", "x,1/4,z        ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 75)= wyck_info_type("P B C N     ", 3,     &
                    (/"0,0,0          ", "0,1/2,0        ", "0,y,1/4        ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 76)= wyck_info_type("P B C A     ", 2,     &
                    (/"0,0,0          ", "0,0,1/2        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 77)= wyck_info_type("P N M A     ", 3,     &
                    (/"0,0,0          ", "0,0,1/2        ", "x,1/4,z        ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 78)= wyck_info_type("C M C M     ", 7,     &
                    (/"0,0,0          ", "0,1/2,0        ", "0,y,1/4        ",    &
                      "1/4,1/4,0      ", "x,0,0          ", "0,y,z          ",    &
                      "x,y,1/4        ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 79)= wyck_info_type("C M C A     ", 6,     &
                    (/"0,0,0          ", "1/2,0,0        ", "1/4,1/4,0      ",    &
                      "x,0,0          ", "1/4,y,1/4      ", "0,y,z          ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 80)= wyck_info_type("C M M M     ",17,     &
                    (/"0,0,0          ", "1/2,0,0        ", "1/2,0,1/2      ",    &
                      "0,0,1/2        ", "1/4,1/4,0      ", "1/4,1/4,1/2    ",    &
                      "x,0,0          ", "x,0,1/2        ", "0,y,0          ",    &
                      "0,y,1/2        ", "0,0,z          ", "0,1/2,z        ",    &
                      "1/4,1/4,z      ", "0,y,z          ", "x,0,z          ",    &
                      "x,y,0          ", "x,y,1/2        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 81)= wyck_info_type("C C C M     ",12,     &
                    (/"0,0,1/4        ", "0,1/2,1/4      ", "0,0,0          ",    &
                      "0,1/2,0        ", "1/4,1/4,0      ", "1/4,3/4,0      ",    &
                      "x,0,1/4        ", "0,y,1/4        ", "0,0,z          ",    &
                      "0,1/2,z        ", "1/4,1/4,z      ", "x,y,0          ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 82)= wyck_info_type("C M M A     ",14,     &
                    (/"1/4,0,0        ", "1/4,0,1/2      ", "0,0,0          ",    &
                      "0,0,1/2        ", "1/4,1/4,0      ", "1/4,1/4,1/2    ",    &
                      "0,1/4,z        ", "x,0,0          ", "x,0,1/2        ",    &
                      "1/4,y,0        ", "1/4,y,1/2      ", "1/4,0,z        ",    &
                      "0,y,z          ", "x,1/4,z        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 83)= wyck_info_type("C C C A:1   ", 8,     &
                    (/"0,0,0          ", "0,0,1/2        ", "1/4,0,1/4      ",    &
                      "0,1/4,1/4      ", "x,0,0          ", "0,y,0          ",    &
                      "0,0,z          ", "1/4,1/4,z      ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 84)= wyck_info_type("C C C A     ", 8,     &
                    (/"0,1/4,1/4      ", "0,1/4,3/2      ", "1/4,3/4,0      ",    &
                      "0,0,0          ", "x,1/4,1/4      ", "0,y,1/4        ",    &
                      "0,1/4,z        ", "1/4,0,z        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 85)= wyck_info_type("F M M M     ",15,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,1/4,1/4      ",    &
                      "1/4,0,1/4      ", "1/4,1/4,0      ", "1/4,1/4,1/4    ",    &
                      "x,0,0          ", "0,y,0          ", "0,0,z          ",    &
                      "1/4,1/4,z      ", "1/4,y,1/4      ", "x,1/4,1/4      ",    &
                      "0,y,z          ", "x,0,z          ", "x,y,0          ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 86)= wyck_info_type("F D D D:1   ", 7,     &
                    (/"0,0,0          ", "0,0,1/2        ", "1/8,1/8,1/8    ",    &
                      "5/8,5/8,5/8    ", "x,0,0          ", "0,y,0          ",    &
                      "0,0,z          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 87)= wyck_info_type("F D D D     ", 7,     &
                    (/"1/8,1/8,1/8    ", "1/8,1/8,5/8    ", "0,0,0          ",    &
                      "1/2,1/2,1/2    ", "x,1/8,1/8      ", "1/8,y,1/8      ",    &
                      "1/8,1/8,z      ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 88)= wyck_info_type("I M M M     ",14,     &
                    (/"0,0,0          ", "0,1/2,1/2      ", "1/2,1/2,0      ",    &
                      "1/2,0,1/2      ", "x,0,0          ", "x,1/2,0        ",    &
                      "0,y,0          ", "0,y,1/2        ", "0,0,z          ",    &
                      "1/2,0,z        ", "1/4,1/4,1/4    ", "0,y,z          ",    &
                      "x,0,z          ", "x,y,0          ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 89)= wyck_info_type("I B A M     ",10,     &
                    (/"0,0,1/4        ", "1/2,0,1/4      ", "0,0,0          ",    &
                      "1/2,0,0        ", "1/4,1/4,1/4    ", "x,0,1/4        ",    &
                      "0,y,1/4        ", "0,0,z          ", "0,1/2,z        ",    &
                      "x,y,0          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 90)= wyck_info_type("I B C A     ", 5,     &
                    (/"0,0,0          ", "1/4,1/4,1/4    ", "x,0,1/4        ",    &
                      "1/4,y,0        ", "0,1/4,z        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 91)= wyck_info_type("I M M A     ", 9,     &
                    (/"0,0,0          ", "0,0,1/2        ", "1/4,1/4,1/4    ",    &
                      "1/4,1/4,3/4    ", "0,1/4,z        ", "x,0,0          ",    &
                      "1/4,y,1/4      ", "0,y,z          ", "x,1/4,z        ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 92)= wyck_info_type("P 4         ", 3,     &
                    (/"0,0,z          ", "1/2,1/2,z      ", "0,1/2,z        ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 93)= wyck_info_type("P 41        ", 0,     &
                    (/"               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 94)= wyck_info_type("P 42        ", 3,     &
                    (/"0,0,z          ", "1/2,1/2,z      ", "0,1/2,z        ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 95)= wyck_info_type("P 43        ", 0,     &
                    (/"               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 96)= wyck_info_type("I 4         ", 2,     &
                    (/"0,0,z          ", "0,1/2,z        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 97)= wyck_info_type("I 41        ", 1,     &
                    (/"0,0,z          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 98)= wyck_info_type("P -4        ", 7,     &
                    (/"0,0,0          ", "0,0,1/2        ", "1/2,1/2,0      ",    &
                      "1/2,1/2,1/2    ", "0,0,z          ", "1/2,1/2,z      ",    &
                      "0,1/2,z        ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info( 99)= wyck_info_type("I -4        ", 6,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,1/2,1/4      ",    &
                      "0,1/2,3/4      ", "0,0,z          ", "0,1/2,z        ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(100)= wyck_info_type("P 4/M       ",11,     &
                    (/"0,0,0          ", "0,0,1/2        ", "1/2,1/2,0      ",    &
                      "1/2,1/2,1/2    ", "0,1/2,0        ", "0,1/2,1/2      ",    &
                      "0,0,z          ", "1/2,1/2,z      ", "0,1/2,z        ",    &
                      "x,y,0          ", "x,y,1/2        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(101)= wyck_info_type("P 42/M      ",10,     &
                    (/"0,0,0          ", "1/2,1/2,0      ", "0,1/2,0        ",    &
                      "0,1/2,1/2      ", "0,0,1/4        ", "1/2,1/2,1/4    ",    &
                      "0,0,z          ", "1/2,1/2,z      ", "0,1/2,z        ",    &
                      "x,y,0          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(102)= wyck_info_type("P 4/N:1     ", 6,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,1/2,z        ",    &
                      "1/4,1/4,0      ", "1/4,1/4,1/2    ", "0,0,z          ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(103)= wyck_info_type("P 4/N       ", 6,     &
                    (/"1/4,3/4,0      ", "1/4,3/4,1/2    ", "1/4,1/4,z      ",    &
                      "0,0,0          ", "0,0,1/2        ", "1/4,3/4,z      ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(104)= wyck_info_type("P 42/N:1    ", 6,     &
                    (/"0,0,0          ", "0,0,1/2        ", "1/4,1/4,1/4    ",    &
                      "1/4,1/4,3/4    ", "0,1/2,z        ", "0,0,z          ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(105)= wyck_info_type("P 42/N      ", 6,     &
                    (/"1/4,1/4,1/4    ", "1/4,1/4,3/4    ", "0,0,0          ",    &
                      "0,0,1/2        ", "3/4,1/4,z      ", "1/4,1/4,z      ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(106)= wyck_info_type("I 4/M       ", 8,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,1/2,0        ",    &
                      "0,1/2,1/4      ", "0,0,z          ", "1/4,1/4,1/4    ",    &
                      "0,1/2,z        ", "x,y,0          ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(107)= wyck_info_type("I 41/A:1    ", 5,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,1/4,1/8      ",    &
                      "0,1/4,5/8      ", "0,0,z          ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(108)= wyck_info_type("I 41/A      ", 5,     &
                    (/"0,1/4,1/8      ", "0,1/4,5/8      ", "0,0,0          ",    &
                      "0,0,1/2        ", "0,1/4,z        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(109)= wyck_info_type("P 4 2 2     ",15,     &
                    (/"0,0,0          ", "0,0,1/2        ", "1/2,1/2,0      ",    &
                      "1/2,1/2,1/2    ", "1/2,0,0        ", "1/2,0,1/2      ",    &
                      "0,0,z          ", "1/2,1/2,z      ", "0,1/2,z        ",    &
                      "x,x,0          ", "x,x,1/2        ", "x,0,0          ",    &
                      "x,1/2,1/2      ", "x,0,1/2        ", "x,1/2,0        ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(110)= wyck_info_type("P 4 21 2    ", 6,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,1/2,z        ",    &
                      "0,0,z          ", "x,x,0          ", "x,x,1/2        ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(111)= wyck_info_type("P 41 2 2    ", 3,     &
                    (/"0,y,0          ", "1/2,y,0        ", "x,x,3/8        ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(112)= wyck_info_type("P 41 21 2   ", 1,     &
                    (/"x,x,0          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(113)= wyck_info_type("P 42 2 2    ",15,     &
                    (/"0,0,0          ", "1/2,1/2,0      ", "0,1/2,0        ",    &
                      "0,1/2,1/2      ", "0,0,1/4        ", "1/2,1/2,1/4    ",    &
                      "0,0,z          ", "1/2,1/2,z      ", "0,1/2,z        ",    &
                      "x,0,0          ", "x,1/2,1/2      ", "x,0,1/2        ",    &
                      "x,1/2,0        ", "x,x,1/4        ", "x,x,3/4        ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(114)= wyck_info_type("P 42 21 2   ", 6,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,0,z          ",    &
                      "0,1/2,z        ", "x,x,0          ", "x,x,1/2        ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(115)= wyck_info_type("P 43 2 2    ", 3,     &
                    (/"0,y,0          ", "1/2,y,0        ", "x,x,5/8        ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(116)= wyck_info_type("P 43 21 2   ", 1,     &
                    (/"x,x,0          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(117)= wyck_info_type("I 4 2 2     ",10,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,1/2,0        ",    &
                      "0,1/2,1/4      ", "0,0,z          ", "0,1/2,z        ",    &
                      "x,x,0          ", "x,0,0          ", "x,0,1/2        ",    &
                      "x,x+1/2,1/4    ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(118)= wyck_info_type("I 41 2 2    ", 6,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,0,z          ",    &
                      "x,x,0          ", "-x,x,0         ", "x,1/4,1/8      ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )

       wyckoff_info(119)= wyck_info_type("P 4 M M     ", 6,     &
                    (/"0,0,z          ", "1/2,1/2,z      ", "1/2,0,z        ",    &
                      "x,x,z          ", "x,0,z          ", "x,1/2,z        ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(120)= wyck_info_type("P 4 B M     ", 3,     &
                    (/"0,0,z          ", "0,1/2,z        ", "x,x+1/2,z      ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(121)= wyck_info_type("P 42 C M    ", 4,     &
                    (/"0,0,z          ", "1/2,1/2,z      ", "0,1/2,z        ",    &
                      "x,x,z          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(122)= wyck_info_type("P 42 N M    ", 3,     &
                    (/"0,0,z          ", "0,1/2,z        ", "x,x,z          ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(123)= wyck_info_type("P 4 C C     ", 3,     &
                    (/"0,0,z          ", "1/2,1/2,z      ", "0,1/2,z        ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(124)= wyck_info_type("P 4 N C     ", 2,     &
                    (/"0,0,z          ", "0,1/2,z        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(125)= wyck_info_type("P 42 M C    ", 5,     &
                    (/"0,0,z          ", "1/2,1/2,z      ", "0,1/2,z        ",    &
                      "x,0,z          ", "x,1/2,z        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(126)= wyck_info_type("P 42 B C    ", 2,     &
                    (/"0,0,z          ", "0,1/2,z        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(127)= wyck_info_type("I 4 M M     ", 4,     &
                    (/"0,0,z          ", "0,1/2,z        ", "x,x,z          ",    &
                      "x,0,z          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(128)= wyck_info_type("I 4 C M     ", 3,     &
                    (/"0,0,z          ", "1/2,0,z        ", "x,x+1/2,z      ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(129)= wyck_info_type("I 41 M D    ", 2,     &
                    (/"0,0,z          ", "0,y,z          ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(130)= wyck_info_type("I 41 C D    ", 1,     &
                    (/"0,0,z          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(131)= wyck_info_type("P -4 2 M    ",14,     &
                    (/"0,0,0          ", "1/2,1/2,1/2    ", "0,0,1/2        ",    &
                      "1/2,1/2,0      ", "1/2,0,0        ", "1/2,0,1/2      ",    &
                      "0,0,z          ", "1/2,1/2,z      ", "x,0,0          ",    &
                      "x,1/2,1/2      ", "x,0,1/2        ", "x,1/2,0        ",    &
                      "0,1/2,z        ", "x,x,z          ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(132)= wyck_info_type("P -4 2 C    ",13,     &
                    (/"0,0,1/4        ", "1/2,0,1/4      ", "1/2,1/2,1/4    ",    &
                      "0,1/2,1/4      ", "0,0,0          ", "1/2,1/2,0      ",    &
                      "x,0,1/4        ", "1/2,y,1/4      ", "x,1/2,1/4      ",    &
                      "0,y,1/4        ", "0,0,z          ", "1/2,1/2,z      ",    &
                      "0,1/2,z        ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(133)= wyck_info_type("P -4 21 M   ", 5,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,1/2,z        ",    &
                      "0,0,z          ", "x,x+1/2,z      ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(134)= wyck_info_type("P -4 21 C   ", 4,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,0,z          ",    &
                      "0,1/2,z        ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(135)= wyck_info_type("P -4 M 2    ",11,     &
                    (/"0,0,0          ", "1/2,1/2,0      ", "1/2,1/2,1/2    ",    &
                      "0,0,1/2        ", "0,0,z          ", "1/2,1/2,z      ",    &
                      "0,1/2,z        ", "x,x,0          ", "x,x,1/2        ",    &
                      "x,0,z          ", "x,1/2,z        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(136)= wyck_info_type("P -4 C 2    ", 9,     &
                    (/"0,0,1/4        ", "1/2,1/2,1/4    ", "0,0,0          ",    &
                      "1/2,1/2,0      ", "x,x,1/4        ", "x,x,3/4        ",    &
                      "0,0,z          ", "1/2,1/2,z      ", "0,1/2,z        ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(137)= wyck_info_type("P -4 B 2    ", 8,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,1/2,0        ",    &
                      "0,1/2,1/2      ", "0,0,z          ", "0,1/2,z        ",    &
                      "x,x+1/2,0      ", "x,x+1/2,1/2    ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(138)= wyck_info_type("P -4 N 2    ", 8,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,1/2,1/4      ",    &
                      "0,1/2,3/4      ", "0,0,z          ", "x,-x+1/2,1/4   ",    &
                      "x,x+1/2,1/4    ", "0,1/2,z        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(139)= wyck_info_type("I -4 M 2    ", 9,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,1/2,1/4      ",    &
                      "0,1/2,3/4      ", "0,0,z          ", "0,1/2,z        ",    &
                      "x,x,0          ", "x,x+1/2,1/4    ", "x,0,z          ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(140)= wyck_info_type("I -4 C 2    ", 8,     &
                    (/"0,0,1/4        ", "0,0,0          ", "0,1/2,1/4      ",    &
                      "0,1/2,0        ", "x,x,1/4        ", "0,0,z          ",    &
                      "0,1/2,z        ", "x,x+1/2,0      ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(141)= wyck_info_type("I -4 2 M    ", 9,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,1/2,0        ",    &
                      "0,1/2,1/4      ", "0,0,z          ", "x,0,0          ",    &
                      "x,0,1/2        ", "0,1/2,z        ", "x,x,z          ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(142)= wyck_info_type("I -4 2 D    ", 4,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,0,z          ",    &
                      "x,1/4,1/8      ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(143)= wyck_info_type("P 4/M M M   ",20,     &
                    (/"0,0,0          ", "0,0,1/2        ", "1/2,1/2,0      ",    &
                      "1/2,1/2,1/2    ", "0,1/2,1/2      ", "0,1/2,0        ",    &
                      "0,0,z          ", "1/2,1/2,z      ", "0,1/2,z        ",    &
                      "x,x,0          ", "x,x,1/2        ", "x,0,0          ",    &
                      "x,0,1/2        ", "x,1/2,0        ", "x,1/2,1/2      ",    &
                      "x,y,0          ", "x,y,1/2        ", "x,x,z          ",    &
                      "x,0,z          ", "x,1/2,z        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(144)= wyck_info_type("P 4/M C C   ",13,     &
                    (/"0,0,1/4        ", "0,0,0          ", "1/2,1/2,1/4    ",    &
                      "1/2,1/2,0      ", "0,1/2,0        ", "0,1/2,1/4      ",    &
                      "0,0,z          ", "1/2,1/2,z      ", "0,1/2,z        ",    &
                      "x,x,1/4        ", "x,0,1/4        ", "x,1/2,1/4      ",    &
                      "x,y,0          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(145)= wyck_info_type("P 4/N B M:1 ",13,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,1/2,0        ",    &
                      "0,1/2,1/2      ", "1/4,1/4,0      ", "1/4,1/4,1/2    ",    &
                      "0,0,z          ", "0,1/2,z        ", "x,x,0          ",    &
                      "x,x,1/2        ", "x,0,0          ", "x,0,1/2        ",    &
                      "x,x+1/2,z      ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(146)= wyck_info_type("P 4/N B M   ",13,     &
                    (/"1/4,1/4,0      ", "1/4,1/4,1/2    ", "3/4,1/4,0      ",    &
                      "3/4,1/4,1/2    ", "0,0,0          ", "0,0,1/2        ",    &
                      "1/4,1/4,z      ", "3/4,1/4,z      ", "x,x,0          ",    &
                      "x,x,1/2        ", "x,1/4,0        ", "x,1/4,1/2      ",    &
                      "x,-x,z         ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(147)= wyck_info_type("P 4/N N C:1 ",10,     &
                    (/"0,0,0          ", "0,0,1/2        ", "1/2,0,0        ",    &
                      "1/2,0,1/4      ", "0,0,z          ", "1/4,1/4,1/4    ",    &
                      "1/2,0,z        ", "x,x,0          ", "x,0,0          ",    &
                      "x,0,1/2        ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(148)= wyck_info_type("P 4/N N C   ",10,     &
                    (/"1/4,1/4,1/4    ", "1/4,1/4,3/4    ", "1/4,3/4,3/4    ",    &
                      "1/4,1/4,0      ", "1/4,1/4,z      ", "0,0,0          ",    &
                      "1/4,3/4,z      ", "x,x,1/4        ", "x,1/4,1/4      ",    &
                      "x,3/4,1/4      ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(149)= wyck_info_type("P 4/M B M   ",11,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,1/2,1/2      ",    &
                      "0,1/2,0        ", "0,0,z          ", "0,1/2,z        ",    &
                      "x,x+1/2,0      ", "x,x+1/2,1/2    ", "x,y,0          ",    &
                      "x,y,1/2        ", "x,x+1/2,z      ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(150)= wyck_info_type("P 4/M N C   ", 8,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,1/2,0        ",    &
                      "0,1/2,1/4      ", "0,0,z          ", "0,1/2,z        ",    &
                      "x,x+1/2,1/4    ", "x,y,0          ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(151)= wyck_info_type("P 4/N M M:1 ",10,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,1/2,z        ",    &
                      "1/4,1/4,0      ", "1/4,1/4,1/2    ", "0,0,z          ",    &
                      "x,x,0          ", "x,x,1/2        ", "0,y,z          ",    &
                      "x,x+1/2,z      ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(152)= wyck_info_type("P 4/N M M   ",10,     &
                    (/"3/4,1/4,0      ", "3/4,1/4,1/2    ", "1/4,1/4,z      ",    &
                      "0,0,0          ", "0,0,1/2        ", "3/4,1/4,z      ",    &
                      "x,-x,0         ", "x,-x,1/2       ", "1/4,y,z        ",    &
                      "x,x,z          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(153)= wyck_info_type("P 4/N C C:1 ", 6,     &
                    (/"0,0,1/4        ", "0,0,0          ", "0,1/2,z        ",    &
                      "1/4,1/4,0      ", "0,0,z          ", "x,x,1/4        ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(154)= wyck_info_type("P 4/N C C   ", 6,     &
                    (/"3/4,1/4,1/4    ", "3/4,1/4,0      ", "1/4,1/4,z      ",    &
                      "0,0,0          ", "3/4,1/4,z      ", "x,-x,1/4       ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(155)= wyck_info_type("P 42/M M C  ",17,     &
                    (/"0,0,0          ", "1/2,1/2,0      ", "0,1/2,0        ",    &
                      "0,1/2,1/2      ", "0,0,1/4        ", "1/2,1/2,1/4    ",    &
                      "0,0,z          ", "1/2,1/2,z      ", "0,1/2,z        ",    &
                      "x,0,0          ", "x,1/2,1/2      ", "x,0,1/2        ",    &
                      "x,1/2,0        ", "x,x,1/4        ", "0,y,z          ",    &
                      "1/2,y,z        ", "x,y,0          ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(156)= wyck_info_type("P 42/M C M  ",15,     &
                    (/"0,0,0          ", "0,0,1/4        ", "1/2,1/2,0      ",    &
                      "1/2,1/2,1/4    ", "0,1/2,1/4      ", "0,1/2,0        ",    &
                      "0,0,z          ", "1/2,1/2,z      ", "x,x,0          ",    &
                      "x,x,1/2        ", "0,1/2,z        ", "x,0,1/4        ",    &
                      "x,1/2,1/4      ", "x,y,0          ", "x,x,z          ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(157)= wyck_info_type("P 42/N B C:1",10,     &
                    (/"0,1/2,1/4      ", "0,0,1/4        ", "0,1/2,0        ",    &
                      "0,0,0          ", "1/4,1/4,1/4    ", "0,1/2,z        ",    &
                      "0,0,z          ", "x,0,1/4        ", "x,0,3/4        ",    &
                      "x,x+1/2,0      ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(158)= wyck_info_type("P 42/N B C  ",10,     &
                    (/"1/4,1/4,0      ", "3/4,1/4,0      ", "1/4,1/4,1/4    ",    &
                      "3/4,1/4,3/4    ", "0,0,0          ", "1/4,1/4,z      ",    &
                      "3/4,1/4,z      ", "x,1/4,0        ", "x,1/4,1/2      ",    &
                      "x,x,1/4        ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(159)= wyck_info_type("P 42/N N M:1",13,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,1/2,0        ",    &
                      "0,1/2,1/4      ", "1/4,1/4,1/4    ", "3/4,3/4,3/4    ",    &
                      "0,0,z          ", "0,1/2,z        ", "x,0,0          ",    &
                      "x,0,1/2        ", "x,x+1/2,1/4    ", "x,x+1/2,3/4    ",    &
                      "x,x,z          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(160)= wyck_info_type("P 42/N N M  ",13,     &
                    (/"1/4,3/4,1/4    ", "3/4,1/4,1/4    ", "1/4,1/4,1/4    ",    &
                      "1/4,1/4,0      ", "0,0,1/2        ", "0,0,0          ",    &
                      "3/4,1/4,z      ", "1/4,1/4,z      ", "x,1/4,3/4      ",    &
                      "x,1/4,1/4      ", "x,x,0          ", "x,x,1/2        ",    &
                      "x,-x,z         ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(161)= wyck_info_type("P 42/M B C  ", 8,     &
                    (/"0,0,0          ", "0,0,1/4        ", "0,1/2,0        ",    &
                      "0,1/2,1/4      ", "0,0,z          ", "0,1/2,z        ",    &
                      "x,x+1/2,1/4    ", "x,y,0          ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(162)= wyck_info_type("P 42/M N M  ",10,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,1/2,0        ",    &
                      "0,1/2,1/4      ", "0,0,z          ", "x,x,0          ",    &
                      "x,-x,0         ", "0,1/2,z        ", "x,y,0          ",    &
                      "x,x,z          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(163)= wyck_info_type("P 42/N M C:1", 7,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,0,z          ",    &
                      "0,1/2,z        ", "1/4,1/4,1/4    ", "x,x,0          ",    &
                      "0,y,z          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(164)= wyck_info_type("P 42/N M C  ", 7,     &
                    (/"3/4,1/4,3/4    ", "3/4,1/4,1/4    ", "3/4,1/4,z      ",    &
                      "1/4,1/4,z      ", "0,0,0          ", "x,-x,1/4       ",    &
                      "1/4,y,z        ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(165)= wyck_info_type("P 42/N C M:1", 9,     &
                    (/"0,0,1/4        ", "0,0,0          ", "1/4,1/4,1/4    ",    &
                      "1/4,1/4,3/4    ", "0,1/2,z        ", "0,0,z          ",    &
                      "x,x,1/4        ", "x,x,3/4        ", "x,x+1/2,z      ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(166)= wyck_info_type("P 42/N C M  ", 9,     &
                    (/"3/4,1/4,0      ", "3/4,1/4,3/4    ", "0,0,1/2        ",    &
                      "0,0,0          ", "1/4,1/4,z      ", "3/4,1/4,z      ",    &
                      "x,-x,1/2       ", "x,-x,0         ", "x,x,z          ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(167)= wyck_info_type("I 4/M M M   ",14,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,1/2,0        ",    &
                      "0,1/2,1/4      ", "0,0,z          ", "1/4,1/4,1/4    ",    &
                      "0,1/2,z        ", "x,x,0          ", "x,0,0          ",    &
                      "x,1/2,0        ", "x,x+1/2,1/4    ", "x,y,0          ",    &
                      "x,x,z          ", "0,y,z          ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(168)= wyck_info_type("I 4/M C M   ",12,     &
                    (/"0,0,1/4        ", "0,1/2,1/4      ", "0,0,0          ",    &
                      "0,1/2,0        ", "1/4,1/4,1/4    ", "0,0,z          ",    &
                      "0,1/2,z        ", "x,x+1/2,0      ", "x,x,1/4        ",    &
                      "x,0,1/4        ", "x,y,0          ", "x,x+1/2,z      ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(169)= wyck_info_type("I 41/A M D:1", 8,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,1/4,1/8      ",    &
                      "0,1/4,5/8      ", "0,0,z          ", "x,1/4,1/8      ",    &
                      "x,x,0          ", "0,y,z          ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(170)= wyck_info_type("I 41/A M D  ", 8,     &
                    (/"0,3/4,1/8      ", "0,1/4,3/8      ", "0,0,0          ",    &
                      "0,0,1/2        ", "0,1/4,z        ", "x,0,0          ",    &
                      "x,x+1/4,7/8    ", "0,y,z          ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(171)= wyck_info_type("I 41/A C D:1", 6,     &
                    (/"0,0,0          ", "0,0,1/4        ", "0,1/4,1/8      ",    &
                      "0,0,z          ", "1/4,y,1/8      ", "x,x,1/4        ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(172)= wyck_info_type("I 41/A C D  ", 6,     &
                    (/"0,1/4,3/8      ", "0,1/4,1/8      ", "0,0,0          ",    &
                      "0,1/4,z        ", "x,0,1/4        ", "x,x+1/4,1/8    ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(173)= wyck_info_type("P 3         ", 3,     &
                    (/"0,0,z          ", "1/3,2/3,z      ", "2/3,1/3,z      ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(174)= wyck_info_type("P 31        ", 0,     &
                    (/"               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(175)= wyck_info_type("P 32        ", 0,     &
                    (/"               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(176)= wyck_info_type("R 3         ", 1,     &
                    (/"0,0,z          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(177)= wyck_info_type("R 3:H       ", 1,     &
                    (/"x,x,x          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(178)= wyck_info_type("P -3        ", 6,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,0,z          ",    &
                      "1/3,2/3,z      ", "1/2,0,0        ", "1/2,0,1/2      ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(179)= wyck_info_type("R -3        ", 5,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,0,z          ",    &
                      "1/2,0,1/2      ", "1/2,0,0        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(180)= wyck_info_type("R -3:H      ", 5,     &
                    (/"0,0,0          ", "1/2,1/2,1/2    ", "x,x,x          ",    &
                      "1/2,0,0        ", "0,1/2,1/2      ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(181)= wyck_info_type("P 3 1 2     ",11,     &
                    (/"0,0,0          ", "0,0,1/2        ", "1/3,2/3,0      ",    &
                      "1/3,2/3,1/2    ", "2/3,1/3,0      ", "2/3,1/3,1/2    ",    &
                      "0,0,z          ", "1/3,2/3,z      ", "2/3,1/3,z      ",    &
                      "x,-x,0         ", "x,-x,1/2       ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(182)= wyck_info_type("P 3 2 1     ", 6,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,0,z          ",    &
                      "1/3,2/3,z      ", "x,0,0          ", "x,0,1/2        ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(183)= wyck_info_type("P 31 1 2    ", 2,     &
                    (/"x,-x,1/3       ", "x,-x,5/6       ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(184)= wyck_info_type("P 31 2 1    ", 2,     &
                    (/"x,0,1/3        ", "x,0,5/6        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(185)= wyck_info_type("P 32 1 2    ", 2,     &
                    (/"x,-x,2/3       ", "x,-x,1/6       ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(186)= wyck_info_type("P 32 2 1    ", 2,     &
                    (/"x,0,2/3        ", "x,0,1/6        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(187)= wyck_info_type("R 3 2       ", 5,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,0,z          ",    &
                      "x,0,0          ", "x,0,1/2        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(188)= wyck_info_type("R 3 2:R     ", 5,     &
                    (/"0,0,0          ", "1/2,1/2,1/2    ", "x,x,x          ",    &
                      "0,y,-y         ", "1/2,y,-y       ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(189)= wyck_info_type("P 3 M 1     ", 4,     &
                    (/"0,0,z          ", "1/3,2/3,z      ", "2/3,1/3,z      ",    &
                      "x,-x,z         ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(190)= wyck_info_type("P 3 1 M     ", 3,     &
                    (/"0,0,z          ", "1/3,2/3,z      ", "x,0,z          ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(191)= wyck_info_type("P 3 C 1     ", 3,     &
                    (/"0,0,z          ", "1/3,2/3,z      ", "2/3,1/3,z      ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(192)= wyck_info_type("P 3 1 C     ", 2,     &
                    (/"0,0,z          ", "1/3,2/3,z      ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(193)= wyck_info_type("R 3 M       ", 2,     &
                    (/"0,0,z          ", "x,-x,z         ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(194)= wyck_info_type("R 3 M:R     ", 2,     &
                    (/"x,x,x          ", "x,x,z          ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(195)= wyck_info_type("R 3 C       ", 1,     &
                    (/"0,0,z          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(196)= wyck_info_type("R 3 C:R     ", 1,     &
                    (/"x,x,x          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(197)= wyck_info_type("P -3 1 M    ",11,     &
                    (/"0,0,0          ", "0,0,1/2        ", "1/3,2/3,0      ",    &
                      "1/3,2/3,1/2    ", "0,0,z          ", "1/2,0,0        ",    &
                      "1/2,0,1/2      ", "1/3,2/3,z      ", "x,-x,0         ",    &
                      "x,-x,1/2       ", "x,0,z          ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(198)= wyck_info_type("P -3 1 C    ", 8,     &
                    (/"0,0,1/4        ", "0,0,0          ", "1/3,2/3,1/4    ",    &
                      "2/3,1/3,1/4    ", "0,0,z          ", "1/3,2/3,z      ",    &
                      "1/2,0,0        ", "x,-x,1/4       ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(199)= wyck_info_type("P -3 M 1    ", 9,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,0,z          ",    &
                      "1/3,2/3,z      ", "1/2,0,0        ", "1/2,0,1/2      ",    &
                      "x,0,0          ", "x,0,1/2        ", "x,-x,z         ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(200)= wyck_info_type("P -3 C 1    ", 6,     &
                    (/"0,0,1/4        ", "0,0,0          ", "0,0,z          ",    &
                      "1/3,2/3,z      ", "1/2,0,0        ", "x,0,1/4        ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(201)= wyck_info_type("R -3 M      ", 8,     &
                    (/"0,0,0          ", "0,0,1/2        ", "0,0,z          ",    &
                      "1/2,0,1/2      ", "1/2,0,0        ", "x,0,0          ",    &
                      "x,0,1/2        ", "x,-x,z         ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(202)= wyck_info_type("R -3 M:R    ", 8,     &
                    (/"0,0,0          ", "1/2,1/2,1/2    ", "x,x,x          ",    &
                      "1/2,0,0        ", "0,1/2,1/2      ", "x,-x,0         ",    &
                      "x,-x,1/2       ", "x,x,z          ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(203)= wyck_info_type("R -3 C      ", 5,     &
                    (/"0,0,1/4        ", "0,0,0          ", "0,0,z          ",    &
                      "1/2,0,0        ", "x,0,1/4        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(204)= wyck_info_type("R -3 C:R    ", 5,     &
                    (/"1/4,1/4,1/4    ", "0,0,0          ", "x,x,x          ",    &
                      "1/2,0,0        ", "x,-x+1/2,1/4   ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(205)= wyck_info_type("P 6         ", 3,     &
                    (/"0,0,z          ", "1/3,2/3,z      ", "1/2,0,z        ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(206)= wyck_info_type("P 61        ", 0,     &
                    (/"               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(207)= wyck_info_type("P 65        ", 0,     &
                    (/"               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(208)= wyck_info_type("P 62        ", 2,     &
                    (/"0,0,z          ", "1/2,1/2,z      ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(209)= wyck_info_type("P 64        ", 2,     &
                    (/"0,0,z          ", "1/2,1/2,z      ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(210)= wyck_info_type("P 63        ", 2,     &
                    (/"0,0,z          ", "1/3,2/3,z      ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(211)= wyck_info_type("P -6        ",11,     &
                    (/"0,0,0          ", "0,0,1/2        ", "1/3,2/3,0      ",    &
                      "1/3,2/3,1/2    ", "2/3,1/3,0      ", "2/3,1/3,1/2    ",    &
                      "0,0,z          ", "1/3,2/3,z      ", "2/3,1/3,z      ",    &
                      "x,y,0          ", "x,y,1/2        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(212)= wyck_info_type("P 6/M       ",11,     &
                    (/"0,0,0          ", "0,0,1/2        ", "1/3,2/3,0      ",    &
                      "1/3,2/3,1/2    ", "0,0,z          ", "1/2,0,0        ",    &
                      "1/2,0,1/2      ", "1/3,2/3,z      ", "1/2,0,z        ",    &
                      "x,y,0          ", "x,y,1/2        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(213)= wyck_info_type("P 63/M      ", 8,     &
                    (/"0,0,1/4        ", "0,0,0          ", "1/3,2/3,1/4    ",    &
                      "2/3,1/3,1/4    ", "0,0,z          ", "1/3,2/3,z      ",    &
                      "1/2,0,0        ", "x,y,1/4        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(214)= wyck_info_type("P 6 2 2     ",13,     &
                    (/"0,0,0          ", "0,0,1/2        ", "1/3,2/3,0      ",    &
                      "1/3,2/3,1/2    ", "0,0,z          ", "1/2,0,0        ",    &
                      "1/2,0,1/2      ", "1/3,2/3,z      ", "1/2,0,z        ",    &
                      "x,0,0          ", "x,0,1/2        ", "x,-x,0         ",    &
                      "x,-x,1/2       ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(215)= wyck_info_type("P 61 2 2    ", 2,     &
                    (/"x,0,0          ", "x,2x,1/4       ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(216)= wyck_info_type("P 65 2 2    ", 2,     &
                    (/"x,0,0          ", "x,2x,3/4       ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(217)= wyck_info_type("P 62 2 2    ",10,     &
                    (/"0,0,0          ", "0,0,1/2        ", "1/2,0,0        ",    &
                      "1/2,0,1/2      ", "0,0,z          ", "1/2,0,z        ",    &
                      "x,0,0          ", "x,0,1/2        ", "x,2x,0         ",    &
                      "x,2x,1/2       ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(218)= wyck_info_type("P 64 2 2    ",10,     &
                    (/"0,0,0          ", "0,0,1/2        ", "1/2,0,0        ",    &
                      "1/2,0,/1,2     ", "0,0,z          ", "1/2,0,z        ",    &
                      "x,0,0          ", "x,0,1/2        ", "x,2x,0         ",    &
                      "x,2x,1/2       ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(219)= wyck_info_type("P 63 2 2    ", 8,     &
                    (/"0,0,0          ", "0,0,1/4        ", "1/3,2/3,1/4    ",    &
                      "1/3,2/3,3/4    ", "0,0,z          ", "1/3,2/3,z      ",    &
                      "x,0,0          ", "x,2x,1/4       ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(220)= wyck_info_type("P 6 M M     ", 5,     &
                    (/"0,0,z          ", "1/3,2/3,z      ", "1/2,0,z        ",    &
                      "x,0,z          ", "x,-x,z         ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(221)= wyck_info_type("P 6 C C     ", 3,     &
                    (/"0,0,z          ", "1/3,2/3,z      ", "1/2,0,z        ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(222)= wyck_info_type("P 63 C M    ", 3,     &
                    (/"0,0,z          ", "1/3,2/3,z      ", "x,0,z          ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(223)= wyck_info_type("P 63 M C    ", 3,     &
                    (/"0,0,z          ", "1/3,2/3,z      ", "x,-x,z         ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(224)= wyck_info_type("P -6 M 2    ",14,     &
                    (/"0,0,0          ", "0,0,1/2        ", "1/3,2/3,0      ",    &
                      "1/3,2/3,1/2    ", "2/3,1/3,0      ", "2/3,1/3,1/2    ",    &
                      "0,0,z          ", "1/3,2/3,z      ", "2/3,1/3,z      ",    &
                      "x,-x,0         ", "x,-x,1/2       ", "x,y,0          ",    &
                      "x,y,1/2        ", "x,-x,z         ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(225)= wyck_info_type("P -6 C 2    ",11,     &
                    (/"0,0,0          ", "0,0,1/4        ", "1/3,2/3,0      ",    &
                      "1/3,2/3,1/4    ", "2/3,1/3,0      ", "2/3,1/3,1/4    ",    &
                      "0,0,z          ", "1/3,2/3,z      ", "2/3,1/3,z      ",    &
                      "x,-x,0         ", "x,y,1/4        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(226)= wyck_info_type("P -6 2 M    ",11,     &
                    (/"0,0,0          ", "0,0,1/2        ", "1/3,2/3,0      ",    &
                      "1/3,2/3,1/2    ", "0,0,z          ", "x,0,0          ",    &
                      "x,0,1/2        ", "1/3,2/3,z      ", "x,0,z          ",    &
                      "x,y,0          ", "x,y,1/2        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(227)= wyck_info_type("P -6 2 C    ", 8,     &
                    (/"0,0,0          ", "0,0,1/4        ", "1/3,2/3,1/4    ",    &
                      "2/3,1/3,1/4    ", "0,0,z          ", "1/3,2/3,z      ",    &
                      "x,0,0          ", "x,y,1/4        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(228)= wyck_info_type("P 6/M M M   ",17,     &
                    (/"0,0,0          ", "0,0,1/2        ", "1/3,2/3,0      ",    &
                      "1/3,2/3,1/2    ", "0,0,z          ", "1/2,0,0        ",    &
                      "1/2,0,1/2      ", "1/3,2/3,z      ", "1/2,0,z        ",    &
                      "x,0,0          ", "x,0,1/2        ", "x,2x,0         ",    &
                      "x,2x,1/2       ", "x,0,z          ", "x,2x,z         ",    &
                      "x,y,0          ", "x,y,1/2        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(229)= wyck_info_type("P 6/M C C   ",12,     &
                    (/"0,0,1/4        ", "0,0,0          ", "1/3,2/3,1/4    ",    &
                      "1/3,2/3,0      ", "0,0,z          ", "1/2,0,1/4      ",    &
                      "1/2,0,0        ", "1/3,2/3,z      ", "1/2,0,z        ",    &
                      "x,0,1/4        ", "x,2x,1/4       ", "x,y,0          ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(230)= wyck_info_type("P 63/M C M  ",11,     &
                    (/"0,0,1/4        ", "0,0,0          ", "1/3,2/3,1/4    ",    &
                      "1/3,2/3,0      ", "0,0,z          ", "1/2,0,0        ",    &
                      "x,0,1/4        ", "1/3,2/3,z      ", "x,2x,0         ",    &
                      "x,y,1/4        ", "x,0,z          ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(231)= wyck_info_type("P 63/M M C  ",11,     &
                    (/"0,0,0          ", "0,0,1/4        ", "1/3,2/3,1/4    ",    &
                      "1/3,2/3,3/4    ", "0,0,z          ", "1/3,2/3,z      ",    &
                      "1/2,0,0        ", "x,2x,1/4       ", "x,0,0          ",    &
                      "x,y,1/4        ", "x,2x,z         ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(232)= wyck_info_type("P 2 3       ", 9,     &
                    (/"0,0,0          ", "1/2,1/2,1/2    ", "0,1/2,1/2      ",    &
                      "1/2,0,0        ", "x,x,x          ", "x,0,0          ",    &
                      "x,0,1/2        ", "x,1/2,0        ", "x,1/2,1/2      ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(233)= wyck_info_type("F 2 3       ", 7,     &
                    (/"0,0,0          ", "1/2,1/2,1/2    ", "1/4,1/4,1/4    ",    &
                      "3/4,3/4,3/4    ", "x,x,x          ", "x,0,0          ",    &
                      "x,1/4,1/4      ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(234)= wyck_info_type("I 2 3       ", 5,     &
                    (/"0,0,0          ", "0,1/2,1/2      ", "x,x,x          ",    &
                      "x,0,0          ", "x,1/2,0        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(235)= wyck_info_type("P 21 3      ", 1,     &
                    (/"x,x,x          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(236)= wyck_info_type("I 21 3      ", 2,     &
                    (/"x,x,x          ", "x,0,1/4        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(237)= wyck_info_type("P M -3      ",11,     &
                    (/"0,0,0          ", "1/2,1/2,1/2    ", "0,1/2,1/2      ",    &
                      "1/2,0,0        ", "x,0,0          ", "x,0,1/2        ",    &
                      "x,1/2,0        ", "x,1/2,1/2      ", "x,x,x          ",    &
                      "0,y,z          ", "1/2,y,z        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(238)= wyck_info_type("P N -3:1    ", 7,     &
                    (/"0,0,0          ", "1/4,1/4,1/4    ", "3/4,3/4,3/4    ",    &
                      "0,1/2,1/2      ", "x,x,x          ", "x,0,0          ",    &
                      "x,1/2,0        ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(239)= wyck_info_type("P N -3      ", 7,     &
                    (/"1/4,1/4,1/4    ", "0,0,0          ", "1/2,1/2,1/2    ",    &
                      "1/4,3/4,3/4    ", "x,x,x          ", "x,1/4,1/4      ",    &
                      "x,3/4,1/4      ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(240)= wyck_info_type("F M -3      ", 8,     &
                    (/"0,0,0          ", "1/2,1/2,1/2    ", "1/4,1/4,1/4    ",    &
                      "0,1/4,1/4      ", "x,0,0          ", "x,x,x          ",    &
                      "x,1/4,1/4      ", "0,y,z          ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(241)= wyck_info_type("F D -3:1    ", 6,     &
                    (/"0,0,0          ", "1/2,1/2,1/2    ", "1/8,1/8,1/8    ",    &
                      "5/8,5/8,5/8    ", "x,x,x          ", "x,0,0          ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(242)= wyck_info_type("F D -3      ", 6,     &
                    (/"1/8,1/8,1/8    ", "5/8,5/8,5/8    ", "0,0,0          ",    &
                      "1/2,1/2,1/2    ", "x,x,x          ", "x,1/8,1/8      ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(243)= wyck_info_type("I M -3      ", 7,     &
                    (/"0,0,0          ", "0,1/2,1/2      ", "1/4,1/4,1/4    ",    &
                      "x,0,0          ", "x,0,1/2        ", "x,x,x          ",    &
                      "0,y,z          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(244)= wyck_info_type("P A -3      ", 3,     &
                    (/"0,0,0          ", "1/2,1/2,1/2    ", "x,x,x          ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(245)= wyck_info_type("I A -3      ", 4,     &
                    (/"0,0,0          ", "1/4,1/4,1/4    ", "x,x,x          ",    &
                      "x,0,1/4        ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(246)= wyck_info_type("P 4 3 2     ",10,     &
                    (/"0,0,0          ", "1/2,1/2,1/2    ", "0,1/2,1/2      ",    &
                      "1/2,0,0        ", "x,0,0          ", "x,1/2,1/2      ",    &
                      "x,x,x          ", "x,1/2,0        ", "0,y,y          ",    &
                      "1/2,y,y        ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(247)= wyck_info_type("P 42 3 2    ",12,     &
                    (/"0,0,0          ", "1/4,1/4,1/4    ", "3/4,3/4,3/4    ",    &
                      "0,1/2,1/2      ", "1/4,0,1/2      ", "1/4,1/2,0      ",    &
                      "x,x,x          ", "x,0,0          ", "x,0,1/2        ",    &
                      "x,1/2,0        ", "1/4,y,-y+1/2   ", "1/4,y,y+1/2    ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(248)= wyck_info_type("F 4 3 2    ", 9,     &
                    (/"0,0,0          ", "1/2,1/2,1/2    ", "1/4,1/4,1/4    ",    &
                      "0,1/4,1/4      ", "x,0,0          ", "x,x,x          ",    &
                      "0,y,y          ", "1/2,y,y        ", "x,1/4,1/4      ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(249)= wyck_info_type("F 41 3 2    ", 7,     &
                    (/"0,0,0          ", "1/2,1/2,1/2    ", "1/8,1/8,1/8    ",    &
                      "5/8,5/8,5/8    ", "x,x,x          ", "x,0,0          ",    &
                      "1/8,y,-y+1/4   ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(250)= wyck_info_type("I 4 3 2     ", 9,     &
                    (/"0,0,0          ", "0,1/2,1/2      ", "1/4,1/4,1/4    ",    &
                      "1/4,1/2,0      ", "x,0,0          ", "x,x,x          ",    &
                      "x,1/2,0        ", "0,y,y          ", "1/4,y,-y+1/2   ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(251)= wyck_info_type("P 43 3 2    ", 4,     &
                    (/"1/8,1/8,1/8    ", "5/8,5/8,5/8    ", "x,x,x          ",    &
                      "1/8,y,-y+1/4   ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(252)= wyck_info_type("P 41 3 2    ", 4,     &
                    (/"3/8,3/8,3/8    ", "7/8,7/8,7/8    ", "x,x,x          ",    &
                      "1/8,y,y+1/4    ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(253)= wyck_info_type("I 41 3 2    ", 8,     &
                    (/"1/8,1/8,1/8    ", "7/8,7/8,7/8    ", "1/8,0,1/4      ",    &
                      "5/8,0,1/4      ", "x,x,x          ", "x,0,1/4        ",    &
                      "1/8,y,y+1/4    ", "1/8,y,-y+1/4   ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(254)= wyck_info_type("P -4 3 M    ", 9,     &
                    (/"0,0,0          ", "1/2,1/2,1/2    ", "0,1/2,1/2      ",    &
                      "1/2,0,0        ", "x,x,x          ", "x,0,0          ",    &
                      "x,1/2,1/2      ", "x,1/2,0        ", "x,x,z          ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(255)= wyck_info_type("F -4 3 M    ", 8,     &
                    (/"0,0,0          ", "1/2,1/2,1/2    ", "1/4,1/4,1/4    ",    &
                      "3/4,3/4,3/4    ", "x,x,x          ", "x,0,0          ",    &
                      "x,1/4,1/4      ", "x,x,z          ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(256)= wyck_info_type("I -4 3 M    ", 7,     &
                    (/"0,0,0          ", "0,1/2,1/2      ", "x,x,x          ",    &
                      "1/4,1/2,0      ", "x,0,0          ", "x,1/2,0        ",    &
                      "x,x,z          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(257)= wyck_info_type("P -4 3 N    ", 8,     &
                    (/"0,0,0          ", "0,1/2,1/2      ", "1/4,1/2,0      ",    &
                      "1/4,0,1/2      ", "x,x,x          ", "x,0,0          ",    &
                      "x,1/2,0        ", "x,0,1/2        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(258)= wyck_info_type("F -4 3 C    ", 7,     &
                    (/"0,0,0          ", "1/4,1/4,1/4    ", "0,1/4,1/4      ",    &
                      "1/4,0,0        ", "x,x,x          ", "x,0,0          ",    &
                      "x,1/4,1/4      ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(259)= wyck_info_type("I -4 3 D    ", 4,     &
                    (/"3/8,0,1/4      ", "7/8,0,1/4      ", "x,x,x          ",    &
                      "x,0,1/4        ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(260)= wyck_info_type("P M -3 M    ",13,     &
                    (/"0,0,0          ", "1/2,1/2,1/2    ", "0,1/2,1/2      ",    &
                      "1/2,0,0        ", "x,0,0          ", "x,1/2,1/2      ",    &
                      "x,x,x          ", "x,1/2,0        ", "0,y,y          ",    &
                      "1/2,y,y        ", "0,y,z          ", "1/2,y,z        ",    &
                      "x,x,z          ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(261)= wyck_info_type("P N -3 N:1  ", 8,     &
                    (/"0,0,0          ", "0,1/2,1/2      ", "1/4,1/4,1/4    ",    &
                      "1/4,0,1/2      ", "x,0,0          ", "x,x,x          ",    &
                      "x,0,1/2        ", "0,y,y          ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(262)= wyck_info_type("P N -3 N    ", 8,     &
                    (/"1/4,1/4,1/4    ", "3/4,1/4,1/4    ", "0,0,0          ",    &
                      "0,3/4,1/4      ", "x,1/4,1/4      ", "x,x,x          ",    &
                      "x,3/4,1/4      ", "1/4,y,y        ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(263)= wyck_info_type("P M -3 N    ",11,     &
                    (/"0,0,0          ", "0,1/2,1/2      ", "1/4,0,1/2      ",    &
                      "1/4,1/2,0      ", "1/4,1/4,1/4    ", "x,0,0          ",    &
                      "x,0,1/2        ", "x,1/2,0        ", "x,x,x          ",    &
                      "1/4,y,y+1/2    ", "0,y,z          ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(264)= wyck_info_type("P N -3 M:1  ",11,     &
                    (/"0,0,0          ", "1/4,1/4,1/4    ", "3/4,3/4,3/4    ",    &
                      "0,1/2,1/2      ", "x,x,x          ", "1/4,0,1/2      ",    &
                      "x,0,0          ", "x,0,1/2        ", "1/4,y,-y+1/2   ",    &
                      "1/4,y,y+1/2    ", "x,x,z          ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(265)= wyck_info_type("P N -3 M    ",11,     &
                    (/"1/4,1/4,1/4    ", "0,0,0          ", "1/2,1/2,1/2    ",    &
                      "1/4,3/4,3/4    ", "x,x,x          ", "1/2,1/4,3/4    ",    &
                      "x,1/4,1/4      ", "x,1/4,3/4      ", "1/2,y,y+1/2    ",    &
                      "1/2,y,-y       ", "x,x,z          ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(266)= wyck_info_type("F M -3 M    ",11,     &
                    (/"0,0,0          ", "1/2,1/2,1/2    ", "1/4,1/4,1/4    ",    &
                      "0,1/4,1/4      ", "x,0,0          ", "x,x,x          ",    &
                      "x,1/4,1/4      ", "0,y,y          ", "1/2,y,y        ",    &
                      "0,y,z          ", "x,x,z          ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(267)= wyck_info_type("F M -3 C    ", 9,     &
                    (/"1/4,1/4,1/4    ", "0,0,0          ", "1/4,0,0        ",    &
                      "0,1/4,1/4      ", "x,0,0          ", "x,1/4,1/4      ",    &
                      "x,x,x          ", "1/4,y,y        ", "0,y,z          ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(268)= wyck_info_type("F D -3 M:1  ", 8,     &
                    (/"0,0,0          ", "1/2,1/2,1/2    ", "1/8,1/8,1/8    ",    &
                      "5/8,5/8,5/8    ", "x,x,x          ", "x,0,0          ",    &
                      "x,x,z          ", "1/8,y,-y+1/4   ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(269)= wyck_info_type("F D -3 M    ", 8,     &
                    (/"1/8,1/8,1/8    ", "3/8,3/8,3/8    ", "0,0,0          ",    &
                      "1/2,1/2,1/2    ", "x,x,x          ", "x,1/8,1/8      ",    &
                      "x,x,z          ", "0,y,-y         ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(270)= wyck_info_type("F D -3 C:1  ", 7,     &
                    (/"0,0,0          ", "1/8,1/8,1/8    ", "3/8,3/8,3/8    ",    &
                      "1/4,0,0        ", "x,x,x          ", "x,0,0          ",    &
                      "1/8,y,-y+1/4   ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(271)= wyck_info_type("F D -3 C    ", 7,     &
                    (/"1/8,1/8,1/8    ", "1/4,1/4,1/4    ", "0,0,0          ",    &
                      "7/8,1/8,1/8    ", "x,x,x          ", "x,1/8,1/8      ",    &
                      "1/4,y,-y       ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
      wyckoff_info(272)= wyck_info_type("I M -3 M    ",11,     &
                    (/"0,0,0          ", "0,1/2,1/2      ", "1/4,1/4,1/4    ",    &
                      "1/4,0,1/2      ", "x,0,0          ", "x,x,x          ",    &
                      "x,0,1/2        ", "0,y,y          ", "1/4,y,-y+1/2   ",    &
                      "0,y,z          ", "x,x,z          ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )
       wyckoff_info(273)= wyck_info_type("I A -3 D    ", 7,     &
                    (/"0,0,0          ", "1/8,1/8,1/8    ", "1/8,0,1/4      ",    &
                      "3/8,0,1/4      ", "x,x,x          ", "x,0,1/4        ",    &
                      "1/8,y,-y+1/4   ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               ", "               ",    &
                      "               ", "               "/) )


       return
    End Subroutine Set_Wyckoff_Info

 End Module CFML_Symmetry_Tables
