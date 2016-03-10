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
!!---- MODULE: CFML_GlobalDeps  (Linux version)
!!----   INFO: Precision for CrysFML library and Operating System information
!!----         All the global variables defined in this module are implicitly public.
!!----
!!---- HISTORY
!!--..    Update: 02/03/2011
!!--..
!!---- VARIABLES
!!--..
!!--..    Operating system
!!--..
!!----    OPS
!!----    OPS_NAME
!!----    OPS_SEP
!!--..
!!--..    Precision Data
!!--..
!!----    SP
!!----    DP
!!----    CP
!!--..
!!--..    Trigonometric
!!--..
!!----    PI
!!----    TO_DEG
!!----    TO_RAD
!!----    TPI
!!--..
!!--..    Numeric
!!--..
!!----    DEPS
!!----    EPS
!!--..
!!---- FUNCTIONS
!!--..
!!----    DIRECTORY_EXISTS
!!----
!!---- SUBROUTINES
!!--..
!!----    WRITE_DATE_TIME
!!----
!!
Module CFML_GlobalDeps

   !---- Variables ----!
   implicit None

   public

   !------------------------------------!
   !---- Operating System variables ----!
   !------------------------------------!

   !!----
   !!---- OPS
   !!----   Integer variable 1: Windows, 2: Linux, 3: MacOs, ....
   !!----   This is a variable set by the user of the library for the case
   !!----   that there is no external library with a procedure for getting
   !!----   the operating system.
   !!----
   !!---- Update: March 2009
   !!
   integer, parameter :: OPS= 2    ! Linux

   !!----
   !!---- OPS_NAME
   !!----   Character variable containing the name of the operating system
   !!----   This is a variable set by the user of the library for the case
   !!----   that there is no external library with a procedure for getting
   !!----   the operating system.
   !!----
   !!---- Update: March 2009
   !!
   character(len=*), parameter :: OPS_NAME="Linux"

   !!----
   !!---- OPS_SEP
   !!----   ASCII code of directory separator character
   !!----   Here it is written explicitly as a character variable
   !!----
   !!---- Update: March 2009
   !!
   character(len=*), parameter :: OPS_SEP="/"

   !------------------------------!
   !---- Precision Parameters ----!
   !------------------------------!

   !!----
   !!---- SP
   !!----    SP: Single precision ( sp = selected_real_kind(6,30) )
   !!----
   !!---- Update: January - 2009
   !!
   integer, parameter :: sp = selected_real_kind(6,30)

   !!----
   !!---- DP
   !!----    DP: Double precision ( dp = selected_real_kind(14,150) )
   !!----
   !!---- Update: January - 2009
   !!
   integer, parameter :: dp = selected_real_kind(14,150)

   !!----
   !!---- CP
   !!----    CP: Current precision
   !!----
   !!---- Update: January - 2009
   !!
   integer, parameter :: cp = sp

   !----------------------------------!
   !---- Trigonometric Parameters ----!
   !----------------------------------!

   !!----
   !!---- PI
   !!----    real(kind=dp), parameter ::  pi = 3.141592653589793238463_dp
   !!----
   !!----    Pi value
   !!----
   !!---- Update: January - 2009
   !!
   real(kind=dp), parameter ::  pi = 3.141592653589793238463_dp

   !!----
   !!---- TO_DEG
   !!----    real(kind=dp), parameter ::  to_DEG = 180.0_dp/pi
   !!----
   !!----    Conversion from Radians to Degrees
   !!----
   !!---- Update: January - 2009
   !!
   real(kind=dp), parameter ::  to_DEG  = 180.0_dp/pi

   !!----
   !!---- TO_RAD
   !!----    real(kind=dp), parameter ::  to_RAD  = pi/180.0_dp
   !!----
   !!----    Conversion from Degrees to Radians
   !!----
   !!---- Update: January - 2009
   !!
   real(kind=dp), parameter ::  to_RAD  = pi/180.0_dp

   !!----
   !!---- TPI
   !!----  real(kind=dp), parameter ::  tpi = 6.283185307179586476925_dp
   !!----
   !!----  2.0*Pi value
   !!----
   !!---- Update: January - 2009
   !!
   real(kind=dp), parameter ::  tpi = 6.283185307179586476925_dp

   !----------------------------!
   !---- Numeric Parameters ----!
   !----------------------------!

   !!----
   !!---- DEPS
   !!----    real(kind=dp), parameter :: deps=0.00000001_dp
   !!----
   !!----    Epsilon value use for comparison of real numbers
   !!----
   !!---- Update: January - 2009
   !!
   real(kind=dp), parameter, public :: deps=0.00000001_dp

   !!----
   !!----  EPS
   !!----     real(kind=cp), public ::  eps=0.00001_cp
   !!----
   !!----     Epsilon value use for comparison of real numbers
   !!----
   !!----  Update: January - 2009
   !!
   real(kind=cp),  parameter, public  ::  eps=0.00001_cp

 Contains

   !-------------------!
   !---- Functions ----!
   !-------------------!

   !!----
   !!---- Function Directory_Exists(Dirname) Result(info)
   !!----    character(len=*), intent(in) :: Dirname
   !!----    logical                      :: info
   !!----
   !!---- Generic function dependent of the compiler that return
   !!---- a logical value if a directory exists or not.
   !!----
   !!---- Update: April - 2009
   !!
   Function Directory_Exists(Dirname) Result(info)
      !---- Argument ----!
      character(len=*), intent(in) :: Dirname
      logical                      :: info

      !---- Local Variables ----!
      character(len=512) :: linea
      integer            :: nlong

      ! Init value
      info=.false.

      linea=trim(dirname)
      nlong=len_trim(linea)
      if (nlong ==0) return

      if (linea(nlong:nlong) /= ops_sep) linea=trim(linea)//ops_sep

      ! All compilers except Intel
      inquire(file=trim(linea)//'.' , exist=info)

      ! Intel
      !inquire(directory=trim(linea), exist=info)

      return
   End Function Directory_Exists

   !---------------------!
   !---- Subroutines ----!
   !---------------------!

   !!----
   !!---- Subroutine Write_Date_Time(lun,dtim)
   !!----  integer,         optional,intent(in) :: lun
   !!----  character(len=*),optional,intent(out):: dtim
   !!----
   !!---- Generic subroutine for writing the date and time
   !!---- in form   Date: Day/Month/Year  Time: hour:minute:second
   !!---- to a file with logical unit = lun. The output argument
   !!---- can be provided to get a string with the same information
   !!----
   !!---- Updated: January - 2014
   !!
   Subroutine Write_Date_Time(lun,dtim)
     integer,         optional,intent(in) :: lun
     character(len=*),optional,intent(out):: dtim
     !--- Local variables ----!
     character (len=10) :: dat
     character (len=10) :: tim
     call date_and_time(date=dat,time=tim)
     if(present(lun)) &
     write(unit=lun,fmt="(/,4a)") &
       " => Date: ",dat(7:8)//"/"//dat(5:6)//"/"//dat(1:4),      &
         "  Time: ",tim(1:2)//":"//tim(3:4)//":"//tim(5:10)
     if(present(dtim)) &
      dtim="#   Date: "//dat(7:8)//"/"//dat(5:6)//"/"//dat(1:4)//      &
            "  Time: "//tim(1:2)//":"//tim(3:4)//":"//tim(5:10)
     return
   End Subroutine Write_Date_Time

End Module CFML_GlobalDeps
