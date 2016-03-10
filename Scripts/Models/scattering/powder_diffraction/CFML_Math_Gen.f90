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
!!---- MODULE: CFML_Math_General
!!----   INFO: Mathematic general utilities for use in Crystallography and
!!----         Solid State Physics and Chemistry.
!!----
!!---- HISTORY
!!----    Updated: 02/03/2011
!!----
!!---- DEPENDENCIES
!!----
!!----    CFML_GlobalDeps
!!----
!!---- VARIABLES
!!--++    EPSS                         [Private]
!!--++    EP_SS                        [Private]
!!----    ERR_MathGen
!!----    ERR_MathGen_Mess
!!----
!!---- PROCEDURES
!!----    Functions:
!!--..
!!--..    Trigonometric Functions
!!----       ACOSD
!!--++       ACOSD_dp                  [Overloaded]
!!--++       ACOSD_sp                  [Overloaded]
!!----       ASIND
!!--++       ASIND_dp                  [Overloaded]
!!--++       ASIND_sp                  [Overloaded]
!!----       ATAN2D
!!--++       ATAN2D_dp                 [Overloaded]
!!--++       ATAN2D_sp                 [Overloaded]
!!----       ATAND
!!--++       ATAND_dp                  [Overloaded]
!!--++       ATAND_sp                  [Overloaded]
!!----       COSD
!!--++       COSD_dp                   [Overloaded]
!!--++       COSD_sp                   [Overloaded]
!!----       SIND
!!--++       SIND_dp                   [Overloaded]
!!--++       SIND_sp                   [Overloaded]
!!----       TAND
!!--++       TAND_dp                   [Overloaded]
!!--++       TAND_sp                   [Overloaded]
!!--..
!!--..    Special Functions
!!----       BESSJ0
!!----       BESSJ1
!!----       BESSJ
!!--..
!!--..    Scalar Functions
!!----       FACTORIAL
!!----       NEGLIGIBLE
!!--++       NEGLIGIBLEC               [Overloaded]
!!--++       NEGLIGIBLER               [Overloaded]
!!----       PGCD
!!----       PPCM
!!----       PYTHAG
!!--++       PYTHAG_dp                 [Overloaded]
!!--++       PYTHAG_sp                 [Overloaded]
!!--..
!!--..    Arrays and Vectors Functions
!!----       CO_LINEAR
!!--++       CO_LINEAR_C               [Overloaded]
!!--++       CO_LINEAR_I               [Overloaded]
!!--++       CO_LINEAR_R               [Overloaded]
!!----       CO_PRIME
!!----       EQUAL_MATRIX
!!--++       EQUAL_MATRIX_I            [Overloaded]
!!--++       EQUAL_MATRIX_R            [Overloaded]
!!----       EQUAL_VECTOR
!!--++       EQUAL_VECTOR_I            [Overloaded]
!!--++       EQUAL_VECTOR_R            [Overloaded]
!!----       EUCLIDEAN_NORM
!!----       IMAXLOC
!!--++       IMAXLOC_I                 [Overloaded]
!!--++       IMAXLOC_R                 [OVerloaded]
!!----       IMINLOC
!!--++       IMINLOC_I                 [Overloaded]
!!--++       IMINLOC_R                 [OVerloaded]
!!----       LOCATE
!!--++       LOCATE_I                  [Overloaded]
!!--++       LOCATE_IB                 [Overloaded]
!!--++       LOCATE_R                  [Overloaded]
!!--++       LOCATE_RB                 [Overloaded]
!!----       LOWER_TRIANGULAR
!!--++       LOWER_TRIANGULAR_I        [Overloaded]
!!--++       LOWER_TRIANGULAR_R        [Overloaded]
!!----       MODULO_LAT
!!----       NORM
!!--++       NORM_I                    [Overloaded]
!!--++       NORM_R                    [Overloaded]
!!----       OUTERPROD
!!--++       OUTERPROD_dp              [Overloaded]
!!--++       OUTERPROD_sp              [Overloaded]
!!----       SCALAR
!!--++       SCALAR_I                  [Overloaded]
!!--++       SCALAR_R                  [Overloaded]
!!----       TRACE
!!--++       TRACE_C                   [Overloaded]
!!--++       TRACE_I                   [Overloaded]
!!--++       TRACE_R                   [Overloaded]
!!----       UPPER_TRIANGULAR
!!--++       UPPER_TRIANGULAR_I        [Overloaded]
!!--++       UPPER_TRIANGULAR_R        [Overloaded]
!!----       ZBELONG
!!--++       ZBELONGM                  [Overloaded]
!!--++       ZBELONGN                  [Overloaded]
!!--++       ZBELONGV                  [Overloaded]
!!--..
!!----
!!----    Subroutines:
!!--..
!!--..    Init Routine
!!----       INIT_ERR_MATHGEN
!!----       SET_EPSG
!!----       SET_EPSG_DEFAULT
!!--..
!!--..    Trigonometric Subroutines
!!----       RTAN
!!--++       RTAN_dp                   [Overloaded]
!!--++       RTAN_sp                   [Overloaded]
!!--..
!!--..    Arrays and Vectors Functions
!!----       CO_PRIME_VECTOR
!!----       DETERMINANT
!!--++       DETERMINANT_C             [Overloaded]
!!--++       DETERMINANT_R             [Overloaded]
!!----       DIAGONALIZE_SH
!!--++       DIAGONALIZE_HERM          [Overloaded]
!!--++       DIAGONALIZE_SYMM          [Overloaded]
!!--++       EIGSRT                    [Private]
!!----       FIRST_DERIVATIVE
!!----       IN_SORT
!!----       INVERT_MATRIX
!!----       LINEAR_DEPENDENT
!!--++       LINEAR_DEPENDENTC         [Overloaded]
!!--++       LINEAR_DEPENDENTI         [Overloaded]
!!--++       LINEAR_DEPENDENTR         [Overloaded]
!!----       LU_BACKSUB
!!----       LU_DECOMP
!!----       MATINV
!!--++       PARTITION                 [Private]
!!----       POINTS_IN_LINE2D
!!----       RANK
!!--++       RANK_dp                   [Overloaded]
!!--++       RANK_sp                   [Overloaded]
!!----       SECOND_DERIVATIVE
!!----       SMOOTHINGVEC
!!----       SORT
!!--++       SORT_I                    [Overloaded]
!!--++       SORT_R                    [Overloaded]
!!----       SORT_STRINGS
!!----       SPLINE
!!----       SPLINT
!!----       SVDCMP
!!--++       SVDCMP_dp                 [Overloaded]
!!--++       SVDCMP_sp                 [Overloaded]
!!----       SWAP
!!--++       SWAP_C                    [Overloaded]
!!--++       SWAP_CM                   [Overloaded]
!!--++       SWAP_CV                   [Overloaded]
!!--++       SWAP_I                    [Overloaded]
!!--++       SWAP_IM                   [Overloaded]
!!--++       SWAP_IV                   [Overloaded]
!!--++       SWAP_R                    [Overloaded]
!!--++       SWAP_RM                   [Overloaded]
!!--++       SWAP_RV                   [Overloaded]
!!--++       MASKED_SWAP_R             [Overloaded]
!!--++       MASKED_SWAP_RM            [Overloaded]
!!--++       MASKED_SWAP_RV            [Overloaded]
!!--++       TQLI1                     [Private]
!!--++       TQLI2                     [Private]
!!--++       TRED1                     [Private]
!!--++       TRED2                     [Private]
!!--++
!!
 Module CFML_Math_General
    !---- Use Modules ----!
    Use CFML_GlobalDeps

    !---- Variables ----!
    implicit none

    private

    !---- List of public functions ----!
    public :: Bessj0, Bessj1, Bessj, Factorial, Pgcd, Ppcm, Modulo_Lat, Co_Prime, &
              Euclidean_Norm,Erf

    !---- List of public overloaded procedures: functions ----!
    public :: Acosd, Asind, Atan2d, Atand, Cosd, Sind, Tand, Negligible, Pythag,   &
              Co_Linear, Equal_Matrix, Equal_Vector, Locate, Outerprod, Trace,     &
              Zbelong, Imaxloc, Iminloc, Norm, Scalar, In_limits, Lower_Triangular,&
              Upper_Triangular

    !---- List of private functions ----!
    private :: Acosd_dp, Acosd_sp, Asind_dp, Asind_sp, Atan2d_dp, Atan2d_sp,       &
               Atand_dp, Atand_sp, Cosd_dp, Cosd_sp, Sind_dp, Sind_sp, Tand_dp,    &
               Tand_sp, Negligiblec, Negligibler, Pythag_dp, Pythag_sp,            &
               Co_linear_C, Co_linear_I, Co_linear_R, Equal_Matrix_I,              &
               Equal_Matrix_R, Equal_Vector_I, Equal_Vector_R, Locate_I, Locate_R, &
               Outerprod_dp, Outerprod_sp, Trace_C, Trace_I, Trace_R, ZbelongM,    &
               ZbelongN, ZbelongV, Imaxloc_I, Imaxloc_R, Iminloc_R, Iminloc_I,     &
               Norm_I, Norm_R, Scalar_I, Scalar_R, Locate_Ib, Locate_Rb,           &
               In_limits_dp, In_limits_sp, In_Limits_int, Lower_Triangular_I,      &
               Lower_Triangular_R, Upper_Triangular_I, Upper_Triangular_R

    !---- List of public subroutines ----!
    public ::  Init_Err_Mathgen, Invert_Matrix, LU_Decomp, LU_Backsub, Matinv,        &
               Sort_Strings, Spline, Splint, Set_Epsg, Set_Epsg_Default,In_Sort,      &
               First_Derivative, Second_Derivative, SmoothingVec, Points_in_Line2D,   &
               Co_Prime_vector

    !---- List of public overloaded procedures: subroutines ----!
    public ::  RTan, Determinant, Diagonalize_Sh, Linear_Dependent, Rank, Sort,   &
               Svdcmp, Swap

    !---- List of private subroutines ----!
    private :: RTan_dp, RTan_sp, Determinant_C,Determinant_R, Diagonalize_Herm,   &
               Diagonalize_Symm, Eigsrt, Linear_DependentC, Linear_DependentI,    &
               Linear_DependentR, Rank_dp, Rank_sp, Sort_I, Sort_R, Svdcmp_dp,    &
               Svdcmp_sp, Swap_C, Swap_Cm, Swap_Cv, Swap_I, Swap_Im, Swap_Iv,     &
               Swap_R, Swap_Rm, Swap_Rv, Masked_Swap_R, Masked_Swap_Rm,           &
               Masked_Swap_Rv, Tqli1, Tqli2, Tred1, Tred2, Partition


    !---- Definitions ----!

    !!--++
    !!--++ EPSS
    !!--++    real(kind=cp)  :: epss=1.0E-5_cp
    !!--++
    !!--++    Internal epsilon value used for comparing reals to integers
    !!--++    in crystallographic applications where the maximum precision in the
    !!--++    measured values is of the order of 10^-5.
    !!--++
    !!--++ Update: April - 2005
    !!
    real(kind=cp),   private :: epss=1.0E-5_cp

    !!--++
    !!--++ EP_SS
    !!--++    real(kind=cp), parameter, private  :: ep_ss=1.0E-12_cp
    !!--++
    !!--++    Internal epsilon value used for comparison in matrix operations
    !!--++
    !!--++ Update: February - 2005
    !!
    real(kind=cp), parameter, private :: ep_ss=1.0E-12_cp

    !!----
    !!---- ERR_MathGen
    !!----    logical :: ERR_MathGen
    !!----
    !!----    Logical Variable indicating an error in CFML_Math_General module
    !!----
    !!---- Update: February - 2005
    !!
    logical, public :: ERR_MathGen

    !!----
    !!---- ERR_MathGen_Mess
    !!----    character(len=150) :: ERR_MathGen_Mess
    !!----
    !!----    String containing information about the last error
    !!----
    !!---- Update: February - 2005
    !!
    character(len=150), public:: ERR_MathGen_Mess

    !!----
    !!---- Primes
    !!----    integer, parameter, dimension(1000), public :: primes
    !!----
    !!----    List of the first 1000 prime numbers.
    !!----    Used by the subroutine Co_Prime_Vector and function Co_Prime
    !!----
    !!----  Created: January - 2011
    !!
    integer, parameter, dimension(1000), public :: primes =                                       &
           (/ 2,      3,      5,      7,     11,     13,     17,     19,     23,     29,  &
             31,     37,     41,     43,     47,     53,     59,     61,     67,     71,  &
             73,     79,     83,     89,     97,    101,    103,    107,    109,    113,  &
            127,    131,    137,    139,    149,    151,    157,    163,    167,    173,  &
            179,    181,    191,    193,    197,    199,    211,    223,    227,    229,  &
            233,    239,    241,    251,    257,    263,    269,    271,    277,    281,  &
            283,    293,    307,    311,    313,    317,    331,    337,    347,    349,  &
            353,    359,    367,    373,    379,    383,    389,    397,    401,    409,  &
            419,    421,    431,    433,    439,    443,    449,    457,    461,    463,  &
            467,    479,    487,    491,    499,    503,    509,    521,    523,    541,  &
            547,    557,    563,    569,    571,    577,    587,    593,    599,    601,  &
            607,    613,    617,    619,    631,    641,    643,    647,    653,    659,  &
            661,    673,    677,    683,    691,    701,    709,    719,    727,    733,  &
            739,    743,    751,    757,    761,    769,    773,    787,    797,    809,  &
            811,    821,    823,    827,    829,    839,    853,    857,    859,    863,  &
            877,    881,    883,    887,    907,    911,    919,    929,    937,    941,  &
            947,    953,    967,    971,    977,    983,    991,    997,   1009,   1013,  &
           1019,   1021,   1031,   1033,   1039,   1049,   1051,   1061,   1063,   1069,  &
           1087,   1091,   1093,   1097,   1103,   1109,   1117,   1123,   1129,   1151,  &
           1153,   1163,   1171,   1181,   1187,   1193,   1201,   1213,   1217,   1223,  &
           1229,   1231,   1237,   1249,   1259,   1277,   1279,   1283,   1289,   1291,  &
           1297,   1301,   1303,   1307,   1319,   1321,   1327,   1361,   1367,   1373,  &
           1381,   1399,   1409,   1423,   1427,   1429,   1433,   1439,   1447,   1451,  &
           1453,   1459,   1471,   1481,   1483,   1487,   1489,   1493,   1499,   1511,  &
           1523,   1531,   1543,   1549,   1553,   1559,   1567,   1571,   1579,   1583,  &
           1597,   1601,   1607,   1609,   1613,   1619,   1621,   1627,   1637,   1657,  &
           1663,   1667,   1669,   1693,   1697,   1699,   1709,   1721,   1723,   1733,  &
           1741,   1747,   1753,   1759,   1777,   1783,   1787,   1789,   1801,   1811,  &
           1823,   1831,   1847,   1861,   1867,   1871,   1873,   1877,   1879,   1889,  &
           1901,   1907,   1913,   1931,   1933,   1949,   1951,   1973,   1979,   1987,  &
           1993,   1997,   1999,   2003,   2011,   2017,   2027,   2029,   2039,   2053,  &
           2063,   2069,   2081,   2083,   2087,   2089,   2099,   2111,   2113,   2129,  &
           2131,   2137,   2141,   2143,   2153,   2161,   2179,   2203,   2207,   2213,  &
           2221,   2237,   2239,   2243,   2251,   2267,   2269,   2273,   2281,   2287,  &
           2293,   2297,   2309,   2311,   2333,   2339,   2341,   2347,   2351,   2357,  &
           2371,   2377,   2381,   2383,   2389,   2393,   2399,   2411,   2417,   2423,  &
           2437,   2441,   2447,   2459,   2467,   2473,   2477,   2503,   2521,   2531,  &
           2539,   2543,   2549,   2551,   2557,   2579,   2591,   2593,   2609,   2617,  &
           2621,   2633,   2647,   2657,   2659,   2663,   2671,   2677,   2683,   2687,  &
           2689,   2693,   2699,   2707,   2711,   2713,   2719,   2729,   2731,   2741,  &
           2749,   2753,   2767,   2777,   2789,   2791,   2797,   2801,   2803,   2819,  &
           2833,   2837,   2843,   2851,   2857,   2861,   2879,   2887,   2897,   2903,  &
           2909,   2917,   2927,   2939,   2953,   2957,   2963,   2969,   2971,   2999,  &
           3001,   3011,   3019,   3023,   3037,   3041,   3049,   3061,   3067,   3079,  &
           3083,   3089,   3109,   3119,   3121,   3137,   3163,   3167,   3169,   3181,  &
           3187,   3191,   3203,   3209,   3217,   3221,   3229,   3251,   3253,   3257,  &
           3259,   3271,   3299,   3301,   3307,   3313,   3319,   3323,   3329,   3331,  &
           3343,   3347,   3359,   3361,   3371,   3373,   3389,   3391,   3407,   3413,  &
           3433,   3449,   3457,   3461,   3463,   3467,   3469,   3491,   3499,   3511,  &
           3517,   3527,   3529,   3533,   3539,   3541,   3547,   3557,   3559,   3571,  &
           3581,   3583,   3593,   3607,   3613,   3617,   3623,   3631,   3637,   3643,  &
           3659,   3671,   3673,   3677,   3691,   3697,   3701,   3709,   3719,   3727,  &
           3733,   3739,   3761,   3767,   3769,   3779,   3793,   3797,   3803,   3821,  &
           3823,   3833,   3847,   3851,   3853,   3863,   3877,   3881,   3889,   3907,  &
           3911,   3917,   3919,   3923,   3929,   3931,   3943,   3947,   3967,   3989,  &
           4001,   4003,   4007,   4013,   4019,   4021,   4027,   4049,   4051,   4057,  &
           4073,   4079,   4091,   4093,   4099,   4111,   4127,   4129,   4133,   4139,  &
           4153,   4157,   4159,   4177,   4201,   4211,   4217,   4219,   4229,   4231,  &
           4241,   4243,   4253,   4259,   4261,   4271,   4273,   4283,   4289,   4297,  &
           4327,   4337,   4339,   4349,   4357,   4363,   4373,   4391,   4397,   4409,  &
           4421,   4423,   4441,   4447,   4451,   4457,   4463,   4481,   4483,   4493,  &
           4507,   4513,   4517,   4519,   4523,   4547,   4549,   4561,   4567,   4583,  &
           4591,   4597,   4603,   4621,   4637,   4639,   4643,   4649,   4651,   4657,  &
           4663,   4673,   4679,   4691,   4703,   4721,   4723,   4729,   4733,   4751,  &
           4759,   4783,   4787,   4789,   4793,   4799,   4801,   4813,   4817,   4831,  &
           4861,   4871,   4877,   4889,   4903,   4909,   4919,   4931,   4933,   4937,  &
           4943,   4951,   4957,   4967,   4969,   4973,   4987,   4993,   4999,   5003,  &
           5009,   5011,   5021,   5023,   5039,   5051,   5059,   5077,   5081,   5087,  &
           5099,   5101,   5107,   5113,   5119,   5147,   5153,   5167,   5171,   5179,  &
           5189,   5197,   5209,   5227,   5231,   5233,   5237,   5261,   5273,   5279,  &
           5281,   5297,   5303,   5309,   5323,   5333,   5347,   5351,   5381,   5387,  &
           5393,   5399,   5407,   5413,   5417,   5419,   5431,   5437,   5441,   5443,  &
           5449,   5471,   5477,   5479,   5483,   5501,   5503,   5507,   5519,   5521,  &
           5527,   5531,   5557,   5563,   5569,   5573,   5581,   5591,   5623,   5639,  &
           5641,   5647,   5651,   5653,   5657,   5659,   5669,   5683,   5689,   5693,  &
           5701,   5711,   5717,   5737,   5741,   5743,   5749,   5779,   5783,   5791,  &
           5801,   5807,   5813,   5821,   5827,   5839,   5843,   5849,   5851,   5857,  &
           5861,   5867,   5869,   5879,   5881,   5897,   5903,   5923,   5927,   5939,  &
           5953,   5981,   5987,   6007,   6011,   6029,   6037,   6043,   6047,   6053,  &
           6067,   6073,   6079,   6089,   6091,   6101,   6113,   6121,   6131,   6133,  &
           6143,   6151,   6163,   6173,   6197,   6199,   6203,   6211,   6217,   6221,  &
           6229,   6247,   6257,   6263,   6269,   6271,   6277,   6287,   6299,   6301,  &
           6311,   6317,   6323,   6329,   6337,   6343,   6353,   6359,   6361,   6367,  &
           6373,   6379,   6389,   6397,   6421,   6427,   6449,   6451,   6469,   6473,  &
           6481,   6491,   6521,   6529,   6547,   6551,   6553,   6563,   6569,   6571,  &
           6577,   6581,   6599,   6607,   6619,   6637,   6653,   6659,   6661,   6673,  &
           6679,   6689,   6691,   6701,   6703,   6709,   6719,   6733,   6737,   6761,  &
           6763,   6779,   6781,   6791,   6793,   6803,   6823,   6827,   6829,   6833,  &
           6841,   6857,   6863,   6869,   6871,   6883,   6899,   6907,   6911,   6917,  &
           6947,   6949,   6959,   6961,   6967,   6971,   6977,   6983,   6991,   6997,  &
           7001,   7013,   7019,   7027,   7039,   7043,   7057,   7069,   7079,   7103,  &
           7109,   7121,   7127,   7129,   7151,   7159,   7177,   7187,   7193,   7207,  &
           7211,   7213,   7219,   7229,   7237,   7243,   7247,   7253,   7283,   7297,  &
           7307,   7309,   7321,   7331,   7333,   7349,   7351,   7369,   7393,   7411,  &
           7417,   7433,   7451,   7457,   7459,   7477,   7481,   7487,   7489,   7499,  &
           7507,   7517,   7523,   7529,   7537,   7541,   7547,   7549,   7559,   7561,  &
           7573,   7577,   7583,   7589,   7591,   7603,   7607,   7621,   7639,   7643,  &
           7649,   7669,   7673,   7681,   7687,   7691,   7699,   7703,   7717,   7723,  &
           7727,   7741,   7753,   7757,   7759,   7789,   7793,   7817,   7823,   7829,  &
           7841,   7853,   7867,   7873,   7877,   7879,   7883,   7901,   7907,   7919 /)

    !---- Interfaces - Overloaded ----!
    Interface  Acosd
       Module Procedure Acosd_dp
       Module Procedure Acosd_sp
    End Interface

    Interface  Asind
       Module Procedure Asind_dp
       Module Procedure Asind_sp
    End Interface

    Interface  Atan2d
       Module Procedure Atan2d_dp
       Module Procedure Atan2d_sp
    End Interface

    Interface  Atand
       Module Procedure Atand_dp
       Module Procedure Atand_sp
    End Interface

    Interface  Cosd
       Module Procedure Cosd_dp
       Module Procedure Cosd_sp
    End Interface

    Interface  Sind
       Module Procedure Sind_dp
       Module Procedure Sind_sp
    End Interface

    Interface  Tand
       Module Procedure Tand_dp
       Module Procedure Tand_sp
    End Interface

    Interface  Negligible
       Module Procedure Negligibler
       Module Procedure Negligiblec
    End Interface

    Interface  Pythag
       Module Procedure Pythag_dp
       Module Procedure Pythag_sp
    End Interface

    Interface  Co_Linear
       Module Procedure Co_linear_C
       Module Procedure Co_linear_I
       Module Procedure Co_linear_R
    End Interface

    Interface  Equal_Matrix
       Module Procedure Equal_Matrix_I
       Module Procedure Equal_Matrix_R
    End Interface

    Interface  Equal_Vector
       Module Procedure Equal_Vector_I
       Module Procedure Equal_Vector_R
    End Interface

    Interface  IMaxloc
       Module Procedure IMaxloc_I
       Module Procedure IMaxloc_R
    End Interface

    Interface  IMinloc
       Module Procedure IMinloc_I
       Module Procedure IMinloc_R
    End Interface

    Interface  Locate
       Module Procedure Locate_I
       Module Procedure Locate_R
       Module Procedure Locate_Ib
       Module Procedure Locate_Rb
    End Interface

    Interface  Lower_Triangular
       Module Procedure Lower_Triangular_I
       Module Procedure Lower_Triangular_R
    End Interface

    Interface Norm
       Module Procedure Norm_I
       Module Procedure Norm_R
    End Interface Norm

    Interface  Outerprod
       Module Procedure Outerprod_dp
       Module Procedure Outerprod_sp
    End Interface

    Interface Scalar
       Module Procedure Scalar_I
       Module Procedure Scalar_R
    End Interface Scalar

    Interface  Trace
       Module Procedure Trace_C
       Module Procedure Trace_I
       Module Procedure Trace_R
    End Interface

    Interface  Upper_Triangular
       Module Procedure Upper_Triangular_I
       Module Procedure Upper_Triangular_R
    End Interface

    Interface  Zbelong
       Module Procedure ZbelongM
       Module Procedure ZbelongN
       Module Procedure ZbelongV
    End Interface

    Interface  Rtan
       Module Procedure Rtan_dp
       Module Procedure Rtan_sp
    End Interface

    Interface  Determinant
       Module Procedure Determinant_c
       Module Procedure Determinant_r
    End Interface

    Interface  Diagonalize_SH
       Module Procedure Diagonalize_HERM
       Module Procedure Diagonalize_SYMM
    End Interface

    Interface In_Limits
       Module Procedure In_Limits_int
       Module Procedure In_Limits_dp
       Module Procedure In_Limits_sp
    End Interface

    Interface  Linear_Dependent
       Module Procedure Linear_Dependentc
       Module Procedure Linear_Dependenti
       Module Procedure Linear_Dependentr
    End Interface

    Interface  Rank
       Module Procedure Rank_dp
       Module Procedure Rank_sp
    End Interface

    Interface  Sort
       Module Procedure Sort_I
       Module Procedure Sort_R
    End Interface

    Interface  Svdcmp
       Module Procedure Svdcmp_dp
       Module Procedure Svdcmp_sp
    End Interface

    Interface Swap
        Module Procedure swap_c
        Module Procedure swap_cm
        Module Procedure swap_cv
        Module Procedure swap_i
        Module Procedure swap_im
        Module Procedure swap_iv
        Module Procedure swap_r
        Module Procedure swap_rm
        Module Procedure swap_rv
        Module Procedure masked_swap_r
        Module Procedure masked_swap_rm
        Module Procedure masked_swap_rv
    End interface

 Contains

    !---- Functions ----!

    !!----
    !!---- Elemental Function Acosd(x) Result(arc_cos)
    !!----    real(kind=sp/dp), intent(in) :: x
    !!----    real(kind=sp/dp)             :: arc_cos
    !!----
    !!----    Inverse cosine function -> output in Degrees
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Elemental Function Acosd_dp(x) Result(arc_cos)
    !!--++    real(kind=dp), intent(in) :: x
    !!--++    real(kind=dp)             :: arc_cos
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse cosine function -> output in Degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Acosd_dp(x) Result(arc_cos)
       !---- Argument ----!
       real(kind=dp), intent(in) :: x
       real(kind=dp)             :: arc_cos

       if (abs(x) > 1.0_dp ) then
          if (x > 0.0_dp)  then
             arc_cos=0.0_dp
          else
             arc_cos=180.0_dp
          end if
       else
          arc_cos=acos(x)*to_DEG
       end if

       return
    End Function Acosd_dp

    !!--++
    !!--++ Elemental Function Acosd_sp(x) Result(arc_cos)
    !!--++    real(kind=sp), intent(in) :: x
    !!--++    real(kind=sp)             :: arc_cos
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse cosine function -> output in Degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Acosd_sp(x) Result(arc_cos)
       !---- Argument ----!
       real(kind=sp), intent(in) :: x
       real(kind=sp)             :: arc_cos

       if (abs(x) > 1.0_sp ) then
          if (x > 0.0_sp)  then
             arc_cos=0.0_sp
          else
             arc_cos=180.0_sp
          end if
       else
          arc_cos=acos(x)*to_DEG
       end if

       return
    End Function Acosd_sp

    !!----
    !!---- Function Asind(x) Result(arc_sin)
    !!----    real(kind=sp/dp), intent(in) :: x
    !!----    real(kind=sp/dp)             :: arc_sin
    !!----
    !!----    Inverse sine function -> output in Degrees
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Elemental Function Asind_dp(x) result(arc_sin)
    !!--++    real(kind=dp), intent(in) :: x
    !!--++    real(kind=dp)             :: arc_sin
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse sine function -> output in Degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Asind_dp(x) Result(arc_sin)
       !---- Argument ----!
       real(kind=dp), intent(in) :: x
       real(kind=dp)             :: arc_sin

       if (abs(x) > 1.0_dp ) then
          if (x > 0.0_dp) then
             arc_sin=90.0_dp
          else
             arc_sin=-90.0_dp
          end if
       else
          arc_sin=asin(x)*to_DEG
       end if

       return
    End Function Asind_dp

    !!--++
    !!--++ Elemental Function Asind_sp(x) result(arc_sin)
    !!--++    real(kind=sp), intent(in) :: x
    !!--++    real(kind=sp)             :: arc_sin
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse sine function -> output in Degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Asind_sp(x) Result(arc_sin)
       !---- Argument ----!
       real(kind=sp), intent(in) :: x
       real(kind=sp)             :: arc_sin

       if (abs(x) > 1.0_sp ) then
          if (x > 0.0_sp) then
             arc_sin=90.0_sp
          else
             arc_sin=-90.0_sp
          end if
       else
          arc_sin=asin(x)*to_DEG
       end if

       return
    End Function Asind_sp

    !!----
    !!---- Elemental Function Atan2d(y,x) Result(atande)
    !!----    real(kind=sp/dp), intent(in) :: y,x
    !!----    real(kind=sp/dp)             :: atande
    !!----
    !!----    Inverse tangent function of y/x
    !!----    y,x have the same units -> output in Degrees
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Elemental Function Atan2d_dp(y,x) Result(atande)
    !!--++    real(kind=dp), intent(in) :: y,x
    !!--++    real(kind=dp)             :: atande
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse tangent function of y/x
    !!--++    y,x have the same units -> output in Degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Atan2d_dp(y,x) Result(atand)
       !---- Argument ----!
       real(kind=dp), intent(in) :: y,x
       real(kind=dp)             :: atand

       atand=atan2(y,x)*to_DEG

       return
    End Function Atan2d_dp

    !!--++
    !!--++ Elemental Function Atan2d_sp(y,x) Result(atande)
    !!--++    real(kind=sp), intent(in) :: y,x
    !!--++    real(kind=sp)             :: atande
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse tangent function of y/x
    !!--++    y,x have the same units -> output in Degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Atan2d_sp(y,x) Result(atande)
       !---- Argument ----!
       real(kind=sp), intent(in) :: y,x
       real(kind=sp)             :: atande

       atande=atan2(y,x)*to_DEG

       return
    End Function Atan2d_sp

    !!----
    !!---- Elemental Function Atand(x) Result(atande)
    !!----    real(kind=sp/dp), intent(in) :: x
    !!----    real(kind=sp/dp)             :: atande
    !!----
    !!----    Inverse tangent function, X no units -> output in Degrees
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Elemental Function Atand_dp(x) result(atande)
    !!--++    real(kind=dp), intent(in) :: x
    !!--++    real(kind=dp)             :: atande
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse tangent function, X no units -> output in Degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Atand_dp(x) Result(atand)
       !---- Argument ----!
       real(kind=dp), intent(in) :: x
       real(kind=dp)             :: atand

       atand=atan(x)*to_DEG

       return
    End Function Atand_dp

    !!--++
    !!--++ Function Atand_sp(x) result(atande)
    !!--++    real(kind=sp), intent(in) :: x
    !!--++    real(kind=sp)             :: atande
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse tangent function, X no units -> output in Degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Atand_sp(x) Result(atande)
       !---- Argument ----!
       real(kind=sp), intent(in) :: x
       real(kind=sp)             :: atande

       atande=atan(x)*to_DEG

       return
    End Function Atand_sp



    !!----
    !!---- Elemental Function Cosd(x) Result(cosine)
    !!----    real(kind=sp/dp), intent(in) :: x
    !!----    real(kind=sp/dp)             :: cosine
    !!----
    !!----    Cosine function, X in degrees
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Elemental Function Cosd_dp(x) Result(cosine)
    !!--++    real(kind=dp), intent(in) :: x
    !!--++    real(kind=dp)             :: cosine
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Cosine function, X in degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Cosd_dp(x) Result(cosine)
       !---- Argument ----!
       real(kind=dp), intent(in) :: x
       real(kind=dp)             :: cosine

       cosine=cos(to_RAD*x)

       return
    End Function Cosd_dp

    !!--++
    !!--++ Elemental Function Cosd_sp(x) Result(cosine)
    !!--++    real(kind=sp), intent(in) :: x
    !!--++    real(kind=sp)             :: cosine
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Cosine function, X in degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Cosd_sp(x) Result(cosine)
       !---- Argument ----!
       real(kind=sp), intent(in) :: x
       real(kind=sp)             :: cosine

       cosine=cos(to_RAD*x)

       return
    End Function Cosd_sp

    !!----
    !!---- Elemental Function Sind(x) Result(sine)
    !!----    real(kind=sp/dp), intent(in) :: x
    !!----    real(kind=sp/dp)             :: sine
    !!----
    !!----    Sine function, X in degrees
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Elemental Function Sind_dp(x) Result(sine)
    !!--++    real(kind=dp), intent(in) :: x
    !!--++    real(kind=dp)             :: sine
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Sine function, X in degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Sind_dp(x) Result(sine)
       !---- Argument ----!
       real(kind=dp), intent(in) :: x
       real(kind=dp)             :: sine

       sine=sin(to_RAD*x)

       return
    End Function Sind_dp

    !!--++
    !!--++ Elemental Function Sind_sp(x) Result(sine)
    !!--++    real(kind=sp), intent(in) :: x
    !!--++    real(kind=sp)             :: sine
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Sine function, X in degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Sind_sp(x) Result(sine)
       !---- Argument ----!
       real(kind=sp), intent(in) :: x
       real(kind=sp)             :: sine

       sine=sin(to_RAD*x)

       return
    End Function Sind_sp

    !!----
    !!---- Elemental Function Tand(x) Result(tande)
    !!----    real(kind=sp/dp), intent(in) :: x
    !!----    real(kind=sp/dp)             :: tande
    !!----
    !!----    Tangent function, X in degrees
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Elemental Function Tand_dp(x) Result(tande)
    !!--++    real(kind=dp), intent(in) :: x
    !!--++    real(kind=dp)             :: tande
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Tangent function, X in degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Tand_dp(x) Result(tand)
       !---- Argument ----!
       real(kind=dp), intent(in) :: x
       real(kind=dp)             :: tand

       tand=tan(to_RAD*x)

       return
    End Function Tand_dp

    !!--++
    !!--++ Elemental Function Tand_sp(x) Result(tande)
    !!--++    real(kind=sp), intent(in) :: x
    !!--++    real(kind=sp)             :: tande
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Tangent function, X in degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Tand_sp(x) Result(tande)
       !---- Argument ----!
       real(kind=sp), intent(in) :: x
       real(kind=sp)             :: tande

       tande=tan(to_RAD*x)

       return
    End Function Tand_sp

    !!----
    !!---- Elemental Function BessJ0(x) Result(bessj_0)
    !!----    real(kind=sp), intent(in) :: x
    !!----    real(kind=sp)             :: bessj_0
    !!----
    !!----    Bessel Fuction J0(x)
    !!----
    !!---- Update: February - 2005
    !!
    Elemental Function BessJ0(x) Result(bessj_0)
       !---- Arguments ----!
       real(kind=cp), intent(in) :: x
       real(kind=cp)             :: bessj_0

       !---- Local variables ----!
       real(kind=dp), parameter :: p1=   1.0_dp
       real(kind=dp), parameter :: p2=  -0.1098628627e-2_dp
       real(kind=dp), parameter :: p3=   0.2734510407e-4_dp
       real(kind=dp), parameter :: p4=  -0.2073370639e-5_dp
       real(kind=dp), parameter :: p5=   0.2093887211e-6_dp
       real(kind=dp), parameter :: q1=  -0.1562499995e-1_dp
       real(kind=dp), parameter :: q2=   0.1430488765e-3_dp
       real(kind=dp), parameter :: q3=  -0.6911147651e-5_dp
       real(kind=dp), parameter :: q4=   0.7621095161e-6_dp
       real(kind=dp), parameter :: q5=  -0.934945152e-7_dp
       real(kind=dp), parameter :: r1=   57568490574.0_dp
       real(kind=dp), parameter :: r2=  -13362590354.0_dp ! corrected by LCC and ADA 16 june 2004
       real(kind=dp), parameter :: r3=     651619640.7_dp
       real(kind=dp), parameter :: r4=     -11214424.18_dp
       real(kind=dp), parameter :: r5=         77392.33017_dp
       real(kind=dp), parameter :: r6=          -184.9052456_dp
       real(kind=dp), parameter :: s1=   57568490411.0_dp
       real(kind=dp), parameter :: s2=    1029532985.0_dp
       real(kind=dp), parameter :: s3=       9494680.718_dp
       real(kind=dp), parameter :: s4=         59272.64853_dp
       real(kind=dp), parameter :: s5=           267.8532712_dp
       real(kind=dp), parameter :: s6=             1.0_dp

       real(kind=dp)            :: y
       real(kind=cp)            :: ax, xx, z

       if (abs(x) < 1.0e-05) then
          bessj_0=1.0
          return
       end if
       if (abs(x) < 8.0)then
          y=x**2
          bessj_0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y*  &
                  (s4+y*(s5+y*s6)))))
       else
          ax=abs(x)
          z=8.0/ax
          y=z**2
          xx=ax-0.785398164
          bessj_0=sqrt(0.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*  &
                  p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
       end if

       return
    End Function BessJ0

    !!----
    !!---- Elemental Function BessJ1(x) Result(bessj_1)
    !!----    real(kind=sp), intent(in) : x
    !!----    real(kind=sp)             : bessj_1
    !!----
    !!----    Bessel Fuction J1(x)
    !!----
    !!---- Update: February - 2005
    !!
    Elemental Function BessJ1(x) Result(bessj_1)
       !---- Arguments ----!
       real(kind=cp), intent(in) :: x
       real(kind=cp)             :: bessj_1

       !---- Local variales ----!
       real(kind=dp), parameter :: p1= 1.0_dp
       real(kind=dp), parameter :: p2=  0.183105e-2_dp
       real(kind=dp), parameter :: p3= -0.3516396496e-4_dp
       real(kind=dp), parameter :: p4=  0.2457520174e-5_dp
       real(kind=dp), parameter :: p5= -0.240337019e-6_dp
       real(kind=dp), parameter :: q1=  0.04687499995_dp
       real(kind=dp), parameter :: q2= -0.2002690873e-3_dp
       real(kind=dp), parameter :: q3=  0.8449199096e-5_dp
       real(kind=dp), parameter :: q4= -0.88228987e-6_dp
       real(kind=dp), parameter :: q5=  0.105787412e-6_dp
       real(kind=dp), parameter :: r1=  72362614232.0_dp
       real(kind=dp), parameter :: r2=  -7895059235.0_dp
       real(kind=dp), parameter :: r3=    242396853.1_dp
       real(kind=dp), parameter :: r4=     -2972611.439_dp
       real(kind=dp), parameter :: r5=        15704.48260_dp
       real(kind=dp), parameter :: r6=          -30.16036606_dp
       real(kind=dp), parameter :: s1= 144725228442.0_dp
       real(kind=dp), parameter :: s2=   2300535178.0_dp
       real(kind=dp), parameter :: s3=     18583304.74_dp
       real(kind=dp), parameter :: s4=        99447.43394_dp
       real(kind=dp), parameter :: s5=          376.9991397_dp
       real(kind=dp), parameter :: s6=            1.0_dp

       real(kind=dp)            :: y
       real(kind=cp)            :: ax,xx,z

       if (abs(x) < 1.0e-05) then
          bessj_1=0.0
          return
       end if
       if (abs(x) < 8.0)then
          y=x**2
          bessj_1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+  &
                  y*(s4+y*(s5+y*s6)))))
       else
          ax=abs(x)
          z=8.0/ax
          y=z**2
          xx=ax-2.356194491
          bessj_1=sqrt(0.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*p5))))  &
                       -z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))*sign(1.0_cp,x)
       end if

       return
    End Function BessJ1

    !!----
    !!---- Function BessJ(n,x) Result (bessj)
    !!----    real(kind=cp), intent(in) : x
    !!----    real(kind=cp)             : bessj
    !!----
    !!----    Bessel Fuction Jn(x)
    !!----    Returns the Bessel function Jn(x) for any real x and n >= 2.
    !!----
    !!----  Update:  June - 2004
    !!
    Function BessJ(n,x) Result(bessj_n)
       !---- Arguments ----!
       integer,        intent(in)  :: n
       real (kind=cp), intent(in)  :: x
       real (kind=cp)              :: bessj_n

       !---- Local Arguments ----!
       integer,    parameter       :: iacc=40
       integer                     :: j,jsum,m
       real (kind=cp), parameter   :: bigno=1.e10,bigni=1.e-10
       real (kind=cp)              :: ax,bj,bjm,suma,tox
       real (kind=cp), save        :: bjp

       if (n==0) then
          bessj_n=Bessj0(x)
          return
       else if (n==1) then
          bessj_n=Bessj1(x)
          return
       end if

       ax=abs(x)
       if (ax==0.0)then
          bessj_n=0.0
          return

       else if (ax > float(n))then ! Upwards recurrence from J0 and J1.
          tox=2./ax
          bjm=bessj0(ax)
          bj=bessj1(ax)
          do j=1,n-1
             bjp=j*tox*bj-bjm
             bjm=bj
             bj=bjp
          end do
          bessj_n=bj
          return

       else ! Downwards recurrence from an even m here computed.
            !Make IACC larger to increase accuracy.
          tox=2./ax
          m=2*((n+int(sqrt(float(IACC*n))))/2)
          bessj_n=0.
          jsum=0       !jsum will alternate between 0 and 1; when it is 1, we
                       !accumulate in sum the even terms in (5.5.16).
          suma=0.
          bjp=0.
          bj=1.0
          do j=m,1,-1 ! The downward recurrence.
             bjm=j*tox*bj-bjp
             bjp=bj
             bj=bjm
             if (abs(bj)>BIGNO) then ! Renormalize to prevent overflows.
                bj=bj*BIGNI
                bjp=bjp*BIGNI
                bessj_n=bessj_n*BIGNI
                suma=suma*BIGNI
             end if
             if (jsum/=0) suma=suma+bj  ! Accumulate the sum.
             jsum=1-jsum                ! Change 0 to 1 or vice versa.
             if (j==n) bessj_n=bjp      ! the unnormalized answer.
          end do
          suma=2.*suma-bj          ! Compute (5.5.16)
          bessj_n=bessj_n/suma     ! and use it to normalize the answer.
       end if
       if ((x<0.0).and.(mod(n,2)==1)) bessj_n=-bessj_n
       return
    End Function BessJ

    !!----
    !!---- Function ERF(X)
    !!----    real(kind=cp), intent(in) : x
    !!----    real(kind=cp)             : Erf
    !!----
    !!----    Error Function
    !!----    Returns the Error Function
    !!----
    !!----  Update:  January - 2016
    !!----
    Elemental Function Erf(X) Result(Fn_Val)
       !---- Arguments ----!
       real (kind=cp), intent(in) :: x
       real (kind=cp)             :: fn_Val

       !---- Local Variables ----!
       real (kind=cp), parameter :: C    = 0.564189583547756_cp
       real (kind=cp), parameter :: ONE  = 1.0_cp
       real (kind=cp), parameter :: HALF = 0.5_cp
       real (kind=cp), parameter :: ZERO = 0.0_cp

       real (kind=dp), parameter :: A(5) = (/0.771058495001320D-04, -0.133733772997339D-02, 0.323076579225834D-01, &
                                             0.479137145607681D-01,  0.128379167095513D+00 /)
       real (kind=dp), parameter :: B(3) = (/0.301048631703895D-02,  0.538971687740286D-01, 0.375795757275549D+00 /)
       real (kind=dp), parameter :: P(8) = (/-1.36864857382717D-07,  5.64195517478974D-01,  7.21175825088309D+00, &
                                              4.31622272220567D+01,  1.52989285046940D+02,  3.39320816734344D+02, &
                                              4.51918953711873D+02,  3.00459261020162D+02 /)
       real (kind=dp), parameter :: Q(8) = (/ 1.00000000000000D+00,  1.27827273196294D+01,  7.70001529352295D+01, &
                                              2.77585444743988D+02,  6.38980264465631D+02,  9.31354094850610D+02, &
                                              7.90950925327898D+02,  3.00459260956983D+02 /)
       real (kind=dp), parameter :: R(5) = (/ 2.10144126479064D+00,  2.62370141675169D+01,  2.13688200555087D+01, &
                                              4.65807828718470D+00,  2.82094791773523D-01 /)
       real (kind=dp), parameter :: S(4) = (/ 9.41537750555460D+01,  1.87114811799590D+02,  9.90191814623914D+01, &
                                              1.80124575948747D+01 /)

       real (kind=cp) :: ax, bot, t, top, x2

       !> Init
       ax = ABS(x)

       if (ax <= half) then
          t = x*x
          top = ((((a(1)*t + a(2))*t + a(3))*t + a(4))*t + a(5)) + one
          bot = ((b(1)*t + b(2))*t + b(3))*t + one
          fn_val = x*(top/bot)
          return
       end if

       if (ax <= 4.0_cp) then
          top = ((((((p(1)*ax + p(2))*ax + p(3))*ax + p(4))*ax + p(5))*ax  &
              + p(6))*ax + p(7))*ax + p(8)
          bot = ((((((q(1)*ax + q(2))*ax + q(3))*ax + q(4))*ax + q(5))*ax  &
              + q(6))*ax + q(7))*ax + q(8)

          fn_val = half + (half - exp(-x*x)*top/bot)
          if (x < zero) fn_val = -fn_val
          return
       end if

       if (ax < 5.8_cp) then
          x2 = x*x
          t = one / x2
          top = (((r(1)*t + r(2))*t + r(3))*t + r(4))*t + r(5)
          bot = (((s(1)*t + s(2))*t + s(3))*t + s(4))*t + one
          fn_val = (c - top/(x2*bot)) / ax
          fn_val = half + (half - exp(-x2)*fn_val)
          if (x < zero) fn_val = -fn_val
          return
       end if

       fn_val = SIGN(one, x)
       return
    End Function Erf

    !!----
    !!---- Elemental Function Factorial(n) Result(fact)
    !!----    integer, intent(in) : n
    !!----
    !!----    Factorial of N
    !!----
    !!---- Update: February - 2005
    !!
    Elemental Function Factorial(n) Result(fact)
       !---- Argument ----!
       integer, intent(in) :: n
       integer             :: fact

       !---- Local variables ----!
       integer   :: nt, np

       if (n ==0) then
          fact=1
       else
          nt=1
          np=abs(n)
          do
             nt=nt*np
             np=np-1
             if(np == 1) exit
          end do
          fact=nt
       end if

       return
    End Function Factorial

    !!----
    !!---- Elemental Function Negligible(v)
    !!----    complex/real(kind=sp),    intent( in) :: v
    !!----
    !!----    Provides the value .TRUE. if the real/complex
    !!----    number V is less than EPS
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Elemental Function Negligiblec(v)
    !!--++    complex, intent( in) :: v
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate if a complex number is negligible
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Negligiblec(v) Result(Neglig)
       !---- Argument ----!
       complex, intent( in) :: v
       logical              :: Neglig

       Neglig=.false.
       if (abs(v) > epss) return
       Neglig=.true.

       return
    End Function Negligiblec

    !!--++
    !!--++ Elemental Function Negligibler(v)
    !!--++    real(kind=cp), intent( in) :: v
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if a real number is negligible (abs < EPSS)
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Negligibler(v) Result(neglig)
       !---- Argument ----!
       real(kind=cp), intent( in) :: v
       logical                    :: Neglig

       Neglig=.false.
       if (abs(v) > epss) return
       Neglig=.true.

       return
    End Function Negligibler

    !!----
    !!---- Function Pgcd(i,j) Result(mcd)
    !!----    integer, intent(in) :: i
    !!----    integer, intent(in) :: j
    !!----    integer             :: mcd
    !!----
    !!----    Function calculating the maximum common divisor of two integers
    !!----
    !!---- Update: February - 2005
    !!
    Function Pgcd(a,b) Result(mcd)
       !---- Arguments ----!
       integer, intent(in) :: a,b
       integer             :: mcd

       !---- Local variables ----!
       integer  :: u,v,m

       u=max(a,b)
       v=min(a,b)
       m=0
       do
          if (m == 1) exit
          m=mod(u,v)
          u=v
          v=m
       end do
       mcd=u

       return
    End Function Pgcd

    !!----
    !!---- Function Ppcm(i,j) result(mcm)
    !!----    integer, intent(in) :: i
    !!----    integer, intent(in) :: j
    !!----    integer             :: mcm
    !!----
    !!----    Function calculating the minimum common multiple of two integers
    !!----
    !!---- Update: February - 2005
    !!
    Function Ppcm(a,b) result(mcm)
       !---- Arguments ----!
       integer, intent(in) :: a,b
       integer             :: mcm

       !---- Local variables ----!
       integer :: u,v,w,i

       u=max(a,b)
       v=min(a,b)
       mcm=1
       if (v <= 1) then
          mcm=u
          return
       end if
       w=int(sqrt(real(u)))+1
       do i=2,w
          do
             if(.not. ((mod(u,i)==0) .or. (mod(v,i)==0)) ) exit
             mcm=mcm*i
             if (modulo(u,i) == 0) u=u/i
             if (modulo(v,i) == 0) v=v/i
          end do
       end do

       return
    End Function Ppcm

    !!----
    !!---- Function Pythag(a,b) Result (c)
    !!----    real(sp/dp),intent(in):: a,b
    !!----    real(sp/dp)           :: c
    !!--<<
    !!----    Computes c=sqrt(a^2 +b^2 ) without destructive underflow or overflow.
    !!----    Adapted from Numerical Recipes.
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Function Pythag_dp(a,b) Result (c)
    !!--++    real(dp),intent(in):: a,b
    !!--++    real(dp)           :: c
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Computes c=sqrt(a^2 +b^2 ) without destructive underflow or overflow.
    !!--++    Adapted from Numerical Recipes.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Pythag_dp(a,b) Result (c)
       !---- Arguments ----!
       real(kind=dp),intent(in):: a,b
       real(kind=dp)           :: c

       !---- Local variables ----!
       real(kind=dp)           :: absa,absb

       absa=abs(a)
       absb=abs(b)
       if (absa >absb)then
          c=absa*sqrt(1.0_dp+(absb/absa)**2)
       else
          if (absb < tiny(1.0_dp))then
             c=0.0
          else
             c=absb*sqrt(1.0_dp+(absa/absb)**2)
          end if
       end if

       return
    End Function Pythag_dp

    !!--++
    !!--++ Function Pythag_sp(a,b) result (c)
    !!--++    real(sp),intent(in):: a,b
    !!--++    real(sp)           :: c
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Computes c=sqrt(a^2 +b^2 ) without destructive underflow or overflow.
    !!--++    Adapted from Numerical Recipes.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Pythag_sp(a,b) Result (c)
       !---- Arguments ----!
       real(kind=sp),intent(in):: a,b
       real(kind=sp)           :: c

       !---- Local variables ----!
       real(kind=sp)           :: absa,absb

       absa=abs(a)
       absb=abs(b)
       if (absa > absb) then
          c=absa*sqrt(1.0_sp+(absb/absa)**2)
       else
          if (absb < tiny(1.0_sp)) then
             c=0.0
          else
             c=absb*sqrt(1.0_sp+(absa/absb)**2)
          end if
       end if

       return
    End Function Pythag_sp

    !!----
    !!---- Logical Function Co_Linear(A,B,N)
    !!----    complex/integer/real(kind=sp), dimension(:), intent(in)  :: a
    !!----    complex/integer/real(kind=sp), dimension(:), intent(in)  :: b
    !!----    integer,                                     intent(in)  :: n
    !!----
    !!----    Provides the value .TRUE. if the vectors A and B are co-linear
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Logical Function Co_Linear_C(A, B, N)
    !!--++    complex, dimension(:), intent(in)  :: a
    !!--++    complex, dimension(:), intent(in)  :: b
    !!--++    integer,               intent(in)  :: n
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two complex vectors are co-linear
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Co_linear_C(a,b,n) Result(co_linear)
       !---- Argument ----!
       complex, dimension(:), intent(in) :: a,b
       integer,               intent(in) :: n
       logical                           :: co_linear

       !---- Local variables ----!
       integer :: i,ia,ib
       complex :: c

       co_linear=.true.
       do i=1,n
          if (abs(a(i)) > epss) then
             ia=i
             exit
          end if
       end do
       do i=1,n
          if (abs(b(i)) > epss) then
             ib=i
             exit
          end if
       end do
       if (ia /= ib) then
          co_linear=.false.
          return
       else
          c=a(ia)/b(ib)
          do i=1,n
             if (abs(a(i)-c*b(i)) > epss) then
                co_linear=.false.
                return
             end if
          end do
       end if

       return
    End Function Co_linear_C

    !!--++
    !!--++ Logical Function Co_Linear_I(A, B, N)
    !!--++    integer, dimension(:), intent(in)  :: a
    !!--++    integer, dimension(:), intent(in)  :: b
    !!--++    integer,               intent(in)  :: n
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two integer vectors are co-linear
    !!--++
    !!--++ Update: October - 2008
    !!
    Function Co_linear_I(a,b,n) Result(co_linear)
       !---- Argument ----!
       integer, dimension(:), intent(in) :: a,b
       integer,               intent(in) :: n
       logical                           :: co_linear

       !---- Local variables ----!
       integer       :: i,ia,ib
       real(kind=cp) :: c

       co_linear=.true.
       do i=1,n
          if (abs(a(i)) > 0) then
             ia=i
             exit
          end if
       end do
       do i=1,n
          if (abs(b(i)) > 0) then
             ib=i
             exit
          end if
       end do
       if (ia /= ib) then
          co_linear=.false.
          return
       else
          c=real(a(ia))/real(b(ib))
          do i=1,n
             if (abs( real(a(i))-c*real(b(i)) ) > epss) then
                co_linear=.false.
                return
             end if
          end do
       end if

       return
    End Function Co_linear_I

    !!--++
    !!--++ Logical Function Co_Linear_R(A, B, N)
    !!--++    real(kind=cp), dimension(:), intent(in)  :: a
    !!--++    real(kind=cp), dimension(:), intent(in)  :: b
    !!--++    integer,                     intent(in)  :: n
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two real vectors are co-linear
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Co_linear_R(a,b,n) Result(co_linear)
       !---- Argument ----!
       real(kind=cp), dimension(:), intent(in) :: a,b
       integer,                     intent(in) :: n
       logical                                 :: co_linear

       !---- Local variables ----!
       integer       :: i,ia,ib
       real(kind=cp) :: c

       co_linear=.true.
       do i=1,n
          if (abs(a(i)) > epss) then
             ia=i
             exit
          end if
       end do
       do i=1,n
          if (abs(b(i)) > epss) then
             ib=i
             exit
          end if
       end do
       if (ia /= ib) then
          co_linear=.false.
          return
       else
          c=a(ia)/b(ib)
          do i=1,n
             if (abs(a(i)-c*b(i)) > epss) then
                co_linear=.false.
                return
             end if
          end do
       end if

       return
    End Function Co_linear_R

    !!----
    !!---- Function Co_Prime(v,imax) result(cop)
    !!----   integer, dimension(:), intent(in) :: v
    !!----   integer,  optional,    intent(in) :: imax !Maximun prime number to be tested
    !!----   Logical                           :: cop
    !!----
    !!---- Provides the value .TRUE. if the array V contains co-prime
    !!---- integers: there is no common divisor for all the integers.
    !!---- Only the first 1000 prime numbers are stored in the module array "primes"
    !!---- imax is the maximum prime number to be tested. It is calculated if not given.
    !!----
    !!---- Created: January - 2011
    !!---- Updated: February - 2012  (JRC) (imax argument made optional, really not needed)
    !!
    Function Co_Prime(v,imax) result(cop)
      integer, dimension(:), intent(in) :: v
      integer, optional,     intent(in) :: imax
      Logical                           :: cop
      !---- Local variables ----!
      integer :: i,j,im,k,dimv,imaxv,maxv

      cop=.true.
      maxv=maxval(abs(v))
      if(present(imax)) then
        imaxv=imax
      else
        imaxv=maxv
      end if
      !---- If the maximum value of the indices is 1 they are coprimes
      if (maxv == 1) return
      if (maxv == 0) then
         cop=.false.
         return
      end if
      !---- Search the maximum prime number to be tested
      if (imaxv > 7919) then
        im=1000
      else
        do i=1,1000
           if(imaxv > primes(i)) cycle
           im=i
           exit
        end do
      end if

      !---- Indices greater than 1
      dimv=size(v)
      do_p: do i=1,im
         k=primes(i)
         do j=1,dimv
            if( mod(v(j),k) /= 0) cycle do_p
         end do
         cop=.false.
         exit
      end do do_p

      return
    End Function Co_Prime

    !!----
    !!---- Logical Function Equal_Matrix(A,B,N)
    !!----    integer/real(kind=cp), dimension(:,:), intent(in)  :: a,b
    !!----    integer,                               intent(in)  :: n
    !!----
    !!----    Provides the value .TRUE. if the array A is equal to array B
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Logical Function Equal_Matrix_I(A, B, N)
    !!--++    integer, dimension(:,:), intent(in)  :: a
    !!--++    integer, dimension(:,:), intent(in)  :: b
    !!--++    integer,                 intent(in)  :: n
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two integer arrays are equal in NxN
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Equal_Matrix_I(a,b,n) result(info)
       !---- Argument ----!
       integer, dimension(:,:), intent(in) :: a,b
       integer                , intent(in) :: n
       logical                             :: info

       !---- Local variables ----!
       integer :: i,j

       info=.false.
       do i=1,n
          do j=1,n
             if (a(i,j) /= b(i,j)) return
          end do
       end do
       info=.true.

       return
    End Function Equal_Matrix_I

    !!--++
    !!--++ Logical Function Equal_Matrix_R(A, B, N)
    !!--++    real(kind=sp), dimension(:,:), intent(in)  :: a
    !!--++    real(kind=sp), dimension(:,:), intent(in)  :: b
    !!--++    integer,                       intent(in)  :: n
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two integer arrays are equal in NxN
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Equal_Matrix_R(a,b,n) result(info)
       !---- Argument ----!
       real(kind=cp), dimension(:,:)   , intent(in) :: a,b
       integer,                          intent(in) :: n
       logical                                      :: info

       !---- Local variables ----!
       integer :: i,j

       info=.false.
       do i=1,n
          do j=1,n
             if (abs(a(i,j) - b(i,j)) > epss ) return
          end do
       end do
       info=.true.

       return
    End Function Equal_Matrix_R

    !!----
    !!---- Logical Function Equal_Vector(A,B,N)
    !!----    integer/real(kind=sp), dimension(:),   intent(in)  :: a,b
    !!----    integer,                               intent(in)  :: n
    !!----
    !!----    Provides the value .TRUE. if the vector A is equal to vector B
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Logical Function Equal_Vector_I(A, B, N)
    !!--++    integer, dimension(:), intent(in)  :: a
    !!--++    integer, dimension(:), intent(in)  :: b
    !!--++    integer,               intent(in)  :: n
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two integer vectors are equal in N
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Equal_Vector_I(a,b,n) result(info)
       !---- Argument ----!
       integer, dimension(:),   intent(in) :: a,b
       integer                , intent(in) :: n
       logical                             :: info

       !---- Local variables ----!
       integer :: i

       info=.false.
       do i=1,n
          if (a(i) /= b(i)) return
       end do
       info=.true.

       return
    End Function Equal_Vector_I

    !!--++
    !!--++ Logical Function Equal_Vector_R(A, B, N)
    !!--++    real(kind=sp), dimension(:), intent(in)  :: a
    !!--++    real(kind=sp), dimension(:), intent(in)  :: b
    !!--++    integer,                     intent(in)  :: n
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two real(kind=sp) vectors are equal in N
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Equal_Vector_R(a,b,n) result(info)
       !---- Argument ----!
       real(kind=cp), dimension(:)   ,   intent(in) :: a,b
       integer,                          intent(in) :: n
       logical                                      :: info

       !---- Local variables ----!
       integer :: i

       info=.false.
       do i=1,n
          if (abs(a(i) - b(i)) > epss ) return
       end do
       info=.true.

       return
    End Function Equal_Vector_R

    !!----
    !!----  Function Euclidean_Norm(n,x) Result(Fn_Val)
    !!----    Integer,                      Intent(In)  :: n
    !!----    Real (Kind=cp), Dimension(:), Intent(In)  :: x
    !!----    Real (Kind=cp)                            :: Fn_Val
    !!----
    !!----  This function calculates safely the Euclidean norm of a vector.
    !!----  Intermediate overflows are avoided using this function. The original
    !!----  name "enorm" from MINPACK has been changed and the subroutine has
    !!----  been translated to Fortran 90.
    !!----
    !!----
    !!--..  Original documentation (from MINPACK):
    !!--..
    !!--..  Function enorm
    !!--..
    !!--..  Given an n-vector x, this function calculates the euclidean norm of x.
    !!--..
    !!--..  The euclidean norm is computed by accumulating the sum of squares in
    !!--..  three different sums.  The sums of squares for the small and large
    !!--..  components are scaled so that no overflows occur.  Non-destructive
    !!--..  underflows are permitted.  Underflows and overflows do not occur in the
    !!--..  computation of the unscaled sum of squares for the intermediate
    !!--..  components.  The definitions of small, intermediate and large components
    !!--..  depend on two constants, rdwarf and rgiant.  The main restrictions on
    !!--..  these constants are that rdwarf**2 not underflow and rgiant**2 not
    !!--..  overflow.  The constants given here are suitable for every known computer.
    !!--..
    !!--..  The function statement is
    !!--..
    !!--..    REAL (kind=cp) function enorm(n,x)
    !!--..
    !!--..  where
    !!--..
    !!--..    n is a positive integer input variable.
    !!--..
    !!--..    x is an input array of length n.
    !!--..
    !!--..  Subprograms called
    !!--..
    !!--..    Fortran-supplied ... ABS,SQRT
    !!--..
    !!--..  Argonne National Laboratory. MINPACK project. march 1980.
    !!--..  Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
    !!----
    !!----  Update: August - 2009
    !!----
    Function Euclidean_Norm(n,x) Result(Fn_Val)
       !---- Arguments ----!
       Integer,                      Intent(In)  :: n
       Real (Kind=cp), Dimension(:), Intent(In)  :: x
       Real (Kind=cp)                            :: Fn_Val

       !--- Local Variables ---!
       Integer                   :: i
       Real (Kind=cp)            :: agiant, floatn, s1, s2, s3, xabs, x1max, x3max
       Real (Kind=cp), Parameter :: one = 1.0_cp, zero = 0.0_cp, rdwarf = 3.834e-20_cp,  &
                                   rgiant = 1.304e+19_cp

       s1 = zero
       s2 = zero
       s3 = zero
       x1max = zero
       x3max = zero
       floatn = n
       agiant = rgiant/floatn
       do i = 1, n
          xabs = Abs(x(i))
          if (.Not. (xabs > rdwarf .AND. xabs < agiant)) then
             ! sum for large components.
             if (xabs > rdwarf) then
                if (xabs > x1max) then
                   s1 = one + s1*(x1max/xabs)**2
                   x1max = xabs
                   cycle
                end if
                s1 = s1 + (xabs/x1max)**2
                cycle
             End If

             ! sum for small components.
             If (xabs > x3max) Then
                s3 = one + s3*(x3max/xabs)**2
                x3max = xabs
                Cycle
             End If

             If (xabs /= zero) s3 = s3 + (xabs/x3max)**2
             Cycle
          End if

          !  sum for intermediate components.
          s2 = s2 + xabs**2
       End Do

       ! calculation of norm.
       If (s1 /= zero) Then
          Fn_Val = x1max*Sqrt(s1 + (s2/x1max)/x1max)
          Return
       End If

       If (s2 /= zero) Then
          If (s2 >= x3max) Fn_Val = Sqrt(s2*(one + (x3max/s2)*(x3max*s3)))
          If (s2 < x3max) Fn_Val = Sqrt(x3max*((s2/x3max) + (x3max*s3)))
          Return
       End If

       Fn_Val = x3max*Sqrt(s3)

       Return
    End Function Euclidean_Norm

    !!----
    !!---- Function Imaxloc(arr) Result(mav)
    !!----  real(kind=sp)/integer, dimension(:), intent(in) :: arr
    !!----  integer                                         :: mav
    !!----
    !!----   Index of maxloc on an array
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Function Imaxloc_I(arr) Result(mav)
    !!--++  integer, dimension(:), intent(in) :: arr
    !!--++  integer                           :: mav
    !!--++
    !!--++   Index of maxloc on an array (from Numerical Recipes)
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Imaxloc_I(iarr) Result(mav)
       !---- Arguments ----!
       integer, dimension(:), intent(in) :: iarr
       integer                           :: mav

       !---- Local variables ----!
       integer, dimension(1) :: imax

       imax=maxloc(iarr(:))
       mav=imax(1)

       return
    End Function Imaxloc_I

    !!--++
    !!--++ Function Imaxloc_R(arr) Result(mav)
    !!--++  real(kind=sp), dimension(:), intent(in) :: arr
    !!--++  integer                                 :: mav
    !!--++
    !!--++   Index of maxloc on an array (from Numerical Recipes)
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Imaxloc_R(arr) Result(mav)
       !---- Arguments ----!
       real(kind=cp), dimension(:), intent(in) :: arr
       integer                                 :: mav

       !---- Local variables ----!
       integer, dimension(1) :: imax

       imax=maxloc(arr(:))
       mav=imax(1)

       return
    End Function Imaxloc_R

    !!----
    !!---- Function Iminloc(arr)  Result(miv)
    !!----  real(kind=sp)/integer, dimension(:), intent(in) :: arr
    !!----  integer                                         :: miv
    !!----
    !!----   Index of minloc on an array  (from Numerical Recipes)
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Function Iminloc_I(arr)  Result(miv)
    !!--++  integer, dimension(:), intent(in) :: arr
    !!--++  integer                           :: miv
    !!--++
    !!--++   Index of minloc on an array (from Numerical Recipes)
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Iminloc_I(arr)  Result(miv)
       !---- Arguments ----!
       integer, dimension(:), intent(in) :: arr
       integer                           :: miv

       !---- Local variables ----!
       integer, dimension(1) :: imin

       imin=minloc(arr(:))
       miv=imin(1)

       return
    End Function Iminloc_I

    !!--++
    !!--++ Function Iminloc_R(arr)  Result(miv)
    !!--++  real(kind=sp), dimension(:), intent(in) :: arr
    !!--++  integer                                 :: miv
    !!--++
    !!--++   Index of minloc on an array (from Numerical Recipes)
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Iminloc_R(arr)  Result(miv)
       !---- Arguments ----!
       real(kind=cp), dimension(:), intent(in) :: arr
       integer                                 :: miv

       !---- Local variables ----!
       integer, dimension(1) :: imin

       imin=minloc(arr(:))
       miv=imin(1)

       return
    End Function Iminloc_R

    !!----
    !!---- Function in_limits(n,limits,vect) result(ok)
    !!----   integer,                      intent(in) :: n
    !!----   integer/real, dimension(:,:), intent(in) :: limits   ! Normally (2,n)
    !!----   integer/real, dimension(n),   intent(in) :: vect
    !!----   logical                                  :: ok
    !!----
    !!----   Logical function that is true if all the components of the vector vect
    !!----   are within the limits:   limits(1,i)  <= vect(i) <=  limits(2,i), for all i.
    !!----
    !!----   Updated: March - 2013
    !!
    !!--++
    !!--++ Function in_limits_int(n,limits,vect) result(ok)
    !!--++   integer,                 intent(in) :: n
    !!--++   integer, dimension(:,:), intent(in) :: limits   ! Normally (2,n)
    !!--++   integer, dimension(n),   intent(in) :: vect
    !!--++   logical                              :: ok
    !!--++
    !!--++   Logical function that is true if all the components of the vector vect
    !!--++   are within the limits:   limits(1,i)  <= vect(i) <=  limits(2,i), for all i.
    !!--++
    !!--++   Updated: March - 2013
    !!
    Function in_limits_int(n,limits,vect) result(ok)
      integer,                 intent(in) :: n
      integer, dimension(:,:), intent(in) :: limits   ! Normally (2,n)
      integer, dimension(n),   intent(in) :: vect
      logical :: ok
      integer :: i
      ok=.true.
      do i=1,n
        if(vect(i) >= limits(1,i) .and. vect(i) <= limits(2,i)) cycle
        ok=.false.
        exit
      end do
      return
    End Function in_limits_int

    !!--++
    !!--++ Function in_limits_dp(n,limits,vect) result(ok)
    !!--++   integer,                       intent(in) :: n
    !!--++   real(kind=dp), dimension(:,:), intent(in) :: limits   ! Normally (2,n)
    !!--++   real(kind=dp), dimension(n),   intent(in) :: vect
    !!--++   logical                                   :: ok
    !!--++
    !!--++   Logical function that is true if all the components of the vector vect
    !!--++   are within the limits:   limits(1,i)  <= vect(i) <=  limits(2,i), for all i.
    !!--++
    !!--++   Updated: March - 2013
    !!
    Function in_limits_dp(n,limits,vect) result(ok)
      integer,                       intent(in) :: n
      real(kind=dp), dimension(:,:), intent(in) :: limits   ! Normally (2,n)
      real(kind=dp), dimension(n),   intent(in) :: vect
      logical :: ok
      integer :: i
      ok=.true.
      do i=1,n
        if(vect(i) >= limits(1,i) .and. vect(i) <= limits(2,i)) cycle
        ok=.false.
        exit
      end do
      return
    End Function in_limits_dp

    !!--++
    !!--++ Function in_limits_sp(n,limits,vect) result(ok)
    !!--++   integer,                       intent(in) :: n
    !!--++   real(kind=sp), dimension(:,:), intent(in) :: limits   ! Normally (2,n)
    !!--++   real(kind=sp), dimension(n),   intent(in) :: vect
    !!--++   logical                                   :: ok
    !!--++
    !!--++   Logical function that is true if all the components of the vector vect
    !!--++   are within the limits:   limits(1,i)  <= vect(i) <=  limits(2,i), for all i.
    !!--++
    !!--++   Updated: March - 2013
    !!
    Function in_limits_sp(n,limits,vect) result(ok)
      integer,                       intent(in) :: n
      real(kind=sp), dimension(:,:), intent(in) :: limits   ! Normally (2,n)
      real(kind=sp), dimension(n),   intent(in) :: vect
      logical :: ok
      integer :: i
      ok=.true.
      do i=1,n
        if(vect(i) >= limits(1,i) .and. vect(i) <= limits(2,i)) cycle
        ok=.false.
        exit
      end do
      return
    End Function in_limits_sp

    !!----
    !!---- Function Locate(xx, n, x) Result(j)
    !!----     or
    !!---- Function Locate(xx,x) Result(j)
    !!----
    !!----    integer/real(kind=sp), dimension(n),intent(in)  :: xx
    !!----    integer ,                           intent(in)  :: n
    !!----    integer/real(kind=sp),              intent(in)  :: x
    !!----    integer ,                           intent(out) :: j
    !!----
    !!----    Function for locating the index J of an array XX(N)
    !!----    satisfying:
    !!--<<
    !!----               XX(J) <= X < XX(J+1)
    !!-->>
    !!----
    !!---- Update: June - 2011
    !!

    !!--++
    !!--++ Function Locate_I(xx, n, x) Result(j)
    !!--++    integer, dimension(:),intent(in)  :: xx
    !!--++    integer ,             intent(in)  :: n
    !!--++    integer,              intent(in)  :: x
    !!--++    integer ,             intent(out) :: j
    !!--++
    !!--++    Subroutine for locating the index J of an array XX(N)
    !!--++    satisfying:
    !!--++
    !!--++               XX(J) <= X < XX(J+1)
    !!--++
    !!--++
    !!--++ Update: June - 2011
    !!
    Function Locate_I(xx,n,x) Result(j)
       !---- Argument ----!
       integer, dimension(:), intent(in):: xx
       integer ,              intent(in):: n
       integer,               intent(in):: x
       integer                          :: j

       !---- Local Variables ----!
       integer :: jl, ju, jm

       if(x <= xx(1)) then
         j=1
         return
       end if
       if(x >= xx(n)) then
         j=n
         return
       end if
       jl=0
       ju=n+1
       do
          if(ju-jl <= 1) exit
          jm=(ju+jl)/2
          if ((xx(n) > xx(1)) .eqv. (x > xx(jm))) then
             jl=jm
          else
             ju=jm
          end if
       end do
       j=jl

       return
    End Function Locate_I
    !!--++
    !!--++ Function Locate_Ib(xx, x) Result(j)
    !!--++    integer, dimension(:),intent(in)  :: xx
    !!--++    integer,              intent(in)  :: x
    !!--++    integer ,             intent(out) :: j
    !!--++
    !!--++    Subroutine for locating the index J of an array XX(:)
    !!--++    satisfying:
    !!--++
    !!--++               XX(J) <= X < XX(J+1)
    !!--++
    !!--++
    !!--++ Update: June - 2011
    Function Locate_Ib(xx,x) Result(j)
       !---- Argument ----!
       integer, dimension(:), intent(in):: xx
       integer,               intent(in):: x
       integer                          :: j

       !---- Local Variables ----!
       integer :: jl, ju, jm, i1,i2
       integer, dimension(1) :: mi

       mi=lbound(xx)
       i1=mi(1)
       mi=ubound(xx)
       i2=mi(1)

       if(x <= xx(i1)) then
         j=i1
         return
       end if
       if(x >= xx(i2)) then
         j=i2
         return
       end if
       jl=i1-1
       ju=i2+1
       do
          if(ju-jl <= 1) exit
          jm=(ju+jl)/2
          if ((xx(i2) > xx(i1)) .eqv. (x > xx(jm))) then
             jl=jm
          else
             ju=jm
          end if
       end do
       j=jl
       return
    End Function Locate_Ib
    !!--++
    !!--++ Function Locate_R(xx, n, x) Result(j)
    !!--++    real(kind=cp), dimension(:),intent(in)  :: xx
    !!--++    integer ,                   intent(in)  :: n
    !!--++    real(kind=cp),              intent(in)  :: x
    !!--++    integer ,                   intent(out) :: j
    !!--++
    !!--++    Function for locating the index J of an array XX(N)
    !!--++    satisfying:
    !!--++
    !!--++               XX(J) <= X < XX(J+1)
    !!--++
    !!--++
    !!--++ Update: June - 2011
    !!
    Function Locate_R(xx,n,x) Result(j)
       !---- Argument ----!
       real(kind=cp), dimension(:), intent(in):: xx
       integer ,                    intent(in):: n
       real(kind=cp),               intent(in):: x
       integer                                :: j

       !---- Local Variables ----!
       integer :: jl, ju, jm

       if(x <= xx(1)) then
         j=1
         return
       end if
       if(x >= xx(n)) then
         j=n
         return
       end if
       jl=0
       ju=n+1
       do
          if(ju-jl <= 1) exit
          jm=(ju+jl)/2
          if ((xx(n) > xx(1)) .eqv. (x > xx(jm))) then
             jl=jm
          else
             ju=jm
          end if
       end do
       j=jl

       return
    End Function Locate_R

    !!--++
    !!--++ Function Locate_Rb(xx, x) Result(j)
    !!--++    real(kind=cp), dimension(:),intent(in)  :: xx
    !!--++    real(kind=cp),              intent(in)  :: x
    !!--++    integer ,                   intent(out) :: j
    !!--++
    !!--++    Function for locating the index J of an array XX(:)
    !!--++    satisfying:
    !!--++
    !!--++               XX(J) <= X < XX(J+1)
    !!--++
    !!--++
    !!--++ Update: June - 2011
    !!
    Function Locate_Rb(xx,x) Result(j)
       !---- Argument ----!
       real(kind=cp), dimension(:), intent(in):: xx
       real(kind=cp),               intent(in):: x
       integer                                :: j

       !---- Local Variables ----!
       integer :: jl, ju, jm, i1,i2
       integer, dimension(1) :: mi

       mi=lbound(xx)
       i1=mi(1)
       mi=ubound(xx)
       i2=mi(1)

       if(x <= xx(i1)) then
         j=i1
         return
       end if
       if(x >= xx(i2)) then
         j=i2
         return
       end if
       jl=i1-1
       ju=i2+1
       do
          if(ju-jl <= 1) exit
          jm=(ju+jl)/2
          if ((xx(i2) > xx(i1)) .eqv. (x > xx(jm))) then
             jl=jm
          else
             ju=jm
          end if
       end do
       j=jl

       return
    End Function Locate_Rb

    !!----
    !!---- Function Lower_Triangular_I(A,n) Result (T)
    !!----   integer, dimension(:,:), intent(in) :: A
    !!----   integer,                 intent(in) :: n
    !!----   integer, dimension(n,n)             :: T
    !!----
    !!----   Updated: October - 2014
    !!----
    Function Lower_Triangular_I(A,n) Result (T)
       !---- Argument ----!
       integer, dimension(:,:), intent(in) :: A
       integer,                 intent(in) :: n
       integer, dimension(n,n)             :: T
       integer :: i,j,p,q,m
       m=n
       p=size(A(:,1)); q=size(A(1,:))
       if(n > p .or. n > q) m=min(p,q)
       T=0
       do j=1,m
         do i=j,m
           T(i,j)=A(i,j)
         end do
       end do
    End Function  Lower_Triangular_I

    !!----
    !!---- Function Lower_Triangular_R(A,n) Result (T)
    !!----   real(kind=cp), dimension(:,:), intent(in) :: A
    !!----   integer,                       intent(in) :: n
    !!----   real(kind=cp), dimension(n,n)             :: T
    !!----
    !!----   Updated: October - 2014
    !!----
    Function Lower_Triangular_R(A,n) Result (T)
       !---- Argument ----!
       real(kind=cp), dimension(:,:), intent(in) :: A
       integer,                       intent(in) :: n
       real(kind=cp), dimension(n,n)             :: T
       integer :: i,j,p,q,m
       m=n
       p=size(A(:,1)); q=size(A(1,:))
       if(n > p .or. n > q) m=min(p,q)
       T=0
       do j=1,m
         do i=j,m
           T(i,j)=A(i,j)
         end do
       end do
    End Function  Lower_Triangular_R

    !!---- Function Modulo_Lat(U)
    !!----    real(kind=cp), dimension(:), intent(in) :: u
    !!----
    !!----    Reduces a real vector to another with components in
    !!----    the interval [0,1)
    !!----
    !!---- Updated: February - 2005
    !!
    Function Modulo_Lat(u) result(v)
       !---- Argument ----!
       real(kind=cp), dimension(:), intent( in) :: u
       real(kind=cp), dimension(1:size(u))      :: v

       v=mod(u+10.0_cp,1.0_cp)

       return
    End Function  Modulo_Lat

    !!----
    !!---- Function Norm(X,G) Result(R)
    !!----    real(kind=cp)/integer, dimension(:),   intent(in) :: x
    !!----    real(kind=cp),         dimension(:,:), intent(in) :: g
    !!----
    !!----    Calculate the Norm of a vector
    !!----
    !!---- Update: April - 2009
    !!

    !!--++
    !!--++ Function Norm_I(X,G) Result(R)
    !!--++    integer,      dimension(:),   intent(in) :: x
    !!--++    real(kind=cp),dimension(:,:), intent(in) :: g
    !!--++
    !!--++    Calculate the Norm of a vector
    !!--++
    !!--++ Update: April - 2009
    !!
    Function Norm_I(X,G) Result(R)
       !---- Arguments ----!
       integer,       dimension(:),   intent(in) :: x
       real(kind=cp), dimension(:,:), intent(in) :: g
       real(kind=cp)                             :: r

       if (size(x)*size(x) /= size(g)) then
          r=tiny(0.0)
       else
          r=sqrt(dot_product(real(x), matmul(g,real(x))))
       end if

       return
    End Function Norm_I

    !!--++
    !!--++ Function Norm_R(X,G) Result(R)
    !!--++    real(kind=cp),dimension(:),   intent(in) :: x
    !!--++    real(kind=cp),dimension(:,:), intent(in) :: g
    !!--++
    !!--++    Calculate the Norm of a vector
    !!--++
    !!--++ Update: April - 2009
    !!
    Function Norm_R(X,G) Result(R)
       !---- Arguments ----!
       real(kind=cp), dimension(:),   intent(in) :: x
       real(kind=cp), dimension(:,:), intent(in) :: g
       real(kind=cp)                             :: r

       if (size(x)*size(x) /= size(g)) then
          r=tiny(0.0)
       else
          r=sqrt(dot_product(x, matmul(g,x)))
       end if

       return
    End Function Norm_R


    !!----
    !!---- Function Outerprod(a,b) Result(c)
    !!----    real(sp/dp),dimension(:),intent(in)    :: a,b
    !!----    real(sp/dp),dimension(size(a),size(b)) :: c
    !!----
    !!----    Computes the outer product (tensorial product) of two
    !!----    vectors to give a tensor (matrix) as the result:
    !!--<<
    !!----                   c(i,j) = a(i)*b(j).
    !!-->>
    !!--..    It uses the intrinsic Fortran 90 function SPREAD.
    !!--..    Function adapted from Numerical Recipes.
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Function Outerprod_dp(a,b) Result(c)
    !!--++    real(dp),dimension(:),intent(in)    :: a,b
    !!--++    real(dp),dimension(size(a),size(b)) :: c
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Computes the outer product (tensorial product) of two
    !!--++    vectors to give a tensor (matrix) as the result:
    !!--++                   c(i,j) = a(i)*b(j).
    !!--++
    !!--++    It uses the intrinsic Fortran 90 function SPREAD.
    !!--++    Taken from Numerical Recipes.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Outerprod_dp(a,b)  Result(c)
       !---- Arguments ----!
       real(kind=dp),dimension(:),intent(in)    :: a,b
       real(kind=dp),dimension(size(a),size(b)) :: c

       c =spread(a,dim=2,ncopies=size(b))*spread(b,dim=1,ncopies=size(a))

       return
    End Function Outerprod_dp

    !!--++
    !!--++ Function Outerprod_sp(a,b) Result(c)
    !!--++    real(sp),dimension(:),intent(in)    :: a,b
    !!--++    real(sp),dimension(size(a),size(b)) :: c
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Computes the outer product (tensorial product) of two
    !!--++    vectors to give a tensor (matrix) as the result:
    !!--++                   c(i,j) = a(i)*b(j).
    !!--++
    !!--++    It uses the intrinsic Fortran 90 function SPREAD.
    !!--++    Taken from Numerical Recipes.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Outerprod_sp(a,b)  Result(c)
       !---- Arguments ----!
       real(kind=sp),dimension(:),intent(in)    :: a,b
       real(kind=sp),dimension(size(a),size(b)) :: c

       c =spread(a,dim=2,ncopies=size(b))*spread(b,dim=1,ncopies=size(a))

       return
    End Function Outerprod_sp

    !!----
    !!---- Function Scalar(X,Y,G) Result(R)
    !!----    integer/real(kind=cp), dimension(:),   intent(in) :: x
    !!----    integer/real(kind=cp), dimension(:),   intent(in) :: y
    !!----    real(kind=cp),         dimension(:,:), intent(in) :: g
    !!----
    !!----    Scalar Product including metrics
    !!----
    !!---- Update: April - 2009
    !!

    !!--++
    !!--++ Function Scalar_R(X,Y,G) Result(R)
    !!--++    integer, dimension(:),   intent(in) :: x
    !!--++    integer, dimension(:),   intent(in) :: y
    !!--++    real(kind=cp), dimension(:,:), intent(in) :: g
    !!--++
    !!--++    Scalar Product including metrics
    !!--++
    !!--++ Update: April - 2009
    !!
    Function Scalar_I(X,Y,G) Result(R)
       !---- Arguments ----!
       integer, dimension(:),   intent(in) :: x
       integer, dimension(:),   intent(in) :: y
       real(kind=cp), dimension(:,:), intent(in) :: g
       real(kind=cp)                             :: r

       if (size(x)/= size(y) .or. size(x)*size(x) /= size(g)) then
          r=tiny(0.0)
       else
          r=dot_product(real(x), matmul(g,real(y)))
       end if

       return
    End Function Scalar_I

    !!--++
    !!--++ Function Scalar_R(X,Y,G) Result(R)
    !!--++    real(kind=cp), dimension(:),   intent(in) :: x
    !!--++    real(kind=cp), dimension(:),   intent(in) :: y
    !!--++    real(kind=cp), dimension(:,:), intent(in) :: g
    !!--++
    !!--++    Scalar Product including metrics
    !!--++
    !!--++ Update: April - 2009
    !!
    Function Scalar_R(X,Y,G) Result(R)
       !---- Arguments ----!
       real(kind=cp), dimension(:),   intent(in) :: x
       real(kind=cp), dimension(:),   intent(in) :: y
       real(kind=cp), dimension(:,:), intent(in) :: g
       real(kind=cp)                             :: r

       if (size(x)/= size(y) .or. size(x)*size(x) /= size(g)) then
          r=tiny(0.0)
       else
          r=dot_product(x, matmul(g,y))
       end if

       return
    End Function Scalar_R

    !!----
    !!---- Function Trace(A)
    !!----    complex/integer/real(kind=cp), dimension(:,:), intent(in)  :: a
    !!----
    !!----    Provides the trace of a complex/real or integer matrix
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Function Trace_C(A)
    !!--++    complex, dimension(:,:), intent(in)  :: a
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the trace of a complex nxn array
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Trace_C(a) Result(b)
       !---- Argument ----!
       complex, dimension(:,:), intent(in) :: a
       complex                             :: b

       !---- Local variables ----!
       integer :: i,imax

       b=(0.0,0.0)
       imax=min(size(a,1),size(a,2))
       do i=1,imax
          b=b+a(i,i)
       end do

       return
    End Function Trace_C

    !!--++
    !!--++ Function Trace_I(A)
    !!--++    integer, dimension(:,:), intent(in)  :: a
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the trace of an integer 3x3 array
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Trace_I(a) Result(b)
       !---- Argument ----!
       integer, dimension(:,:), intent(in) :: a
       integer                             :: b

       !---- Local variables ----!
       integer :: i,imax

       b=0
       imax=min(size(a,1),size(a,2))
       do i=1,imax
          b=b+a(i,i)
       end do

       return
    End Function Trace_I

    !!--++
    !!--++ Function Trace_R(A)
    !!--++    real(kind=cp), dimension(:,:), intent(in)  :: a
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the trace of a real 3x3 array
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Trace_R(a) Result(b)
       !---- Argument ----!
       real(kind=cp), dimension(:,:), intent(in) :: a
       real(kind=cp)                             :: b

       !---- Local variables ----!
       integer :: i,imax

       b=0.0
       imax=min(size(a,1),size(a,2))
       do i=1,imax
          b=b+a(i,i)
       end do

       return
    End Function Trace_R


    !!----
    !!---- Function Upper_Triangular_I(A,n) Result (T)
    !!----   integer, dimension(:,:), intent(in) :: A
    !!----   integer,                 intent(in) :: n
    !!----   integer, dimension(n,n)             :: T
    !!----
    !!----   Updated: October - 2014
    !!----
    Function Upper_Triangular_I(A,n) Result (T)
       !---- Argument ----!
       integer, dimension(:,:), intent(in) :: A
       integer,                 intent(in) :: n
       integer, dimension(n,n)             :: T
       integer :: i,j,p,q,m
       m=n
       p=size(A(:,1)); q=size(A(1,:))
       if(n > p .or. n > q) m=min(p,q)
       T=0
       do j=1,m
         do i=1,j
           T(i,j)=A(i,j)
         end do
       end do
    End Function  Upper_Triangular_I

    !!----
    !!---- Function Upper_Triangular_R(A,n) Result (T)
    !!----   real(kind=cp), dimension(:,:), intent(in) :: A
    !!----   integer,                       intent(in) :: n
    !!----   real(kind=cp), dimension(n,n)             :: T
    !!----
    !!----   Updated: October - 2014
    !!----
    Function Upper_Triangular_R(A,n) Result (T)
       !---- Argument ----!
       real(kind=cp), dimension(:,:), intent(in) :: A
       integer,                       intent(in) :: n
       real(kind=cp), dimension(n,n)             :: T
       integer :: i,j,p,q,m
       m=n
       p=size(A(:,1)); q=size(A(1,:))
       if(n > p .or. n > q) m=min(p,q)
       T=0
       do j=1,m
         do i=1,j
           T(i,j)=A(i,j)
         end do
       end do
    End Function  Upper_Triangular_R

    !!----
    !!---- Logical Function Zbelong(V)
    !!----    real(kind=cp),   dimension(:,:), intent( in) :: v
    !!----                      or
    !!----    real(kind=cp),   dimension(:),   intent( in) :: v
    !!----                      or
    !!----    real(kind=cp),                   intent( in) :: v
    !!----
    !!----    Provides the value .TRUE. if the real number (or array) V is close enough
    !!----    (whithin EPS) to an integer.
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Logical Function ZbelongM(V)
    !!--++    real(kind=cp),   dimension(:,:), intent( in) :: v
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if a real array is an Integer matrix
    !!--++
    !!--++ Update: February - 2005
    !!
    Function ZbelongM(v) Result(belong)
       !---- Argument ----!
       real(kind=cp),   dimension(:,:), intent( in) :: v
       logical                                      :: belong

       !---- Local variables ----!
       real(kind=cp),   dimension(size(v,1),size(v,2)) :: vec

       vec= abs(real(nint (v))-v)
       belong=.not. ANY(vec > epss)

       return
    End Function ZbelongM

    !!--++
    !!--++ Logical Function ZbelongN(A)
    !!--++    real(kind=cp),  intent(in)  :: a
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if a real number is an Integer
    !!--++
    !!--++ Update: February - 2005
   !!
    Function ZbelongN(a) Result(belong)
       !---- Argument ----!
       real(kind=cp), intent( in) :: a
       logical                    :: belong

       belong=.false.
       if (abs(real(nint (a))-a) > epss) return
       belong=.true.

       return
    End Function ZbelongN

    !!--++
    !!--++ Logical Function ZbelongV(V)
    !!--++    real(kind=sp),   dimension(:), intent( in) :: v
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if a real vector is an Integer vector
    !!--++
    !!--++ Update: February - 2005
    !!
    Function ZbelongV(v) Result(belong)
       !---- Argument ----!
       real(kind=cp),   dimension(:), intent( in) :: v
       logical                                    :: belong

       !---- Local variables ----!
       integer                             :: i
       real(kind=cp),   dimension(size(v)) :: vec

       belong=.false.
       vec= abs(real(nint (v))-v)
       do i=1,size(v)
          if (vec(i) > epss) return
       end do
       belong=.true.

       return
    End Function ZbelongV

    !---------------------!
    !---- Subroutines ----!
    !---------------------!

    !!----
    !!---- Subroutine Init_Err_Mathgen()
    !!----
    !!----    Initialize the errors flags in CFML_Math_General
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Init_Err_MathGen()

       ERR_MathGen=.false.
       ERR_MathGen_Mess=" "

       return
    End Subroutine Init_Err_MathGen

    !!----
    !!---- Subroutine Set_Epsg(Neweps)
    !!----    real(kind=cp), intent( in) :: neweps
    !!----
    !!----    Sets global EPSS to the value "neweps"
    !!----
    !!---- Update: April - 2005
    !!
    Subroutine Set_Epsg(Neweps)
       !---- Arguments ----!
       real(kind=cp), intent( in) :: neweps

       epss=neweps

       return
    End Subroutine Set_Epsg

    !!----
    !!---- Subroutine Set_Epsg_Default()
    !!----
    !!----    Sets global EPSS to the default value: epss=1.0E-5_sp
    !!----
    !!---- Update: April - 2005
    !!
    Subroutine Set_Epsg_Default()

       epss=1.0E-5_sp

       return
    End Subroutine Set_Epsg_Default

    !!----
    !!---- Subroutine Rtan(y,x,ang,deg)
    !!----    real(sp/dp),               intent( in) :: x,y
    !!----    real(sp/dp),               intent(out) :: ang
    !!----    character(len=*),optional, intent( in) :: deg
    !!----
    !!----    Returns ang=arctan(y/x) in the quadrant where the signs sin(ang) and
    !!----    cos(ang) are those of y and x. If deg is present, return ang in degrees.
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine Rtan_dp(y,x,ang,deg)
    !!--++    real(dp),                  intent( in) :: x,y
    !!--++    real(dp),                  intent(out) :: ang
    !!--++    character(len=*),optional, intent( in) :: deg
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Returns ang=arctan(y/x) in the quadrant where the signs sin(ang) and
    !!--++    cos(ang) are those of y and x. If deg is present, return ang in degrees.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Rtan_dp(y,x,ang,deg)
       !---- Arguments ----!
       real(kind=dp),              Intent( In)   :: x,y
       real(kind=dp),              Intent(Out)   :: ang
       character(len=*), optional, Intent( In)   :: deg

       !---- Local variables ----!
       real(kind=dp):: abx,aby

       abx=abs(x)
       aby=abs(y)
       if ((abx < eps) .and. (aby < eps)) then
          ang = 0.0_dp
          return
       else if(abx < eps) then
          ang = pi/2.0_dp
       else if(aby < abx) then
          ang = atan(aby/abx)
          if(x < 0.0_dp) ang = pi-ang
       else
          ang = pi/2.0_dp - atan(abx/aby)
          if(x < 0.0_dp) ang = pi-ang
       end if
       if (y < 0.0_dp) ang = -ang
       if (present(deg)) ang = ang*to_deg

       return
    End Subroutine Rtan_dp

    !!--++
    !!--++ Subroutine Rtan_sp(x,y,ang,deg)
    !!--++    real(sp),                  intent( in) :: x,y
    !!--++    real(sp),                  intent(out) :: ang
    !!--++    character(len=*),optional, intent( in) :: deg
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Returns ang=arctan(y/x) in the quadrant where the signs sin(ang) and
    !!--++    cos(ang) are those of y and x. If deg is present, return ang in degrees.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Rtan_sp(y,x,ang,deg)
       !---- Arguments ----!
       real(kind=sp),              Intent( In)   :: x,y
       real(kind=sp),              Intent(Out)   :: ang
       character(len=*), optional, Intent( In)   :: deg

       !---- local variables ----!
       real(kind=sp):: abx,aby

       abx=abs(x)
       aby=abs(y)
       if ((abx < eps) .and. (aby < eps)) then
          ang = 0.0_sp
          return
       else if(abx < eps) then
          ang = pi/2.0_sp
       else if(aby < abx) then
          ang = atan(aby/abx)
          if(x < 0.0_sp) ang = pi-ang
       else
          ang = pi/2.0_sp - atan(abx/aby)
          if(x < 0.0_sp) ang = pi-ang
       end if
       if(y < 0.0_sp) ang = -ang
       if (present(deg)) ang = ang*to_deg

       return
    End Subroutine Rtan_sp

    !!----
    !!----  Subroutine Co_Prime_Vector(V,Cop,F)
    !!----     integer, dimension(:), intent(in)  :: v      !input integer vector
    !!----     integer, dimension(:), intent(out) :: cop    !Output co-prime vector
    !!----     integer,  optional,    intent(out) :: f      !Common multiplicative factor
    !!----
    !!----     Calculates the co-prime vector (cop) parallel to the input vector (v)
    !!----     It uses the list of the first thousand prime numbers.
    !!----
    !!----   Updated: January 2012 (JRC), copied from Nodal_Indices (Laue_Mod) in July 2013 (JRC)
    !!----
    Subroutine Co_Prime_Vector(V,Cop,f)
       !---- Arguments ----!
       integer, dimension(:), intent(in)  :: v
       integer, dimension(:), intent(out) :: cop
       integer,  optional,    intent(out) :: f

       !---- Local variables ----!
       integer                     :: i,j,max_ind,k,im,dimv,n

       cop=v
       n=1
       if (present(f)) f=1
       max_ind=maxval(abs(cop))
       !---- If the maximum value of the indices is 1 they are already coprimes
       if (max_ind <= 1) return
       !---- Indices greater than 1
       dimv=size(v)
       im=0
       do i=1,size(primes)
          if(primes(i) > max_ind) then  !primes is an array within this module
             im=i
             exit
          end if
       end do
       if(im == 0) return
       do_p: do i=1,im
         k=primes(i)
         do
           do j=1,dimv
              if( mod(cop(j),k) /= 0) cycle do_p
           end do
           n=n*k
           cop=cop/k
         end do
       end do do_p

       if (present(f)) f=n

       return
    End Subroutine Co_Prime_vector

    !!----
    !!---- Subroutine Determinant(A,n,determ)
    !!----    complex/real(sp), dimension(:,:), intent( in) :: A      !input square matrix (n,n)
    !!----    integer,                          intent( in) :: n      !actual dimension of A
    !!----    real(kind=sp),                    intent(out) :: determ !det(A) if real
    !!----                                                             det(AR)^2 + det(AI)^2 if complex
    !!----
    !!----    Calculates the determinant of a real square matrix.
    !!----    Calculates the pseudo-determinant of a complex square matrix.
    !!----    The calculated value is only useful for linear dependency purposes.
    !!----    It tell us if the complex matrix is singular or not.
    !!--..
    !!--..    Calculates the determinant of a complex square matrix selected from a rectangular
    !!--..    matrix A, n x m, where m >= n. determ=determinant_of_A(1:n,icol:icol+n-1)
    !!--..    If icol is absent, the calculation is performed as if icol=1.
    !!--..    If icol+n-1 > m, or m < n, determ is set to 0.0 and an error message is generated.
    !!----
    !!--..    P R O V I S I O N A L (The determinant of A is not calculated at present)
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine Determinant_C(A,n,determ)
    !!--++    complex,          dimension(:,:), intent( in) :: A      !input square matrix (n,n)
    !!--++    integer,                          intent( in) :: n      !actual dimension of A
    !!--++    real(kind=cp),                    intent(out) :: determ !det(A) if real
    !!--++                                                             det(AR)^2 + det(AI)^2 if complex
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the determinant of a real square matrix.
    !!--++    Calculates the pseudo-determinant of a complex square matrix.
    !!--++    The calculated value is only useful for linear dependency purposes.
    !!--++    It tell us if the complex matrix is singular or not.
    !!--++
    !!--++    P R O V I S I O N A L (The determinant of A is not calculated at present)
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Determinant_C(A,n,determ)
       !---- Arguments ----!
       complex, dimension(:,:), intent( in) :: A
       integer,                 intent( in) :: n
       real(kind=cp),           intent(out) :: determ

       !---- local variables ----!
       real(kind=cp),    dimension(2*n,2*n) :: AC   !real square matrix
       real(kind=cp)                        :: d
       integer                              :: i,nn
       logical                              :: singular

       nn=2*n
       AC(  1:n ,  1:n ) =  real(A(1:n ,1:n))
       AC(n+1:nn,  1:n ) = aimag(A(1:n ,1:n))
       AC(n+1:nn,n+1:nn) =    AC(  1:n ,1:n)
       AC(  1:n ,n+1:nn) =   -AC(n+1:nn,1:n)

       call lu_decomp(ac(1:nn,1:nn),d,singular)

       if (singular) then
          determ=0.0
       else
          determ=0.0
          do i=1,nn
             d=d*sign(1.0_cp,ac(i,i))
             determ=determ+ log(abs(ac(i,i)))
          end do
          determ=d*exp(determ)
       end if

       return
    End Subroutine Determinant_C

    !!--++
    !!--++ Subroutine Determinant_R(A,n,determ)
    !!--++    real(kind=cp), dimension(:,:),intent( in) :: A   (input square matrix (n,n))
    !!--++    integer,                      intent( in) :: n   (actual dimension of A)
    !!--++    real(kind=cp),                intent(out) :: determ  (determinant )
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the determinant of a real square matrix.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Determinant_R(A,n,determ)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:), intent( in) :: A
       integer,                       intent( in) :: n
       real(kind=cp),                 intent(out) :: determ

       !---- local variables ----!
       real(kind=cp),    dimension(n,n)  :: AC
       real(kind=cp)                     :: d
       integer                           :: i
       logical                           :: singular

       ac=A(1:n,1:n)
       call lu_decomp(ac,d,singular)

       if (singular) then
          determ=0.0
       else
          determ=0.0
          do i=1,n
             d=d*sign(1.0_cp,ac(i,i))
             determ=determ + log(abs(ac(i,i)))
          end do
          determ=d*exp(determ)
       end if

       return
    End Subroutine Determinant_R

    !!----
    !!---- Subroutine Diagonalize_SH(A,N,E_val,E_vect)
    !!----    complex/real,      dimension(:,:), intent( in)  :: A
    !!----    integer,                           intent( in)  :: n
    !!----    real(kind=cp),     dimension(:),   intent(out)  :: E_val
    !!----    complex, optional, dimension(:,:), intent(out)  :: E_vect
    !!----
    !!----    Diagonalize Symmetric/Hermitian matrices.
    !!----    The eigen_values E_val are sorted in descending order. The columns
    !!----    of E_vect are the corresponding eigenvectors.
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine Diagonalize_Herm(a,n,e_val,e_vect)
    !!--++    complex,           dimension(:,:), intent( in)  :: A
    !!--++    integer,                           intent( in)  :: n
    !!--++    real(kind=cp),     dimension(:),   intent(out)  :: E_val
    !!--++    complex, optional, dimension(:,:), intent(out)  :: E_vect
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Diagonalize Hermitian matrices.
    !!--++    The eigen_values E_val are sorted in descending order. The columns
    !!--++    of E_vect are the corresponding eigenvectors.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Diagonalize_Herm(a,n,e_val,e_vect)
       !---- Arguments ----!
       complex,           dimension(:,:), intent( in)  :: A
       integer,                           intent( in)  :: n
       real(kind=cp),     dimension(:),   intent(out)  :: E_val
       complex, optional, dimension(:,:), intent(out)  :: E_vect

       !---- Local variables ----!
       real(kind=cp),        dimension(2*n,2*n)   :: aux
       real(kind=cp),        dimension(2*n)       :: e,d
       integer :: nn

       e_val=0.0
       call init_err_mathgen()
       if (n > size(A,1) .or. n > size(A,2)) then
          ERR_MathGen=.true.
          ERR_MathGen_Mess=" Diagonalize_HERM: Error in dimension of input matrix: A(m,m) with m < n "
          return
       end if

       nn=2*n
       aux(  1:n ,  1:n ) =  real(a(1:n ,1:n))   !      (  U   V )
       aux(n+1:nn,n+1:nn) =  real(a(1:n ,1:n))   !   M=(          ),   A = U + i V
       aux(n+1:nn,  1:n ) = aimag(a(1:n ,1:n))   !      ( -V   U )
       aux(  1:n ,n+1:nn) =-aimag(a(1:n ,1:n))   !

       if (present(E_vect)) then
          call tred2(aux,nn,d,e)
          call tqli2(d,e,nn,aux)
          call eigsrt(d,aux,nn,1)
          e_vect(1:n,1:n)=cmplx(aux(1:n,1:nn:2),aux(n+1:nn,1:nn:2))
       else
          call tred1(aux,nn,d,e)
          call tqli1(d,e,nn)
          call eigsrt(d,aux,nn,0)
       end if
       e_val(1:n)=d(1:nn:2)

       return
    End Subroutine Diagonalize_Herm

    !!--++
    !!--++ Subroutine Diagonalize_Symm(a,n,e_val,e_vect)
    !!--++    real(kind=cp)            dimension(:,:),intent( in)  :: A      (input matrix with)
    !!--++    integer,                                intent( in)  :: n      (actual dimension)
    !!--++    real(kind=cp),           dimension(:),  intent(out)  :: E_val  (eigenvalues)
    !!--++    real(kind=cp), optional, dimension(:,:),intent(out)  :: E_vect (eigenvectors)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Diagonalize symmetric matrices
    !!--++    The eigen_values E_val are sorted in descending order. The columns
    !!--++    of E_vect are the corresponding eigenvectors.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Diagonalize_Symm(A,n,E_Val,E_vect)
       !---- Arguments ----!
       real(kind=cp),           dimension(:,:), intent( in)  :: A
       integer,                                 intent( in)  :: n
       real(kind=cp),           dimension(:),   intent(out)  :: E_val
       real(kind=cp), optional, dimension(:,:), intent(out)  :: E_vect

       !---- Local variables ----!
       real(kind=cp),        dimension(n,n)   :: aux
       real(kind=cp),        dimension(n)     :: e

       e_val=0.0
       call init_err_mathgen()
       if (n > size(A,1) .or. n > size(A,2)) then
          ERR_MathGen=.true.
          ERR_MathGen_Mess=" Diagonalize_SYMM: Error in dimension of input matrix: A(m,m) with m < n "
          return
       end if

       aux=a(1:n,1:n)
       if (present(E_vect)) then
          call tred2(aux,n,E_val,e)
          call tqli2(E_val,e,n,aux)
          call eigsrt(E_val,aux,n,1)
          e_vect(1:n,1:n)=aux
       else
          call tred1(aux,n,E_val,e)
          call tqli1(E_val,e,n)
          call eigsrt(E_val,aux,n,0)
       end if

       return
    End Subroutine Diagonalize_Symm

    !!--++
    !!--++ Subroutine Eigsrt(d,v,n,io)
    !!--++    real(kind=cp), dimension(:),   intent(in out) :: d
    !!--++    real(kind=cp), dimension(:,:), intent(in out) :: v
    !!--++    integer,                       intent (in)    :: n
    !!--++    integer,                       intent (in)    :: io
    !!--++
    !!--++    (PRIVATE)
    !!--++    Subroutine for sorting eigenvalues in d(n) and eigenvectors
    !!--++    in columns of v(n,n). Sorts d(n) in descending order and
    !!--++    rearranges v(n,n) correspondingly. The method is the straight
    !!--++    insertion. If io=0 order  only the eigenvalues are treated.
    !!--++    Adapted from Numerical Recipes. Valid for hermitian matrices
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Eigsrt(d,v,n,io)
       !---- Arguments ----!
       real(kind=cp), dimension(:),   intent(in out) :: d
       real(kind=cp), dimension(:,:), intent(in out) :: v
       integer,                       intent(in)     :: n
       integer,                       intent(in)     :: io

       !---- Local Variables ----!
       integer          :: i,j,k
       real(kind=cp)    :: p

       do i=1,n-1
          k=i
          p=d(i)
          do j=i+1,n
             if (d(j) >= p) then
                k=j
                p=d(j)
             end if
          end do
          if (k /= i) then
             d(k)=d(i)
             d(i)=p
             if (io == 1) then
                do j=1,n
                   p=v(j,i)
                   v(j,i)=v(j,k)
                   v(j,k)=p
                end do
             end if
          end if
       end do

       return
    End Subroutine Eigsrt

    !!----
    !!---- Subroutine First_Derivative(x, y, n, d2y, d1y)
    !!----    real(kind=cp),    intent(in),     dimension(:) :: x     !  In -> Array X
    !!----    real(kind=cp),    intent(in),     dimension(:) :: y     !  In -> Array Yi=F(Xi)
    !!----    integer ,         intent(in)                   :: n     !  In -> Dimension of X, Y
    !!----    real(kind=cp),    intent(in),     dimension(:) :: d2y   !  In -> array containing second derivatives
    !!----                                                                     at the given points
    !!----    real(kind=cp),    intent(out),    dimension(:) :: d1y   ! Out -> array containing first derivatives
    !!----                                                                     at the given points
    !!----
    !!----    Calculate the First derivate values of the N points
    !!----
    !!---- Update: January - 2006
    !!
    Subroutine First_Derivative(x,y,n,d2y,d1y)
       !---- Arguments ----!
       real(kind=cp), dimension(:), intent(in)  :: x
       real(kind=cp), dimension(:), intent(in)  :: y
       integer ,                    intent(in)  :: n
       real(kind=cp), dimension(:), intent(in)  :: d2y
       real(kind=cp), dimension(:), intent(out) :: d1y

       !---- Local Variables ----!
       integer       :: i
       real(kind=cp) :: step, x0, y0, y1, y2

       do i=1,n
         if (i /= n) then
           step = x(i+1)-x(i)
         end if
         x0 = x(i) - step/2.0
         call splint(x,y, d2y, n, x0, y0)
         y1 = y0
         x0 = x(i) + step/2
         call splint(x,y, d2y, n, x0, y0)
         y2 = y0
         d1y(i) = (y2 - y1) / step
       end do

       return
    End Subroutine First_Derivative

    !!----
    !!---- Subroutine In_Sort(id,n,p,q)
    !!----    integer, dimension(:), intent(in) :: id  !Integer array to be sorted
    !!----    integer,               intent(in) :: n   !Number items in the array
    !!----    integer, dimension(:), intent(in) :: p   !Initial pointer from a previous related call
    !!----    integer, dimension(:), intent(out):: q   !Final pointer doing the sort of id
    !!--<<
    !!----    Subroutine to order in ascending mode the integer array "id".
    !!----    The input value "n" is the number of items to be ordered in "id".
    !!----    The array "p" is the initial pointer to "id" (coming from a previous call)
    !!----    The final pointer holding the order of items.
    !!-->>
    !!----
    !!---- Update: November - 2008
    !!
    Subroutine In_Sort(id,n,p,q)
       !---- Arguments ----!
       integer, dimension(:), intent(in) :: id  !Integer array to be sorted
       integer,               intent(in) :: n   !Number items in the array
       integer, dimension(:), intent(in) :: p   !Initial pointer from a previous related call
       integer, dimension(:), intent(out):: q   !Final pointer doing the sort of id

       !--- Local Variables ----!
       integer :: i,j,k,l,m
       integer, dimension(:),allocatable :: it

       l=minval(id)
       m=maxval(id)
       l=l-1
       m=m-l
       allocate(it(m))
       it(1:m)=0
       do i=1,n
          j=id(p(i))-l
          it(j)=it(j)+1
       end do
       j=0
       do i=1,m
          k=j
          j=j+it(i)
          it(i)=k
       end do
       do i=1,n
          j=id(p(i))-l
          it(j)=it(j)+1
          j=it(j)
          q(j)=p(i)
       end do

       return
    End Subroutine In_Sort

    !!----
    !!---- Subroutine Invert_Matrix(a,b,singular,perm)
    !!----    real(kind=cp), dimension(:,:),  intent( in) :: a
    !!----    real(kind=cp), dimension(:,:),  intent(out) :: b
    !!----    LOGICAL,                        intent(out) :: singular
    !!----    integer, dimension(:),optional, intent(out) :: perm
    !!--<<
    !!----    Subroutine to invert a real matrix using LU decomposition.
    !!----    In case of singular matrix (singular=.true.) instead of the inverse
    !!----    matrix, the subroutine provides the LU decomposed matrix as used
    !!----    in Numerical Recipes.
    !!----    The input matrix is preserved and its inverse (or its LU decomposition)
    !!----    is provided in "b". The optional argument "perm" holds the row permutation
    !!----    performed to obtain the LU decomposition.
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Invert_Matrix(a,b,singular,perm)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:),  intent(in ) :: a
       real(kind=cp), dimension(:,:),  intent(out) :: b
       logical,                        intent(out) :: singular
       integer, dimension(:),optional, intent(out) :: perm

       !---- Local variables ----!
       integer                                       :: i,n
       integer,       dimension(size(a,1))           :: indx
       real(kind=cp)                                 :: d, det
       real(kind=cp), dimension(size(a,1),size(a,1)) :: lu

       n=size(a,1)
       lu=a(1:n,1:n)

       call LU_Decomp(lu,d,singular,indx)
       if (present(perm)) perm(1:n)=indx(1:n)

       if (singular) then
          b=lu
          return
       else
          det=0.0
          do i=1,n
             d=d*sign(1.0_cp,lu(i,i))
             det=det + log(abs(lu(i,i)))
          end do
          det=d*exp(det)
          if (abs(det) <= 1.0e-36) then
             singular=.true.
             b=lu
             return
          end if
       end if

       b=0.0
       do i=1,n
          b(i,i)=1.0
          call LU_backsub(lu,indx,b(:,i))
       end do

       return
    End Subroutine Invert_Matrix

    !!----
    !!---- Subroutine Linear_Dependent(a,na,b,nb,mb,info)
    !!----    complex/integer/real(kind=cp), dimension(:),   intent(in)  :: a
    !!----    complex/integer/real(kind=cp), dimension(:,:), intent(in)  :: b
    !!----    integer,                                       intent(in)  :: na,nb,mb
    !!----    logical,                                       intent(out) :: info
    !!--<<
    !!----    Provides the value .TRUE. if the vector A is linear dependent of the
    !!----    vectors constituting the rows (columns) of the matrix B. In input nb & mb
    !!----    are the number of rows and columns of B to be considered. The actual
    !!----    dimension of vector a should be na=max(nb,mb).
    !!----    The problem is equivalent to determine the rank (in algebraic sense)
    !!----    of the composite matrix C(nb+1,mb)=(B/A) or C(nb,mb+1)=(B|A). In the first
    !!----    case it is supposed that na = mb and in the second na = nb.
    !!----    and the rank of B is min(nb, mb). If na /= nb and na /= mb an error condition
    !!----    is generated. The function uses floating arithmetic for all types.
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine Linear_DependentC(a,na,b,nb,mb,info)
    !!--++    complex, dimension(:),   intent(in)  :: a
    !!--++    complex, dimension(:,:), intent(in)  :: b
    !!--++    integer,                 intent(in)  :: na,nb,mb
    !!--++    logical,                 intent(out) :: info
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Provides the value .TRUE. if the vector A is linear dependent of the
    !!--++    vectors constituting the rows (columns) of the matrix B. In input nb & mb
    !!--++    are the number of rows and columns of B to be considered. The actual
    !!--++    dimension of vector a should be na=max(nb,mb).
    !!--++    The problem is equivalent to determine the rank (in algebraic sense)
    !!--++    of the composite matrix C(nb+1,mb)=(B/A) or C(nb,mb+1)=(B|A). In the first
    !!--++    case it is supposed that na = mb and in the second na = nb.
    !!--++    and the rank of B is min(nb, mb). If na /= nb and na /= mb an error condition
    !!--++    is generated
    !!--++
    !!--++    For the case of complex vectors in Cn the problem can be reduced to real vectors
    !!--++    of dimension R2n. Each complex vector contributes as two real vectors of dimension
    !!--++    2n: (R,I) and (-I,R). A complex vector V is linearly dependent on n complex vectors
    !!--++    if V can be written as: V = Sigma{j=1,n}(Cj.Vj), with Cj complex numbers and Vj
    !!--++    having n complex components. One may write:
    !!--++
    !!--++     V = Sigma{j=1,n}(Cj.Vj)
    !!--++     (R,I) = Sigma{j=1,n} (Cjr Vj + i Cji Vj) = Sigma{j=1,n} (Cjr (Rj,Ij) +  Cji (-Ij,Rj) )
    !!--++     (R,I) = Sigma{j=1,n} (aj (Rj,Ij) + bj (-Ij,Rj) )  = Sigma{j=1,2n} (Aj.Uj)
    !!--++     Were Uj=(Rj,Ij) and U(j+1)= (-Ij,Rj)
    !!--++
    !!--++    The function uses floating arithmetic for all types.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Linear_DependentC(A,na,B,nb,mb,info)
       !---- Arguments ----!
       complex, dimension(:),   intent(in)  :: a
       complex, dimension(:,:), intent(in)  :: b
       integer,                 intent(in)  :: na,nb,mb
       logical,                 intent(out) :: info

       !---- Local variables ----!
       integer                                                     :: r,n1
       real(kind=dp), parameter                                    :: tol= 100.0_dp*deps
       real(kind=dp), dimension(2*max(nb+1,mb+1),2*max(nb+1,mb+1)) :: c

       c=0.0
       call init_err_mathgen()
       info=.true.
       if (nb > size(b,1) .or. mb > size(b,2) .or. na > size(a) ) then
          ERR_MathGen=.true.
          ERR_MathGen_Mess=" Linear_DependentC: Error in dimension of input matrix or vector"
          return
       end if

       if ( na == mb) then
          n1=2*nb+1
          if(n1+1 > 2*mb) return !the vector is linear dependent
          c(1:nb,           1:mb) =  real(b(1:nb,1:mb))
          c(1:nb,     mb+1:mb+na) = aimag(b(1:nb,1:mb))
          c(nb+1:2*nb,      1:mb) =-aimag(b(1:nb,1:mb))
          c(nb+1:2*nb,mb+1:mb+na) =  real(b(1:nb,1:mb))
          c(n1,             1:mb) =  real(a(1:na))
          c(n1,      mb+1:mb+na ) = aimag(a(1:na))
          c(n1+1,           1:mb) =-aimag(a(1:na))
          c(n1+1,    mb+1:mb+na ) =  real(a(1:na))
          call rank(c,tol,r)
          if(r == min(n1+1,2*mb)) info=.false.
       else if( na == nb) then
          n1=2*mb+1
          if(n1+1 > 2*nb) return !the vector is linear dependent
          c(1:nb,           1:mb) =  real(b(1:nb,1:mb))
          c(nb+1:nb+na,     1:mb) = aimag(b(1:nb,1:mb))
          c(1:nb,      mb+1:2*mb) =-aimag(b(1:nb,1:mb))
          c(nb+1:nb+na,mb+1:2*mb) =  real(b(1:nb,1:mb))
          c(1:na,             n1) =  real(a(1:na))
          c(nb+1:nb+na,       n1) = aimag(a(1:na))
          c(1:na,           1+n1) =-aimag(a(1:na))
          c(nb+1:nb+na,     1+n1) =  real(a(1:na))
          call rank(c,tol,r)
          if(r == min(n1+1,2*nb)) info=.false.
       else
          ERR_MathGen=.true.
          ERR_MathGen_Mess=" Linear_DependentC: input dimension of vector incompatible with matrix"
       end if

       return
    End Subroutine Linear_DependentC

    !!--++
    !!--++ Subroutine Linear_DependentI(a,na,b,nb,mb,info)
    !!--++    integer, dimension(:),   intent(in)  :: a
    !!--++    integer, dimension(:,:), intent(in)  :: b
    !!--++    integer,                 intent(in)  :: na,nb,mb
    !!--++    logical,                 intent(out) :: info
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Provides the value .TRUE. if the vector A is linear dependent of the
    !!--++    vectors constituting the rows (columns) of the matrix B. In input nb & mb
    !!--++    are the number of rows and columns of B to be considered. The actual
    !!--++    dimension of vector a should be na=max(nb,mb).
    !!--++    The problem is equivalent to determine the rank (in algebraic sense)
    !!--++    of the composite matrix C(nb+1,mb)=(B/A) or C(nb,mb+1)=(B|A). In the first
    !!--++    case it is supposed that na = mb and in the second na = nb.
    !!--++    and the rank of B is min(nb, mb). If na /= nb and na /= mb an error condition
    !!--++    is generated
    !!--++    The function uses floating arithmetic for all types.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Linear_DependentI(A,na,B,nb,mb,info)
       !---- Arguments ----!
       integer, dimension(:),   intent(in)  :: a
       integer, dimension(:,:), intent(in)  :: b
       integer,                 intent(in)  :: na,nb,mb
       logical,                 intent(out) :: info

       !---- Local variables ----!
       integer                                                 :: r,n1
       real(kind=dp), parameter                                :: tol= 100.0_dp*deps
       real(kind=dp), dimension(max(nb+1,mb+1),max(nb+1,mb+1)) :: c

       c=0.0
       call init_err_mathgen()
       info=.true.
       if (nb > size(b,1) .or. mb > size(b,2) .or. na > size(a) ) then
          ERR_MathGen=.true.
          ERR_MathGen_Mess=" Linear_DependentI: Error in dimension of input matrix or vector"
          return
       end if

       if ( na == mb) then
          n1=nb+1
          if(n1 > mb) return !the vector is linear dependent
          c(1:nb,1:mb)=real(b(1:nb,1:mb))
          c(n1,  1:mb)=real(a(1:na))      !C(nb+1,mb)
          call rank(c,tol,r)
          if(r == min(n1,mb)) info=.false.
       else if( na == nb) then
          n1=mb+1
          if(n1 > nb) return !the vector is linear dependent
          c(1:nb,1:mb)=real(b(1:nb,1:mb))
          c(1:nb,  n1)=real(a(1:na))     !C(nb,mb+1)
          call rank(c,tol,r)
          if(r == min(n1,nb)) info=.false.
       else
          ERR_MathGen=.true.
          ERR_MathGen_Mess=" Linear_DependentI: input dimension of vector incompatible with matrix"
       end if

       return
    End Subroutine Linear_DependentI

    !!--++
    !!--++ Subroutine Linear_DependentR(a,na,b,nb,mb,info)
    !!--++    real(kind=cp), dimension(:),   intent(in)  :: a
    !!--++    real(kind=cp), dimension(:,:), intent(in)  :: b
    !!--++    integer,                       intent(in)  :: na,nb,mb
    !!--++    logical,                       intent(out) :: info
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Provides the value .TRUE. if the vector A is linear dependent of the
    !!--++    vectors constituting the rows (columns) of the matrix B. In input nb & mb
    !!--++    are the number of rows and columns of B to be considered. The actual
    !!--++    dimension of vector a should be na=max(nb,mb).
    !!--++    The problem is equivalent to determine the rank (in algebraic sense)
    !!--++    of the composite matrix C(nb+1,mb)=(B/A) or C(nb,mb+1)=(B|A). In the first
    !!--++    case it is supposed that na = mb and in the second na = nb.
    !!--++    and the rank of B is min(nb, mb). If na /= nb and na /= mb an error condition
    !!--++    is generated
    !!--++    The function uses floating arithmetic for all types.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Linear_DependentR(A,na,B,nb,mb,info)
       !---- Arguments ----!
       real(kind=cp), dimension(:),   intent(in)  :: a
       real(kind=cp), dimension(:,:), intent(in)  :: b
       integer,                       intent(in)  :: na,nb,mb
       logical,                       intent(out) :: info

       !---- Local Variables ----!
       integer                                                 :: r,n1
       real(kind=dp), parameter                                :: tol= 100.0_dp*deps
       real(kind=dp), dimension(max(nb+1,mb+1),max(nb+1,mb+1)) :: c

       c=0.0
       call init_err_mathgen()
       info=.true.
       if (nb > size(b,1) .or. mb > size(b,2) .or. na > size(a) ) then
          ERR_MathGen=.true.
          ERR_MathGen_Mess=" Linear_DependentR: Error in dimension of input matrix or vector"
          return
       end if

       if ( na == mb) then    !Vector added as an additional row
          n1=nb+1
          if(n1 > mb) return !the vector is linear dependent
          c(1:nb,1:mb)=b(1:nb,1:mb)
          c(n1,  1:mb)=a(1:na)      !C(nb+1,mb)
          call rank(c,tol,r)
          if(r == min(n1,mb)) info=.false.
       else if( na == nb) then   !Vector added as an additional column
          n1=mb+1
          if(n1 > nb) return !the vector is linear dependent
          c(1:nb,1:mb)=b(1:nb,1:mb)
          c(1:nb,  n1)=a(1:na)     !C(nb,mb+1)
          call rank(c,tol,r)
          if(r == min(n1,nb)) info=.false.
       else
          ERR_MathGen=.true.
          ERR_MathGen_Mess=" Linear_DependentR: input dimension of vector incompatible with matrix"
       end if

       return
    End Subroutine Linear_DependentR

    !!----
    !!---- Subroutine LU_Backsub(a,indx,b)
    !!----    real(kind=cp),    dimension(:,:),intent(in)     :: a
    !!----    integer,          dimension(:),  intent(in)     :: indx
    !!----    real(kind=cp),    dimension(:),  intent(in out) :: b
    !!--<<
    !!----    Adapted from Numerical Recipes.
    !!----    Solves the set of N linear equations A  X = B. Here the N x N matrix A is input,
    !!----    not as the original matrix A, but rather as its LU decomposition, determined
    !!----    by the routine LU_DECOMP. INDX is input as the permutation vector of length N
    !!----    returned by LU_DECOMP. B is input as the right-hand-side vector B,
    !!----    also of length N, and returns with the solution vector X.
    !!----    A and INDX are not modified by this routine and can be left in place for successive calls
    !!----    with different right-hand sides B. This routine takes into account the possibility that B will
    !!----    begin with many zero elements, so it is efficient for use in matrix inversion.
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine LU_Backsub(a,indx,b)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:), intent(in)     :: a
       integer,         dimension(:), intent(in)     :: indx
       real(kind=cp),   dimension(:), intent(in out) :: b

       !---- Local Variables ----!
       integer       :: i,ii,ll,n
       real(kind=cp) :: summ

       n=size(a,1)
       ii=0              !When ii is set to a positive value, it will become the index
       do i=1,n          !of the first nonvanishing element of b. We now do
          ll=indx(i)     !the forward substitution. The only new wrinkle is to
          summ=b(ll)     !unscramble the permutation as we go.
          b(ll)=b(i)
          if (ii /= 0) then
             summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
          else if(summ /= 0.0) then   !A nonzero element was encountered, so from now on
             ii=i                       !we will have to do the dot product above.
          end if
          b(i)=summ
       end do

       do i=n,1,-1       !Now we do the backsubstitution
          b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
       end do

       return
    End Subroutine LU_Backsub

    !!----
    !!---- Subroutine LU_Decomp(a,d,singular,indx)
    !!----    real(kind=cp),    dimension(:,:),intent(in out) :: a
    !!----    real(kind=cp),                   intent(out)    :: d
    !!----    logical,                         intent(out)    :: singular
    !!----    integer, dimension(:), optional, intent(out)    :: indx
    !!--<<
    !!----    Subroutine to make the LU decomposition of an input matrix A.
    !!----    The input matrix is destroyed and replaced by a matrix containing
    !!----    in its upper triangular part (plus diagonal) the matrix U. The
    !!----    lower triangular part contains the nontrivial part (Lii=1) of matrix L.
    !!----    The output is rowwise permutation of the initial matrix. The vector INDX
    !!----    recording the row permutation. D is output as +/-1 depending on whether
    !!----    the number of row interchanges was even or odd, respectively.
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine LU_Decomp(a,d,singular,indx)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:), intent(in out) :: a
       real(kind=cp),                 intent(out)    :: d
       logical,                       intent(out)    :: singular
       integer,  dimension(:), intent(out), optional :: indx

       !---- Local variables ----!
       real(kind=cp), dimension(size(a,1)):: vv  !vv stores the implicit scaling of each row.
       real(kind=cp), parameter           :: vtiny = 1.0e-20_sp !A small number.
       integer                            :: j,imax,n

       singular=.false.
       n=size(a,1)
       d=1.0                      !No row interchanges yet.
       vv=maxval(abs(a),dim=2)    !Loop over rows to get the implicit scaling information.
       if (any(abs(vv) <= vtiny)) then   !There is a row of zeros.
          singular=.true.
          return
       end if
       vv=1.0_sp/vv     !Save the scaling.
       do j=1,n
          imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j)))   !Find the pivot row.
          if (j /= imax) then                         !Do we need to interchange rows?
             call swap(a(imax,:),a(j,:))              !Yes, do so...
             d=-d                                     !...and change the parity of d.
             vv(imax)=vv(j)                           !Also interchange the scale factor.
          end if
          if (present(indx)) indx(j)=imax
          if (abs(a(j,j)) <= vtiny) then !If the pivot element is zero the matrix is singular.
             a(j,j)=vtiny                !(at least to the precision of the algorithm)
             singular=.true.             !For some applications on singular matrices,
             return                      !it is desirable to substitute vtiny for zero.
          end if                         !This is actually the present case
          a(j+1:n,j)=a(j+1:n,j)/a(j,j)                                    !Divide by the pivot element.
          a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))  !Reduce remaining submatrix.
       end do

       return
    End Subroutine LU_Decomp

    !!----
    !!---- Subroutine Matinv(a,n)
    !!----    real(kind=cp), dimension(:,:),intent(in out) :: a
    !!----    integer     ,                 intent(in)     :: n
    !!----
    !!----  Subroutine for inverting a real square matrix.
    !!----  The input matrix is replaced in output with its inverse.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Matinv(a,n)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:), intent(in out) :: a
       integer     ,                  intent(in)     :: n

       !---- Local variables ----!
       real(kind=cp)                 :: amax,savec
       integer, dimension(size(a,1)) :: ik,jk
       integer                       :: i,j,k,l

       !---- Subroutine to invert a real matrix ----!
       do k=1,n
          amax=0.0
          do
             do
                do i=k,n
                   do j=k,n
                      if (abs(amax)-abs(a(i,j)) > 0.0) cycle
                      amax=a(i,j)
                      ik(k)=i
                      jk(k)=j
                   end do
                end do
                i=ik(k)
                if (i-k < 0) cycle
                exit
             end do

             if (i-k /= 0) then
                do j=1,n
                   savec=a(k,j)
                   a(k,j)=a(i,j)
                   a(i,j)=-savec
                end do
             end if

             j=jk(k)
             if (j-k < 0) cycle
             exit
          end do

          if (j-k /= 0) then
             do i=1,n
                savec=a(i,k)
                a(i,k)=a(i,j)
                a(i,j)=-savec
             end do
          end if

          do i=1,n
             if (i-k /= 0)  then
                a(i,k)=-a(i,k)/amax
             end if
          end do
          do i=1,n
             do j=1,n
                if (i-k == 0 .or. j-k == 0) cycle
                a(i,j)=a(i,j)+a(i,k)*a(k,j)
             end do
          end do
          do j=1,n
             if (j-k == 0)   cycle
             a(k,j)=a(k,j)/amax
          end do
          a(k,k)=1.0/amax
       end do     !k

       do l=1,n
          k=n-l+1
          j=ik(k)
          if (j-k > 0) then
             do i=1,n
                savec=a(i,k)
                a(i,k)=-a(i,j)
                a(i,j)=savec
             end do
          end if
          i=jk(k)
          if (i-k > 0) then
             do j=1,n
                savec=a(k,j)
                a(k,j)=-a(i,j)
                a(i,j)=savec
             end do
          end if
       end do

       return
    End Subroutine Matinv

    !!--++
    !!--++ Subroutine Partition(A, marker)
    !!--++    character(len=*), dimension(:), intent(in out) :: A
    !!--++    integer,                        intent(out)    :: marker
    !!--++
    !!--++    (Private)
    !!--++    Utilised by Sort_Strings.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Partition(A, Marker)
       !---- Arguments ----!
       character(len=*), dimension(:), intent(in out) :: A
       integer,                        intent(   out) :: marker

       !---- Local variables ----!
       integer                  :: i, j
       character(len=len(A(1))) :: temp
       character(len=len(A(1))) :: x      ! pivot point

       x = A(1)
       i= 0
       j= size(A) + 1

       do
          j = j-1
          do
             if (A(j) <= x) exit
             j = j-1
          end do
          i = i+1
          do
             if (A(i) >= x) exit
             i = i+1
          end do
          if (i < j) then
             !---- exchange A(i) and A(j)
             temp = A(i)
             A(i) = A(j)
             A(j) = temp
          else if (i == j) then
             marker = i+1
             return
          else
             marker = i
             return
          end if
       end do

       return
    End Subroutine Partition

    !!----
    !!---- Subroutine Points_In_Line2D(X1, XN, N, XP)
    !!----    real(kind=cp), dimension(2),   intent(in)  :: X1   ! Point1 in 2D
    !!----    real(kind=cp), dimension(2),   intent(in)  :: XN   ! PointN in 2D
    !!----    integer,                       intent(in)  :: N    ! Number of Total points
    !!----    real(kind=cp), dimension(:,:), intent(out) :: XP   ! List of points
    !!----
    !!----    The routine calculate N points belonging to the line defined
    !!----    by X1 and Xn with equal distance between them. XP contains
    !!----    X1,X2,.....,XN points.
    !!----
    !!---- Update: April 2008
    !!
    Subroutine Points_In_Line2D(X1, XN, N, XP)
       !---- Arguments ----!
       real(kind=cp), dimension(2),   intent(in)  :: X1   ! Point1 in 2D
       real(kind=cp), dimension(2),   intent(in)  :: XN   ! PointN in 2D
       integer,                       intent(in)  :: N    ! Number of Total points
       real(kind=cp), dimension(:,:), intent(out) :: XP   ! List of points

       !---- Local Variables ----!
       integer :: i
       real(kind=cp)    :: ml,bl,dl,t
       real(kind=cp)    :: a,b,c,d
       real(kind=cp)    :: xa,xb

       xp=0.0

       if (n <= 1) return

       !---- Calculating the distance between two points to
       !---- eliminate rare considerations as the same point
       dl=sqrt( (xn(1)-x1(1))**2 + (xn(2)-x1(2))**2 )
       if (dl <= 0.0001) return

       !---- When N=2 is trivial case ----!
       if (n == 2) then
          xp(:,1)=x1
          xp(:,2)=xn
          return
       end if

       !---- Case 1: Y=cte ----!
       !Xn(2) and X1(2) are equal, then we have a line  with Y=cte
       if (abs(xn(2)-x1(2)) <= 0.0001) then
          dl=abs(xn(1)-x1(1))
          d=dl/real(n-1)
          xp(:,1)=x1
          if (xn(1) > x1(1)) then
             do i=2,n-1
                xp(1,i)=xp(1,i-1)+d
                xp(2,i)=xp(2,1)
             end do
          else
             do i=2,n-1
                xp(1,i)=xp(1,i-1)-d
                xp(2,i)=xp(2,1)
             end do
          end if
          xp(:,n)=xn

          return
       end if

       !---- Case 2: X=cte ----!
       !Xn(1) - X1(1) are equal, then we have a line with X=cte
       if (abs(xn(1)-x1(1)) <= 0.0001) then
          dl=abs(xn(2)-x1(2))
          d=dl/real(n-1)
          xp(:,1)=x1
          if (xn(2) > x1(2)) then
             do i=2,n-1
                xp(1,i)=xp(1,1)
                xp(2,i)=xp(2,i-1)+d
             end do
          else
             do i=2,n-1
                xp(1,i)=xp(1,1)
                xp(2,i)=xp(2,i-1)-d
             end do
          end if
          xp(:,n)=xn

          return
       end if

       !---- Case 3: General case ----!
       ml=(x1(2)-xn(2))/(x1(1)-xn(1))
       bl=x1(2) - (ml * x1(1))

       !---- Distance between X1 and XN ----!
       dl=sqrt( (xn(1)-x1(1))**2 + (xn(2)-x1(2))**2 )

       !---- Creating the list ----!
       a=ml**2 + 1.0
       b=2.0 *( ml*(bl-x1(2)) -x1(1) )

       xp(:,1)=x1
       do i=2,n-1
          t=(dl**2)*((real(i-1)/real(n-1))**2)
          c=(x1(2)-bl)**2 + x1(1)**2 - t

          xa=(-b + sqrt(b**2 - 4.0*a*c))/(2.0*a)
          xb=(-b - sqrt(b**2 - 4.0*a*c))/(2.0*a)
          if (x1(1) <= xa .and. xa <= xn(1)) then
             xp(1,i)=xa
             xp(2,i)=ml*xa+bl
          else
             xp(1,i)=xb
             xp(2,i)=ml*xb+bl
          end if
       end do
       xp(:,n)=xn

       return
    End Subroutine Points_In_Line2D

    !!----
    !!---- Subroutine Rank(a,tol,r)
    !!----    real(sp/dp), dimension(:,:), intent( in) :: a
    !!----    real(sp/dp),                 intent( in) :: tol
    !!----    integer,                     intent(out) :: r
    !!----
    !!----    Computes the rank (in algebraic sense) of the rectangular matrix A.
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine Rank_dp(a,tol,r)
    !!--++    real(dp), dimension(:,:), intent( in) :: a
    !!--++    real(dp),                 intent( in) :: tol
    !!--++    integer,                  intent(out) :: r
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Computes the rank (in algebraic sense) of the rectangular matrix A.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Rank_dp(a,tol,r)
       !---- Arguments ----!
       real(kind=dp), dimension(:,:),intent( in)      :: a
       real(kind=dp),                intent( in)      :: tol
       integer,                      intent(out)      :: r

       !---- Arguments ----!
       real(kind=dp), dimension(size(a,1),size(a,2))  :: u
       real(kind=dp), dimension(size(a,2))            :: w
       real(kind=dp), dimension(size(a,2),size(a,2))  :: v
       integer                                        :: i

       u=a
       call svdcmp(u,w,v)
       if (ERR_MathGen) then
          r=0
       else
          r=0
          do i=1,size(a,2)
             if(w(i) > tol) r=r+1
          end do
       end if

       return
    End Subroutine Rank_dp

    !!--++
    !!--++ Subroutine Rank_sp(a,tol,r)
    !!--++    real(sp), dimension(:,:), intent( in) :: a
    !!--++    real(sp),                 intent( in) :: tol
    !!--++    integer,                  intent(out) :: r
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Computes the rank (in algebraic sense) of the rectangular matrix A.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Rank_sp(a,tol,r)
       !---- Arguments ----!
       real(kind=sp), dimension(:,:),intent( in)      :: a
       real(kind=sp),                intent( in)      :: tol
       integer,                      intent(out)      :: r

       !---- Local variables ----!
       real(kind=sp), dimension(size(a,1),size(a,2))  :: u
       real(kind=sp), dimension(size(a,2))            :: w
       real(kind=sp), dimension(size(a,2),size(a,2))  :: v
       integer :: i

       u=a
       call svdcmp(u,w,v)
       if (ERR_MathGen) then
          r=0
       else
          r=0
          do i=1,size(a,2)
             if(w(i) > tol) r=r+1
          end do
       end if

       return
    End Subroutine Rank_sp

    !!----
    !!---- Subroutine Second_Derivative(x, y, n, d2y)
    !!----    real(kind=cp),    intent(in),     dimension(n) :: x     !  In -> Array X
    !!----    real(kind=cp),    intent(in),     dimension(n) :: y     !  In -> Array Yi=F(Xi)
    !!----    integer ,         intent(in)                   :: n     !  In -> Dimension of X, Y
    !!----    real(kind=cp),    intent(out),    dimension(n) :: d2y    ! Out -> array containing second derivatives
    !!----                                                                     at the given points
    !!----
    !!----    Calculate the second derivate of N Points
    !!----
    !!---- Update: January - 2006
    !!
    Subroutine Second_Derivative(x,y,n,d2y)
       !---- Arguments ----!
       real(kind=cp), dimension(:), intent(in)  :: x
       real(kind=cp), dimension(:), intent(in)  :: y
       integer ,                    intent(in)  :: n
       real(kind=cp), dimension(:), intent(out) :: d2y

       !---- Local Variables ----!
       integer                     :: i, k
       real(kind=cp), dimension(n) :: u
       real(kind=cp)               :: yp1, ypn, sig, p, qn, un

       yp1=(y(2) - y(1))   / (x(2) - x(1))     ! derivative at point 1
       ypn=(y(n) - y(n-1)) / (x(n) - x(n-1))   ! derivative at point n

       d2y(1)=-0.5
       u(1)=(3.0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)

       do i=2,n-1
          sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
          p=sig*d2y(i-1)+2.0
          d2y(i)=(sig-1.0)/p
          u(i)=(6.0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))  &
               /(x(i+1)-x(i-1))-sig*u(i-1))/p
       end do

       qn=0.5
       un=(3.0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
       d2y(n)=(un-qn*u(n-1))/(qn*d2y(n-1)+1.0)
       do k=n-1,1,-1
          d2y(k)=d2y(k)*d2y(k+1)+u(k)
       end do

       return
    End Subroutine Second_Derivative

    !!----
    !!---- Subroutine SmoothingVec(Y, N, NIter, Ys)
    !!----    real(kind=cp),    dimension(:),           intent(in out) :: Y      !  In Out-> Array to be smoothed
    !!----    integer,                                  intent(in)     :: N      !  In -> Number of points
    !!----    integer,                                  intent(in)     :: NIter  !  In -> Number of iterations
    !!----    real(kind=cp),    dimension(:), optional, intent(out)    :: datY   !  Out-> Array smoothed
    !!----
    !!----    Procedure to smooth the array values
    !!----
    !!---- Update: January - 2006
    !!
    Subroutine SmoothingVec(Y, N, Niter, Ys)
       !---- Arguments ----!
       real(kind=cp),dimension(:),            intent(in out) :: Y
       integer,                               intent(in)     :: n
       integer,                               intent(in)     :: niter
       real(kind=cp),dimension(:), optional,  intent(out)    :: Ys

       !---- Local Variables ----!
       integer                     :: n1, n2
       integer                     :: i, iter
       real(kind=cp), dimension (n):: datYs


       n1 = 4
       n2 = n-3

       do iter = 1 ,niter
          datYs(n1-1)=((Y(n1-2)+Y(n1))*10.0+(Y(n1-3)+Y(n1+1))*5.0+Y(n1+2))/31.0
          datYs(n1-2)=((Y(n1-3)+Y(n1-1))*10.0+Y(n1)*5.0+Y(n1+1))/26.0
          datYs(n1-3)=(Y(n1-2)*10.0+Y(n1-1)*5.0+Y(n1))/16.0

          do i=n1,n2
             datYs(i)=(Y(i-3)+Y(i+3)+5.0*(Y(i-2)+Y(i+2))+10.0*(Y(i-1)+Y(i+1)))/ 32.0
          end do

          datYs(n2+1)=((Y(n2+2)+Y(n2))*10.0+(Y(n2+3)+Y(n2-1))*5.0+Y(n2-2))/31.0
          datYs(n2+2)=((Y(n2+3)+Y(n2+1))*10.0+Y(n2)*5.0+Y(n2-1))/26.0
          datYs(n2+3)=(Y(n2+2)*10.0+Y(n2+1)*5.0+Y(n2))/16.0

          if(present(Ys)) then
             Ys(1:n) = datYs(1:n)
          else
             Y(1:n) = datYs(1:n)
          end if
       end do

       return
    End Subroutine SmoothingVec

    !!---
    !!---- Subroutine Sort(a,n,indx)
    !!----    integer/real(kind=cp)  dimension(:), intent( in) :: a
    !!----    integer,                             intent( in) :: n
    !!----    integer,               dimension(:), intent(out) :: indx
    !!----
    !!----    Sort an array such the a(indx(j)) is in ascending
    !!----    order for j=1,2,...,N.
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine Sort_I(Arr,N,Indx)
    !!--++    integer, dimension(:), intent( in) :: arr
    !!--++    integer,               intent( in) :: n
    !!--++    integer, dimension(:), intent(out) :: indx
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Sort an array such the arr(indx(j)) is in ascending
    !!--++    order for j=1,2,...,N.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Sort_I(arr,n,indx)
       !---- Arguments ----!
       integer, dimension(:), intent(in ) :: arr
       integer              , intent(in ) :: n
       integer, dimension(:), intent(out) :: indx

       !---- Local Variables ----!
       integer, parameter           :: m=7
       integer, parameter           :: nstack=50  !nstack=2log2(n)
       integer, dimension(nstack)   :: istack
       integer                      :: i,indxt,ir,itemp,j,jstack,k,l
       integer                      :: a

       call init_Err_MathGen()
       do j=1,n
          indx(j)=j
       end do

       istack=0
       jstack=0
       l=1
       ir=n
       do
          if (ir-l < m) then
             doext: do j=l+1,ir
                indxt=indx(j)
                a=arr(indxt)
                do i=j-1,1,-1
                   if (arr(indx(i)) <= a)  then
                      indx(i+1)=indxt
                      cycle doext
                   end if
                   indx(i+1)=indx(i)
                end do
                i=0
                indx(i+1)=indxt
             end do doext

             if (jstack == 0) exit
             ir=istack(jstack)
             l=istack(jstack-1)
             jstack=jstack-2
          else
             k=(l+ir)/2
             itemp=indx(k)
             indx(k)=indx(l+1)
             indx(l+1)=itemp
             if (arr(indx(l+1)) > arr(indx(ir)))then
                itemp=indx(l+1)
                indx(l+1)=indx(ir)
                indx(ir)=itemp
             end if
             if (arr(indx(l)) > arr(indx(ir)))then
                itemp=indx(l)
                indx(l)=indx(ir)
                indx(ir)=itemp
             end if
             if (arr(indx(l+1)) > arr(indx(l)))then
                itemp=indx(l+1)
                indx(l+1)=indx(l)
                indx(l)=itemp
             end if
             i=l+1
             j=ir
             indxt=indx(l)
             a=arr(indxt)
             do
                i=i+1
                if (arr(indx(i)) < a)  cycle
                do
                   j=j-1
                   if (arr(indx(j)) > a) cycle
                   exit
                end do
                if (j < i) exit
                itemp=indx(i)
                indx(i)=indx(j)
                indx(j)=itemp
             end do
             indx(l)=indx(j)
             indx(j)=indxt
             jstack=jstack+2
             if (jstack > nstack) then
                ERR_MathGen=.true.
                ERR_MathGen_Mess=" NSTACK too small in SORT"
                return
             end if
             if (ir-i+1 >= j-l) then
                istack(jstack)=ir
                istack(jstack-1)=i
                ir=j-1
             else
                istack(jstack)=j-1
                istack(jstack-1)=l
                l=i
             end if
          end if
       end do

       return
    End Subroutine Sort_I

    !!--++
    !!--++ Subroutine Sort_R(arr,n,indx)
    !!--++    real(kind=cp),dimension(:), intent( in) :: arr
    !!--++    integer,                    intent( in) :: n
    !!--++    integer,      dimension(:), intent(out) :: indx
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Sort an array such the arr(indx(j)) is in ascending
    !!--++    order for j=1,2,...,N.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Sort_R(arr,n,indx)
       !---- Arguments ----!
       real(kind=cp),dimension(:), intent(in) :: arr
       integer,                    intent(in) :: n
       integer,      dimension(:), intent(out):: indx

       !---- Local Variables ----!
       integer, parameter           :: m=7
       integer, parameter           :: nstack=50  !nstack=2log2(n)
       integer, dimension(nstack)   :: istack
       integer :: i,indxt,ir,itemp,j,jstack,k,l
       real(kind=cp)    :: a

       call init_Err_MathGen()
       do j=1,n
          indx(j)=j
       end do

       istack=0
       jstack=0
       l=1
       ir=n
       do
          if (ir-l < m) then
             doext: do j=l+1,ir
                indxt=indx(j)
                a=arr(indxt)
                do i=j-1,1,-1
                   if (arr(indx(i)) <= a)  then
                      indx(i+1)=indxt
                      cycle doext
                   end if
                   indx(i+1)=indx(i)
                end do
                i=0
                indx(i+1)=indxt
             end do doext

             if (jstack == 0) exit
             ir=istack(jstack)
             l=istack(jstack-1)
             jstack=jstack-2
          else
             k=(l+ir)/2
             itemp=indx(k)
             indx(k)=indx(l+1)
             indx(l+1)=itemp
             if (arr(indx(l+1)) > arr(indx(ir)))then
                itemp=indx(l+1)
                indx(l+1)=indx(ir)
                indx(ir)=itemp
             end if
             if (arr(indx(l)) > arr(indx(ir)))then
                itemp=indx(l)
                indx(l)=indx(ir)
                indx(ir)=itemp
             end if
             if (arr(indx(l+1)) > arr(indx(l)))then
                itemp=indx(l+1)
                indx(l+1)=indx(l)
                indx(l)=itemp
             end if
             i=l+1
             j=ir
             indxt=indx(l)
             a=arr(indxt)
             do
                i=i+1
                if (arr(indx(i)) < a)  cycle
                do
                   j=j-1
                   if (arr(indx(j)) > a) cycle
                   exit
                end do
                if (j < i) exit
                itemp=indx(i)
                indx(i)=indx(j)
                indx(j)=itemp
             end do
             indx(l)=indx(j)
             indx(j)=indxt
             jstack=jstack+2
             if (jstack > nstack) then
                ERR_MathGen=.true.
                ERR_MathGen_Mess=" NSTACK too small in SORT"
                return
             end if
             if (ir-i+1 >= j-l) then
                istack(jstack)=ir
                istack(jstack-1)=i
                ir=j-1
             else
                istack(jstack)=j-1
                istack(jstack-1)=l
                l=i
             end if
          end if
       end do

       return
    End Subroutine Sort_R

    !!---
    !!---- Subroutine Sort_Strings(arr)
    !!----    character(len=*), dimension(:), intent(in out) :: arr
    !!----
    !!----    Sort an array of string
    !!----
    !!---- Update: March - 2005
    !!
    Recursive Subroutine Sort_Strings(Arr)
       !---- Argument ----!
       character(len=*), dimension(:), intent(in out) :: Arr

       !---- Local variables ----!
       integer :: iq

       if (size(Arr) > 1) then
          call Partition(Arr, iq)
          call Sort_Strings(Arr(:iq-1))
          call Sort_Strings(Arr(iq:))
       end if

       return
    End Subroutine Sort_Strings

    !!----
    !!---- Subroutine Spline(x, y, n, yp1, ypn, y2)
    !!----    real(kind=cp),    intent(in),     dimension(n) :: x     !  In -> Array X
    !!----    real(kind=cp),    intent(in),     dimension(n) :: y     !  In -> Array Yi=F(Xi)
    !!----    integer ,         intent(in)                   :: n     !  In -> Dimension of X, Y
    !!----    real(kind=cp),    intent(in)                   :: yp1   !  In -> Derivate of Point 1
    !!----    real(kind=cp),    intent(in)                   :: ypn   !  In -> Derivate of Point N
    !!----    real(kind=cp),    intent(out),    dimension(n) :: y2    ! Out -> array containing second derivatives
    !!----                                                                     at the given points
    !!----    Spline  N points
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Spline(x,y,n,yp1,ypn,y2)
       !---- Arguments ----!
       real(kind=cp), dimension(:), intent(in)  :: x
       real(kind=cp), dimension(:), intent(in)  :: y
       integer ,                    intent(in)  :: n
       real(kind=cp),               intent(in)  :: yp1
       real(kind=cp),               intent(in)  :: ypn
       real(kind=cp), dimension(:), intent(out) :: y2

       !---- Local Variables ----!
       integer                     :: i, k
       real(kind=cp), dimension(n) :: u
       real(kind=cp)               :: sig, p, qn, un

       if (yp1 > 1.0e+30) then
          y2(1)=0.0
          u(1)=0.0
       else
          y2(1)=-0.5
          u(1)=(3.0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
       end if

       do i=2,n-1
          sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
          p=sig*y2(i-1)+2.0
          y2(i)=(sig-1.0)/p
          u(i)=(6.0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))  &
               /(x(i+1)-x(i-1))-sig*u(i-1))/p
       end do
       if (ypn > 1.0e+30) then
          qn=0.0
          un=0.0
       else
          qn=0.5
          un=(3.0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
       end if
       y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0)
       do k=n-1,1,-1
          y2(k)=y2(k)*y2(k+1)+u(k)
       end do

       return
    End Subroutine Spline

    !!----
    !!---- Subroutine Splint(x, y, y2, n, xp, yp)
    !!----    real(kind=cp),    intent(in), dimension(n) :: x  !  In -> Array X
    !!----    real(kind=cp),    intent(in), dimension(n) :: y  !  In -> Array Y=F(X)
    !!----    real(kind=cp),    intent(in), dimension(n) :: y2 !  In -> Array Second Derivatives in X
    !!----    integer ,         intent(in)               :: n  !  In -> Dimension of XA,YA,Y2A
    !!----    real(kind=cp),    intent(in)               :: xp !  In -> Point to evaluate
    !!----    real(kind=cp),    intent(out),             :: yp ! Out -> Interpoled value
    !!----
    !!----    Spline Interpolation
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Splint(xa,ya,y2a,n,x,y)
       !---- Arguments ----!
       real(kind=cp), dimension(:), intent(in)  :: xa
       real(kind=cp), dimension(:), intent(in)  :: ya
       real(kind=cp), dimension(:), intent(in)  :: y2a
       integer ,                    intent(in)  :: n
       real(kind=cp),               intent(in)  :: x
       real(kind=cp),               intent(out) :: y

       !---- Local Variables ----!
       integer          :: klo, khi, k
       real(kind=cp)    :: h, a, b

       klo=1
       khi=n
       do
          if (khi-klo > 1) then
             k=(khi+klo)/2
             if (xa(k) > x) then
                khi=k
             else
                klo=k
             end if
             cycle
          else
             exit
          end if
       end do

       h=xa(khi)-xa(klo)
       a=(xa(khi)-x)/h
       b=(x-xa(klo))/h
       y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)* y2a(khi))*(h**2)/6.0

       return
    End Subroutine Splint

    !!----
    !!---- Subroutine Svdcmp(a,w,v)
    !!----    real(sp/dp),dimension(:,:),intent(in out) :: a  !A(m,n)
    !!----    real(sp/dp),dimension(:),  intent(   out) :: w  !W(n)
    !!----    real(sp/dp),dimension(:,:),intent(   out) :: v  !V(n,n)
    !!--<<
    !!----    Given an M�N matrix A ,this routine computes its singular value decomposition,
    !!----    A = U �W �VT . The matrix U replaces A on output. The diagonal matrix of
    !!----    singular values W is output as the N-dimensional vector w. The N�N matrix V
    !!----    (not the transpose VT )is output as v .
    !!----    Adapted from Numerical Recipes. Valid for arbitrary real matrices
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine Svdcmp_dp(a,w,v)
    !!--++    real(dp),dimension(:,:),intent(in out) :: a  !A(m,n)
    !!--++    real(dp),dimension(:),  intent(   out) :: w  !W(n)
    !!--++    real(dp),dimension(:,:),intent(   out) :: v  !V(n,n)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Given an M �N matrix A ,this routine computes its singular value decomposition,
    !!--++    A = U �W �VT . The matrix U replaces A on output. The diagonal matrix of
    !!--++    singular values W is output as the N-dimensional vector w. The N�N matrix V
    !!--++    (not the transpose VT )is output as v .
    !!--++    Adapted from Numerical Recipes. Valid for arbitrary real matrices
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Svdcmp_dp(a,w,v)
       !---- Arguments ----!
       real(kind=dp),dimension(:,:),intent(in out) ::a
       real(kind=dp),dimension(:),  intent(   out) ::w
       real(kind=dp),dimension(:,:),intent(   out) ::v

       !---- Local variables ----!
       integer, parameter                          :: num_its=500
       integer                                     ::i,its,j,k,l,m,n,nm
       real(kind=dp)                               ::anorm,c,f,g,h,s,scal,x,y,z
       real(kind=dp),dimension(size(a,1))          ::tempm
       real(kind=dp),dimension(size(a,2))          ::rv1,tempn

       m=size(a,1)
       n=size(a,2)
       call init_err_mathgen()
       if ( .not. (size(v,1) == n .and. size(v,2) == n .and. size(w) == n)) then
          ERR_MathGen = .true.
          ERR_MathGen_Mess = " => Physical dimensions of arguments in SVDcmp_dp are not compatible "
          return
       end if
       g=0.0_dp
       scal=0.0_dp
       do i=1,n
          l=i+1
          rv1(i)=scal*g
          g=0.0_dp
          scal=0.0_dp
          if (i <=m)then
             scal=sum(abs(a(i:m,i)))
             if ( abs(scal) > tiny(1.0_dp) ) then
                a(i:m,i)=a(i:m,i)/scal
                s=dot_product(a(i:m,i),a(i:m,i))
                f=a(i,i)
                g=-sign(sqrt(s),f)
                h=f*g-s
                a(i,i)=f-g
                tempn(l:n)=matmul(a(i:m,i),a(i:m,l:n))/h
                a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
                a(i:m,i)=scal*a(i:m,i)
             end if
          end if
          w(i)=scal*g
          g=0.0_dp
          scal=0.0_dp
          if ((i <=m).and.(i /=n))then
             scal=sum(abs(a(i,l:n)))
             if ( abs(scal) > tiny(1.0_dp) ) then
                a(i,l:n)=a(i,l:n)/scal
                s=dot_product(a(i,l:n),a(i,l:n))
                f=a(i,l)
                g=-sign(sqrt(s),f)
                h=f*g-s
                a(i,l)=f-g
                rv1(l:n)=a(i,l:n)/h
                tempm(l:m)=matmul(a(l:m,l:n),a(i,l:n))
                a(l:m,l:n)=a(l:m,l:n)+outerprod(tempm(l:m),rv1(l:n))
                a(i,l:n)=scal*a(i,l:n)
             end if
          end if
       end do
       anorm=maxval(abs(w)+abs(rv1))
       do i=n,1,-1
          if (i <n) then
             if ( abs(g) > tiny(1.0_dp) ) then
                v(l:n,i)=(a(i,l:n)/a(i,l))/g
                tempn(l:n)=matmul(a(i,l:n),v(l:n,l:n))
                v(l:n,l:n)=v(l:n,l:n)+outerprod(v(l:n,i),tempn(l:n))
             end if
             v(i,l:n)=0.0_dp
             v(l:n,i)=0.0_dp
          end if
          v(i,i)=1.0_dp
          g=rv1(i)
          l=i
       end do
       do i=min(m,n),1,-1
          l=i+1
          g=w(i)
          a(i,l:n)=0.0_dp
          if ( abs(g) > tiny(1.0_dp) ) then
             g=1.0_dp/g
             tempn(l:n)=(matmul(a(l:m,i),a(l:m,l:n))/a(i,i))*g
             a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
             a(i:m,i)=a(i:m,i)*g
          else
             a(i:m,i)=0.0_dp
          end if
          a(i,i)=a(i,i)+1.0_dp
       end do
       do k=n,1,-1
          do its=1,num_its
             do l=k,1,-1
                nm=l-1
                if ((abs(rv1(l))+anorm)==anorm) exit
                if ((abs(w(nm))+anorm)==anorm) then
                   c=0.0_dp
                   s=1.0_dp
                   do i=l,k
                      f=s*rv1(i)
                      rv1(i)=c*rv1(i)
                      if ((abs(f)+anorm)==anorm)exit
                      g=w(i)
                      h=pythag(f,g)
                      w(i)=h
                      h=1.0_dp/h
                      c=(g*h)
                      s=-(f*h)
                      tempm(1:m)=a(1:m,nm)
                      a(1:m,nm)=a(1:m,nm)*c+a(1:m,i)*s
                      a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
                   end do
                   exit
                end if
             end do
             z=w(k)
             if (l ==k)then
                if (z <0.0_dp)then
                   w(k)=-z
                   v(1:n,k)=-v(1:n,k)
                end if
                exit
             end if
             if (its == num_its) then
                ERR_MathGen = .true.
                ERR_MathGen_Mess = " => SVDcmp_dp: convergence not reached ! "
                return
             end if
             x=w(l)
             nm=k-1
             y=w(nm)
             g=rv1(nm)
             h=rv1(k)
             f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0_dp*h*y)
             g=pythag(f,1.0_dp)
             f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
             c=1.0_dp
             s=1.0_dp
             do j=l,nm
                i=j+1
                g=rv1(i)
                y=w(i)
                h=s*g
                g=c*g
                z=pythag(f,h)
                rv1(j)=z
                c=f/z
                s=h/z
                f=(x*c)+(g*s)
                g=-(x*s)+(g*c)
                h=y*s
                y=y*c
                tempn(1:n)=v(1:n,j)
                v(1:n,j)=v(1:n,j)*c+v(1:n,i)*s
                v(1:n,i)=-tempn(1:n)*s+v(1:n,i)*c
                z=pythag(f,h)
                w(j)=z
                if ( abs(z) > tiny(1.0_dp) ) then
                   z=1.0_dp/z
                   c=f*z
                   s=h*z
                end if
                f=(c*g)+(s*y)
                x=-(s*g)+(c*y)
                tempm(1:m)=a(1:m,j)
                a(1:m,j)=a(1:m,j)*c+a(1:m,i)*s
                a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
             end do
             rv1(l)=0.0_dp
             rv1(k)=f
             w(k)=x
          end do
       end do

       return
    End Subroutine Svdcmp_dp

    !!--++
    !!--++ Subroutine Svdcmp_sp(a,w,v)
    !!--++    real(sp),dimension(:,:),intent(in out) :: a  !A(m,n)
    !!--++    real(sp),dimension(:),  intent(   out) :: w  !W(n)
    !!--++    real(sp),dimension(:,:),intent(   out) :: v  !V(n,n)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Given an M �N matrix A ,this routine computes its singular value decomposition,
    !!--++    A = U �W �VT . The matrix U replaces A on output. The diagonal matrix of
    !!--++    singular values W is output as the N-dimensional vector w. The N�N matrix V
    !!--++    (not the transpose VT )is output as v .
    !!--++    Adapted from Numerical Recipes. Valid for arbitrary real matrices
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Svdcmp_sp(a,w,v)
       !---- Arguments ----!
       real(kind=sp),dimension(:,:),intent(in out) :: a
       real(kind=sp),dimension(:),  intent(   out) :: w
       real(kind=sp),dimension(:,:),intent(   out) :: v

       !---- Local variables ----!
       integer, parameter                          :: num_its=500
       integer                                     ::i,its,j,k,l,m,n,nm
       real(kind=sp)                               ::anorm,c,f,g,h,s,scala,x,y,z
       real(kind=sp),dimension(size(a,1))          ::tempm
       real(kind=sp),dimension(size(a,2))          ::rv1,tempn


       m=size(a,1)
       n=size(a,2)
       call init_err_mathgen()
       if ( .not. (size(v,1) == n .and. size(v,2) == n .and. size(w) == n)) then
          ERR_MathGen = .true.
          ERR_MathGen_Mess = " => Physical dimensions of arguments in SVDcmp_sp are not compatible "
          return
       end if
       g=0.0
       scala=0.0
       do i=1,n                        !Householder reduction to bidiagonal form.
          l=i+1
          rv1(i)=scala*g
          g=0.0
          scala=0.0
          if (i <=m)then
             scala=sum(abs(a(i:m,i)))
             if (abs(scala) > tiny(1.0_sp))then
                a(i:m,i)=a(i:m,i)/scala
                s=dot_product(a(i:m,i),a(i:m,i))
                f=a(i,i)
                g=-sign(sqrt(s),f)
                h=f*g-s
                a(i,i)=f-g
                tempn(l:n)=matmul(a(i:m,i),a(i:m,l:n))/h
                a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
                a(i:m,i)=scala*a(i:m,i)
             end if
          end if
          w(i)=scala*g
          g=0.0
          scala=0.0
          if ((i <=m).and.(i /=n))then
             scala=sum(abs(a(i,l:n)))
             if (abs(scala) > tiny(1.0_sp))then
                a(i,l:n)=a(i,l:n)/scala
                s=dot_product(a(i,l:n),a(i,l:n))
                f=a(i,l)
                g=-sign(sqrt(s),f)
                h=f*g-s
                a(i,l)=f-g
                rv1(l:n)=a(i,l:n)/h
                tempm(l:m)=matmul(a(l:m,l:n),a(i,l:n))
                a(l:m,l:n)=a(l:m,l:n)+outerprod(tempm(l:m),rv1(l:n))
                a(i,l:n)=scala*a(i,l:n)
             end if
          end if
       end do
       anorm=maxval(abs(w)+abs(rv1))
       do i=n,1,-1                    ! Accumulation of right-hand transformations.
          if (i <n)then
             if (abs(g) > tiny(1.0_sp))then
                v(l:n,i)=(a(i,l:n)/a(i,l))/g   !Double division to avoid possible underflow.
                tempn(l:n)=matmul(a(i,l:n),v(l:n,l:n))
                v(l:n,l:n)=v(l:n,l:n)+outerprod(v(l:n,i),tempn(l:n))
             end if
             v(i,l:n)=0.0
             v(l:n,i)=0.0
          end if
          v(i,i)=1.0
          g=rv1(i)
          l=i
       end do
       do i=min(m,n),1,-1  !Accumulation of left-hand transformations.
          l=i+1
          g=w(i)
          a(i,l:n)=0.0
          if (abs(g) > tiny(1.0_sp))then
             g=1.0_sp/g
             tempn(l:n)=(matmul(a(l:m,i),a(l:m,l:n))/a(i,i))*g
             a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
             a(i:m,i)=a(i:m,i)*g
          else
             a(i:m,i)=0.0
          end if
          a(i,i)=a(i,i)+1.0_sp
       end do
       do k=n,1,-1           !Diagonalization of the idiagonal form:Loop over
          do its=1,num_its    !singular values,and over allowed iterations.
             do l=k,1,-1      !Test for splitting.
                nm=l-1        !Note that rv1(1)is always zero,so can never fall through bottom of loop.
                if ((abs(rv1(l))+anorm)==anorm) exit
                if ((abs(w(nm))+anorm)==anorm) then
                   c=0.0       ! Cancellation of rv1(l),if l >1 .
                   s=1.0
                   do i=l,k
                      f=s*rv1(i)
                      rv1(i)=c*rv1(i)
                      if ((abs(f)+anorm)==anorm)exit
                      g=w(i)
                      h=pythag(f,g)
                      w(i)=h
                      h=1.0_sp/h
                      c=(g*h)
                      s=-(f*h)
                      tempm(1:m)=a(1:m,nm)
                      a(1:m,nm)=a(1:m,nm)*c+a(1:m,i)*s
                      a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
                   end do
                   exit
                end if
             end do
             z=w(k)
             if (l ==k) then    !Convergence.
                if (z <0.0)then !Singular value is made nonnegative.
                   w(k)=-z
                   v(1:n,k)=-v(1:n,k)
                end if
                exit
             end if
             if (its == num_its) then
                ERR_MathGen = .true.
                ERR_MathGen_Mess = " => SVDcmp_sp: convergence not reached ! "
                return
             end if
             x=w(l)             !Shift from ottom 2-y-2 minor.
             nm=k-1
             y=w(nm)
             g=rv1(nm)
             h=rv1(k)
             f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0_sp*h*y)
             g=pythag(f,1.0_sp)
             f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
             c=1.0  ! Next QR transformation:
             s=1.0
             do j=l,nm
                i=j+1
                g=rv1(i)
                y=w(i)
                h=s*g
                g=c*g
                z=pythag(f,h)
                rv1(j)=z
                c=f/z
                s=h/z
                f=(x*c)+(g*s)
                g=-(x*s)+(g*c)
                h=y*s
                y=y*c
                tempn(1:n)=v(1:n,j)
                v(1:n,j)=v(1:n,j)*c+v(1:n,i)*s
                v(1:n,i)=-tempn(1:n)*s+v(1:n,i)*c
                z=pythag(f,h)
                w(j)=z                 !Rotation can e arbitrary if z =0 .
                if (abs(z) > tiny(1.0_sp) )then
                   z=1.0_sp/z
                   c=f*z
                   s=h*z
                end if
                f=(c*g)+(s*y)
                x=-(s*g)+(c*y)
                tempm(1:m)=a(1:m,j)
                a(1:m,j)=a(1:m,j)*c+a(1:m,i)*s
                a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
             end do
             rv1(l)=0.0
             rv1(k)=f
             w(k)=x
          end do
       end do

       return
    End Subroutine Svdcmp_sp

    !!----
    !!---- Subroutine Swap(a,b) or Swap(a,b,mask)
    !!----    integer,real(cp),complex, intent( in out) :: a, b
    !!----      or
    !!----    integer,real(cp),complex, dimension(:), intent( in out) :: a, b
    !!----      or
    !!----    integer,real(cp),complex, dimension(:,:), intent( in out) :: a, b
    !!----      or
    !!----    real(kind=cp),  intent(in out) :: a,b
    !!----    logical,        intent(in)     :: mask
    !!----      or
    !!----    real(kind=cp), dimension(:), intent(in out) :: a,b
    !!----    logical,       dimension(:), intent(in)     :: mask
    !!----      or
    !!----    real(kind=cp), dimension(:,:), intent(in out) :: a,b
    !!----    logical,       dimension(:,:), intent(in)     :: mask
    !!----
    !!----    Swap the contents of a and b, when mask (if given) is true.
    !!----
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine Swap_C(a,b)
    !!--++    complex, intent(in out) :: a,b
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Swap_C(a,b)
       !---- Arguments ----!
       complex, intent(in out) :: a
       complex, intent(in out) :: b

       !---- Local variables ----!
       complex :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_C

    !!--++
    !!--++ Subroutine Swap_Cm(A,B)
    !!--++    complex, dimension(:,:), intent(in out) :: a,b
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Swap_Cm(a,b)
       !---- Arguments ----!
       complex, dimension(:,:), intent(in out) :: a
       complex, dimension(:,:), intent(in out) :: b

       !---- Local variables ----!
       complex, dimension(size(a,1),size(a,2)) :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_Cm

    !!--++
    !!--++ Subroutine Swap_Cv(a,b)
    !!--++    complex, dimension(:), intent(in out) :: a,b
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Swap_Cv(a,b)
       !---- Arguments ----!
       complex, dimension(:), intent(in out) :: a
       complex, dimension(:), intent(in out) :: b

       !---- Local variables ----!
       complex, dimension(size(a)) :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_Cv

    !!--++
    !!--++ Subroutine Swap_I(A,B)
    !!--++    integer , intent(in out) :: a,b
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Swap_I(A,B)
       !---- Arguments ----!
       integer , intent(in out) :: a
       integer , intent(in out) :: b

       !---- Local variables ----!
       integer  :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_I

    !!--++
    !!--++ Subroutine Swap_Im(A,B)
    !!--++    integer, dimension(:,:), intent(in out) :: a,b
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Swap_Im(A,B)
       !---- Arguments ----!
       integer, dimension(:,:), intent(in out) :: a
       integer, dimension(:,:), intent(in out) :: b

       !---- Local Variables ----!
       integer, dimension(size(a,1),size(a,2)) :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_Im

    !!--++
    !!--++ Subroutine Swap_Iv(A,B)
    !!--++    integer, dimension(:), intent(in out) :: a,b
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Swap_Iv(A,B)
       !---- Arguments ----!
       integer, dimension(:), intent(in out) :: a
       integer, dimension(:), intent(in out) :: b

       !---- Local Variables ----!
       integer, dimension(size(a)) :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_Iv

    !!--++
    !!--++ Subroutine Swap_R(A,B)
    !!--++    real(kind=cp) , intent(in out) :: a,b
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Swap_R(A,B)
       !---- Arguments ----!
       real(kind=cp), intent(in out) :: a
       real(kind=cp), intent(in out) :: b

       !---- Local variables ----!
       real(kind=cp) :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_R

    !!--++
    !!--++ Subroutine Swap_Rm(A,B)
    !!--++    real(kind=cp), dimension(:,:), intent(in out) :: a,b
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Swap_Rm(A,B)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:), intent(in out) :: a
       real(kind=cp), dimension(:,:), intent(in out) :: b

       !---- Local variables ----!
       real(kind=cp), dimension(size(a,1),size(a,2)) :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_Rm

    !!--++
    !!--++ Subroutine Swap_Rv(A,B)
    !!--++    real(kind=cp), dimension(:), intent(in out) :: a,b
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Swap_Rv(A,B)
       !---- Arguments ----!
       real(kind=cp), dimension(:), intent(in out) :: a
       real(kind=cp), dimension(:), intent(in out) :: b

       !---- Local variables ----!
       real(kind=cp), dimension(size(a)) :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_Rv

    !!--++
    !!--++ Subroutine Masked_Swap_R(A,B,Mask)
    !!--++    real(kind=cp), intent(in out) :: a,b
    !!--++    logical,           intent(in) :: mask
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b if mask=.true.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Masked_Swap_R(A,B,Mask)
       !---- Arguments ----!
       real(kind=cp), intent(in out) :: a
       real(kind=cp), intent(in out) :: b
       logical,           intent(in) :: mask

       !---- Local Variables ----!
       real(kind=cp) :: swp

       if (mask) then
          swp=a
          a=b
          b=swp
       end if

       return
    End Subroutine Masked_Swap_R

    !!--++
    !!--++ Subroutine Masked_Swap_Rm(A,B,Mask)
    !!--++    real(kind=cp), dimension(:,:),intent(in out) :: a,b
    !!--++    logical,       dimension(:,:),    intent(in) :: mask
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b where mask=.true.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Masked_Swap_Rm(A,B,Mask)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:), intent(in out) :: a
       real(kind=cp), dimension(:,:), intent(in out) :: b
       logical,       dimension(:,:), intent(in)     :: mask

       !---- Local variables ----!
       real(kind=cp), dimension(size(a,1),size(a,2)) :: swp

       where (mask)
          swp=a
          a=b
          b=swp
       end where

       return
    End Subroutine Masked_Swap_Rm

    !!--++
    !!--++ Subroutine Masked_Swap_Rv(A,B,Mask)
    !!--++    real(kind=cp), dimension(:),intent(in out) :: a,b
    !!--++    logical,       dimension(:),    intent(in) :: mask
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b where mask=.true.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Masked_Swap_Rv(A,B,Mask)
       !---- Arguments ----!
       real(kind=cp), dimension(:), intent(in out) :: a
       real(kind=cp), dimension(:), intent(in out) :: b
       logical,       dimension(:), intent(in)     :: mask

       !---- Local variables ----!
       real(kind=cp), dimension(size(a))           :: swp

       where (mask)
          swp=a
          a=b
          b=swp
       end where

       return
    End Subroutine Masked_Swap_Rv

    !!--++
    !!--++ Subroutine Tqli1(d,e,n)
    !!--++    real(kind=cp), dimension(:), intent (in out):: d
    !!--++    real(kind=cp), dimension(:), intent (in out):: e
    !!--++    integer,                     intent (in)    :: n
    !!--++
    !!--++    (PRIVATE)
    !!--++    QL-algorithm with implicit shifts, to determine the eigenvalues
    !!--++    and eigenvectors of a real tridiagonal symmetric matrix, or of
    !!--++    a real symmetric matrix previously reduced by tred. D is a vector
    !!--++    with the diagonal elements of the tridiagonal matrix. on output
    !!--++    it returns the eigenvalues. the vector e inputs the subdiagonal
    !!--++    elements of the tridiagonal matrix, with E(1) arbitrary. on
    !!--++    output e is destroyed.
    !!--++    In TLQ1 only the eigenvalues are calculated
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Tqli1(d,e,n)
       !---- Arguments ----!
       real(kind=cp), dimension(:), intent(in out):: d, e ! d(np),e(np)
       integer,                     intent(in )   :: n

       !---- Local variables ----!
       integer      :: i, iter, l, m, mv
       real(kind=cp):: b, c, dd, f, g, p, r, s, comp

       call init_Err_MathGen()
       do i=2,n
          e(i-1)=e(i)
       end do
       e(n)=0.0
       do l=1,n
          iter=0
          do_g : do
             mv=n
             do m=l,n-1
                dd=abs(d(m))+abs(d(m+1))
                comp= abs(e(m))+dd
                if (abs(comp-dd) <= ep_ss) then
                   mv=m
                   exit
                end if
             end do
             m=mv

             if (m /= l) then
                if (iter == 40) then
                   ERR_MathGen=.true.
                   ERR_MathGen_Mess=" Too many iterations in TQLI1"
                   exit
                end if

                iter=iter+1
                g=(d(l+1)-d(l))/(2.0*e(l))
                r=sqrt(g*g+1.0)
                g=d(m)-d(l)+e(l)/(g+sign(r,g))
                s=1.0
                c=1.0
                p=0.0
                do i=m-1,l,-1
                   f=s*e(i)
                   b=c*e(i)
                   r=sqrt(f*f+g*g)
                   e(i+1)=r
                   if (abs(r)  <= ep_ss) then
                      d(i+1)=d(i+1)-p
                      e(m)=0.0
                      cycle do_g
                   end if
                   s=f/r
                   c=g/r
                   g=d(i+1)-p
                   r=(d(i)-g)*s+2.0*c*b
                   p=s*r
                   d(i+1)=g+p
                   g=c*r-b
                end do
                d(l)=d(l)-p
                e(l)=g
                e(m)=0.0
                cycle do_g
             end if
             exit
          end do do_g
       end do

       return
    End Subroutine Tqli1

    !!--++
    !!--++ Subroutine Tqli2(d,e,n,z)
    !!--++    real(kind=cp), dimension(:)  , intent (in out):: d
    !!--++    real(kind=cp), dimension(:)  , intent (in out):: e
    !!--++    integer,                       intent (in)    :: n
    !!--++    real(kind=cp), dimension(:,:), intent (in out):: z
    !!--++
    !!--++    (PRIVATE)
    !!--++    QL-algorithm with implicit shifts, to determine the eigenvalues
    !!--++    and eigenvectors of a real tridiagonal symmetric matrix, or of
    !!--++    a real symmetric matrix previously reduced by tred. D is a vector
    !!--++    with the diagonal elements of the tridiagonal matrix. on output
    !!--++    it returns the eigenvalues. the vector e inputs the subdiagonal
    !!--++    elements of the tridiagonal matrix, with E(1) arbitrary. on
    !!--++    output e is destroyed.
    !!--++    The eigenvectors of the tridiagonal matrix are calculated in TLQ2
    !!--++    by providing the matrix Z  as the identity matrix on input. if the
    !!--++    eigenvectors of the matrix reduced by tred are required, then Z
    !!--++    is input as the matrix output of tred. in either cased, the k-th
    !!--++    column of Z returns the mormalized eigenvector corresponding to
    !!--++    D(k).
    !!--++
    !!--++  Update: February - 2005
    !!
    Subroutine Tqli2(d,e,n,z)
       !---- Arguments ----!
       real(kind=cp), dimension(:),   intent(in out) :: d, e ! d(np),e(np)
       integer,                       intent(in )    :: n
       real(kind=cp), dimension(:,:), intent(in out) :: z    ! z(np,np)

       !---- Local Variables ----!
       integer       :: i, iter, k, l, m, mv
       real(kind=cp) :: b, c, dd, f, g, p, r, s, comp

       call init_Err_MathGen()
       do i=2,n
          e(i-1)=e(i)
       end do

       e(n)=0.0
       do l=1,n
          iter=0
          do_g: do
             mv=n
             do m=l,n-1
                dd=abs(d(m))+abs(d(m+1))
                comp= abs(e(m))+dd
                if (abs(comp-dd) <= ep_ss) then
                   mv=m
                   exit
                end if
             end do
             m=mv
             if (m /= l) then
                if (iter == 40) then
                   ERR_MathGen=.true.
                   ERR_MathGen_Mess=" Too many iterations in TQLI2"
                   exit
                end if

                iter=iter+1
                g=(d(l+1)-d(l))/(2.0*e(l))
                r=sqrt(g*g+1.0)
                g=d(m)-d(l)+e(l)/(g+sign(r,g))
                s=1.0
                c=1.0
                p=0.0
                do i=m-1,l,-1
                   f=s*e(i)
                   b=c*e(i)
                   r=sqrt(f*f+g*g)
                   e(i+1)=r
                   if (abs(r) <= ep_ss) then
                      d(i+1)=d(i+1)-p
                      e(m)=0.0
                      cycle do_g
                   end if
                   s=f/r
                   c=g/r
                   g=d(i+1)-p
                   r=(d(i)-g)*s+2.0*c*b
                   p=s*r
                   d(i+1)=g+p
                   g=c*r-b

                   !---- omit lines from here ...
                   do k=1,n
                      f=z(k,i+1)
                      z(k,i+1)=s*z(k,i)+c*f
                      z(k,i)=c*z(k,i)-s*f
                   end do

                   !---- ... to here when finding only eigenvalues.
                end do
                d(l)=d(l)-p
                e(l)=g
                e(m)=0.0
                cycle do_g
             end if
             exit
          end do do_g
       end do

       return
    End Subroutine Tqli2

    !!--++
    !!--++ Subroutine Tred1(a,n,d,e)
    !!--++    real(kind=cp), dimension(:,:), intent (in out):: a
    !!--++    integer,                       intent (in)    :: n
    !!--++    real(kind=cp), dimension(:)  , intent (in out):: d
    !!--++    real(kind=cp), dimension(:)  , intent (in out):: e
    !!--++
    !!--++    (PRIVATE)
    !!--++    Subroutine for preparing the matrix to find only eigenvalues
    !!--++    Householder reduction of a real symetric nxn matrix A.
    !!--++    On output A is replaced by the orthogonal matrix Q effecting
    !!--++    the transformation. D returns the diagonal elements of the tri-
    !!--++    diagonal matrix and E the off-diagonal elements with E(1)=0.
    !!--++    In tred1 several lines have been deleted and A contains no
    !!--++    useful information on output.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Tred1(a,n,d,e)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:), intent(in out) :: a    ! a(np,np)
       integer,                       intent(in)     :: n
       real(kind=cp), dimension(:),   intent(in out) :: d, e ! d(np),e(np)

       !---- Local Variables ----!
       integer :: i, j, k, l
       real(kind=cp)    :: f, g, h, hh, scala

       do i=n,2,-1
          l=i-1
          h=0.0
          scala=0.0
          if (l > 1)then
             do k=1,l
                scala=scala+abs(a(i,k))
             end do
             if (abs(scala) <= ep_ss) then
                e(i)=a(i,l)
             else
                do k=1,l
                   a(i,k)=a(i,k)/scala
                   h=h+a(i,k)**2
                end do
                f=a(i,l)
                g=-sign(sqrt(h),f)
                e(i)=scala*g
                h=h-f*g
                a(i,l)=f-g
                f=0.0
                do j=1,l
                   g=0.0
                   do k=1,j
                      g=g+a(j,k)*a(i,k)
                   end do
                   do k=j+1,l
                      g=g+a(k,j)*a(i,k)
                   end do
                   e(j)=g/h
                   f=f+e(j)*a(i,j)
                end do
                hh=f/(h+h)
                do j=1,l
                   f=a(i,j)
                   g=e(j)-hh*f
                   e(j)=g
                   do k=1,j
                      a(j,k)=a(j,k)-f*e(k)-g*a(i,k)
                   end do
                end do
             end if
          else
             e(i)=a(i,l)
          end if
          d(i)=h
       end do

       e(1)=0.0
       do i=1,n
          d(i)=a(i,i)
       end do

       return
    End Subroutine Tred1

    !!--++
    !!--++ Subroutine Tred2(a,n,d,e)
    !!--++    real(kind=cp), dimension(:,:), intent (in out) :: a
    !!--++    integer,                       intent (in)     :: n
    !!--++    real(kind=cp), dimension(:)  , intent (in out) :: d
    !!--++    real(kind=cp), dimension(:)  , intent (in out) :: e
    !!--++
    !!--++    (PRIVATE)
    !!--++    Subroutine for preparing the matrix to find the complete set
    !!--++    of eigenvectors.
    !!--++    Householder reduction of a real symetric nxn matrix A.
    !!--++    On output A is replaced by the orthogonal matrix Q effecting
    !!--++    the transformation. D returns the diagonal elements of the tri-
    !!--++    diagonal matrix and E the off-diagonal elements with E(1)=0.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Tred2(a,n,d,e)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:), intent(in out) :: a    ! a(np,np)
       integer,                       intent(in)     :: n
       real(kind=cp), dimension(:),   intent(in out) :: d, e ! d(np),e(np)

       !---- Local variables ----!
       integer :: i, j, k, l
       real(kind=cp)    :: f, g, h, hh, scala

       do i=n,2,-1
          l=i-1
          h=0.0
          scala=0.0
          if (l > 1)then
             do k=1,l
                scala=scala+abs(a(i,k))
             end do
             if (abs(scala) <= ep_ss) then
                e(i)=a(i,l)
             else
                do k=1,l
                   a(i,k)=a(i,k)/scala
                   h=h+a(i,k)**2
                end do
                f=a(i,l)
                g=-sign(sqrt(h),f)
                e(i)=scala*g
                h=h-f*g
                a(i,l)=f-g
                f=0.0
                do j=1,l
                   !---- omit following line if finding only eigenvalues
                   a(j,i)=a(i,j)/h
                   g=0.0
                   do k=1,j
                      g=g+a(j,k)*a(i,k)
                   end do
                   do k=j+1,l
                      g=g+a(k,j)*a(i,k)
                   end do
                   e(j)=g/h
                   f=f+e(j)*a(i,j)
                end do
               hh=f/(h+h)
                do j=1,l
                   f=a(i,j)
                   g=e(j)-hh*f
                   e(j)=g
                   do k=1,j
                      a(j,k)=a(j,k)-f*e(k)-g*a(i,k)
                   end do
                end do
             end if
          else
             e(i)=a(i,l)
          end if
          d(i)=h
       end do

       !---- omit following line if finding only eigenvalues.
       d(1)=0.0
       e(1)=0.0
       do i=1,n
          !---- delete lines from here ...
          l=i-1
          if (abs(d(i)) > ep_ss)then
             do j=1,l
                g=0.0
                do k=1,l
                   g=g+a(i,k)*a(k,j)
                end do
                do k=1,l
                   a(k,j)=a(k,j)-g*a(k,i)
                end do
             end do
          end if
          !---- ... to here when finding only eigenvalues.
          d(i)=a(i,i)
          !---- also delete lines from here ...
          a(i,i)=1.0
          do j=1,l
             a(i,j)=0.0
             a(j,i)=0.0
          end do
          !---- ... to here when finding only eigenvalues.
       end do

       return
    End Subroutine Tred2

 End Module CFML_Math_General

