                       
                       
Welcome to the iFit/iData package
=================================

                              <ifit.mccode.org>
                              
                        E. Farhi, ILL/CS <farhi@ill.fr>
                           Version @IFIT_VERSION@ - @IFIT_DATE@

Purpose
-------

  This library aims at providing basic functionality to achieve some of the 
  general tasks needed for scientific data analysis:
    Load, Plot, Save, Fit, Math operations, define and use Models
    
  It also includes specific 'applications' for neutron scattering:
  * instrument simulation and optimisation using McStas <http://www.mcstas.org>
  * neutron scattering triple-axis spectrometer (TAS) resolution calculation 'ResLibCal'
  * neutron scattering time-of-flight spectroscopy S(q,w) analysis
  * lattice dynamics 4D S(q,w) computation using ASE and DFT codes, which can be coupled to ResLibCal.

Requirements
------------

**Requirements for standalone (binary) package: NONE**
  Stand-alone versions do not require a Matlab license, and have no dependency
  except the Matlab Compiler Runtime (MCR). This latter is included in some of the 
  packages, or have to be installed separately.
  You can get the MCR installer at 
    http://ifit.mccode.org/Downloads/binary

**Requirements for Source package**
  Matlab (any version from 7.x, R2007a), possibly a C compiler for the looktxt  
  and the cbf_uncompress MeX.

Installation
------------

**Installation from Source package**
  Copy the library directories where-ever you want or in MALTAB/toolbox/local:
    /home/joe/Matlab/iFit
  or
    /opt/MATLAB/R2010a/toolbox/local/iFit
  
  Then start Matlab and type in, e.g.:
    >> addpath(genpath('/home/joe/Matlab/iFit'))
  or
    >> addpath(genpath('/opt/MATLAB/R2010a/toolbox/local/iFit'))
  
**Installation for standalone (binary) package**
  Install the MCR (see above, prefer /opt/MATLAB location on Linux systems), then 
  extract the iFit binary package and launch 'ifit'.

**Quick start**
  type at Matlab prompt:
    >> addpath(genpath('/path/to/iFit'))
    >> doc(iData)

  Then refer to the Quick Start tutorial (in iFit/Docs/QuickStart.html).
  
Useful links
------------
 * Homepage: http://ifit.mccode.org
 * Download: http://ifit.mccode.org/Install.html
 * Asking for help: http://mail.mccode.org/cgi-bin/mailman/listinfo/ifit-users
 * Issue tracking: https://github.com/McStasMcXtrace/iFit/issues
 * Developer site: https://github.com/McStasMcXtrace/iFit
 * Project statistics: https://www.openhub.net/p/ifit

Contacts
--------
  You can register to the iFit mailing list at <http://mail.mccode.org/cgi-bin/mailman/listinfo/ifit-users>
  Send messages to the <ifit-users@mail.mccode.org>.
  Help pages are available at <http://ifit.mccode.org>
  
License: EUPL
-------------

  Basically this is open-source. Use it if you find it useful, and enrich it.
  If you do produce new methods, please send them back to me so that they are 
  added in the software and thus benefit to the community.

  In short, you can use, copy, distribute and modify the Code. However, a number 
  of restrictions apply, especially when producing derived work (that is modify 
  and redistribute the code in other products). In particular, the derived work 
  must be licensed under the EUPL or a Compatible License, label all modifications 
  explicitly, distribute the Source Code, and cite the Original work.
  The Source code of iFit is freely available at <http://ifit.mccode.org>

  A number of additions, included in the software, where obtained from the Matlab 
  Central contributions, and are BSD licensed.
  
  Matlab is a registered trademark of The Mathworks Inc.

  Contributions are listed in the Credits page at <http://ifit.mccode.org/Credits.html>

Disclaimer
----------
  This is not a professional tool, and there is no Dev team to actively take care
  of it. Expect occasional failures and bugs. However, I try my best to make the 
  software efficient and reliable. 

--------------------------------------------------------------------------------
$Revision$ - $Date$
