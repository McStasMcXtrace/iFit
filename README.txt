                       Welcome to the iFit/iData package
                              <ifit.mccode.org>
                              
                        E. Farhi, ILL/CS <farhi@ill.fr>
                           Version @IFIT_VERSION@ - @IFIT_DATE@

** Purpose:
This small library aims at providing basic functionality to achieve some of the 
general tasks needed for scienfic data analysis:
    Load, Plot, Save, Fit, Math operations

** License: GPL 2
Basically this is open-source. Use it if you find it usefull, and enrich it.
If you do produce new methods, please send them back to me so that they are 
added in the software and thus benefit to the community.

Any commercial use of this library, that is re-use of the code for commercial 
purposes, is prevented (this is GPL 2). Any inclusion into an other software must 
properly include source code, and cite authors.

A number of additions, included in the software, where obtained from the Matlab 
Central contributions, and are BSD licensed.

Contributions are listed in the iFit/Docs/Credits.html page.

** Disclaimer:
This is not a professional tool, and there is no Dev team to actively take care
of it. Expect occasional failures and bugs. However, I try my best to make the 
software efficient and reliable. 

** Requirements:
Matlab (any version from 6.x), possibly a C compiler for the looktxt  and the
cbf_uncompress MeX.

** Installation:
Copy the library directories where-ever you want or in MALTAB/toolbox/local:
  /home/joe/Matlab/iFit
or
  /usr/local/matlab/toolbox/local/iFit
  
Then start Matlab and type in, e.g.:
  >> addpath(genpath('/home/joe/Matlab/iFit'))
or
  >> addpath(genpath('/usr/local/matlab/toolbox/local/iFit'))

** Quick start:
type in:
  >> addpath(genpath('/path/to/iFit'))
  >> doc(iData)

Then refer to the Quick Start tutorial (in iFit/Docs/QuickStart.html).




