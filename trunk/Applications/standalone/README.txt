HOWTO create standalone applications for iFit:

Requires to have Matlab compiler installed (mcc).

Linux (Debian/Ubuntu, works on amd64 and i386):
-----------------------------------------------
  1- navigate into svn 'trunk'
  2- launch: % ./mkdist <major.minor>
will create a src.zip source, a <arch>.zip linux binary and a <arch>.deb

Mac OSX (prefer i386):
----------------------
  1- copy the ifit-src.zip file for the distribution to MacOSX
  2- launch matlab and navigate to ifit-src directory
  3- type: addpath(genpath(pwd))
  4- launch: >> ifitdeploy
will create a 'ifit-maci' binary distribution.

To create a MacOSX App:
----------------------
  1- launch Platypus, set app title to 'iFit'
  2- use the ifit-src/Applications/standalone/macosx/app-script.sh script
  
#!/bin/sh 
# Script to open a terminal which launches iFit

open -a terminal standalone/ifit &

osascript  <<EOF
tell app "Terminal"
  set custom title of front window to "iFit (c) ILL <ifit.mccode.org>"
  set normal text color of front window to "blue"
end tell
EOF

  3- set the icon to ifit-src/Docs/images/iFit-logo.png (drag-n-drop from iFit src directory)
  4- set output to None, version number, allow drag-n-drop, uncheck 'remain running'
  5- drag-n-drop the 'ifit-maci' binary distribution directory renamed as 'standalone'
  6- click create
an 'iFit.app' is created.

To create a MacOSX installer:
-----------------------------
  1- open /osx/Developer/Applications/Utilities/PackageMaker
  2- in Tab configuration, enter Title 'iFit', and a description as 
  ifit-src/Applications/standalone/macosx/description.txt
  
** iFit generic data analysis and fitting to models **
Simple methods to be used for complex data analysis

The iFit program provides a set of methods to load, analyze, plot, fit and optimize models,
 and export results. iFit is based on Matlab, but stand-alone version does not require a
 Matlab license to run. Any text file can be imported straight away, and a set of binary files are
 supported. Any data dimensionality can be handled, including event based data sets.
 .
 The spirit of the software is to include simple object definitions for Data sets and Models, with   
 a set of methods that provide all the means to perform the usual data analysis procedures.
 - iData objects to hold data sets ; Import with:          iData('filename')
 - iFunc objects to hold models ;    Create new ones with: iFunc('expression')
 - fit model to data with: fits(data, model)
 - documentation is available
 
 Main functionalities are: [ iData Load Plot Math Fit Save Optimization iFunc Models ]
 
 To use this software, you need to install the Matlab Compiler Runtime from the DMG available in the ifit.mccode.org website / Download binary section.
 
 To start iFit, start it from the Applications folder.
 
 Refer to <ifit.mccode.org> for on-line documentation.
 Matlab is a registered trademark of The Mathworks Inc.
  
  3- in the Requirements Tab, add a 
    If: file exists on System: /Applications/MATLAB/MATLAB_Compiler_Runtime
    and the failure message: 'MATLAB Compiler Runtime is not installed' 
    'The MATLAB Compiler Runtime is not installed, but it is required for iFit to launch. 
    get it at <http://ifit.mccode.org/Downloads/binary/mac32/MCRInstaller32.dmg>'
    NOTE: this may lead to a JavaScript error. Then remove the requirement.
  4- drag-n-drop the iFit.app in the Contents left panel.
  5- in the Configuration tab, update version, check require admin, uncheck custom location (/Applications)
  6- in Scripts tab, add path to ifit-src/Applications/standalone/macosx/pkg-postscript.sh
  
#!/bin/sh -e
if [ ! -f "`which matlab`"  ]; then
  ln -s /Applications/iFit.app/Contents/MacOS/iFit /usr/bin/matlab
fi
 
  
