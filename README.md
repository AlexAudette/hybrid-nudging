# Documentation on what files to modify to use the hybrid nudging method
Author: Alexandre Audette


## Source code modifications

CESM versions tested : CESM1.2.2 (others might work, but were not tested)

All the changes to the code are included in the src.cice directory. All you have to do is copy them into the SourceMods/src.cice/ directory in your CESM case directory.

## User namelist files

Add the lines from the user_nl_cice files to the files with the same name in your case directory. This will activate the nudging.

## Build namelist files

In the build_namelists directory in this repo is the namelist definitions for the nudging. CESM doesn't allow (to my knowledge) for namelist definitions for individual cases, so you have to change the source file with that one. The source file is located in

cesm1_2_2/models/ice/cice/bld/namelist_files/

I recommend strongly keeping a copy of the original file for your other cases.

Have fun! :)