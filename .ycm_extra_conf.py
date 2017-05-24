## ---------------------------------------------------------------------
##
## Copyright (C) 2016-17 by the lazyten authors
##
## This file is part of lazyten.
##
## lazyten is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published
## by the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## lazyten is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with lazyten. If not, see <http://www.gnu.org/licenses/>.
##
## ---------------------------------------------------------------------

import os
import yaml

# The builddir (relative to the path of this file),
# which we use to find some generated files:
builddir = "build"

# This file is loosely based upon the file
# cpp/ycm/.ycm_extra_conf.py from the youcompleteme daemon process
# available on github:
# https://github.com/Valloric/ycmd/blob/master/cpp/ycm/.ycm_extra_conf.py

# Static flags:
flags = [
    # Warnings: For a very detailed discussion about this
    # see the following stackexchange post:
    # https://programmers.stackexchange.com/questions/122608#124574
    '-Wall',
    '-Wextra',
    '-Wnon-virtual-dtor',
    '-Woverloaded-virtual',
    '-Wold-style-cast',
    '-Wcast-align',
    '-Wconversion',
    '-Wsign-conversion',
    '-pedantic',
    '-Werror',
    # Generate unwind information
    '-fexceptions',
    # Compile debug code as well
    '-DDEBUG=DEBUG',
    # Compile as c++14
    '-std=c++14',
    # Treat .h header files as c++:
    '-x', 'c++',
    # C++14 code blocks in krims
    '-DKRIMS_HAVE_CXX14',
    # Include other libraries and show errors and
    # warnings within them
    # To suppress errors shown here, use "-isystem"
    # instead of "-I"
    '-I', 'src',
    '-I', builddir + '/src',
    #
    '-isystem', 'external/krims/src',
    '-isystem', builddir+'/external/krims/src',
    '-isystem', '../krims/src',
    '-isystem', '../krims/'+builddir+'/src',
    #
    '-isystem', 'external/rapidcheck/include',
    '-isystem', 'external/rapidcheck/ext/catch/include',
    '-isystem', '../rapidcheck/ext/catch/include',
    '-isystem', '../rapidcheck/include',
    #
    '-isystem', '/usr/lib/ycmd/clang_includes',
]

def ParseExtraIncludes():
  thisdir = os.path.dirname( os.path.abspath( __file__ ) )
  extrafile = thisdir + "/"+builddir+"/ycm_extra_includes.yaml"
  if os.path.exists(extrafile):
    with open(extrafile,"r") as f:
      paths = yaml.safe_load(f)
      return [ f for path in paths
               for f in [ "-isystem", path ] ]
  return []

def MakeRelativePathsInFlagsAbsolute( flags, working_directory ):
  if not working_directory:
    return list( flags )
  new_flags = []
  make_next_absolute = False
  path_flags = [ '-isystem', '-I', '-iquote', '--sysroot=' ]
  for flag in flags:
    new_flag = flag

    if make_next_absolute:
      make_next_absolute = False
      if not flag.startswith( '/' ):
        new_flag = os.path.join( working_directory, flag )

    for path_flag in path_flags:
      if flag == path_flag:
        make_next_absolute = True
        break

      if flag.startswith( path_flag ):
        path = flag[ len( path_flag ): ]
        new_flag = path_flag + os.path.join( working_directory, path )
        break

    if new_flag:
      new_flags.append( new_flag )
  return new_flags


#SOURCE_EXTENSIONS = [ '.cpp', '.cxx', '.cc', '.c', '.C' ]
#def IsHeaderFile( filename ):
#  extension = os.path.splitext( filename )[ 1 ]
#  return extension in [ '.h', '.hxx', '.hpp', '.hh' ]

def FlagsForFile( filename, **kwargs ):
  relative_to = os.path.dirname( os.path.abspath( __file__ ) )
  flags_in = flags + ParseExtraIncludes()
  final_flags = MakeRelativePathsInFlagsAbsolute( flags_in, relative_to )

  return {
    'flags': final_flags,
    'do_cache': True
  }
