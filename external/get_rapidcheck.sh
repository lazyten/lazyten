#!/bin/sh
## ---------------------------------------------------------------------
##
## Copyright (C) 2016 by the linalgwrap authors
##
## This file is part of linalgwrap.
##
## linalgwrap is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published
## by the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## linalgwrap is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with linalgwrap. If not, see <http://www.gnu.org/licenses/>.
##
## ---------------------------------------------------------------------

## DO NOT EDIT
## This file is automatically generated from a file in the repository "krims".
## Edit the original and call the script "update_from_sister_repos.sh" instead.

#
# Settings
#
# The git repo to checkout
FROM="https://github.com/emil-e/rapidcheck"

# Folder to check it out to
WHAT="rapidcheck"

# Interval: How often to update:
INTERVAL="1 hour"

# File to use in order to test a successful checkout
CHECKFILE="CMakeLists.txt"

# Branch to checkout (empty for default)
BRANCH=""

if [ ! -f "$PWD/get.lib.sh" ]; then
	echo "PWD needs to be the location of the get.lib.sh file."
	exit 1
fi

. "$PWD/get.lib.sh" || exit 1
update_repo "$FROM" "$WHAT" "$CHECKFILE" "$INTERVAL" "$BRANCH"
exit $?