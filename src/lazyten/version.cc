//
// Copyright (C) 2016-17 by the lazyten authors
//
// This file is part of lazyten.
//
// lazyten is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// lazyten is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with lazyten. If not, see <http://www.gnu.org/licenses/>.
//

#include "lazyten/version.hh"
#include <sstream>

namespace lazyten {

const std::map<std::string, bool> version::feature_availability =
      detail::build_features_map();

std::string version::version_string() {
  std::stringstream ss;
  ss << major << "." << minor << "." << patch;
  return ss.str();
}

bool version::has_feature(const std::string& feature) {
  auto it = feature_availability.find(feature);
  if (it == std::end(feature_availability)) return false;
  return it->second;
}

std::string version::feature_string() {
  std::stringstream ss;
  for (auto it = feature_availability.begin(); it != feature_availability.end(); ++it) {
    if (it != feature_availability.begin()) ss << ' ';
    ss << (it->second ? '+' : '-') << it->first;
  }
  return ss.str();
}

}  // namespace lazyten
