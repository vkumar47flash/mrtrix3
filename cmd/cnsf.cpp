/* Copyright (c) 2008-2023 the MRtrix3 contributors.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * Covered Software is provided under this License on an "as is"
 * basis, without warranty of any kind, either expressed, implied, or
 * statutory, including, without limitation, warranties that the
 * Covered Software is free of defects, merchantable, fit for a
 * particular purpose or non-infringing.
 * See the Mozilla Public License v. 2.0 for more details.
 *
 * For more details, see http://www.mrtrix.org/.
 */

#include "command.h"
#include "debug.h"
#include "exception.h"
#include "image.h"
#include "image_helpers.h"
#include "math/ZSH.h"
#include <Eigen/src/Core/Matrix.h>

using namespace MR;
using namespace App;

void usage() {
  AUTHOR = "Vashisht Kumar (vashishtkumar2525@yahoo.com) and J-Donald Tournier (jdtournier@gmail.com)";

  SYNOPSIS = "blah";

  DESCRIPTION + "Blah";

  ARGUMENTS
  +Argument("zsh",
            "the input image consisting of zonal harmonics (ZSH) "
            "coefficients of the dMRI signal.")
          .type_image_in() +
      Argument("odf", "the input image of SH coefficients of the ODF for each tissue type.")
          .type_image_in()
          .allow_multiple();
}

using value_type = float;

void run() {
  auto zsh_image = Image<value_type>::open(argument[0]);
  if (zsh_image.ndim() < 5)
    throw Exception("ZSH image is expected to have at least 5 dimensions");

  std::vector<Image<value_type>> odf_images;
  for (size_t n = 1; n < argument.size(); ++n) {
    odf_images.push_back(Image<value_type>::open(argument[n]));
    check_dimensions(zsh_image, odf_images.back(), 0, 3);
  }

  size_t nvoxels = 0;
  for (auto l = Loop(zsh_image, 0, 3)(zsh_image); l; ++l) {
    bool notzero = false;
    for (auto l2 = Loop(zsh_image, 3, 5)(zsh_image); l2; ++l2) {
      if (zsh_image.value() != 0) {
        notzero = true;
        break;
      }
    }
    if (notzero)
      nvoxels++;
  }

  const size_t ntissues = odf_images.size();
  const size_t nbvalues = zsh_image.size(4);
  const size_t nl_zsh = zsh_image.size(3);

  INFO("found " + str(nvoxels) + " voxels in ZSH image");
  INFO("assuming " + str(ntissues) + " tissue types");

  std::vector<std::vector<Eigen::MatrixXf>> zsh(nbvalues);
  for (auto &z : zsh)
    z = std::vector<Eigen::MatrixXf>(nl_zsh, Eigen::MatrixXf(nvoxels, ntissues));
}
