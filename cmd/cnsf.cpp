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
#include "file/matrix.h"
#include "math/SH.h"
#include "math/ZSH.h"
#include <Eigen/src/Core/Matrix.h>

using namespace MR;
using namespace App;

void usage() {
  AUTHOR = "Vashisht Kumar (vashishtkumar2525@yahoo.com) and J-Donald Tournier (jdtournier@gmail.com)";

  SYNOPSIS = "blah";

  DESCRIPTION + "Blah";

  ARGUMENTS
  + Argument("zsh",
            "the input image consisting of zonal harmonics (ZSH) "
            "coefficients of the dMRI signal.")
          .type_image_in()

  + Argument("response odf", "a pair of arguments, consisting of the input response function "
      "and its matching ODF for a given tissue type. Multiple such pairs can be provided to "
      "handle multiple tissue types.")
          .allow_multiple();
}






using value_type = float;



void run() {
  auto zsh_image = Image<value_type>::open(argument[0]);
  if (zsh_image.ndim() < 5)
    throw Exception("ZSH image is expected to have at least 5 dimensions");

  const size_t nbvalues = zsh_image.size(4);
  const size_t nl_zsh = zsh_image.size(3);

  INFO ("input ZSH image contains " + str(nbvalues) + " b-value shells");

  if ((argument.size()-1)&1U)
    throw Exception ("expected responses and ODFs to be provided in matching pairs");

  std::vector<Eigen::MatrixXf> responses;
  std::vector<Image<value_type>> odf_images;
  for (size_t n = 1; n < argument.size(); n+=2) {
    responses.push_back (File::Matrix::load_matrix<float> (argument[n]));
    if (responses.back().rows() != nbvalues)
      throw Exception ("number of b-value shells in response function does not match number in ZSH input");
    odf_images.push_back(Image<value_type>::open(argument[n+1]));
    check_dimensions(zsh_image, odf_images.back(), 0, 3);
  }

  size_t nvoxels = 0;
  for (auto l = Loop(zsh_image, 0, 3)(zsh_image); l; ++l) {
    bool notzero = false;
    for (auto l2 = Loop(zsh_image, 3, 5)(zsh_image); l2; ++l2) {
      if (std::isfinite (zsh_image.value()) && zsh_image.value() != 0) {
        notzero = true;
        break;
      }
    }
    if (notzero)
      nvoxels++;
  }

  const size_t ntissues = odf_images.size();

  INFO("found " + str(nvoxels) + " voxels in ZSH image");
  INFO("assuming " + str(ntissues) + " tissue types:");
  for (size_t n = 0; n < responses.size(); ++n) {
    const int lmax = Math::SH::LforN (odf_images[n].size(3));
    const int nlmax = lmax/2+1;
    INFO("  " + odf_images[n].name() + ": lmax = " + str(lmax));
    if (responses[n].cols() < nlmax)
      throw Exception ("too few harmonic coefficients in response for \"" + odf_images[n].name() + "\"");
    if (responses[n].cols() > nlmax) {
      responses[n].conservativeResize(responses[n].rows(), nlmax);
      DEBUG ("resized response for \"" + odf_images[n].name() + "\" to " + str(responses[n].cols()) + " to match ODF");
    }
  }

  size_t nl = 0;
  for (const auto& R : responses)
    nl = std::max (nl, size_t(R.cols()));

  INFO ("maximum harmonic order is " + str(2*(nl-1)));

  std::vector<Eigen::MatrixXf> zsh (nl, Eigen::MatrixXf::Zero (nvoxels,nbvalues));
  for (size_t n = 0; n < nl; ++n) {
    zsh_image.index(3) = n;

    size_t v = 0;
    Eigen::VectorXf vals;
    for (auto l = Loop ("loading ZSH coefficients for L=" + str(2*n), zsh_image, 0, 3) (zsh_image); l; ++l) {
      vals = zsh_image.row(4);
      if (vals.allFinite() && vals.any()) {
        zsh[n].row(v++) = vals;
      }
    }
  }

}
