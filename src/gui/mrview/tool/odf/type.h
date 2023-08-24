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

#ifndef __gui_mrview_tool_odf_type_h__
#define __gui_mrview_tool_odf_type_h__

#include "gui/dwi/renderer.h"

namespace MR {

class Header;

namespace GUI {
namespace MRView {
namespace Tool {

// TODO Remove
using odf_type_t = GUI::DWI::Renderer::mode_t;

} // namespace Tool
} // namespace MRView
} // namespace GUI
} // namespace MR

#endif
