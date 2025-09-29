#SPectral ANalysis software (SPAN).
#Written by Daniele Gasparri#

"""
    Copyright (C) 2020-2025, Daniele Gasparri

    E-mail: daniele.gasparri@gmail.com

    SPAN is a GUI software that allows to modify and analyze 1D astronomical spectra.

    1. This software is licensed for non-commercial, academic and personal use only.
    2. The source code may be used and modified for research and educational purposes, 
    but any modifications must remain for private use unless explicitly authorized 
    in writing by the original author.
    3. Redistribution of the software in its original, unmodified form is permitted 
    for non-commercial purposes, provided that this license notice is always included.
    4. Redistribution or public release of modified versions of the source code 
    is prohibited without prior written permission from the author.
    5. Any user of this software must properly attribute the original author 
    in any academic work, research, or derivative project.
    6. Commercial use of this software is strictly prohibited without prior 
    written permission from the author.

    DISCLAIMER:
    THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES, OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT, OR OTHERWISE, ARISING FROM, OUT OF, OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""
# Import local modules needed by SPAN

import importlib
import sys

# Import GUI module
try:
    from FreeSimpleGUI_local import FreeSimpleGUI as sg
except ModuleNotFoundError:
    from span.FreeSimpleGUI_local import FreeSimpleGUI as sg

# ------------------------------------------------------------------
# Dynamic imports (with alias)
# ------------------------------------------------------------------
modules = {
    # Core functions
    "stm": "span_functions.system_span",
    "uti": "span_functions.utilities",
    "spman": "span_functions.spec_manipul",
    "spmt": "span_functions.spec_math",
    "ls": "span_functions.linestrength",
    "span": "span_functions.spec_analysis",
    "cubextr": "span_functions.cube_extract",

    # Modules
    "layouts": "span_modules.layouts",
    "misc": "span_modules.misc",
    "sub_programs": "span_modules.sub_programs",
    "spec_manipulation": "span_modules.spec_manipulation",
    "param_windows": "span_modules.param_windows",
    "files_setup": "span_modules.files_setup",
    "utility_tasks": "span_modules.utility_tasks",
    "apply_spec_tasks": "span_modules.apply_spec_tasks",
    "apply_analysis_tasks": "span_modules.apply_analysis_tasks",
    "check_spec": "span_modules.check_spec",
    "settings": "span_modules.settings",
    "file_writer": "span_modules.file_writer",

    # New imports (with alias for clarity)
    "zoom": "span_modules.ui_zoom",
    "preview_tools": "span_modules.preview_tools",
}

# Import dynamically with fallback
for alias, module in modules.items():
    try:
        imported_module = importlib.import_module(module)
    except ModuleNotFoundError:
        imported_module = importlib.import_module(f"span.{module}")

    globals()[alias] = imported_module

# ------------------------------------------------------------------
# Direct imports (no alias, use full module name)
# ------------------------------------------------------------------
from span_modules import listbox_events
# from span_modules import preview_tools_extra

# ------------------------------------------------------------------
# Data classes / definitions
# ------------------------------------------------------------------
from params import SpectraParams



