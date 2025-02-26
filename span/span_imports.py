#SPectral ANalysis software (SPAN).
#Written by Daniele Gasparri#

"""
    Copyright (C) 2020-2025, Daniele Gasparri

    E-mail: daniele.gasparri@gmail.com

    SPAN is a GUI interface that allows to modify and analyse 1D astronomical spectra.

    1. This software is licensed **for non-commercial use only**.
    2. The source code may be **freely redistributed**, but this license notice must always be included.
    3. Any user who redistributes or uses this software **must properly attribute the original author**.
    4. The source code **may be modified** for non-commercial purposes, but any modifications must be clearly documented.
    5. **Commercial use is strictly prohibited** without prior written permission from the author.

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

# List of modules to import dynamically
modules = {
    "stm": "span_functions.system_span",
    "uti": "span_functions.utilities",
    "spman": "span_functions.spec_manipul",
    "spmt": "span_functions.spec_math",
    "ls": "span_functions.linestrength",
    "span": "span_functions.spec_analysis",
    "cubextr": "span_functions.cube_extract",
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
}

# Try importing modules dynamically
for alias, module in modules.items():
    try:
        imported_module = importlib.import_module(module)
    except ModuleNotFoundError:
        imported_module = importlib.import_module(f"span.{module}")

    globals()[alias] = imported_module
