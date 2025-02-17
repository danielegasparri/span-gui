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
