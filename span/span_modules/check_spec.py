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

# Functions to check and validate the spectra loaded by the user, either in a spectra list or as a single spectrum

try: #try local import if executed as script
    #GUI import
    from FreeSimpleGUI_local import FreeSimpleGUI as sg
    #SPAN functions import
    from span_functions import system_span as stm
    from params import SpectraParams

except ModuleNotFoundError: #local import if executed as package
    #GUI import
    from span.FreeSimpleGUI_local import FreeSimpleGUI as sg
    #SPAN functions import
    from span.span_functions import system_span as stm
    from .params import SpectraParams

import os
import numpy as np
import pandas as pd

# from params import SpectraParams
from dataclasses import replace

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(CURRENT_DIR)



def load_and_validate_spectra(spectra_list, lambda_units, window):
    """
    Loads and validates a list of spectra, ensuring they exist and are readable.
    """
    fatal_condition = 0  # Flag to indicate critical failure

    # Check if the spectra list file exists
    if not os.path.isfile(spectra_list):
        sg.popup("Error: The spectra file list does not exist. Try again...")
        return 0, [], [], 1  # Fatal condition

    try:
        # Open the file and read all lines as complete strings without splitting by spaces
        with open(spectra_list, "r", encoding="utf-8") as f:
            lines = [line.strip() for line in f if line.strip() and not line.startswith("#")]
    except Exception:
        sg.popup("Cannot read the spectra list. Ensure the format is correct.")
        return 0, [], [], 1  # Fatal condition

    spectra_number = len(lines)

    print(f"You want to load {spectra_number} spectra")
    print("Now checking if they really exist and are valid...")
    print("")

    # Resolve absolute/relative paths
    spec_names = []
    for line in lines:
        if not os.path.isabs(line):
            spec_names.append(os.path.join(BASE_DIR, line))
        else:
            spec_names.append(line)

    # Normalize paths to avoid inconsistencies
    spec_names = [os.path.normpath(spectrum) for spectrum in spec_names]
    spec_names_nopath = [os.path.basename(spectrum) for spectrum in spec_names]

    # Update the listbox in the GUI
    window["-LIST-"].Update(spec_names_nopath)

    # Check if spectra files exist
    spec_not_exist = [s for s in spec_names if not os.path.isfile(s)]
    if spec_not_exist:
        if len(spec_not_exist) == len(spec_names):
            sg.popup("None of the spectra exist. Please check the file list and reload.")
            print("Loaded 0 valid spectra.")
            return 0, [], [], 1  # Fatal condition
        else:
            sg.popup("Some spectra do not exist. Updating the list to include only existing spectra.")
            print(f"The following spectra do not exist: {spec_not_exist}")
            spec_names = np.setdiff1d(spec_names, spec_not_exist).tolist()
            spec_names_nopath = [os.path.basename(s) for s in spec_names]
            window["-LIST-"].Update(spec_names_nopath)
            spectra_number = len(spec_names)

    # Validate the spectra content
    spec_not_readable = []
    spec_not_valid = []

    for i, spectrum in enumerate(spec_names):
        # Show a progress meter if more than 5 spectra
        if len(spec_names) > 5:
            if not sg.OneLineProgressMeter("Reading spectra", i + 1, len(spec_names), "Processing spectra:",
                                           orientation="h", button_color=("white", "red")):
                print("*** CANCELLED ***\n")
                break

        # Check FITS files
        if spectrum.lower().endswith(".fits"):
            valid, message = stm.is_valid_spectrum(spectrum)
            if not valid:
                spec_not_valid.append(spectrum)
        else:
            # Attempt to read non-FITS spectra
            try:
                stm.read_spec(spectrum, lambda_units)
            except Exception as e:
                print(f"Error reading spectrum: {spectrum} | {e}")
                spec_not_readable.append(spectrum)

    # Handle unreadable spectra
    if spec_not_readable:
        if len(spec_not_readable) == len(spec_names):
            sg.popup("None of the spectra are readable. Check the file list and reload.")
            print("Loaded 0 valid spectra.")
            return 0, [], [], 1  # Fatal condition
        else:
            sg.popup("Some spectra are unreadable. Updating the list to exclude them.")
            print(f"Unreadable spectra: {spec_not_readable}")
            spec_names = np.setdiff1d(spec_names, spec_not_readable).tolist()
            spec_names_nopath = [os.path.basename(s) for s in spec_names]
            window["-LIST-"].Update(spec_names_nopath)
            spectra_number = len(spec_names)

    # Handle invalid FITS spectra
    if spec_not_valid:
        if len(spec_not_valid) == len([s for s in spec_names if s.lower().endswith(".fits")]):
            sg.popup("None of the FITS spectra are valid. Check the file list and reload.")
            print("Loaded 0 valid FITS spectra.")
            return 0, [], [], 1  # Fatal condition
        else:
            sg.popup("Some FITS files are invalid. Updating the list to exclude them.")
            print(f"Invalid FITS spectra: {spec_not_valid}")
            spec_names = np.setdiff1d(spec_names, spec_not_valid).tolist()
            spec_names_nopath = [os.path.basename(s) for s in spec_names]
            window["-LIST-"].Update(spec_names_nopath)
            spectra_number = len(spec_names)

    print("")
    print(f"Successfully loaded {spectra_number} valid spectra.")
    print("")
    print("*** Please check the wavelength range using the plot button before proceeding ***")
    print("")

    return spectra_number, spec_names, spec_names_nopath, fatal_condition


def validate_and_load_spectrum(params, window):

    """
    Loads and validates the single spectrum, ensuring it exists and is readable

    """

    # Check if the file exists
    cond00 = os.path.isfile(params.spectra_list)
    params = replace(params, spectra_number=1)
    valid_spec = False
    print('Guessing the type of the spectrum. Is it correct?')

    if not cond00:
        sg.popup("We don't start well: the spectrum does not exist. Try again...")
        return params, valid_spec  # Return params without modifications

    # Assign the selected spectrum
    params = replace(params, prev_spec=params.spectra_list)
    params = replace(params, prev_spec_nopath=os.path.splitext(os.path.basename(params.prev_spec))[0])

    try:
        # Check if the file is a valid FITS spectrum
        if params.prev_spec.lower().endswith('.fits'):
            valid, message = stm.is_valid_spectrum(params.prev_spec)
            print(f"{params.prev_spec}: {message}")  # Log for debugging
            if not valid:
                params = replace(params, spec_not_valid=params.prev_spec)
                sg.popup("Your FITS file does not seem to be a spectrum. Please, load a valid spectrum.")
                params = replace(params, prev_spec='')
                return params, valid_spec  # Exit function with updated params

        # Read the spectrum and update parameters
        wavelength, flux = stm.read_spec(params.prev_spec, params.lambda_units)[:2]
        params = replace(params, wavelength=wavelength, flux=flux)

        # Update filename without path
        params = replace(params, prev_spec_nopath=os.path.splitext(os.path.basename(params.prev_spec))[0])

        # Update GUI Listbox
        spec_name = [params.prev_spec_nopath, ' ']
        window['-LIST-'].Update(spec_name)
        valid_spec = True

    except Exception:
        sg.popup("Oops! Cannot read the spectrum. Are you sure it's a spectrum?")
        print("Please, reload a valid spectrum or a list.")
        valid_spec = False

    return params, valid_spec  # Return updated params
