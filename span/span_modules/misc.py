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

#Miscellaneous collection of routines used by the GUI

try: #try local import if executed as script
    #GUI import
    from FreeSimpleGUI_local import FreeSimpleGUI as sg
    from span_modules import layouts

except ModuleNotFoundError: #local import if executed as package
    #GUI import
    from span.FreeSimpleGUI_local import FreeSimpleGUI as sg
    from . import layouts

#Python imports
import json
import os
import numpy as np
import urllib.request
import zipfile
import ssl
import certifi


CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(CURRENT_DIR)



def save_config_folder(data, config_file):
    """Save the result path to a JSON file."""
    with open(config_file, 'w') as f:
        json.dump(data, f)


def load_config(config_file):
    """ Load the result path folder from the JSON file."""
    if os.path.exists(config_file):
        with open(config_file, 'r') as f:
            return json.load(f)
    return {}


def ask_user_for_result_path():
    """Ask user to select a folder to store the results of SPAN."""
    layout, scale_win, fontsize, default_size = get_layout()
    sg.theme('LightBlue')
    layout = [
        [sg.Text("Select the path to store the SPAN_results folder:")],
        [sg.InputText(), sg.FolderBrowse()],
        [sg.Button("Confirm"), sg.Button("Cancel")]
    ]
    window = sg.Window("Select folder", layout)

    while True:
        event, values = window.read()
        if event in (sg.WINDOW_CLOSED, "Cancel"):
            window.close()
            return None
        if event == "Confirm" and values[0]:
            window.close()
            return values[0]
        sg.popup("Please, select a valid path")


def create_result_structure(base_path):
    """Creating the directory structure"""
    result_data = os.path.join(base_path, 'SPAN_results')
    subdirectories = [
        'processed_spectra', 'SNR', 'planck_black_body_fitting', 'cross-correlation', 'velocity_dispersion',
        'line-strength_analysis', 'line_fitting', 'stars_and_gas_kinematics', 'stellar_populations_and_sfh',
        'line-strength_sigma_coefficients', 'plots'
    ]
    os.makedirs(result_data, exist_ok=True)
    for subdir in subdirectories:
        os.makedirs(os.path.join(result_data, subdir), exist_ok=True)
    return result_data


def change_result_path(config_folder, config_file):
    """Function to allow the user to change the result directory directly in the GUI"""
    new_path = ask_user_for_result_path()
    if new_path:
        config_folder["result_path"] = new_path
        save_config_folder(config_folder, config_file)
        create_result_structure(new_path)  # Assicura che la struttura sia ricreata
        sg.popup(f"The new SPAN-result folder now is in: {new_path}")
    else:
        sg.popup("Path of the SPAN_result folder has not changed")


def get_layout():
    """Function to select the layout based on the OS"""
    current_os = os.name  # 'posix' for Linux/Mac, 'nt' for Windows
    if current_os == "nt":
        # Adapting to scaling factors on windows, both for the GUI and Matplotlib
        import ctypes
        import matplotlib

        # DPI awareness for windows
        ctypes.windll.shcore.SetProcessDpiAwareness(2)
        ctypes.windll.user32.SetProcessDPIAware()

        # Set the scaling for windows
        dpi_scale = ctypes.windll.shcore.GetScaleFactorForDevice(0) / 100.0

        #with scaling <1.5, do not scale automatically
        if dpi_scale < 1.5:
            scale_win = 1.5
            matplotlib.rcParams['figure.dpi'] = 100
        else:
            scale_win = dpi_scale
            matplotlib.rcParams['figure.dpi'] = 100 * dpi_scale
        fontsize = sg.set_options(font=("Helvetica", 11))
        default_size = 11

        return layouts.layout_windows, scale_win, fontsize, default_size

    elif current_os == "posix":
        # Further check between Linux e macOS
        if "ANDROID_BOOTLOGO" in os.environ:  # Check for Android
            scale_win = 2.25
            fontsize = sg.set_options(font=("Helvetica", 10))
            default_size = 10
            return layouts.layout_android, scale_win, fontsize, default_size
        elif os.uname().sysname == "Darwin":  # Check for macOS
            scale_win = 1 #for macos the scaling does not work, so I set to 1
            fontsize = sg.set_options(font=("Helvetica", 14))
            default_size = 14
            return layouts.layout_macos, scale_win, fontsize, default_size
        else:  # Linux
            scale_win = 1.5
            fontsize = sg.set_options(font=("Helvetica", 10))
            default_size = 10
            return layouts.layout_linux, scale_win, fontsize, default_size
    else:
        print ("Operating system not recognised. Using Linux layout") # In case the system is not recognized, using the Linux layout
        scale_win = 1.5
        fontsize = sg.set_options(font=("Helvetica", 10))
        return layouts.layout_linux, scale_win, fontsize, default_size


#Function to check if the spectralTemplates folder is available
SPECTRAL_TEMPLATES_DIR = os.path.join(BASE_DIR, "spectralTemplates")

# Link to my website to download the spectralTemplates folder
DOWNLOAD_URL = "https://www.danielegasparri.com/spectralTemplates.zip"

# Temporary path to save the zipped file
TEMP_ZIP_PATH = os.path.join(BASE_DIR, "spectralTemplates.zip")

def download_with_progress(url, dest):
    """Download a file with a progress bar"""

    context = ssl.create_default_context(cafile=certifi.where())

    # Get the file size
    response = urllib.request.urlopen(url, context=context)
    total_size = int(response.getheader('Content-Length', 0))

    # Setting up the progress bar
    layout = [[sg.Text("Downloading spectral templates...")],
            [sg.ProgressBar(total_size, orientation='h', size=(50, 10), key='PROG_BAR')],
            [sg.Cancel()]]
    window = sg.Window("SPAN Download", layout, finalize=True, keep_on_top=True)

    # Opening file
    with open(dest, 'wb') as f:
        downloaded = 0
        block_size = 8192  # 8 KB for block

        while True:
            buffer = response.read(block_size)
            if not buffer:
                break  # end download

            f.write(buffer)
            downloaded += len(buffer)

            # Update the bar
            window['PROG_BAR'].update(downloaded)

            # If 'Cancel' is pressed, I stop the download
            event, _ = window.read(timeout=10)
            if event == "Cancel":
                window.close()
                os.remove(dest)  # and delete the incomplete file downloaded
                sg.popup("Download cancelled!", title="SPAN Error", keep_on_top=True)
                return False

    window.close()
    return True


def download_and_extract_files():
    """Download and extract the file"""
    try:
        # starting the download
        success = download_with_progress(DOWNLOAD_URL, TEMP_ZIP_PATH)
        if not success:
            return  # stop the download if cancelled

        # Extract the downloaded zip file
        with zipfile.ZipFile(TEMP_ZIP_PATH, "r") as zip_ref:
            zip_ref.extractall(BASE_DIR)

        # delete the zip file
        os.remove(TEMP_ZIP_PATH)

        # yeeee success!
        sg.popup("Download completed! SPAN is now ready to use.", title="SPAN Success", keep_on_top=True)

    except Exception as e:
        sg.popup_error(f"Error downloading auxiliary files:\n{str(e)}", title="SPAN Error", keep_on_top=True)


# function to check if the folder 'spectralTemplates' exists in the root folder of SPAN
def check_and_download_spectral_templates():
    sg.theme('LightBlue')
    layout, scale_win, fontsize, default_size = get_layout()
    """Checking if the spectralTemplates/ exists."""
    if not os.path.exists(SPECTRAL_TEMPLATES_DIR):
        # If spectralTemplates does not exist, I should download it, if the user agrees
        choice = sg.popup_yes_no(
            "SPAN must download and extract the spectralTemplates folder to work properly. Do you want to continue? Size = 250MB. This might take a while...\n \nYou can also download the folder here: https://www.danielegasparri.com/spectralTemplates.zip , unzip the folder and put in the root folder of span",
            title="SPAN Missing Files",
            keep_on_top=True
        )

        if choice == "Yes":
            download_and_extract_files()
        else:
            sg.popup(
                "Without the required files, SPAN functionalities are limited, but you can still perform some tasks.",
                title="SPAN Warning",
                keep_on_top=True)
