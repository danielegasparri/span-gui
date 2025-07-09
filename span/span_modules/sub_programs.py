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

# Sub-programs definition and implementation routines

try: #try local import if executed as script
    #GUI import
    from FreeSimpleGUI_local import FreeSimpleGUI as sg
    #SPAN functions import
    from span_functions import system_span as stm
    from span_functions import cube_extract as cubextr
    from span_modules import misc
    from span_modules import layouts
    from params import SpectraParams

except ModuleNotFoundError: #local import if executed as package
    #GUI import
    from span.FreeSimpleGUI_local import FreeSimpleGUI as sg
    #SPAN functions import
    from span.span_functions import system_span as stm
    from span.span_functions import cube_extract as cubextr
    from . import misc
    from . import layouts
    from .params import SpectraParams

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.backends.backend_tkagg
from matplotlib.ticker import MultipleLocator
from matplotlib.backend_bases import MouseButton
from matplotlib.widgets import Slider
from skimage.measure import label, regionprops
import os
import numpy as np
from astropy.io import fits
from astropy.table import Table

from dataclasses import replace

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(CURRENT_DIR)


# 1) PLOT DATA
def plot_data_window(BASE_DIR, layout):

    # Example file to plot
    file_to_plot = os.path.join(BASE_DIR, "example_files", "results", "NGC5320_populations.dat")

    layout, scale_win, fontsize, default_size = misc.get_layout()
    sg.theme('DarkBlue3')
    markers = ['red', 'green', 'yellow', 'blue', 'purple', 'black', 'orange']
    plot_layout = [
        [sg.Text('File to plot:', font=("Helvetica", 15, 'bold'), text_color = 'lightgreen'), sg.InputText(file_to_plot,key='-FILE-', readonly=True, size=(35, 1)), sg.FileBrowse(file_types=(('Text Files', '*.*'),),tooltip='Browse an ASCII file with space or tab soaced columns and first line containing the names of columns'), sg.Button('Load', button_color=('black','light green'), size=(7, 1), font=("Helvetica", 15), tooltip='After you browse for your data file, click here to load it')],
        [sg.Text(' x-axis data:', font=("Helvetica", 15, 'bold')), sg.Listbox(values=[], select_mode=sg.LISTBOX_SELECT_MODE_EXTENDED, key='-X_COLUMN-', enable_events=True, size=(20, 5)), sg.Text('y-axis data:', font=("Helvetica", 15, 'bold')), sg.Listbox(values=[], select_mode=sg.LISTBOX_SELECT_MODE_EXTENDED, key='-Y_COLUMNS-', size=(20, 5))],
        [sg.Checkbox('x errbars:',default = False, key = '-XERRBARS-', font=("Helvetica", 15)), sg.Listbox(values=[], select_mode=sg.LISTBOX_SELECT_MODE_EXTENDED, key='-X_ERR-', enable_events=True, size=(20, 5)), sg.Checkbox('y errbars:',default = False, key = '-YERRBARS-', font=("Helvetica", 15)), sg.Listbox(values=[], select_mode=sg.LISTBOX_SELECT_MODE_EXTENDED, key='-Y_ERR-', size=(20, 5))],
        [sg.Checkbox('Linear Fit of the data', default=False, key='-LINEAR_FIT-'), sg.Push(), sg.Checkbox('X log scale', default=False, key='-X_LOG-', font=("Helvetica", 10)),sg.Checkbox('Y log scale', default=False, key='-Y_LOG-', font=("Helvetica", 10))],
        [sg.HorizontalSeparator()],
        [sg.Text('X-axis Range:'), sg.InputText(key='-X_RANGE_MIN-', size=(8, 1)), sg.Text(' to '), sg.InputText(key='-X_RANGE_MAX-', size=(8, 1)),
        sg.Text('Y-axis Range:'), sg.InputText(key='-Y_RANGE_MIN-', size=(8, 1)), sg.Text(' to '), sg.InputText(key='-Y_RANGE_MAX-', size=(8, 1))],
        [sg.Text('X-axis label:'), sg.InputText(key='-X_LABEL-', size=(20, 1)),sg.Text('Y-axis label:'), sg.InputText(key='-Y_LABEL-', size=(20, 1))],
        [sg.Text('X-axis label size:'), sg.InputText(default_text='14', key='-X_LABEL_SIZE-', size=(5, 1)), sg.Text('Y-axis label size:'), sg.InputText(default_text='14', key='-Y_LABEL_SIZE-', size=(5, 1))],
        [sg.Text('X-ticks size:'), sg.InputText(default_text='14', key='-X_TICKS_SIZE-', size=(5, 1)), sg.Text('Y-ticks size:'), sg.InputText(default_text='14', key='-Y_TICKS_SIZE-', size=(5, 1))],
        [sg.Text('Marker color:'), sg.InputCombo(markers, key='-MARKER_COLOR-',default_value=markers[0],readonly=True), sg.Text('Marker size:'), sg.Slider(range=(1, 100), orientation='h', default_value=40, key='-MARKER_SIZE-')],
        [sg.Text('Plot size (inches):'), sg.InputText('8, 6', key='-PLOT_SIZE-', size=(5, 1)), sg.Checkbox ('Show legend', default = True, key = '-LEGEND-')],
        [sg.Button('Plot!', button_color=('white','orange'), size=(15, 1), font=("Helvetica", 15)), sg.Button('Save img', button_color=('black','light gray'), size=(15, 1), font=("Helvetica", 15)), sg.Button("Help", size=(12, 1),button_color=('black','orange')), sg.Push(), sg.Button('Exit', size=(15, 1))]
    ]

    print ('*** Plotting window open. The main panel will be inactive until you close the window ***')

    plot_window = sg.Window('Data Plotter', plot_layout)

    while True:
        plot_event, plot_values = plot_window.read()

        if plot_event == sg.WIN_CLOSED or plot_event == 'Exit':
            print('Plot window closed. This main panel is now active again')
            print('')
            break

        try:
            file_to_plot = plot_values['-FILE-']
        except Exception:
            sg.popup ('Cannot read the file to plot')
            continue


        if plot_event == 'Load':
            file_path = plot_values['-FILE-']
            if file_path:
                column_names = stm.get_column_names(file_path)
                plot_window['-X_COLUMN-'].update(values=column_names)
                plot_window['-Y_COLUMNS-'].update(values=column_names)
                plot_window['-X_ERR-'].update(values=column_names)
                plot_window['-Y_ERR-'].update(values=column_names)

        elif plot_event == 'Plot!' or plot_event == 'Save img':
            file_path = plot_values['-FILE-']
            x_column = plot_values['-X_COLUMN-']
            y_columns = plot_values['-Y_COLUMNS-']
            x_err = plot_values['-X_ERR-']
            y_err = plot_values['-Y_ERR-']
            x_label = plot_values['-X_LABEL-']
            y_label = plot_values['-Y_LABEL-']
            legend = plot_values['-LEGEND-']
            marker_color = plot_values['-MARKER_COLOR-']
            add_error_bars_x = plot_values['-XERRBARS-']
            add_error_bars_y = plot_values['-YERRBARS-']
            marker_size = int(plot_values['-MARKER_SIZE-'])
            plot_size = tuple(map(float, plot_values['-PLOT_SIZE-'].split(',')))
            x_label_size = int(plot_values['-X_LABEL_SIZE-'])
            y_label_size = int(plot_values['-Y_LABEL_SIZE-'])
            x_tick_size = int(plot_values['-X_TICKS_SIZE-'])
            y_tick_size = int(plot_values['-Y_TICKS_SIZE-'])
            enable_linear_fit = plot_values['-LINEAR_FIT-']
            x_log_scale = plot_values['-X_LOG-']
            y_log_scale = plot_values['-Y_LOG-']

            try:
                x_range_min = float(plot_values['-X_RANGE_MIN-']) if plot_values['-X_RANGE_MIN-'] else None
                x_range_max = float(plot_values['-X_RANGE_MAX-']) if plot_values['-X_RANGE_MAX-'] else None
                y_range_min = float(plot_values['-Y_RANGE_MIN-']) if plot_values['-Y_RANGE_MIN-'] else None
                y_range_max = float(plot_values['-Y_RANGE_MAX-']) if plot_values['-Y_RANGE_MAX-'] else None
            except ValueError:
                sg.popup ('Range values not valid!')
                continue

            if plot_event == 'Plot!':
                stm.plot_data(file_path, x_column, y_columns, x_label, y_label, marker_color, marker_size, plot_size, x_label_size, y_label_size, x_tick_size, y_tick_size, legend, add_error_bars_x, add_error_bars_y, x_err, y_err, False, enable_linear_fit, x_log_scale, y_log_scale, x_range_min, x_range_max, y_range_min, y_range_max)

            if plot_event == 'Save img':
                stm.plot_data(file_path, x_column, y_columns, x_label, y_label, marker_color, marker_size, plot_size, x_label_size, y_label_size, x_tick_size, y_tick_size, legend, add_error_bars_x, add_error_bars_y, x_err, y_err, True, enable_linear_fit, x_log_scale, y_log_scale, x_range_min, x_range_max, y_range_min, y_range_max)

        if plot_event == 'Help':
            f = open(os.path.join(BASE_DIR, "help_files", "help_me_plot.txt"), 'r')
            file_contents = f.read()
            if layout == layouts.layout_android:
                sg.popup_scrolled(file_contents, size=(120, 30))
            else:
                sg.popup_scrolled(file_contents, size=(100, 40))

    plot_window.close()



# 2) PLOT MAPS
def plot_maps_window(BASE_DIR, layout):
    layout, scale_win, fontsize, default_size = misc.get_layout()
    sg.theme("DarkBlue3")
    map_layout = [
        [sg.Text("1. Select the FITS file (*_table.fits) with spaxel and bin info", font=("Helvetica", 14))],
        [sg.Input(key="-FITS-", size=(46, 1), font=("Helvetica", 12)), sg.FileBrowse(file_types=(("FITS files", "*.fits"),), font=("Helvetica", 12))],
        [sg.Text("2. Select the text file with spectral analysis results", font=("Helvetica", 14))],
        [sg.Input(key="-TXT-", size=(46, 1), font=("Helvetica", 12)), sg.FileBrowse(file_types=(("Text files", "*.txt *.dat"),), font=("Helvetica", 12))],
        [sg.Button("Load Files", font=("Helvetica", 14), button_color=('black','light green')), sg.Push(), sg.Button('Help', size=(9, 1), font=("Helvetica", 14), button_color=('black','orange'))],
        [sg.HorizontalSeparator()],
        [sg.Text("Select the quantity to plot:", font=("Helvetica", 14))],
        [sg.Listbox(values=[], size=(44, 10), key="-LIST-", enable_events=True, font=("Helvetica", 14))],
        [sg.Text("Colormap:", font=("Helvetica", 14)), sg.Combo(
            values=["inferno", "viridis", "plasma", "magma", "cividis"],
            default_value="inferno",
            key="-CMAP-",
            readonly=True,
            font=("Helvetica", 12)
        )],
        [sg.Button("Plot Map", size=(9, 1), font=("Helvetica", 14), button_color=('white','orange')), sg.Button("Save selected", size=(13, 1), font=("Helvetica", 14), button_color=('black','light gray')), sg.Button("Save ALL", size=(9, 1), font=("Helvetica", 14), button_color=('black','gray')), sg.Button("Exit", size=(9, 1), font=("Helvetica", 14))]

    ]

    map_window = sg.Window("2D Map Viewer", map_layout)

    x, y, bin_id = None, None, None
    result_df = None
    selected_quantity = None

    while True:
        map_event, map_values = map_window.read()
        if map_event in (sg.WIN_CLOSED, "Exit"):
            break
        elif map_event == "Load Files":
            try:
                fits_path = map_values["-FITS-"]
                txt_path = map_values["-TXT-"]

                # Extract RUN_NAME from the file
                fits_run = fits_path.split("/")[-1].split("_table")[0]
                txt_run = txt_path.split("/")[-1].split("_bins_list")[0]

                if fits_run != txt_run:
                    proceed = sg.popup_yes_no(
                        "WARNING: The selected files appear to belong to different runs.\n\n"
                        f"FITS table run: {fits_run}\nTXT run: {txt_run}\n\n"
                        "Proceed anyway?",
                        title="Run name mismatch"
                    )
                    if proceed != "Yes":
                        sg.popup("File loading aborted.")
                        continue

                x, y, bin_id = stm.load_fits_data(fits_path)
                result_df = stm.load_analysis_results(txt_path)
                col_names = list(result_df.columns)[1:]  # exclude first column (filename)
                map_window["-LIST-"].update(values=col_names)
                sg.popup("Files loaded successfully.")

            except Exception as e:
                sg.popup_error(f"Error loading files:\n{e}")

        elif map_event == "-LIST-":
            if map_values["-LIST-"]:
                selected_quantity = map_values["-LIST-"][0]
        elif map_event == "Plot Map":
            if x is not None and result_df is not None and selected_quantity:
                stm.plot_voronoi_map_clickable(x, y, bin_id, result_df, selected_quantity, cmap=map_values["-CMAP-"])
            else:
                sg.popup("Please load files and select a quantity.")

        elif map_event == "Save selected":
            if x is not None and result_df is not None and selected_quantity:
                save_path = sg.popup_get_file(
                    "Save PNG file",
                    save_as=True,
                    no_window=True,
                    file_types=(("PNG Files", "*.png"),),
                    default_extension=".png"
                )
                if save_path:
                    stm.plot_voronoi_map(x, y, bin_id, result_df, selected_quantity, save_path=save_path, cmap=map_values["-CMAP-"])
                    sg.popup("Image saved successfully.")
            else:
                sg.popup("Please load files and select a quantity before saving.")

        elif map_event == "Save ALL":
            if x is not None and result_df is not None:
                folder = sg.popup_get_folder("Select output folder for PNGs", no_window=True)
                if folder:
                    for quantity in result_df.columns[1:]:  # skip first col (filename)
                        filename = f"{folder}/{stm.sanitize_filename(quantity)}.png"
                        stm.plot_voronoi_map(x, y, bin_id, result_df, quantity, save_path=filename, cmap=map_values["-CMAP-"])
                    sg.popup("All maps saved successfully.")
            else:
                sg.popup("Please load files before saving.")

        elif map_event == 'Help':
            f = open(os.path.join(BASE_DIR, "help_files", "help_maps..txt"), 'r')
            file_contents = f.read()
            if layout == layouts.layout_android:
                sg.popup_scrolled(file_contents, size=(120, 30))
            else:
                sg.popup_scrolled(file_contents, size=(100, 40))

    map_window.close()



#3) TEXT EDITOR
def text_editor_window(layout):
    print('***** Text editor open. The main panel will be inactive until you close the editor *****')

    layout, scale_win, fontsize, default_size = misc.get_layout()
    sg.theme('DarkBlue3')

    if layout == layouts.layout_android:
        editor_layout = [
            [sg.Multiline(size=(120, 20), key='-TEXT-', font=('Helvetica', 12))],
            [sg.Button('Open', key='-OPEN-', size = (15,1), font=('Helvetica', 12), button_color=('black','light green')), sg.Button('Save', key='-SAVE-', size = (15,1), font=('Helvetica', 12)),  sg.Button('Find', size = (15,1), font=('Helvetica', 12)), sg.Button('Undo', size = (15,1), font=('Helvetica', 12)) ],
            [sg.Button('Find/Replace', size = (15,1), font=('Helvetica', 12)), sg.Button('Match rows', size = (15,1), font=('Helvetica', 12)), sg.Button('Create New Column', size = (15,1), font=('Helvetica', 12)), sg.Button('Delete Columns', size = (15,1), font=('Helvetica', 12)), sg.Push(), sg.Button('Close', button_color=('white','orange'), size = (15,1), font=('Helvetica', 12, 'bold'))]
        ]
    else:
        editor_layout = [
            [sg.Multiline(size=(90, 30), key='-TEXT-', font=('Helvetica', 12))],
            [sg.Button('Open', key='-OPEN-', size = (15,1), font=('Helvetica', 12), button_color=('black','light green')), sg.Button('Save', key='-SAVE-', size = (15,1), font=('Helvetica', 12)),  sg.Button('Find', size = (15,1), font=('Helvetica', 12)), sg.Button('Undo', size = (15,1), font=('Helvetica', 12)) ],
            [sg.Button('Find/Replace', size = (15,1), font=('Helvetica', 12)), sg.Button('Match rows', size = (15,1), font=('Helvetica', 12)), sg.Button('Create New Column', size = (15,1), font=('Helvetica', 12)), sg.Button('Delete Columns', size = (15,1), font=('Helvetica', 12)), sg.Push(), sg.Button('Close', button_color=('white','orange'), size = (15,1), font=('Helvetica', 12, 'bold'))]
        ]

    window_editor = sg.Window('Text editor', editor_layout)
    file_modified = False
    text_backup = ""

    while True:
        editor_event, editor_values = window_editor.read()

        if editor_event == sg.WIN_CLOSED or editor_event == 'Close':
            if file_modified:
                confirm_close = sg.popup_yes_no(
                    'Changes have not been saved. Are you sure you want to close?', 'Close')
                if confirm_close == 'Yes':
                    print('Text editor closed. This main panel is now active again')
                    print('')
                    break
            else:
                print('Text editor closed. This main panel is now active again')
                print('')
                break
        elif editor_event == '-SAVE-':
            text = editor_values['-TEXT-']
            filename = sg.popup_get_file('Save the file', save_as=True, default_extension=".txt")
            if filename:
                stm.save_file(filename, text)
                sg.popup(f'The file {filename} has been saved.')
                file_modified = False  # Reset the modification flag
        elif editor_event == '-OPEN-':
            if file_modified:
                confirm_open = sg.popup_yes_no(
                    'Changes have not been saved. Are you sure you want to open a new file?', 'Open')
                if confirm_open == 'No':
                    continue
            filename = sg.popup_get_file('Open file', default_extension=".txt")
            if filename:
                with open(filename, 'r') as file:
                    try:
                        text = file.read()
                        text_backup = text  # Backup the original text
                    except UnicodeDecodeError:
                        sg.popup('Invalid file. I just read simple txt files!')
                        continue
                window_editor['-TEXT-'].update(text)
                file_modified = False  # Reset the modification flag
        elif editor_event == 'Find':
            find_text = sg.popup_get_text('Enter text to find:')
            if find_text:
                text_to_search = editor_values['-TEXT-']
                if find_text in text_to_search:
                    sg.popup(f'Text found at position: {text_to_search.find(find_text)}')
                else:
                    sg.popup('Text not found.')
        elif editor_event == 'Find/Replace':
            find_text = sg.popup_get_text('Enter text to find:')
            if find_text:
                replace_text = sg.popup_get_text('Enter text to replace with:')
                if replace_text is not None:
                    replace_all = sg.popup_yes_no('Replace all occurrences?', 'Replace All')
                    if replace_all == 'Yes':
                        text_to_search = editor_values['-TEXT-']
                        updated_text = stm.find_replace(text_to_search, find_text, replace_text, replace_all)
                        window_editor['-TEXT-'].update(updated_text)
                        file_modified = True  # Set the modification flag
        elif editor_event == 'Match rows':
            sg.theme('DarkBlue3')

            match_layout = [
                [sg.Text('Select the first file:'), sg.InputText(key='-FILE1-', readonly=True),
                    sg.FileBrowse(file_types=(("Text Files", "*.*"),))],
                [sg.Text('Select the second file:'), sg.InputText(key='-FILE2-', readonly=True),
                    sg.FileBrowse(file_types=(("Text Files", "*.*"),))],
                [sg.Text('Select the common column:'), sg.InputText(key='-COMMON_COLUMN-')],
                [sg.Button('Merge'), sg.Button('Exit')]
            ]

            match_window = sg.Window('Match and merge rows', match_layout)

            while True:
                match_event, match_values = match_window.read()

                if match_event == sg.WIN_CLOSED or match_event == 'Exit':
                    break
                elif match_event == 'Merge':
                    file1_path = match_values['-FILE1-']
                    file2_path = match_values['-FILE2-']
                    common_column = match_values['-COMMON_COLUMN-']

                    if file1_path and file2_path and common_column:
                        stm.merge_files(file1_path, file2_path, common_column)
                    else:
                        sg.popup_error('Please select both files and enter a common column.')

            match_window.close()
        elif editor_event == 'Create New Column':
            new_column_name = sg.popup_get_text('Enter the name for the new column:')
            if new_column_name:
                col1_name = sg.popup_get_text('Enter the name of the first column:')
                col2_name = sg.popup_get_text('Enter the name of the second column:')
                expression = sg.popup_get_text('Enter the algebraic expression (e.g., col1 + col2):')
                if col1_name and col2_name and expression:
                    try:
                        df = pd.read_csv(io.StringIO(editor_values['-TEXT-']), delimiter=r'\s+')

                        col1_name_clean = col1_name.replace(' ', '_')
                        col2_name_clean = col2_name.replace(' ', '_')

                        df = stm.create_new_column(df, new_column_name, col1_name_clean, col2_name_clean, expression)

                        if df is not None:
                            window_editor['-TEXT-'].update(df.to_csv(index=False, sep=' ', na_rep=''))
                            file_modified = True  # Set the modification flag
                    except Exception as e:
                        sg.popup_error(f'Error creating the new column: {str(e)}')
                else:
                    sg.popup_error('Please enter names for both columns and the algebraic expression.')

        elif editor_event == 'Delete Columns':
            columns_to_delete = sg.popup_get_text(
                'Enter column names to delete (comma-separated):')
            if columns_to_delete:
                try:
                    df = pd.read_csv(io.StringIO(editor_values['-TEXT-']), delimiter=r'\s+')
                    columns_to_delete_list = [col.strip() for col in columns_to_delete.split(',')]
                    df = df.drop(columns=columns_to_delete_list, errors='ignore')
                    window_editor['-TEXT-'].update(df.to_csv(index=False, sep=' ', na_rep=''))
                    file_modified = True  # Set the modification flag
                except Exception as e:
                    sg.popup_error(f'Error deleting columns: {str(e)}')


        elif editor_event == 'Undo':
            # Ripristina il testo alla sua versione precedente
            window_editor['-TEXT-'].update(text_backup)
            file_modified = True  # Imposta il flag di modifica
        elif editor_event == '-TEXT-':
            # Aggiorna il backup del testo quando viene modificato
            text_backup = editor_values['-TEXT-']
            file_modified = True  # Imposta il flag di modifica


    window_editor.close()



# 4) FITS HEADER EDITOR
def fits_header_window():

    layout, scale_win, fontsize, default_size = misc.get_layout()
    sg.theme('DarkBlue3')
    fitsheader_layout = [
        [sg.Text('Please, select what operation you want to perform on the header of the fits file(s)')],
        [sg.Button('Single fits header editor',button_color= ('black','orange'), size = (22,2), key ='hdr_single_file'), sg.Button('List of fits header editor',button_color= ('black','orange'), size = (22,2), key ='hdr_list_file'), sg.Button('Extract keyword from a list',button_color= ('black','orange'), size = (22,2), key ='extract_keyword')],
        [sg.Button('Close')]
        ]

    print ('*** Fits header editor open. The main panel will be inactive until you close the window ***')
    fitsheader_window = sg.Window('Fits header editor', fitsheader_layout)

    while True:

        fitsheader_event, fitsheader_values = fitsheader_window.read()

        if fitsheader_event == sg.WIN_CLOSED or fitsheader_event == 'Close':
            print ('Fits editor closed. This main panel is now active again')
            print ('')
            break

        #modify/add/delete keyword for a single fits
        if fitsheader_event == 'hdr_single_file':
            sg.theme('DarkBlue3')
            subfitsheader_layout = [
                [sg.Text("Select a FITS file")],
                [sg.Input(key='-FILE-', enable_events=True), sg.FileBrowse()],
                [sg.Multiline(key='-HEADER-', size=(60, 15), disabled=True)],
                [sg.Text("Modify/Add Keyword"),
                sg.Input(key='-KEY-', size=(20, 1)),
                sg.Input(key='-VALUE-', size=(20, 1)),
                sg.Checkbox("Numerical value", key='-NUMERIC-', enable_events=True),
                sg.Button("Add/Modify"), sg.Button("Delete")],  # Added "Delete" button
                [sg.Button("Save Header"), sg.Button("Exit")]
            ]

            subfitsheader_window = sg.Window("Single FITS header editor", subfitsheader_layout)

            while True:
                subfitsheader_event, subfitsheader_values = subfitsheader_window.read()

                if subfitsheader_event == sg.WINDOW_CLOSED or subfitsheader_event == 'Exit':
                    break

                if subfitsheader_event == '-FILE-':
                    file_path = subfitsheader_values['-FILE-']
                    header = stm.read_fits_header(file_path)
                    subfitsheader_window['-HEADER-'].update(repr(header))

                if subfitsheader_event == 'Add/Modify':
                    key = subfitsheader_values['-KEY-']
                    value = subfitsheader_values['-VALUE-']
                    is_numeric = subfitsheader_values['-NUMERIC-']

                    if key:
                        if is_numeric:
                            try:
                                value = float(value)
                            except ValueError:
                                sg.popup_error("Vale must be a number")
                                continue

                        header[key] = value
                        subfitsheader_window['-HEADER-'].update(repr(header))

                if subfitsheader_event == 'Delete':
                    key_to_delete = subfitsheader_values['-KEY-']
                    if key_to_delete:
                        delete_result = stm.delete_keyword(header, key_to_delete)
                        if delete_result is True:
                            subfitsheader_window['-HEADER-'].update(repr(header))
                        else:
                            sg.popup_error(f"Error during deletion: {delete_result}")


                if subfitsheader_event == 'Save Header':
                    if 'file_path' in locals():
                        save_result = stm.save_fits_header(file_path, header)
                        if save_result is True:
                            sg.popup("Header saved with succes!")
                        else:
                            sg.popup_error(f"Error during the saving of the header: {save_result}")
                    else:
                        sg.popup_error("First select a FITS file.")

            subfitsheader_window.close()


        if fitsheader_event == 'hdr_list_file':
            sg.theme('DarkBlue3')

            hdr_list_layout = [
                [sg.Text("Select a list containing only FITS files")],
                [sg.Input(key='-FILELIST-', enable_events=True), sg.FileBrowse()],
                [sg.Text("Select a list containing the keyword you want to change (format: key=value=type)"), sg.Input(key='-KEYFILE-', enable_events=True), sg.FileBrowse()],
                [sg.Multiline(key='-HEADER-', size=(60, 15), disabled=True)],
                [sg.Button("Add/Modify"), sg.Button("Delete"), sg.Button("Exit")]
            ]

            hdr_list_window = sg.Window("FITS header editor", hdr_list_layout)

            while True:
                hdr_list_event, hdr_list_values = hdr_list_window.read()

                if hdr_list_event == sg.WINDOW_CLOSED or hdr_list_event == 'Exit':
                    break

                if hdr_list_event == '-FILELIST-':
                    file_list_path = hdr_list_values['-FILELIST-']

                if hdr_list_event == '-KEYFILE-':
                    key_file_path = hdr_list_values['-KEYFILE-']
                try:
                    if hdr_list_event == 'Add/Modify':
                        if not file_list_path or not key_file_path:
                            sg.popup_error("Something is wrong. Check the files")
                            continue

                        file_paths = stm.read_file_list(file_list_path)
                        if not file_paths:
                            sg.popup_error("Something is wrong. Check the files")
                            continue

                        try:
                            key_paths = stm.read_file_list(key_file_path)
                            key_name, key_value = stm.read_keyword_values_from_file(key_file_path)
                        except ValueError:
                            sg.popup_error('The keyword file is not correct. Chek it!')
                            continue

                        if len(file_paths) != len(key_paths):
                            sg.popup ('The length of the fits file is different from the length of the key file, or you just loaded wrong files. Try again')
                            continue

                        cond = 0

                        try:
                            for i in range (len(file_paths)):
                                with fits.open(file_paths[i], mode='update') as hdul:
                                # Put the keyword in the first HDU)
                                    hdul[0].header[key_name[i]] = (key_value[i])

                                # Save
                                    hdul.flush()

                                cond = cond +1
                                hdr = hdul[0].header

                                hdr_list_window['-HEADER-'].update(repr(hdr))

                        except (AttributeError, ValueError):
                            sg.popup_error('Something is wrong. Check and try again')

                        sg.popup('Successfully modified headers: ', cond, '/', len(file_paths))
                except Exception:
                    sg.popup('File missing')
                    continue

                try:
                    if hdr_list_event == 'Delete':
                        file_paths = stm.read_file_list(file_list_path)
                        key_to_delete = sg.popup_get_text("Enter the keyword to delete:")
                        if key_to_delete:
                            cond = 0
                            try:
                                for i in range(len(file_paths)):
                                    try:
                                        with fits.open(file_paths[i], mode='update') as hdul:
                                            # Aggiungi la keyword all'header del primo HDU (Header Data Unit)
                                            header = hdul[0].header
                                            header.remove(key_to_delete)
                                            # Salva le modifiche
                                            hdul.flush()

                                        cond = cond + 1
                                        hdr = hdul[0].header
                                    except KeyError:
                                        print ('Keyword not found')
                                        continue

                                    hdr_list_window['-HEADER-'].update(repr(hdr))
                            except FileNotFoundError:
                                sg.popup ('Incorrect file or missing')
                                continue
                            if cond > 0:
                                sg.popup(f'Successfully deleted keyword "{key_to_delete}" from headers: {cond}/{len(file_paths)}')
                            else:
                                sg.popup ('Keyword not found')
                except NameError:
                    sg.popup ('No file to process!')
                    continue
            hdr_list_window.close()


        if fitsheader_event == 'extract_keyword':

            sg.theme('DarkBlue3')

            ext_key_layout = [
                [sg.Text("Select a list of FITS files")],
                [sg.Input(key='-FILELIST-', enable_events=True), sg.FileBrowse()],
                [sg.Text("Insert the keyword to extract (case insensitive)"), sg.Input(key='-KEYWORD-')],
                [sg.Multiline(key='-OUTPUT-', size=(60, 15), disabled=True)],
                [sg.Button("Extract and Save"), sg.Button("Exit")]
            ]

            ext_key_window = sg.Window("Extract and Save Keyword", ext_key_layout)

            while True:
                ext_key_event, ext_key_values = ext_key_window.read()

                if ext_key_event == sg.WINDOW_CLOSED or ext_key_event == 'Exit':
                    break

                if ext_key_event == '-FILELIST-':
                    file_list_path = ext_key_values['-FILELIST-']

                try:
                    if ext_key_event == 'Extract and Save':
                        if not file_list_path:
                            sg.popup_error("Select a list of FITS files with relative path included.")
                            continue

                        keyword = ext_key_values['-KEYWORD-'].strip()
                        if not keyword:
                            sg.popup_error("Insert the keyword you want to extract")
                            continue
                        try:
                            file_paths = [line.strip() for line in open(file_list_path) if not line.startswith('#')]
                            data = []
                        except Exception:
                            sg.popup ('Problems with the fits list file. Chek it, please')
                            continue

                        for file_path in file_paths:
                            file_path = file_path.strip()
                            value = stm.extract_keyword(file_path, keyword)
                            data.append({'file': file_path, 'keyword': keyword, 'value': value})

                        ext_key_window['-OUTPUT-'].update('')
                        for entry in data:
                            ext_key_window['-OUTPUT-'].print(f"{entry['file']} - {entry['keyword']}: {entry['value']}")

                        output_file = sg.popup_get_file('Save on file', save_as=True, file_types=(("Text Files", "*.txt"),))
                        if output_file:
                            stm.save_to_text_file(data, output_file)
                            sg.popup(f"Results saved on: '{output_file}'")
                except NameError:
                    sg.popup ('File not found!')
                    continue


            ext_key_window.close()

    fitsheader_window.close()



# 5) LONG-SLIT EXTRACTION
def long_slit_extraction(BASE_DIR, layout, params):

    file_path_spec_extr = params.file_path_spec_extr
    trace_y_range_str = params.trace_y_range_str
    poly_degree_str = params.poly_degree_str
    extract_y_range_str = params.extract_y_range_str
    snr_threshold_str = params.snr_threshold_str
    pixel_scale_str = params.pixel_scale_str
    result_data = params.result_data
    result_long_slit_extract = params.result_long_slit_extract

    layout, scale_win, fontsize, default_size = misc.get_layout()
    sg.theme('DarkBlue3')
    x_axis = np.array([])
    # Define FreeSimpleGUI layout
    spec_extr_layout = [
        [sg.Text("Select FITS File:", font=("Helvetica", 15, 'bold'), text_color = 'light blue'), sg.InputText(default_text = file_path_spec_extr, key="file_path", size = (38,1)), sg.FileBrowse()],
        [sg.Text("Select Y Range for Trace Fitting:"), sg.InputText(default_text= trace_y_range_str, key="trace_y_range",size=(12, 1)), sg.Text("Degree to fit:"), sg.InputText(default_text= poly_degree_str, key="poly_degree",size=(5, 1))],
        [sg.Button("1) Open 2D Spectrum", button_color=('black','light blue'), size=(18, 1), tooltip='First we have to load and visualise the spectrum'), sg.Button("2) Fit Trace",button_color=('black','light green'), size=(19, 1), tooltip='Now we find the trace of the spectrum along the dispersion axis'), sg.Button("3) Correct Spectrum",button_color=('black','orange'), size=(18, 1), tooltip='Finally we correct the distortion of the spectrum before the extraction')],
        [sg.HorizontalSeparator()],
        [sg.Button("Extract 1D Spectrum",size=(20, 1), tooltip='Extract 1D spectrum within the selected Y range. Useful for point sources'), sg.Text("Y Range for Extract 1D Spectrum:"), sg.InputText(default_text=extract_y_range_str, key="extract_y_range",size=(12, 1))],
        [sg.Button("Extract SNR bin Spectra",size=(20, 1), tooltip='Extract n bins with the selected SNR Threshold. Useful for extended sources'), sg.Text('SNR Threshold:'), sg.InputText(key='snr',size=(4, 1), default_text=snr_threshold_str), sg.Text('Pix scale ("/pix):'), sg.InputText(key='pix_scale',size=(5, 1), default_text=pixel_scale_str)],
        [sg.Canvas(key="-CANVAS-")],
        [sg.Button("Help", size=(12, 1),button_color=('black','orange')), sg.Push(), sg.Button("Exit", size=(12, 1))]
    ]

    print ('*** 2D spectra extraction open. The main panel will be inactive until you close the window ***')

    spec_extr_window = sg.Window("2D spectra extraction", spec_extr_layout, finalize=True)
    canvas_elem = spec_extr_window["-CANVAS-"]
    canvas = canvas_elem.Widget

    trace_model = None

    # Event loop
    while True:
        spec_extr_event, spec_extr_values = spec_extr_window.read()

        if spec_extr_event == (sg.WIN_CLOSED):
            print('2D spec window closed. This main panel is now active again')
            print('')
            break

        file_path_spec_extr= spec_extr_values['file_path']
        trace_y_range_str= spec_extr_values['trace_y_range']
        poly_degree_str= spec_extr_values['poly_degree']
        extract_y_range_str= spec_extr_values['extract_y_range']
        snr_threshold_str= spec_extr_values['snr']
        pixel_scale_str= spec_extr_values['pix_scale']


        if spec_extr_event == ('Exit'):
            print('2D spec window closed. This main panel is now active again')
            print('')
            break

        if spec_extr_event == "1) Open 2D Spectrum":

            try:
                trace_model = None
                spectrum, header = stm.open_fits(file_path_spec_extr)
                x_axis = np.arange(len(spectrum[0]))
                plt.imshow(spectrum, cmap="viridis", norm=LogNorm())
                plt.title("2D Spectrum")
                plt.show()
                plt.close()
            except Exception:
                sg.popup ('Spectrum not valid. Must be a 2D fits image!')

        if spec_extr_event == "2) Fit Trace":

            try:
                trace_y_range = eval(trace_y_range_str)
                poly_degree = int(poly_degree_str)

                if poly_degree < 1 or poly_degree > 5:
                    sg.popup ('The polynomial degree should be between 1 and 5')
                    continue

                trace_model = stm.find_and_fit_spectroscopic_trace(spectrum, trace_y_range, poly_degree, True, True)

            except Exception as e:
                sg.popup_error(f"Error: {str(e)}")

        if spec_extr_event == "3) Correct Spectrum":
            trace_y_range = eval(trace_y_range_str)

            if trace_model is not None:
                corrected_spectrum = stm.correct_distortion_slope(spectrum, trace_model, trace_y_range)
            else:
                sg.popup_error("Please find and fit the spectroscopic trace first.")

        if spec_extr_event == "Extract 1D Spectrum":
            if trace_model is not None:
                try:
                    extract_y_range = eval(extract_y_range_str)
                    extracted_filename = os.path.splitext(os.path.basename(file_path_spec_extr))[0]
                    #creating the directory
                    result_longslit_extraction = result_data+'/longslit_extracted/'+extracted_filename+'/'
                    os.makedirs(result_longslit_extraction, exist_ok=True)
                    stm.extract_1d_spectrum(corrected_spectrum, extract_y_range, header, x_axis, output_fits_path= (result_longslit_extraction + f"{extracted_filename}_extracted_.fits"))
                    sg.popup ('1D spectrum saved in the working directory')
                except Exception:
                    sg.popup('Extraction parameters not valid. Spectrum not extracted')
                    continue
            else:
                sg.popup_error("Please find and fit the spectroscopic trace first.")

        if spec_extr_event == "Extract SNR bin Spectra":
            if trace_model is not None and 'corrected_spectrum' in locals():

                try:
                    snr_threshold = float(snr_threshold_str)
                    pixel_scale = float(pixel_scale_str)
                except Exception:
                    sg.popup ('SNR or pixel scale values not valid!')
                    continue

                if snr_threshold < 0 or pixel_scale < 0:
                    sg.popup ('SNR and pixel scale values canoot be negative')
                    continue

                y_correction_trace_position = trace_y_range[0]
                stm.extract_and_save_snr_spectra(corrected_spectrum, trace_model, header, x_axis, snr_threshold, pixel_scale, file_path_spec_extr, y_correction_trace_position, result_long_slit_extract)

                #create spectra list of the bins to use with SPAN
                extracted_filename = os.path.splitext(os.path.basename(file_path_spec_extr))[0]
                result_longslit_extraction_bins = result_data+'/longslit_extracted/'+extracted_filename+'/bins/'
                os.makedirs(result_longslit_extraction_bins, exist_ok=True)
                file_list = stm.get_files_in_folder(result_longslit_extraction_bins)
                output_file = extracted_filename + '_bins_list.txt'
                stm.save_to_text_file(file_list, output_file)
                sg.Popup('Spectra file list of the bins saved in the working directory', output_file, 'You can now browse and load this list file')

            else:
                sg.popup_error("Please correct the spectrum and find the spectroscopic trace first.")

        if spec_extr_event == 'Help':
            f = open(os.path.join(BASE_DIR, "help_files", "help_2d_spec.txt"), 'r')
            file_contents = f.read()
            if layout == layouts.layout_android:
                sg.popup_scrolled(file_contents, size=(120, 30))
            else:
                sg.popup_scrolled(file_contents, size=(100, 40))

    spec_extr_window.close()


    #updating the params
    params = replace(params,
                    file_path_spec_extr = file_path_spec_extr,
                    trace_y_range_str = trace_y_range_str,
                    poly_degree_str = poly_degree_str,
                    extract_y_range_str = extract_y_range_str,
                    snr_threshold_str = snr_threshold_str,
                    pixel_scale_str = pixel_scale_str,
                     )


    return params



# 6) DATACUBE EXTRACTION
def datacube_extraction(params):

    result_data = params.result_data
    ifs_run_id = params.ifs_run_id
    ifs_input = params.ifs_input
    ifs_redshift = params.ifs_redshift
    ifs_lfs_data_default = params.ifs_lfs_data_default
    ifs_ow_config = params.ifs_ow_config
    ifs_ow_output = params.ifs_ow_output
    ifs_lmin_tot = params.ifs_lmin_tot
    ifs_lmax_tot = params.ifs_lmax_tot
    ifs_preloaded_routine = params.ifs_preloaded_routine
    ifs_min_snr_mask = params.ifs_min_snr_mask
    ifs_target_snr = params.ifs_target_snr
    ifs_routine_read = params.ifs_routine_read
    ifs_routine_read_default = params.ifs_routine_read_default
    ifs_user_routine = params.ifs_user_routine
    ifs_user_routine_file = params.ifs_user_routine_file
    ifs_origin = params.ifs_origin
    ifs_mask = params.ifs_mask
    ifs_output = params.ifs_output
    ifs_lmin_snr_default = params.ifs_lmin_tot #params.ifs_lmin_snr_default
    ifs_lmax_snr_default = params.ifs_lmax_tot #params.ifs_lmax_snr_default
    ifs_manual_bin = params.ifs_manual_bin
    ifs_voronoi = params.ifs_voronoi
    ifs_existing_bin = params.ifs_existing_bin
    ifs_existing_bin_folder = params.ifs_existing_bin_folder
    ifs_bin_method = params.ifs_bin_method
    ifs_covariance = params.ifs_covariance
    ifs_prepare_method = params.ifs_prepare_method



    layout, scale_win, fontsize, default_size = misc.get_layout()
    sg.theme('LightBlue1')

    cube_ifs_layout = [
        [sg.Text('Select a fits cube:', font = ('', default_size, 'bold'), tooltip='Select a datacube WITHIN the inputData folder'), sg.InputText(ifs_input, size=(30, 1), key = 'ifs_input'), sg.FileBrowse(file_types=(('fits file', '*.fits *.fit'),)), sg.Button('View datacube', button_color=('black','light blue'), size = (18,1), tooltip='Take a look at the datacube, it may be useful')],
        [sg.Text('Name of the run:', tooltip='Just give a name for this session'), sg.InputText(ifs_run_id, size = (15,1), key = 'ifs_run_id'), sg.Text('z:', tooltip='Redshift estimation. Put zero to not correct for redshift'), sg.InputText(ifs_redshift, size = (8,1), key = 'ifs_redshift'), sg.Text('Wave to extract (A):', tooltip='Wavelength range you want to extract. Look at the datacube if you do not know'), sg.InputText(ifs_lmin_tot, size = (6,1), key = 'ifs_lmin_tot'), sg.Text('-'), sg.InputText(ifs_lmax_tot, size = (6,1), key = 'ifs_lmax_tot')],

        [sg.HorizontalSeparator()],

        [sg.Radio('Using a pre-loaded routine for extraction:', "RADIOCUBEROUTINE", default = ifs_preloaded_routine, key = 'ifs_preloaded_routine', font = ('', default_size, 'bold'), tooltip='These are pre-loaded routines for reading the most commin datacubes'), sg.InputCombo(ifs_routine_read,key='ifs_routine_read',default_value=ifs_routine_read_default, readonly=True, size = (18,1))],
        [sg.Radio('Using a user defined routine for extraction:', "RADIOCUBEROUTINE", default = ifs_user_routine, key = 'ifs_user_routine', font = ('', default_size, 'bold'), tooltip='If you have your .py datacube read routine, load it here'), sg.InputText(ifs_user_routine_file, size=(15, 1), key = 'ifs_user_routine_file'), sg.FileBrowse(file_types=(('py file', '*.py'),))],
        [sg.Text('Origin (in pixel) of the coordinate systems:', tooltip='Pixel coordinates of the centre or the object you want to study. Look at the datacube to know it'), sg.InputText(ifs_origin, size = (9,1), key = 'ifs_origin')],

        [sg.HorizontalSeparator()],

        [sg.Text('Select a mask:', font = ('', default_size, 'bold'), tooltip='If you do not have a mask file, you can create with the button on the right'), sg.InputText(ifs_mask, size=(20, 1), key = 'ifs_mask'), sg.FileBrowse(file_types=(('fits file', '*.fits'),)), sg.Button('Generate mask',button_color=('black','light blue'), size = (18,1), tooltip='Generate a mask, even without really masking any spaxel'), sg.Text('Min S/N to mask:', tooltip='Masking all the spaxels with low S/N'), sg.InputText(ifs_min_snr_mask, size = (6,1), key = 'ifs_min_snr_mask')],

        [sg.HorizontalSeparator()],

        [sg.Radio('Voronoi bin', "RADIOVOR", default=ifs_voronoi, key='ifs_voronoi', tooltip='Automatic voronoi rebinning'), sg.Text('S/N:', tooltip='Select the S/N treshold of the binned spaxels. A good starting value is 30-50'), sg.InputText(ifs_target_snr, size = (4,1), key = 'ifs_target_snr'), sg.Radio('Manual bin', "RADIOVOR", default=ifs_manual_bin, key='ifs_manual_bin', tooltip='Select region(s) to bin. Masking is not applied here'), sg.Button('Manual binning'), sg.Radio('Existing bins', "RADIOVOR", default=ifs_existing_bin, key='ifs_existing_bin', tooltip='Use already available mask and bin info'), sg.Input(ifs_existing_bin_folder, key='ifs_existing_bin_folder', size = (11,1)), sg.FolderBrowse(tooltip='Browse the folder where your *_table.fits and *_mask.fits are located')],

        [sg.HorizontalSeparator()],

        [sg.Button('Preview bins',button_color=('black','light green'), size = (18,1)), sg.Button('Extract!',button_color= ('white','black'), size = (18,1)), sg.Push(), sg.Button('I need help',button_color=('black','orange'), size = (12,1)), sg.Exit(size=(18, 1))],
    ]

    print ('*** Cube extraction routine open. The main panel will be inactive until you close the window ***')
    cube_ifs_window = sg.Window('Cube extraction using GIST standard', cube_ifs_layout)

    while True:

        cube_ifs_event, cube_ifs_values = cube_ifs_window.read()

        if cube_ifs_event == (sg.WIN_CLOSED):
            print ('Cube extraction routine closed. This main panel is now active again')
            print ('')
            break

        #assigning user values
        ifs_run_id = cube_ifs_values['ifs_run_id']
        ifs_input = cube_ifs_values['ifs_input']
        ifs_routine_read_default = cube_ifs_values['ifs_routine_read']
        ifs_routine_selected  = os.path.join(BASE_DIR, "span_functions", "cube_extract_functions", f"{ifs_routine_read_default}.py")

        ifs_origin = cube_ifs_values['ifs_origin']
        ifs_mask = cube_ifs_values['ifs_mask']
        ifs_output_dir = ifs_output + ifs_run_id


        ifs_preloaded_routine = cube_ifs_values['ifs_preloaded_routine']
        ifs_user_routine = cube_ifs_values['ifs_user_routine']
        ifs_user_routine_file = cube_ifs_values['ifs_user_routine_file']

        ifs_manual_bin = cube_ifs_values['ifs_manual_bin']
        ifs_voronoi = cube_ifs_values['ifs_voronoi']

        ifs_existing_bin = cube_ifs_values['ifs_existing_bin']
        if ifs_existing_bin:
            ifs_existing_bin_folder = cube_ifs_values['ifs_existing_bin_folder']

        try:
            ifs_redshift = float(cube_ifs_values['ifs_redshift'])
            ifs_lmin_tot = float(cube_ifs_values['ifs_lmin_tot'])
            ifs_lmax_tot = float(cube_ifs_values['ifs_lmax_tot'])
            ifs_min_snr_mask = float(cube_ifs_values['ifs_min_snr_mask'])
            ifs_target_snr = float(cube_ifs_values['ifs_target_snr'])

            #converting to A
            ifs_lmin_snr = ifs_lmin_tot
            ifs_lmax_snr = ifs_lmax_tot

        except Exception:
            sg.popup ('Invalid input parameters!')
            continue


        if ifs_user_routine:
            ifs_routine_selected = ifs_user_routine_file


        if cube_ifs_event == ('Exit'):
            print ('Cube extraction routine closed. This main panel is now active again')
            print ('')
            break

        #routine to view the datacube
        if cube_ifs_event == 'View datacube':

            try:
                # Load the datacube
                data, wave = stm.read_datacube(ifs_input)

                # Create the figure and axes with adjusted layout
                fig, ax_img = plt.subplots(figsize=(10, 7))
                plt.subplots_adjust(bottom=0.15,top=0.95)  # Less space for the slider

                # Display the initial image of the datacube
                img = ax_img.imshow(data[0], cmap="gray", origin='lower')
                ax_img.set_title(f"Wavelength Index: {0}", pad=15)

                # Position the slider just below the main plot
                ax_slider = plt.axes([0.15, 0.05, 0.7, 0.03], facecolor='lightgoldenrodyellow')
                slider = Slider(ax_slider, 'Wavelength', wave[0], wave[-1], valinit=wave[0], valfmt='%0.0f')

                # Update function for the slider
                def update(val):
                    wave_cube = slider.val
                    closest_index = (np.abs(wave - wave_cube)).argmin()
                    img.set_data(data[closest_index])
                    ax_img.set_title(f'Wavelength: {wave_cube:.2f} Å')
                    fig.canvas.draw_idle()

                slider.on_changed(update)

                plt.show()
                plt.close()

            except Exception as e:
                sg.popup('Datacube not found or not valid')
                continue


        #routine for generating a mask
        if cube_ifs_event == 'Generate mask':

            # Considering non mobile devices with a keyboard to perform the masking. Using ctrl+click to mask
            if layout != layouts.layout_android:

                try:
                    ifs_input_mask = ifs_input
                    data, wave = stm.read_datacube(ifs_input_mask)

                    # Preparing the mask and show the first slice of the cube
                    mask = np.zeros(data.shape[1:], dtype=int)
                    fig, (ax_img, ax_slider) = plt.subplots(
                        2, 1, gridspec_kw={'height_ratios': [9, 1]}, figsize=(12, 8)
                    )
                    plt.subplots_adjust(hspace=0.3)  # Slider-image spacing

                    # # Showing the image of the first wavelength value
                    img = ax_img.imshow(data[0], cmap="gray", origin='lower')
                    mask_overlay = ax_img.imshow(
                        np.ma.masked_where(mask == 0, mask),
                        cmap='Reds', alpha=0.5, origin='lower'
                    )
                    ax_img.set_title(f"Wavelength Index: {0}")

                    # Slider for wavelength selection
                    slider = plt.Slider(ax_slider, 'Wavelength Index', 0, data.shape[0] - 1, valinit=0, valfmt='%0.0f')

                    # Instructions on the plot
                    instructions = (
                        "Ctrl + left click: mask a pixel\n"
                        "Ctrl + right click: unmask a pixel\n"
                        "Ctrl + left click and drag: mask an area\n"
                        "Ctrl + right click and drag: unmask an area\n"
                        "Close this window to save the mask"
                    )
                    ax_img.text(1.05, 0.5, instructions, transform=ax_img.transAxes,
                                fontsize=10, va='center', ha='left', color='blue')

                    # Variables for tracing the selection
                    # global start_point, dragging, deselecting
                    start_point = None
                    dragging = False
                    deselecting = False

                    # Updating the slider
                    def update(val):
                        wavelength_index = int(slider.val)
                        img.set_data(data[wavelength_index])
                        ax_img.set_title(f"Wavelength Index: {wavelength_index}")
                        fig.canvas.draw_idle()

                    slider.on_changed(update)

                    # Functions for selecting and de-selecting the spaxels to mask
                    def on_click(event):
                        global start_point, dragging, deselecting
                        # Ignoring the mask if using the slider
                        if event.inaxes != ax_img:
                            return

                        # Check if Crtl or Control is pressed
                        if event.key is None or (('control' not in event.key.lower()) and ('ctrl' not in event.key.lower())):
                            return

                        if event.button == MouseButton.LEFT:  # Ctrl + left click to mask
                            start_point = (int(event.xdata), int(event.ydata))
                            dragging = True
                            deselecting = False
                        elif event.button == MouseButton.RIGHT:  # Ctrl + right click to unmask
                            start_point = (int(event.xdata), int(event.ydata))
                            dragging = False
                            deselecting = True

                    # Dragging functions
                    def on_release(event):
                        global start_point, dragging, deselecting
                        if event.inaxes == ax_img and start_point:
                            end_point = (int(event.xdata), int(event.ydata))
                            x0, y0 = start_point
                            x1, y1 = end_point

                            if dragging:  # Masking
                                mask[min(y0, y1):max(y0, y1)+1, min(x0, x1):max(x0, x1)+1] = 1
                            elif deselecting:  # Unmasking
                                mask[min(y0, y1):max(y0, y1)+1, min(x0, x1):max(x0, x1)+1] = 0

                            # Update the figure
                            mask_overlay.set_data(np.ma.masked_where(mask == 0, mask))
                            fig.canvas.draw()
                            start_point = None
                            dragging = False
                            deselecting = False

                    fig.canvas.mpl_connect("button_press_event", on_click)
                    fig.canvas.mpl_connect("button_release_event", on_release)
                    plt.show()
                    plt.close()

                    mask_name = "mask_" + ifs_run_id + "_.fits"
                    stm.save_mask_as_fits(mask, result_data + "/" + mask_name)
                    mask_path = result_data + '/' + mask_name
                    sg.popup("Mask saved as ", mask_name, "in ", result_data, " folder and loaded.")

                    #Updating the mask path in the GUI
                    cube_ifs_window['ifs_mask'].update(mask_path)

                except Exception as e:
                    sg.popup('Fits datacube not valid.')
                    continue



            # Considering mobile devices without a keyboard to perform the masking. Using simple tap and drag to mask
            if layout == layouts.layout_android:
                try:
                    ifs_input_mask = ifs_input
                    data, wave = stm.read_datacube(ifs_input_mask)

                    # Preparing the mask and show the first slice of the cube
                    mask = np.zeros(data.shape[1:], dtype=int)
                    fig, (ax_img, ax_slider) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [9, 1]}, figsize=(12, 8))
                    plt.subplots_adjust(hspace=0.3)  # Slider-image spacing

                    # Showing the image of the first wavelength value
                    img = ax_img.imshow(data[0], cmap="gray", origin='lower')
                    mask_overlay = ax_img.imshow(np.ma.masked_where(mask == 0, mask), cmap='Reds', alpha=0.5, origin='lower')
                    ax_img.set_title(f"Wavelength Index: {0}")

                    # Slider for wavelength selection
                    slider = plt.Slider(ax_slider, 'Wavelength Index', 0, data.shape[0] - 1, valinit=0, valfmt='%0.0f')

                    # Instructions on the plot
                    instructions = (
                        "Left click: mask a pixel\n"
                        "Right click: unmask pixel\n"
                        "Left click and drag: mask an area\n"
                        "Right click and drag: unmask an area\n"
                        "Close this window to save the mask"
                    )
                    ax_img.text(1.05, 0.5, instructions, transform=ax_img.transAxes, fontsize=10, va='center', ha='left', color='blue')

                    # Variables for tracing the selection
                    start_point = None
                    dragging = False
                    deselecting = False

                    # Updating the slider
                    def update(val):
                        wavelength_index = int(slider.val)
                        img.set_data(data[wavelength_index])
                        ax_img.set_title(f"Wavelength Index: {wavelength_index}")
                        fig.canvas.draw_idle()

                    slider.on_changed(update)

                    # Functions for selecting and de-selecting the spaxels to mask
                    def on_click(event):
                        global start_point, dragging, deselecting
                        # Ignoring the mask if using the slider
                        if event.inaxes == ax_img:
                            if event.button == MouseButton.LEFT:  # Left click for masking
                                start_point = (int(event.xdata), int(event.ydata))
                                dragging = True
                                deselecting = False
                            elif event.button == MouseButton.RIGHT:  # RIght click for unmasking
                                start_point = (int(event.xdata), int(event.ydata))
                                dragging = False
                                deselecting = True

                    # Dragging functions
                    def on_release(event):
                        global start_point, dragging, deselecting
                        if event.inaxes == ax_img and start_point:
                            end_point = (int(event.xdata), int(event.ydata))
                            x0, y0 = start_point
                            x1, y1 = end_point

                            if dragging:  # dragging with the left mouse click
                                mask[min(y0, y1):max(y0, y1)+1, min(x0, x1):max(x0, x1)+1] = 1
                            elif deselecting:  # Dragging with the right mouse click
                                mask[min(y0, y1):max(y0, y1)+1, min(x0, x1):max(x0, x1)+1] = 0

                            # Update mask and plot
                            mask_overlay.set_data(np.ma.masked_where(mask == 0, mask))
                            fig.canvas.draw()
                            start_point = None
                            dragging = False
                            deselecting = False

                    fig.canvas.mpl_connect("button_press_event", on_click)
                    fig.canvas.mpl_connect("button_release_event", on_release)
                    plt.show()
                    plt.close()

                    mask_name = "mask_"+ifs_run_id+"_.fits"
                    stm.save_mask_as_fits(mask, result_data+"/"+mask_name)
                    mask_path = result_data + '/'+ mask_name
                    sg.popup("Mask saved as ", mask_name, "in ", result_data, " folder and loaded.")

                    #Updating the mask path in the GUI
                    cube_ifs_window['ifs_mask'].update(mask_path)
                except Exception as e:
                    sg.popup ('Fits datacube not valid.')
                    continue


        if ifs_existing_bin:
            cubextr.handle_existing_bin_files(ifs_existing_bin_folder, ifs_output_dir, ifs_run_id)

        # preview mode for voronoi binning
        if cube_ifs_event == 'Preview bins' and not ifs_manual_bin:
            voronoi = True
            preview = True

            # Creating the disctionary to be passed to the cube_extract module
            config = cubextr.buildConfigFromGUI(
                ifs_run_id, ifs_input, ifs_output_dir, ifs_redshift,
                ifs_lfs_data_default, ifs_ow_config,
                ifs_ow_output, ifs_routine_selected, ifs_origin,
                ifs_lmin_tot, ifs_lmax_tot, ifs_lmin_snr, ifs_lmax_snr,
                ifs_min_snr_mask, ifs_mask, ifs_bin_method, ifs_target_snr,
                ifs_covariance, ifs_prepare_method)
            try:
                cubextr.extract(config, preview, voronoi, ifs_manual_bin, ifs_existing_bin)
            except Exception as e:
                sg.popup("Error showing the bins:", str(e))
                continue



        # Performing manual binning by the user, by selecting one or multiple regions in a matplotlib iterative window
        # Using a modified version of the mask routine above to select the manual binning regions. Then inverting the mask to consider ONLY the selected spaxels.
        if cube_ifs_event == 'Manual binning':

            #activating the Manual bin option, if not activated
            cube_ifs_window['ifs_manual_bin'].update(True)
            ifs_manual_bin = cube_ifs_values['ifs_manual_bin']

            try:


            # Considering non mobile devices with a keyboard to perform the masking. Using ctrl+click to mask
                if layout != layouts.layout_android:

                    ifs_input_mask = ifs_input
                    data, wave = stm.read_datacube(ifs_input_mask)

                    # Preparing the mask and show the first slice of the cube
                    bin_mask = np.zeros(data.shape[1:], dtype=int)
                    fig, (ax_img, ax_slider) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [9, 1]}, figsize=(12, 8))
                    plt.subplots_adjust(hspace=0.3)  # Slider-image spacing

                    # Showing the image of the first wavelength value
                    img = ax_img.imshow(data[0], cmap="gray", origin='lower')
                    bin_mask_overlay = ax_img.imshow(np.ma.masked_where(bin_mask == 0, bin_mask), cmap='Reds', alpha=0.5, origin='lower')
                    ax_img.set_title(f"Wavelength Index: {0}")

                    # Slider for wavelength selection
                    slider = plt.Slider(ax_slider, 'Wavelength Index', 0, data.shape[0] - 1, valinit=0, valfmt='%0.0f')

                    # Instructions on the plot
                    instructions = (
                        "Ctrl + Left click: Select a spaxel\n"
                        "Ctrl + Right click: Deselect a spaxel\n"
                        "Ctrl + Left click and drag: Select an area\n"
                        "Ctrl + Right click and drag: Deselect an area\n"
                        "Ctrl + Close this window to save the work"
                    )
                    ax_img.text(1.05, 0.5, instructions, transform=ax_img.transAxes,
                                fontsize=10, va='center', ha='left', color='blue')


                    # Variables for tracing the selection
                    start_point = None
                    dragging = False
                    deselecting = False

                    # Updating the slider
                    def update(val):
                        wavelength_index = int(slider.val)
                        img.set_data(data[wavelength_index])
                        ax_img.set_title(f"Wavelength Index: {wavelength_index}")
                        fig.canvas.draw_idle()

                    slider.on_changed(update)

                    # Functions for selecting and de-selecting the spaxels to mask
                    def on_click(event):
                        global start_point, dragging, deselecting
                        # Ignoring the mask if using the slider
                        if event.inaxes != ax_img:
                            return

                        # Check if Crtl or Control is pressed
                        if event.key is None or (('control' not in event.key.lower()) and ('ctrl' not in event.key.lower())):
                            return

                        if event.button == MouseButton.LEFT:  # Ctrl + left click to mask
                            start_point = (int(event.xdata), int(event.ydata))
                            dragging = True
                            deselecting = False
                        elif event.button == MouseButton.RIGHT:  # Ctrl + right click to unmask
                            start_point = (int(event.xdata), int(event.ydata))
                            dragging = False
                            deselecting = True

                    # Dragging functions
                    def on_release(event):
                        global start_point, dragging, deselecting
                        if event.inaxes == ax_img and start_point:
                            end_point = (int(event.xdata), int(event.ydata))
                            x0, y0 = start_point
                            x1, y1 = end_point

                            if dragging:  # Masking
                                bin_mask[min(y0, y1):max(y0, y1)+1, min(x0, x1):max(x0, x1)+1] = 1
                            elif deselecting:  # Unmasking
                                bin_mask[min(y0, y1):max(y0, y1)+1, min(x0, x1):max(x0, x1)+1] = 0

                            # Update mask and plot
                            bin_mask_overlay.set_data(np.ma.masked_where(bin_mask == 0, bin_mask))
                            fig.canvas.draw()
                            start_point = None
                            dragging = False
                            deselecting = False

                    fig.canvas.mpl_connect("button_press_event", on_click)
                    fig.canvas.mpl_connect("button_release_event", on_release)
                    plt.show()
                    plt.close()



                if layout == layouts.layout_android:
                    ifs_input_mask = ifs_input
                    data, wave = stm.read_datacube(ifs_input_mask)

                    # Preparing the mask and show the first slice of the cube
                    bin_mask = np.zeros(data.shape[1:], dtype=int)
                    fig, (ax_img, ax_slider) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [9, 1]}, figsize=(12, 8))
                    plt.subplots_adjust(hspace=0.3)  # Slider-image spacing

                    # Showing the image of the first wavelength value
                    img = ax_img.imshow(data[0], cmap="gray", origin='lower')
                    bin_mask_overlay = ax_img.imshow(np.ma.masked_where(bin_mask == 0, bin_mask), cmap='Reds', alpha=0.5, origin='lower')
                    ax_img.set_title(f"Wavelength Index: {0}")

                    # Slider for wavelength selection
                    slider = plt.Slider(ax_slider, 'Wavelength Index', 0, data.shape[0] - 1, valinit=0, valfmt='%0.0f')

                    # Instructions on the plot
                    instructions = (
                        "Left click: Select a spaxel\n"
                        "Right click: Deselect a spaxel\n"
                        "Left click and drag: Select an area\n"
                        "Right click and drag: Deselect an area\n"
                        "Close this window to save the work"
                    )
                    ax_img.text(1.05, 0.5, instructions, transform=ax_img.transAxes, fontsize=10, va='center', ha='left', color='blue')

                    # Variables for tracing the selection
                    start_point = None
                    dragging = False
                    deselecting = False

                    # Updating the slider
                    def update(val):
                        wavelength_index = int(slider.val)
                        img.set_data(data[wavelength_index])
                        ax_img.set_title(f"Wavelength Index: {wavelength_index}")
                        fig.canvas.draw_idle()

                    slider.on_changed(update)

                    # Functions for selecting and de-selecting the spaxels to mask
                    def on_click(event):
                        global start_point, dragging, deselecting
                        # Ignoring the mask if using the slider
                        if event.inaxes == ax_img:
                            if event.button == MouseButton.LEFT:  # Left click for masking
                                start_point = (int(event.xdata), int(event.ydata))
                                dragging = True
                                deselecting = False
                            elif event.button == MouseButton.RIGHT:  # RIght click for unmasking
                                start_point = (int(event.xdata), int(event.ydata))
                                dragging = False
                                deselecting = True

                    def on_release(event):
                        global start_point, dragging, deselecting
                        if event.inaxes == ax_img and start_point:
                            end_point = (int(event.xdata), int(event.ydata))
                            x0, y0 = start_point
                            x1, y1 = end_point

                            if dragging:  # dragging with the left mouse click
                                bin_mask[min(y0, y1):max(y0, y1)+1, min(x0, x1):max(x0, x1)+1] = 1
                            elif deselecting:  # Dragging with the right mouse click
                                bin_mask[min(y0, y1):max(y0, y1)+1, min(x0, x1):max(x0, x1)+1] = 0

                            # Update mask and plot
                            bin_mask_overlay.set_data(np.ma.masked_where(bin_mask == 0, bin_mask))
                            fig.canvas.draw()
                            start_point = None
                            dragging = False
                            deselecting = False

                    fig.canvas.mpl_connect("button_press_event", on_click)
                    fig.canvas.mpl_connect("button_release_event", on_release)
                    plt.show()
                    plt.close()


                #DOING THE MAGIC: Finding all the contigous selected spaxels, assign an integer flag for each region selected by the user.
                labeled_mask = label(bin_mask, connectivity=1)

                #inverting the mask: the selected regions will be the active regions to bin
                bin_mask = 1 - bin_mask

                #saving the masked regions
                bin_mask_name = "bin_mask2D_"+ifs_run_id+"_.fits"
                stm.save_mask_as_fits(bin_mask, result_data+"/"+bin_mask_name)
                bin_mask_path = result_data + '/'+ bin_mask_name

                # Generating the text file with all spaxel info:
                txt_filename = f"mask_regions_{ifs_run_id}.txt"
                mask_labels = save_mask_regions_txt(labeled_mask, result_data + "/" + txt_filename)

                sg.Popup("Manual binned regions computed. Now click 'Extract!' to extract the binned spectra")

            except Exception as e:
                sg.Popup('Sorry, manual binning has failed and I do not know why. Maybe the datacube does not exist?')
                continue


        # preview mode for the manual binning: reading the datacube and showing the S/N of the spaxels contained in the selected regions
        if cube_ifs_event == 'Preview bins' and ifs_manual_bin:
            ifs_min_snr_mask_bin = 0 #No S/N cut
            ifs_bin_method_manual = 'False' #No voronoi binning
            voronoi_bin = False #No voronoi binning
            ifs_target_snr_manual = 0
            ifs_covariance_manual = 0

            # 2) extracting the label column of all the spaxels with the bin_info of the selected manual regions
            try:

                # Creating the dictionary to be passed to the cube_extract module
                config_manual = cubextr.buildConfigFromGUI(
                    ifs_run_id, ifs_input, ifs_output_dir, ifs_redshift,
                    ifs_lfs_data_default, ifs_ow_config,
                    ifs_ow_output, ifs_routine_selected, ifs_origin,
                    ifs_lmin_tot, ifs_lmax_tot, ifs_lmin_snr, ifs_lmax_snr,
                    ifs_min_snr_mask_bin, bin_mask_path, ifs_bin_method_manual, ifs_target_snr_manual,
                    ifs_covariance_manual, ifs_prepare_method)
                try:
                    cubextr.extract(config_manual, True, voronoi_bin, ifs_manual_bin, ifs_existing_bin)
                except Exception as e:
                    sg.popup("Sorry, Error.", str(e))
                    continue

            except Exception as e:
                sg.Popup('You first need to define the regions to be binned!\nOtherwise Select the Voronoi rebinning to automatically rebin the data')
                continue


        # NOW WE HAVE THE MAP WITH THE LABELED SPAXELS. Negative labels means spaxels not selected, therefore not considered. Positive labels identify the spaxels to consider for binning. Contiguous regions are marked with the same identifier (e.g. 1). This map has been stretched to 1D following the same order that the cubextr stores the spaxel infos in the _table.fit file. Now we need to generate the _table.fit file without any rebinning in order to have the BIN_ID of each spaxel, then we replace the BIN_ID array of the file with the bin info stored in the third component of the mask_labels array.

        # now we apply the manual bin by running the cube_extract_module in two steps
        if cube_ifs_event == 'Extract!': #and ifs_manual_bin:

            if ifs_manual_bin:
            # 1) RUNNING cubextract in preview mode without any rebinning to extract the info of the spaxels stored in the _table.fits file.
                ifs_min_snr_mask_bin = 0 #No S/N cut
                ifs_bin_method_manual = 'False'
                voronoi_bin = False #No voronoi binning
                ifs_target_snr_manual = 0
                ifs_covariance_manual = 0

            # 2) extracting the label column of all the spaxels with the bin_info of the selected manual regions
                try:
                    region_labels = mask_labels[:, 2].copy() #creating a copy, otherwise if exectuted more than one time is erodes the bin number
                    #Starting from BIN_ID zero and not one!
                    region_labels[region_labels > 0] -= 1
                except Exception as e:
                    sg.Popup('You first need to define the regions to be binned!\nOtherwise Select the Voronoi rebinning to automatically rebin the data')
                    continue

                # Creating the dictionary to be passed to the cube_extract module
                config_manual = cubextr.buildConfigFromGUI(
                    ifs_run_id, ifs_input, ifs_output_dir, ifs_redshift,
                    ifs_lfs_data_default, ifs_ow_config,
                    ifs_ow_output, ifs_routine_selected, ifs_origin,
                    ifs_lmin_tot, ifs_lmax_tot, ifs_lmin_snr, ifs_lmax_snr,
                    ifs_min_snr_mask_bin, bin_mask_path, ifs_bin_method_manual, ifs_target_snr_manual,
                    ifs_covariance_manual, ifs_prepare_method)
                try:
                    #running the cubextract module to produce the spaxel and BIN_ID map
                    cubextr.extract(config_manual, True, voronoi_bin, ifs_manual_bin, ifs_existing_bin)
                except Exception as e:
                    sg.popup("Error! Cannot show the bins", str(e))
                    continue


                # #3) REPLACE THE BIN_INFO IN THE _TABLE.FITS WITH THE LABELLED VALUES STORED IN region_labels
                fits_table_path = result_data + '/' + ifs_run_id + '/' + ifs_run_id + '_table.fits'

                # Opening fits
                with fits.open(fits_table_path, mode="update") as hdul:
                    data_table = hdul[1].data

                    # Checking
                    if len(data_table) != len(region_labels):
                        raise ValueError("Mismatch size between the labelled spaxels list and the actual spaxel list in the _table.fits file")

                    # Updating
                    data_table['BIN_ID'] = region_labels
                    hdul.flush()

                # 4) Now calculate the mean position and SNR of the spaxels to be binned
                with fits.open(fits_table_path, mode="update") as hdul:
                    data_hdu = hdul[1]
                    tbl = Table(data_hdu.data)  # Convert to Astropy Table for better handling

                    # Spaxels selected for binning
                    valid_mask = (tbl['BIN_ID'] >= 0)
                    unique_bins = np.unique(tbl['BIN_ID'][valid_mask])

                    # Calculating the positions
                    for b in unique_bins:
                        region_mask = (tbl['BIN_ID'] == b)

                        # Mean (NOT weighted) position of the bins
                        mean_x = np.mean(tbl['X'][region_mask])
                        mean_y = np.mean(tbl['Y'][region_mask])

                        # Spaxel number to be binned
                        n_spax = np.count_nonzero(region_mask)

                        # Calculate the S/N of the bins
                        flux_i = tbl['FLUX'][region_mask]
                        sn_i = tbl['SNR'][region_mask]
                        S_total = np.sum(flux_i)
                        noise_i = flux_i / sn_i
                        noise_quad_sum = np.sum(noise_i**2)
                        SNR_bin = S_total / np.sqrt(noise_quad_sum)

                        # Updating the values
                        tbl['XBIN'][region_mask] = mean_x
                        tbl['YBIN'][region_mask] = mean_y
                        tbl['NSPAX'][region_mask] = n_spax
                        tbl['SNRBIN'][region_mask] = SNR_bin

                    # updating
                    hdul[1].data = tbl.as_array()
                    hdul.flush()

                # 5) Run cubextract again with the new bin configuration
                try:
                    mock_voronoi = True # Fake voronoi bin required
                    cubextr.extract(config_manual, False, mock_voronoi, ifs_manual_bin, ifs_existing_bin)
                except Exception as e:
                    sg.Popup("ERROR performing the extraction")


            # With voronoi rebinning things are easier:
            if not ifs_manual_bin:

            # Creating the dictionary to be passed to the cube_extract module
                config = cubextr.buildConfigFromGUI(
                    ifs_run_id, ifs_input, ifs_output_dir, ifs_redshift,
                    ifs_lfs_data_default, ifs_ow_config,
                    ifs_ow_output, ifs_routine_selected, ifs_origin,
                    ifs_lmin_tot, ifs_lmax_tot, ifs_lmin_snr, ifs_lmax_snr,
                    ifs_min_snr_mask, ifs_mask, ifs_bin_method, ifs_target_snr,
                    ifs_covariance, ifs_prepare_method
    )


                print ('This might take a while. Please, relax...')

                try:
                    voronoi = True
                    preview = False
                    #calling the cube_extraction routine
                    cubextr.extract(config, preview, voronoi, ifs_manual_bin, ifs_existing_bin)
                except Exception as e:
                    sg.Popup ('ERROR performing the extraction')
                    continue


            #extracting the bin positions infos and saving in a txt file and in lists
            root_spectra_file_bin_info = result_data+'/'+ifs_run_id+'/'+ifs_run_id+'_table.fits'
            output_file_bin_data = result_data+'/'+ifs_run_id+'/'+ifs_run_id+'_bin_info.txt'

            try:
                with fits.open(root_spectra_file_bin_info) as hdul:
                    tbl = Table(hdul[1].data)
            except Exception as e:
                sg.Popup('Cannot read the datacube')
                continue

            # Select only binned spaxels
            valid_mask = (tbl['BIN_ID'] >= 0)
            unique_bins = np.unique(tbl['BIN_ID'][valid_mask])

            bin_id_array = []
            bin_x_array = []
            bin_y_array = []

            with open(output_file_bin_data, "w") as f:
                # Header
                f.write("#BIN_ID BIN_NUMBER XBIN YBIN SNRBIN NSPAX\n")

                # For all the bins
                for b in unique_bins:
                    region_mask = (tbl['BIN_ID'] == b)

                    # Taking the first index of the binned regions
                    idx_first = np.where(region_mask)[0][0]

                    bin_number = b + 1

                    # Extracting the values
                    bin_id   = b
                    bin_x    = tbl['XBIN'][idx_first]
                    bin_y    = tbl['YBIN'][idx_first]
                    bin_snr  = tbl['SNRBIN'][idx_first]
                    bin_nspx = tbl['NSPAX'][idx_first]

                    #Storing the interesting values in a list to be used later
                    bin_id_array.append(bin_id)
                    bin_x_array.append(bin_x)
                    bin_y_array.append(bin_y)

                    # writing to a file
                    f.write(f"{bin_id} {bin_number} {bin_x} {bin_y} {bin_snr} {bin_nspx}\n")

            print("Text file written with BIN info:", output_file_bin_data)


            #saving the extracted spectra also in single fits files SPAN-ready
            try:
                root_spectra_file = result_data+'/'+ifs_run_id+'/'+ifs_run_id+'_BinSpectra_linear.fits'
                hdul = fits.open(root_spectra_file)
                data_flux = hdul[1].data['SPEC']
                data_variance = hdul[1].data['ESPEC']
                wavelengths = hdul[2].data['WAVE']

                #creating the subdirectoy to store the single bins spectra
                single_bins_dir = result_data+'/'+ifs_run_id+'/bins'
                os.makedirs(single_bins_dir, exist_ok=True)
            except Exception:
                sg.popup('Wavelength interval not covered by the spectra!')
                continue

            #saving the spectra
            # writing the single spectra with 'BIN_ID' in the filename
            try:
                for i in range(data_flux.shape[0]):
                    flux = data_flux[i]
                    variance = data_variance[i]
                    t = Table([wavelengths, flux, variance], names=('wavelength', 'flux', 'variance'))

                    # Primary HDU and keywords
                    primary_hdu = fits.PrimaryHDU()
                    primary_hdu.header['BIN_ID'] = bin_id_array[i]
                    primary_hdu.header['X'] = bin_x_array[i]
                    primary_hdu.header['Y'] = bin_y_array[i]

                    # Creating bintable to store the data
                    table_hdu = fits.BinTableHDU(t)

                    # Craring the HDU
                    hdulist = fits.HDUList([primary_hdu, table_hdu])
                    filename = f"{single_bins_dir}/{ifs_run_id}_bin_id_{i:04}.fits"
                    hdulist.writeto(filename, overwrite=True)
            except Exception as e:
                sg.Popup('Results already present in the folder. Please, change the run_id name and try again')
                continue

            print('Single binned spectra saved in:', single_bins_dir, 'Wavelength units: A')

            #closing the fits file _BinSpectra_linear.
            hdul.close()

            #create spectra list of the bins to use with SPAN
            folder_path = single_bins_dir
            if folder_path:
                file_list = stm.get_files_in_folder(folder_path)
                output_file = ifs_run_id + '_bins_list.txt'
                stm.save_to_text_file(file_list, output_file)
                sg.Popup('Spectra file list of the bins saved in the working directory', output_file, 'You can now browse and load this list file\n\nWARNING: change the name of the run to process again')


        #showing the help file
        if cube_ifs_event == 'I need help':
            f = open(os.path.join(BASE_DIR, "help_files", "help_3d_spec.txt"), 'r')
            file_contents = f.read()
            if layout == layouts.layout_android:
                sg.popup_scrolled(file_contents, size=(120, 30))
            else:
                sg.popup_scrolled(file_contents, size=(100, 40))

    cube_ifs_window.close()

    #updating the parameters
    params = replace(params,
                    ifs_run_id = ifs_run_id,
                    ifs_input = ifs_input,
                    ifs_redshift = ifs_redshift,
                    ifs_lfs_data_default = ifs_lfs_data_default,
                    ifs_ow_config = ifs_ow_config,
                    ifs_ow_output = ifs_ow_output,
                    ifs_lmin_tot = ifs_lmin_tot,
                    ifs_lmax_tot = ifs_lmax_tot,
                    ifs_preloaded_routine = ifs_preloaded_routine,
                    ifs_min_snr_mask = ifs_min_snr_mask,
                    ifs_target_snr = ifs_target_snr,
                    ifs_routine_read = ifs_routine_read,
                    ifs_routine_read_default = ifs_routine_read_default,
                    ifs_user_routine = ifs_user_routine,
                    ifs_user_routine_file = ifs_user_routine_file,
                    ifs_origin = ifs_origin,
                    ifs_mask = ifs_mask,
                    ifs_lmin_snr_default = ifs_lmin_snr_default,
                    ifs_lmax_snr_default = ifs_lmax_snr_default,
                    ifs_manual_bin = ifs_manual_bin,
                    ifs_voronoi = ifs_voronoi,
                    ifs_existing_bin = ifs_existing_bin,
                    ifs_existing_bin_folder = ifs_existing_bin_folder,
                    ifs_bin_method = ifs_bin_method,
                    ifs_covariance = ifs_covariance,
                    ifs_prepare_method = ifs_prepare_method
                     )

    return params


# saving the spaxels for the Cube extract panel and manual bin info in a txt file and store in the array.
def save_mask_regions_txt(labeled_mask, output_filename):
    """
    - If labeled_mask[y,x] == 0 → label = -1 (not selected)
    - Otherwise label = labeled_mask[y,x]
    """
    rows, cols = labeled_mask.shape

    # Prepare a list to store the values (y, x, label)
    mask_labels_list = []

    with open(output_filename, "w") as f:
        f.write("# y\tx\tregion_label\n")
        for y in range(rows):
            for x in range(cols):
                lbl = labeled_mask[y, x]
                region_label = lbl if lbl != 0 else -1
                # Saving the text file
                f.write(f"{y}\t{x}\t{region_label}\n")
                # Fill the list
                mask_labels_list.append([y, x, region_label])

    # Converting the list to numpy and return it
    mask_labels = np.array(mask_labels_list, dtype=int)
    return mask_labels
