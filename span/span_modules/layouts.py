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

# Layout definitions for the main GUI for different OS environments. Modify this code to change the aspect of the main GUI.

try: #try local import if executed as script
    #GUI import
    from FreeSimpleGUI_local import FreeSimpleGUI as sg
except ModuleNotFoundError: #local import if executed as package
    #GUI import
    from span.FreeSimpleGUI_local import FreeSimpleGUI as sg

import os

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(CURRENT_DIR)

################ FreeSimpleGUI User Interface construction ####################
listbox1 = ['Load a spectra file list and click Load!']
default_spectra_list = os.path.join(BASE_DIR, "example_files", "xshooter_vis_sample_list_spectra.dat")

sg.SetOptions(tooltip_time=1000) #tooltip time after mouse over
sg.theme('DarkBlue3')
menu_def = ['Unused', ['↑ Move Up', '↓ Move Down', 'Remove', 'Undo', '---', 'Compare spectra']] # right click event menu for the listbox
#************************************************************************************
#************************************************************************************

#Layout optimized for Windows systems
layout_windows = [
            [sg.Menu([
                ['&File', ['&Load!', '&Save parameters...', 'Load parameters...', 'Restore default parameters', 'E&xit']],
                ['&Edit', ['Clear all tas&ks', 'Clean output', 'Show result folder', 'Change result folder...', 'Save current spectra list...']],
                ['&Window', ['Long-slit extraction', 'DataCube extraction', 'Text editor', 'FITS header editor', 'Spectra manipulation', 'Utilities']],
                ['P&rocess',['Pl&ot', 'Pre&view spec.']],
                ['&Analysis', ['Preview res&ult', 'Proc&ess selected', 'Process a&ll']],
                ['&Plotting', ['Plot data', 'Plot maps']],
                ['&View', ['Zoom In', 'Zoom Out', 'Reset Zoom']],
                ['&Help', ['&Quick start', '&Read me', 'Tips and tricks', 'SPAN Manual']],
                ['&About', ['About SPAN', 'Version', 'Read me']]
                ])
            ],

            [sg.Frame('Prepare and load spectra', [
            [sg.Text('1. Extract 1D spectra from 2D or 3D FITS images', font = ('', 11 ,'bold'))],
            [sg.Button('Long-slit extraction', tooltip='Stand alone program to extract 1D spectra from 2D fits',button_color= ('black','light blue')), sg.Button('DataCube extraction', tooltip='Stand alone program to extract 1D spectra from data cubes',button_color= ('black','light blue'))],
            [sg.HorizontalSeparator(pad=(0, 8))],
            [sg.Text('2. Generate a spectra list containing 1D spectra', font = ('', 11 ,'bold'))],
            [sg.Button('Generate spectra list containing 1D spectra', key = 'listfile',tooltip='If you do not have a spectra file list, you can generate here',size = (37,2))],
            [sg.HorizontalSeparator(pad=(0, 8))],
            [sg.Text('3. Browse the spectra list or just one spectrum', font = ('', 11 ,'bold'))],
            [sg.InputText(default_spectra_list, size=(34, 1), key='spec_list' ), sg.FileBrowse(tooltip='Load an ascii file list of spectra or a single (fits, txt) spectrum')],
            [sg.Checkbox('I browsed a single spectrum', font = ('Helvetica', 11, 'bold'), key='one_spec',tooltip='Check this if you want to load just one spectrum instead a text file containing the names of the spectra')],
            [sg.Text('Wavelength units:',tooltip='Set the correct wavelength units of your spectra: Angstrom, nm, mu', font = ('', 12)), sg.Radio('nm', "RADIO2", default=True, key = 'wave_units_nm' ), sg.Radio('A', "RADIO2", key = 'wave_units_a'), sg.Radio('mu', "RADIO2" , key = 'wave_units_mu')],
            [sg.HorizontalSeparator(pad=(0, 8))],
            [sg.Text('4. Finally load the browsed spectra to SPAN', font = ('', 11 ,'bold'))],
            [sg.Button('Load!', font = ("Helvetica", 11, 'bold'),button_color=('black','light green'), size = (11,1)), sg.Push(), sg.Button('Plot',button_color=('black','light gray'), size = (10,1))],
            ], font=("Helvetica", 14, 'bold'), title_color = 'orange'),

            sg.Frame('Loaded spectra', [[
            sg.Listbox(values=listbox1,size=(46, 20),key='-LIST-',horizontal_scroll=True,enable_events=True,right_click_menu=menu_def,select_mode=sg.LISTBOX_SELECT_MODE_EXTENDED)]], font=("Helvetica", 12, 'bold'), title_color='white'),
            
            sg.Frame('Preview', [
            [sg.Canvas(key='-CANVAS-', size=(580, 320))],
            [sg.Button('Header',   key='-SHOWHDR-',  button_color=('black','light gray'), size=(10,1)), sg.Button('Sampling', key='-SHOWSTEP-', button_color=('black','light gray'), size=(10,1)), sg.Push(), sg.Button('Spectra manipulation', size=(18,1), button_color=('black','light blue')),
            sg.Button('Preview spec.', button_color=('black','light gray'), size=(13,1))],
            ], font=("Helvetica", 12, 'bold'))],

            #Spectral analysis frame
            [sg.Frame('Spectral analysis', [

            #1) Black-body fitting
            [sg.Checkbox('Planck Blackbody fitting', font = ('Helvetica', 12, 'bold'), key = 'bb_fitting',tooltip='Blackdoby Planck function fitting. Works fine for stellar spectra and wide wavelength range'),sg.Push(), sg.Button('Blackbody parameters',button_color= ('black','light blue'), size = (22,1))],

            #2) Cross-correlation
            [sg.Checkbox('Cross-correlation', font = ('Helvetica', 12, 'bold'), key = 'xcorr',tooltip='Cross-correlating a band with a template. Use Stars and gas Kinematics to refine the value found'),sg.Push(), sg.Button('Cross-corr parameters',button_color= ('black','light blue'), size = (22,1))],

            #3) Velocity disperion measurement
            [sg.Checkbox('Velocity dispersion', font = ('Helvetica', 12, 'bold'), key = 'sigma_measurement',tooltip='Fitting a band with a template. Rough but fast. Use Kinematics for accurate science results'),sg.Push(), sg.Button('Sigma parameters',button_color= ('black','light blue'), size = (22,1))],

            #4) Line fitting
            [sg.Checkbox('Line(s) fitting', font = ('Helvetica', 12, 'bold'), key = 'line_fitting',tooltip='User line or automatic CaT band fitting with gaussian functions'),sg.Push(), sg.Button('Line fitting parameters',button_color= ('black','light blue'), size = (22,1))],
            
            #5) Line-strength
            [sg.Checkbox('Line-strength analysis', font = ('Helvetica', 12, 'bold'), key = 'ew_measurement',tooltip='Equivalent width measurement for a list of indices, a single user defined index and Lick/IDS indices'),sg.Push(), sg.Button('Line-strength parameters',button_color= ('black','light blue'), size = (22,1))],

            #6) Kinematics with ppxf
            [sg.Checkbox('Stars and gas kinematics', font = ('Helvetica', 12, 'bold'), key = 'ppxf_kin',tooltip='Perform the fitting of a spectral region and gives the kinematics'),sg.Push(), sg.Button('Kinematics parameters',button_color= ('black','light blue'), size = (22,1))  ],

            #7) Stellar populations with ppxf
            [sg.Checkbox('Stellar populations and SFH', font = ('Helvetica', 12, 'bold'), key = 'ppxf_pop',tooltip='Perform the fitting of a spectral region and gives the properties of the stellar populations'),sg.Push(), sg.Button('Population parameters',button_color= ('black','light blue'), size = (22,1))  ],
            ], font=("Helvetica", 14, 'bold'), title_color='yellow'),

            # Buttons to perform the spectral analysis actions
            sg.Frame('An. Actions',[
            [sg.Button('Help me',button_color=('black','orange'), size = (12,1))],
            [sg.Text('')],
            [sg.Button('Preview result',button_color=('black','light gray'),tooltip='Preview all the results of the Spectral analysis frame', size = (12,2))],
            [sg.Text('')],
            [sg.Text('')],
            [sg.Text('')],
            [sg.Text('')],
            [sg.Text('', font = ('Helvetica',16))],
            ],font=("Helvetica", 9, 'bold')),

            #COMMENT THE FOLLOWING THREE LINES TO HAVE THE EXTERNAL OUTPUT
            sg.Frame('Output', [
            [sg.Output(size=(94, 14), key='-OUTPUT-', font = ('Helvetica', 11))],
            ] ,font=("Helvetica", 12, 'bold')),

            ],

            #General buttons at the end of the panel
            [sg.Button('Process selected', button_color=('white','orange'), size=(15, 1),tooltip='Process the selected spectrum by performing all the enabled tasks'), sg.Button('Process all', button_color=('white','red'), size=(15, 1), tooltip='Process all the loaded spectra by performing all the enabled tasks'), sg.Checkbox('Save spectral analysis plots', default = False, text_color='yellow', key = 'save_plots', tooltip='To save all the plots generated by the Spectral Analysis tasks activated and the Process All mode', font = ("Helvetica", 10, 'bold')), sg.Push(), sg.Exit(size=(15, 1),tooltip='See you soon!')]

                ]


#************************************************************************************
#************************************************************************************
#Layout optimized for Linux systems
layout_linux = [
            [sg.Menu([
                ['&File', ['&Load!', '&Save parameters...', 'Load parameters...', 'Restore default parameters', 'E&xit']],
                ['&Edit', ['Clear all tas&ks', 'Clean output', 'Show result folder', 'Change result folder...', 'Save current spectra list...']],
                ['&Window', ['Long-slit extraction', 'DataCube extraction', 'Text editor', 'FITS header editor', 'Spectra manipulation', 'Utilities']],
                ['P&rocess',['Pl&ot', 'Pre&view spec.']],
                ['&Analysis', ['Preview res&ult', 'Proc&ess selected', 'Process a&ll']],
                ['&Plotting', ['Plot data', 'Plot maps']],
                ['&View', ['Zoom In', 'Zoom Out', 'Reset Zoom']],
                ['&Help', ['&Quick start', '&Read me', 'Tips and tricks', 'SPAN Manual']],
                ['&About', ['About SPAN', 'Version', 'Read me']]
                ])
            ],

            [sg.Frame('Prepare and load spectra', [
            [sg.Text('1. Extract spectra from 2D or 3D FITS', font = ('', 11 ,'bold'))],
            [sg.Button('Long-slit extraction', tooltip='Stand alone program to extract 1D spectra from 2D fits',button_color= ('black','light blue')), sg.Button('DataCube extraction', tooltip='Stand alone program to extract 1D spectra from data cubes',button_color= ('black','light blue'))],
            #[sg.Text('', font = ('', 1))],
            [sg.HorizontalSeparator(pad=(0, 10))],
            [sg.Text('2. Generate a spectra list with 1D spectra', font = ('', 11 ,'bold'))],
            [sg.Button('Generate spectra list containing 1D spectra', key = 'listfile',tooltip='If you do not have a spectra file list, you can generate here',size = (37,2))],
            [sg.HorizontalSeparator(pad=(0, 9))],
            [sg.Text('3. Browse the spectra list or one spectrum', font = ('', 11 ,'bold'))],
            [sg.InputText(default_spectra_list, size=(36, 1), key='spec_list' ), sg.FileBrowse(tooltip='Load an ascii file list of spectra or a single (fits, txt) spectrum')],
            [sg.Checkbox('I browsed a single spectrum', font = ('Helvetica', 11, 'bold'), key='one_spec',tooltip='Check this if you want to load just one spectrum instead a text file containing the names of the spectra')],
            [sg.Text('Wavelength units:', font = ('', 11),tooltip='Set the correct wavelength units of your spectra: Angstrom, nm, mu'), sg.Radio('nm', "RADIO2", default=True, key = 'wave_units_nm' ), sg.Radio('A', "RADIO2", key = 'wave_units_a'), sg.Radio('mu', "RADIO2" , key = 'wave_units_mu')],
            [sg.HorizontalSeparator(pad=(0, 10))],
            [sg.Text('4. Finally load the spectra to SPAN', font = ('', 11 ,'bold'))],
            [sg.Button('Load!', font = ("Helvetica", 11, 'bold'),button_color=('black','light green'), size = (11,1)), sg.Push(), sg.Button('Plot',button_color=('black','light gray'), size = (10,1))],
            ], font=("Helvetica", 14, 'bold'), title_color = 'orange'),

            sg.Frame('Loaded spectra', [[
            sg.Listbox(values=listbox1,size=(50, 24),key='-LIST-',horizontal_scroll=True,enable_events=True,right_click_menu=menu_def,select_mode=sg.LISTBOX_SELECT_MODE_EXTENDED)]], font=("Helvetica", 12, 'bold'), title_color='white'),

            sg.Frame('Preview', [
            [sg.Canvas(key='-CANVAS-', size=(680, 345))],
            [sg.Button('Header',   key='-SHOWHDR-',  button_color=('black','light gray'), size=(10,1)), sg.Button('Sampling', key='-SHOWSTEP-', button_color=('black','light gray'), size=(10,1)), sg.Push(), sg.Button('Spectra manipulation', size=(18,1), button_color=('black','light blue')),
            sg.Button('Preview spec.', button_color=('black','light gray'), size=(13,1))],
            ], font=("Helvetica", 12, 'bold'))],
            
            #Spectral analysis frame
            [sg.Frame('Spectral analysis', [

            #1) Black-body fitting
            [sg.Checkbox('Planck Blackbody fitting', font = ('Helvetica', 12, 'bold'), key = 'bb_fitting',tooltip='Blackdoby Planck function fitting. Works fine for stellar spectra and wide wavelength range'),sg.Push(), sg.Button('Blackbody parameters',button_color= ('black','light blue'), size = (22,1))],

            #2) Cross-correlation
            [sg.Checkbox('Cross-correlation', font = ('Helvetica', 12, 'bold'), key = 'xcorr',tooltip='Cross-correlating a band with a template. Use Stars and gas Kinematics to refine the value found'),sg.Push(), sg.Button('Cross-corr parameters',button_color= ('black','light blue'), size = (22,1))],

            #3) Velocity disperion measurement
            [sg.Checkbox('Velocity dispersion', font = ('Helvetica', 12, 'bold'), key = 'sigma_measurement',tooltip='Fitting a band with a template. Rough but fast. Use Kinematics for accurate science results'),sg.Push(), sg.Button('Sigma parameters',button_color= ('black','light blue'), size = (22,1))],

            #4) Line fitting
            [sg.Checkbox('Line(s) fitting', font = ('Helvetica', 12, 'bold'), key = 'line_fitting',tooltip='User line or automatic CaT band fitting with gaussian functions'),sg.Push(), sg.Button('Line fitting parameters',button_color= ('black','light blue'), size = (22,1))],
            
            #5) Line-strength
            [sg.Checkbox('Line-strength analysis', font = ('Helvetica', 12, 'bold'), key = 'ew_measurement',tooltip='Equivalent width measurement for a list of indices, a single user defined index and Lick/IDS indices'),sg.Push(), sg.Button('Line-strength parameters',button_color= ('black','light blue'), size = (22,1))],

            #6) Kinematics with ppxf
            [sg.Checkbox('Stars and gas kinematics', font = ('Helvetica', 12, 'bold'), key = 'ppxf_kin',tooltip='Perform the fitting of a spectral region and gives the kinematics'),sg.Push(), sg.Button('Kinematics parameters',button_color= ('black','light blue'), size = (22,1))  ],

            #7) Stellar populations with ppxf
            [sg.Checkbox('Stellar populations and SFH', font = ('Helvetica', 12, 'bold'), key = 'ppxf_pop',tooltip='Perform the fitting of a spectral region and gives the properties of the stellar populations'),sg.Push(), sg.Button('Population parameters',button_color= ('black','light blue'), size = (22,1))  ],
            ], font=("Helvetica", 14, 'bold'), title_color='yellow'),

            # Buttons to perform the spectral analysis actions
            sg.Frame('An. Actions',[
            [sg.Button('Help me',button_color=('black','orange'), size = (12,1))],
            [sg.Text('')],
            [sg.Button('Preview result',button_color=('black','light gray'),tooltip='Preview all the results of the Spectral analysis frame', size = (12,2))],
            [sg.Text('')],
            [sg.Text('')],
            [sg.Text('')],
            [sg.Text('')],
            [sg.Text('', font = ('Helvetica',24))],

            ],font=("Helvetica", 9, 'bold')),

            #COMMENT THE FOLLOWING THREE LINES TO HAVE THE EXTERNAL OUTPUT
            sg.Frame('Output', [
            [sg.Output(size=(104, 16), key='-OUTPUT-' , font=('Helvetica', 11))],
            ] ,font=("Helvetica", 12, 'bold')),

            ],

            #General buttons at the end of the panel
            [sg.Button('Process selected', button_color=('white','orange'), size=(15, 1),tooltip='Process the selected spectrum by performing all the enabled tasks'), sg.Button('Process all', button_color=('white','red'), size=(15, 1), tooltip='Process all the loaded spectra by performing all the enabled tasks'), sg.Checkbox('Save spectral analysis plots', default = False, text_color='yellow', key = 'save_plots', tooltip='To save all the plots generated by the Spectral Analysis tasks activated and the Process All method', font = ("Helvetica", 10, 'bold')), sg.Push(), sg.Exit(size=(15, 1),tooltip='See you soon!')]

                ]


#************************************************************************************
#************************************************************************************
#Layout optimized for MacOS systems
layout_macos = [
            [sg.Menu([
                ['&File', ['&Load!', '&Save parameters...', 'Load parameters...', 'Restore default parameters', 'E&xit']],
                ['&Edit', ['Clear all tas&ks', 'Show result folder', 'Change result folder...', 'Save current spectra list...']],
                ['&Window', ['Long-slit extraction', 'DataCube extraction', 'Text editor', 'FITS header editor', 'Spectra manipulation', 'Utilities']],
                ['P&rocess',['Pl&ot', 'Pre&view spec.']],
                ['&Analysis', ['Preview res&ult', 'Proc&ess selected', 'Process a&ll']],
                ['&Plotting', ['Plot data', 'Plot maps']],
                ['&View', ['Zoom In', 'Zoom Out', 'Reset Zoom']],
                ['&Help', ['&Quick start', '&Read me', 'Tips and tricks', 'SPAN Manual']],
                ['&About', ['About SPAN', 'Version', 'Read me']]
                ])
            ],

            [sg.Frame('Prepare and load spectra', [
            [sg.Text('1. Extract spectra from 2D or 3D FITS', font = ('', 14 ,'bold'))],
            [sg.Button('Long-slit extraction', tooltip='Stand alone program to extract 1D spectra from 2D fits',button_color= ('black','light blue')), sg.Button('DataCube extraction', tooltip='Stand alone program to extract 1D spectra from data cubes',button_color= ('black','light blue'))],
            [sg.HorizontalSeparator()],
            [sg.Text('2. Generate a 1D spectra list', font = ('', 14 ,'bold'))],
            [sg.Button('Generate spectra list containing 1D spectra', key = 'listfile',tooltip='If you do not have a spectra file list, you can generate here')],
            [sg.HorizontalSeparator()],
            [sg.Text('3. Browse the list or one spectrum', font = ('', 14 ,'bold'))],
            [sg.InputText(default_spectra_list, size=(34, 1), key='spec_list' ), sg.FileBrowse(tooltip='Load an ascii file list of spectra or a single (fits, txt) spectrum')],
            [sg.Checkbox('I browsed a single spectrum', key='one_spec',tooltip='Check this if you want to load just one spectrum instead a text file containing the names of the spectra')],
            [sg.Text('Wavelength of the spectra is in:',tooltip='Set the correct wavelength units of your spectra: Angstrom, nm, mu', font = ('Helvetica', 14)), sg.Radio('nm', "RADIO2", default=True, key = 'wave_units_nm' ), sg.Radio('A', "RADIO2", key = 'wave_units_a'), sg.Radio('mu', "RADIO2" , key = 'wave_units_mu')],
            [sg.HorizontalSeparator()],
            [sg.Text('4. Finally load the spectra to SPAN', font = ('', 14 ,'bold'))],
            [sg.Button('Load!', font = ("Helvetica", 14, 'bold'),button_color=('black','light green'), size = (11,1)), sg.Push(), sg.Button('Plot',button_color=('black','light gray'), size = (10,1))],
            ], font=("Helvetica", 18, 'bold'), title_color = 'orange'),

            sg.Frame('Loaded spectra', [[
            sg.Listbox(values=listbox1,size=(45, 19),key='-LIST-',horizontal_scroll=True,enable_events=True,right_click_menu=menu_def,select_mode=sg.LISTBOX_SELECT_MODE_EXTENDED)]], font=("Helvetica", 14, 'bold'), title_color='white'),

            sg.Frame('Preview', [
            [sg.Canvas(key='-CANVAS-', size=(660, 310))],
            [sg.Button('Header',   key='-SHOWHDR-',  button_color=('black','light gray'), size=(10,1)), sg.Button('Sampling', key='-SHOWSTEP-', button_color=('black','light gray'), size=(10,1)), sg.Push(), sg.Button('Spectra manipulation', size=(18,1), button_color=('black','light blue')),
            sg.Button('Preview spec.', button_color=('black','light gray'), size=(13,1))],
            ], font=("Helvetica", 14, 'bold'))],
            
            #Spectral analysis frame
            [sg.Frame('Spectral analysis', [

            #1) Black-body fitting
            [sg.Checkbox('Planck Blackbody fitting ', font = ('Helvetica', 16, 'bold'), key = 'bb_fitting',tooltip='Blackdoby Planck function fitting. Works fine for stellar spectra and wide wavelength range'), sg.Text('    '), sg.Button('Blackbody parameters',button_color= ('black','light blue'), size = (22,1)), sg.Text('          '), sg.Checkbox('Cross-correlation', font = ('Helvetica', 16, 'bold'), key = 'xcorr',tooltip='Cross-correlating a band with a template. Use Stars and gas Kinematics to refine the value found'),sg.Push(), sg.Button('Cross-corr parameters',button_color= ('black','light blue'), size = (22,1))],

            #2) Velocity dispersion measurement
            [sg.Checkbox('Velocity dispersion   ', font = ('Helvetica', 16, 'bold'), key = 'sigma_measurement',tooltip='Fitting a band with a template. Rough but fast. Use Kinematics for accurate science results'),sg.Text('              '), sg.Button('Sigma parameters',button_color= ('black','light blue'), size = (22,1)), sg.Text('          '), sg.Checkbox('Line(s) fitting', font = ('Helvetica', 16, 'bold'), key = 'line_fitting',tooltip='User line or automatic CaT band fitting with gaussian functions'),sg.Push(), sg.Button('Line fitting parameters',button_color= ('black','light blue'), size = (22,1))],

            #3) Line-strength
            [sg.Checkbox('Line-strength analysis  ', font = ('Helvetica', 16, 'bold'), key = 'ew_measurement',tooltip='Equivalent width measurement for a list of indices, a single user defined index and Lick/IDS indices'), sg.Text('        '), sg.Button('Line-strength parameters',button_color= ('black','light blue'), size = (22,1)), sg.Text('          '),sg.Checkbox('Kinematics', font = ('Helvetica', 16, 'bold'), key = 'ppxf_kin',tooltip='Perform the fitting of a spectral region and gives the kinematics'),sg.Push(), sg.Button('Kinematics parameters',button_color= ('black','light blue'), size = (22,1))  ],

            #4) Stellar populations with ppxf
            [sg.Checkbox('Stellar populations and SFH ', font = ('Helvetica', 16, 'bold'), key = 'ppxf_pop',tooltip='Perform the fitting of a spectral region and gives the properties of the stellar populations'), sg.Button('Population parameters',button_color= ('black','light blue'), size = (22,1))  ],
            ], font=("Helvetica", 18, 'bold'), title_color='yellow'),

            # Buttons to perform the spectral analysis actions
            sg.Frame('An. Actions',[
            [sg.Button('Help me',button_color=('black','orange'), size = (12,1))],
            [sg.Text('')],
            [sg.Button('Preview result',button_color=('black','light gray'),tooltip='Preview all the results of the Spectral analysis frame', size = (12,1), font=("Helvetica", 14, 'bold'))],

            ],font=("Helvetica", 14, 'bold')),


            ],
            [sg.HorizontalSeparator()],

            #General buttons at the end of the panel
            [sg.Button('Process selected', button_color=('white','orange'), size=(15, 1),tooltip='Process the selected spectrum by performing all the enabled tasks'), sg.Button('Process all', button_color=('white','red'), size=(15, 1), tooltip='Process all the loaded spectra by performing all the enabled tasks'), sg.Checkbox('Save spectral analysis plots', default = False, text_color='yellow', key = 'save_plots', tooltip='To save all the plots generated by the Spectral Analysis tasks activated and the Process All method', font = ("Helvetica", 14, 'bold')), sg.Push(), sg.Exit(size=(15, 1),tooltip='See you soon!')]

                ]


#************************************************************************************
#Layout optimized for Android systems
layout_android = [
            [sg.Button('Read me', button_color=('black','orange'), tooltip='Open the SPAN readme'), sg.Button('Quick start', button_color=('black','orange'), tooltip='A fast guide to begin using SPAN'), sg.Button('Tips and tricks', button_color=('black','orange'), tooltip='Some tricks to master SPAN'), sg.Push(), sg.Button('Change result folder...', button_color=('black','light blue')), sg.Button('Save parameters...', button_color=('black','light blue'), tooltip='Save the current parameters in a json file'), sg.Button('Load parameters...', button_color=('black','light blue'), tooltip='Load the saved parameters'), sg.Button('Restore default parameters', button_color=('black','light blue'), tooltip='Restore the default parameters'), sg.Button('Clear all tasks', button_color=('black','light blue'), tooltip='De-activate all the tasks, including from the spectral manipulation panel'), sg.Button('Clean output', button_color=('black','light blue'), tooltip='Delete the output window')],
            [sg.HorizontalSeparator()],

            [sg.Frame('Prepare and load spectra', [
            [sg.Text('1. Extract 1D spectra and/or generate a spectra list', font = ('Helvetica', 11, 'bold'))],
            [sg.Button('Long-slit extraction', tooltip='Stand alone program to extract 1D spectra from 2D fits',button_color= ('black','light blue'), size=(13, 2)), sg.Button('DataCube extraction', tooltip='Stand alone program to extract 1D spectra from data cubes',button_color= ('black','light blue'), size=(11, 2)), sg.Button('Gen. spectra list', key = 'listfile',tooltip='If you do not have a spectra file list, you can generate here', size=(14, 2))],
            [sg.Text('', font = ("Helvetica", 1))],
            [sg.Text('2. Browse the spectra list or just one spectrum', font = ('Helvetica', 11, 'bold'))],
            [sg.InputText(default_spectra_list, size=(39, 1), key='spec_list' ), sg.FileBrowse(tooltip='Load an ascii file list of spectra or a single (fits, txt) spectrum')],
            [sg.Checkbox('I browsed a single spectrum', font = ('Helvetica', 10, 'bold'), key='one_spec',tooltip='Check this if you want to load just one spectrum instead a text file containing the names of the spectra')],
            [sg.Text('W. scale:',tooltip='Set the correct wavelength units of your spectra: Angstrom, nm, mu'), sg.Radio('nm', "RADIO2", default=True, key = 'wave_units_nm' ), sg.Radio('A', "RADIO2", key = 'wave_units_a'), sg.Radio('mu', "RADIO2" , key = 'wave_units_mu'), sg.Push(), sg.Button('Load!', font = ('Helvetica', 11, 'bold'),button_color=('black','light green'), size = (6,1)), sg.Button('Plot',button_color=('black','light gray'), size = (4,1))],
            ], font=("Helvetica", 14, 'bold'), title_color = 'orange'),

            sg.Frame('Loaded spectra', [[
            sg.Listbox(values=listbox1,size=(42, 11),key='-LIST-',horizontal_scroll=True,enable_events=True,select_mode=sg.LISTBOX_SELECT_MODE_EXTENDED)]], font=("Helvetica", 12, 'bold'), title_color='white'),

            sg.Frame('Preview', [
            [sg.Canvas(key='-CANVAS-', size=(840, 275))],
            [sg.Button('Header',   key='-SHOWHDR-',  button_color=('black','light gray'), size=(8,1)), sg.Button('Sampling', key='-SHOWSTEP-', button_color=('black','light gray'), size=(8,1)), sg.Button('Utilities', button_color=('black','light gray'), size=(8,1))],
            ], font=("Helvetica", 12, 'bold'))],

            [sg.Frame('Spectral analysis', [

            #1) Black-body fitting
            [sg.Checkbox('Planck Blackbody fitting', font = ('Helvetica', 12, 'bold'), key = 'bb_fitting',tooltip='Blackdoby Planck function fitting. Works fine for stellar spectra and wide wavelength range'),sg.Push(), sg.Button('Blackbody parameters',button_color= ('black','light blue'), size = (22,1))],

            #2) Cross-correlation
            [sg.Checkbox('Cross-correlation', font = ('Helvetica', 12, 'bold'), key = 'xcorr',tooltip='Cross-correlating a band with a template. Use Stars and gas Kinematics to refine the value found'),sg.Push(), sg.Button('Cross-corr parameters',button_color= ('black','light blue'), size = (22,1))],

            #3) Velocity disperion measurement
            [sg.Checkbox('Velocity dispersion', font = ('Helvetica', 12, 'bold'), key = 'sigma_measurement',tooltip='Fitting a band with a template. Rough but fast. Use Kinematics for accurate science results'),sg.Push(), sg.Button('Sigma parameters',button_color= ('black','light blue'), size = (22,1))],

            #4) Line fitting
            [sg.Checkbox('Line(s) fitting', font = ('Helvetica', 12, 'bold'), key = 'line_fitting',tooltip='User line or automatic CaT band fitting with gaussian functions'),sg.Push(), sg.Button('Line fitting parameters',button_color= ('black','light blue'), size = (22,1))],
            
            #4) Line-strength
            [sg.Checkbox('Line-strength analysis', font = ('Helvetica', 12, 'bold'), key = 'ew_measurement',tooltip='Equivalent width measurement for a list of indices, a single user defined index and Lick/IDS indices'),sg.Push(), sg.Button('Line-strength parameters',button_color= ('black','light blue'), size = (22,1))],

            #6) Kinematics with ppxf
            [sg.Checkbox('Stars and gas kinematics', font = ('Helvetica', 12, 'bold'), key = 'ppxf_kin',tooltip='Perform the fitting of a spectral region and gives the kinematics'),sg.Push(), sg.Button('Kinematics parameters',button_color= ('black','light blue'), size = (22,1))  ],

            #7) Stellar populations with ppxf
            [sg.Checkbox('Stellar populations and SFH', font = ('Helvetica', 12, 'bold'), key = 'ppxf_pop',tooltip='Perform the fitting of a spectral region and gives the properties of the stellar populations'),sg.Push(), sg.Button('Population parameters',button_color= ('black','light blue'), size = (22,1))  ],
            ], font=("Helvetica", 14, 'bold'), title_color='yellow'),

            # Buttons to open the spectral manipulation panel and perform the spectral analysis actions
            sg.Frame('Actions',[
            [sg.Button('Spectra manipulation', size = (12,2),button_color= ('black','light blue'), font=("Helvetica", 10, 'bold'), tooltip='Open the spectra manipulation panel to modify the spectra', key = 'Spectra manipulation')],
            [sg.Button('Preview spec.',button_color=('black','light gray'), size = (12,2), font=("Helvetica", 10, 'bold'),tooltip='View the modified version of the selected spectrum')],
            [sg.Text('', font = ('Helvetica',1))],
            [sg.HorizontalSeparator()],
            [sg.Text('', font = ('Helvetica',1))],
            [sg.Button('Preview result',button_color=('black','light gray'),tooltip='Preview the results of the Spectral analysis frame', size = (12,2), font=("Helvetica", 10, 'bold'))],
            [sg.Text('')],
            [sg.Button('Help me',button_color=('black','orange'), size = (12,1),tooltip='Getting help for the spectral analysis')],

            ],font=("Helvetica", 9, 'bold')),

            #COMMENT THE FOLLOWING THREE LINES TO HAVE THE EXTERNAL OUTPUT
            sg.Frame('Output', [
            [sg.Output(size=(79, 12), key='-OUTPUT-' , font=('Helvetica', 11))],
            ] ,font=("Helvetica", 12, 'bold')),

            ],

            [sg.Button('Process selected', button_color=('white','orange'), size=(15, 2),tooltip='Process the selected spectrum by performing all the enabled tasks'), sg.Button('Process all', button_color=('white','red'), size=(15, 2), tooltip='Process all the loaded spectra by performing all the enabled tasks'), sg.Checkbox('Save spectral analysis plots', default = False, text_color='yellow', key = 'save_plots', tooltip='To save all the plots generated by the Spectral Analysis tasks activated and the Process All method', font = ("Helvetica", 10, 'bold')), sg.Push(), sg.Button('Text editor', tooltip='Stand alone simple text editor',button_color= ('black','light blue'),size =(8,1)),sg.Button('FITS header editor', tooltip='Stand alone FITS header editor',button_color= ('black','light blue'), size = (13,1)), sg.Button('Plot data', tooltip='Stand alone data plotter. ASCII files with spaced rows',button_color= ('black','light blue'), size = (7,1)), sg.Button('Plot maps', tooltip='Stand alone datacube maps plotter',button_color= ('black','light blue'), size =(8,1)), sg.Exit(size=(15, 2),tooltip='See you soon!')]

                ]
