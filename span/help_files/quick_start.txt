SPAN: SPectral ANalysis software V6.6
Daniele Gasparri, June 2025


### Quick start ###

To get started, press the "Load!" button in the "Prepare and Load Spectra" frame to load the example files.

The loaded spectra will appear in the central frame on the right (the white window). Simply select a spectrum with the mouse, then click "Plot" to display it. Close the plot window to re-enable the main panel.
You can now:
- Analyse the selected spectrum by activating any spectral analysis tasks.
- Modify the spectrum by opening the "Spectra Manipulation" panel, located in the "Utilities" frame on the right.


Applying a Spectral Manipulation Task
Let's try modifying a spectrum:
- Open the "Spectra Manipulation" panel.
- Activate one of the available tasks, e.g., "Add noise".
- Confirm the selection by pressing "Confirm".
- Press "Preview Spec." in the main panel to view the result.
- If satisfied, save the modified spectrum by clicking "Process Selected" (ensure the plot window is closed first!).
- To apply the selected task to all loaded spectra, click "Process All". 
The results will be saved in the "SPAN_results/spec" folder, located within the directory you selected when first opening SPAN.

Tip: The output window provides real-time updates on all tasks performed by the program.


Performing a Spectral Analysis Task
Now, let's try a spectral analysis task:
- In the "Spectral Analysis" frame, activate "Line-Strength Analysis".
- Click "Line-Strength Parameters" to view the available settings.
- Select the "Single Index" option and confirm by clicking "Confirm". The Line-Strength Analysis window will close automatically, returning you to the main panel.
- Click "Preview Result" to see the outcome.

Important:
If the spectrum looks strange, check whether the "Add Noise" task is still active.
The Spectral Analysis frame processes spectra after applying any selected Spectra Manipulation tasks.
If multiple tasks are active, the processed spectrum will reflect all applied modifications.
If no tasks are selected, the original spectrum is used.


Resetting and Managing Parameters
If too many tasks are activated and the program becomes difficult to manage, don't worry!
- Click "Edit → Clear All Tasks" to reset all active tasks.
- To restore default parameters, select "File → Restore Default Parameters".
- To save your custom settings, go to "File → Save Parameters", allowing you to reload them later.
- To change the SPAN_results folder location, use "Edit → Change Result Folder...".


Exploring Additional Sample Spectra
You can experiment with other preloaded 1D spectra lists provided with SPAN.
For example:
- In "3. Browse the spectra list or just one spectrum", load the spectra list "ngc5806_bins.dat" in the "example_files" folder. This file contains 39 optical spectra of the galaxy NGC 5806.
- Set the Wavelength units to "A" in the "Wavelength of the spectra is in:" option
- Press the "Load!" button to load the 39 1D spectra in the listbox. 
- Have fun!