SPAN: SPectral ANalysis software V6.6
Daniele Gasparri, June 2025

The Spectra Manipulation panel provides several essential operations for processing spectra, organized into three sections:
- Spectra Pre-Processing
- Spectra Processing
- Spectra Math


Activating and Reordering Tasks
To apply a task, enable the corresponding checkbox.
Tasks are executed sequentially, following the order in which they appear, starting from the leftmost frame.
If you need to change the execution order, first activate the desired tasks, then click "Reorder tasks" to arrange them in the window that opens.
Important: If you activate or deactivate tasks after reordering, you must reorder again, or tasks will be executed in their default panel order.

Multiple tasks can be selected simultaneously (e.g., rebinning, Doppler correction, adding noise) without restrictions.
The "Preview Spec." button in the main panel allows you to visualise the spectrum after you activated the task(s) and CLOSED the Spectra Manipulation panel.
The preview displays the final processed spectrum, combining all activated tasks.
Intermediate previews are not available. If you apply three tasks, the preview reflects their cumulative effect.
If no task is selected, the original spectrum is shown.


Mathematical Operations in "Spectra Math"
The four operations in the "Spectra Math" frame:
- Average All
- Normalize and Average All
- Sum All
- Normalize and Sum All
operate on all original spectra, regardless of any individual task selected.
If "Use for Spec. Analysis" is enabled, the program will prioritise the combined result for spectral analysis, overriding any single-spectrum tasks.
Use this option carefully: A terminal message will confirm when combined spectra are set for analysis.


Impact on Spectral Analysis
Any task performed in these frames modifies the selected spectrum and will affect the results in the "Spectral Analysis" panel if you intend to use it.
If at least one task has been activated in the Spectra Manipulation panel, the button color in the main GUI will change to a very annoying red hue. This will remember you that you have active tasks in the Spectra Manipulation panel that are modifying your spectra.
