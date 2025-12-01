SPAN: SPectral ANalysis Software V7.3
Daniele Gasparri, December 2025

# Known issues

SPAN tasks work on all supported operating systems. However, some GUI quirks can appear depending on how Tkinter, the display server (Wayland/X11), and OS scaling interact. The following are minor visual issues currently known that I was not able to solve.

- **Reset zoom on Windows:** After zooming in or out, the reset zoom (View → Reset Zoom) might not reset to the default zoom level all the sub-windows. 

- **Screen scaling on Linux:** On some Linux distributions (especially Ubuntu on Wayland), when desktop scaling > 100% the main GUI elements can lose their vertical alignment. This is cosmetic only: functionality is unaffected. The behavior stems from how Tk handles HiDPI/fractional scaling on certain Linux setups. **Workaround:** if you want a perfect symmetry in the main GUI, work with the screen scaling set to 100% (this is a general advice, since some Linux distributions do not work well with fractional scaling). If the GUI panel appears small, use the zooming option of SPAN: View → zoom in. If you find a reliable fix for your environment (e.g., specific Tk/OS settings), please let me know: I’m happy to test and incorporate improvements.

- **Zooming on macOS:** Zooming functions on macOS systems do not resize the buttons and their text. I could not test extensively SPAN on macOS systems to fix this issue. **Workaround:** if possible, avoid the zooming. 

- **Output window on macOS:** The embedded terminal window of SPAN on macOS systems has been disabled since it slowed-down all the tasks by a factor 4-5. 

- **Other GUI aspect differences on... guess where... yes, macOS**: Depending on where you installed Python, from the official Python.org distribution or via Homebrew, something magical happens: the GUI appearance changes. More precisely, with the Homebrew build of Python, buttons, text, and listboxes appear oversized, while the preview shrinks as if by design. This is caused by differences in the bundled Tk version. For the best (and most consistent) appearance, install Python from the official Python.org distribution. 

- **Real time protection in Windows:** Well, this is a known and general behavior of Windows systems that slows down anything. SPAN is no exception. As an example, when real time Windows protection is active, the initial loading and checking of the spectra are 3 times slower. If you work with large datasets, consider temporarily deactivating the real time protection. By doing this, the performances increase and will be comparable to Linux systems. 
