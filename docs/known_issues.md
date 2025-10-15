SPAN: SPectral ANalysis Software V7.1
Daniele Gasparri, October 2025

# Known issues

SPAN’s features work on all supported operating systems. However, some GUI quirks can appear depending on how Tkinter, the display server (Wayland/X11), and OS scaling interact. The following are minor visual issues currently known that I was not able to solve.

- **Reset zoom on Windows:** After zooming in/out, if you open and close one or more sub-windows and then choose View → Reset Zoom, the main panel may not fully return to its initial geometry. **Workaround:** Set your preferred zoom level at the start of a session and keep it unchanged.

- **Screen scaling on Linux:** On some Linux distributions (especially Ubuntu on Wayland), when desktop scaling > 100% the main GUI elements can lose perfect alignment. This is cosmetic only; functionality is unaffected. The behavior stems from how Tk handles HiDPI/fractional scaling on certain Linux setups. **Workaround:** work always with the screen scaling set to 100% (this is a general advice, since some Linux distributions do not work well with fractional scaling). If the GUI panel appears small, use the zooming option of SPAN: View → zoom in. If you find a reliable fix for your environment (e.g., specific Tk/OS settings), please let me know: I’m happy to test and incorporate improvements.

- **Output window on macOS:** The embedded terminal window of SPAN on macOS systems has been disables since it slowed-down all the tasks by a factor 4-5. 
