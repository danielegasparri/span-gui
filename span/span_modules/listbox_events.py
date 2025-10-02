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

#Functions to handle the listbox and preview dynamic events

try: #try local import if executed as script
    #GUI import
    from FreeSimpleGUI_local import FreeSimpleGUI as sg
except ModuleNotFoundError: #local import if executed as package
    #GUI import
    from span.FreeSimpleGUI_local import FreeSimpleGUI as sg

import os
import numpy as np
from dataclasses import replace
import matplotlib.ticker as mticker
from scipy.signal import medfilt

# --------------------------------------------------------------------
# Handle selection in the listbox (update preview and params)
# --------------------------------------------------------------------
def handle_list_select(event, values, window, params,
                       _plot_line, ax, fig, preview_interactor, redshift_shifter, stm,
                       sel_full_override=None):
    try:
        if sel_full_override is not None:
            sel_full = sel_full_override
        else:
            selected_base = values['-LIST-'][0] if values['-LIST-'] else ""
            sel_full = ""
            if selected_base:
                if hasattr(window, "metadata") and window.metadata:
                    sel_full = window.metadata.get(selected_base, "")
                else:
                    sel_full = next(
                        (s for s in params.spec_names
                         if isinstance(s, str) and os.path.basename(s) == str(selected_base)),
                        ""
                    )
    except Exception:
        sel_full = ""

    if sel_full and getattr(params, 'spectra_number', 0) > 0:
        wl, fl, *_ = stm.read_spec(sel_full, params.lambda_units)
        wl = np.asarray(wl, dtype=float)
        fl = np.asarray(fl, dtype=float)

        wl_log = np.log10(wl)
        _plot_line.set_data(wl_log, fl)
        ax._last_xydata = (wl, fl)
        ax._x_increasing = bool(wl_log[0] <= wl_log[-1])
        ax.set_xlim(np.nanmin(wl_log), np.nanmax(wl_log))
        ax.set_ylim(np.nanmin(fl), np.nanmax(fl))
        ax.margins(x=0.02, y=0.05)
        ax.relim(); ax.autoscale_view()
        preview_interactor.update_home()
        fig.canvas.draw_idle()

        redshift_shifter._xdata_orig = wl_log.copy()
        redshift_shifter._ydata_orig = fl.copy()
        redshift_shifter.refresh_labels()
        redshift_shifter._cumulative_dx = 0.0
        if redshift_shifter.hud_text:
            redshift_shifter.hud_text.set_text("Estimated z reset")

        def log_to_lin(x, pos):
            return f"{10**x:.0f}"
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(log_to_lin))

        params = replace(params,
                         prev_spec=sel_full,
                         prev_spec_nopath=os.path.basename(sel_full))

    # --- Clean comparison labels if present (back to single spectrum) ---
    if hasattr(ax, "_compare_labels"):
        for t in ax._compare_labels:
            try:
                t.remove()
            except Exception:
                pass
        ax._compare_labels = []
        fig.canvas.draw_idle()
        
    return params


# --------------------------------------------------------------------
# Handle delete with Canc key
# --------------------------------------------------------------------
def handle_list_delete(event, values, window, params,
                       _plot_line, ax, fig, preview_interactor, redshift_shifter):
    last_state = (list(params.spec_names),
                  dict(window.metadata) if hasattr(window, "metadata") else {})

    selected_items = values['-LIST-']
    if not selected_items:
        return params, last_state

    # Retrieve paths
    if hasattr(window, "metadata") and window.metadata:
        selected_paths = [window.metadata[item] for item in selected_items if item in window.metadata]
    else:
        selected_paths = [s for s in params.spec_names if os.path.basename(s) in selected_items]

    new_spec_names = [s for s in params.spec_names if s not in selected_paths]

    if not new_spec_names:
        params = replace(params, spec_names=[0], spectra_number=0,
                         prev_spec="", prev_spec_nopath="")
        window['-LIST-'].update(values=[])
        sg.popup("The spectra list is now empty.\nPlease load new spectra to continue.")

        # Clear preview
        _plot_line.set_data([], [])
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        fig.canvas.draw_idle()

        redshift_shifter._xdata_orig = np.array([])
        redshift_shifter._ydata_orig = np.array([])
        redshift_shifter._cumulative_dx = 0.0
        if redshift_shifter.hud_text:
            redshift_shifter.hud_text.set_text("")
        ax._last_xydata = None
        preview_interactor.update_home()
    else:
        # Update names with LOW SNR flag
        display_names = []
        for s in new_spec_names:
            base = os.path.basename(s)
            if hasattr(window, "metadata") and window.metadata:
                for dn, path in window.metadata.items():
                    if path == s and "[LOW SNR]" in dn:
                        base += "  [LOW SNR]"
                        break
            display_names.append(base)

        new_spec_names_nopath = [os.path.basename(spectrum) for spectrum in new_spec_names]
        
        params = replace(params, spec_names=new_spec_names,
                         spectra_number=len(new_spec_names),
                         prev_spec="", prev_spec_nopath="", spec_names_nopath = new_spec_names_nopath)
        window['-LIST-'].update(values=display_names)
        window.metadata = {dn: s for dn, s in zip(display_names, new_spec_names)}

    return params, last_state


# --------------------------------------------------------------------
# Handle right-click menu (Move Up, Move Down, Remove)
# --------------------------------------------------------------------
def handle_list_menu(event, values, window, params,
                     _plot_line, ax, fig, preview_interactor, redshift_shifter, last_state):
    last_state = (list(params.spec_names),
                  dict(window.metadata) if hasattr(window, "metadata") else {})

    selected_items = values['-LIST-']
    if not selected_items:
        return params, last_state

    if hasattr(window, "metadata") and window.metadata:
        selected_paths = [window.metadata.get(item, "") for item in selected_items]
        selected_paths = [p for p in selected_paths if p]
    else:
        selected_paths = []

    if not selected_paths:
        return params, last_state

    # --- MOVE UP ---
    if event == '↑ Move Up':
        indices = [params.spec_names.index(p) for p in selected_paths]
        if min(indices) > 0:
            for i in sorted(indices):
                params.spec_names[i-1], params.spec_names[i] = params.spec_names[i], params.spec_names[i-1]

            params.spec_names_nopath = [os.path.basename(spectrum) for spectrum in params.spec_names]
            
    # --- MOVE DOWN ---
    elif event == '↓ Move Down':
        indices = [params.spec_names.index(p) for p in selected_paths]
        if max(indices) < len(params.spec_names) - 1:
            for i in sorted(indices, reverse=True):
                params.spec_names[i+1], params.spec_names[i] = params.spec_names[i], params.spec_names[i+1]
                
            params.spec_names_nopath = [os.path.basename(spectrum) for spectrum in params.spec_names]
             

    # --- REMOVE ---
    elif event == 'Remove':
        new_spec_names = [s for s in params.spec_names if s not in selected_paths]
        if not new_spec_names:
            params = replace(params, spec_names=[0], spectra_number=0,
                             prev_spec="", prev_spec_nopath="")
            window['-LIST-'].update(values=[])
            sg.popup("The spectra list is now empty.\nPlease load new spectra to continue.")

            # Clear preview
            _plot_line.set_data([], [])
            ax.set_xlim(0, 1); ax.set_ylim(0, 1)
            fig.canvas.draw_idle()

            redshift_shifter._xdata_orig = np.array([])
            redshift_shifter._ydata_orig = np.array([])
            redshift_shifter._cumulative_dx = 0.0
            if redshift_shifter.hud_text:
                redshift_shifter.hud_text.set_text("")
            ax._last_xydata = None
            preview_interactor.update_home()
            return params, last_state
        else:
            new_spec_names_nopath = [os.path.basename(spectrum) for spectrum in new_spec_names]
            params = replace(params, spec_names=new_spec_names,
                             spectra_number=len(new_spec_names),
                             prev_spec="", prev_spec_nopath="", spec_names_nopath = new_spec_names_nopath)

    # Update listbox names
    display_names = []
    for s in params.spec_names:
        base = os.path.basename(s)
        if hasattr(window, "metadata") and window.metadata:
            for dn, path in window.metadata.items():
                if path == s and "[LOW SNR]" in dn:
                    base += "  [LOW SNR]"
                    break
        display_names.append(base)

    window['-LIST-'].update(values=display_names)
    window.metadata = {dn: s for dn, s in zip(display_names, params.spec_names)}

    # Re-select moved items
    if event in ('↑ Move Up', '↓ Move Down'):
        new_indices = [params.spec_names.index(p) for p in selected_paths]
        window['-LIST-'].update(set_to_index=new_indices)
        listbox_widget = window['-LIST-'].Widget
        for ni in new_indices:
            listbox_widget.activate(ni)
            listbox_widget.see(ni)

    return params, last_state


# --------------------------------------------------------------------
# Handle Undo
# --------------------------------------------------------------------
def handle_undo(event, values, window, params, last_state):
    if last_state is None:
        return params

    old_spec_names, old_metadata = last_state
    old_spec_names_nopath = [os.path.basename(spectrum) for spectrum in old_spec_names]
    
    params = replace(params, spec_names=list(old_spec_names),
                     spectra_number=len(old_spec_names),
                     prev_spec="", prev_spec_nopath="", spec_names_nopath = old_spec_names_nopath)

    display_names = list(old_metadata.keys())
    window['-LIST-'].update(values=display_names)
    window.metadata = dict(old_metadata)
    print("Last action undone.")
    return params


# --------------------------------------------------------------------
# Handle save list
# --------------------------------------------------------------------
def handle_save_list(event, values, window, params):
    try:
        if not params.spec_names or params.spec_names[0] == 0:
            sg.popup("No spectra loaded. Nothing to save.")
        else:
            save_path = sg.popup_get_file(
                "Save spectra list as...",
                save_as=True,
                no_window=True,
                default_extension=".txt",
                file_types=(("Text Files", "*.txt"), ("All Files", "*.*"))
            )
            if save_path:
                with open(save_path, "w") as f:
                    f.write("#Spectrum\n")
                    for spec in params.spec_names:
                        if isinstance(spec, str):
                            f.write(spec + "\n")
                sg.popup(f"Spectra list saved to:\n{save_path}")
    except Exception as e:
        sg.popup(f"Error while saving spectra list:\n You have loaded just one spectrum or not a valid spectra list")



# --------------------------------------------------------------------
# Handle compare spectra (two selected spectra shown together)
# --------------------------------------------------------------------
# --------------------------------------------------------------------
# Handle compare spectra (two selected spectra shown together)
# --------------------------------------------------------------------
def handle_compare_spectra(event, values, window, params,
                           _plot_line, _plot_line2, ax, fig,
                           preview_interactor, redshift_shifter, stm):
    selected_items = values['-LIST-']
    if len(selected_items) != 2:
        sg.popup("Please select exactly two spectra to compare.")
        return params

    # Get paths
    if hasattr(window, "metadata") and window.metadata:
        selected_paths = [window.metadata.get(item, "") for item in selected_items]
        selected_paths = [p for p in selected_paths if p]
    else:
        selected_paths = []

    if len(selected_paths) != 2:
        sg.popup("Could not retrieve both spectra paths for comparison.")
        return params

    # --- First spectrum (blue, main one) ---
    wl, fl, *_ = stm.read_spec(selected_paths[0], params.lambda_units)
    wl = np.asarray(wl, dtype=float)
    fl = np.asarray(fl, dtype=float)
    wl_log = np.log10(wl)

    _plot_line.set_data(wl_log, fl)
    ax._last_xydata = (wl, fl)  # HUD legato al primo
    ax._x_increasing = bool(wl_log[0] <= wl_log[-1])
    ax.set_xlim(np.nanmin(wl_log), np.nanmax(wl_log))
    ax.set_ylim(np.nanmin(fl), np.nanmax(fl))
    ax.margins(x=0.02, y=0.05)
    ax.relim(); ax.autoscale_view()
    preview_interactor.update_home()

    redshift_shifter._xdata_orig = wl_log.copy()
    redshift_shifter._ydata_orig = fl.copy()
    redshift_shifter.refresh_labels()
    redshift_shifter._cumulative_dx = 0.0
    if redshift_shifter.hud_text:
        redshift_shifter.hud_text.set_text("Estimated z reset")

    # --- Second spectrum (red, static reference) ---
    wl2, fl2, *_ = stm.read_spec(selected_paths[1], params.lambda_units)
    wl2 = np.asarray(wl2, dtype=float)
    fl2 = np.asarray(fl2, dtype=float)
    wl_log2 = np.log10(wl2)
    _plot_line2.set_data(wl_log2, fl2)

    # --- Manage labels ---
    # Remove old labels if present
    if hasattr(ax, "_compare_labels"):
        for t in ax._compare_labels:
            try:
                t.remove()
            except Exception:
                pass

    # Add new labels (file basename)
    base1 = os.path.basename(selected_paths[0])
    base2 = os.path.basename(selected_paths[1])
    ax._compare_labels = [
        ax.text(0.01, 0.97, f"{base1}", transform=ax.transAxes,
                ha="left", va="top", color="blue", fontsize=8, weight="bold"),
        ax.text(0.01, 0.92, f"{base2}", transform=ax.transAxes,
                ha="left", va="top", color="red", fontsize=8, weight="bold")
    ]

    # Redraw
    fig.canvas.draw_idle()

    # Axis ticks formatter
    def log_to_lin(x, pos):
        return f"{10**x:.0f}"
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(log_to_lin))

    # Update params only for first spectrum
    params = replace(params,
                     prev_spec=selected_paths[0],
                     prev_spec_nopath=os.path.basename(selected_paths[0]))
    return params


def handle_list_doubleclick(values, window, params, stm):
    selected_items = values['-LIST-']
    if not selected_items:
        return

    selected_base = selected_items[0]
    if hasattr(window, "metadata") and window.metadata:
        sel_full = window.metadata.get(selected_base, "")
    else:
        sel_full = next(
            (s for s in params.spec_names
             if isinstance(s, str) and os.path.basename(s) == str(selected_base)),
            ""
        )

    if not sel_full:
        return

    try:
        wl, fl, *_ = stm.read_spec(sel_full, params.lambda_units)
        wl = np.asarray(wl, dtype=float)
        fl = np.asarray(fl, dtype=float)

        wl_min, wl_max = np.nanmin(wl), np.nanmax(wl)
        n_points = len(wl)
        step = (wl_max - wl_min) / (n_points - 1) if n_points > 1 else np.nan

        flux_mean = np.nanmean(fl)
        flux_median = np.nanmedian(fl)
        flux_std = np.nanstd(fl)

        continuum = medfilt(fl, kernel_size=101)
        residual = fl - continuum

        noise_std = np.nanstd(residual)
        signal = np.nanmedian(fl)

        snr_global_pix = signal / noise_std if noise_std > 0 else np.nan
        snr_global_ang = snr_global_pix / step if step > 0 else np.nan

        flag = ""
        if "[LOW SNR]" in selected_base:
            flag = "FLAG: LOW SNR\n"

        info = f"""Spectrum info

File: {os.path.basename(sel_full)}

Wavelength range: {wl_min:.1f} - {wl_max:.1f} Å
Number of pixels: {n_points}
Step: {step:.2f} Å

Flux mean: {flux_mean:.3g}
Flux median: {flux_median:.3g}
Flux std: {flux_std:.3g}

Global S/N (per pixel): {snr_global_pix:.2f}
Global S/N (per Å): {snr_global_ang:.2f}
{flag}"""

        # set theme and popup
        sg.theme('LightBlue1')
        sg.popup_scrolled(info,
                          title="Spectrum quick info",
                          size=(40, 16),
                          font=("Helvetica", 12))

    except Exception as e:
        sg.popup(f"Error reading spectrum:\n{e}")
