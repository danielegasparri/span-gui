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

# --- PreviewInteractor: pan/zoom/readout for the Matplotlib preview ----------------
from dataclasses import dataclass
import time
import numpy as np
from typing import Optional, Callable, Any
from scipy.optimize import curve_fit
import matplotlib.patches as patches
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import matplotlib.transforms as mtransforms
import matplotlib.pyplot as plt
import os
from datetime import datetime
import math
import platform
from scipy.signal import find_peaks
import tkinter as tk

        
try:
    from matplotlib.backend_bases import MouseButton
except Exception:
    MouseButton = None


@dataclass
class PreviewInteractor:
    ax: Any                                # Matplotlib Axes
    status_setter: Callable[[str], None]   # fn(str) -> None; updates a status label
    get_snr: Optional[Callable[[float], Optional[float]]] = None
    zoom_step: float = 0.9
    throttle_ms: int = 25
    hud_text: Any = None
    snr_mode: str = "points"
    snr_halfwin_pts: int = 25
    snr_halfwin_A: float = 20.0
    spec_name: Optional[str] = None  
    results_dir: Optional[str] = None 

    def __post_init__(self):
        self._is_linux = platform.system().lower() == "linux"
        self._home_xlim = self.ax.get_xlim()
        self._home_ylim = self.ax.get_ylim()
        self._is_panning = False
        self._is_selecting_line = False
        self._alt_pressed = False
        self._press_pixel = None
        self._press_lambda = None
        self._sel_rect = None
        self._fit_lines = []  # store fit lines for later removal
        self._last_move_ts = 0.0
        self._last_hud_msg = ""
        fig = self.ax.figure
        
        # === EW selection state ===
        self._e_pressed = False
        self._is_selecting_ew_quick = False
        self._ew_quick_rect = None
        self._ew_press_lambda = None

        # === S/N selection state ===
        self._s_pressed = False
        self._is_selecting_sn = False
        self._sn_rect = None
        self._sn_press_lambda = None
        self._sn_lines = []

        # === Delta-lambda/velocity mode ===
        self._d_pressed = False
        self._d_points = []
        self._d_lines = []
        
        # === Flux mode ===
        self._intflux_press_lambda = None
        self._intflux_rect = None
        self._intflux_lines = []
        
        self._linefinder_lines = []
        
        self._active_mode = None 
        self._mode_label = None 

        self._cids = [
            fig.canvas.mpl_connect('scroll_event', self._on_scroll),
            fig.canvas.mpl_connect('button_press_event', self._on_press),
            fig.canvas.mpl_connect('button_release_event', self._on_release),
            fig.canvas.mpl_connect('motion_notify_event', self._on_move),
            fig.canvas.mpl_connect('key_press_event', self._on_key_press),
        ]

    # === Gaussian model ===
    @staticmethod
    def _gauss(x, amp, mu, sigma, c0, c1):
        """Gaussian plus linear continuum."""
        return amp * np.exp(-0.5 * ((x - mu) / sigma) ** 2) + c0 + c1 * x

    @staticmethod
    def _zoom_around(lims, centre, scale):
        a, b = lims
        return (centre + (a - centre) * scale,
                centre + (b - centre) * scale)


    def _draw_rect(self, rect_attr_name, x0_log, x1_log, color):
        y0, y1 = self.ax.get_ylim()
        rect = getattr(self, rect_attr_name, None)
        if rect is None:
            rect = patches.Rectangle((min(x0_log, x1_log), y0),
                                     abs(x1_log - x0_log), y1 - y0,
                                     color=color, alpha=0.25, lw=0)
            setattr(self, rect_attr_name, rect)
            self.ax.add_patch(rect)
        else:
            rect.set_x(min(x0_log, x1_log))
            rect.set_width(abs(x1_log - x0_log))
            rect.set_y(y0)
            rect.set_height(y1 - y0)
        self.ax.figure.canvas.draw_idle()


    def reset_view(self, clear_fits=True):
        """Reset zoom/pan and remove all fit and EW curves if requested."""
        self.ax.set_xlim(self._home_xlim)
        self.ax.set_ylim(self._home_ylim)
        if clear_fits:
            for line in getattr(self, "_fit_lines", []):
                try:
                    line.remove()
                except Exception:
                    pass
            self._fit_lines.clear()

            for obj in getattr(self, "_ew_lines", []):
                try:
                    obj.remove()
                except Exception:
                    pass
            if hasattr(self, "_ew_lines"):
                self._ew_lines.clear()

            for obj in getattr(self, "_sn_lines", []):
                try: obj.remove()
                except Exception: pass
            if hasattr(self, "_sn_lines"):
                self._sn_lines.clear()

            for obj in getattr(self, "_d_lines", []):
                try: obj.remove()
                except Exception: pass
            if hasattr(self, "_d_lines"):
                self._d_lines.clear()

            for obj in getattr(self, "_intflux_lines", []):
                try: obj.remove()
                except Exception: pass
            if hasattr(self, "_intflux_lines"):
                self._intflux_lines.clear()

            for obj_pair in getattr(self, "_linefinder_lines", []):
                for obj in obj_pair:
                    try: obj.remove()
                    except Exception: pass
            if hasattr(self, "_linefinder_lines"):
                self._linefinder_lines.clear()


        self.ax.figure.canvas.draw_idle()

        if self._mode_label is not None:
            try:
                self._mode_label.remove()
            except Exception:
                pass
            self._mode_label = None


    def update_home(self):
        self._home_xlim = self.ax.get_xlim()
        self._home_ylim = self.ax.get_ylim()

    def _on_key_press(self, event):
        """Toggle interactive modes (fit, EW, S/N, Δλ/Δv, Flux) with transient HUD label."""
        key = event.key.lower() if event.key else ""

        # === Snapshot (P) ===
        if key == 'p':
            self._save_preview_snapshot()
            return
        
        # Toggle logic
        if key == 'f':
            self._active_mode = None if self._active_mode == 'fit' else 'fit'
        elif key == 'e':
            self._active_mode = None if self._active_mode == 'ew' else 'ew'
        elif key == 's':
            self._active_mode = None if self._active_mode == 'sn' else 'sn'
        elif key == 'd':
            self._active_mode = None if self._active_mode == 'dv' else 'dv'
        elif key == 'i':
            self._active_mode = None if self._active_mode == 'intflux' else 'intflux'
        elif key == 'l':
            self._compute_linefinder()
            if self.hud_text:
                self.hud_text.set_text("Line peaks identified (press C or double-click to clear)")
            self.ax.figure.canvas.draw_idle()
            return

        elif key == 'c':
            self.clear_overlays()
            return
        else:
            return

        # Remove old label if exists
        if self._mode_label is not None:
            try:
                self._mode_label.remove()
            except Exception:
                pass
            self._mode_label = None

        # Create label only if a mode is active
        if self._active_mode is not None:
            text = f"Mode: {self._active_mode.upper()}"
            self._mode_label = self.ax.text(
                0.02, 0.98, text,
                transform=self.ax.transAxes,
                fontsize=10, fontweight='bold',
                color='#F2E9DC',
                ha='left', va='top',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='black', alpha=0.5, edgecolor='none'),
                zorder=10,
            )

        # Update main HUD message (bottom-left)
        mode_msg = {
            None: "Preview idle",
            'fit': "Line fit mode active (press F again to exit)",
            'ew': "EW mode active (press E again to exit)",
            'sn': "S/N mode active (press S again to exit)",
            'dv': "Δλ/Δv mode active (press D again to exit)",
            'intflux': "Flux mode active (press I again to exit)",
            'linefinder': "Line finder mode active (press L again to exit)",
            'clear': "Cleared all overlays (press C)"
        }[self._active_mode]

        if self.hud_text:
            self.hud_text.set_text(mode_msg)
        self.ax.figure.canvas.draw_idle()


    def _on_scroll(self, event):
        if event.inaxes != self.ax:
            return
        key = (event.key or "").lower()
        only_x = 'control' in key and 'shift' not in key
        only_y = 'shift' in key and 'control' not in key
        scale = self.zoom_step if event.button == 'up' else (1.0 / self.zoom_step)
        cx, cy = event.xdata, event.ydata
        if not only_y:
            self.ax.set_xlim(self._zoom_around(self.ax.get_xlim(), cx, scale))
        if not only_x:
            self.ax.set_ylim(self._zoom_around(self.ax.get_ylim(), cy, scale))
        self.ax.figure.canvas.draw_idle()

    def _on_press(self, event):
        if event.inaxes != self.ax:
            return
        is_left = (MouseButton and event.button == MouseButton.LEFT) or (event.button == 1)

        # === Double click left → full reset ===
        if getattr(event, 'dblclick', False) and is_left:
            self.reset_view(clear_fits=True)
            self._active_mode = None

            # Clear textual HUD
            if self.hud_text:
                self.hud_text.set_text("Preview idle")

            # Remove transient mode label (if present)
            if getattr(self, "_mode_label", None) is not None:
                try:
                    self._mode_label.remove()
                except Exception:
                    pass
                self._mode_label = None
                self.ax.figure.canvas.draw_idle()

            return


        # === FIT mode ===
        if self._active_mode == 'fit' and is_left and event.xdata is not None:
            self._is_selecting_line = True
            self._press_lambda = 10 ** float(event.xdata)
            if self._sel_rect:
                self._sel_rect.remove()
            y0, y1 = self.ax.get_ylim()
            self._sel_rect = patches.Rectangle(
                (event.xdata, y0), 0, y1 - y0, color='orange', alpha=0.25, lw=0
            )
            self.ax.add_patch(self._sel_rect)
            self.ax.figure.canvas.draw_idle()
            return

        # === EW mode ===
        if self._active_mode == 'ew' and is_left and event.xdata is not None:
            self._is_selecting_ew_quick = True
            self._ew_press_lambda = 10 ** float(event.xdata)
            self._draw_rect('_ew_quick_rect', event.xdata, event.xdata, color='purple')
            return

        # === S/N mode ===
        if self._active_mode == 'sn' and is_left and event.xdata is not None:
            self._is_selecting_sn = True
            self._sn_press_lambda = 10 ** float(event.xdata)
            self._draw_rect('_sn_rect', event.xdata, event.xdata, color='teal')
            return

        # === Δλ/Δv mode ===
        if self._active_mode == 'dv' and is_left and event.xdata is not None:
            lam_click = 10 ** float(event.xdata)
            flux_click = float(event.ydata) if event.ydata is not None else np.nan
            self._d_points.append((lam_click, flux_click))
            if len(self._d_points) == 2:
                (lam1, _), (lam2, _) = self._d_points
                self._compute_delta_lambda(lam1, lam2)
                self._d_points.clear()
            return

        # === Integrated Flux mode (I + left drag) ===
        if self._active_mode == 'intflux' and is_left and event.xdata is not None:
            self._intflux_press_lambda = 10 ** float(event.xdata)
            self._draw_rect('_intflux_rect', event.xdata, event.xdata, color='lime')
            return

        # === Normal panning (only when idle) ===
        if is_left and self._active_mode is None:
            self._is_panning = True
            self._press_pixel = (event.x, event.y)
            self._pan_xlim0 = self.ax.get_xlim()
            self._pan_ylim0 = self.ax.get_ylim()
            try:
                self.ax.set_cursor(1)
            except Exception:
                pass


    def _on_release(self, event):
        if event.inaxes != self.ax or event.xdata is None:
            return

        # === FIT mode ===
        if self._active_mode == 'fit' and self._is_selecting_line:
            lam1 = self._press_lambda
            lam2 = 10 ** float(event.xdata)
            if lam1 is not None and abs(lam2 - lam1) > 1.0:
                self._fit_feature_region(min(lam1, lam2), max(lam1, lam2))
            if self._sel_rect:
                self._sel_rect.remove()
                self._sel_rect = None
            self._is_selecting_line = False
            self.ax.figure.canvas.draw_idle()
            return

        # === EW mode ===
        if self._active_mode == 'ew' and self._is_selecting_ew_quick:
            lam1 = self._ew_press_lambda
            lam2 = 10 ** float(event.xdata)
            self._is_selecting_ew_quick = False
            self._ew_press_lambda = None
            if self._ew_quick_rect:
                self._ew_quick_rect.remove()
                self._ew_quick_rect = None
            self.ax.figure.canvas.draw_idle()
            if lam1 is not None and abs(lam2 - lam1) > 1.0:
                self._compute_ew_quick(min(lam1, lam2), max(lam1, lam2))
            return

        # === S/N mode ===
        if self._active_mode == 'sn' and self._is_selecting_sn:
            lam1 = self._sn_press_lambda
            lam2 = 10 ** float(event.xdata)
            self._is_selecting_sn = False
            self._sn_press_lambda = None
            if self._sn_rect:
                try:
                    self._sn_rect.remove()
                except Exception:
                    pass
                self._sn_rect = None
            self.ax.figure.canvas.draw_idle()
            if lam1 is not None and abs(lam2 - lam1) > 1.0:
                self._compute_sn_region(min(lam1, lam2), max(lam1, lam2))
            return


        # === Integrated Flux mode ===
        if self._active_mode == 'intflux' and event.inaxes == self.ax and event.xdata is not None:
            lam1 = self._intflux_press_lambda
            lam2 = 10 ** float(event.xdata)
            self._intflux_press_lambda = None

            # remove the rectangle
            if getattr(self, "_intflux_rect", None) is not None:
                try:
                    self._intflux_rect.remove()
                except Exception:
                    pass
                self._intflux_rect = None
                self.ax.figure.canvas.draw_idle()

            if lam1 is not None and abs(lam2 - lam1) > 1.0:
                self._compute_integrated_flux(min(lam1, lam2), max(lam1, lam2))
            return

        # === End of panning ===
        if self._is_panning:
            self._is_panning = False
            self._press_pixel = None
            try:
                self.ax.set_cursor(0)
            except Exception:
                pass



    def _fit_feature_region(self, lam1, lam2):
        """Perform a Gaussian fit in the selected region (works for emission and absorption)."""

        data = getattr(self.ax, "_last_xydata", None)
        if not data:
            return
        x, y = np.asarray(data[0]), np.asarray(data[1])
        ok = np.isfinite(x) & np.isfinite(y)
        lam = x[ok]
        flux = y[ok]
        sel = (lam >= lam1) & (lam <= lam2)
        if np.sum(sel) < 5:
            return

        lam_sel, flux_sel = lam[sel], flux[sel]

        # === Initial guess ===
        median_flux = np.median(flux_sel)
        amp_guess = flux_sel[np.argmax(np.abs(flux_sel - median_flux))] - median_flux
        if np.mean(flux_sel) < median_flux:
            amp_guess = -abs(amp_guess)

        mu_guess = lam_sel[np.argmax(np.abs(flux_sel - median_flux))]
        sigma_guess = (lam2 - lam1) / 6
        p0 = [amp_guess, mu_guess, sigma_guess, median_flux, 0]

        try:
            popt, pcov = curve_fit(self._gauss, lam_sel, flux_sel, p0=p0)
            amp, mu, sigma, c0, c1 = popt
            perr = np.sqrt(np.diag(pcov)) if pcov is not None else np.full_like(popt, np.nan)
            amp_err, mu_err, sigma_err, _, _ = perr

            fwhm = 2.355 * abs(sigma)
            fwhm_err = 2.355 * abs(sigma_err)

            eqw = np.trapz(1 - (flux_sel / (c0 + c1 * lam_sel)), lam_sel)
            residuals = flux_sel - self._gauss(lam_sel, *popt)
            eqw_err = np.std(residuals / (c0 + c1 * lam_sel)) * (lam2 - lam1) / np.sqrt(len(lam_sel))

            flux_int = abs(amp) * abs(sigma) * np.sqrt(2 * np.pi)
            flux_err = np.sqrt(
                (np.sqrt(2 * np.pi) * sigma * amp_err) ** 2 +
                (np.sqrt(2 * np.pi) * amp * sigma_err) ** 2
            )

            sn_line = flux_int / flux_err if np.isfinite(flux_err) and flux_err > 0 else np.nan
            line_type = "emission" if amp > 0 else "absorption"

            # Plot fit
            fit_curve, = self.ax.plot(np.log10(lam_sel), self._gauss(lam_sel, *popt),
                                      color='orange', lw=1.5, alpha=0.8)
            self._fit_lines.append(fit_curve)
            self.ax.figure.canvas.draw_idle()

            # === HUD message ===
            msg = (f"λ₀ = {mu:.2f} Å · FWHM = {fwhm:.2f} Å · "
                   f"EW ≈ {eqw:.3f} Å · Flux = {flux_int:.3e}")
            if self.status_setter:
                self.status_setter(msg)
            if self.hud_text:
                self.hud_text.set_text(msg)
                self.ax.figure.canvas.draw_idle()

            # === Terminal print (detailed) ===
            print(f"\n[SPAN-Preview] Line fit results:")
            print(f"  Type: {line_type}")
            print(f"  Range: {lam1:.2f} – {lam2:.2f} Å")
            print(f"  λ₀ (centre): {mu:.4f} ± {mu_err:.4f} Å")
            print(f"  FWHM: {fwhm:.4f} ± {fwhm_err:.4f} Å")
            print(f"  EW: {eqw:.4f} ± {eqw_err:.4f} Å")
            print(f"  Flux: {flux_int:.4e} ± {flux_err:.4e}")
            print(f"  Amplitude: {amp:.4e} ± {amp_err:.4e}")
            print(f"  σ (dispersion): {sigma:.4f} ± {sigma_err:.4f} Å")
            print(f"  S/N (Flux / Flux_err): {sn_line:.2f}\n")

            # === Append to log file ===
            if self.results_dir is not None and os.path.isdir(self.results_dir):
                out_dir = self.results_dir
            else:
                # fallback: crea SPAN_results se non esiste
                out_dir = os.path.join(os.getcwd(), "SPAN_results")
                os.makedirs(out_dir, exist_ok=True)

            # Saving
            log_path = os.path.join(out_dir, "span_preview_fits.log")

            header = ("# Timestamp,Spectrum,Type,lambda0_A,FWHM_A,EW_A,"
                      "Flux,Flux_err,SN,Range_start_A,Range_end_A\n")

            spec_label = self.spec_name if self.spec_name is not None else "Unknown"

            line = (f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')},"
                    f"{spec_label},{line_type},{mu:.4f},{fwhm:.4f},{eqw:.4f},"
                    f"{flux_int:.4e},{flux_err:.4e},{sn_line:.2f},{lam1:.2f},{lam2:.2f}\n")

            if not os.path.exists(log_path):
                with open(log_path, "w") as f:
                    f.write(header)
                    f.write(line)
            else:
                with open(log_path, "a") as f:
                    f.write(line)
                    
        except Exception as e:
            print("Fit failed:", e)
            

    def _compute_ew_quick(self, lam1, lam2):
        """Compute approximate EW from a single selected region."""
        data = getattr(self.ax, "_last_xydata", None)
        if not data:
            return
        x, y = np.asarray(data[0]), np.asarray(data[1])
        ok = np.isfinite(x) & np.isfinite(y)
        lam = x[ok]
        flux = y[ok]
        sel = (lam >= lam1) & (lam <= lam2)
        if np.sum(sel) < 5:
            return

        lam_sel = lam[sel]
        flux_sel = flux[sel]

        # Use 15% of edges as pseudo-sidebands for continuum
        n = len(lam_sel)
        k = max(3, int(0.15 * n))
        lam_blue, flux_blue = lam_sel[:k], flux_sel[:k]
        lam_red, flux_red = lam_sel[-k:], flux_sel[-k:]

        # Linear continuum from edges
        fb, fr = np.median(flux_blue), np.median(flux_red)
        lb, lr = np.median(lam_blue), np.median(lam_red)
        if not np.isfinite([fb, fr, lb, lr]).all() or abs(lr - lb) < 1e-6:
            print("[SPAN-Preview] EW quick: unable to estimate continuum.")
            return
        a = (fr - fb) / (lr - lb)
        b = fb - a * lb

        Fc = a * lam_sel + b
        ew = np.trapz(1 - (flux_sel / Fc), lam_sel)

        # Rough error from local noise
        noise = np.std((flux_sel - Fc) / Fc)
        ew_err = abs(ew) * noise

        # === Plot: continuum line + fill area ===
        cont_line, = self.ax.plot(np.log10(lam_sel), Fc, color='purple', lw=1.5, alpha=0.8, ls='--')
        fill_poly = self.ax.fill_between(np.log10(lam_sel), flux_sel, Fc,
                                         color='purple', alpha=0.25, lw=0)
        # save for removal on reset
        if not hasattr(self, "_ew_lines"):
            self._ew_lines = []
        self._ew_lines.extend([cont_line, fill_poly])
        self.ax.figure.canvas.draw_idle()

        # === HUD + Terminal ===
        msg = f"EW ≈ {ew:.3f} ± {ew_err:.3f} Å (quick)"
        if self.status_setter:
            self.status_setter(msg)
        if self.hud_text:
            self.hud_text.set_text(msg)
            self.ax.figure.canvas.draw_idle()

        print(f"[SPAN-Preview] EW quick: range={lam1:.2f}–{lam2:.2f} Å, EW={ew:.4f} ± {ew_err:.4f} Å")

        # === Log ===
        if self.results_dir is not None and os.path.isdir(self.results_dir):
            out_dir = self.results_dir
        else:
            out_dir = os.path.join(os.getcwd(), "SPAN_results")
            os.makedirs(out_dir, exist_ok=True)

        log_path = os.path.join(out_dir, "span_preview_EW.log")
        header = "# Timestamp,Spectrum,Mode,EW_A,EW_err_A,Range_start_A,Range_end_A\n"
        spec_label = self.spec_name if self.spec_name else "Unknown"
        line = (f"{datetime.now():%Y-%m-%d %H:%M:%S},{spec_label},quick,"
                f"{ew:.4f},{ew_err:.4f},{lam1:.2f},{lam2:.2f}\n")
        if not os.path.exists(log_path):
            with open(log_path, "w") as f:
                f.write(header)
                f.write(line)
        else:
            with open(log_path, "a") as f:
                f.write(line)
        # print(f"[SPAN-Preview] EW logged to: {log_path}")

    def _compute_sn_region(self, lam1, lam2):
        """Measure ⟨S⟩, σ (MAD-based) and S/N"""
        data = getattr(self.ax, "_last_xydata", None)
        if not data:
            return
        x, y = np.asarray(data[0]), np.asarray(data[1])
        ok = np.isfinite(x) & np.isfinite(y)
        lam = x[ok]
        flux = y[ok]

        sel = (lam >= lam1) & (lam <= lam2)
        if np.sum(sel) < 5:
            return

        lam_sel = lam[sel]
        flux_sel = flux[sel]

        S_med = float(np.median(flux_sel))
        # Noise: sigma ≈ 1.4826 * MAD
        mad = float(np.median(np.abs(flux_sel - S_med)))
        sigma = 1.4826 * mad if np.isfinite(mad) else np.nan

        sn_val = (S_med / sigma) if (np.isfinite(sigma) and sigma > 0) else np.nan

        # === Draw ===
        yline = np.full_like(lam_sel, S_med, dtype=float)
        sn_line_obj, = self.ax.plot(np.log10(lam_sel), yline, color='teal', lw=1.2, ls='--', alpha=0.9)
        self._sn_lines.append(sn_line_obj)
        self.ax.figure.canvas.draw_idle()

        # === HUD ===
        msg = f"S/N ≈ {sn_val:.1f}  ·  ⟨Flux⟩={S_med:.3e}  ·  σ={sigma:.3e}  ·  [{lam1:.1f}-{lam2:.1f}] Å"
        if self.status_setter:
            self.status_setter(msg)
        if self.hud_text:
            self.hud_text.set_text(msg)
            self.ax.figure.canvas.draw_idle()

        # === Terminal ===
        print("[SPAN-Preview] S/N region:")
        print(f"  Range: {lam1:.2f} – {lam2:.2f} Å  (N={lam_sel.size})")
        print(f"  <Flux> (median): {S_med:.6e}")
        print(f"  sigma (1.4826*MAD): {sigma:.6e}")
        print(f"  S/N: {sn_val:.2f}\n")

        # === Log  ===
        if self.results_dir is not None and os.path.isdir(self.results_dir):
            out_dir = self.results_dir
        else:
            out_dir = os.path.join(os.getcwd(), "SPAN_results")
            os.makedirs(out_dir, exist_ok=True)

        log_path = os.path.join(out_dir, "span_preview_SNR.log")
        header = "# Timestamp,Spectrum,Flux_median,Noise_sigma,SN,Npts,Range_start_A,Range_end_A\n"
        spec_label = self.spec_name if self.spec_name else "Unknown"
        line = (f"{datetime.now():%Y-%m-%d %H:%M:%S},{spec_label},"
                f"{S_med:.6e},{sigma:.6e},{sn_val:.3f},{lam_sel.size},"
                f"{lam1:.2f},{lam2:.2f}\n")
        if not os.path.exists(log_path):
            with open(log_path, "w") as f:
                f.write(header); f.write(line)
        else:
            with open(log_path, "a") as f:
                f.write(line)
        # print(f"[SPAN-Preview] S/N logged to: {log_path}")


    def _compute_delta_lambda(self, lam1, lam2):
        """Calculate Δλ e Δv between two points."""

        c = 299792.458  # km/s
        lam1, lam2 = float(lam1), float(lam2)
        dlam = lam2 - lam1
        lmean = 0.5 * (lam1 + lam2)
        dv = c * dlam / lmean

        # === Draw ===
        xline = np.log10([lam1, lam2])
        y0, y1 = self.ax.get_ylim()
        line, = self.ax.plot(xline, [y0 + 0.1*(y1 - y0)]*2, color='deepskyblue', lw=1.8, alpha=0.8)
        self._d_lines.append(line)
        self.ax.figure.canvas.draw_idle()

        # === HUD message ===
        msg = f"Δλ = {dlam:.2f} Å  ·  Δv = {dv:.1f} km/s"
        if self.status_setter:
            self.status_setter(msg)
        if self.hud_text:
            self.hud_text.set_text(msg)
            self.ax.figure.canvas.draw_idle()

        # === Terminal ===
        print(f"[SPAN-Preview] Δλ/Δv measurement:")
        print(f"  λ1 = {lam1:.2f} Å,  λ2 = {lam2:.2f} Å")
        print(f"  Δλ = {dlam:.3f} Å")
        print(f"  Δv = {dv:.2f} km/s\n")

        # === Log ===

        if self.results_dir is not None and os.path.isdir(self.results_dir):
            out_dir = self.results_dir
        else:
            out_dir = os.path.join(os.getcwd(), "SPAN_results")
            os.makedirs(out_dir, exist_ok=True)

        log_path = os.path.join(out_dir, "span_preview_deltav.log")
        header = "# Timestamp,Spectrum,lambda1_A,lambda2_A,delta_lambda_A,delta_v_kms\n"
        spec_label = self.spec_name if self.spec_name else "Unknown"
        line = (f"{datetime.now():%Y-%m-%d %H:%M:%S},{spec_label},"
                f"{lam1:.3f},{lam2:.3f},{dlam:.3f},{dv:.2f}\n")
        if not os.path.exists(log_path):
            with open(log_path, "w") as f:
                f.write(header); f.write(line)
        else:
            with open(log_path, "a") as f:
                f.write(line)
        # print(f"[SPAN-Preview] Δλ/Δv logged to: {log_path}")


    def _compute_integrated_flux(self, lam1, lam2):
        """Compute integrated flux and approximate EW in a user-selected region."""

        data = getattr(self.ax, "_last_xydata", None)
        if not data:
            return

        x, y = np.asarray(data[0]), np.asarray(data[1])
        ok = np.isfinite(x) & np.isfinite(y)
        lam, flux = x[ok], y[ok]
        sel = (lam >= lam1) & (lam <= lam2)
        if np.sum(sel) < 5:
            return

        lam_sel, flux_sel = lam[sel], flux[sel]

        # === Continuum line between edges ===
        cont_slope = (flux_sel[-1] - flux_sel[0]) / (lam_sel[-1] - lam_sel[0])
        cont_intercept = flux_sel[0] - cont_slope * lam_sel[0]
        cont = cont_slope * lam_sel + cont_intercept

        # === Integrated flux and EW ===
        flux_int = np.trapz(flux_sel - cont, lam_sel)
        ew_num = np.trapz(1 - (flux_sel / cont), lam_sel)

        # === Visual overlay ===
        band = self.ax.fill_between(np.log10(lam_sel), cont, flux_sel,
                                    color='lime', alpha=0.3)
        self._intflux_lines.append(band)
        self.ax.figure.canvas.draw_idle()

        # === HUD message ===
        msg = f"Integrated flux = {flux_int:.3e}  ·  EW ≈ {ew_num:.2f} Å"
        if self.status_setter:
            self.status_setter(msg)
        if self.hud_text:
            self.hud_text.set_text(msg)
            self.ax.figure.canvas.draw_idle()

        # === Terminal print ===
        print(f"[SPAN-Preview] Integrated flux measurement:")
        print(f"  Range: {lam1:.2f} – {lam2:.2f} Å")
        print(f"  Integrated flux: {flux_int:.4e}")
        print(f"  EW (numerical): {ew_num:.3f} Å\n")

        # === Log to file ===
        if self.results_dir is not None and os.path.isdir(self.results_dir):
            out_dir = self.results_dir
        else:
            out_dir = os.path.join(os.getcwd(), "SPAN_results")
            os.makedirs(out_dir, exist_ok=True)

        log_path = os.path.join(out_dir, "span_preview_intflux.log")
        header = "# Timestamp,Spectrum,lambda1_A,lambda2_A,Integrated_flux,EW_num_A\n"
        spec_label = self.spec_name if self.spec_name else "Unknown"
        line = (f"{datetime.now():%Y-%m-%d %H:%M:%S},{spec_label},"
                f"{lam1:.3f},{lam2:.3f},{flux_int:.4e},{ew_num:.3f}\n")

        if not os.path.exists(log_path):
            with open(log_path, "w") as f:
                f.write(header)
                f.write(line)
        else:
            with open(log_path, "a") as f:
                f.write(line)


    def _save_preview_snapshot(self):
        """Save only the visible spectrum plot region (with small margins)."""

        # Determine output directory
        if self.results_dir and os.path.isdir(self.results_dir):
            out_dir = os.path.join(self.results_dir, "preview_snapshots")
        else:
            out_dir = os.path.join(os.getcwd(), "SPAN_results", "preview_snapshots")
        os.makedirs(out_dir, exist_ok=True)

        spec_label = self.spec_name if getattr(self, "spec_name", None) else "snapshot"
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        path = os.path.join(out_dir, f"snapshot_{spec_label}_{timestamp}.png")

        # === Prepare canvas ===
        self.ax.figure.canvas.draw()
        renderer = self.ax.figure.canvas.get_renderer()

        # Compute tight bounding box for the Axes
        bbox = self.ax.get_window_extent(renderer=renderer)

        # Add a small margin (in pixels)
        margin = 70  # adjust if needed
        expanded_bbox = mtransforms.Bbox.from_extents(
            bbox.x0 - margin, bbox.y0 - margin, bbox.x1 + margin, bbox.y1 + margin
        )

        # Convert to inches
        bbox_inches = expanded_bbox.transformed(self.ax.figure.dpi_scale_trans.inverted())

        # Save only that cropped region
        self.ax.figure.savefig(path, dpi=300, bbox_inches=bbox_inches, facecolor="white")

        # === Feedback ===
        msg = f"Snapshot saved → {os.path.basename(path)}"
        print(f"[SPAN-Preview] {msg}")
        if self.hud_text:
            self.hud_text.set_text(msg)
        if self.status_setter:
            self.status_setter(msg)
        self.ax.figure.canvas.draw_idle()



    def _compute_linefinder(self):
        """Detect and mark spectral lines automatically using simple peak detection."""

        data = getattr(self.ax, "_last_xydata", None)
        if not data:
            return

        lam, flux = np.asarray(data[0]), np.asarray(data[1])
        ok = np.isfinite(lam) & np.isfinite(flux)
        lam, flux = lam[ok], flux[ok]

        # === Peak detection ===
        median = np.nanmedian(flux)
        std = np.nanstd(flux)
        threshold = median + 2 * std
        peaks, props = find_peaks(flux, height=threshold, distance=5)

        self._linefinder_peaks = lam[peaks]
        
        y0, y1 = self.ax.get_ylim()

        for p in peaks:
            lam_p = lam[p]
            line = self.ax.axvline(np.log10(lam_p), color="lime", ls="--", lw=1.2, alpha=0.7)
            txt = self.ax.text(
                np.log10(lam_p),
                y1 - 0.05 * (y1 - y0),
                f"{lam_p:.0f}",
                color="lime",
                fontsize=8,
                ha="center",
                va="top",
            )
            self._linefinder_lines.append((line, txt))

        self.ax.figure.canvas.draw_idle()

        msg = f"Found {len(peaks)} lines"
        if self.hud_text:
            self.hud_text.set_text(msg)
        if self.status_setter:
            self.status_setter(msg)

        print(f"[SPAN-Preview] Line finder: {len(peaks)} lines detected.")



    def update_linefinder_positions(self, delta_loglam):
        """Shift existing line-finder markers to follow spectrum shifts."""
        if not hasattr(self, "_linefinder_lines") or not self._linefinder_lines:
            return
        if not hasattr(self, "_linefinder_peaks") or not np.any(self._linefinder_peaks):
            return

        for (line, txt), lam0 in zip(self._linefinder_lines, self._linefinder_peaks):
            new_x = np.log10(lam0) + delta_loglam
            line.set_xdata([new_x, new_x])
            txt.set_x(new_x)

        self.ax.figure.canvas.draw_idle()


    def clear_overlays(self):
        """Remove all graphical overlays (fit, EW, S/N, Δλ/Δv, flux, linefinder) but keep zoom/pan."""
        groups = [
            "_fit_lines", "_ew_lines", "_sn_lines", "_d_lines",
            "_intflux_lines", "_linefinder_lines"
        ]

        for group in groups:
            if hasattr(self, group):
                objs = getattr(self, group)
                for obj in objs:
                    # Gestisce sia linee che fill_between (PolyCollections)
                    try:
                        if hasattr(obj, "remove"):
                            obj.remove()
                        elif isinstance(obj, tuple):
                            for o in obj:
                                try: o.remove()
                                except Exception: pass
                    except Exception:
                        pass
                setattr(self, group, [])

        # Aggiorna canvas
        self.ax.figure.canvas.draw_idle()

        # Aggiorna HUD
        if self.hud_text:
            self.hud_text.set_text("All overlays cleared (zoom/pan preserved)")


    def _on_move(self, event):
        now = time.time()
        if (now - self._last_move_ts) * 1000.0 >= self.throttle_ms:
            self._last_move_ts = now

            # === Update selection rectangle ===
            if self._is_selecting_line and self._sel_rect is not None and event.xdata is not None:
                x0 = np.log10(self._press_lambda)
                x1 = event.xdata
                self._sel_rect.set_x(min(x0, x1))
                self._sel_rect.set_width(abs(x1 - x0))
                self.ax.figure.canvas.draw_idle()
                return

            # === Update EW quick rectangle ===
            if self._is_selecting_ew_quick and event.inaxes == self.ax and event.xdata is not None and self._ew_quick_rect is not None:
                x0 = np.log10(self._ew_press_lambda)
                x1 = event.xdata
                self._draw_rect('_ew_quick_rect', x0, x1, color='purple')
                return

            # === Update S/N selection rectangle ===
            if self._is_selecting_sn and event.inaxes == self.ax and event.xdata is not None and self._sn_rect is not None:
                x0 = np.log10(self._sn_press_lambda)
                x1 = event.xdata
                self._draw_rect('_sn_rect', x0, x1, color='teal')
                return

            # === Update Integrated Flux selection rectangle ===
            if self._active_mode == 'intflux' and self._intflux_rect is not None and event.xdata is not None:
                x0 = np.log10(self._intflux_press_lambda)
                x1 = event.xdata
                self._intflux_rect.set_x(min(x0, x1))
                self._intflux_rect.set_width(abs(x1 - x0))
                self.ax.figure.canvas.draw_idle()
                return

            # === HUD update ===
            if not self._is_panning and not self._is_selecting_line and event.inaxes == self.ax and event.xdata is not None:
                lam_log = float(event.xdata)
                lam = 10 ** lam_log
                data = getattr(self.ax, "_last_xydata", None)
                flux = float(event.ydata) if event.ydata is not None else float("nan")
                if data:
                    x, y = data
                    x = np.asarray(x, dtype=float)
                    y = np.asarray(y, dtype=float)
                    ok = np.isfinite(x) & np.isfinite(y)
                    if np.any(ok):
                        x = x[ok]
                        y = y[ok]
                        if x.size >= 2:
                            if x[0] <= x[-1]:
                                flux = float(np.interp(lam, x, y))
                            else:
                                flux = float(np.interp(lam[::-1], y[::-1]))

                snr_txt = " · SNR n/a"
                if self.get_snr is not None:
                    try:
                        snr_val = self.get_snr(lam)
                        if snr_val is not None and np.isfinite(snr_val):
                            snr_txt = f" · SNR≈{snr_val:.1f}"
                    except Exception:
                        snr_txt = " · SNR n/a"

                msg = f"λ = {lam:.2f} Å · Flux = {flux:.4g}{snr_txt} (50 pts)"
                if self.hud_text is not None and msg != self._last_hud_msg:
                    self.hud_text.set_text(msg)
                    self._last_hud_msg = msg
                    self.ax.figure.canvas.draw_idle()
                if self.status_setter is not None:
                    self.status_setter(msg)

            #=== Panning ===
            if self._is_panning and event.inaxes == self.ax and self._press_pixel is not None:
                dx_pix = event.x - self._press_pixel[0]
                dy_pix = event.y - self._press_pixel[1]
                inv = self.ax.transData.inverted()
                x0_data, y0_data = inv.transform(self._press_pixel)
                x1_data, y1_data = inv.transform((event.x, event.y))
                dx = x1_data - x0_data
                dy = y1_data - y0_data
                self.ax.set_xlim(self._pan_xlim0[0] - dx, self._pan_xlim0[1] - dx)
                self.ax.set_ylim(self._pan_ylim0[0] - dy, self._pan_ylim0[1] - dy)
                self.ax.figure.canvas.draw_idle()



# Class for estimating the redshift by manual shifting the spectrum on thr Preview window
class SpectrumShifterInteractor:
    """
    Allow shifting the spectrum horizontally (Right Click + Drag) 
    to align with fixed rest-frame lines and estimate redshift.
    Right double-click resets the spectrum to its original position.
    """

    def __init__(self, ax, line, hud_text=None, parent=None):
        self.ax = ax
        self.line = line  # the Line2D object of the spectrum
        self.hud_text = hud_text
        self.parent = parent 
        self._press_event = None
        self._xdata0 = None
        self._last_valid_x = None
        
        # Save original data
        self._xdata_orig = line.get_xdata().copy()
        self._ydata_orig = line.get_ydata().copy()
        self._cumulative_dx = 0.0
        self._labels = []
        
        self._xdata_orig = np.log10(line.get_xdata().copy())
        self._ydata_orig = line.get_ydata().copy()
        line.set_xdata(self._xdata_orig)

        # Rest-frame reference lines (estese per high-z)
        self.ref_lines = {
            "[O II]": 3727.0,
            "Hβ": 4861.0,
            "[O III]": 5007.0,
            "Hα": 6563.0,
            "Ca II 8498": 8498.0,
            "Ca II 8542": 8542.0,
            "Ca II 8662": 8662.0,
        }

        self._draw_markers()

        # Connect events
        fig = ax.figure
        fig.canvas.mpl_connect("button_press_event", self.on_press)
        fig.canvas.mpl_connect("button_release_event", self.on_release)
        fig.canvas.mpl_connect("motion_notify_event", self.on_motion)

    def _draw_markers(self):
        # Clear old labels
        for t in self._labels:
            try: 
                t.remove()
            except Exception: 
                pass
        self._labels = []

        for name, lam in self.ref_lines.items():
            lam_log = np.log10(lam)
            self.ax.axvline(lam_log, color="orange", ls="--", lw=0.8, alpha=0.7, zorder=5)
            t = self.ax.text(
                lam_log, self.ax.get_ylim()[1]*0.95, f"{name} ({lam:.0f} Å)",
                rotation=90, va="top", ha="center",
                fontsize=8, color="darkred", zorder=6
            )
            self._labels.append(t)


    def on_press(self, event):
        if event.inaxes != self.ax:
            return
        if event.button == 3:
            if getattr(event, "dblclick", False):
                self.line.set_data(self._xdata_orig.copy(), self._ydata_orig.copy())
                self._cumulative_dx = 0.0
                self.ax._last_xydata = (10**self._xdata_orig, self._ydata_orig)
                if self.hud_text:
                    self.hud_text.set_text("Estimated z reset")
                self.ax.figure.canvas.draw_idle()
                            
                if hasattr(self, "parent") and hasattr(self.parent, "update_linefinder_positions"):
                    self.parent.update_linefinder_positions(0.0)
                return
            # Start drag
            self._press_event = event.xdata
            self._xdata0 = self.line.get_xdata().copy()
            self._last_valid_x = event.xdata

    def on_motion(self, event):
        if self._press_event is None or event.inaxes != self.ax:
            return
        if event.xdata is None:
            return  # ignore if out of plot margins

        dx = event.xdata - self._press_event
        self.line.set_xdata(self._xdata0 + dx)
        self.ax._last_xydata = (10**self.line.get_xdata(), self.line.get_ydata())

        self._cumulative_dx += dx
        self._press_event = event.xdata 
        self._xdata0 = self.line.get_xdata().copy()
        
        delta_log = self._cumulative_dx 
        z = (10**(-delta_log) - 1)


        if self.hud_text:
            self.hud_text.set_text(f"Estimated z ≈ {z:.3f}")
        else:
            self.ax.set_title(f"Estimated z ≈ {z:.3f}", fontsize=10)

        self.ax.figure.canvas.draw_idle()
        
        # === Notify parent (if linefinder lines exist) ===
        if self.parent is not None and hasattr(self.parent, "update_linefinder_positions"):
            try:
                self.parent.update_linefinder_positions(self._cumulative_dx)
            except Exception:
                pass
            
        self._last_valid_x = event.xdata

    def on_release(self, event):
        if event.button == 3 and self._press_event is not None:
            self._press_event = None
            self._xdata0 = None

    def refresh_labels(self):
        self._draw_markers()
        self.ax.figure.canvas.draw_idle()


# --------------------------------------------------------------
# Create preview figure and canvas
# --------------------------------------------------------------

def create_preview(layout, window, preview_key='-CANVAS-'):
    """
    Matplotlib-in-Tk Canvas with:
      - proper embedding via create_window (no clipping, correct origin),
      - responsive sizing on <Configure>,
      - dynamic typography scaling (labels, ticks, HUD, line width).
    """
    if layout == 'linux' or layout =='windows': # Linux and Windows systems deserve a special treatment for handling the scaling factor of the screen (if applied)
        tk_canvas = window[preview_key].TKCanvas  # real tkinter.Canvas
        try:
            tk_canvas.configure(highlightthickness=0, bd=0, bg=window.TKroot['bg'])
        except Exception:
            try:
                tk_canvas.configure(highlightthickness=0, bd=0)
            except Exception:
                pass

        # --- Base visual parameters (tuned for your plot) ---
        dpi = 100
        BASE_W_PX, BASE_H_PX = 900, 420     # reference pixels used for scaling
        BASE_LABEL = 12                     # axis label font (pt) at reference size
        BASE_TICK  = 12                      # tick labels font (pt)
        BASE_HUD   = 11                     # HUD font (pt)
        BASE_LW    = 0.9                    # plotted line width (points)

        # --- Figure & axes ---
        fig = Figure(figsize=(6.0, 3.0), dpi=dpi, constrained_layout=False)
        ax = fig.add_subplot(111)
        fig.subplots_adjust(left=0.11, right=0.98, top=0.93, bottom=0.16)

        (plot_line,) = ax.plot([], [], lw=BASE_LW)
        ax.set_xlabel("Wavelength [Å]")
        ax.set_ylabel("Flux")
        hud_text = ax.text(
            0.995, 0.01, "", transform=ax.transAxes,
            ha='right', va='bottom', fontsize=BASE_HUD,
            color='black', zorder=50, clip_on=False
        )

        # --- Embed MPL widget into the Tk Canvas via create_window ---
        mplt_canvas = FigureCanvasTkAgg(fig, master=tk_canvas)
        mplt_widget = mplt_canvas.get_tk_widget()
        window_id = tk_canvas.create_window(0, 0, window=mplt_widget, anchor="nw")

        def _apply_style_scale(scale: float):
            """Scale fonts and strokes according to the current canvas size."""
            # clamp to a sensible range to avoid extremes on tiny/huge windows
            s = max(0.7, min(scale, 2.0))

            ax.xaxis.label.set_size(BASE_LABEL * s)
            ax.yaxis.label.set_size(BASE_LABEL * s)
            ax.tick_params(axis='both', labelsize=BASE_TICK * s, width=0.75 * s, length=3.5 * s)

            # scale spines and plotted line
            for spine in ax.spines.values():
                spine.set_linewidth(0.8 * s)
            plot_line.set_linewidth(max(0.5, BASE_LW * s))

            # scale HUD text
            hud_text.set_fontsize(BASE_HUD * s)

        def _sync_to_canvas(width_px: int, height_px: int):
            """Resize embedded widget, figure (inches) and style to the Canvas pixel size."""
            tk_canvas.coords(window_id, 0, 0)
            tk_canvas.itemconfig(window_id, width=width_px, height=height_px)
            tk_canvas.configure(scrollregion=(0, 0, width_px, height_px))

            # match figure pixel area
            fig.set_size_inches(max(width_px, 1) / dpi, max(height_px, 1) / dpi, forward=True)

            # compute scale factor relative to a reference pixel box
            scale = min(width_px / BASE_W_PX, height_px / BASE_H_PX)
            _apply_style_scale(scale)

            try:
                mplt_canvas.draw_idle()
            except tk.TclError:
                pass

        def _on_canvas_configure(event):
            # ignore spurious tiny events during initial layout
            if event.width < 150 or event.height < 120:
                return
            _sync_to_canvas(event.width, event.height)

        tk_canvas.bind("<Configure>", _on_canvas_configure)

        # Initial sync once layout is ready
        try:
            window.refresh()
        except Exception:
            pass
        w = max(tk_canvas.winfo_width(), 1)
        h = max(tk_canvas.winfo_height(), 1)
        _sync_to_canvas(w, h)
        return fig, ax, plot_line, hud_text, mplt_canvas
    
    else:
        if layout == 'macos':
            fig = Figure(figsize=(6.6, 3.05), dpi=100)
            ax = fig.add_subplot(111)
            fig.subplots_adjust(left=0.11, right=0.98, top=0.93, bottom=0.16)
        elif layout == 'android':
            fig = Figure(figsize=(8.3, 2.9), dpi=100)
            ax = fig.add_subplot(111)
            fig.subplots_adjust(left=0.11, right=0.98, top=0.92, bottom=0.16)
        else:
            fig = Figure(figsize=(6.0, 3.0), dpi=100)
            ax = fig.add_subplot(111)

        (_plot_line,) = ax.plot([], [], lw=0.8)
        ax.set_xlabel("Wavelength [Å]")
        ax.set_ylabel("Flux")
        hud_text = ax.text(0.995, 0.01, "", transform=ax.transAxes,
                        ha='right', va='bottom', fontsize=9,
                        color='black', zorder=50, clip_on=False)

        _preview_canvas = FigureCanvasTkAgg(fig, window[preview_key].TKCanvas)
        widget = _preview_canvas.get_tk_widget()
        widget.pack(side='top', fill='both', expand=0)
        _preview_canvas.draw()

        return fig, ax, _plot_line, hud_text, _preview_canvas


# --------------------------------------------------------------
# Real-time SNR provider for the preview
# --------------------------------------------------------------
def snr_provider(lam_x: float, ax, preview_interactor):
    data = getattr(ax, "_last_xydata", None)
    if not data:
        return None
    x, y = data

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if not np.any(ok):
        return None
    x = x[ok]; y = y[ok]
    n = x.size
    if n < 9:
        return None

    # Intercepting the nearest indices
    ascending = bool(x[0] <= x[-1])
    if ascending:
        idx = int(np.searchsorted(x, lam_x))
        if idx <= 0: idx = 0
        elif idx >= n: idx = n - 1
    else:
        idx = int(np.argmin(np.abs(x - lam_x)))

    # Window size
    mode = getattr(preview_interactor, "snr_mode", "points")
    if mode == "angstrom":
        dx = np.median(np.abs(np.diff(x))) if n > 1 else np.nan
        if not np.isfinite(dx) or dx <= 0:
            halfwin_pts = int(max(8, getattr(preview_interactor, "snr_halfwin_pts", 20)))
        else:
            half_A = float(getattr(preview_interactor, "snr_halfwin_A", 10.0))
            halfwin_pts = int(max(8, min(300, round(half_A / dx))))
    else:
        halfwin_pts = int(max(8, min(300, getattr(preview_interactor, "snr_halfwin_pts", 20))))

    lo = max(0, idx - halfwin_pts)
    hi = min(n, idx + halfwin_pts + 1)
    if hi - lo < 9:
        return None

    xs = x[lo:hi]; ys = y[lo:hi]
    w_ok = np.isfinite(xs) & np.isfinite(ys)
    xs = xs[w_ok]; ys = ys[w_ok]
    if xs.size < 9:
        return None

    # detrend
    try:
        p = np.polyfit(xs, ys, 1)
        baseline = p[0]*xs + p[1]
        resid = ys - baseline
    except Exception:
        med = float(np.nanmedian(ys))
        resid = ys - med

    resid = resid[np.isfinite(resid)]
    if resid.size < 5:
        return None

    # sigma
    mad = float(np.nanmedian(np.abs(resid)))
    sigma = mad * 1.4826 if mad > 0 else float(np.nanstd(resid))
    if not np.isfinite(sigma) or sigma <= 0:
        amp = float(np.nanmax(np.abs(ys))) if np.any(np.isfinite(ys)) else 1.0
        sigma = max(amp * 1e-12, 1e-12)

    # signal
    if ascending:
        signal = float(np.interp(lam_x, x, y))
    else:
        signal = float(np.interp(lam_x, x[::-1], y[::-1]))

    return abs(signal) / sigma
