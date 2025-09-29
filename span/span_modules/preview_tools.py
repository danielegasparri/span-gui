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

import matplotlib.ticker as mticker
from dataclasses import replace
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

try: #try local import if executed as script
    #GUI import
    from FreeSimpleGUI_local import FreeSimpleGUI as sg

except ModuleNotFoundError: #local import if executed as package
    #GUI import
    from span.FreeSimpleGUI_local import FreeSimpleGUI as sg

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
    snr_mode: str = "points"     # "points" or "angstrom" mode for SNR. I set points because is more general.
    snr_halfwin_pts: int = 25      #+-25 points for SNR calculation
    snr_halfwin_A: float = 20.0    #+- 20 A for SNR calculation

    def __post_init__(self):
        self._home_xlim = self.ax.get_xlim()
        self._home_ylim = self.ax.get_ylim()
        self._is_panning = False
        self._last_xy = (None, None)
        self._last_move_ts = 0.0
        self._last_hud_msg = "" 
        fig = self.ax.figure
        self._cids = [
            fig.canvas.mpl_connect('scroll_event', self._on_scroll),
            fig.canvas.mpl_connect('button_press_event', self._on_press),
            fig.canvas.mpl_connect('button_release_event', self._on_release),
            fig.canvas.mpl_connect('motion_notify_event', self._on_move),
        ]

    @staticmethod
    def _zoom_around(lims, centre, scale):
        a, b = lims
        return (centre + (a - centre) * scale,
                centre + (b - centre) * scale)

    def reset_view(self):
        self.ax.set_xlim(self._home_xlim)
        self.ax.set_ylim(self._home_ylim)
        self.ax.figure.canvas.draw_idle()

    def update_home(self):
        self._home_xlim = self.ax.get_xlim()
        self._home_ylim = self.ax.get_ylim()

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

        if getattr(event, 'dblclick', False):
            if is_left:
                self.reset_view()
            return

        if is_left:
            self._is_panning = True
            # 🔧 Salva posizione iniziale in pixel
            self._press_pixel = (event.x, event.y)
            self._pan_xlim0 = self.ax.get_xlim()
            self._pan_ylim0 = self.ax.get_ylim()
            try:
                self.ax.set_cursor(1)
            except Exception:
                pass


    def _on_release(self, event):
        if self._is_panning:
            self._is_panning = False
            self._press_pixel = None
            try:
                self.ax.set_cursor(0)
            except Exception:
                pass


    def _on_move(self, event):
        now = time.time()
        if (now - self._last_move_ts) * 1000.0 >= self.throttle_ms:
            self._last_move_ts = now

            # === HUD update (solo se NON stai facendo pan) ===
            if not self._is_panning and event.inaxes == self.ax and event.xdata is not None:
                lam_log = float(event.xdata)
                lam = 10**lam_log
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

            # === PAN update (solo se stai trascinando) ===
            if self._is_panning and event.inaxes == self.ax and self._press_pixel is not None:
                dx_pix = event.x - self._press_pixel[0]
                dy_pix = event.y - self._press_pixel[1]

                # Converti spostamento pixel → dati
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

    def __init__(self, ax, line, hud_text=None):
        self.ax = ax
        self.line = line  # the Line2D object of the spectrum
        self.hud_text = hud_text
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
    if layout == 'windows':
        fig = Figure(figsize=(5.2, 3.2), dpi=100)
        ax = fig.add_subplot(111)
        fig.subplots_adjust(left=0.13, right=0.98, top=0.93, bottom=0.16)
    elif layout == 'linux':
        fig = Figure(figsize=(6.6, 3.05), dpi=100)
        ax = fig.add_subplot(111)
        fig.subplots_adjust(left=0.11, right=0.98, top=0.93, bottom=0.16)
    elif layout == 'macos':
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
