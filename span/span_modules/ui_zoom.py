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

# Class and functions to handle the zooming events in the GUI.

from dataclasses import dataclass
from typing import Dict, Optional
import tkinter.ttk as ttk 

try:
    from FreeSimpleGUI_local import FreeSimpleGUI as sg
except ModuleNotFoundError:
    from span.FreeSimpleGUI_local import FreeSimpleGUI as sg

import sys
import tkinter as tk
import tkinter.font as tkfont
import matplotlib

_TK_FONT_NAMES = (
    "TkDefaultFont",
    "TkTextFont",
    "TkFixedFont",
    "TkMenuFont",
    "TkHeadingFont",
    "TkIconFont",
    "TkTooltipFont",
)

@dataclass
class _FontRecord:
    size: int

class ZoomManager:
    _instance: Optional["ZoomManager"] = None

    def __init__(self):
        self._scale = 1.0
        self._base_fonts = {}
        self._attached_roots = {}
        self._base_dpi = None
        self._base_tk_scaling = None 
        
    @classmethod
    def get(cls) -> "ZoomManager":
        if cls._instance is None:
            cls._instance = ZoomManager()
        return cls._instance

    def attach_window(self, window):
        root = window.TKroot
        self._attached_roots[id(root)] = root

        if self._base_tk_scaling is None:
            try:
                self._base_tk_scaling = float(root.tk.call('tk', 'scaling'))
            except Exception:
                self._base_tk_scaling = 1.0

        if not self._base_fonts:
            for name in _TK_FONT_NAMES:
                try:
                    f = tkfont.nametofont(name)
                    self._base_fonts[name] = _FontRecord(size=f.cget("size"))
                except tk.TclError:
                    pass
        if self._base_dpi is None:
            self._base_dpi = matplotlib.rcParams.get("figure.dpi", 100.0)

        try:
            self._apply_to_root(root)
        except Exception:
            pass
        
    def set_scale(self, scale: float):
        scale = max(0.6, min(2.5, float(scale)))
        self._scale = scale         
        for root in list(self._attached_roots.values()):
            self._apply_to_root(root)
        if self._base_dpi:
            matplotlib.rcParams["figure.dpi"] = self._base_dpi * self._scale

        # When resetting zooming, keep the original settings
        for fn in getattr(self, "_redraw_hooks", []):
            try: fn()
            except Exception: pass
            
            
    def zoom_in(self, step: float = 0.1):
        self.set_scale(self._scale + step)

    def zoom_out(self, step: float = 0.1):
        self.set_scale(self._scale - step)

    def current_scale(self) -> float:
        return self._scale

    def _apply_to_root(self, root):
            s = self._scale
            base = self._base_tk_scaling or 1.0
            desired = base * s

            try:
                root.tk.call('tk', 'scaling', desired)
            except tk.TclError:
                pass

            # 1) Named fonts (TkDefaultFont, ecc.)
            for name, rec in self._base_fonts.items():
                try:
                    f = tkfont.nametofont(name)
                    if abs(s - 1.0) < 1e-6:
                        new_size = rec.size
                    else:
                        new_size = max(6, int(round(rec.size * s)))
                    if f.cget("size") != new_size:
                        f.configure(size=new_size)
                except tk.TclError:
                    continue
            
            # 2) Explicit fonts (sg.Text, sg.Frame, ecc.)
            try:
                if not hasattr(self, "_widget_base_fonts"):
                    self._widget_base_fonts = {}

                if id(root) not in self._widget_base_fonts:
                    self._widget_base_fonts[id(root)] = {}

                    def iter_widgets(w):
                        yield w
                        for c in w.winfo_children():
                            yield from iter_widgets(c)

                    for w in iter_widgets(root):
                        try:
                            current_font = w.cget("font")
                            if not current_font:
                                continue
                            f = tkfont.Font(root=root, font=current_font)
                            self._widget_base_fonts[id(root)][str(w)] = {
                                "family": f.cget("family"),
                                "size": int(f.cget("size")),
                                "weight": "bold" if f.cget("weight") == "bold" else "normal",
                                "slant": "italic" if f.cget("slant") == "italic" else "roman",
                                "underline": int(f.cget("underline")),
                                "overstrike": int(f.cget("overstrike")),
                                "original": current_font,
                            }
                        except Exception:
                            continue

                for widget_id, bf in self._widget_base_fonts[id(root)].items():
                    try:
                        w = root.nametowidget(widget_id)
                    except Exception:
                        continue

                    if abs(s - 1.0) < 1e-6:
                        # Reset: restore original font
                        try:
                            w.configure(font=bf["original"])
                            continue
                        except Exception:
                            pass

                    # Otherwise apply scaling
                    new_size = max(6, int(round(bf["size"] * s)))
                    try:
                        w.configure(font=(bf["family"], new_size, bf["weight"], bf["slant"]))
                    except Exception:
                        try:
                            w.configure(font=(bf["family"], new_size))
                        except Exception:
                            pass
            except Exception:
                pass


            # 3) Widget Canvas/Frame size
            try:
                if not hasattr(self, "_widget_base_sizes"):
                    self._widget_base_sizes = {}

                def iter_widgets(w):
                    yield w
                    for c in w.winfo_children():
                        yield from iter_widgets(c)

                for w in iter_widgets(root):
                    cls_name = w.winfo_class()
                    if cls_name in ("Canvas", "Frame"):
                        widget_id = str(w)
                        if widget_id not in self._widget_base_sizes:
                            base_w, base_h = w.winfo_reqwidth(), w.winfo_reqheight()
                            self._widget_base_sizes[widget_id] = (base_w, base_h)

                        base_w, base_h = self._widget_base_sizes[widget_id]

                        if sys.platform == "win32" or sys.platform == "darwin":
                            # Windows & macOS: Tk behaves consistently
                            if abs(s - 1.0) < 1e-6:
                                new_w, new_h = base_w, base_h
                            else:
                                new_w = max(50, int(base_w * s))
                                new_h = max(30, int(base_h * s))

                        else: #For Linux and maybe MacOS?
                            if abs(s - 1.0) < 1e-6:
                                new_w, new_h = base_w, base_h
                            
                                # FIX for Linux: Tk stay out of this!
                                if cls_name == "Canvas":
                                    try:
                                        w.pack_propagate(False)
                                        w.config(width=base_w, height=base_h)
                                    except Exception:
                                        pass
                            else:
                                new_w = max(50, int(base_w * s))
                                new_h = max(30, int(base_h * s))


                        try:
                            w.config(width=new_w, height=new_h)
                        except Exception:
                            pass
            except Exception:
                pass

            # 4) Padding
            pad = max(1, int(round(4 * s)))
            PAD_PATTERNS = (
                "*Button.padX", "*Button.padY",
                "*Label.padX",  "*Label.padY",
                "*Entry.padX",  "*Entry.padY",
            )
            try:
                if abs(s - 1.0) > 1e-6:
                    for pat in PAD_PATTERNS:
                        root.tk.call('option', 'add', pat, pad, 'interactive')
                    self._pads_applied = True
                else:
                    if getattr(self, "_pads_applied", False):
                        try:
                            root.tk.call('option', 'clear')
                        except tk.TclError:
                            pass
                        self._pads_applied = False
            except tk.TclError:
                pass

            # 5) Re-layout + force geometry
            try:
                root.update_idletasks()
                req_w, req_h = root.winfo_reqwidth(), root.winfo_reqheight()
                x, y = root.winfo_x(), root.winfo_y()
                root.geometry(f"{req_w}x{req_h}+{x}+{y}")
            except tk.TclError:
                pass


    def current_tk_scaling(self) -> float:
        base = self._base_tk_scaling or 1.0
        return base * self._scale

    def register_redraw(self, fn):
        if not hasattr(self, "_redraw_hooks"):
            self._redraw_hooks = []
        self._redraw_hooks.append(fn)


# wrapper for the subwindows scaling 
def open_subwindow(title, layout_or_factory, *, zm=None, **window_kwargs):
    layout = layout_or_factory() if callable(layout_or_factory) else layout_or_factory
    zm = zm or ZoomManager.get()

    # default
    window_kwargs.setdefault("finalize", True)
    window_kwargs.setdefault("resizable", True)
    window_kwargs.setdefault("modal", False)

    win = sg.Window(title, layout, **window_kwargs)
    zm.attach_window(win)           
    return win


def zoom_to(target: float, step: float = 0.1):
    """Zooming"""
    zm = ZoomManager.get()
    cur = zm.current_scale()
    # Clamp for safety
    target = max(0.6, min(2.5, float(target)))

    if abs(cur - target) < 1e-3:
        return

    # Direction
    if target > cur:
        while cur + 1e-9 < target:
            zm.zoom_in(step)
            cur = zm.current_scale()
    else:
        while cur - 1e-9 > target:
            zm.zoom_out(step)
            cur = zm.current_scale()

    zm.set_scale(target)
