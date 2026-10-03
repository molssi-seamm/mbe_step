# -*- coding: utf-8 -*-

"""The graphical part of a MBE step"""

import pprint  # noqa: F401
import tkinter as tk
import tkinter.ttk as ttk

import seamm
from seamm_util import ureg, Q_, units_class  # noqa: F401
import seamm_widgets as sw

#: The parameters of each group, in order
GROUPS = {
    "levels": (
        "Levels of theory",
        ("high level", "molecular low level", "periodic low level", "cell low level"),
    ),
    "fragments": (
        "Fragments",
        (
            "maximum order",
            "distance criterion",
            "cutoffs",
            "pair cutoff",
            "pair cutoff table",
            "triple rule",
            "triple cutoff",
            "triple cutoff table",
        ),
    ),
    "assignment": (
        "Periodic low level",
        ("periodic monomers", "periodic pair cutoff"),
    ),
    "labels": ("Labels", ("energy offsets", "extxyz file")),
    "execution": (
        "Execution",
        (
            "molecular ranks",
            "molecular memory",
            "molecular bundle",
            "periodic ranks",
            "periodic memory",
            "periodic bundle",
            "cell ranks",
            "cell memory",
            "bundle walltime",
            "archive",
        ),
    ),
}

#: Parameters whose value changes the layout
DRIVERS = ("maximum order", "cutoffs", "triple rule", "periodic low level")


class TkMbe(seamm.TkNode):
    """
    The graphical part of a MBE step in a flowchart.

    The dialog shows only what applies (Paul's rule: invalid combinations are
    not offered): the triple controls only for order 3, the pair controls only
    from order 2, single cutoffs or the type-pair tables but not both, and the
    periodic-assignment controls only with a periodic low level. The order is
    capped at 3 here; 4-body terms are designed for in seamm_mbe but not yet
    validated. The level fields offer the installed model chemistries and
    accept a typed one, or a $variable.
    """

    def __init__(
        self,
        tk_flowchart=None,
        node=None,
        canvas=None,
        x=None,
        y=None,
        w=200,
        h=50,
    ):
        """Initialize a graphical node."""
        self.dialog = None
        super().__init__(
            tk_flowchart=tk_flowchart,
            node=node,
            canvas=canvas,
            x=x,
            y=y,
            w=w,
            h=h,
        )

    def create_dialog(self):
        """Create the dialog: a Parameters tab and the standard Results tab."""
        super().create_dialog(title="MBE", widget="notebook", results_tab=True)
        P = self.node.parameters
        frame = self["frame"]

        self["structures frame"] = ttk.LabelFrame(
            frame, borderwidth=4, relief="sunken", text="Structures", padding=10
        )
        self.create_structure_selection_widgets(self["structures frame"])

        for group, (title, keys) in GROUPS.items():
            subframe = self[f"{group} frame"] = ttk.LabelFrame(
                frame, borderwidth=4, relief="sunken", text=title, padding=10
            )
            for key in keys:
                self[key] = P[key].widget(subframe)
        for key in DRIVERS:
            self[key].bind("<<ComboboxSelected>>", self.reset_dialog)
            self[key].bind("<Return>", self.reset_dialog)
            self[key].bind("<FocusOut>", self.reset_dialog)
        for key in ("pair cutoff table", "triple cutoff table", "energy offsets"):
            if hasattr(self[key], "entry"):
                self[key].entry.configure(width=40)

        self._fill_levels()
        self.reset_dialog()

    def _fill_levels(self):
        """Offer the installed model chemistries in the level fields."""
        try:
            from model_chemistry_step.model_chemistry import (
                discover_model_chemistries,
            )

            molecular = sorted(discover_model_chemistries())
            periodic = sorted(discover_model_chemistries(periodic_only=True))
        except Exception:
            molecular = periodic = []
        values = {
            "high level": ["current model chemistry"] + molecular,
            "molecular low level": molecular,
            "periodic low level": ["none"] + molecular,
            "cell low level": ["automatic"] + periodic + molecular,
        }
        for key, choices in values.items():
            if hasattr(self[key], "combobox"):
                self[key].combobox.configure(values=choices, width=50)

    def reset_dialog(self, widget=None):
        """Lay out the dialog for the current choices."""
        frame = self["frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        row = 0
        self["structures frame"].grid(row=row, column=0, sticky=tk.EW, pady=5)
        sframe = self["structures frame"]
        for slave in sframe.grid_slaves():
            slave.grid_forget()
        _, widgets = self.layout_structure_selection(row=0)
        sw.align_labels(widgets, sticky=tk.E)
        row += 1

        shown = self.shown()
        for group in GROUPS:
            subframe = self[f"{group} frame"]
            for slave in subframe.grid_slaves():
                slave.grid_forget()
            keys = [k for k in GROUPS[group][1] if k in shown]
            if not keys:
                continue
            subframe.grid(row=row, column=0, sticky=tk.EW, pady=5)
            row += 1
            for r, key in enumerate(keys):
                self[key].grid(row=r, column=0, sticky=tk.EW)
            sw.align_labels([self[k] for k in keys], sticky=tk.E)
            subframe.columnconfigure(0, weight=1)
        frame.columnconfigure(0, weight=1)
        return row

    def shown(self):
        """The parameters that apply to the current choices."""
        try:
            order = int(self._value("maximum order"))
        except (TypeError, ValueError):
            order = 3  # e.g. a $variable: show everything that might apply
        table = self._value("cutoffs") != "single value"
        triples = order >= 3 and self._value("triple rule") != "none"
        periodic = self._value("periodic low level") != "none"

        shown = set(GROUPS["levels"][1]) | set(GROUPS["labels"][1])
        shown |= {"maximum order", "distance criterion"}
        shown |= {"molecular ranks", "molecular memory", "molecular bundle"}
        shown |= {"cell ranks", "cell memory", "bundle walltime", "archive"}
        if order >= 2:
            shown.add("cutoffs")
            shown.add("pair cutoff table" if table else "pair cutoff")
        if order >= 3:
            shown.add("triple rule")
        if triples:
            shown.add("triple cutoff table" if table else "triple cutoff")
        if periodic:
            shown |= set(GROUPS["assignment"][1])
            shown |= {"periodic ranks", "periodic memory", "periodic bundle"}
        return shown

    def _value(self, key):
        """The current value of a widget, falling back to the parameter's."""
        try:
            return self[key].get()
        except Exception:
            return self.node.parameters[key].value

    def right_click(self, event):
        """Handle a right-click: add the Edit... item and post the menu."""
        super().right_click(event)
        self.popup_menu.add_command(label="Edit..", command=self.edit)
        self.popup_menu.tk_popup(event.x_root, event.y_root, 0)
