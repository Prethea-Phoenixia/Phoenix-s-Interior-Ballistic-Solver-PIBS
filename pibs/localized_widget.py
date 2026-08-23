"""Localized tkinter widgets for PIBS.

Provides a hierarchy of localized, state-manageable widgets for the interior
ballistics UI. All widgets support runtime language switching and can be
programmatically enabled/disabled or locked during computation.

Class hierarchy:

    Localizable (ABC)
        Abstract base for anything with a localize() method.

    LocalizedLabel
        Helper (not a widget): manages a tooltip StringVar and localization
        keys for any target widget. Shared across widget types.

    LocalizableWidget (Localizable)
        Base for data-bearing widgets. Owns a tkinter.Variable, handles
        registration in the localization system, and tracks enable/disable
        state. inhibit()/disinhibit() are no-ops here — use ComputableWidget
        for widgets that must be locked during computation.

    ComputableWidget (LocalizableWidget, ABC)
        Extends LocalizableWidget with actual tkinter state application.
        Used by input widgets (entries, dropdowns, checkboxes) that need
        to be disabled while a calculation runs.

    Descriptive (ABC)
        Mixin for widgets whose values are serialized to JSON. Requires
        get(), get_descriptive() (the JSON key), and reset().

Concrete widgets:

    Loc2LineDisp   — Read-only: label + single entry (2 rows).
    Loc3LineDisp   — Read-only: label + two entries (3 rows).
    Loc2Input      — Editable: label + entry (1 row). Lockable.
    Loc3Input      — Editable: label + entry + unit label. Lockable.
    LocDropdown    — Combobox with localized options. Lockable.
    LocLabelCheck  — Checkbutton with localized label. Lockable.
    LocLabelFrame  — LabelFrame with localized title + tooltip.

Helpers:

    RowBuilder     — Auto-increments row index while creating widgets.
                     Delegates to LocalizedFrame factories.

    LocalizedFrame — Top-level container managing localization dict,
                     language switching, and widget registry. Can be nested:
                     child frames share the parent's language variable and
                     widget list, so a single language change propagates
                     through the entire tree.
"""

from __future__ import annotations

import tkinter
import warnings
from abc import ABC, abstractmethod
from tkinter import BooleanVar, Event, Frame, Menu, StringVar, Tk, Toplevel, ttk
from tkinter.font import Font
from typing import (
    Any,
    Callable,
    Literal,
    Union,
)

from .misc import format_float_input
from .tip import create_tool_tip

# Type aliases for localization functions
localize_function_type = Union[
    Callable[[str, bool], str],
    Callable[[str], str],
]
"""Type alias for localization functions that can optionally force default language."""

placeholder_loc_func: localize_function_type = lambda _: ""


def warn(func):
    def wrapped(*args, font=None, loc_func=None, all_localized=None, **kwargs):
        if loc_func or all_localized or font:
            warnings.warn(
                "LocalizedFrame sets the font, loc_fun and all_localized parameters, supplied arguments will be ignored."
            )
        return func(*args, **kwargs)

    return wrapped


class Localizable(ABC):
    @abstractmethod
    def localize(self, *args: Any, **kwargs: Any) -> None: ...


class LocalizedLabel:
    """Manages a localized label text and optional tooltip on a target widget."""

    def __init__(
        self,
        target_widget: ttk.Widget,
        loc_func: localize_function_type,
        label_key: str,
        tooltip_key: str = "",
        font: Font | None = None,
    ):
        self.loc_func = loc_func
        self.label_key = label_key
        self.tooltip_key = tooltip_key
        self.tooltip_var = StringVar(value=loc_func(tooltip_key))
        if tooltip_key:
            create_tool_tip(target_widget, self.tooltip_var, font=font)

    def localize(self, new_key: str = "", new_tooltip_key: str = "") -> None:
        if new_key:
            self.label_key = new_key
        if new_tooltip_key:
            self.tooltip_key = new_tooltip_key
        self.tooltip_var.set(self.loc_func(self.tooltip_key))


class LocalizableWidget(Localizable):
    """Base for widgets that participate in localization and state management.

    inhibit()/disinhibit() are no-ops by default — override via ComputableWidget
    for widgets that need to be locked during computation.
    """

    def __init__(
        self,
        *args: Any,
        var: tkinter.Variable,
        loc_func: localize_function_type,
        all_localized: list[LocalizableWidget] | None = None,
        nominal_state: Literal["normal", "disabled", "readonly"] = "readonly",
        **kwargs: Any,
    ) -> None:
        self.loc_func = loc_func
        if isinstance(all_localized, list):
            all_localized.append(self)
        self.nominal_state = nominal_state
        self.target_state = nominal_state
        self.var = var

    @abstractmethod
    def localize(self, *args: Any, **kwargs: Any) -> None:
        """Update the widget's text to the current language."""
        ...

    def trace_add(self, *args: Any) -> None:
        self.var.trace_add(*args)

    def set(self, val: Any) -> None: ...

    def get(self) -> Any:
        return self.var.get()

    def disable(self) -> None:
        self.target_state = "disabled"

    def enable(self) -> None:
        self.target_state = self.nominal_state

    def inhibit(self) -> None:
        pass

    def disinhibit(self) -> None:
        pass


class ComputableWidget(LocalizableWidget, ABC):
    """Extends LocalizableWidget with actual tkinter state application.

    Use this for input widgets that must be locked during computation.
    """

    def __init__(
        self,
        state_target: ttk.Widget,
        state_key: str = "state",
        **kwargs: Any,
    ):
        super().__init__(**kwargs)
        self._state_target = state_target
        self._state_key = state_key

    def _apply_tk_state(self, state: str) -> None:
        self._state_target.config(**{self._state_key: state})

    def disable(self) -> None:
        self._apply_tk_state("disabled")
        super().disable()

    def enable(self) -> None:
        self._apply_tk_state(self.nominal_state)
        super().enable()

    def inhibit(self) -> None:
        self._apply_tk_state("disabled")

    def disinhibit(self) -> None:
        self._apply_tk_state(self.target_state)


class Descriptive(ABC):
    @abstractmethod
    def get(self): ...

    @abstractmethod
    def get_descriptive(self) -> str:
        """Returns the JSON key to this widget's value."""
        ...

    @abstractmethod
    def reset(self, *args: Any, **kwargs) -> None:
        """Resets the value of this widget to the default value."""
        ...


class Loc2LineDisp(LocalizableWidget):
    def __init__(
        self,
        parent: ttk.Widget,
        loc_func: localize_function_type = placeholder_loc_func,
        font: Font | None = None,
        row: int = 0,
        col: int = 0,
        label_loc_key: str = "",
        default: str = "",
        entry_width: int = 10,
        justify: Literal["left", "right", "center"] = "center",
        tooltip_loc_key: str = "",
        all_localized: list[LocalizableWidget] | None = None,
    ) -> None:
        e = StringVar()
        e.set(default)
        super().__init__(loc_func=loc_func, all_localized=all_localized, var=e)
        self.label_widget = ttk.Label(parent, text=loc_func(label_loc_key))
        self.label_widget.grid(row=row, column=col, sticky="nsew", padx=2, pady=2)

        self.entry_widget = ttk.Entry(parent, textvariable=e, width=entry_width, state="disabled", justify=justify)
        self.entry_widget.grid(row=row + 1, column=0, sticky="nsew", padx=2, pady=2)

        self._loc = LocalizedLabel(self.label_widget, loc_func, label_loc_key, tooltip_loc_key, font)
        self.default = default

    def localize(self, new_loc_key: str = "", new_tooltip_key: str = "") -> None:
        self._loc.localize(new_loc_key, new_tooltip_key)
        self.label_widget.config(text=self.loc_func(self._loc.label_key))

    def set(self, val: str) -> None:
        self.var.set(val)

    def reset(self) -> None:
        self.set(self.default)

    def remove(self) -> None:
        self.label_widget.grid_remove()
        self.entry_widget.grid_remove()

    def restore(self) -> None:
        self.label_widget.grid()
        self.entry_widget.grid()


class Loc3LineDisp(Loc2LineDisp):
    def __init__(
        self,
        parent: ttk.Widget,
        loc_func: localize_function_type = placeholder_loc_func,
        font: Font | None = None,
        row: int = 0,
        col: int = 0,
        label_loc_key: str = "",
        default_up: str = "",
        default_dn: str = "",
        entry_width: int = 10,
        justify_up: Literal["left", "right", "center"] = "center",
        justify_dn: Literal["left", "right", "center"] = "center",
        tooltip_loc_key: str = "",
        all_localized: list[LocalizableWidget] | None = None,
    ) -> None:
        super().__init__(
            parent=parent,
            row=row,
            col=col,
            font=font,
            label_loc_key=label_loc_key,
            default=default_up,
            entry_width=entry_width,
            justify=justify_up,
            tooltip_loc_key=tooltip_loc_key,
            loc_func=loc_func,
            all_localized=all_localized,
        )
        self.aux_entry_var = StringVar()
        self.aux_entry_var.set(default_dn)
        self.aux_entry_widget = ttk.Entry(
            parent, textvariable=self.aux_entry_var, width=entry_width, state="disabled", justify=justify_dn
        )
        self.aux_entry_widget.grid(row=row + 2, column=0, sticky="nsew", padx=2, pady=2)
        self.aux_default = default_dn

    def set(self, val: tuple[str, str] | str) -> None:
        assert isinstance(val, tuple)
        val_1, val_2 = val
        super().set(val_1)
        self.aux_entry_var.set(val_2)

    def reset(self) -> None:
        self.set((self.default, self.aux_default))

    def remove(self) -> None:
        super().remove()
        self.aux_entry_widget.grid_remove()

    def restore(self) -> None:
        super().restore()
        self.aux_entry_widget.grid()


class Loc2Input(ComputableWidget, Descriptive):
    def __init__(
        self,
        parent: ttk.Widget,
        loc_func: localize_function_type = placeholder_loc_func,
        font: Font | None = None,
        default: str = "",
        row: int = 0,
        col: int = 0,
        label_loc_key: str = "",
        desc_label_key: str | None = "",
        label_width: int = 20,
        entry_width: int = 10,
        formatter: Callable[[Event, StringVar], None] = format_float_input,
        color: str = "",
        tooltip_loc_key: str = "",
        anchor: Literal["nw", "n", "ne", "w", "center", "e", "sw", "s", "se"] = "w",
        reverse: bool = False,
        all_localized: list[LocalizableWidget] | None = None,
        dtype: Callable[[str], Any] = lambda v: int(float(v)),
    ) -> None:
        e = StringVar(value=default)

        lb = ttk.Label(parent, text=loc_func(label_loc_key), width=label_width, anchor=anchor)
        lb.grid(row=row, column=col + (1 if reverse else 0), sticky="nsew", padx=2, pady=2)

        parent.rowconfigure(row, weight=0)

        en = ttk.Entry(
            parent,
            textvariable=e,
            width=entry_width,
            foreground=color,
            justify="center",
        )
        en.grid(row=row, column=col + (0 if reverse else 1), sticky="nsew", padx=2, pady=2)
        en.bind("<FocusOut>", lambda event: formatter(event, e))

        super().__init__(
            state_target=en,
            var=e,
            loc_func=loc_func,
            all_localized=all_localized,
            nominal_state="normal",
        )

        self._loc = LocalizedLabel(lb, loc_func, label_loc_key, tooltip_loc_key, font)

        self.default = default
        self.label_widget = lb
        self.input_widget = en
        self.row = row
        self.desc_label_key = desc_label_key
        self.dtype = dtype

    def localize(self, new_loc_key: str = "", new_tooltip_key: str = "") -> None:
        self._loc.localize(new_loc_key, new_tooltip_key)
        self.label_widget.config(text=self.loc_func(self._loc.label_key))

    def remove(self) -> None:
        self.label_widget.grid_remove()
        self.input_widget.grid_remove()

    def restore(self) -> None:
        self.label_widget.grid()
        self.input_widget.grid()

    def reset(self, *_) -> None:
        self.var.set(self.default)

    def get(self) -> Any:
        return self.dtype(self.var.get())

    def set(self, val: Any) -> None:
        self.var.set(val)

    def get_descriptive(self) -> str:
        if self.desc_label_key is None:
            return ""
        elif self.desc_label_key:
            return self.loc_func(self.desc_label_key, True)
        else:
            return self.loc_func(self._loc.label_key, True)


class Loc3Input(Loc2Input):
    def __init__(
        self,
        parent: ttk.Widget,
        loc_func: localize_function_type = placeholder_loc_func,
        font: Font | None = None,
        row: int = 0,
        col: int = 0,
        label_loc_key: str = "",
        desc_label_key: str | None = "",
        unit_text: str = "",
        default: str = "",
        label_width: int = 20,
        entry_width: int = 10,
        formatter: Callable[[Event, StringVar], None] = format_float_input,
        color: str = "",
        tooltip_loc_key: str = "",
        anchor: Literal["w", "e", "s", "n"] = "w",
        reverse: bool = False,
        all_localized: list[LocalizableWidget] | None = None,
        dtype: Callable[[str], Any] = lambda v: int(float(v)),
    ) -> None:
        super().__init__(
            parent=parent,
            font=font,
            row=row,
            col=col,
            label_loc_key=label_loc_key,
            desc_label_key=desc_label_key,
            default=default,
            label_width=label_width,
            entry_width=entry_width,
            formatter=formatter,
            color=color,
            tooltip_loc_key=tooltip_loc_key,
            anchor=anchor,
            reverse=reverse,
            loc_func=loc_func,
            all_localized=all_localized,
            dtype=dtype,
        )

        ulb = ttk.Label(parent, text=unit_text)
        ulb.grid(row=row, column=col + 2, sticky="nsew", padx=2, pady=2)
        self.unit_text = unit_text
        self.unit_label = ulb

    def remove(self) -> None:
        super().remove()
        self.unit_label.grid_remove()

    def restore(self) -> None:
        super().restore()
        self.unit_label.grid()

    def get_descriptive(self) -> str:
        if self.desc_label_key is None:
            return ""
        else:
            return super().get_descriptive() + (f" ({self.unit_text})" if self.unit_text else "")


class LocDropdown(ComputableWidget, Descriptive):
    def __init__(
        self,
        parent: ttk.Widget,
        loc_func: localize_function_type = placeholder_loc_func,
        font: Font | None = None,
        str_obj_dict: dict[str, object] | None = None,
        all_localized: list[LocalizableWidget] | None = None,
        desc_label_key: str = "",
        tooltip_loc_key: str = "",
    ) -> None:
        self.str_obj_dict: dict[str, object] = str_obj_dict if str_obj_dict else {"": ""}
        self.loc_str_obj_dict = self._update_loc_str_obj_dict(loc_func)
        var = StringVar()
        widget = ttk.Combobox(
            parent,
            textvariable=var,
            values=tuple(self.loc_str_obj_dict.keys()),
            justify="center",
            state="readonly",
        )
        widget.option_add("*TCombobox*Listbox.Justify", "center")
        widget.current(0)

        super().__init__(
            state_target=widget,
            var=var,
            loc_func=loc_func,
            all_localized=all_localized,
            nominal_state="readonly",
        )

        self.widget = widget
        self.desc_label_key = desc_label_key
        self._loc = LocalizedLabel(widget, loc_func, "", tooltip_loc_key, font)

    def localize(self) -> None:
        index = self.widget["values"].index(self.var.get())
        self.loc_str_obj_dict = self._update_loc_str_obj_dict(self.loc_func)
        self.widget.config(values=tuple(self.loc_str_obj_dict.keys()))
        self.widget.current(index)
        self._loc.localize()

    def _update_loc_str_obj_dict(self, loc_func: localize_function_type) -> dict[str, object]:
        return {loc_func(k): v for k, v in self.str_obj_dict.items()}

    def get(self) -> str:
        return self.get_obj().__str__()

    def get_obj(self) -> object:
        return self.loc_str_obj_dict[self.var.get()]

    def set_by_str(self, string: str) -> None:
        self.widget.set(self.widget["values"][list(self.str_obj_dict.keys()).index(string)])

    def set_by_obj(self, obj: object) -> None:
        index = list(self.str_obj_dict.values()).index(obj)
        self.widget.current(index)

    def set(self, val: str) -> None:
        self.set_by_str(val)

    def grid(self, **kwargs: Any) -> None:
        self.widget.grid(**kwargs)

    def get_descriptive(self) -> str:
        return self.loc_func(self.desc_label_key, True)

    def reset(self, str_obj_dict: dict[str, object] | None = None, overwrite: bool = True) -> None:
        if str_obj_dict is not None:
            self.str_obj_dict = str_obj_dict
            self.loc_str_obj_dict = self._update_loc_str_obj_dict(self.loc_func)
        self.widget["values"] = tuple(self.loc_str_obj_dict.keys())

        if overwrite:
            self.widget.current(0)
        elif self.widget.get() not in self.widget["values"]:
            self.widget.current(0)


class LocLabelFrame(ttk.LabelFrame, Localizable):
    def __init__(
        self,
        *args: Any,
        loc_func: localize_function_type = placeholder_loc_func,
        font: Font | None = None,
        label_loc_key: str = "",
        tooltip_loc_key: str = "",
        all_localized: list[LocalizableWidget] | None = None,
        **kwargs: Any,
    ) -> None:
        super().__init__(*args, text=loc_func(label_loc_key), **kwargs)
        self.loc_key = label_loc_key
        self.loc_func = loc_func

        self._loc = LocalizedLabel(self, loc_func, label_loc_key, tooltip_loc_key, font)

        if isinstance(all_localized, list):
            all_localized.append(self)

    def localize(self) -> None:
        self._loc.localize()
        self.config(text=self.loc_func(self.loc_key))


class LocLabelCheck(ComputableWidget, Descriptive):
    def __init__(
        self,
        parent: ttk.Widget,
        loc_func: localize_function_type = placeholder_loc_func,
        font: Font | None = None,
        default: bool = True,
        skip_grid: bool = False,
        row: int = 0,
        col: int = 0,
        columnspan: int = 1,
        label_loc_key: str = "",
        desc_label_key: str | None = "",
        tooltip_loc_key: str = "",
        width: int = 0,
        all_localized: list[LocalizableWidget] | None = None,
    ) -> None:
        self.default = default

        var = BooleanVar(value=default)
        check_widget = ttk.Checkbutton(parent, text=loc_func(label_loc_key), variable=var, width=width)
        if not skip_grid:
            check_widget.grid(row=row, column=col, sticky="nsew", columnspan=columnspan, padx=2, pady=2)

        super().__init__(
            state_target=check_widget,
            var=var,
            loc_func=loc_func,
            all_localized=all_localized,
            nominal_state="normal",
        )

        self.check_widget = check_widget
        self.desc_label_key = desc_label_key
        self._loc = LocalizedLabel(check_widget, loc_func, label_loc_key, tooltip_loc_key, font)

    def localize(self, new_loc_key: str = "") -> None:
        self._loc.localize(new_loc_key)
        self.check_widget.config(text=self.loc_func(self._loc.label_key))

    def remove(self) -> None:
        self.check_widget.grid_remove()

    def restore(self) -> None:
        self.check_widget.grid()

    def set(self, val: bool) -> None:
        self.var.set(val)

    def trace_add(self, *args: Any) -> None:
        self.var.trace_add(*args)

    def get_descriptive(self) -> str:
        if self.desc_label_key is None:
            return ""
        elif self.desc_label_key:
            return self.loc_func(self.desc_label_key, True)
        else:
            return self.loc_func(self._loc.label_key, True)

    def reset(self, *args) -> None:
        self.var.set(self.default)


class RowBuilder:
    """Helper to build vertically stacked widgets with automatic row management.

    Tracks both the parent frame (for widget placement) and the current row.
    Create a new builder for each parent frame.

    Usage:
        b = RowBuilder(localized_frame, parent_frame)
        self.cal_mm = b.input_3(label_loc_key="calLabel", unit_text="mm", default="50.0")
    """

    def __init__(self, localized_frame: LocalizedFrame, parent: ttk.Widget, start_row: int = 0) -> None:
        self.localized_frame = localized_frame
        self.parent = parent
        self.row = start_row

    @property
    def current_row(self) -> int:
        return self.row

    def reset(self, start_row: int = 0) -> None:
        self.row = start_row

    def next(self) -> int:
        self.row += 1
        return self.row

    def input_3(self, **kwargs) -> Loc3Input:
        kwargs.setdefault("parent", self.parent)
        kwargs.setdefault("row", self.row)
        widget = self.localized_frame.add_localized_3_input(**kwargs)
        self.row += 1
        return widget

    def input_2(self, **kwargs) -> Loc2Input:
        kwargs.setdefault("parent", self.parent)
        kwargs.setdefault("row", self.row)
        widget = self.localized_frame.add_localized_2_input(**kwargs)
        self.row += 1
        return widget

    def dropdown(self, grid_kwargs: dict[str, Any] | None = None, **kwargs) -> LocDropdown:
        kwargs.setdefault("parent", self.parent)
        widget = self.localized_frame.add_localized_dropdown(**kwargs)
        widget.grid(row=self.row, sticky="nsew", padx=2, pady=2, **(grid_kwargs or {}))
        self.row += 1
        return widget

    def label_frame(self, **kwargs) -> LocLabelFrame:
        kwargs.setdefault("parent", self.parent)
        widget = self.localized_frame.add_localized_label_frame(**kwargs)
        self.row += 1
        return widget

    def check(self, **kwargs) -> LocLabelCheck:
        kwargs.setdefault("parent", self.parent)
        kwargs.setdefault("row", self.row)
        widget = self.localized_frame.add_localized_label_check(**kwargs)
        self.row += 1
        return widget

    def display_2(self, **kwargs) -> Loc2LineDisp:
        kwargs.setdefault("parent", self.parent)
        kwargs.setdefault("row", self.row)
        widget = self.localized_frame.add_localized_2_line_display(**kwargs)
        self.row += 2
        return widget

    def display_3(self, **kwargs) -> Loc3LineDisp:
        kwargs.setdefault("parent", self.parent)
        kwargs.setdefault("row", self.row)
        widget = self.localized_frame.add_localized_3_line_display(**kwargs)
        self.row += 3
        return widget


class LocalizedFrame(Frame):
    def __init__(
        self,
        master: Tk | Toplevel | LocalizedFrame,
        *args: Any,
        font: Font | None = None,
        localization_dict: dict[str, dict[str, str]],
        default_lang: str,
        menubar: Menu | None = None,
        lang_var: StringVar | None = None,
        **kwargs: Any,
    ) -> None:
        super().__init__(master, *args, **kwargs)
        self.localization_dict = localization_dict

        if isinstance(master, LocalizedFrame):
            self.localized_widgets = master.localized_widgets
        else:
            self.localized_widgets: list[LocalizableWidget] = []

        self.font = font

        if isinstance(master, (Toplevel, Tk)):

            self.lang_var = StringVar(
                value=(
                    default_lang
                    if default_lang in self.localization_dict.keys()
                    else list(self.localization_dict.keys())[0]
                )
            )
            if menubar:
                lang_menu = Menu(menubar)
                menubar.add_cascade(label="Lang 语言", menu=lang_menu)

                for lang in localization_dict.keys():
                    lang_menu.add_radiobutton(label=lang, variable=self.lang_var, value=lang, command=self.change_lang)
        elif isinstance(master, LocalizedFrame):
            self.lang_var = master.lang_var
        else:
            if lang_var is None:
                raise ValueError("lang_var must be provided when master is not a Toplevel, Tk, or LocalizedFrame")
            self.lang_var = lang_var

    def change_lang(self, *args, **kwargs) -> None:
        for loc_widget in self.localized_widgets:
            loc_widget.localize()

    def _make_localized(self, cls, *args, **kwargs):
        return cls(*args, loc_func=self.get_loc_str, all_localized=self.localized_widgets, font=self.font, **kwargs)

    @warn
    def add_localized_2_line_display(self, *args, **kwargs):
        return self._make_localized(Loc2LineDisp, *args, **kwargs)

    @warn
    def add_localized_3_line_display(self, *args, **kwargs):
        return self._make_localized(Loc3LineDisp, *args, **kwargs)

    @warn
    def add_localized_label_check(self, *args, **kwargs):
        return self._make_localized(LocLabelCheck, *args, **kwargs)

    @warn
    def add_localized_2_input(self, *args, **kwargs):
        return self._make_localized(Loc2Input, *args, **kwargs)

    @warn
    def add_localized_3_input(self, *args, **kwargs):
        return self._make_localized(Loc3Input, *args, **kwargs)

    @warn
    def add_localized_dropdown(self, *args, **kwargs):
        return self._make_localized(LocDropdown, *args, **kwargs)

    @warn
    def add_localized_label_frame(self, *args, **kwargs):
        return self._make_localized(LocLabelFrame, *args, **kwargs)

    def get_loc_str(self, name: str, force_default: bool = False) -> str:
        try:
            lang = "English" if force_default else self.lang_var.get()
            return self.localization_dict[lang][name]
        except KeyError:
            try:
                return self.localization_dict["English"][name]
            except KeyError:
                return name
