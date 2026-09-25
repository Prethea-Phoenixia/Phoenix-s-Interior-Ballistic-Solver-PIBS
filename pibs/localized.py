"""Localized tkinter widgets for PIBS.

All widgets support runtime language switching and can be enabled/disabled or
locked during computation.

Class hierarchy:

    Localizable (ABC)      — anything with a localize() method.
    LocalizedLabel         — helper (not a widget): tooltip StringVar +
                             localization keys for a target widget.
    LocalizableWidget      — base for data-bearing widgets; owns a tkinter
                             Variable and tracks enable/disable state.
    ComputableWidget (ABC) — adds actual tkinter state application for
                             widgets locked during computation.
    Descriptive (ABC)      — mixin for widgets serialized to JSON
                             (get()/reset()/get_descriptive()).

Concrete widgets:

    Loc2LineDisp / Loc3LineDisp — read-only label + entry display.
    Loc2Input / Loc3Input       — editable label + entry (+ unit). Lockable.
    LocDropdown                 — combobox with localized options. Lockable.
    LocLabelCheck               — checkbutton with localized label. Lockable.
    LocLabelFrame               — LabelFrame with localized title + tooltip.
    ScrollableText              — Text with auto-attached scrollbars.

Helpers: RowBuilder (auto-increments row index) and LocalizedFrame (top-level
container managing the localization dict, language switching, and widget
registry; child frames share the parent's lang variable and widget list).
"""

from __future__ import annotations

import tkinter
import tkinter.font
from abc import ABC, abstractmethod
from tkinter import BooleanVar, Event, Frame, Menu, StringVar, Text, Tk, Toplevel, ttk
from tkinter.font import Font
from typing import Any, Callable, Literal, TypeVar, Union, overload

from . import FONTNAME, FONTSIZE
from .misc import format_float_input
from .tip import create_tool_tip

# Type aliases for localization functions
localize_function_type = Union[
    Callable[[str, bool], str],
    Callable[[str], str],
]
"""Type alias for localization functions that can optionally force default language."""

placeholder_loc_func: localize_function_type = lambda _: ""

# localization.json stores mostly strings, but structured entries exist
# (e.g. "columnList" maps gun-type values to lists of column names).
localization_dict_type = dict[str, dict[str, str | dict[str, list[str]]]]


class Localizable(ABC):
    @abstractmethod
    def localize(self, *args: Any, **kwargs: Any) -> None: ...

    def inhibit(self) -> None:
        pass

    def disinhibit(self) -> None:
        pass


class LocalizedLabel:
    """Manages a localized label text and optional tooltip on a target widget."""

    def __init__(
        self,
        target_widget: ttk.Widget,
        loc_func: localize_function_type,
        label_key: str,
        font: Font,
        tooltip_key: str = "",
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
        loc: LocalizedLabel,
        var: tkinter.Variable,
        loc_func: localize_function_type,
        all_localized: list[LocalizableWidget] | None = None,
        nominal_state: Literal["normal", "disabled", "readonly"] = "readonly",
        **kwargs: Any,
    ) -> None:
        self.loc = loc
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

    def set(self, val: Any) -> None:
        self.var.set(val)

    def get(self) -> Any:
        return self.var.get()

    def disable(self) -> None:
        self.target_state = "disabled"

    def enable(self) -> None:
        self.target_state = self.nominal_state


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


class Descriptive(LocalizableWidget):

    def __init__(self, *args, desc_label_key, **kwargs):
        super().__init__(*args, **kwargs)
        self.desc_label_key = desc_label_key

    def get_descriptive(self) -> str:
        """Returns the JSON key to this widget's value."""
        if self.desc_label_key is None:
            return ""
        if self.desc_label_key:
            return self.loc_func(self.desc_label_key, True)
        return self.loc_func(self.loc.label_key, True)

    @abstractmethod
    def reset(self, *args: Any, **kwargs) -> None:
        """Resets the value of this widget to the default value."""
        ...


class Loc2LineDisp(LocalizableWidget):
    def __init__(
        self,
        parent: ttk.Widget,
        font: Font,
        loc_func: localize_function_type = placeholder_loc_func,
        row: int = 0,
        col: int = 0,
        label_loc_key: str = "",
        default: str = "",
        entry_width: int = 16,
        justify: Literal["left", "right", "center"] = "center",
        tooltip_loc_key: str = "",
        all_localized: list[LocalizableWidget] | None = None,
    ) -> None:
        e = StringVar()
        e.set(default)

        self.label_widget = ttk.Label(parent, text=loc_func(label_loc_key))
        self.label_widget.grid(row=row, column=col, sticky="nsew", padx=2, pady=2)

        self.entry_widget = ttk.Entry(parent, textvariable=e, width=entry_width, state="disabled", justify=justify)
        self.entry_widget.grid(row=row + 1, column=0, sticky="nsew", padx=2, pady=2)

        loc = LocalizedLabel(
            target_widget=self.label_widget,
            loc_func=loc_func,
            label_key=label_loc_key,
            tooltip_key=tooltip_loc_key,
            font=font,
        )
        self.default = default

        super().__init__(loc_func=loc_func, all_localized=all_localized, var=e, loc=loc)

    def localize(self, new_loc_key: str = "", new_tooltip_key: str = "") -> None:
        self.loc.localize(new_loc_key, new_tooltip_key)
        self.label_widget.config(text=self.loc_func(self.loc.label_key))

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
        font: Font,
        loc_func: localize_function_type = placeholder_loc_func,
        row: int = 0,
        col: int = 0,
        label_loc_key: str = "",
        default_up: str = "",
        default_dn: str = "",
        entry_width: int = 16,
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
        font: Font,
        loc_func: localize_function_type = placeholder_loc_func,
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

        loc = LocalizedLabel(
            target_widget=lb, loc_func=loc_func, label_key=label_loc_key, tooltip_key=tooltip_loc_key, font=font
        )

        self.default = default
        self.label_widget = lb
        self.input_widget = en
        self.row = row
        self.dtype = dtype

        super().__init__(
            state_target=en,
            var=e,
            loc=loc,
            loc_func=loc_func,
            all_localized=all_localized,
            desc_label_key=desc_label_key,
            nominal_state="normal",
        )

    def localize(self, new_loc_key: str = "", new_tooltip_key: str = "") -> None:
        self.loc.localize(new_loc_key, new_tooltip_key)
        self.label_widget.config(text=self.loc_func(self.loc.label_key))

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


class Loc3Input(Loc2Input):
    def __init__(
        self,
        parent: ttk.Widget,
        font: Font,
        loc_func: localize_function_type = placeholder_loc_func,
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


T = TypeVar("T")


class LocDropdown(ComputableWidget, Descriptive):
    def __init__(
        self,
        parent: ttk.Widget,
        font: tkinter.font.Font,
        loc_func: localize_function_type = placeholder_loc_func,
        str_obj_dict: dict[str, T] | None = None,
        all_localized: list[LocalizableWidget] | None = None,
        desc_label_key: str = "",
        tooltip_loc_key: str = "",
    ) -> None:
        self.str_obj_dict: dict[str, T] = str_obj_dict if str_obj_dict else {"": ""}
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

        self.widget = widget
        loc = LocalizedLabel(
            target_widget=widget, loc_func=loc_func, label_key="", tooltip_key=tooltip_loc_key, font=font
        )

        super().__init__(
            state_target=widget,
            var=var,
            loc_func=loc_func,
            all_localized=all_localized,
            desc_label_key=desc_label_key,
            nominal_state="readonly",
            loc=loc,
        )

    def localize(self) -> None:
        index = self.widget["values"].index(self.var.get())
        self.loc_str_obj_dict = self._update_loc_str_obj_dict(self.loc_func)
        self.widget.config(values=tuple(self.loc_str_obj_dict.keys()))
        self.widget.current(index)
        self.loc.localize()

    def _update_loc_str_obj_dict(self, loc_func: localize_function_type) -> dict[str, object]:
        return {loc_func(k): v for k, v in self.str_obj_dict.items()}

    def get(self) -> str:
        return self.get_obj().__str__()

    def get_obj(self) -> T:
        return self.loc_str_obj_dict[self.var.get()]

    def set_by_str(self, string: str) -> None:
        self.widget.set(self.widget["values"][list(self.str_obj_dict.keys()).index(string)])

    def set_by_obj(self, obj: T) -> None:
        index = list(self.str_obj_dict.values()).index(obj)
        self.widget.current(index)

    def set(self, val: str) -> None:
        self.set_by_str(val)

    def grid(self, **kwargs: Any) -> None:
        self.widget.grid(**kwargs)

    def reset(self, str_obj_dict: dict[str, T] | None = None, overwrite: bool = True) -> None:
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
        font: Font,
        loc_func: localize_function_type = placeholder_loc_func,
        label_loc_key: str = "",
        tooltip_loc_key: str = "",
        all_localized: list[Localizable] | None = None,
        **kwargs: Any,
    ) -> None:
        super().__init__(*args, text=loc_func(label_loc_key), **kwargs)
        self.loc_key = label_loc_key
        self.loc_func = loc_func
        self.loc = LocalizedLabel(
            target_widget=self, loc_func=loc_func, label_key=label_loc_key, tooltip_key=tooltip_loc_key, font=font
        )

        if isinstance(all_localized, list):
            all_localized.append(self)

    def localize(self) -> None:
        self.loc.localize()
        self.config(text=self.loc_func(self.loc_key))


class LocLabelCheck(ComputableWidget, Descriptive):
    def __init__(
        self,
        parent: ttk.Widget,
        font: Font,
        loc_func: localize_function_type = placeholder_loc_func,
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

        self.check_widget = check_widget
        loc = LocalizedLabel(
            target_widget=check_widget,
            loc_func=loc_func,
            label_key=label_loc_key,
            tooltip_key=tooltip_loc_key,
            font=font,
        )

        super().__init__(
            state_target=check_widget,
            var=var,
            loc=loc,
            loc_func=loc_func,
            all_localized=all_localized,
            desc_label_key=desc_label_key,
            nominal_state="normal",
        )

    def localize(self, new_loc_key: str = "") -> None:
        self.loc.localize(new_loc_key)
        self.check_widget.config(text=self.loc_func(self.loc.label_key))

    def remove(self) -> None:
        self.check_widget.grid_remove()

    def restore(self) -> None:
        self.check_widget.grid()

    def set(self, val: bool) -> None:
        self.var.set(val)

    def reset(self, *args) -> None:
        self.var.set(self.default)


class ScrollableText(Text):
    """Text widget with attached scrollbars; self-grids at (row, col).

    Vertical scrollbar at (row, col+1); optional horizontal scrollbar at
    (row+1, col), which extends the vertical scrollbar's rowspan to 2.
    Callers consuming extra rows (e.g. RowBuilder) must advance past them.
    """

    def __init__(self, parent: ttk.Widget, row: int, col: int = 0, hscroll: bool = False, **text_kwargs: Any) -> None:
        text_kwargs.setdefault("wrap", "word")
        text_kwargs.setdefault("height", 0)
        text_kwargs.setdefault("width", 0)
        text_kwargs.setdefault("font", (FONTNAME, FONTSIZE))
        super().__init__(parent, **text_kwargs)

        vscroll = ttk.Scrollbar(parent, orient="vertical", command=self.yview)
        self.configure(yscrollcommand=vscroll.set)
        self.grid(row=row, column=col, sticky="nsew")
        vscroll.grid(row=row, column=col + 1, rowspan=2 if hscroll else 1, sticky="nsew")
        if hscroll:
            h_scroll = ttk.Scrollbar(parent, orient="horizontal", command=self.xview)
            self.configure(xscrollcommand=h_scroll.set)
            h_scroll.grid(row=row + 1, column=col, sticky="nsew")


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

    def label_frame(self, grid_kwargs: dict[str, Any] | None = None, **kwargs) -> LocLabelFrame:
        widget = self.localized_frame.add_localized_label_frame(self.parent, **kwargs)
        widget.grid(row=self.row, column=0, sticky="nsew", padx=2, pady=2, **(grid_kwargs or {}))
        widget.columnconfigure(0, weight=1)
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
        font: Font,
        localization_dict: localization_dict_type,
        default_lang: str,
        *args: Any,
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

    def add_localized_2_line_display(self, *args, **kwargs):
        return self._make_localized(Loc2LineDisp, *args, **kwargs)

    def add_localized_3_line_display(self, *args, **kwargs):
        return self._make_localized(Loc3LineDisp, *args, **kwargs)

    def add_localized_label_check(self, *args, **kwargs):
        return self._make_localized(LocLabelCheck, *args, **kwargs)

    def add_localized_2_input(self, *args, **kwargs):
        return self._make_localized(Loc2Input, *args, **kwargs)

    def add_localized_3_input(self, *args, **kwargs):
        return self._make_localized(Loc3Input, *args, **kwargs)

    def add_localized_dropdown(self, *args, **kwargs):
        return self._make_localized(LocDropdown, *args, **kwargs)

    def add_localized_label_frame(self, *args, **kwargs):
        return self._make_localized(LocLabelFrame, *args, **kwargs)

    @overload
    def get_loc_str(self, name: Literal["columnList"], force_default: bool = False) -> dict[str, list[str]]: ...

    @overload
    def get_loc_str(self, name: str, force_default: bool = False) -> str: ...

    def get_loc_str(self, name: str, force_default: bool = False) -> str | dict[str, list[str]]:
        """Return the localized entry for ``name``.

        Falls back to English when ``name`` is absent from the current
        language, then to ``name`` itself when absent from English too — the
        ``str`` fallback only occurs if the key is missing from the
        localization file entirely, which the overloads do not model.
        """
        try:
            lang = "English" if force_default else self.lang_var.get()
            return self.localization_dict[lang][name]
        except KeyError:
            try:
                return self.localization_dict["English"][name]
            except KeyError:
                return name
