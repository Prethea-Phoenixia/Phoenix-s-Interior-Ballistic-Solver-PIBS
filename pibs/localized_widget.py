from __future__ import annotations

import tkinter
import warnings
from abc import ABC, abstractmethod
from functools import wraps
from tkinter import BooleanVar, Event, Frame, Menu, StringVar, ttk, Tk, Toplevel
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
                "LocalizedFrame sets the font, loc_fun and all_localized parameters to "
                "instance of itself, supplied arguments will be ignored."
            )
        return func(*args, **kwargs)

    return wrapped


class Localizable(ABC):
    @abstractmethod
    def localize(self, *args: Any, **kwargs: Any) -> None: ...


class LocalizableWidget(Localizable):
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

    def get(self) -> Any:
        return self.var.get()

    def disable(self) -> None:
        self.target_state = "disabled"

    def enable(self) -> None:
        self.target_state = self.nominal_state

    @abstractmethod
    def inhibit(self) -> None: ...

    @abstractmethod
    def disinhibit(self) -> None: ...


class Descriptive(ABC):
    @abstractmethod
    def get_descriptive(self) -> str:
        """
        Returns the JSON key to this widget's value.
        """
        ...

    @abstractmethod
    def reset(self, *args: Any) -> None:
        """
        Resets the value of this widget to the default value.
        """
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
        entry_width: int = 20,
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

        self.loc_tooltip_var = StringVar(value=loc_func(tooltip_loc_key))
        create_tool_tip(self.label_widget, self.loc_tooltip_var, font=font)
        self.label_loc_key = label_loc_key
        self.tooltip_loc_key = tooltip_loc_key
        self.default = default

    def localize(self, new_loc_key: str = "", new_tooltip_key: str = "") -> None:
        if new_loc_key:
            self.label_loc_key = new_loc_key
        self.label_widget.config(text=self.loc_func(self.label_loc_key))

        if new_tooltip_key:
            self.tooltip_loc_key = new_tooltip_key
        self.loc_tooltip_var.set(self.loc_func(self.tooltip_loc_key))

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

    def disable(self) -> None:
        pass

    def enable(self) -> None:
        pass

    def inhibit(self) -> None:
        pass

    def disinhibit(self) -> None:
        pass


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
        entry_width: int = 20,
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

    def set(self, val: tuple[str, str]) -> None:
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


class Loc2Input(LocalizableWidget, Descriptive):
    def __init__(
        self,
        parent: ttk.Widget,
        loc_func: localize_function_type = placeholder_loc_func,
        font: Font | None = None,
        default: str = "",
        row: int = 0,
        col: int = 0,
        label_loc_key: str = "",
        desc_label_key: str | None = "",  # None can be specified to ignore some inputs.
        validation: str = "",
        label_width: int = 25,
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

        super().__init__(loc_func=loc_func, all_localized=all_localized, var=e, nominal_state="normal")
        lb = ttk.Label(parent, text=loc_func(label_loc_key), width=label_width, anchor=anchor)
        lb.grid(row=row, column=col + (1 if reverse else 0), sticky="nsew", padx=2, pady=2)

        self.loc_tooltip_var = StringVar(value=loc_func(tooltip_loc_key))
        create_tool_tip(lb, self.loc_tooltip_var, font=font)

        parent.rowconfigure(row, weight=0)

        en = ttk.Entry(
            parent,
            textvariable=e,
            validate="key",
            validatecommand=(validation, "%P") if validation else "",
            width=entry_width,
            foreground=color,
            justify="center",
        )
        en.grid(row=row, column=col + (0 if reverse else 1), sticky="nsew", padx=2, pady=2)
        en.bind("<FocusOut>", lambda event: formatter(event, e))

        self.default = default
        self.label_widget = lb
        self.input_widget = en

        self.row = row
        self.label_loc_key = label_loc_key
        self.desc_label_key = desc_label_key
        self.tooltip_loc_key = tooltip_loc_key

        self.dtype = dtype

    def localize(self, new_loc_key: str = "", new_tooltip_key: str = "") -> None:
        if new_loc_key:
            self.label_loc_key = new_loc_key
        self.label_widget.config(text=self.loc_func(self.label_loc_key))
        if new_tooltip_key:
            self.tooltip_loc_key = new_tooltip_key
        self.loc_tooltip_var.set(self.loc_func(self.tooltip_loc_key))

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

    def disable(self) -> None:
        self.input_widget.config(state="disabled")
        super().disable()

    def enable(self) -> None:
        self.input_widget.config(state=self.nominal_state)
        super().enable()

    def inhibit(self) -> None:
        self.input_widget.config(state="disabled")
        super().inhibit()

    def disinhibit(self) -> None:
        self.input_widget.config(state=self.target_state)
        super().disinhibit()

    def get_descriptive(self) -> str:
        if self.desc_label_key is None:  # None, explicit exclude
            return ""
        elif self.desc_label_key:  # str, include with label.
            return self.loc_func(self.desc_label_key, True)
        else:  # "", defers to label
            return self.loc_func(self.label_loc_key, True)


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
        validation: str = "",
        label_width: int = 25,
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
            validation=validation,
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
        if self.desc_label_key is None:  # None, explicit exclude
            return ""
        else:
            return super().get_descriptive() + (f" ({self.unit_text:})" if self.unit_text else "")


class LocDropdown(LocalizableWidget, Descriptive):
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
        super().__init__(loc_func=loc_func, all_localized=all_localized, var=StringVar(), nominal_state="readonly")
        self.nominal_state = "readonly"
        self.loc_func = loc_func
        self.desc_label_key = desc_label_key

        self.str_obj_dict: dict[str, object] = str_obj_dict if str_obj_dict else {"": ""}
        self.loc_str_obj_dict = self.update_loc_str_obj_dict()
        self.widget = ttk.Combobox(
            parent,
            textvariable=self.var,
            values=tuple(self.loc_str_obj_dict.keys()),
            justify="center",
            state=self.nominal_state,
        )

        self.tooltip_loc_key = tooltip_loc_key

        self.loc_tooltip_var = StringVar(value=loc_func(tooltip_loc_key))
        create_tool_tip(self.widget, self.loc_tooltip_var, font=font)

        self.widget.option_add("*TCombobox*Listbox.Justify", "center")
        self.widget.current(0)

    def localize(self) -> None:
        index = self.widget["values"].index(self.var.get())
        self.loc_str_obj_dict = self.update_loc_str_obj_dict()
        self.widget.config(values=tuple(self.loc_str_obj_dict.keys()))
        self.widget.current(index)
        self.loc_tooltip_var.set(self.loc_func(self.tooltip_loc_key))

    def update_loc_str_obj_dict(self) -> dict[str, object]:
        return {self.loc_func(k): v for k, v in self.str_obj_dict.items()}

    def get(self) -> str:
        return self.get_obj().__str__()

    def get_obj(self) -> object:
        return self.loc_str_obj_dict[self.var.get()]

    def set_by_str(self, string: str) -> None:
        """
        Given an unlocalized string / localization key, set the drop-down menu
        to the correct position.
        """
        self.widget.set(self.widget["values"][list(self.str_obj_dict.keys()).index(string)])

    def set_by_obj(self, obj: object) -> None:
        index = list(self.str_obj_dict.values()).index(obj)
        self.widget.current(index)

    def set(self, string: str) -> None:
        self.set_by_str(string)

    def grid(self, **kwargs: Any) -> None:
        self.widget.grid(**kwargs)

    def disable(self) -> None:
        self.widget.configure(state="disabled")
        super().disable()

    def enable(self) -> None:
        self.widget.configure(state=self.nominal_state)
        super().enable()

    def inhibit(self) -> None:
        self.widget.config(state="disabled")
        super().inhibit()

    def disinhibit(self) -> None:
        self.widget.config(state=self.target_state)
        super().disinhibit()

    def get_descriptive(self) -> str:
        return self.loc_func(self.desc_label_key, True)

    def reset(self, str_obj_dict: dict[str, object] | None = None, overwrite: bool = True) -> None:
        if str_obj_dict is not None:
            self.str_obj_dict = str_obj_dict
            self.loc_str_obj_dict = self.update_loc_str_obj_dict()
        self.widget["values"] = tuple(self.loc_str_obj_dict.keys())

        if overwrite:
            self.widget.current(0)
        else:

            if self.widget.get() not in self.widget["values"]:
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

        self.loc_tooltip_var = StringVar(value=loc_func(tooltip_loc_key))
        create_tool_tip(self, self.loc_tooltip_var, font=font)

        self.tooltip_loc_key = tooltip_loc_key
        if isinstance(all_localized, list):
            all_localized.append(self)

    def localize(self) -> None:
        self.config(text=self.loc_func(self.loc_key))
        self.loc_tooltip_var.set(self.loc_func(self.tooltip_loc_key))


class LocLabelCheck(LocalizableWidget, Descriptive):
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
        super().__init__(
            loc_func=loc_func, all_localized=all_localized, var=BooleanVar(value=default), nominal_state="normal"
        )
        self.default = default

        self.loc_func = loc_func
        self.check_widget = ttk.Checkbutton(parent, text=loc_func(label_loc_key), variable=self.var, width=width)
        if not skip_grid:
            self.check_widget.grid(row=row, column=col, sticky="nsew", columnspan=columnspan, padx=2, pady=2)

        self.loc_tooltip_var = StringVar(value=loc_func(tooltip_loc_key))
        create_tool_tip(self.check_widget, self.loc_tooltip_var, font=font)

        self.label_loc_key = label_loc_key
        self.desc_label_key = desc_label_key
        self.tooltip_loc_key = tooltip_loc_key

    def localize(self, new_loc_key: str = "") -> None:
        if new_loc_key:
            self.label_loc_key = new_loc_key
        self.check_widget.config(text=self.loc_func(self.label_loc_key))
        self.loc_tooltip_var.set(self.loc_func(self.tooltip_loc_key))

    def remove(self) -> None:
        self.check_widget.grid_remove()

    def restore(self) -> None:
        self.check_widget.grid()

    def set(self, value: bool) -> None:
        self.var.set(value)

    def trace_add(self, *args: Any) -> None:
        self.var.trace_add(*args)

    def disable(self):
        self.check_widget.config(state="disabled")
        super().disable()

    def enable(self) -> None:
        self.check_widget.config(state=self.nominal_state)
        super().enable()

    def inhibit(self) -> None:
        self.check_widget.config(state="disabled")
        super().inhibit()

    def disinhibit(self) -> None:
        self.check_widget.config(state=self.target_state)
        super().disinhibit()

    def get_descriptive(self) -> str:
        if self.desc_label_key is None:  # None, explicit exclude
            return ""
        elif self.desc_label_key:  # str, include with label.
            return self.loc_func(self.desc_label_key, True)
        else:  # "", defers to label
            return self.loc_func(self.label_loc_key, True)

    def reset(self) -> None:
        self.var.set(self.default)


class LocalizedFrame(Frame):
    def __init__(
        self,
        master: Tk | Toplevel | LocalizedFrame,
        *args: Any,
        font: Font | None = None,
        localization_dict: dict[str, dict[str, str]],
        menubar: Menu,
        default_lang: str,
        lang_var: StringVar | None = None,
        **kwargs: Any,
    ) -> None:
        super().__init__(master, *args, **kwargs)
        self.localization_dict = localization_dict

        self.localized_widgets: list[LocalizableWidget] = []
        self.font = font
        self.lang_var: StringVar

        if isinstance(master, (Toplevel, Tk)):  # if this LocalizedFrame acts as the top-level frame:

            self.lang_var = StringVar(
                value=(
                    default_lang
                    if default_lang in self.localization_dict.keys()
                    else list(self.localization_dict.keys())[0]
                )
            )

            lang_menu = Menu(menubar)
            menubar.add_cascade(label="Lang 语言", menu=lang_menu)

            for lang in localization_dict.keys():
                lang_menu.add_radiobutton(label=lang, variable=self.lang_var, value=lang, command=self.change_lang)
        elif isinstance(master, LocalizedFrame):  # otherwise, piggyback on to the parent's top-level frame
            self.lang_var = master.lang_var
        else:  # otherwise, let the children decide how to handle this.
            if lang_var is None:
                raise ValueError("lang_var must be provided when master is not a Toplevel, Tk, or LocalizedFrame")
            self.lang_var = lang_var

    def change_lang(self, lang_var: tkinter.StringVar | None = None) -> None:
        for loc_widget in self.localized_widgets:
            loc_widget.localize()

    @warn
    @wraps(Loc2LineDisp.__init__)
    def add_localized_2_line_display(self, *args: Any, **kwargs: Any) -> Loc2LineDisp:
        return Loc2LineDisp(
            *args, loc_func=self.get_loc_str, all_localized=self.localized_widgets, font=self.font, **kwargs
        )

    @warn
    @wraps(Loc3LineDisp.__init__)
    def add_localized_3_line_display(self, *args: Any, **kwargs: Any) -> Loc3LineDisp:
        return Loc3LineDisp(
            *args, loc_func=self.get_loc_str, all_localized=self.localized_widgets, font=self.font, **kwargs
        )

    @warn
    @wraps(LocLabelCheck.__init__)
    def add_localized_label_check(self, *args: Any, **kwargs: Any) -> LocLabelCheck:
        return LocLabelCheck(
            *args, loc_func=self.get_loc_str, all_localized=self.localized_widgets, font=self.font, **kwargs
        )

    @warn
    @wraps(Loc2Input.__init__)
    def add_localized_2_input(self, *args: Any, **kwargs: Any) -> Loc2Input:
        return Loc2Input(
            *args, loc_func=self.get_loc_str, all_localized=self.localized_widgets, font=self.font, **kwargs
        )

    @warn
    @wraps(Loc3Input.__init__)
    def add_localized_3_input(self, *args: Any, **kwargs: Any) -> Loc3Input:
        return Loc3Input(
            *args, loc_func=self.get_loc_str, all_localized=self.localized_widgets, font=self.font, **kwargs
        )

    @warn
    @wraps(LocDropdown.__init__)
    def add_localized_dropdown(self, *args: Any, **kwargs: Any) -> LocDropdown:
        return LocDropdown(
            *args, loc_func=self.get_loc_str, all_localized=self.localized_widgets, font=self.font, **kwargs
        )

    @warn
    @wraps(LocLabelFrame.__init__)
    def add_localized_label_frame(self, *args: Any, **kwargs: Any) -> LocLabelFrame:
        return LocLabelFrame(
            *args, loc_func=self.get_loc_str, all_localized=self.localized_widgets, font=self.font, **kwargs
        )

    def get_loc_str(self, name: str, force_default: bool = False) -> str:
        try:
            lang = "English" if force_default else self.lang_var.get()  # type: ignore[union-attr]
            return self.localization_dict[lang][name]
        except KeyError:
            try:
                return self.localization_dict["English"][name]
            except KeyError:
                return name
