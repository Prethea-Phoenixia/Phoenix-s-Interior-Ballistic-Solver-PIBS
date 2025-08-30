from __future__ import annotations

import warnings
from abc import ABC, abstractmethod
from functools import wraps
from tkinter import BooleanVar, Event, Frame, Menu, StringVar, ttk
from tkinter.font import Font
from typing import Any, Callable, Literal

from .misc import format_float_input
from .tip import create_tool_tip

localize_function_type = Callable[[str, bool], str] | Callable[[str], str]
placeholder_loc_func = lambda _: ""


class Localizable(ABC):
    @abstractmethod
    def localize(self, *args: Any, **kwargs: Any) -> None: ...


class LocalizableWidget(Localizable):
    def __init__(
        self,
        *args,
        var: Any,
        loc_func: localize_function_type,
        all_localized: list[LocalizableWidget] | None = None,
        **kwargs,
    ):
        self.loc_func = loc_func
        if isinstance(all_localized, list):
            all_localized.append(self)

        self.var = var

    @abstractmethod
    def localize(self, *args: Any, **kwargs: Any) -> None: ...


class Descriptive(ABC):
    @abstractmethod
    def get_descriptive(self) -> str:
        """
        returns the json key to this widget's value.
        """
        ...

    @abstractmethod
    def get(self) -> Any:
        """
        returns the object containing the data this widget holds.
        """
        ...

    @abstractmethod
    def reset(self, *_) -> None:
        """
        resets the value of this widget to the default value.
        """
        ...


class Loc12Disp(LocalizableWidget):
    def __init__(
        self,
        parent,
        loc_func: localize_function_type = placeholder_loc_func,
        font: None | Font = None,
        row: int = 0,
        col: int = 0,
        label_loc_key: str = "",
        unit_text: str = "",
        default: str = "",
        entry_width: int = 20,
        justify: Literal["left", "right", "center"] = "center",
        tooltip_loc_key: str = "",
        reverse: bool = False,
        all_localized: list | None = None,
    ):

        e = StringVar()
        e.set(default)
        super().__init__(loc_func=loc_func, all_localized=all_localized, var=e)
        lb = ttk.Label(parent, text=loc_func(label_loc_key))
        lb.grid(row=row, column=col, columnspan=2, sticky="nsew", padx=2, pady=2)

        parent.rowconfigure(row, weight=0)
        en = ttk.Entry(parent, textvariable=e, width=entry_width, state="disabled", justify=justify)
        # noinspection PyUnresolvedReferences
        en.grid(row=row + 1, column=col + (1 if reverse else 0), sticky="nsew", padx=2, pady=2)
        # noinspection PyUnresolvedReferences
        ttk.Label(parent, text=unit_text).grid(
            row=row + 1, column=col + (0 if reverse else 1), sticky="nsew", padx=2, pady=2
        )
        self.loc_tooltip_var = StringVar(value=loc_func(tooltip_loc_key))
        create_tool_tip(lb, self.loc_tooltip_var, font=font)

        self.label_loc_key = label_loc_key
        self.label_widget = lb

        self.tooltip_loc_key = tooltip_loc_key
        self.entry_widget = en

        self.default = default

    def localize(self, new_loc_key: str = "", new_tooltip_key: str = "") -> None:
        if new_loc_key:
            self.label_loc_key = new_loc_key
        self.label_widget.config(text=self.loc_func(self.label_loc_key))

        if new_tooltip_key:
            self.tooltip_loc_key = new_tooltip_key
        self.loc_tooltip_var.set(self.loc_func(self.tooltip_loc_key))

    def set(self, val) -> None:
        self.var.set(val)

    def reset(self) -> None:
        self.var.set(self.default)


class Loc122Disp(Loc12Disp):
    def __init__(
        self,
        parent,
        loc_func: localize_function_type = placeholder_loc_func,
        font: Font | None = None,
        row: int = 0,
        col: int = 0,
        label_loc_key: str = "",
        unit_text_up: str = "",
        unit_text_dn: str = "",
        default_up: str = "",
        default_dn: str = "",
        entry_width: int = 20,
        justify_up: Literal["left", "right", "center"] = "center",
        justify_dn: Literal["left", "right", "center"] = "center",
        tooltip_loc_key: str = "",
        reverse: bool = False,
        all_localized: None | list = None,
    ):
        super().__init__(
            parent=parent,
            row=row,
            col=col,
            font=font,
            label_loc_key=label_loc_key,
            unit_text=unit_text_up,
            default=default_up,
            entry_width=entry_width,
            justify=justify_up,
            tooltip_loc_key=tooltip_loc_key,
            reverse=reverse,
            loc_func=loc_func,
            all_localized=all_localized,
        )
        e2 = StringVar()
        e2.set(default_dn)
        parent.rowconfigure(row, weight=0)
        en2 = ttk.Entry(parent, textvariable=e2, width=entry_width, state="disabled", justify=justify_dn)
        en2.grid(row=row + 2, column=col + (1 if reverse else 0), sticky="nsew", padx=2, pady=2)
        ttk.Label(parent, text=unit_text_dn).grid(
            row=row + 2, column=col + (0 if reverse else 1), sticky="nsew", padx=2, pady=2
        )

        self.aux_entry_var = e2
        self.aux_entry_widget = en2
        self.aux_default = default_dn

    def set(self, val):
        val_1, val_2 = val
        self.var.set(val_1)
        self.aux_entry_var.set(val_2)

    def reset(self):
        self.aux_entry_var.set(self.aux_default)
        super().reset()


class Loc2Input(LocalizableWidget, Descriptive):
    def __init__(
        self,
        parent,
        loc_func: localize_function_type = placeholder_loc_func,
        font: None | Font = None,
        default: str = "",
        row: int = 0,
        col: int = 0,
        label_loc_key: str = "",
        desc_label_key: None | str = "",  # None can be specified to ignore some inputs.
        validation: str | None = None,
        label_width: int = 20,
        entry_width: int = 10,
        formatter: Callable[[Event, StringVar], None] = format_float_input,
        color: str = "",
        tooltip_loc_key: str = "",
        anchor: Literal["nw", "n", "ne", "w", "center", "e", "sw", "s", "se"] = "w",
        reverse: bool = False,
        all_localized: None | list = None,
        dtype: Callable[[str], any] = lambda v: int(float(v)),
    ):

        e = StringVar(value=default)

        super().__init__(loc_func=loc_func, all_localized=all_localized, var=e)
        lb = ttk.Label(parent, text=loc_func(label_loc_key), width=label_width, anchor=anchor)

        lb.grid(row=row, column=col + (1 if reverse else 0), sticky="nsew", padx=2, pady=2)

        self.loc_tooltip_var = StringVar(value=loc_func(tooltip_loc_key))
        create_tool_tip(lb, self.loc_tooltip_var, font=font)

        parent.rowconfigure(row, weight=0)

        # noinspection PyTypeChecker
        en = ttk.Entry(
            parent,
            textvariable=e,
            validate="key",
            validatecommand=(validation, "%P"),  # PyTypeChecker is wrong.
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
        self.nominalState = "normal"

        self.dtype = dtype

    def localize(self, new_loc_key: str = "", new_tooltip_key: str = ""):
        if new_loc_key:
            self.label_loc_key = new_loc_key
        self.label_widget.config(text=self.loc_func(self.label_loc_key))
        if new_tooltip_key:
            self.tooltip_loc_key = new_tooltip_key
        self.loc_tooltip_var.set(self.loc_func(self.tooltip_loc_key))

    def remove(self):
        self.label_widget.grid_remove()
        self.input_widget.grid_remove()

    def restore(self):
        self.label_widget.grid()
        self.input_widget.grid()

    def reset(self, *_) -> None:
        self.var.set(self.default)

    def get(self) -> Any:
        return self.dtype(self.var.get())

    def set(self, val) -> None:
        self.var.set(val)

    def disable(self) -> None:
        self.input_widget.config(state="disabled")

    def enable(self) -> None:
        self.input_widget.config(state="normal")

    def inhibit(self) -> None:
        self.nominalState = self.input_widget.cget("state")
        self.input_widget.config(state="disabled")

    def disinhibit(self) -> None:
        self.input_widget.config(state=self.nominalState)

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
        parent,
        loc_func: localize_function_type = placeholder_loc_func,
        font: None | Font = None,
        row: int = 0,
        col: int = 0,
        label_loc_key: str = "",
        desc_label_key: None | str = "",
        unit_text: str = "",
        default: str = "",
        validation: str | None = None,
        label_width: int = 20,
        entry_width: int = 10,
        formatter: Callable[[Event, StringVar], None] = format_float_input,
        color: None | str = None,
        tooltip_loc_key: str = "",
        anchor: Literal["w", "e", "s", "n"] = "w",
        reverse: bool = False,
        all_localized: None | list = None,
        dtype: Callable[[str], any] = lambda v: int(float(v)),
    ):
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

    def remove(self):
        super().remove()
        self.unit_label.grid_remove()

    def restore(self):
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
        parent,
        loc_func: localize_function_type = placeholder_loc_func,
        font: None | Font = None,
        str_obj_dict: dict[str, object] | None = None,
        all_localized: list | None = None,
        desc_label_key: str = "",
        tooltip_loc_key: str = "",
    ):
        """
        localized key of string type: underlying object
        """

        str_obj_dict = str_obj_dict if str_obj_dict else {"": ""}

        super().__init__(loc_func=loc_func, all_localized=all_localized, var=StringVar())
        self.nominal_state = "readonly"
        self.loc_func = loc_func
        self.desc_label_key = desc_label_key

        self.str_obj_dict = str_obj_dict
        self.loc_str_obj_dict = {self.loc_func(k): v for k, v in str_obj_dict.items()}

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

    def localize(self):
        index = self.widget["values"].index(self.var.get())
        self.loc_str_obj_dict = {self.loc_func(k): v for k, v in self.str_obj_dict.items()}
        self.widget.config(values=tuple(self.loc_str_obj_dict.keys()))
        self.widget.current(index)
        self.loc_tooltip_var.set(self.loc_func(self.tooltip_loc_key))

    def get(self) -> str:
        return str(self.get_obj())

    def get_obj(self):
        return self.loc_str_obj_dict[self.var.get()]

    def set_by_str(self, string):
        """
        Given an unlocalized string / localization key, set the drop-down menu
        to the correct position.
        """
        self.widget.set(self.widget["values"][list(self.str_obj_dict.keys()).index(string)])

    def set_by_obj(self, obj):
        index = list(self.str_obj_dict.values()).index(obj)
        self.widget.current(index)

    def set(self, string: str):
        self.set_by_str(string)

    def grid(self, **kwargs):
        self.widget.grid(**kwargs)

    def trace_add(self, *args):
        self.var.trace_add(*args)

    def disable(self):
        self.widget.configure(state="disabled")

    def enable(self):
        self.widget.configure(state="readonly")

    def inhibit(self):
        self.nominal_state = self.widget.cget("state")
        self.widget.config(state="disabled")

    def disinhibit(self):
        self.widget.config(state=self.nominal_state)

    def get_descriptive(self) -> str:
        return self.loc_func(self.desc_label_key, True)

    def reset(self, str_obj_dict=None):
        if str_obj_dict:
            self.str_obj_dict = str_obj_dict
        self.loc_str_obj_dict = {self.loc_func(k): v for k, v in self.str_obj_dict.items()}
        self.widget["values"] = tuple(self.loc_str_obj_dict.keys())
        self.widget.current(0)


class LocLabelFrame(ttk.LabelFrame, Localizable):
    def __init__(
        self,
        *args,
        loc_func: localize_function_type = placeholder_loc_func,
        font: Font | None = None,
        label_loc_key: str = "",
        tooltip_loc_key: str = "",
        all_localized: list | None = None,
        **kwargs,
    ):
        super().__init__(*args, text=loc_func(label_loc_key), **kwargs)
        self.loc_key = label_loc_key
        self.loc_func = loc_func

        self.loc_tooltip_var = StringVar(value=loc_func(tooltip_loc_key))
        create_tool_tip(self, self.loc_tooltip_var, font=font)

        self.tooltip_loc_key = tooltip_loc_key
        if isinstance(all_localized, list):
            all_localized.append(self)

    def localize(self):
        self.config(text=self.loc_func(self.loc_key))
        self.loc_tooltip_var.set(self.loc_func(self.tooltip_loc_key))


class LocLabelCheck(LocalizableWidget, Descriptive):
    def __init__(
        self,
        parent,
        loc_func: localize_function_type = placeholder_loc_func,
        font: Font | None = None,
        default: bool = True,
        row: int = 0,
        col: int = 0,
        columnspan: int = 1,
        label_loc_key: str = "",
        desc_label_key: str | None = "",
        tooltip_loc_key: str = "",
        width: int = 25,
        all_localized: list | None = None,
    ):
        super().__init__(loc_func=loc_func, all_localized=all_localized, var=BooleanVar(value=default))
        self.default = default
        self.nominal_state = "normal"
        self.loc_func = loc_func
        self.check_widget = ttk.Checkbutton(parent, text=loc_func(label_loc_key), variable=self.var, width=width)
        self.check_widget.grid(row=row, column=col, sticky="nsew", columnspan=columnspan, padx=2, pady=2)

        self.loc_tooltip_var = StringVar(value=loc_func(tooltip_loc_key))
        create_tool_tip(self.check_widget, self.loc_tooltip_var, font=font)

        self.label_loc_key = label_loc_key
        self.desc_label_key = desc_label_key
        self.tooltip_loc_key = tooltip_loc_key

    def localize(self, new_loc_key: str = ""):
        if new_loc_key:
            self.label_loc_key = new_loc_key
        self.check_widget.config(text=self.loc_func(self.label_loc_key))
        self.loc_tooltip_var.set(self.loc_func(self.tooltip_loc_key))

    def disable(self):
        self.check_widget.config(state="disabled")

    def enable(self):
        self.check_widget.config(state="normal")

    def remove(self):
        self.check_widget.grid_remove()

    def restore(self):
        self.check_widget.grid()

    def get(self) -> bool:
        return bool(self.var.get())

    def set(self, value):
        self.var.set(value)

    def trace_add(self, *args):
        self.var.trace_add(*args)

    def inhibit(self):
        self.nominal_state = self.check_widget.cget("state")
        self.check_widget.config(state="disabled")

    def disinhibit(self):
        self.check_widget.config(state=self.nominal_state)

    def get_descriptive(self) -> str:
        if self.desc_label_key is None:  # None, explicit exclude
            return ""
        elif self.desc_label_key:  # str, include with label.
            return self.loc_func(self.desc_label_key, True)
        else:  # "", defers to label
            return self.loc_func(self.label_loc_key, True)

    def reset(self) -> None:
        self.var.set(self.default)


def warn(func):

    def wrapped(*args, font=None, loc_func=None, all_localized=None, **kwargs):
        if loc_func or all_localized or font:
            warnings.warn(
                "LocalizedFrame sets the font, loc_fun and all_localized parameters to instance of itself, supplied arguments will be ignored."
            )
        return func(*args, **kwargs)

    return wrapped


class LocalizedFrame(Frame):
    def __init__(
        self,
        *args,
        font: Font | None = None,
        localization_dict: dict[str, dict[str, str]],
        menubar: Menu,
        default_lang: str,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.localization_dict = localization_dict
        self.langVar = StringVar(
            value=(
                default_lang
                if default_lang in self.localization_dict.keys()
                else list(self.localization_dict.keys())[0]
            )
        )
        self.locs = []
        self.font = font

        lang_menu = Menu(menubar)
        menubar.add_cascade(
            label="Lang 语言",
            menu=lang_menu,
        )

        for lang in localization_dict.keys():
            lang_menu.add_radiobutton(label=lang, variable=self.langVar, value=lang, command=self.change_lang)

    def change_lang(self):
        pass

    @warn
    @wraps(Loc12Disp.__init__)
    def add_localized_12_display(self, *args, **kwargs) -> Loc12Disp:
        # logging.warn("")
        return Loc12Disp(*args, loc_func=self.get_loc_str, all_localized=self.locs, font=self.font, **kwargs)

    @warn
    @wraps(Loc122Disp.__init__)
    def add_localized_122_display(self, *args, **kwargs) -> Loc122Disp:
        return Loc122Disp(*args, loc_func=self.get_loc_str, all_localized=self.locs, font=self.font, **kwargs)

    @warn
    @wraps(LocLabelCheck.__init__)
    def add_localized_label_check(self, *args, **kwargs) -> LocLabelCheck:
        return LocLabelCheck(*args, loc_func=self.get_loc_str, all_localized=self.locs, font=self.font, **kwargs)

    @warn
    @wraps(Loc2Input.__init__)
    def add_localized_2_input(self, *args, **kwargs) -> Loc2Input:
        return Loc2Input(*args, loc_func=self.get_loc_str, all_localized=self.locs, font=self.font, **kwargs)

    @warn
    @wraps(Loc3Input.__init__)
    def add_localized_3_input(self, *args, **kwargs) -> Loc3Input:
        return Loc3Input(*args, loc_func=self.get_loc_str, all_localized=self.locs, font=self.font, **kwargs)

    @warn
    @wraps(LocDropdown.__init__)
    def add_localized_dropdown(self, *args, **kwargs) -> LocDropdown:
        return LocDropdown(*args, loc_func=self.get_loc_str, all_localized=self.locs, font=self.font, **kwargs)

    @warn
    @wraps(LocLabelFrame.__init__)
    def add_localized_label_frame(self, *args, **kwargs) -> LocLabelFrame:
        return LocLabelFrame(*args, loc_func=self.get_loc_str, all_localized=self.locs, font=self.font, **kwargs)

    def get_loc_str(self, name, force_default: bool = False):
        try:
            return self.localization_dict["English" if force_default else self.langVar.get()][name]
        except KeyError:
            try:
                return self.localization_dict["English"][name]
            except KeyError:
                return name
