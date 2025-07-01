from tkinter import IntVar, StringVar, ttk
from typing import Literal

from .misc import formatFloatInput
from .tip import CreateToolTip


class Loc12Disp:
    def __init__(
        self,
        parent,
        row=0,
        col=0,
        labelLocKey="",
        unitText="",
        default="",
        entryWidth=20,
        justify: Literal["left", "right", "center"] = "center",
        tooltipLocKey=None,
        reverse=False,
        locFunc=None,
        allDisps=None,
    ):
        lb = ttk.Label(parent, text=locFunc(labelLocKey))
        lb.grid(row=row, column=col, columnspan=2, sticky="nsew", padx=2, pady=2)
        e = StringVar(parent)
        e.default = default
        e.set(default)
        parent.rowconfigure(row, weight=0)
        en = ttk.Entry(parent, textvariable=e, width=entryWidth, state="disabled", justify=justify)
        en.grid(row=row + 1, column=col + (1 if reverse else 0), sticky="nsew", padx=2, pady=2)
        ttk.Label(parent, text=unitText).grid(
            row=row + 1, column=col + (0 if reverse else 1), sticky="nsew", padx=2, pady=2
        )
        if tooltipLocKey is not None:
            self.locTooltipVar = StringVar(value=locFunc(tooltipLocKey))
            CreateToolTip(lb, self.locTooltipVar)
        else:
            self.locTooltipVar = None

        self.locFunc = locFunc
        self.labelLocKey = labelLocKey
        self.labelWidget = lb

        self.tooltipLocKey = tooltipLocKey

        self.entryVar = e
        self.entryWidget = en

        if allDisps is not None:
            allDisps.append(self)

    def reLocalize(self, newLocKey=None, newTooltipKey=None):
        if newLocKey is not None:
            self.labelLocKey = newLocKey
        self.labelWidget.config(text=self.locFunc(self.labelLocKey))

        if self.locTooltipVar is not None:
            if newTooltipKey is not None:
                self.tooltipLocKey = newTooltipKey
            self.locTooltipVar.set(self.locFunc(self.tooltipLocKey))

    def set(self, val):
        self.entryVar.set(val)


class Loc122Disp(Loc12Disp):
    def __init__(
        self,
        parent,
        row=0,
        col=0,
        labelLocKey="",
        unitText_up="",
        unitText_dn="",
        default_up="",
        default_dn="",
        entryWidth=20,
        justify_up: Literal["left", "right", "center"] = "center",
        justify_dn: Literal["left", "right", "center"] = "center",
        tooltipLocKey=None,
        reverse=False,
        locFunc=None,
        allDisps=None,
    ):
        super().__init__(
            parent=parent,
            row=row,
            col=col,
            labelLocKey=labelLocKey,
            unitText=unitText_up,
            default=default_up,
            entryWidth=entryWidth,
            justify=justify_up,
            tooltipLocKey=tooltipLocKey,
            reverse=reverse,
            locFunc=locFunc,
        )
        e2 = StringVar(parent)
        e2.default = default_dn
        e2.set(default_dn)
        parent.rowconfigure(row, weight=0)
        en2 = ttk.Entry(parent, textvariable=e2, width=entryWidth, state="disabled", justify=justify_dn)
        en2.grid(row=row + 2, column=col + (1 if reverse else 0), sticky="nsew", padx=2, pady=2)
        ttk.Label(parent, text=unitText_dn).grid(
            row=row + 2, column=col + (0 if reverse else 1), sticky="nsew", padx=2, pady=2
        )

        self.auxEntryVar = e2
        self.auxEntryWidget = en2

        if allDisps is not None:
            allDisps.append(self)

    def set(self, val):
        val_1, val_2 = val
        self.entryVar.set(val_1)
        self.auxEntryVar.set(val_2)


class Loc2Input:
    def __init__(
        self,
        parent,
        row=0,
        col=0,
        labelLocKey="",
        descLabelKey=None,
        default="",
        validation=None,
        labelWidth=20,
        entryWidth=10,
        formatter=formatFloatInput,
        color=None,
        tooltipLocKey=None,
        anchor="w",
        reverse=False,
        locFunc=None,
        allInputs=None,
    ):
        # noinspection PyTypeChecker
        lb = ttk.Label(parent, text=locFunc(labelLocKey), width=labelWidth, anchor=anchor)

        lb.grid(
            row=row,
            column=col + (1 if reverse else 0),
            sticky="nsew",
            padx=2,
            pady=2,
        )
        if tooltipLocKey is not None:
            self.locTooltipVar = StringVar(parent, value=locFunc(tooltipLocKey))
            CreateToolTip(lb, self.locTooltipVar)
        else:
            self.locTooltipVar = None

        parent.rowconfigure(row, weight=0)
        e = StringVar(parent)
        e.set(default)
        en = ttk.Entry(
            parent,
            textvariable=e,
            validate="key",
            validatecommand=(validation, "%P"),
            width=entryWidth,
            foreground=color,
            justify="center",
        )
        en.default = default
        en.grid(
            row=row,
            column=col + (0 if reverse else 1),
            sticky="nsew",
            padx=2,
            pady=2,
        )
        en.bind("<FocusOut>", lambda event: formatter(event, e))

        self.labelWidget = lb
        self.inputVar = e
        self.inputWidget = en

        self.row = row
        self.labelLocKey = labelLocKey
        self.descLabelKey = descLabelKey
        self.tooltipLocKey = tooltipLocKey
        self.locFunc = locFunc
        self.nominalState = "normal"

        if allInputs is not None:
            allInputs.append(self)

    def reLocalize(self, newLocKey=None, newTooltipKey=None):
        if newLocKey is not None:
            self.labelLocKey = newLocKey
        self.labelWidget.config(text=self.locFunc(self.labelLocKey))

        if self.locTooltipVar is not None:
            if newTooltipKey is not None:
                self.tooltipLocKey = newTooltipKey
            self.locTooltipVar.set(self.locFunc(self.tooltipLocKey))

    def remove(self):
        self.labelWidget.grid_remove()
        self.inputWidget.grid_remove()

    def restore(self):
        self.labelWidget.grid()
        self.inputWidget.grid()

    def get(self):
        return self.inputVar.get()

    def set(self, val):
        self.inputVar.set(val)

    def trace_add(self, *args):
        self.inputVar.trace_add(*args)

    def disable(self):
        self.inputWidget.config(state="disabled")

    def enable(self):
        self.inputWidget.config(state="normal")

    def inhibit(self):
        self.nominalState = self.inputWidget.cget("state")
        self.inputWidget.config(state="disabled")

    def disinhibit(self):
        self.inputWidget.config(state=self.nominalState)

    def getDescriptive(self):
        return self.locFunc(
            self.labelLocKey if self.descLabelKey is None else self.descLabelKey,
            forceDefault=True,
        )


class Loc3Input(Loc2Input):
    def __init__(
        self,
        parent,
        row=0,
        col=0,
        labelLocKey="",
        descLabelKey=None,
        unitText="",
        default="",
        validation=None,
        labelWidth=20,
        entryWidth=10,
        formatter=formatFloatInput,
        color=None,
        tooltipLocKey=None,
        anchor="w",
        reverse=False,
        locFunc=None,
        allInputs=None,
    ):
        super().__init__(
            parent=parent,
            row=row,
            col=col,
            labelLocKey=labelLocKey,
            descLabelKey=descLabelKey,
            default=default,
            validation=validation,
            labelWidth=labelWidth,
            entryWidth=entryWidth,
            formatter=formatter,
            color=color,
            tooltipLocKey=tooltipLocKey,
            anchor=anchor,
            reverse=reverse,
            locFunc=locFunc,
        )

        ulb = ttk.Label(parent, text=unitText)
        ulb.grid(row=row, column=col + 2, sticky="nsew", padx=2, pady=2)
        self.unitText = unitText
        self.unitLabel = ulb

        if allInputs is not None:
            allInputs.append(self)

    def remove(self):
        super().remove()
        self.unitLabel.grid_remove()

    def restore(self):
        super().restore()
        self.unitLabel.grid()

    def getDescriptive(self):
        return self.locFunc(
            self.labelLocKey if self.descLabelKey is None else self.descLabelKey,
            forceDefault=True,
        ) + (f" ({self.unitText:})" if self.unitText != "" else "")


class LocDropdown:
    def __init__(self, parent, strObjDict, locFunc, dropdowns=None, descLabelKey=""):
        """
        localized key of string type: underlying object
        """
        self.nominalState = "readonly"
        self.textVar = StringVar(parent)

        self.locFunc = locFunc
        self.descLabelKey = descLabelKey

        self.strObjDict = strObjDict
        self.locStrObjDict = {self.locFunc(k): v for k, v in strObjDict.items()}

        self.widget = ttk.Combobox(
            parent,
            textvariable=self.textVar,
            values=tuple(self.locStrObjDict.keys()),
            justify="center",
            state="readonly",
        )
        self.widget.option_add("*TCombobox*Listbox.Justify", "center")
        self.widget.current(0)

        if dropdowns is not None:
            dropdowns.append(self)

    def reLocalize(self):
        index = self.widget["values"].index(self.textVar.get())
        self.locStrObjDict = {self.locFunc(k): v for k, v in self.strObjDict.items()}
        self.widget.config(values=tuple(self.locStrObjDict.keys()))
        self.widget.current(index)

    def get(self):
        return self.getObj()

    def getObj(self):
        return self.locStrObjDict[self.textVar.get()]

    def setByStr(self, string):
        """
        Given an unlocalized string / localization key, set the drop-down menu
        to the correct position.
        """
        self.widget.set(self.widget["values"][list(self.strObjDict.keys()).index(string)])

    def setByObj(self, obj):
        index = list(self.strObjDict.values()).index(obj)
        self.widget.current(index)

    def set(self, str):
        self.setByStr(str)

    def grid(self, **kwargs):
        self.widget.grid(**kwargs)

    def trace_add(self, *args):
        self.textVar.trace_add(*args)

    def disable(self):
        self.widget.configure(state="disabled")

    def enable(self):
        self.widget.configure(state="readonly")

    def inhibit(self):
        self.nominalState = self.widget.cget("state")
        self.widget.config(state="disabled")

    def disinhibit(self):
        self.widget.config(state=self.nominalState)

    def getDescriptive(self):
        return self.locFunc(self.descLabelKey, forceDefault=True)

    def reset(self, strObjDict):
        currentObject = self.getObj()
        self.strObjDict = strObjDict
        self.locStrObjDict = {self.locFunc(k): v for k, v in strObjDict.items()}
        self.widget["values"] = tuple(self.locStrObjDict.keys())
        try:
            self.setByObj(currentObject)
        except ValueError:
            self.widget.current(0)


class LocLabelFrame(ttk.LabelFrame):
    def __init__(self, *args, locKey="", tooltipLocKey=None, locFunc=None, allLLF=None, **kwargs):
        self.locKey = locKey
        self.locFunc = locFunc
        super().__init__(*args, text=locFunc(locKey), **kwargs)

        if tooltipLocKey is not None:
            self.locTooltipVar = StringVar(value=locFunc(tooltipLocKey))
            CreateToolTip(self, self.locTooltipVar)
        else:
            self.locTooltipVar = None

        self.tooltipLocKey = tooltipLocKey

        if allLLF is not None:
            allLLF.append(self)

    def reLocalize(self):
        self.config(text=self.locFunc(self.locKey))
        if self.locTooltipVar is not None:
            self.locTooltipVar.set(self.locFunc(self.tooltipLocKey))


class LocLabelCheck:
    def __init__(
        self,
        parent,
        default=1,
        row=0,
        col=0,
        columnspan=None,
        labelLocKey="",
        descLabelKey=None,
        tooltipLocKey=None,
        locFunc=None,
        width=25,
        allLC=None,
    ):
        self.nominalState = "normal"
        self.checkVar = IntVar(value=default)
        self.locFunc = locFunc
        self.checkWidget = ttk.Checkbutton(
            parent,
            text=locFunc(labelLocKey),
            variable=self.checkVar,
            width=width,
        )
        self.checkWidget.grid(
            row=row,
            column=col,
            sticky="nsew",
            columnspan=columnspan,
            padx=2,
            pady=2,
        )

        if tooltipLocKey is not None:
            self.locTooltipVar = StringVar(value=locFunc(tooltipLocKey))
            CreateToolTip(self.checkWidget, self.locTooltipVar)
        else:
            self.locTooltipVar = None

        self.labelLocKey = labelLocKey
        self.descLabelKey = descLabelKey
        self.tooltipLocKey = tooltipLocKey

        if allLC is not None:
            allLC.append(self)

    def reLocalize(self, newLocKey=None):
        if newLocKey is not None:
            self.labelLocKey = newLocKey
        self.checkWidget.config(text=self.locFunc(self.labelLocKey))
        if self.locTooltipVar is not None:
            self.locTooltipVar.set(self.locFunc(self.tooltipLocKey))

    def disable(self):
        self.checkWidget.config(state="disabled")

    def enable(self):
        self.checkWidget.config(state="normal")

    def remove(self):
        self.checkWidget.grid_remove()

    def restore(self):
        self.checkWidget.grid()

    def get(self):
        return self.checkVar.get()

    def set(self, value):
        self.checkVar.set(value)

    def trace_add(self, *args):
        self.checkVar.trace_add(*args)

    def inhibit(self):
        self.nominalState = self.checkWidget.cget("state")
        self.checkWidget.config(state="disabled")

    def disinhibit(self):
        self.checkWidget.config(state=self.nominalState)

    def getDescriptive(self):
        return self.locFunc(self.labelLocKey if self.descLabelKey is None else self.descLabelKey, forceDefault=True)
