from __future__ import annotations

import csv
from tkinter import ttk, Frame, filedialog, messagebox

from .ballistics import CONVENTIONAL, RECOILLESS
from .ballistics import (
    POINT_BURNOUT,
    POINT_EXIT,
    POINT_FRACTURE,
    POINT_PEAK_AVG,
    POINT_PEAK_BREECH,
    POINT_PEAK_SHOT,
    POINT_PEAK_STAG,
    POINT_START,
    COMPUTE,
)
from .ballistics.gun import GunResult
from .ballistics.recoilless import RecoillessResult
from .localized_widget import LocalizedFrame


class TableFrame(LocalizedFrame):
    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        table_label_frame = self.add_localized_label_frame(self, label_loc_key="tblFrmLabel")
        table_label_frame.grid(row=0, column=0, sticky="nsew")

        table_label_frame.rowconfigure(0, weight=1)
        table_label_frame.columnconfigure(0, weight=1)

        v_scroll = ttk.Scrollbar(table_label_frame, orient="vertical")  # create a scrollbar
        v_scroll.grid(row=0, rowspan=2, column=1, sticky="nsew")

        h_scroll = ttk.Scrollbar(table_label_frame, orient="horizontal")
        h_scroll.grid(row=1, column=0, sticky="nsew")

        # trick: intentionally mix .grid and .place geometry managers to prevent table from blowing up the entire frame on update.
        table_place_frame = Frame(table_label_frame)
        table_place_frame.grid(row=0, column=0, sticky="nsew")
        self.tv = ttk.Treeview(
            table_place_frame, selectmode="browse", yscrollcommand=v_scroll.set, xscrollcommand=h_scroll.set
        )
        self.tv.place(relwidth=1, relheight=1)

        v_scroll.configure(command=self.tv.yview)
        h_scroll.configure(command=self.tv.xview)

    def format_table(
        self, gun_result: GunResult | RecoillessResult, acc_exp: int
    ) -> list[
        tuple[str, str, str, str, str, str, str, str, str]
        | tuple[str, str, str, str, str, str, str, str, str, str, str, str]
    ]:
        table_data, use_sn, units = [], (), ()

        sigfig = acc_exp

        if isinstance(gun_result, GunResult):
            table_data = [
                (
                    self.get_loc_str(entry.tag),
                    f"{entry.time * 1e3:.{sigfig}g}",
                    f"{entry.travel:.{sigfig}g}",
                    f"{entry.burnup:.0%}",
                    f"{entry.velocity:.{sigfig}g}",
                    f"{entry.breech_pressure * 1e-6:.{sigfig}g}",
                    f"{entry.avg_pressure * 1e-6:.{sigfig}g}",
                    f"{entry.shot_pressure * 1e-6:.{sigfig}g}",
                    f"{entry.temperature:.{sigfig}g}",
                )
                for entry in gun_result.table_data
            ]

        elif isinstance(gun_result, RecoillessResult):
            table_data = [
                (
                    self.get_loc_str(entry.tag),
                    f"{entry.time * 1e3:.{sigfig}g}",
                    f"{entry.travel:.{sigfig}g}",
                    f"{entry.burnup:.0%}",
                    f"{entry.outflow_fraction:.0%}",
                    f"{entry.velocity:.{sigfig}g}",
                    f"{entry.outflow_velocity:.{sigfig}g}",
                    f"{entry.breech_pressure * 1e-6:.{sigfig}g}",
                    f"{entry.stag_pressure * 1e-6:.{sigfig}g}",
                    f"{entry.avg_pressure * 1e-6:.{sigfig}g}",
                    f"{entry.shot_pressure * 1e-6:.{sigfig}g}",
                    f"{entry.temperature:.{sigfig}g}",
                )
                for entry in gun_result.table_data
            ]

        return table_data

    def update_table(self, gun_result: GunResult | RecoillessResult, acc_exp: int):
        self.tv.delete(*self.tv.get_children())
        if not gun_result:
            return

        loc_table_data = self.format_table(gun_result=gun_result, acc_exp=acc_exp)

        if isinstance(gun_result, GunResult):
            column_list = self.get_loc_str("columnList")[CONVENTIONAL]
        elif isinstance(gun_result, RecoillessResult):
            column_list = self.get_loc_str("columnList")[RECOILLESS]
        else:
            raise ValueError

        self.tv["columns"] = column_list
        self.tv["show"] = "headings"

        self.tv.tag_configure(self.get_loc_str(POINT_PEAK_STAG), foreground="#2e8b57")
        self.tv.tag_configure(self.get_loc_str(POINT_PEAK_AVG), foreground="#2ca02c")
        self.tv.tag_configure(self.get_loc_str(POINT_PEAK_BREECH), foreground="orange")
        self.tv.tag_configure(self.get_loc_str(POINT_PEAK_SHOT), foreground="yellow green")
        self.tv.tag_configure(self.get_loc_str(POINT_BURNOUT), foreground="red")
        self.tv.tag_configure(self.get_loc_str(POINT_FRACTURE), foreground="brown")
        self.tv.tag_configure(self.get_loc_str(POINT_EXIT), foreground="steel blue")
        self.tv.tag_configure(self.get_loc_str(POINT_START), foreground="steel blue")
        self.tv.tag_configure(self.get_loc_str(COMPUTE), foreground="tan")

        self.tv.tag_configure("monospace", font=self.font)
        self.tv.tag_configure("error", font=self.font, foreground="dim gray")

        # we use a fixed width font so any char will do
        font_width, _ = self.font.measure("m"), self.font.metrics("linespace")

        win_width = self.tv.winfo_width()

        for i, column in enumerate(column_list):  # foreach column
            self.tv.heading(i, text=column, anchor="center")  # let the column heading = column name
            min_width = font_width * 14
            ini_width = max(win_width // len(self.tv["columns"]), min_width)
            self.tv.column(column, stretch=True, width=ini_width, minwidth=min_width, anchor="center")

        for i, row in enumerate(loc_table_data):
            self.tv.insert("", "end", str(i + 1), values=row, tags=(row[0].strip(), "monospace"))

    def export_table(self, gun_result: GunResult | RecoillessResult, acc_exp: int, normalized_filename: str):
        if gun_result is None:
            messagebox.showinfo(self.get_loc_str("excTitle"), self.get_loc_str("noDataMsg"))
            return
        file_name = filedialog.asksaveasfilename(
            title=self.get_loc_str("exportLabel"),
            filetypes=(("Comma Separated File", "*.csv"),),
            defaultextension=".csv",
            initialfile=normalized_filename,
        )

        if isinstance(gun_result, GunResult):
            gun_type = CONVENTIONAL
        elif isinstance(gun_result, RecoillessResult):
            gun_type = RECOILLESS
        else:
            raise ValueError

        column_list = self.get_loc_str("columnList")[gun_type]

        with open(file_name, "w", encoding="utf-8", newline="") as csvFile:
            # noinspection PyTypeChecker
            csv_writer = csv.writer(csvFile, delimiter=",", quoting=csv.QUOTE_MINIMAL)
            csv_writer.writerow(column_list)
            table_data = self.format_table(gun_result=gun_result, acc_exp=acc_exp)
            for line in table_data:
                csv_writer.writerow(line)

        messagebox.showinfo(self.get_loc_str("sucTitle"), self.get_loc_str("savedLocMsg") + " {:}".format(file_name))
