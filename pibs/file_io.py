from __future__ import annotations

import json
from pathlib import Path
from tkinter import filedialog, messagebox
from typing import TYPE_CHECKING, Literal

from .misc import filenameize

if TYPE_CHECKING:
    from .ib_frame import InteriorBallisticsFrame


class FileIOManager:
    """Manages file I/O operations for the PIBS application."""

    def __init__(self, frame: InteriorBallisticsFrame):
        self.frame = frame

    def save(self):
        """Save gun design to JSON file."""
        frame = self.frame

        if not frame.has_data():
            messagebox.showinfo(frame.get_loc_str("excTitle"), frame.get_loc_str("noDataMsg"))
            return

        file_name = filedialog.asksaveasfilename(
            title=frame.get_loc_str("saveLabel"),
            filetypes=(("JSON file", "*.json"),),
            defaultextension=".json",
            initialfile=filenameize(frame.get_normalized_name()),
        )

        if not file_name:
            return

        data = frame.get_save_data()
        with open(file_name, "w", encoding="utf-8") as file:
            json.dump(data, file, indent="\t", ensure_ascii=False, sort_keys=True)

        messagebox.showinfo(frame.get_loc_str("sucTitle"), frame.get_loc_str("savedLocMsg") + f" {file_name}")

    def load_gun(self, initial_dir: str | None = None, file_path: str | None = None):
        """Load gun design from JSON file."""
        frame = self.frame

        file_name = file_path or filedialog.askopenfilename(
            title=frame.get_loc_str("loadLabel"),
            filetypes=(("JSON File", "*.json"),),
            defaultextension=".json",
            initialdir=initial_dir,
        )

        if not file_name:
            return

        frame.reset_entries()
        frame.set_name(Path(file_name).stem)

        with open(file_name, "r", encoding="utf-8") as file:
            data = json.load(file)

        frame.apply_loaded_data(data)
        frame.on_calculate()

    def load_propellant(self):
        """Load propellant from CSV file."""
        frame = self.frame

        file_name = filedialog.askopenfilename(
            title=frame.get_loc_str("loadLabel"), filetypes=(("Comma Separated Values File", "*.csv"),)
        )

        if file_name:
            frame.set_propellant_options(file_name)

    def export_table(self):
        """Export table data."""
        frame = self.frame
        frame.export_table_data()

    def export_graph(self, save: Literal["main", "aux", "geom", "guide"]):
        """Export plot as PNG."""
        frame = self.frame

        file_name = filedialog.asksaveasfilename(
            title=frame.get_loc_str("exportGraphLabel"),
            filetypes=(("Portable Network Graphics", "*.png"),),
            defaultextension=".png",
            initialfile=filenameize(f"{frame.get_normalized_name()}_{save}"),
        )

        if not file_name:
            return

        fig = frame.get_figure(save)
        if fig is None:
            return

        with frame.get_export_context():
            fig.savefig(file_name, transparent=True, dpi=300)

        messagebox.showinfo(frame.get_loc_str("sucTitle"), frame.get_loc_str("savedLocMsg") + f" {file_name}")
