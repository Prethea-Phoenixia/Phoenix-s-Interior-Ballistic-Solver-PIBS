import tkinter.font as tkFont
from tkinter import Label, StringVar, Toplevel


class ToolTip(object):
    """solution found in
    https://stackoverflow.com/questions/20399243/display-message
    -when-hovering-over-something-with-mouse-cursor-in-python
    """

    def __init__(self, widget, font):
        self.widget = widget
        self.tipwindow = None
        self.id = None
        self.x = self.y = 0
        self.font = font
        self.text = None

    def showtip(self, text):
        """Display text in tooltip window"""
        if isinstance(text, StringVar):
            self.text = text.get()
        else:
            self.text = text
        if self.tipwindow or not self.text:
            return

        self.tipwindow = tw = Toplevel(self.widget)
        tw.wm_overrideredirect(1)
        root = self.widget.winfo_toplevel()

        # t_Font = tkFont.Font(family="Sarasa Fixed SC", size=8)

        t_font = self.font
        # we use a fixed width font so any char will do
        column_width = 60
        # apparently this doesn't work correctly with CJK fonts.....
        width, height = t_font.measure(" "), t_font.metrics("linespace")

        x, y, _, _ = self.widget.bbox("insert")
        # rx, ry, crx, cry = root.bbox()
        # bounding box coordinate is in relation to origin of widget/window

        if x + self.widget.winfo_rootx() > root.winfo_rootx() + 0.5 * root.winfo_width():
            x = x + self.widget.winfo_rootx() - width * column_width
            y = y + self.widget.winfo_rooty()
        else:
            x = x + self.widget.winfo_rootx() + self.widget.winfo_width()
            y = y + self.widget.winfo_rooty()

        # wrappedText = textwrap.fill(self.text, width=columnWidth - 10)

        label = Label(
            tw,
            text=self.text,
            justify="left",
            background="#ffffe0",
            wraplength=(column_width - 2) * width,  # 2 to account for the padding.
            relief="solid",
            borderwidth=0,
            font=t_font,
        )
        label.pack(ipadx=width, ipady=height * 0.25, anchor="nw", fill="both")

        # tw.update_idletasks()
        wheight = tw.winfo_height()

        margin = (
            y + wheight - root.winfo_rooty() - root.winfo_height()
        )  # ensure that tooltip does not overrun the main window

        if margin > 0:
            y -= margin

        tw.wm_geometry("+%d+%d" % (x, y))
        """
        tw.bind(
            "<Configure>",
            lambda event: label.configure(
                wraplength=label.winfo_width() * 0.75
            ),
        )
        """

    def hidetip(self):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.after(100, lambda: tw.destroy())


def create_tool_tip(widget, text, font):
    tool_tip = ToolTip(widget, font)

    def enter(_):
        tool_tip.showtip(text)

    def leave(_):
        tool_tip.hidetip()

    widget.bind("<Enter>", enter)
    widget.bind("<Leave>", leave)
