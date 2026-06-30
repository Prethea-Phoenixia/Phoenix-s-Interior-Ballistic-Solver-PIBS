from tkinter import Label, StringVar, Toplevel


class ToolTip(object):
    """solution found in
    https://stackoverflow.com/questions/20399243/display-message-when-hovering-over-something-with-mouse-cursor-in-python
    """

    def __init__(self, widget, font):
        self.widget = widget
        self.tipwindow = None
        self.font = font
        self.text = None
        self.initialized: bool = False

    def showtip(self, text):
        """Display text in tooltip window"""
        if isinstance(text, StringVar):
            self.text = text.get()
        else:
            self.text = text

        if self.tipwindow or not self.text:
            return

        self.tipwindow = tw = Toplevel(self.widget)
        self.initialized = False
        tw.wm_overrideredirect(True)
        tw.attributes("-alpha", 0.0)
        t_font = self.font

        # we use a fixed width font so any char will do
        column_width = 60
        # apparently this doesn't work correctly with CJK fonts.....
        width, height = t_font.measure(" "), t_font.metrics("linespace")

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

        root_window = self.widget.winfo_toplevel()
        x, y = root_window.winfo_pointerxy()

        # ensure the tip shows up on the correct side:

        tw.wm_geometry("+%d+%d" % (x + 1, y + 1))
        tw.attributes("-alpha", 1.0)

        self.initialized = True

    def hidetip(self):
        tw = self.tipwindow
        if tw:
            if self.initialized:
                self.tipwindow = None
                tw.destroy()
            else:
                self.tipwindow.after(10, self.hidetip)


def create_tool_tip(widget, text, font):
    tool_tip = ToolTip(widget, font)

    def enter(_):
        tool_tip.showtip(text)

    def leave(_):
        tool_tip.hidetip()

    widget.bind("<Enter>", enter)
    widget.bind("<Leave>", leave)
