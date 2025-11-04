import os
import PyInstaller.__main__
from pibs import __version__

if __name__ == "__main__":
    name = "PIBSv" + __version__

    sep = os.pathsep  # ':' on POSIX, ';' on Windows

    options = [
        "run_pibs.py",
        "--clean",
        "--noconfirm",
        "--windowed",
        "--icon=pibs/ui/logo.ico",
        "--name=" + name,
        "--debug=all",
        f"--add-data=pibs/ballistics/resource{sep}ballistics/resource/",
        f"--add-data=pibs/ui{sep}ui/",
        f"--add-data=pibs/examples{sep}examples/",
        "--exclude-module=pibs.ballistics.therm",
        "--onefile",
    ]

    PyInstaller.__main__.run(options)

    with open(name + ".spec", "r") as f:
        content = f.readlines()

    i = content.index("pyz = PYZ(a.pure)\n")
    """
    The following code injects code defined in `dll_exclusion` to the specs file, used by pyinstaller. 
        a = Analysis...
        >here<
        pyz = PYZ...
    excluded keywords corresponds to:
    - api-ms-win, ucrtbase: MikTex (LaTeX) compiler on Windows machines
    These can be safely excluded since the functionality of PIBS does not depend on them.
    
    """
    dll_exclusion = """
# exclude excessive DLL collected by pyinstaller
key_words = ['api-ms-win', 'ucrtbase']
new_binaries = []
excluded = []
for item in a.binaries:
    name, _, _ = item
    to_include = True
    for key_word in key_words:
        if key_word in name:
            to_include = False
    if to_include:
        new_binaries.append(item)
    else:
        excluded.append(item)

a.binaries = new_binaries

# Filter out all entries in matplotlib/mpl-data/fonts
mpl_font_path = os.path.join('matplotlib', 'mpl-data', 'fonts')
a.datas = [entry for entry in a.datas if not entry[0].startswith(mpl_font_path)]

mpl_data_path = os.path.join('matplotlib', 'mpl-data', 'sample_data')
a.datas = [entry for entry in a.datas if not entry[0].startswith(mpl_data_path)]
"""

    for line in dll_exclusion.split("\n"):
        content.insert(i, line + "\n")
        i += 1

    with open(name + ".spec", "w") as f:
        f.writelines(content)

    os.system("pyinstaller " + name + ".spec" + " --noconfirm")
