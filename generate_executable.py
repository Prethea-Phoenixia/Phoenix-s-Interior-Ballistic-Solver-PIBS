import os

import PyInstaller.__main__

if __name__ == "__main__":

    options = [
        "run_pibs.py",
        "--clean",
        "--noconfirm",
        "--windowed",
        "--icon",
        "pibs/ui/logo.ico",
        "--name",
        "pibs",
        "--debug",
        "all",
        "--add-data",
        "pibs/ballistics/resource;ballistics/resource/",
        "--add-data",
        "pibs/ui;ui/",
        "--exclude-module",
        "pyinstaller",
    ]
    # use this for dist
    PyInstaller.__main__.run(options + ["--onefile"])

    with open("pibs.spec", "r") as f:
        content = f.readlines()

    i = content.index("pyz = PYZ(a.pure)\n")
    """
    move this to spec file for production:
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

    with open("pibs.spec", "w") as f:
        f.writelines(content)

    os.system("pyinstaller pibs.spec --noconfirm")
