import os
import PyInstaller.__main__

from pibs import __version__


def generate_executables(mult_file: bool = False):
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
        "--noupx",
        f"--add-data=pibs/ballistics/resource{sep}ballistics/resource/",
        f"--add-data=pibs/ui{sep}ui/",
        f"--add-data=pibs/examples{sep}examples/",
    ] + (["--onefile"] if not mult_file else [])

    PyInstaller.__main__.run(options)

    with open(name + ".spec", "r") as f:
        content = f.readlines()

    with open(name + ".spec", "w") as f:
        f.writelines(content)

    i = content.index("pyz = PYZ(a.pure)\n")

    dll_exclusion = """# exclude excessive DLL collected by pyinstaller
key_words = ['api-ms-win']
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
"""

    for line in dll_exclusion.split("\n"):
        content.insert(i, line + "\n")
        i += 1

    with open(name + ".spec", "w") as f:
        f.writelines(content)

    import sys
    os.system(f'"{sys.executable}" -m PyInstaller ' + name + ".spec" + " --noconfirm")


if __name__ == "__main__":
    generate_executables(False)
    generate_executables(True)
