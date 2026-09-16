import os
import shutil
import subprocess

from pibs import __version__


def generate_executables(mult_file: bool = False):
    """Build executable using PyInstaller locally."""

    name = "PIBSv" + __version__

    sep = os.pathsep  # ':' on POSIX, ';' on Windows

    options = [
        "run_pibs.py",
        "--windowed",
        "--icon=pibs/ui/logo.ico",
        "--name=" + name,
        "--optimize=2",
        "--noupx",
        f"--add-data=pibs/ballistics/resource{sep}ballistics/resource/",
        f"--add-data=pibs/ui{sep}ui/",
        f"--add-data=pibs/examples{sep}examples/",
    ] + (["--onefile"] if not mult_file else [])

    subprocess.run(["pyi-makespec"] + options, check=True)

    new_hook = """{
        "matplotlib": {
            "backends": "TkAgg",
        },
    }"""

    with open(name + ".spec", "r") as f:
        stuff = f.read()
        stuff = stuff.replace("hooksconfig={}", f"hooksconfig={new_hook}")

    with open(name + ".spec", "w") as f:
        f.write(stuff)

#     with open(name + ".spec", "r") as f:
#         content = f.readlines()
#
#     i = content.index("pyz = PYZ(a.pure)\n")
#
#     dll_exclusion = """# exclude excessive DLL collected by pyinstaller
# key_words = ['api-ms-win', 'ucrtbase']
# new_binaries = []
# excluded = []
# for item in a.binaries:
#     name, _, _ = item
#     to_include = True
#     for key_word in key_words:
#         if key_word in name:
#             to_include = False
#     if to_include:
#         new_binaries.append(item)
#     else:
#         excluded.append(item)
#
# a.binaries = new_binaries"""
#
#     for line in dll_exclusion.split("\n"):
#         content.insert(i, line + "\n")
#         i += 1
#
#     with open(name + ".spec", "w") as f:
#         f.writelines(content)

    subprocess.run(["pyinstaller", f"{name}.spec", "--noconfirm", "--clean"], check=True)


if __name__ == "__main__":
    shutil.rmtree(os.path.join(str(os.path.dirname(__file__)), "dist"), ignore_errors=True)

    generate_executables(mult_file=True)
    generate_executables(mult_file=False)
