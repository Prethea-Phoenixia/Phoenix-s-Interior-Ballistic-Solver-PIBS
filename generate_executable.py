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

    subprocess.run(["pyinstaller", f"{name}.spec", "--noconfirm", "--clean"], check=True)



if __name__ == "__main__":
    shutil.rmtree(os.path.join(str(os.path.dirname(__file__)), "dist"), ignore_errors=True)

    generate_executables(mult_file=True)
    generate_executables(mult_file=False)
