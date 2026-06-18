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

    os.system("pyinstaller " + name + ".spec" + " --noconfirm")


if __name__ == "__main__":
    # parser = argparse.ArgumentParser(description="script to generate PIBS executable")
    # parser.add_argument("--mult_file", action="store_true", default=False)
    #
    # args = parser.parse_args()

    generate_executables(False)
    generate_executables(True)
