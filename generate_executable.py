import PyInstaller.__main__

if __name__ == "__main__":
    PyInstaller.__main__.run(
        [
            "run_pibs.py",
            "--clean",
            "--noconfirm",
            "--onefile",
            # "--onedir",
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
            "--splash",
            "splash.png",
        ],
    )


"""
# move this to spec file for production:
# a = Analysis...
# >here<
# pyz = PYZ...

key_words = ['api-ms-win', 'ucrtbase']
new_binaries = []
for item in a.binaries:
    name, _, _ = item
    to_include = True
    for key_word in key_words:
        if key_word in name:
            to_include = False
    if to_include:
        new_binaries.append(item)

a.binaries = new_binaries
"""
