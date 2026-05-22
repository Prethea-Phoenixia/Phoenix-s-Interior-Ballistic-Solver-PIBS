# -*- mode: python ; coding: utf-8 -*-


a = Analysis(
    ['run_pibs.py'],
    pathex=[],
    binaries=[],
    datas=[('pibs/ballistics/resource', 'ballistics/resource/'), ('pibs/ui', 'ui/'), ('pibs/examples', 'examples/')],
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=True,
    optimize=0,
)
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

# Filter out all 
mpl_data_path = os.path.join('matplotlib', 'mpl-data', 'sample_data')
a.datas = [entry for entry in a.datas if not entry[0].startswith(mpl_data_path)]

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [('v', None, 'OPTION')],
    exclude_binaries=True,
    name='PIBSv0.5.5',
    debug=True,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=['pibs\\ui\\logo.ico'],
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='PIBSv0.5.5',
)
