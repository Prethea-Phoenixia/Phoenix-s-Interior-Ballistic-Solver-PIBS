# Ideas

These are some outstanding aspect that I would much appreciate help in.

## Linux Compatibility
It *should be* relatively straightforward to run PIBS on various Linux distribution, although the maintainer does not possess a Linux device, functionality cannot be tested nor guaranteed. Specifically, the following aspects would need work:

- Factoring out Windows specific code in graphic interface.
- Ensuring file menu operation under Mac OSX due to limitation of tk/tkinter.
- (Potentially) refactor multiprocessing due to differences in how Python implements multiprocessing across Windows and Linux.

## Designs and Propellants
Curated example designs are stored in the [examples folder](/examples). Custom propellant presets are defined in the [propellants.csv](/pibs/ballistics/resource/propellants.csv) file. If possible, include citation or sources in the description for future reference.


