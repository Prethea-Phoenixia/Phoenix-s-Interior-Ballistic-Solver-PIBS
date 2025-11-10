
## Phoenix's Interior Ballistics Solver PIBS
Application for solving the interior ballistics problems of closed-breech and recoilless weapon, with provisions for constrained design and certain optimization calculations.    

## Features
### Modelling
The interior ballistics problem is formulated after the system named after M.E. Serebryakov, widely used in the Soviet Union and People's Republic of China. The calculation is done in the reduced form for conventional and recoilless guns.

Supports shaped propellant described by cubic form function, including sphere, strip, cylinder, and multi-perforated cylinder/prism/rosette grains.

### Operation
Supported operations includes calculation from parameters (forward calculation), and solving parameters required to hit performance targets. 

Additionally, solves the minimum bore volume problem for a given performance criteria.

### Data
A preset of propellant is supplied, with definition for many propellants, sourced from public literature. New propellant definition may be hot-loaded.

To facilitate interoperability with design processes, data export and design save/load is provided. Interior ballistic configuration for known guns are supplied as presets for reference.

### Technical
Compromising between ease of development, distribution, and runtime speed, the following algorithms were implemented in Python for this application:

* Numerical integration up to user specified precision using high-order, adaptive [Runge Kutta Fehlberg](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method) method.
* Root-finding is implemented using the Dekker method, a variant of Brent's method that opportunistically employ polynomial and bisection to speedup root finding.
* Maximal values are found using the Gold Section Search (GSS) method

## Instructions
### For Use as Application:
For Windows, a single-file executable is provided. This package the Python runtime, dependent libraries and data, and the code in a single self-extracting archive that is extracted to the local temporary folder upon use.
- download the latest executable at [latest release](https://github.com/octo-org/octo-repo/releases/latest)

### For Development:

- Install Python (>=3.9)
- Setup the virtual environment:`python -m venv .venv`
- Activate the virtual environment venv 
  - on Windows: `.venv\Scripts\activate.bat`
  - on Linux: `source .venv/bin/activate`
- Install development dependencies:`pip install .[dev]`
- Run the main script: `python run_pibs.py`
- Generate executables: `python generate_executable.py`
      

## Contribution
Your contribution is welcomed! Please feel free to raise issues or propose pull requests regarding features or problems encountered. To get started, see the [community documentation](CONTRIBUTING.md)


## Resources Used
* tcl/tk themes used include "awdark" & "awlight" from awthemes: https://wiki.tcl-lang.org/page/awthemes
* monospaced, CJK compatible font: Sarasa-Gothic, https://github.com/be5invis/Sarasa-Gothic

