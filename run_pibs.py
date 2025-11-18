import sys

from pibs.interior_ballistics import main

if __name__ == "__main__":
    if getattr(sys, "frozen", False):  # if PIBS is run as frozen installation, start without debug flag.
        main(debug=False)
    else:  # if PIBS is run as a script, start with debug flag.
        main(debug=True)
