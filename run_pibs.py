from pibs.interior_ballistics import main
import sys

if __name__ == "__main__":
    # if PIBS is run as a script, start with debug flag.
    if getattr(sys, "frozen", False):
        main(debug=False)
    else:
        main(debug=True)
