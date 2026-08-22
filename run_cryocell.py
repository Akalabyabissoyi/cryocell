#!/usr/bin/env python3
"""Launch CryoCell.   python3 run_cryocell.py"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cryocell.app import main
if __name__ == "__main__":
    main()
