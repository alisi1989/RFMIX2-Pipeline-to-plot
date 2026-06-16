#!/usr/bin/env python3
"""Compatibility wrapper for the GAP v2 plotter.

The historical AncestryGrapher command used `Scripts/GAP.py` to draw the
global ancestry bar plot. The current plotting implementation lives in
`GAP_Plot.py`; this wrapper preserves the old entry point.
"""

from GAP_Plot import main


if __name__ == "__main__":
    raise SystemExit(main())
