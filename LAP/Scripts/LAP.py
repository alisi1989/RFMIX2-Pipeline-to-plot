#!/usr/bin/env python3
"""Compatibility wrapper for the LAP v4 plotter.

The historical AncestryGrapher command used `Scripts/LAP.py` for plotting.
The current implementation lives in `Plot_LAP_v4.py`; this wrapper preserves
the old entry point while exposing all current plotting options.
"""

from Plot_LAP_v4 import main


if __name__ == "__main__":
    raise SystemExit(main())
