"""Import command-line helper modules without installing a site package."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts"))
