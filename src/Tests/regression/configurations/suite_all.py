"""suite_all regression configurations (aggregate of all group modules; new group).

Sibling-style imports: the loader (``suite.config_loader.load_configurations``)
puts this module's own directory on ``sys.path`` before executing it, so the
group modules import as plain top-level modules (loader contract; the package
``__init__`` is deliberately import-free).
"""

from wcsph_dbc import configurations as _dbc
from wcsph_bi import configurations as _bi
from rsph import configurations as _rsph
from defaults import configurations as _defaults
from multiresolution import configurations as _mr

configurations = _dbc + _bi + _rsph + _defaults + _mr


def groups() -> dict:
    """Return the aggregate's source lists keyed by group name."""
    return {
        "wcsph_dbc": _dbc,
        "wcsph_bi": _bi,
        "rsph": _rsph,
        "defaults": _defaults,
        "multiresolution": _mr,
    }
