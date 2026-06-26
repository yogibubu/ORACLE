from pathlib import Path

from gui.xyz_reader import read_xyz
from .elements import atomic_number
from .pipeline import build_topology_objects
from .topology_writer import write_topology_to_xyzin
from .topology_reporting import print_topology_report


def build_topology(xyzin_path: Path, *, force_aromatic=False):
    """
    Build molecular topology from the current xyzin file.

    This function:
    - reads geometry from xyzin
    - constructs continuous and discrete topology
    - updates xyzin via topology_writer
    - writes topology.report via topology_reporting

    Parameters
    ----------
    xyzin_path : Path
        Path to the xyzin file in working directory.
    force_aromatic : bool, optional
        Force aromaticity perception.
    """

    xyzin_path = Path(xyzin_path)

    # --------------------------------------------------
    # Read geometry
    # --------------------------------------------------
    atoms, coords = read_xyz(xyzin_path)

    Z = [atomic_number(sym) for sym in atoms]
    if any(z is None for z in Z):
        raise ValueError("Unknown atomic symbol in XYZ input")

    # --------------------------------------------------
    # Build topology objects
    # --------------------------------------------------
    cg, dg, ringset, synthons, aromaticity = build_topology_objects(
        coords,
        Z,
        force_aromatic=force_aromatic,
    )

    # --------------------------------------------------
    # Update xyzin (state file)
    # --------------------------------------------------
    write_topology_to_xyzin(
        xyzin_path.name,
        cg=cg,
        dg=dg,
        ringset=ringset,
        aromaticity=aromaticity,
        synthons=synthons,
    )

    # --------------------------------------------------
    # Write topology report
    # --------------------------------------------------
    print_topology_report(
        cg,
        dg,
        synthons,
        arom=aromaticity,
        ringset=ringset,
    )

    return cg, dg, ringset, synthons, aromaticity
