import gemmi


def write_map(path: str, grid: gemmi.FloatGrid):
    ccp4 = gemmi.Ccp4Map()
    ccp4.grid = grid
    ccp4.update_ccp4_header()
    ccp4.write_ccp4_map(path)
