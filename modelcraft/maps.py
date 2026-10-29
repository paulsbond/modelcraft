import gemmi


def read_map(path: str) -> gemmi.Ccp4Map:
    ccp4 = gemmi.read_ccp4_map(path, setup=True)
    if any(abs(ccp4.header_float(w)) > 0 for w in (50, 51, 52)):
        print(f"Non-zero origin in {path} is not supported and will be ignored.")
    for word_number in (50, 51, 52):
        ccp4.set_header_float(word_number, 0.0)
    return ccp4
