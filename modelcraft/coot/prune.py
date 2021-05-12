import gtk
import math

# CONSTS

atomic_numbers = {
    "H": 1,
    "HE": 2,
    "LI": 3,
    "BE": 4,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "NE": 10,
    "NA": 11,
    "MG": 12,
    "AL": 13,
    "SI": 14,
    "P": 15,
    "S": 16,
    "CL": 17,
    "AR": 18,
    "K": 19,
    "CA": 20,
    "SC": 21,
    "TI": 22,
    "V": 23,
    "CR": 24,
    "MN": 25,
    "FE": 26,
    "CO": 27,
    "NI": 28,
    "CU": 29,
    "ZN": 30,
    "GA": 31,
    "GE": 32,
    "AS": 33,
    "SE": 34,
    "BR": 35,
    "KR": 36,
    "RB": 37,
    "SR": 38,
    "Y": 39,
    "ZR": 40,
    "NB": 41,
    "MO": 42,
    "TC": 43,
    "RU": 44,
    "RH": 45,
    "PD": 46,
    "AG": 47,
    "CD": 48,
    "IN": 49,
    "SN": 50,
    "SB": 51,
    "TE": 52,
    "I": 53,
    "XE": 54,
    "CS": 55,
    "BA": 56,
    "LA": 57,
    "CE": 58,
    "PR": 59,
    "ND": 60,
    "PM": 61,
    "SM": 62,
    "EU": 63,
    "GD": 64,
    "TB": 65,
    "DY": 66,
    "HO": 67,
    "ER": 68,
    "TM": 69,
    "YB": 70,
    "LU": 71,
    "HF": 72,
    "TA": 73,
    "W": 74,
    "RE": 75,
    "OS": 76,
    "IR": 77,
    "PT": 78,
    "AU": 79,
    "HG": 80,
    "TL": 81,
    "PB": 82,
    "BI": 83,
    "PO": 84,
    "AT": 85,
    "RN": 86,
    "FR": 87,
    "RA": 88,
    "AC": 89,
    "TH": 90,
    "PA": 91,
    "U": 92,
    "NP": 93,
    "PU": 94,
    "AM": 95,
    "CM": 96,
    "BK": 97,
    "CF": 98,
    "ES": 99,
    "FM": 100,
    "MD": 101,
    "NO": 102,
    "LR": 103,
    "RF": 104,
    "DB": 105,
    "SG": 106,
    "BH": 107,
    "HS": 108,
    "MT": 109,
    "DS": 110,
    "RG": 111,
    "CN": 112,
    "NH": 113,
    "FL": 114,
    "MC": 115,
    "LV": 116,
    "TS": 117,
    "OG": 118,
}

protein_residues = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "MSE",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "UNK",
    "VAL",
}

main_chain_atoms = {" N  ", " CA ", " C  ", " O  ", " CB "}

atoms = {
    "ALA": {" N  ", " CA ", " C  ", " O  ", " CB "},
    "ARG": {" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " CD ", " NE ", " CZ ", " NH1", " NH2",},
    "ASN": {" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " OD1", " ND2"},
    "ASP": {" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " OD1", " OD2"},
    "CYS": {" N  ", " CA ", " C  ", " O  ", " CB ", " SG "},
    "GLN": {" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " CD ", " OE1", " NE2"},
    "GLU": {" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " CD ", " OE1", " OE2"},
    "GLY": {" N  ", " CA ", " C  ", " O  "},
    "HIS": {" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " ND1", " CD2", " CE1", " NE2",},
    "ILE": {" N  ", " CA ", " C  ", " O  ", " CB ", " CG1", " CG2", " CD1"},
    "LEU": {" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " CD1", " CD2"},
    "LYS": {" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " CD ", " CE ", " NZ "},
    "MET": {" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " SD ", " CE "},
    "MSE": {" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", "SE  ", " CE "},
    "PHE": {" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " CD1", " CD2", " CE1", " CE2", " CZ ",},
    "PRO": {" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " CD "},
    "SER": {" N  ", " CA ", " C  ", " O  ", " CB ", " OG "},
    "THR": {" N  ", " CA ", " C  ", " O  ", " CB ", " OG1", " CG2"},
    "TRP": {
        " N  ",
        " CA ",
        " C  ",
        " O  ",
        " CB ",
        " CG ",
        " CD1",
        " CD2",
        " NE1",
        " CE2",
        " CE3",
        " CZ2",
        " CZ3",
        " CH2",
    },
    "TYR": {" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " CD1", " CD2", " CE1", " CE2", " CZ ", " OH ",},
    "VAL": {" N  ", " CA ", " C  ", " O  ", " CB ", " CG1", " CG2"},
}

bonded = {
    " C  ": {" CA ", " O  ", " OXT"},
    " CA ": {" C  ", " CB ", " N  "},
    " CB ": {" CA ", " CG ", " CG1", " CG2", " OG ", " OG1", " SG "},
    " CD ": {" CE ", " CG ", " NE ", " NE2", " OE1", " OE2"},
    " CD1": {" CE1", " CG ", " CG1", " NE1"},
    " CD2": {" CE2", " CE3", " CG ", " NE2"},
    " CE ": {" CD ", " NZ ", " SD ", "SE  "},
    " CE1": {" CD1", " CZ ", " ND1", " NE2"},
    " CE2": {" CD2", " CZ ", " CZ2", " NE1"},
    " CE3": {" CD2", " CZ3"},
    " CG ": {" CB ", " CD ", " CD1", " CD2", " ND1", " ND2", " OD1", " OD2", " SD ", "SE  ",},
    " CG1": {" CB ", " CD1"},
    " CG2": {" CB "},
    " CH2": {" CZ2", " CZ3"},
    " CZ ": {" CE1", " CE2", " NE ", " NH1", " NH2", " OH "},
    " CZ2": {" CE2", " CH2"},
    " CZ3": {" CE3", " CH2"},
    " N  ": {" CA "},
    " ND1": {" CE1", " CG "},
    " ND2": {" CG "},
    " NE ": {" CD ", " CZ "},
    " NE1": {" CD1", " CE2"},
    " NE2": {" CD ", " CD2", " CE1"},
    " NH1": {" CZ "},
    " NH2": {" CZ "},
    " NZ ": {" CE "},
    " O  ": {" C  "},
    " OD1": {" CG "},
    " OD2": {" CG "},
    " OE1": {" CD "},
    " OE2": {" CD "},
    " OG ": {" CB "},
    " OG1": {" CB "},
    " OH ": {" CZ "},
    " OXT": {" C  "},
    " SD ": {" CE ", " CG "},
    " SG ": {" CB "},
    "SE  ": {" CE ", " CG "},
}

# UTILS


def mean(values):
    return float(sum(values)) / len(values)


def median(values):
    n = len(values)
    if n < 1:
        return None
    i = n // 2
    if n % 2 == 1:
        return sorted(values)[i]
    return sum(sorted(values)[i - 1 : i + 1]) / 2.0


def median_absolute_deviation(values):
    median_value = median(values)
    abs_deviations = [abs(value - median_value) for value in values]
    return median(abs_deviations)


def dot(xyz1, xyz2):
    return xyz1[0] * xyz2[0] + xyz1[1] * xyz2[1] + xyz1[2] * xyz2[2]


def cross(xyz1, xyz2):
    return [
        xyz1[1] * xyz2[2] - xyz1[2] * xyz2[1],
        xyz1[2] * xyz2[0] - xyz1[0] * xyz2[2],
        xyz1[0] * xyz2[1] - xyz1[1] * xyz2[0],
    ]


def magnitude(xyz):
    return (xyz[0] ** 2 + xyz[1] ** 2 + xyz[2] ** 2) ** 0.5


def unit(xyz):
    length = magnitude(xyz)
    return [xyz[0] / length, xyz[1] / length, xyz[2] / length]


def subtract(xyz1, xyz2):
    return [xyz1[0] - xyz2[0], xyz1[1] - xyz2[1], xyz1[2] - xyz2[2]]


def distance(xyz1, xyz2):
    v = subtract(xyz1, xyz2)
    return magnitude(v)


def angle(xyz1, xyz2, xyz3):
    v1 = subtract(xyz2, xyz1)
    v2 = subtract(xyz2, xyz3)
    angle = math.acos(dot(v1, v2) / (magnitude(v1) * magnitude(v2)))
    return math.degrees(angle)


def torsion(xyz1, xyz2, xyz3, xyz4):
    b1 = subtract(xyz2, xyz1)
    b2 = subtract(xyz3, xyz2)
    b3 = subtract(xyz4, xyz3)
    n1 = cross(b1, b2)
    n2 = cross(b2, b3)
    m1 = cross(n1, n2)
    y = dot(m1, unit(b2))
    x = dot(n1, n2)
    angle = math.degrees(math.atan2(y, x))
    if angle > 180:
        angle -= 360
    return angle


def halfway(xyz1, xyz2):
    x = (xyz1[0] + xyz2[0]) / 2
    y = (xyz1[1] + xyz2[1]) / 2
    z = (xyz1[2] + xyz2[2]) / 2
    return [x, y, z]


def is_protein(imol, name, res_spec):
    if name in protein_residues:
        return True
    atom_names = []
    atom_points = []
    for info in residue_info_py(imol, *res_spec):
        atom_names.append(info[0][0].strip())
        atom_points.append(info[2])
    if any(x not in atom_names for x in ("N", "CA", "C")):
        return False
    natoms = len(atom_names)
    n = [atom_points[i] for i in range(natoms) if atom_names[i] == "N"][0]
    ca = [atom_points[i] for i in range(natoms) if atom_names[i] == "CA"][0]
    c = [atom_points[i] for i in range(natoms) if atom_names[i] == "C"][0]
    return distance(n, ca) < 1.8 and distance(ca, c) < 1.8


def attached_atoms(atom, residue):
    if atom.name == " N  " and residue.prev is not None:
        yield residue.prev.atoms[" C  "]
    if atom.name == " C  " and residue.next is not None:
        yield residue.next.atoms[" N  "]
    for other in bonded[atom.name]:
        if other in residue.atoms:
            yield residue.atoms[other]


# MACHINE LEARNING

training_data = {
    "main": {
        "medians": [
            0.8093730211257935,
            0.0,
            0.3797053494361412,
            0.2948472222222224,
            1.1603445955372607,
            0.07967734307025365,
            -1.04405097291337,
            -0.7666229846405788,
            0.09301174988214152,
            0.2074185971626753,
            1.9998366832733154,
            7.473087887125843,
        ],
        "scaler.mean_": [
            0.7635879979927991,
            0.0119850162929945,
            0.47070085149410484,
            0.8779728019083975,
            1.4010333429538593,
            -0.07809440994354622,
            -1.1634877907354524,
            -0.9231595065864178,
            0.06630610233126723,
            0.7376180259790496,
            1.9333955416970983,
            9.00267277208097,
        ],
        "scaler.scale_": [
            0.15510631635324296,
            0.10881808524983969,
            0.49424003745942746,
            2.1942551160797694,
            1.298509344613297,
            0.8836945495954351,
            0.9345174344061477,
            0.9874691652931548,
            1.2037765763286472,
            1.1729051296356117,
            0.42936735360011963,
            6.84475636137809,
        ],
        "regressor.coefs_": [
            [
                [
                    -0.13040818548223823,
                    -0.22295218610758016,
                    -0.16258080220962234,
                    0.5988690897292985,
                    0.29886699802544087,
                    -0.2138536188789283,
                    0.17640972642170277,
                    0.8323631686040229,
                    0.8080535193667213,
                    -0.64199269078319,
                ],
                [
                    0.3371557445052779,
                    -0.021990355769565446,
                    0.3393433644319566,
                    0.7031031102737114,
                    -1.5952677749631456,
                    -0.6291226390140611,
                    -0.16470505242375907,
                    0.8959152619464158,
                    0.8165087694148113,
                    1.3192349998470725,
                ],
                [
                    -0.0019874340519655917,
                    0.04016420712930254,
                    -0.0016519722530829112,
                    -0.05753723016959708,
                    -1.0970634376604886,
                    0.13781915021293673,
                    -0.043726103959200804,
                    0.010119141446418834,
                    -0.09415567058366356,
                    0.10675615201134665,
                ],
                [
                    -0.15553908141483186,
                    0.7815309718136879,
                    -0.58301702342861,
                    -0.3074194496122189,
                    -0.708058122322275,
                    0.01629181598284882,
                    -0.17713645375637574,
                    -0.5780953291256867,
                    -0.08330249685202995,
                    -0.046797521299526336,
                ],
                [
                    -0.16382881259387705,
                    0.1590019158219002,
                    -0.15271138566338957,
                    0.004631275745636504,
                    0.14085902736133887,
                    0.044487235857714004,
                    -0.22845244910340506,
                    0.03361790942261084,
                    -0.09431223857684323,
                    0.1894309630932748,
                ],
                [
                    0.4600605680220205,
                    1.0475854981548196,
                    0.28528470975660497,
                    -0.07122257966372195,
                    -0.17193073732931968,
                    -0.563441753208365,
                    -0.6552700463525175,
                    -0.08483854601112313,
                    -0.17645291495517423,
                    -0.44659761042818563,
                ],
                [
                    -0.32437585347959863,
                    -0.40621672156363603,
                    -0.22786021141673057,
                    -0.006696637890890252,
                    0.5191491920623607,
                    0.16595743564965426,
                    0.27911102055393494,
                    0.2640726988201516,
                    -0.0016347547075830314,
                    0.24770754766817038,
                ],
                [
                    0.3728665898401538,
                    0.07700714533797552,
                    0.2311708633268853,
                    -0.1401187037969537,
                    0.30234045086128297,
                    -0.24458025060957192,
                    0.4302767405048504,
                    -0.3703955591237981,
                    0.3484976393563264,
                    -0.5594717372377868,
                ],
                [
                    0.00937493962029786,
                    -0.04694394839682417,
                    -0.028893614058871878,
                    0.27141499045973605,
                    -0.07941767649526653,
                    -0.03911716545463574,
                    -0.12103549705527984,
                    0.2148437699515711,
                    0.08720611510332744,
                    -0.13302591593197644,
                ],
                [
                    0.16013205898626096,
                    0.35382721838123166,
                    0.05844089574619886,
                    -0.2090590614354888,
                    0.06453171825149609,
                    -1.172220172769279,
                    -0.14500654761975515,
                    -0.35695026385344086,
                    0.7073129963007073,
                    -1.1972037303520695,
                ],
                [
                    0.12999852919812566,
                    -0.37003217528080845,
                    0.3596837609653873,
                    0.30447584223950475,
                    -0.5844046567044845,
                    1.099321801469266,
                    0.28973479295322985,
                    0.6123382299120337,
                    -0.8254996440849733,
                    0.3860328709153097,
                ],
                [
                    -0.013864646575082658,
                    0.33883378654775376,
                    -0.0432656081821551,
                    -0.03831855907564628,
                    -0.17009170626798892,
                    0.09755681267845129,
                    -0.14149706478955423,
                    0.039883119397616046,
                    -0.11814265819588236,
                    0.13867097327343803,
                ],
            ],
            [
                [0.4113990249648804],
                [-0.1635521743014414],
                [-0.3233742954955824],
                [-0.2832881263154175],
                [0.10184600411275878],
                [0.26297418395412697],
                [-0.38667783069289213],
                [0.18806411489930078],
                [0.2650096089001783],
                [-0.3432819700082653],
            ],
        ],
        "regressor.intercepts_": [
            [
                0.05426873271484922,
                -0.17185888723743994,
                0.49670571655441914,
                -0.26895026535845296,
                1.2489563308304976,
                -0.21430979896458663,
                -1.0229946977640407,
                0.43986781805151526,
                0.32517336900734384,
                -1.1981611190470425,
            ],
            [0.2009473880274795],
        ],
    },
    "side": {
        "medians": [
            0.7128432393074036,
            0.3930308684784101,
            0.12138492565055738,
            0.7697625388342151,
            -0.01155913420195338,
            -0.7577415671212863,
            -0.5053093547058974,
            1.9498454332351685,
            36.654449462890625,
        ],
        "scaler.mean_": [
            0.6699939922392967,
            0.4861278093279958,
            0.5598238115764408,
            1.0098633759640239,
            -0.09092168866565035,
            -0.8734938716432699,
            -0.5939212054114583,
            1.903840522636135,
            41.034051026184834,
        ],
        "scaler.scale_": [
            0.1922931323119439,
            0.5148420550064197,
            1.547688727501242,
            1.6216833496523257,
            0.8535355005798227,
            1.063733997556021,
            0.946219729526465,
            0.42543155678705524,
            29.99679803225501,
        ],
        "regressor.coefs_": [
            [
                [
                    0.056302872646141214,
                    -0.42725673599695374,
                    0.23328615715807235,
                    -0.25692317829492684,
                    0.116828236758043,
                    1.2384242768052554,
                    0.6748323891517939,
                    0.14291234817220463,
                    0.21380028409294824,
                    0.07679070847654368,
                ],
                [
                    -0.46153688659142184,
                    -0.004544981320490466,
                    -0.2964031182985578,
                    0.09456605767558321,
                    -0.06435341411291445,
                    0.12673028609583875,
                    -0.2719367005541086,
                    0.29679974643773216,
                    0.11486700383265894,
                    0.23185475607178951,
                ],
                [
                    0.33853298205315957,
                    0.9079813683093967,
                    0.09346553254580311,
                    -0.025784490716986762,
                    -0.3475803251546684,
                    -0.04639625178653906,
                    -0.500881871296081,
                    0.7809152729815941,
                    0.406203973798232,
                    -0.3284891104603179,
                ],
                [
                    -0.2088543439928225,
                    0.47556925150407825,
                    -0.19039046477157073,
                    0.20421072806892146,
                    -0.23471331429141476,
                    -0.08725373460798075,
                    -0.0811783537165059,
                    0.4593053196175982,
                    0.4351302763377739,
                    0.11837686827471676,
                ],
                [
                    -0.32138852798617806,
                    0.11322678209797299,
                    -0.07514559826050492,
                    0.13497722508180546,
                    0.2308583983764714,
                    0.11180021860598335,
                    -0.6957881978609368,
                    0.0005190646203973872,
                    -0.2880602141275108,
                    -0.13675732785562053,
                ],
                [
                    0.06423107339757392,
                    0.08120627709915425,
                    0.29375027114694363,
                    -0.15903515362475384,
                    -0.038018112519918414,
                    0.2893167384125542,
                    0.14378196316887357,
                    0.04488068580666965,
                    0.1342358964548304,
                    0.00855401074520736,
                ],
                [
                    -0.06837070927792858,
                    -0.45192674269454314,
                    0.059932235100582196,
                    -0.15614805983481528,
                    -0.3859760254395101,
                    -0.3429573286705649,
                    0.6838364312078055,
                    -0.347976240289938,
                    0.049432609557722565,
                    -0.28055309729723454,
                ],
                [
                    0.6431046107009092,
                    0.49740363903294654,
                    0.2748436188374371,
                    -0.18326862373951083,
                    -0.2069737876487165,
                    -0.2999102167712396,
                    -0.1590891264716532,
                    0.13233313109972625,
                    -0.030493374720953297,
                    -0.324131016869104,
                ],
                [
                    -0.44478515534248564,
                    -0.22066143889313278,
                    -0.31968023354915465,
                    0.025341827810811048,
                    -0.08457818649830318,
                    0.16910787319016476,
                    0.38203483276154415,
                    -0.10511631297976766,
                    0.08750071031831075,
                    0.2747264545695668,
                ],
            ],
            [
                [0.10508025998686506],
                [-0.2691991244124654],
                [0.21154928243402793],
                [0.6059699335423194],
                [-0.3520197843417094],
                [0.1473370411925756],
                [0.3586311702546054],
                [0.1428546809034874],
                [-0.2921184613487701],
                [0.20088114454313125],
            ],
        ],
        "regressor.intercepts_": [
            [
                0.5652480230183061,
                -0.14461022167946497,
                0.0734880322671333,
                0.774381965871222,
                -0.5810353861392875,
                -0.6447073648626959,
                0.6469258828055443,
                -0.8684424314182713,
                0.42510685510971963,
                -0.0666156988876174,
            ],
            [0.12043980294087006],
        ],
    },
}


def main_features(model, res):
    return [
        res.main_chain_correlation,
        1 if res.has_pepflip_peak else 0,
        max(atom.max_overlap for atom in res.main_chain_atoms),
        max(atom.bfactor_zscore for atom in res.main_chain_atoms),
        max(atom.bchange_zscore for atom in res.main_chain_atoms),
        mean([atom.density_zscore for atom in res.main_chain_atoms]),
        min(atom.diff_zscore for atom in res.main_chain_atoms),
        min(atom.density_zscore for atom in res.main_chain_atoms),
        None if res.next is None else res.next.atoms[" CA "].diff_zscore,
        res.ramachandran_score,
        model.resolution,
        res.twistedness,
    ]


def side_features(model, res):
    return [
        res.side_chain_correlation,
        max(atom.max_overlap for atom in res.side_chain_atoms),
        max(atom.bfactor_zscore for atom in res.side_chain_atoms),
        max(atom.bchange_zscore for atom in res.side_chain_atoms),
        mean([atom.density_zscore for atom in res.side_chain_atoms]),
        min(atom.diff_zscore for atom in res.side_chain_atoms),
        min(atom.density_zscore for atom in res.side_chain_atoms),
        model.resolution,
        res.rotamer_score,
    ]


def add_medians(X, medians):
    return [medians[i] if X[i] is None else X[i] for i in range(len(X))]


def scale(X, means, scales):
    return [(X[i] - means[i]) / scales[i] for i in range(len(X))]


def predict(X, coefs, intercepts):
    n_layers = len(coefs) + 1
    values = [None] * n_layers
    values[0] = list(X)
    for i in range(1, n_layers):
        values[i] = list(intercepts[i - 1])
        for j in range(len(intercepts[i - 1])):
            for k in range(len(values[i - 1])):
                values[i][j] += values[i - 1][k] * coefs[i - 1][k][j]
            if i < n_layers - 1:
                values[i][j] = math.tanh(values[i][j])
    return values[-1][0]


# ATOM / RESIDUE / CHAIN / MODEL


class Atom:
    def __init__(self, model, residue, atom_info):
        self.name = atom_info[0][0]
        self.alt_conf = atom_info[0][1]
        self.occupancy = atom_info[1][0]
        self.bfactor = atom_info[1][1]
        self.element = atom_info[1][2].strip()
        self.atomic_number = atomic_numbers[self.element]
        self.is_main_chain = self.name in main_chain_atoms
        self.is_side_chain = not self.is_main_chain
        self.point = atom_info[2]
        self.density = density_at_point(model.imap, *self.point)
        self.density_norm = self.density / self.atomic_number
        self.diff_density = density_at_point(model.imap_diff, *self.point)
        self.diff_norm = self.diff_density / self.atomic_number
        self.max_overlap = 0


class Residue:
    def __init__(self, model, spec, name):
        self.spec = spec
        self.name = name
        self.chain = spec[0]
        self.resno = spec[1]
        self.ins_code = spec[2]
        self.next = None
        self.prev = None
        self.twistedness = None
        self.ramachandran_score = None
        self.rotamer_score = rotamer_score(model.imol, spec[0], spec[1], spec[2], "")
        self.atoms = {}
        for atom_info in residue_info_py(model.imol, *spec):
            atom = Atom(model, self, atom_info)
            model.atoms.append(atom)
            self.atoms[atom.name] = atom
        self.main_chain_atoms = [atom for atom in self.atoms.values() if atom.is_main_chain]
        self.side_chain_atoms = [atom for atom in self.atoms.values() if atom.is_side_chain]
        self.truncatable = len(self.side_chain_atoms) > 0


class Chain:
    def __init__(self):
        self.residues = []
        self.correctnesses = []


class Model:
    def __init__(self, imol, imap, imap_diff):
        self.imol = imol
        self.imap = imap
        self.imap_diff = imap_diff
        self.resolution = data_resolution(self.imap)
        self.residues = []
        self.residue_spec_dict = {}
        self.atoms = []
        for ichain in range(n_chains(imol)):
            chain_id = chain_id_py(imol, ichain)
            for serial_num in range(chain_n_residues(chain_id, imol)):
                resno = seqnum_from_serial_number(imol, chain_id, serial_num)
                ins_code = insertion_code_from_serial_number(imol, chain_id, serial_num)
                spec = [chain_id, resno, ins_code]
                name = residue_name(imol, *spec)
                if is_protein(imol, name, spec):
                    residue = Residue(self, spec, name)
                    self.residues.append(residue)
                    self.residue_spec_dict[tuple(spec)] = residue
        self.set_connections()
        self.set_correlations()
        self.set_bond_changes()
        self.set_zscores()
        self.set_ramachandran()
        self.set_overlaps()
        self.set_pepflip_peaks()
        self.set_correctness()
        self.set_chains()

    def set_connections(self):
        for i in range(len(self.residues) - 1):
            res1 = self.residues[i]
            res2 = self.residues[i + 1]
            if " CA " not in res1.atoms or " C  " not in res1.atoms:
                continue
            if " N  " not in res2.atoms or " CA " not in res2.atoms:
                continue
            point1 = res1.atoms[" C  "].point
            point2 = res2.atoms[" N  "].point
            if distance(point1, point2) < 1.7:
                res1.next = res2
                res2.prev = res1
                xyz1 = res1.atoms[" CA "].point
                xyz2 = res1.atoms[" C  "].point
                xyz3 = res2.atoms[" N  "].point
                xyz4 = res2.atoms[" CA "].point
                omega = abs(torsion(xyz1, xyz2, xyz3, xyz4))
                twistedness = min(omega, 180 - omega)
                if twistedness < 0:
                    raise Exception("Negative twistedness between %s and %s" % (res1.spec, res2.spec))
                res1.twistedness = max(twistedness, 0 if res1.twistedness is None else res1.twistedness)
                res2.twistedness = max(twistedness, 0 if res2.twistedness is None else res2.twistedness)

    def set_correlations(self):
        self.set_correlation("main_chain_correlation", 1)
        self.set_correlation("side_chain_correlation", 3)

    def set_bond_changes(self):
        for residue in self.residues:
            for atom in residue.atoms.values():
                atom.bchange = -99999
                for attached in attached_atoms(atom, residue):
                    bchange = (atom.bfactor - attached.bfactor) / attached.bfactor
                    atom.bchange = max(bchange, atom.bchange)
                if atom.bchange == -99999:
                    raise Exception("Residue %s atom %s is not bonded to anything" % (residue.spec, atom.name))

    def set_correlation(self, attr, mask):
        specs = [residue.spec for residue in self.residues]
        for correlation in map_to_model_correlation_per_residue_py(self.imol, specs, mask, self.imap):
            residue = self.residue_spec_dict[tuple(correlation[0][1:])]
            setattr(residue, attr, correlation[1])

    # Uses a modified Z-score calculated using medians
    # Boris Iglewicz and David Hoaglin
    # Volume 16: How to Detect and Handle Outliers
    # The ASQC Basic References in Quality Control: Statistical Techniques
    # Edited by Edward F. Mykytka, Ph.D.
    # Page 11
    # 1993
    def set_zscores(self):
        def set_zscore(attr, zscore_attr):
            main_values = [getattr(atom, attr) for atom in self.atoms if atom.is_main_chain]
            side_values = [getattr(atom, attr) for atom in self.atoms if atom.is_side_chain]
            main_median = median(main_values)
            side_median = median(side_values)
            main_mad = median_absolute_deviation(main_values)
            side_mad = median_absolute_deviation(side_values)
            for residue in self.residues:
                for atom in residue.atoms.values():
                    value = getattr(atom, attr)
                    if atom.is_main_chain:
                        zscore = 0 if main_mad == 0 else 0.6745 * (value - main_median) / main_mad
                    else:
                        zscore = 0 if side_mad == 0 else 0.6745 * (value - side_median) / side_mad
                    setattr(atom, zscore_attr, zscore)

        set_zscore("density_norm", "density_zscore")
        set_zscore("diff_norm", "diff_zscore")
        set_zscore("bfactor", "bfactor_zscore")
        set_zscore("bchange", "bchange_zscore")

    def set_ramachandran(self):
        for item in all_molecule_ramachandran_score_py(self.imol)[5]:
            spec = item[1][1:]
            if tuple(spec) in self.residue_spec_dict:
                residue = self.residue_spec_dict[tuple(spec)]
                residue.ramachandran_score = item[2]

    def set_overlaps(self):
        for overlap in molecule_atom_overlaps_py(self.imol):
            for n in (1, 2):
                spec = tuple(overlap["atom-%d-spec" % n][1:4])
                if spec in self.residue_spec_dict:
                    atom_name = overlap["atom-%d-spec" % n][4]
                    atom = self.residue_spec_dict[spec].atoms[atom_name]
                    atom.max_overlap = max(atom.max_overlap, overlap["overlap-volume"])

    def set_pepflip_peaks(self):
        peaks = map_peaks_around_molecule_py(self.imap_diff, 4.62567528, False, self.imol)

        def has_o_moving_peak(residue):
            if residue.next is None:
                return False
            if " O  " not in residue.atoms:
                return False
            ca1 = residue.atoms[" CA "].point
            c = residue.atoms[" C  "].point
            o = residue.atoms[" O  "].point
            ca2 = residue.next.atoms[" CA "].point
            for peak in peaks:
                point = peak[1]
                cd = distance(point, c)
                if cd < 1.09769614 or cd > 3.22217608:
                    continue
                if angle(o, c, point) < 75.45064459:
                    continue
                ca1d = distance(point, ca1)
                if ca1d < 0.66596873 or ca1d > 2.8910501:
                    continue
                ca2d = distance(point, ca2)
                if ca2d < 2.23908466 or ca2d > 3.84813939:
                    continue
                return True
            return False

        def has_n_moving_peak(residue):
            if residue.prev is None or not residue.prev.has_pepflip_peak:
                return False
            if " O  " not in residue.prev.atoms:
                return False
            o = residue.prev.atoms[" O  "].point
            ca = residue.atoms[" CA "].point
            for peak in peaks:
                point = peak[1]
                if distance(point, o) > 0.52516368:
                    continue
                if distance(point, ca) > 2.76397302:
                    continue
                return True
            return False

        for residue in self.residues:
            residue.has_pepflip_peak = has_o_moving_peak(residue) or has_n_moving_peak(residue)

    def set_correctness(self):
        if len(training_data) == 0:
            return
        for res in self.residues:
            X_main = main_features(self, res)
            X_main = add_medians(X_main, training_data["main"]["medians"])
            X_main = scale(X_main, training_data["main"]["scaler.mean_"], training_data["main"]["scaler.scale_"],)
            res.main_chain_correctness = predict(
                X_main, training_data["main"]["regressor.coefs_"], training_data["main"]["regressor.intercepts_"],
            )
            if res.truncatable:
                X_side = side_features(self, res)
                X_side = add_medians(X_side, training_data["side"]["medians"])
                X_side = scale(X_side, training_data["side"]["scaler.mean_"], training_data["side"]["scaler.scale_"],)
                res.side_chain_correctness = predict(
                    X_side, training_data["side"]["regressor.coefs_"], training_data["side"]["regressor.intercepts_"],
                )

    def set_chains(self):
        self.chains = {}
        for residue in self.residues:
            chain_id = residue.chain
            if chain_id not in self.chains:
                self.chains[chain_id] = Chain()
            chain = self.chains[chain_id]
            chain.residues.append(residue)
            if hasattr(residue, "main_chain_correctness"):
                chain.correctnesses.append(residue.main_chain_correctness)
        for chain in self.chains.values():
            if len(chain.correctnesses) > 0:
                chain.correctness = mean(chain.correctnesses)


# GUI


class TodoListWindow(gtk.Window):
    instance = None

    def __init__(self, title, categories, height=600):
        gtk.Window.__init__(self, gtk.WINDOW_TOPLEVEL)
        self.check_input(categories)
        self.set_title(title)
        self.selected_todo = None
        if TodoListWindow.instance is not None:
            TodoListWindow.instance.destroy()
        TodoListWindow.instance = self
        self.init_content(categories)
        self.show_all()
        self.set_size(height)

    def check_input(self, categories):
        for category in categories:
            assert "title" in category
            assert "columns" in category
            assert len(category["columns"]) > 0
            assert "todos" in category
            for todo in category["todos"]:
                assert all(column in todo for column in category["columns"])
                assert "point" in todo

    def init_content(self, categories):
        self.outside_vbox = gtk.VBox(spacing=2)
        self.outside_vbox.set_border_width(2)
        self.add(self.outside_vbox)
        self.scrolled_window = gtk.ScrolledWindow()
        self.scrolled_window.set_policy(gtk.POLICY_NEVER, gtk.POLICY_NEVER)
        self.outside_vbox.pack_start(self.scrolled_window, True, True, 0)
        self.inside_vbox = gtk.VBox()
        self.scrolled_window.add_with_viewport(self.inside_vbox)
        self.frames = []
        for category in categories:
            frame = TodoListFrame(self, category)
            self.frames.append(frame)
            self.inside_vbox.pack_start(frame, False, False, 5)
        self.close_button = gtk.Button("Close")
        self.close_button.connect("clicked", lambda b: self.destroy())
        self.outside_vbox.pack_end(self.close_button, False, False, 0)

    def set_size(self, height):
        max_frame_height = height - 150
        self.scrolled_window.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        size = self.get_size()
        width = size[0] + 50
        height = min(size[1] + 2, height)
        self.resize(width, height)
        for frame in self.frames:
            if len(frame.todos) == 0:
                continue
            frame.scrolled_window.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
            height = 30 + 26 * len(frame.todos)
            if height > max_frame_height:
                height = max_frame_height
            frame.scrolled_window.set_size_request(width - 50, height)
        self.show_all()


class TodoListFrame(gtk.Frame):
    def __init__(self, window, category):
        gtk.Frame.__init__(self)
        self.parent_window = window
        self.title = category["title"]
        self.todos = category["todos"]
        self.columns = category["columns"]
        self.expanded = True
        self.set_border_width(6)

        title_event_box = gtk.EventBox()
        self.title_label = gtk.Label()
        self.set_title_label()
        title_event_box.add(self.title_label)
        title_event_box.connect("button-press-event", self.expand_contract)
        self.set_label_widget(title_event_box)

        self.contents = gtk.VBox()
        if "description" in category:
            description = gtk.Label(category["description"])
            description.set_line_wrap(True)
            self.contents.pack_start(description, False, False, 5)
        if len(self.todos) > 0:
            if "widget" in category:
                self.contents.pack_start(category["widget"], False, False, 5)
            self.init_list_store()
            self.init_tree_view()
            self.add_tree_view_with_scrolled_window()
        self.add(self.contents)

    def expand_contract(self, event_box, event):
        self.expanded = not self.expanded
        if self.expanded:
            self.contents.show()
        else:
            self.contents.hide()
        self.set_title_label()

    def set_title_label(self):
        button = "-" if self.expanded else "+"
        done = 0 if not hasattr(self, "list_store") else [row[self.done_col] for row in self.list_store].count(True)
        text = "%s %s (%d/%d)" % (button, self.title, done, len(self.todos))
        self.title_label.set_label(text)

    def init_list_store(self):
        types = [type(self.todos[0][name]) for name in self.columns] + [
            bool,
            object,
            object,
        ]
        self.done_col = len(self.columns)
        self.xyz_col = len(self.columns) + 1
        self.todo_col = len(self.columns) + 2
        self.list_store = gtk.ListStore(*types)
        for todo in self.todos:
            todo["category"] = self.title
            row = [todo[name] for name in self.columns] + [False, todo["point"], todo]
            self.list_store.append(row)

    def init_tree_view(self):
        self.tree_view = gtk.TreeView(self.list_store)
        for i, column_name in enumerate(self.columns):
            renderer = gtk.CellRendererText()
            column = gtk.TreeViewColumn(column_name, renderer)
            column.add_attribute(renderer, "text", i)
            column.add_attribute(renderer, "strikethrough", self.done_col)
            column.set_expand(True)
            column.set_sort_column_id(i)
            self.tree_view.append_column(column)
        renderer = gtk.CellRendererToggle()
        renderer.connect("toggled", self.done_toggled)
        column = gtk.TreeViewColumn("Done", renderer)
        column.add_attribute(renderer, "active", self.done_col)
        column.set_expand(True)
        column.set_sort_column_id(self.done_col)
        self.tree_view.append_column(column)
        selection = self.tree_view.get_selection()
        selection.set_mode(gtk.SELECTION_SINGLE)
        self.tree_view.connect("cursor-changed", self.on_tree_view_selection)
        self.tree_view.connect("button-release-event", self.button_released)
        self.done_toggled_this_click = False

    def add_tree_view_with_scrolled_window(self):
        self.scrolled_window = gtk.ScrolledWindow()
        self.scrolled_window.set_policy(gtk.POLICY_NEVER, gtk.POLICY_NEVER)
        self.scrolled_window.add(self.tree_view)
        self.contents.pack_start(self.scrolled_window, False, False, 5)

    def done_toggled(self, widget, path):
        self.list_store[path][self.done_col] = not self.list_store[path][self.done_col]
        self.set_title_label()
        self.done_toggled_this_click = True

    def button_released(self, *args):
        self.done_toggled_this_click = False

    def on_tree_view_selection(self, tree_view):
        if self.done_toggled_this_click:
            return
        selection = tree_view.get_selection()
        model, treeiter = selection.get_selected()
        if treeiter is not None:
            xyz = model[treeiter][self.xyz_col]
            todo = model[treeiter][self.todo_col]
            set_rotation_centre(*xyz)
            for frame in self.parent_window.frames:
                if self == frame or len(frame.todos) == 0:
                    continue
                frame.tree_view.get_selection().unselect_all()
            self.parent_window.selected_todo = todo


def chain_todos(model):
    todos = []
    for chain_id in sorted(model.chains):
        chain = model.chains[chain_id]
        todos.append(
            {
                "Chain": chain_id,
                "Length": len(chain.residues),
                "UNKs": len([r for r in chain.residues if r.name == "UNK"]),
                "Correctness": chain.correctness,
                "point": model.chains[chain_id].residues[0].atoms[" CA "].point,
            }
        )
    todos = sorted(todos, key=lambda todo: todo["Correctness"])
    return todos


def mainchain_todos(model):
    todos = []
    for residue in model.residues:
        todos.append(
            {
                "Chain": residue.chain,
                "Residue": residue.resno,
                "Ins": residue.ins_code,
                "Name": residue.name,
                "Correctness": residue.main_chain_correctness,
                "point": residue.atoms[" CA "].point,
            }
        )
    todos = sorted(todos, key=lambda todo: todo["Correctness"])
    return todos


def sidechain_todos(model):
    todos = []
    for residue in model.residues:
        if residue.truncatable:
            todos.append(
                {
                    "Chain": residue.chain,
                    "Residue": residue.resno,
                    "Ins": residue.ins_code,
                    "Name": residue.name,
                    "Correctness": residue.side_chain_correctness,
                    "point": residue.atoms[" CB "].point,
                }
            )
    todos = sorted(todos, key=lambda todo: todo["Correctness"])
    return todos


def fill_option_menu_with_mol_options(menu, mol_nos):
    for mol_no in mol_nos:
        label = "%s %s" % (mol_no, molecule_name(mol_no))
        menu.append_text(label)
        menu.set_active(0)


def imolmapdiff_chooser(callback):
    imols = []
    imaps = []
    imap_diffs = []
    for i in molecule_number_list():
        if valid_model_molecule_qm(i):
            imols.append(i)
        elif valid_map_molecule_qm(i):
            if map_is_difference_map(i):
                imap_diffs.append(i)
            else:
                imaps.append(i)
    if len(imols) == 0 or len(imaps) == 0 or len(imap_diffs) == 0:
        print("ML_CORRECTNESS: Can't use ML Correctness: " "must have a valid model, map and difference map")
        return

    if len(imols) == 1 and len(imaps) == 1 and len(imap_diffs) == 1:
        callback(imols[0], imaps[0], imap_diffs[0])
        return

    def destroy_window(*args):
        window.destroy()
        return False

    def on_ok_clicked(*args):
        imol = int(get_option_menu_active_molecule(imol_options, imols))
        imap = int(get_option_menu_active_molecule(imap_options, imaps))
        imap_diff = int(get_option_menu_active_molecule(imap_diff_options, imap_diffs))
        callback(imol, imap, imap_diff)
        destroy_window()

    window = gtk.Window(gtk.WINDOW_TOPLEVEL)
    imol_options = gtk.combo_box_new_text()
    fill_option_menu_with_mol_options(imol_options, imols)
    imap_options = gtk.combo_box_new_text()
    fill_option_menu_with_mol_options(imap_options, imaps)
    imap_diff_options = gtk.combo_box_new_text()
    fill_option_menu_with_mol_options(imap_diff_options, imap_diffs)
    ok_button = gtk.Button("  OK  ")
    ok_button.connect("clicked", on_ok_clicked)
    cancel_button = gtk.Button(" Cancel ")
    cancel_button.connect("clicked", destroy_window)

    window.set_title("Coot")
    vbox = gtk.VBox(False, 6)
    window.add(vbox)
    vbox.pack_start(gtk.Label("Model"), False, False, 5)
    vbox.pack_start(imol_options, True, True, 0)
    vbox.pack_start(gtk.Label("Map"), False, False, 5)
    vbox.pack_start(imap_options, True, True, 0)
    vbox.pack_start(gtk.Label("Difference Map"), False, False, 5)
    vbox.pack_start(imap_diff_options, True, True, 0)
    vbox.pack_start(gtk.HSeparator(), True, False, 2)
    hbox_buttons = gtk.HBox(False, 5)
    vbox.pack_start(hbox_buttons, False, False, 5)
    hbox_buttons.pack_start(ok_button, True, False, 5)
    hbox_buttons.pack_start(cancel_button, True, False, 5)
    window.show_all()


def show_gui(imol, imap, imap_diff):
    model = Model(imol, imap, imap_diff)
    TodoListWindow(
        "ML Correctness",
        [
            {"title": "Chains", "columns": ["Chain", "Length", "UNKs", "Correctness"], "todos": chain_todos(model),},
            {
                "title": "Main chains",
                "columns": ["Chain", "Residue", "Ins", "Name", "Correctness"],
                "todos": mainchain_todos(model),
            },
            {
                "title": "Side chains",
                "columns": ["Chain", "Residue", "Ins", "Name", "Correctness"],
                "todos": sidechain_todos(model),
            },
        ],
    )


if "coot_menubar_menu" in locals():
    menu = coot_menubar_menu("ML Correctness")
    add_simple_coot_menu_menuitem(menu, "GUI", lambda x: imolmapdiff_chooser(show_gui))
    add_simple_coot_menu_menuitem(menu, "Prune", lambda x: imolmapdiff_chooser(prune))


# SCRIPTING


def prune(
    imol,
    imap,
    imap_diff,
    chains=True,
    chain_threshold="auto",
    max_chain_fraction=0.2,
    max_chain_length=20,
    residues=True,
    residue_threshold="auto",
    max_residue_fraction=0.2,
    remove_isolated_residues=True,
    sidechains=True,
    sidechain_threshold="auto",
    max_sidechain_fraction=0.2,
):

    model = Model(imol, imap, imap_diff)

    if len(model.chains) < 1:
        return

    if chains:
        main_median = median([r.main_chain_correctness for r in model.residues])
        if chain_threshold == "auto":
            chain_threshold = main_median * 0.2
        print(
            "ML_CORRECTNESS: Deleting chains (up to %d residues long) with scores < %.3f"
            % (max_chain_length, chain_threshold)
        )
        print("ML_CORRECTNESS: Up to %.0f%% of residues will be deleted" % (max_chain_fraction * 100))
        max_deleted = len(model.residues) * max_chain_fraction
        deleted = 0
        remaining = []
        chain_ids = model.chains.keys()
        chain_ids.sort(key=lambda cid: model.chains[cid].correctness)
        for chain_id in sorted(model.chains,):
            chain = model.chains[chain_id]
            if (
                chain.correctness < chain_threshold
                and len(chain.residues) <= max_chain_length
                and deleted + len(chain.residues) <= max_deleted
            ):
                deleted += len(chain.residues)
                delete_chain(imol, chain_id)
            else:
                remaining.extend([r for r in chain.residues])
        print("ML_CORRECTNESS: Deleted %.0f%% of residues" % (float(deleted) / len(model.residues) * 100))
    else:
        remaining = [r for r in model.residues]

    if len(remaining) < 1:
        return

    if residues:
        main_median = median([r.main_chain_correctness for r in remaining])
        if residue_threshold == "auto":
            residue_threshold = main_median * 0.5
        print("ML_CORRECTNESS: Deleting residues with scores < %.3f" % residue_threshold)
        print("ML_CORRECTNESS: Up to %.0f%% of residues will be deleted" % (max_residue_fraction * 100))
        max_deleted = len(remaining) * max_residue_fraction
        deleted = 0
        remaining.sort(key=lambda r: r.main_chain_correctness)
        for residue in remaining:
            residue.delete = False
            if residue.main_chain_correctness < residue_threshold and deleted + 1 <= max_deleted:
                deleted += 1
                residue.delete = True
        if remove_isolated_residues:
            for residue in remaining:
                if residue.prev is not None and not residue.prev.delete:
                    continue
                if residue.next is not None and not residue.next.delete:
                    continue
                residue.delete = True
        for residue in remaining:
            if residue.delete:
                delete_residue(imol, residue.chain, residue.resno, residue.ins_code)
        print("ML_CORRECTNESS: Deleted %.0f%% of residues" % (float(deleted) / len(remaining) * 100))
        remaining = [r for r in remaining if not r.delete]

    remaining = [r for r in remaining if r.truncatable]
    if len(remaining) < 1:
        return

    if sidechains:
        side_median = median([r.side_chain_correctness for r in model.residues if r.truncatable])
        if sidechain_threshold == "auto":
            sidechain_threshold = side_median * 0.5
        print("ML_CORRECTNESS: Deleting sidechains with scores < %.3f" % sidechain_threshold)
        print("ML_CORRECTNESS: Up to %.0f%% of sidechains will be deleted" % (max_sidechain_fraction * 100))
        max_deleted = len(remaining) * max_sidechain_fraction
        deleted = 0
        remaining.sort(key=lambda r: r.side_chain_correctness)
        for residue in remaining:
            if residue.side_chain_correctness < sidechain_threshold and deleted + 1 < max_deleted:
                deleted += 1
                delete_residue_sidechain(imol, residue.chain, residue.resno, residue.ins_code, 0)
        print("ML_CORRECTNESS: Deleted %.0f%% of sidechains" % (float(deleted) / len(remaining) * 100))
