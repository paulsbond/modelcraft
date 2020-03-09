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
    " C  ": {" CA ", " O  "},
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
            0.8139716982841492,
            0.0,
            0.376086983983169,
            0.2961565371112583,
            1.160100205235751,
            0.08482286142947047,
            -1.0453192991086182,
            -0.7553121220334733,
            0.09134544293313959,
            0.2097645763533924,
            1.9984469413757324,
            7.492323662130829,
        ],
        "scaler.mean_": [
            0.7682692677233491,
            0.018148262073208243,
            0.4627909950409583,
            0.9245502611767549,
            1.396284058678193,
            -0.08631770961879058,
            -1.1642906477599513,
            -0.9168867599718515,
            0.06389955815807204,
            0.7460009577500388,
            1.9125921660594094,
            8.994721429841155,
        ],
        "scaler.scale_": [
            0.15444164528082477,
            0.1334874625460613,
            0.48816667792871743,
            2.2783615459089113,
            1.3486772253140358,
            0.8937534379089717,
            0.9239649914154019,
            0.9883603740887346,
            1.2025388394669327,
            1.1802774937250036,
            0.42224754623855043,
            6.809482541249054,
        ],
        "regressor.coefs_": [
            [
                [
                    0.6143797225161832,
                    -0.4595911301115578,
                    -0.02930740287226514,
                    -0.22043430485424026,
                    0.6925544323073416,
                    -0.1240644004159779,
                    0.35133445333519747,
                    -0.5310353304851705,
                    -0.020760906270481615,
                    0.07105901733613063,
                ],
                [
                    -0.19893530700279707,
                    -1.1199067651432602,
                    0.1886375723433218,
                    1.7142328018778785,
                    -0.0639971071085352,
                    1.1116986422124884,
                    -0.8349332063732631,
                    -0.7767166481879676,
                    -0.730715880162232,
                    0.6897588964360177,
                ],
                [
                    -0.34030084051504816,
                    0.0031671266162438573,
                    0.1101730714920676,
                    0.03567228543679168,
                    -0.07233438064206893,
                    -0.7540676053702117,
                    -0.07616736907184114,
                    0.04900731148923605,
                    0.013127437681358647,
                    -0.27217307215835657,
                ],
                [
                    -0.20975295105804911,
                    0.46918100222045,
                    0.015372781420030101,
                    -1.8322384275746633,
                    -0.46163312919940136,
                    -0.018522110001506335,
                    -0.16642043565064057,
                    0.0743249719990909,
                    -0.7715748925000463,
                    0.30689399390613087,
                ],
                [
                    -0.10464175668790789,
                    0.06419287381523431,
                    0.4335222050835248,
                    -0.16942663370103975,
                    -0.05055983875688424,
                    0.03683511080108421,
                    0.16786586236183793,
                    -0.11200247200801353,
                    -0.25132157890694945,
                    -0.12363521385055494,
                ],
                [
                    0.6501423456690897,
                    0.47664657871792687,
                    0.5103686095958361,
                    0.020861321958729973,
                    -0.05955097514669694,
                    -0.1566780965099033,
                    0.5038785808366654,
                    -0.0894798702307193,
                    0.20520681798291154,
                    -0.22582766589015196,
                ],
                [
                    -0.28837541278466655,
                    -0.22033946038341554,
                    -0.20250922212493583,
                    -0.16551447958569931,
                    0.3431004812283165,
                    0.4056541697142576,
                    -0.015808464020352297,
                    -0.042568713459129943,
                    -0.1433957182883911,
                    -0.16405947628808087,
                ],
                [
                    0.2902764745189481,
                    0.19132290301121183,
                    -0.9924721524548433,
                    0.3008648018166311,
                    0.032584037954498234,
                    -0.033641569015730975,
                    -0.1362323766300239,
                    -0.033570880319530046,
                    0.27392029112924,
                    0.25929830981232416,
                ],
                [
                    0.10789591340142685,
                    -0.01413157365996538,
                    0.15121010328795068,
                    -0.07054351489003059,
                    0.1303962686385081,
                    0.10015053630408653,
                    0.041271069010533075,
                    -0.04387673145656446,
                    -0.008145821346669947,
                    0.4945857079704406,
                ],
                [
                    1.2377128494505736,
                    0.024961170838804177,
                    -0.1678390549338977,
                    0.07298633716304263,
                    -0.04206106343619906,
                    -0.26772945250038577,
                    0.09801041147674049,
                    -0.12275255833851273,
                    0.12325247402531916,
                    -0.09880239229665681,
                ],
                [
                    -0.5934746292165817,
                    -0.576961278525674,
                    0.21566938934920554,
                    0.4139406979054004,
                    -0.014123124282201308,
                    0.4119043723112296,
                    -0.5839758598335055,
                    0.4677504262620488,
                    -0.017057276635181653,
                    0.28312861599199796,
                ],
                [
                    -0.13380497642314856,
                    0.002127976889487308,
                    0.23894067942569522,
                    -0.1578335548621471,
                    -0.12137082335145756,
                    -0.10407232275799869,
                    -0.040141602986115375,
                    0.03753865802667223,
                    -0.13877330959772186,
                    -0.046939567351882844,
                ],
            ],
            [
                [0.2582327040700272],
                [0.5359040278244972],
                [0.19131851495722435],
                [-0.12527964181856077],
                [0.2642007690847596],
                [0.0807647408516691],
                [-0.6927601013843842],
                [-0.6300596409527701],
                [0.35332831997657466],
                [-0.08756671193143149],
            ],
        ],
        "regressor.intercepts_": [
            [
                1.5518483724369763,
                1.1052285111036375,
                0.3270333365902817,
                0.6294059472035622,
                0.9312192333310708,
                0.17965543050772695,
                -0.33221202370473135,
                0.045338838976720405,
                -0.22189212726349586,
                -0.20430867462391178,
            ],
            [0.02932469458860332],
        ],
    },
    "side": {
        "medians": [
            0.7164410948753357,
            0.38994837045904174,
            0.12586781857451385,
            0.7630352303851361,
            -0.011975465875031385,
            -0.7592286410664458,
            -0.4993887613407643,
            1.939817190170288,
            37.424095153808594,
        ],
        "scaler.mean_": [
            0.673412019631825,
            0.48089143726991673,
            0.5799702558160009,
            1.0143532611843056,
            -0.0942555396903983,
            -0.8750788889815784,
            -0.5910028207303396,
            1.890088834582728,
            41.473247390709204,
        ],
        "scaler.scale_": [
            0.1906889449197939,
            0.5125132538979486,
            1.5721257494386607,
            1.612785396142665,
            0.8565831079891278,
            1.061990724043148,
            0.9495461490877636,
            0.4218492745313629,
            29.93087152763042,
        ],
        "regressor.coefs_": [
            [
                [
                    0.03961428978473987,
                    0.2873234170392773,
                    0.1809070597580511,
                    -0.20559134559419323,
                    0.7496460978584658,
                    0.25093073518918985,
                    -0.09235170237724591,
                    0.09444160834378032,
                    -0.21072106804793264,
                    -0.11185320439905218,
                ],
                [
                    0.0228075211308879,
                    0.08889494354648758,
                    -0.008493270151176265,
                    0.3177943717900223,
                    -0.13355332982773763,
                    0.018307658993131404,
                    -0.20793890282336028,
                    -0.019026844012020018,
                    0.015933914847746197,
                    0.04946836199564234,
                ],
                [
                    0.29626262813956095,
                    0.616573842486694,
                    -0.19688956127956922,
                    -0.013137276121410331,
                    -0.8350741584004115,
                    0.127626247329841,
                    0.44329062361276783,
                    0.22421546470226758,
                    0.5541506833104448,
                    0.34050195359364827,
                ],
                [
                    -0.1800042790080773,
                    0.3694559663725238,
                    -0.14862565646989984,
                    -0.06825959299702655,
                    -0.22078147887506677,
                    0.17896025956530875,
                    -0.1251541708425913,
                    -0.08522341375747354,
                    -0.152347043746884,
                    0.11104038982523914,
                ],
                [
                    -0.0034331711633664584,
                    -0.27932806135809596,
                    -0.15607868808539002,
                    -0.10495217223737632,
                    -0.5117721494069597,
                    -0.38409920278602394,
                    0.13711189122185688,
                    -0.06064382384505531,
                    0.2665817092212631,
                    0.1346268601868371,
                ],
                [
                    0.13502361437191668,
                    0.21924538563236448,
                    -0.04363224539330135,
                    -0.02058677765553592,
                    0.12305897990753915,
                    0.05653349670836186,
                    -0.18316250091849964,
                    0.02665674276088463,
                    0.050221774186613374,
                    0.03409505540888545,
                ],
                [
                    0.17315481036230357,
                    -0.10506312341139304,
                    0.09818022226456279,
                    0.5597050106203582,
                    0.5416523997580623,
                    0.10760004710641499,
                    -0.07296380147467701,
                    0.12232100562595896,
                    0.021382828678610338,
                    -0.05552024142428309,
                ],
                [
                    0.07493031477495755,
                    0.11229919917752604,
                    -0.09049093431715949,
                    0.2831429766241456,
                    -0.23045942509055142,
                    0.46473417641509623,
                    -0.09686553125722495,
                    0.16561548876851917,
                    0.034398478009819135,
                    -0.11298031818118966,
                ],
                [
                    0.04958093989745634,
                    -0.054652224582197456,
                    -0.5953898290996007,
                    -0.190337322605491,
                    0.31714922174700594,
                    -0.17952610440780983,
                    -0.20064575191072623,
                    -0.012640205995514405,
                    -0.11390468629134977,
                    -0.12346418496893617,
                ],
            ],
            [
                [0.7231591970641398],
                [-0.3748327992519657],
                [-0.144522702759284],
                [-0.08626418815345917],
                [0.44029785379560354],
                [0.3668058253213245],
                [0.1408607921916091],
                [-0.43900584998226844],
                [-0.5158909840566932],
                [0.6449625117582403],
            ],
        ],
        "regressor.intercepts_": [
            [
                0.5484044591733884,
                0.5585838464276648,
                -0.836270022521484,
                0.38594900804597293,
                0.40515818757604916,
                0.644181609423902,
                -0.012399826608033246,
                -0.561202984549227,
                -0.31099363666935415,
                0.1823655237716331,
            ],
            [-0.36766606249515704],
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
                    raise Exception("Atom %s is not bonded to anything" % atom.name)

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
        peaks = map_peaks_around_molecule_py(self.imap_diff, 4.44722675, False, self.imol)

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
                if cd < 0.8867185 or cd > 2.75353879:
                    continue
                if angle(o, c, point) < 60.94925295:
                    continue
                ca1d = distance(point, ca1)
                if ca1d < 1.01077907 or ca1d > 3.70564479:
                    continue
                ca2d = distance(point, ca2)
                if ca2d < 1.84499095 or ca2d > 3.87123115:
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
                if distance(point, o) > 1.45724775:
                    continue
                if distance(point, ca) > 2.09323668:
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
