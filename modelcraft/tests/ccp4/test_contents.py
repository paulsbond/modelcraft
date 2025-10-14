from pytest import approx

from modelcraft.contents import AsuContents, Polymer, PolymerType


def _test_contents(entry: str, expected_json: list, selenomet: bool):
    contents = AsuContents.from_pdbe(entry)
    assert contents.to_json() == expected_json
    assert contents.is_selenomet() == selenomet
    return contents


def test_1o6a():
    expected = {
        "copies": 2,
        "proteins": [
            {
                "sequence": (
                    "SETRKTEVPSDKLELLLDIPLKVTVELGRTRMTLKRVLEMIHGSIIELDKLTGEPVDILV"
                    "NGKLIARGEVVVIDENFGVRITEIVSPKERLELLNE"
                ),
                "stoichiometry": 1,
                "modifications": ["M->MSE"],
            }
        ],
        "rnas": [],
        "dnas": [],
        "carbs": [],
        "ligands": [],
        "buffers": [],
    }
    _test_contents("1o6a", expected, selenomet=True)


def test_4gxy():
    expected = {
        "copies": 1,
        "proteins": [],
        "rnas": [
            {
                "sequence": (
                    "GGCGGCAGGUGCUCCCGACCCUGCGGUCGGGAGUUAAAAGGGAAGCCGGUGCAAGUCCGG"
                    "CACGGUCCCGCCACUGUGACGGGGAGUCGCCCCUCGGGAUGUGCCACUGGCCCGAAGGCC"
                    "GGGAAGGCGGAGGGGCGGCGAGGAUCCGGAGUCAGGAAACCUGCCUGCCGUC"
                ),
                "stoichiometry": 1,
                "modifications": ["1->GTP", "172->CCC"],
            }
        ],
        "dnas": [],
        "carbs": [],
        "ligands": [
            {"code": "B1Z", "stoichiometry": 2},
            {"code": "IRI", "stoichiometry": 7},
        ],
        "buffers": ["MG"],
    }
    _test_contents("4gxy", expected, selenomet=False)


def test_6as7():
    expected = {
        "copies": 1,
        "proteins": [
            {
                "sequence": (
                    "DEEQVFHFYWLDAYEDQYNQPGVVFLFGKVWIESAETHVSCCVMVKNIERTLYFLPREMK"
                    "IDLNTGKETGTPISMKDVYEEFDEKIATKYKIMKFKSKPVEKNYAFEIPDVPEKSEYLEV"
                    "KYSAEMPQLPQDLKGETFSHVFGTNTSSLELFLMNRKIKGPCWLEVKSPQLLNQPVSWCK"
                    "AEAMALKPDLVNVIKDVSPPPLVVMAFSMKTMQNAKNHQNEIIAMAALVHHSFALDKAAPK"
                    "PPFQSHFCVVSKPKDCIFPYAFKEVIEKKNVKVEVAATERTLLGFFLAKVHKIDPDIIVGH"
                    "NIYGFELEVLLQRINVCKAPHWSKIGRLKRSNMPKLGGRSGFGERNATCGRMICDVEISAK"
                    "ELIRCKSYHLSELVQQILKTERVVIPMENIQNMYSESSQLLYLLEHTWKDAKFILQIMCEL"
                    "NVLPLALQITNIAGNIMSRTLMGGRSERNEFLLLHAFYENNYIVPDKQIFRKPQQKLGDED"
                    "EEIDGDTNKYKKGRKKAAYAGGLVLDPKVGFYDKFILLLDFNSLYPSIIQEFNICFTTVQR"
                    "VASEAQKVTEDGEQEQIPELPDPSLEMGILPREIRKLVERRKQVKQLMKQQDLNPDLILQY"
                    "DIRQKALKLTANSMYGCLGFSYSRFYAKPLAALVTYKGREILMHTKEMVQKMNLEVIYGDT"
                    "DSIMINTNSTNLEEVFKLGNKVKSEVNKLYKLLEIDIDGVFKSLLLLKKKKYAALVVEPTS"
                    "DGNYVTKQELKGLDIVRRDWCDLAKDTGNFVIGQILSDQSRDTIVENIQKRLIEIGENVLN"
                    "GSVPVSQFEINKALTKDPQDYPDKKSLPHVHVALWINSQGGRKVKAGDTVSYVICQDGSNL"
                    "TASQRAYAPEQLQKQDNLTIDTQYYLAQQIHPVVARICEPIDGIDAVLIATWLGLDPTQFR"
                    "VHHYHKDEEN"
                ),
                "stoichiometry": 1,
                "modifications": [],
            }
        ],
        "rnas": [],
        "dnas": [
            {
                "sequence": "GCCTGGAGCGC",
                "stoichiometry": 1,
                "modifications": [],
            },
            {
                "sequence": "AGGCGCTCCAGGC",
                "stoichiometry": 1,
                "modifications": [],
            },
        ],
        "carbs": [],
        "ligands": [{"code": "DCP", "stoichiometry": 1}],
        "buffers": ["MG", "CO"],
    }
    _test_contents("6as7", expected, selenomet=False)


def test_4aqd():
    expected = {
        "copies": 1,
        "proteins": [
            {
                "sequence": (
                    "RSEDDIIIATKNGKVRGMNLTVFGGTVTAFLGIPYAQPPLGRLRFKKPQSLTKWSDIWNA"
                    "TKYANSCCQNIDQSFPGFHGSEMWNPNTDLSEDCLYLNVWIPAPKPKNATVLIWIYGGGF"
                    "QTGTSSLHVYDGKFLARVERVIVVSMNYRVGALGFLALPGNPEAPGNMGLFDQQLALQWV"
                    "QKNIAAFGGNPKSVTLFGESAGAASVSLHLLSPGSHSLFTRAILQSGSFNAPWAVTSLYE"
                    "ARNRTLNLAKLTGCSRENETEIIKCLRNKDPQEILLNEAFVVPYGTPLSVNFGPTVDGDF"
                    "LTDMPDILLELGQFKKTQILVGVNKDEGTAFLVYGAPGFSKDNNSIITRKEFQEGLKIFF"
                    "PGVSEFGKESILFHYTDWVDDQRPENYREALGDVVGDYNFICPALEFTKKFSEWGNNAFF"
                    "YYFEHRSSKLPWPEWMGVMHGYEIEFVFGLPLERRDNYTKAEEILSRSIVKRWANFAKYG"
                    "NPNETQNNSTSWPVFKSTEQKYLTLNTESTRIMTKLRAQQCRFWTSFFPKV"
                ),
                "stoichiometry": 2,
                "modifications": [],
            }
        ],
        "rnas": [],
        "dnas": [],
        "carbs": [
            {"codes": {"NAG": 2}, "stoichiometry": 2},
            {"codes": {"FUL": 1, "NAG": 2}, "stoichiometry": 6},
            {"codes": {"MAN": 1, "NAG": 2}, "stoichiometry": 1},
        ],
        "ligands": [
            {"code": "BAL", "stoichiometry": 2},
            {"code": "NAG", "stoichiometry": 6},
            {"code": "PG4", "stoichiometry": 2},
            {"code": "PEG", "stoichiometry": 2},
        ],
        "buffers": ["EDO", "UNX", "CL", "GLY"],
    }
    _test_contents("4aqd", expected, selenomet=False)


def test_1vjr():
    expected = {
        "copies": 1,
        "proteins": [
            {
                "sequence": (
                    "MGSDKIHHHHHHVLDKIELFILDMDGTFYLDDSLLPGSLEFLETLKEKNKRFVFFTNNSS"
                    "LGAQDYVRKLRNMGVDVPDDAVVTSGEITAEHMLKRFGRCRIFLLGTPQLKKVFEAYGHV"
                    "IDEENPDFVVLGFDKTLTYERLKKACILLRKGKFYIATHPDINCPSKEGPVPDAGSIMAA"
                    "IEASTGRKPDLIAGKPNPLVVDVISEKFGVPKERMAMVGDRLYTDVKLGKNAGIVSILVL"
                    "TGETTPEDLERAETKPDFVFKNLGELAKAVQ"
                ),
                "stoichiometry": 1,
                "modifications": ["M->MSE"],
            }
        ],
        "rnas": [],
        "dnas": [],
        "carbs": [],
        "ligands": [],
        "buffers": ["NI", "CL"],
    }
    _test_contents("1vjr", expected, selenomet=True)


def test_1cag():
    expected = {
        "copies": 3,
        "proteins": [
            {
                "sequence": "PPGPPGPPGPPGPPAPPGPPGPPGPPGPPG",
                "stoichiometry": 1,
                "modifications": [
                    "2->HYP",
                    "5->HYP",
                    "8->HYP",
                    "11->HYP",
                    "14->HYP",
                    "17->HYP",
                    "20->HYP",
                    "23->HYP",
                    "26->HYP",
                    "29->HYP",
                ],
            }
        ],
        "rnas": [],
        "dnas": [],
        "carbs": [],
        "ligands": [],
        "buffers": ["ACY"],
    }
    contents = _test_contents("1cag", expected, selenomet=False)
    polymer = contents.proteins[0]
    codes = polymer.residue_codes()
    ppg = ["PRO", "HYP", "GLY"]
    ppa = ["PRO", "HYP", "ALA"]
    assert codes == ppg * 4 + ppa + ppg * 5


def test_1iha():
    expected = {
        "copies": 2,
        "proteins": [],
        "rnas": [
            {
                "sequence": "GCUUCGGCU",
                "stoichiometry": 1,
                "modifications": ["9->BRU"],
            }
        ],
        "dnas": [],
        "carbs": [],
        "ligands": [{"code": "RHD", "stoichiometry": 1}],
        "buffers": ["CL"],
    }
    _test_contents("1iha", expected, selenomet=False)


def test_3ue7():
    expected = {
        "copies": 1,
        "proteins": [
            {
                "sequence": "TTCCPSIVARSNFNACRLPGTPEALCATYTGCIIIPGATCPGDYAN",
                "stoichiometry": 1,
                "modifications": [
                    "T->DTH",
                    "C->DCY",
                    "P->DPR",
                    "S->DSN",
                    "I->DIL",
                    "V->DVA",
                    "A->DAL",
                    "R->DAR",
                    "N->DSG",
                    "F->DPN",
                    "L->DLE",
                    "E->DGL",
                    "Y->DTY",
                    "D->DAS",
                ],
            },
            {
                "sequence": "TTCCPSIVAKSNFNACRLPGTPEALCATYTGCIIIPGATCPGDYAN",
                "stoichiometry": 1,
                "modifications": [],
            },
        ],
        "rnas": [],
        "dnas": [],
        "carbs": [],
        "ligands": [],
        "buffers": [],
    }
    _test_contents("3ue7", expected, selenomet=False)


def test_polymer_weight():
    polymer = Polymer("GG", polymer_type=PolymerType.PROTEIN)
    assert polymer.weight() == approx(132.12, abs=0.01)
