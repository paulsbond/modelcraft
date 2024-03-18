from pytest import approx
from modelcraft.contents import AsuContents, Polymer, PolymerType


def _test_contents(entry: str, expected_json: list, selenomet: bool):
    contents = AsuContents.from_pdbe(entry)
    assert contents.to_json() == expected_json
    assert contents.is_selenomet() == selenomet
    return contents


def test_1o6a():
    expected = {
        "proteins": [
            {
                "sequence": "SETRKTEVPSDKLELLLDIPLKVTVELGRTRMTLKRVLEMIHGSIIELDKLTGEPVDILVNGKLIARGEVVVIDENFGVRITEIVSPKERLELLNE",
                "copies": 2,
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
        "proteins": [],
        "rnas": [
            {
                "sequence": "GGCGGCAGGUGCUCCCGACCCUGCGGUCGGGAGUUAAAAGGGAAGCCGGUGCAAGUCCGGCACGGUCCCGCCACUGUGACGGGGAGUCGCCCCUCGGGAUGUGCCACUGGCCCGAAGGCCGGGAAGGCGGAGGGGCGGCGAGGAUCCGGAGUCAGGAAACCUGCCUGCCGUC",
                "copies": 1,
                "modifications": ["1->GTP", "172->CCC"],
            }
        ],
        "dnas": [],
        "carbs": [],
        "ligands": [
            {"code": "B1Z", "copies": 2},
            {"code": "IRI", "copies": 7},
        ],
        "buffers": ["MG"],
    }
    _test_contents("4gxy", expected, selenomet=False)


def test_6as7():
    expected = {
        "proteins": [
            {
                "sequence": "DEEQVFHFYWLDAYEDQYNQPGVVFLFGKVWIESAETHVSCCVMVKNIERTLYFLPREMKIDLNTGKETGTPISMKDVYEEFDEKIATKYKIMKFKSKPVEKNYAFEIPDVPEKSEYLEVKYSAEMPQLPQDLKGETFSHVFGTNTSSLELFLMNRKIKGPCWLEVKSPQLLNQPVSWCKAEAMALKPDLVNVIKDVSPPPLVVMAFSMKTMQNAKNHQNEIIAMAALVHHSFALDKAAPKPPFQSHFCVVSKPKDCIFPYAFKEVIEKKNVKVEVAATERTLLGFFLAKVHKIDPDIIVGHNIYGFELEVLLQRINVCKAPHWSKIGRLKRSNMPKLGGRSGFGERNATCGRMICDVEISAKELIRCKSYHLSELVQQILKTERVVIPMENIQNMYSESSQLLYLLEHTWKDAKFILQIMCELNVLPLALQITNIAGNIMSRTLMGGRSERNEFLLLHAFYENNYIVPDKQIFRKPQQKLGDEDEEIDGDTNKYKKGRKKAAYAGGLVLDPKVGFYDKFILLLDFNSLYPSIIQEFNICFTTVQRVASEAQKVTEDGEQEQIPELPDPSLEMGILPREIRKLVERRKQVKQLMKQQDLNPDLILQYDIRQKALKLTANSMYGCLGFSYSRFYAKPLAALVTYKGREILMHTKEMVQKMNLEVIYGDTDSIMINTNSTNLEEVFKLGNKVKSEVNKLYKLLEIDIDGVFKSLLLLKKKKYAALVVEPTSDGNYVTKQELKGLDIVRRDWCDLAKDTGNFVIGQILSDQSRDTIVENIQKRLIEIGENVLNGSVPVSQFEINKALTKDPQDYPDKKSLPHVHVALWINSQGGRKVKAGDTVSYVICQDGSNLTASQRAYAPEQLQKQDNLTIDTQYYLAQQIHPVVARICEPIDGIDAVLIATWLGLDPTQFRVHHYHKDEEN",
                "copies": 1,
                "modifications": [],
            }
        ],
        "rnas": [],
        "dnas": [
            {"sequence": "GCCTGGAGCGC", "copies": 1, "modifications": []},
            {"sequence": "AGGCGCTCCAGGC", "copies": 1, "modifications": []},
        ],
        "carbs": [],
        "ligands": [{"code": "DCP", "copies": 1}],
        "buffers": ["MG", "CO"],
    }
    _test_contents("6as7", expected, selenomet=False)


def test_4aqd():
    expected = {
        "proteins": [
            {
                "sequence": "RSEDDIIIATKNGKVRGMNLTVFGGTVTAFLGIPYAQPPLGRLRFKKPQSLTKWSDIWNATKYANSCCQNIDQSFPGFHGSEMWNPNTDLSEDCLYLNVWIPAPKPKNATVLIWIYGGGFQTGTSSLHVYDGKFLARVERVIVVSMNYRVGALGFLALPGNPEAPGNMGLFDQQLALQWVQKNIAAFGGNPKSVTLFGESAGAASVSLHLLSPGSHSLFTRAILQSGSFNAPWAVTSLYEARNRTLNLAKLTGCSRENETEIIKCLRNKDPQEILLNEAFVVPYGTPLSVNFGPTVDGDFLTDMPDILLELGQFKKTQILVGVNKDEGTAFLVYGAPGFSKDNNSIITRKEFQEGLKIFFPGVSEFGKESILFHYTDWVDDQRPENYREALGDVVGDYNFICPALEFTKKFSEWGNNAFFYYFEHRSSKLPWPEWMGVMHGYEIEFVFGLPLERRDNYTKAEEILSRSIVKRWANFAKYGNPNETQNNSTSWPVFKSTEQKYLTLNTESTRIMTKLRAQQCRFWTSFFPKV",
                "copies": 2,
                "modifications": [],
            }
        ],
        "rnas": [],
        "dnas": [],
        "carbs": [
            {"codes": {"NAG": 2}, "copies": 2},
            {"codes": {"FUL": 1, "NAG": 2}, "copies": 6},
            {"codes": {"MAN": 1, "NAG": 2}, "copies": 1},
        ],
        "ligands": [
            {"code": "BAL", "copies": 2},
            {"code": "NAG", "copies": 6},
            {"code": "PG4", "copies": 2},
            {"code": "PEG", "copies": 2},
        ],
        "buffers": ["EDO", "CL", "GLY"],
    }
    _test_contents("4aqd", expected, selenomet=False)


def test_1vjr():
    expected = {
        "proteins": [
            {
                "sequence": "MGSDKIHHHHHHVLDKIELFILDMDGTFYLDDSLLPGSLEFLETLKEKNKRFVFFTNNSSLGAQDYVRKLRNMGVDVPDDAVVTSGEITAEHMLKRFGRCRIFLLGTPQLKKVFEAYGHVIDEENPDFVVLGFDKTLTYERLKKACILLRKGKFYIATHPDINCPSKEGPVPDAGSIMAAIEASTGRKPDLIAGKPNPLVVDVISEKFGVPKERMAMVGDRLYTDVKLGKNAGIVSILVLTGETTPEDLERAETKPDFVFKNLGELAKAVQ",
                "copies": 1,
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
        "proteins": [
            {
                "sequence": "PPGPPGPPGPPGPPAPPGPPGPPGPPGPPG",
                "copies": 3,
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
        "proteins": [],
        "rnas": [],
        "dnas": [{"sequence": "GCUUCGGCU", "copies": 2, "modifications": ["9->BRU"]}],
        "carbs": [],
        "ligands": [{"code": "RHD", "copies": 2}],
        "buffers": ["CL"],
    }
    _test_contents("1iha", expected, selenomet=False)


def test_3ue7():
    expected = {
        "proteins": [
            {
                "sequence": "TTCCPSIVARSNFNACRLPGTPEALCATYTGCIIIPGATCPGDYAN",
                "copies": 1,
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
                "copies": 1,
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
    polymer = Polymer("GG", 1, polymer_type=PolymerType.PROTEIN)
    assert polymer.weight() == approx(132.12, abs=0.01)
