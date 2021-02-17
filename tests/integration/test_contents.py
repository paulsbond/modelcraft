from modelcraft.contents import AsuContents


def _test_contents(pdbid: str, expected: list):
    contents = AsuContents()
    contents.add_from_pdbid(pdbid)
    actual = contents.components_json()
    assert actual == expected


def test_1o6a():
    expected = [
        {
            "sequence": "SETRKTEVPSDKLELLLDIPLKVTVELGRTRMTLKRVLEMIHGSIIELDKLTGEPVDILVNGKLIARGEVVVIDENFGVRITEIVSPKERLELLNE",
            "type": "PROTEIN",
            "start": 1,
            "copies": 2,
            "modifications": ["MSE"],
        }
    ]
    _test_contents("1o6a", expected)


def test_4gxy():
    expected = [
        {
            "sequence": "GGCGGCAGGUGCUCCCGACCCUGCGGUCGGGAGUUAAAAGGGAAGCCGGUGCAAGUCCGGCACGGUCCCGCCACUGUGACGGGGAGUCGCCCCUCGGGAUGUGCCACUGGCCCGAAGGCCGGGAAGGCGGAGGGGCGGCGAGGAUCCGGAGUCAGGAAACCUGCCUGCCGUC",
            "type": "RNA",
            "start": 1,
            "copies": 1,
            "modifications": ["GTP1", "CCC172"],
        },
        {"code": "B1Z", "copies": 2},
        {"code": "IRI", "copies": 7},
        {"code": "MG", "copies": 2},
    ]
    _test_contents("4gxy", expected)


def test_6as7():
    expected = [
        {
            "sequence": "DEEQVFHFYWLDAYEDQYNQPGVVFLFGKVWIESAETHVSCCVMVKNIERTLYFLPREMKIDLNTGKETGTPISMKDVYEEFDEKIATKYKIMKFKSKPVEKNYAFEIPDVPEKSEYLEVKYSAEMPQLPQDLKGETFSHVFGTNTSSLELFLMNRKIKGPCWLEVKSPQLLNQPVSWCKAEAMALKPDLVNVIKDVSPPPLVVMAFSMKTMQNAKNHQNEIIAMAALVHHSFALDKAAPKPPFQSHFCVVSKPKDCIFPYAFKEVIEKKNVKVEVAATERTLLGFFLAKVHKIDPDIIVGHNIYGFELEVLLQRINVCKAPHWSKIGRLKRSNMPKLGGRSGFGERNATCGRMICDVEISAKELIRCKSYHLSELVQQILKTERVVIPMENIQNMYSESSQLLYLLEHTWKDAKFILQIMCELNVLPLALQITNIAGNIMSRTLMGGRSERNEFLLLHAFYENNYIVPDKQIFRKPQQKLGDEDEEIDGDTNKYKKGRKKAAYAGGLVLDPKVGFYDKFILLLDFNSLYPSIIQEFNICFTTVQRVASEAQKVTEDGEQEQIPELPDPSLEMGILPREIRKLVERRKQVKQLMKQQDLNPDLILQYDIRQKALKLTANSMYGCLGFSYSRFYAKPLAALVTYKGREILMHTKEMVQKMNLEVIYGDTDSIMINTNSTNLEEVFKLGNKVKSEVNKLYKLLEIDIDGVFKSLLLLKKKKYAALVVEPTSDGNYVTKQELKGLDIVRRDWCDLAKDTGNFVIGQILSDQSRDTIVENIQKRLIEIGENVLNGSVPVSQFEINKALTKDPQDYPDKKSLPHVHVALWINSQGGRKVKAGDTVSYVICQDGSNLTASQRAYAPEQLQKQDNLTIDTQYYLAQQIHPVVARICEPIDGIDAVLIATWLGLDPTQFRVHHYHKDEEN",
            "type": "PROTEIN",
            "start": 1,
            "copies": 1,
            "modifications": [],
        },
        {
            "sequence": "GCCTGGAGCGC",
            "type": "DNA",
            "start": 1,
            "copies": 1,
            "modifications": [],
        },
        {
            "sequence": "AGGCGCTCCAGGC",
            "type": "DNA",
            "start": 1,
            "copies": 1,
            "modifications": [],
        },
        {"code": "DCP", "copies": 1},
        {"code": "MG", "copies": 2},
        {"code": "CO", "copies": 2},
    ]
    _test_contents("6as7", expected)


def test_4aqd():
    expected = [
        {
            "sequence": "RSEDDIIIATKNGKVRGMNLTVFGGTVTAFLGIPYAQPPLGRLRFKKPQSLTKWSDIWNATKYANSCCQNIDQSFPGFHGSEMWNPNTDLSEDCLYLNVWIPAPKPKNATVLIWIYGGGFQTGTSSLHVYDGKFLARVERVIVVSMNYRVGALGFLALPGNPEAPGNMGLFDQQLALQWVQKNIAAFGGNPKSVTLFGESAGAASVSLHLLSPGSHSLFTRAILQSGSFNAPWAVTSLYEARNRTLNLAKLTGCSRENETEIIKCLRNKDPQEILLNEAFVVPYGTPLSVNFGPTVDGDFLTDMPDILLELGQFKKTQILVGVNKDEGTAFLVYGAPGFSKDNNSIITRKEFQEGLKIFFPGVSEFGKESILFHYTDWVDDQRPENYREALGDVVGDYNFICPALEFTKKFSEWGNNAFFYYFEHRSSKLPWPEWMGVMHGYEIEFVFGLPLERRDNYTKAEEILSRSIVKRWANFAKYGNPNETQNNSTSWPVFKSTEQKYLTLNTESTRIMTKLRAQQCRFWTSFFPKV",
            "type": "PROTEIN",
            "start": 1,
            "copies": 2,
            "modifications": [],
        },
        {"codes": ["NAG"], "length": 2, "copies": 2},
        {"codes": ["FUL", "NAG"], "length": 3, "copies": 6},
        {"codes": ["MAN", "NAG"], "length": 3, "copies": 1},
        {"code": "BAL", "copies": 2},
        {"code": "NAG", "copies": 6},
        {"code": "PG4", "copies": 2},
        {"code": "EDO", "copies": 19},
        {"code": "CL", "copies": 7},
        {"code": "GLY", "copies": 2},
        {"code": "PEG", "copies": 2},
    ]
    _test_contents("4aqd", expected)
