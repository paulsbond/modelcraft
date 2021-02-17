from modelcraft.contents import AsuContents


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
    contents = AsuContents()
    contents.add_from_pdbid("1o6a")
    actual = contents.components_json()
    assert actual == expected


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
    contents = AsuContents()
    contents.add_from_pdbid("4gxy")
    actual = contents.components_json()
    assert actual == expected


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
    contents = AsuContents()
    contents.add_from_pdbid("6as7")
    actual = contents.components_json()
    assert actual == expected
