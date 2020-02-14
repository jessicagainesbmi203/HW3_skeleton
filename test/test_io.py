from smith_waterman.io import read_sequence

def test_read_sequence():
    seq = read_sequence('sequences/prot-0004.fa')
    correct_seq = 'SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM'
    assert seq == correct_seq