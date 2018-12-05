from nose.tools import eq_

from mhctools import NetMHCII
from mhcnames import normalize_allele_name

def test_netmhcii_DRB():
    alleles = [normalize_allele_name("HLA-DRB1*01:01")]
    ii_predictor = NetMHCII(
        alleles=alleles)
    fasta_dictionary = {
        "SMAD4-001": "PAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGT",
        "TP53-001": "SQAMDDLMLSPDDIEQWFTED"
    }
    binding_predictions = ii_predictor.predict_subsequences(
        sequence_dict=fasta_dictionary,
        peptide_lengths=[15, 16])

    unique_lengths = {x.length for x in binding_predictions}
    eq_(unique_lengths, {15, 16})

    unique_alleles = {x.allele for x in binding_predictions}
    eq_(unique_alleles, {"HLA-DRA1*01:01-DRB1*01:01"})

    # length of "PAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGT" is 34
    # length of "SQAMDDLMLSPDDIEQWFTED" is 21
    # Expect (34-15+1) + (34-16+1) + (21-15+1) + (21-16+1) = 52 entries
    assert len(binding_predictions) == 52, \
        "Expected 52 epitopes from %s" % (binding_predictions,)

def test_netmhcii_alpha_beta():
    alleles = [normalize_allele_name("HLA-DPA1*01:03-DPB1*02:01")]
    ii_predictor = NetMHCII(
        alleles=alleles)
    fasta_dictionary = {
        "SMAD4-001": "PAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGT",
        "TP53-001": "SQAMDDLMLSPDDIEQWFTED"
    }

    binding_predictions = ii_predictor.predict_subsequences(
        sequence_dict=fasta_dictionary,
        peptide_lengths=[15, 16])
    unique_lengths = {x.length for x in binding_predictions}
    eq_(unique_lengths, {15, 16})

    unique_alleles = {x.allele for x in binding_predictions}
    eq_(unique_alleles, {"HLA-DPA1*01:03-DPB1*02:01"})

    # length of "PAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGT" is 34
    # length of "SQAMDDLMLSPDDIEQWFTED" is 21
    # Expect (34-15+1) + (34-16+1) + (21-15+1) + (21-16+1) = 52 entries
    assert len(binding_predictions) == 52, \
        "Expected 52 epitopes from %s" % (binding_predictions,)

def test_netmhcii_multiple_alleles():
    alleles = [
        normalize_allele_name("HLA-DPA1*01:03-DPB1*02:01"),
        normalize_allele_name("HLA-DQA1*01:04-DQB1*05:03")
    ]
    ii_predictor = NetMHCII(
        alleles=alleles)
    fasta_dictionary = {
        "TP53-001": "SQAMDDLMLSPDDIEQWFTED"
    }
    binding_predictions = ii_predictor.predict_subsequences(
        sequence_dict=fasta_dictionary,
        peptide_lengths=[15, 16])

    unique_lengths = {x.length for x in binding_predictions}
    eq_(unique_lengths, {15, 16})

    unique_alleles = {x.allele for x in binding_predictions}
    eq_(unique_alleles, {
        "HLA-DPA1*01:03-DPB1*02:01",
        "HLA-DQA1*01:04-DQB1*05:03"
    })

    # length of "SQAMDDLMLSPDDIEQWFTED" is 21
    # Expect 2 * ((21-15+1) + (21-16+1)) = 26 entries
    assert len(binding_predictions) == 26, \
        "Expected 39 epitopes from %s" % (binding_predictions,)

def test_netmhcii_mouse():
    alleles = [normalize_allele_name("H2-IAb")]
    ii_predictor = NetMHCII(alleles=alleles)
    fasta_dictionary = {
        "SMAD4-001": "PAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGT",
        "TP53-001": "SQAMDDLMLSPDDIEQWFTED"
    }
    binding_predictions = ii_predictor.predict_subsequences(
        sequence_dict=fasta_dictionary,
        peptide_lengths=[15, 16])

    unique_lengths = {x.length for x in binding_predictions}
    eq_(unique_lengths, {15, 16})

    unique_alleles = {x.allele for x in binding_predictions}
    eq_(unique_alleles, {"H-2-IAb"})

    # length of "PAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGT" is 34
    # length of "SQAMDDLMLSPDDIEQWFTED" is 21
    # Expect (34-15+1) + (34-16+1) + (21-15+1) + (21-16+1) = 52 entries
    assert len(binding_predictions) == 52, \
        "Expected 52 epitopes from %s" % (binding_predictions,)
