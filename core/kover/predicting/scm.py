"""
    Kover: Learn interpretable computational phenotyping models from k-merized genomic data
    Copyright (C) 2015  Alexandre Drouin

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
import logging
import os

from ..utils import _fasta_to_sequences


def load_model(model_fasta_path):
    """
    Load a SCM model from its fasta file and return it

    """
    # Load the model from the fasta format
    kmers, headers = _fasta_to_sequences(model_fasta_path, return_headers=True)

    # Determine the type of the model (conjunction or disjunction)
    model_type = "disjunction" if "disjunction" in headers[0].lower() else "conjunction"

    # Load the rules and determine their type
    model = zip(kmers, ["presence" if "presence" in header.lower() else "absence" for header in headers])

    # Check that the k-mers are all of the same length and save that length
    kmer_len = len(kmers[0])
    for k in kmers:
        if len(k) != kmer_len:
            raise ValueError("The k-mers in the model don't have the same length.")

    return model_type, model, kmer_len


def predict_on_genome(model_type, model_rules, model_kmer_len, contig_sequences, contig_headers):
    """
    The model is a model object and contigs is a list of contigs

    We return which k-mer was found in which contig alon with the final prediction

    """
    # Find out which model k-mers are in the genome
    kmer_in_genome = {k: [] for k in zip(*model_rules)[0]}
    for header, sequence in zip(contig_headers, contig_sequences):
        for i in xrange(len(sequence) - model_kmer_len):
            kmer = sequence[i : i + model_kmer_len]
            if kmer_in_genome.has_key(kmer):
                kmer_in_genome[kmer].append(dict(contig=header, start=i, end=i + model_kmer_len))

    # Now that we know which of the model k-mers are in the genome, do the logic
    rule_outcomes = []
    rule_kmer_hits = []
    for kmer, rule_type in model_rules:
        if len(kmer_in_genome[kmer]) > 0:
            rule_outcomes.append(True if rule_type == "presence" else False)
            rule_kmer_hits.append(kmer_in_genome[kmer])

    # TODO: As a minimum viable product, we will just return the prediction,
    #       but future versions will include explanations derived from the model
    #       For example:
    #       If it is a disjunction, return which rules were True as the explanation
    #       If it is a conjunction, return either "All rules satisfied, or rules X violated"
    #       Return where the k-mers were found (which contig and position)
    # I believe that producing an interpretable explanation is the job of the interface
    # Make sure it is like text and maybe summary tables so it will look super cool.
    if model_type == "conjunction":
        return int(all(rule_outcomes)), None
    else:
        return int(any(rule_outcomes)), None


def parse_genomic_data_file(path):
    """
    Reads the TSV file and extract the path to all genomes for which predictions
    must be computed.

    Returns a list of pairs of (genome_id, path)

    """
    genome_data = [l.strip().split("\t") for l in open(path, "r")]
    for genome_id, genome_path in genome_data:
        if not os.path.exists(genome_path):
            raise RuntimeError("File not found {}".format(genome_path))
    return genome_data


def predict(genomic_data_file, genomic_data_type, model_fasta_file, progress_callback):
    """
    Compute predictions on a set of genomes using a model

    """
    if progress_callback is None:
        progress_callback = lambda x, y: None

    if genomic_data_type != "fasta":
        raise ValueError("Only fasta is supported at the moment.")

    logging.debug("Parsing model")
    model_type, model_rules, model_kmer_len = load_model(model_fasta_file)

    results = {}
    logging.debug("Parsing the genomic data file")
    genome_info = parse_genomic_data_file(genomic_data_file)
    logging.debug("Found {} genomes".format(len(genome_info)))
    for i, (genome_id, genome_path) in enumerate(genome_info):
        logging.debug("Reading genome {} from path {}".format(genome_id, genome_path))
        contig_sequences, contig_headers = _fasta_to_sequences(genome_path, return_headers=True)  # TODO: modify so that it returns contig names
        prediction, explanation = predict_on_genome(model_type, model_rules, model_kmer_len, contig_sequences, contig_headers)
        results[genome_id] = dict(prediction=prediction, explanation=explanation)
        progress_callback("Predicting on genomes", 1.0 * (i + 1) / len(genome_info))
    return results
