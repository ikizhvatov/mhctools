# Copyright (c) 2016-2017. Mount Sinai School of Medicine
# Copyright (c) 2018. Ilya Kizhvatov, Wahsington University School of Medicine, ilya.kizhvatov@wustl.edu
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function, division, absolute_import

from mhcnames import parse_classi_or_classii_allele_name, normalize_allele_name

from .base_commandline_predictor import BaseCommandlinePredictor
from .binding_prediction import BindingPrediction

MIXMHC2PRED_HEADER = "Peptide\tBestAllele\tScore\tScore_perL\tBestCore\tBest_s"

class MixMHC2pred(BaseCommandlinePredictor):
    def __init__(
            self,
            alleles,
            default_peptide_lengths=[9],
            program_name="MixMHC2pred",
            process_limit=-1,
            extra_flags=[]):

        # Usage: MixMHC2pred --input INPUT --output OUTPUT --alleles ALLELES [ALLELES...] [--no_Nterm] 
        #                    [--no_Cterm] [--flat_ws]

        # The --output /dev/stdout is for redirecting the output to console
        BaseCommandlinePredictor.__init__(
            self,
            program_name=program_name,
            alleles=alleles,
            default_peptide_lengths=default_peptide_lengths,
            parse_output_fn=self.parse_output,
            supported_alleles_flag=None,
            peptide_mode_flags=None,
            tempdir_flag=None,
            input_file_flag="--input",
            length_flag=None,
            allele_flag="--alleles",
            extra_flags= ["--output", "/dev/stdout"] + extra_flags,
            process_limit=1)

        # 
    def predict_peptides(self, peptides, peptide_mode=False):
        """ Override the base class function to set the peptide_mode
            flagn _build_command to False """

        return BaseCommandlinePredictor.predict_peptides(
            self,
            peptides=peptides,
            peptide_mode=peptide_mode)

    def parse_output(
        self,
        stdout,
        prediction_method_name="mixmhc2pred",
        sequence_key_mapping=None):
        """ MixMHC2pred format is so different from netMHC family
            that is make sense to have a fully custom parser """

        binding_predictions = []
        seen_header = False
        for l in stdout.split("\n"):
            l = l.strip()

            # ignore empty lines
            if not l:
                continue
            
            # step through lines till the header
            if l.startswith(MIXMHC2PRED_HEADER):
                seen_header = True
                continue
            if not seen_header:
                continue

            # after having seen the header, parse every line
            fields = l.split("\t")
            binding_predictions.append(BindingPrediction(
                peptide=fields[0],
                allele=normalize_allele_name(fields[1]),
                affinity=0,
                percentile_rank=float(fields[2]),
                prediction_method_name=prediction_method_name))

        return binding_predictions

    def _prepare_drb_allele_name(self, parsed_beta_allele):
        """
        Format defined by MixMHC2pred's Alleles_list.txt
        """
        if "DRB" not in parsed_beta_allele.gene:
            raise ValueError("Unexpected allele %s" % parsed_beta_allele)
        return "%s_%s_%s" % (
            parsed_beta_allele.gene,
            parsed_beta_allele.allele_family,
            parsed_beta_allele.allele_code)

    def prepare_allele_name(self, allele_name):
        """
        MixMHC2pred allele format requirements as defined
        by its Alleles_list.txt
         - DRB1_01_01 (for non-alpha/beta pairs)
         - DPA1_01_03__DPB1_02_01 (for alpha and beta pairs)
         TODOL
         - DPA1_01_03__DPA1_02_01__DPB1_03_01 (for some alpha and beta pairs)

        Only human class II alleles are supported
        """
        parsed_alleles = parse_classi_or_classii_allele_name(allele_name)
        if len(parsed_alleles) == 1:
            allele = parsed_alleles[0]
            return self._prepare_drb_allele_name(allele)

        else:
            alpha, beta = parsed_alleles
            if "DRA" in alpha.gene:
                return self._prepare_drb_allele_name(beta)
            return "%s_%s_%s__%s_%s_%s" % (
                alpha.gene,
                alpha.allele_family,
                alpha.allele_code,
                beta.gene,
                beta.allele_family,
                beta.allele_code)
