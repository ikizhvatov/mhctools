# mhctools
Python interface to running command-line and web-based MHC binding predictors. 

## Example

```python
from mhctools import NetMHCpan
# Run NetMHCpan for alleles HLA-A*01:01 and HLA-A*02:01
predictor = NetMHCpan(alleles=["A*02:01", "hla-a0101"])

# scan the short proteins 1L2Y and 1L3Y for epitopes
protein_sequences = (
  {"1L2Y": "NLYIQWLKDGGPSSGRPPPS",
  "1L3Y": "ECDTINCERYNGQVCGGPGRGLCFCGKCRCHPGFEGSACQA"
}

epitope_collection = predictor.predict(protein_sequences)

# flatten binding predictions into a Pandas DataFrame
df = epitope_collection.dataframe()

# epitope collection is sorted by percentile rank
# of binding predictions
strongest_predicted_binder = epitope_collection[0]
```
## API

The following models are available in `mhctools`: 
* `NetMHCpan`: requires locally installed version of [NetMHCpan](http://www.cbs.dtu.dk/services/NetMHCpan/)
* `NetMHCcons`: requires locally installed version of [NetMHCcons](http://www.cbs.dtu.dk/services/NetMHCcons/)
* `IedbMhcClass1`: Uses IEDB's REST API for class I binding predictions.
* `IedbMhcClass2`: Uses IEDB's REST API for class II binding predictions.
* `RandomBindingPredictor`: Creates binding predictions with random IC50 and percentile rank values.

Every model is constructed with an `alleles` argument specifying the HLA type for which to make predictions. Predictions are generated by calling the `predict` method with a dictionary mapping sequence IDs or names to amino acid sequences.
