"""crispr_correct — forward-looking import alias for the `crispr_ambiguous_mapping` package.

§4.11: the project has three names historically — PyPI `crispr-ambiguous-mapping`,
repo `CRISPR-Correct`, import `crispr_ambiguous_mapping`. This shim unifies
future imports under a single name (`crispr_correct`) matching the repo.

Use this for new code:

```python
import crispr_correct as cc
result = cc.map_fastq(library, fastq_r1_fns=["R1.fq.gz"], config=cc.ParsingConfig(...))
```

The existing `crispr_ambiguous_mapping` package stays available through 0.1.x
for existing drivers; it'll be deprecated in 0.2.0 and removed in 0.3.0.
"""
from crispr_ambiguous_mapping import mapping, utility, models, processing, visualization, quality_control, postprocessing  # noqa: F401
from crispr_ambiguous_mapping.api import map_fastq, count, alleles, ParsingConfig, save, load  # noqa: F401
from crispr_ambiguous_mapping.models.mapping_models import MatchTier  # noqa: F401
