"""
CPCA: FDA's Carcinogenic Potency Categorization Approach for nitrosamine drug
substance-related impurities (NDSRIs).

    >>> from cpca import score_smiles, format_report
    >>> result = score_smiles("CN(N=O)C1(c2ccccc2Cl)CCCCC1=O")
    >>> result.category, result.ai_ng_day
    (5, 1500.0)
    >>> print(format_report(result, "N-nitroso-ketamine"))
"""

from .policy import DEFAULT_RULESET, PRESET_NAMES, PRESETS, Ruleset
from .report import format_report
from .scorer import CPCAResult, Factor, SiteResult, score_smiles

__all__ = ["score_smiles", "format_report", "CPCAResult", "SiteResult", "Factor",
           "Ruleset", "PRESETS", "PRESET_NAMES", "DEFAULT_RULESET"]
