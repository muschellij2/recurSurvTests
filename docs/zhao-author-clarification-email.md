# Draft email to the Zhao et al. authors

**Subject:** Clarification of Equation 6 variance convention and Table 2 scaling in Zhao et al. (2020)

Dear Drs. Zhao, Zhang, LaValley, Massaro, Lunetta, and Chang,

I am preparing an open-source R implementation and reproducibility report for
the extended recurrent-event rank tests in your paper, “Extended Rank Tests for
Analyzing Recurrent Event Data” (*Statistics in Biopharmaceutical Research*,
2020).

Could you please clarify two details that affect numerical reproduction?

1. **Equation 6 residual denominator.** In our reading, the residual appears to
   use a risk-set denominator specific to the subject's treatment group. A
   pooled-risk implementation is also natural because the score uses the pooled
   Wang--Chang weighted risk set. Which denominator did you use in the
   simulations and published Tables 1–3? If the group-specific denominator is
   intended, should it be evaluated before or after applying the LR/GB/PP time
   weight? Any source code or pseudocode for the subject residual and robust
   variance would be extremely helpful.

2. **Table 2 scaling.** The table caption appears to indicate “(*100),” but the
   printed values (for example 27, 25, and 26) seem consistent with percentages
   reported as “(*1000),” giving probabilities 0.027, 0.025, and 0.026. Is
   “*100” a typographical error, or should the printed values be divided by
   100? We currently retain the printed values and treat them as thousandths
   only as an explicitly flagged sensitivity interpretation.

For context, our implementation keeps the pooled-risk and Equation 6 paths as
separate options and reports Monte Carlo uncertainty rather than claiming exact
published-table reproduction. We would be glad to update the manuscript and
software to match your clarification and to acknowledge any recommended code
or numerical benchmark.

Thank you for your time and for developing these methods.

Best regards,

John Muschelli
