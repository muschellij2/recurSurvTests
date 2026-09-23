# Published benchmark provenance

`zhao-2020-published-tables.csv` transcribes the **New method** columns
(LR, PP, GB) of Tables 1–3, pages 93–95, in Zhao et al. (2020),
*Extended Rank Tests for Analyzing Recurrent Event Data*,
[doi:10.1080/19466315.2019.1601596](https://doi.org/10.1080/19466315.2019.1601596).
Source: the locally supplied publisher PDF, `Extended Rank Tests for Analyzing
Recurrent Event Data.pdf`. Table 1 is the empirical SD of the standardized
statistic; Tables 2–3 are rejection probabilities. OLR and Jung–Jeong columns
are not benchmarks for the implementation assessed here.

**Table 2 has a printed scale inconsistency.** The PDF caption says `(*100)`
but the New method cells print integers such as `27`, `26`, and `25`, without
decimal points (confirmed visually on page 94). Section 3.2 specifies a
one-sided 0.025 level and describes these results as controlling type I error.
We interpret 27 as 0.027 (2.7%), i.e. divide the printed integers by 1000.
This is an explicit inference from the text, not an author-confirmed correction.
Both the printed value and interpreted probability are retained in the CSV.
Literal division by 100 would give 27%, inconsistent with the accompanying text.
Table 2 comparisons are conditional on this interpretation.

`log_time_shift` is 0 under the null and 0.25 or 0.5 for Table 3. For the
Weibull with shape 2, the printed changes in `-log(lambda)` are 0.5 and 1,
which correspond to these same log-time shifts. Other distributions use
changes of 0.25 and 0.5 in `-log(lambda)` or `mu`.

All benchmark values have effective precision 0.001. The report allows a
half rounding unit (0.0005) in comparisons. Paper Monte Carlo uncertainty
uses 100,000 replicates for Tables 1–2 and 10,000 for Table 3.
