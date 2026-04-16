# Algorithm Description and Mathematical Formulation

## Introduction

DiffPrimer addresses the "Exclusive Region Identification Problem" in comparative genomics. The objective is to identify genomic subsequences $S$ present in a target taxonomic sets $\mathcal{T}$ that satisfy a uniqueness constraint with respect to a background set $\mathcal{B}$, and subsequently design Polymerase Chain Reaction (PCR) primers $P$ that maximize specificity.

The pipeline integrates three computational paradigms:  
1.  **Alignment-Free Sequence Comparison**: Using $k$-mer multisets for region discovery.
2.  **Bit-Parallel Approximate String Matching**: Myers' algorithm for global homology filtering.
3.  **Semiglobal Dynamic Programming**: For precise local primer verification.

---

## Exclusive Region Identification (Alignment-Free)

Let $\Sigma = \{A, C, G, T\}$ be the alphabet of DNA nucleotides. A genome $G$ is a sequence $G \in \Sigma^*$.
Let $w$ denote a $k$-mer (substring of length $k$). The multiset of $k$-mers in $G$ is defined as:

$$
\mathcal{K}_G = \{ (w, c_w) \mid w \in \Sigma^k \}
$$

where $c_w$ is the count of occurrences of $w$ in $G$.

### Target and Background Sets

Let the Reference Genome be $R$. We define the set of **Target Unique $k$-mers**, $\mathcal{U}_R$, subject to an abundance threshold $\tau$ (experimentally set to $\tau=1$ for single-copy markers):

$$
\mathcal{U}_R = \{ w \in \Sigma^k \mid (w, c_w) \in \mathcal{K}_R \land c_w \le \tau \}
$$

Let the set of Background Genomes be $\mathcal{B} = \{B_1, B_2, \dots, B_N\}$. The **Background Universe** $\Omega_B$ is the union of all background $k$-mers:

$$
\Omega_B = \bigcup_{i=1}^N \{ w \mid (w, c_w) \in \mathcal{K}_{B_i} \}
$$

### Diagnostic Marker Set

The set of candidate diagnostic markers $\mathcal{D}$ is the set difference:

$$
\mathcal{D} = \mathcal{U}_R \setminus \Omega_B
$$

**Contig Assembly:**
Continuous regions are assembled from $\mathcal{D}$. Two $k$-mers $w_i, w_j \in \mathcal{D}$ (at positions $p_i, p_j$ in $R$) are merged into a region $S_{region}$ if $p_j = p_i + 1$. This yields a set of candidate amplicons $\mathcal{S} = \{S_1, S_2, \dots\}$.

---

## Specificity Verification Algorithm

For each candidate region $S \in \mathcal{S}$ and designed primer pair $(P_{fwd}, P_{rev})$, we verify specificity against every background genome $B \in \mathcal{B}$. This is a two-stage process.

### Stage I: Global Homology Filter (Myers' Bit-Vector Algorithm)

To efficiently identify if $S$ appears in $B$ with edit distance $d < d_{max}$, we employ the **Myers Bit-Vector Algorithm** (Myers, 1999). This algorithm exploits the bit-parallelism of modern processors to compute a column of the Dynamic Programming (DP) matrix in $O(m/w)$ operations, where $m=|S|$ and $w=64$ (machine word size).

**State Variables:**
Let $D_{i,j}$ be the edit distance between $S[1..i]$ and $B[1..j]$.
We define delta variables for the vertical differences in the DP matrix:

$$
\Delta V_{i,j} = D_{i,j} - D_{i-1,j} \in \{-1, 0, +1\}
$$

These are encoded using two bit-vectors, $VP$ (Vertical Positive) and $VN$ (Vertical Negative), such that the $i$-th bit is set if $\Delta V_{i,j} = +1$ or $-1$ respectively.

**Complexity:** $O(\lceil \frac{m}{64} \rceil \cdot n)$, where $n=|B|$. This provides a $64\times$ speedup over standard Levenshtein calculation.

### Stage II: Local Verification (Semiglobal Alignment)

When Stage I identifies a "Hit" (substring $B'$ in $B$ such that $\text{Lev}(S, B') \le d_{max}$), we must verify if the **primers** bind to $B'$. Due to insertions/deletions (indels), the primers' positions in $B'$ are shifted relative to $S$.

We solve this using **Semiglobal Alignment** (Needleman-Wunsch variant with free end-gaps).
Let $S$ be the query (Region) and $W$ be the target window extracted from $B$ around the hit.

**Score Optimization:**
Maximized Score matrix $H$:

$$
H_{i,j} = \max \begin{cases} H_{i-1, j-1} + \sigma(S[i], W[j]) \\ H_{i-1, j} + \gamma_{\text{gap}} \\ H_{i, j-1} + \gamma_{\text{gap}} \end{cases}
$$

where $\sigma(a, b) = 1$ if $a=b$ else $-1$, and $\gamma_{\text{gap}} = -1$.
Initialization: $H_{0,j} = 0, H_{i,0} = i \cdot \gamma$ (Semiglobal in target, global in query).

**Coordinate Mapping:**
An optimal alignment trace $\mathcal{T}$ is computed. We map the primer interval $[x_{start}, x_{end}]$ from $S$ to $[y_{start}, y_{end}]$ in $W$ via the trace function $\tau$:

$$
P_{mapped} = \tau(\mathcal{T}, P_{original})
$$

This accounts for all indels between the primer site and the region boundaries.

**Final Specificity Test (Positional Mismatch Penalty):**
Instead of a simple local edit distance, we apply a **positional weight** $\omega(x)$ to each mismatch, insertion, or deletion relative to the mapped primer sequence $P_{seq}$ of length $L_P$, favoring high specificity at the critical 3' end.
Let $x \in [1, L_P]$ be the nucleotide position along the primer (5' to 3').

$$
\omega(x) = \begin{cases} 3.0 & \text{if } x \ge L_P - 5 \text{ (local 3' region)} \\ 1.0 & \text{otherwise} \end{cases}
$$

The positional mismatch score $\delta_{pos}$ is accumulated during the alignment traversal over all editing operations:

$$
\delta_{pos} = \sum_{\text{op} \in \{\text{Subst}, \text{Ins}, \text{Del}\}} \omega(\text{pos}(\text{op}))
$$

If $\delta_{pos} < T_{mismatch}$ (where $T_{mismatch}$ is a user-configurable threshold, default 7.0) for **both** forward and reverse primers, the amplicon is flagged as **NonSpecific** (potential false positive). If the mismatch score exceeds the threshold for either primer, it is safely filtered as **Specific**.

---

## Statistical Scoring

Homology is reported as:

$$
\text{Similarity Score} = \left( 1 - \frac{d_{min}}{L_{amplicon}} \right) \times 100
$$

where $d_{min}$ is the minimum Levenshtein distance found across all background genomes.

This methodology ensures that DiffPrimer minimizes both False Positives (non-specific binding) and simplifies the computational burden of scanning gigobase-sized genomes.
